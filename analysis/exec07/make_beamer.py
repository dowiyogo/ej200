#!/usr/bin/env python3
"""Generate the integral EXEC_09 full/key Beamer decks from analysis products."""

from __future__ import annotations

import argparse
import math
import pathlib
import subprocess

import pandas as pd

N_EVENTS = 2000
KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)


def frame(title: str, body: str) -> str:
    return f"\\begin{{frame}}{{{title}}}\n{body}\n\\end{{frame}}\n"


def image_frame(title: str, filename: str, note: str = "", height: float = 0.78) -> str:
    return frame(
        title,
        rf"\centering\includegraphics[width=\textwidth,height={height}\textheight,keepaspectratio]{{{filename}}}"
        + (rf"\par\scriptsize {note}" if note else ""),
    )


def stats(summary: pd.DataFrame, x_mm: int, group: str) -> pd.Series:
    return summary[(summary.x_beam_mm == x_mm) & (summary.grupo == group)].iloc[0]


def npe_text(row: pd.Series) -> str:
    return f"${row.npe_mean:.1f}\\pm{row.npe_rms / math.sqrt(N_EVENTS):.1f}$"


def time_text(row: pd.Series, field: str, rms_field: str, samples: float) -> str:
    error = row[rms_field] / math.sqrt(samples)
    precision = 5 if error < 0.0005 else 4 if error < 0.005 else 3
    return f"${row[field]:.{precision}f}\\pm{error:.{precision}f}$ ns"


def position_frame(position: pd.Series, summary: pd.DataFrame, localization: pd.DataFrame) -> str:
    x_mm = int(position.x_beam_mm)
    left = stats(summary, x_mm, "end_left_all")
    right = stats(summary, x_mm, "end_right_all")
    top = stats(summary, x_mm, "nearest_top")
    gate = localization.loc[x_mm]
    status = {
        "strict maximum at nearest": "exact geometry",
        "known window-track dip <=15%": "window-track exception",
        "strict statistical tie <=1 sigma": "statistical tie: two nearest",
    }[gate.reason]
    body = rf"""
\begin{{columns}}[T]
  \column{{0.52\textwidth}}
  \includegraphics[width=\textwidth]{{figs/muon_{x_mm}mm_geometry.png}}
  \column{{0.46\textwidth}}
  \includegraphics[width=\textwidth]{{figs/muon_{x_mm}mm_top_profile.png}}
  \tiny
  \begin{{tabular}}{{lr}}
    \toprule Quantity & Value \\ \midrule
    End L $N_{{pe}}$ & {npe_text(left)} \\
    End R $N_{{pe}}$ & {npe_text(right)} \\
    Nearest Top ID & {int(position.nearest_top_id)} ({status}) \\
    Nearest Top $N_{{pe}}$ & {npe_text(top)} \\
    Nearest Top $\langle t\rangle$ &
      {time_text(top, "t_mean_ns", "t_rms_ns", top.npe_mean * N_EVENTS)} \\
    Nearest Top FPT & {time_text(top, "fpt_mean_ns", "fpt_rms_ns", N_EVENTS)} \\
    \bottomrule
  \end{{tabular}}
\end{{columns}}
"""
    return frame(f"Position x={x_mm} mm", body)


def dip_table(profiles: pd.DataFrame) -> str:
    definitions = [
        ("existing_x-650", "Existing -650", 19, 18),
        ("run_A_x-652", "A -652 center", 19, 18),
        ("run_B_x-642", "B -642 midpoint", 19, 18),
        ("run_C1_x-648", "C1 -648 (+4 mm)", 19, 18),
        ("run_C2_x-654_exact_mirror", "C2 -654 mirror", 17, 18),
    ]
    lines = []
    for run, label, high_id, low_id in definitions:
        selected = profiles[profiles.run == run].set_index("channel")
        high, low = selected.loc[high_id], selected.loc[low_id]
        difference = high.npe_mean - low.npe_mean
        error = math.hypot(high.npe_sem, low.npe_sem)
        lines.append(
            rf"{label} & ${difference:.2f}\pm{error:.2f}$ & ${difference/error:.2f}$ \\"
        )
    return "\n".join(lines)


def timing_table(gate: pd.DataFrame) -> str:
    lines = []
    for _, row in gate.iterrows():
        lines.append(
            rf"{int(row.x_beam_mm):+d} & {row.side[0].upper()} & "
            rf"${row.sigma_endtop_ps:.1f}\pm{row.sigma_endtop_error_ps:.1f}$ & "
            rf"${row.sigma_end_only_ps:.1f}\pm{row.sigma_end_only_error_ps:.1f}$ & "
            rf"${row.observed_sigma_ratio:.3f}\pm{row.observed_sigma_ratio_error:.3f}$ \\"
        )
    return "\n".join(lines)


def top_estimate_table(summary: pd.DataFrame) -> str:
    lines = []
    for x_mm in KEY_POSITIONS:
        row = stats(summary, x_mm, "nearest_top")
        npe_sem = row.npe_rms / math.sqrt(N_EVENTS)
        sigma_error = 0.5 * row.sigma_est_ps * npe_sem / row.npe_mean
        lines.append(
            rf"{x_mm:+d} & {npe_text(row)} & ${row.sigma_est_ps:.2f}\pm{sigma_error:.2f}$ \\"
        )
    return "\n".join(lines)


def landau_table(positions: pd.DataFrame, landau: pd.DataFrame) -> str:
    lines = []
    for x_mm in KEY_POSITIONS:
        channel = int(positions.loc[x_mm, "nearest_top_id"])
        row = landau[(landau.x_beam_mm == x_mm) & (landau.channel == channel)].iloc[0]
        lines.append(
            rf"{x_mm:+d} & {channel} & ${row.mpv:.1f}\pm{row.mpv_error:.1f}$ & "
            rf"${row.width:.1f}\pm{row.width_error:.1f}$ & {row.chi2_ndf:.2f} \\"
        )
    return "\n".join(lines)


def velocity_table(velocity_fits: pd.DataFrame) -> str:
    labels = {
        "mean_ns": r"all-photon mean",
        "t50_ns": r"all-photon $t_{50}$",
        "fpt_mean_ns": "mean FPT",
    }
    lines = []
    for estimator in ("mean_ns", "t50_ns", "fpt_mean_ns"):
        row = velocity_fits[(velocity_fits.scope == "combined") & (velocity_fits.estimator == estimator)].iloc[0]
        lines.append(
            rf"{labels[estimator]} & ${row.velocity_cm_ns:.3f}\pm{row.velocity_error_cm_ns:.3f}$ & {row.chi2_ndf:.1f} \\"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Exec-12 pedagogical helpers (T2, T3, T4, T5, T6)
# ---------------------------------------------------------------------------

def readingbox(why: str, axes: str, takeaway: str) -> str:
    """Inline LaTeX block: how to read this plot."""
    return (
        r"\begin{block}{\footnotesize How to read this plot}" + "\n"
        + rf"\scriptsize\textbf{{Why:}}~{why}" + "\n"
        + rf"\par\textbf{{Axes:}}~{axes}" + "\n"
        + rf"\par\textbf{{Takeaway:}}~{takeaway}" + "\n"
        + r"\end{block}" + "\n"
    )


def image_frame_wb(
    title: str,
    filename: str,
    note: str = "",
    why: str = "",
    axes: str = "",
    takeaway: str = "",
    height: float = 0.63,
) -> str:
    """image_frame with an optional readingbox appended."""
    body = (
        rf"\centering\includegraphics[width=\textwidth,height={height}\textheight,"
        rf"keepaspectratio]{{{filename}}}"
        + (rf"\par\scriptsize {note}\vspace{{2pt}}" if note else "")
        + (readingbox(why, axes, takeaway) if why else "")
    )
    return frame(title, body)


def window_track_intro_slide() -> str:
    """T2: one-slide definition of the window-track collection effect."""
    return frame(
        "Window-track collection effect: mechanism",
        r"""
\begin{columns}[T]
\column{0.52\textwidth}
\begin{itemize}\small
  \item Top SiPMs couple through index-matched windows
        ($6\times6$ mm$^2$) in the +Y Mylar face.
  \item First-pass photon collection depends on the solid angle
        subtended by the window from the muon track (geometric
        working hypothesis, consistent with five-run pattern).
  \item Empirically: a track at $\pm2$ mm from the window centre
        collects $\approx9.8$\% more $N_{pe}$ ($\approx12\,\sigma$);
        null at centre and at the midpoint; symmetric in mirror (C2).
  \item Scan positions with $x\equiv10\pmod{20}$ mm fall at
        $\pm2$ mm from a window centre --- hence the
        ``(window-track exception)'' label in the tables.
  \item Effect is local, symmetric, characterised and bounded:
        localization diagonal and channel map remain valid.
\end{itemize}
\column{0.46\textwidth}
\begin{tikzpicture}[scale=0.67,>=stealth,font=\scriptsize]
  \fill[cyan!12] (-3,0) rectangle (3,3.2);
  \node at (0,2.9) {EJ-204 bar};
  \draw[line width=2.5pt,gray!60] (-3,0) -- (-1.2,0);
  \fill[yellow!40,draw=gray!70,line width=0.7pt] (-1.2,-0.35) rectangle (0,0);
  \draw[line width=2.5pt,gray!60] (0,0) -- (1.0,0);
  \fill[yellow!40,draw=gray!70,line width=0.7pt] (1.0,-0.35) rectangle (2.2,0);
  \draw[line width=2.5pt,gray!60] (2.2,0) -- (3,0);
  \node[below] at (-0.6,-0.35) {W$_1$};
  \node[below] at (1.6,-0.35) {W$_2$};
  % Track A: centred
  \draw[ultra thick,blue!80] (-0.6,3.2) -- (-0.6,2.1)
        node[right,blue!80] {\textbf{A} (centred)};
  \draw[->,thick,blue!50,dashed] (-0.6,2.1) -- (-0.6,0.02);
  % Track B: +2 mm
  \draw[ultra thick,red!80] (0.0,3.2) -- (0.0,2.1)
        node[right,red!80] {\textbf{B} ($+2$ mm)};
  \draw[->,thick,red!50,dashed] (0.0,2.1) -- (0.0,-0.15);
  % Track C: midpoint
  \draw[ultra thick,green!55!black] (0.5,3.2) -- (0.5,2.1)
        node[right,green!55!black] {\textbf{C} (midpoint)};
  \draw[->,thick,green!45!black,dashed] (0.5,2.1) -- (0.5,1.8);
  \node[below,align=center] at (-3.5,-1.2) {};
\end{tikzpicture}
\end{columns}
\begin{alertblock}{Consequence for the scan}
Track B subtends a larger first-pass solid angle to W$_1$; fewer multiple
reflections $\Rightarrow$ higher $N_{pe}$ at the nearest channel.
At midpoint (C) both windows contribute equally and the effect is null.
\end{alertblock}
""",
    )


def fano_motivation_slide(c: float, c_err: float) -> str:
    """T5a: why the Fano factor grows linearly with mean N_pe."""
    sqrt_c = math.sqrt(c)
    return frame(
        r"Why the Fano factor grows with $\langle N_{pe}\rangle$",
        rf"""
\begin{{columns}}[T]\column{{0.56\textwidth}}
\begin{{enumerate}}\small
  \item \textbf{{Pure Poisson (fixed $\Delta E$):}}
        $\mathrm{{Var}}(N_{{pe}})=\langle N_{{pe}}\rangle$, hence $F\equiv1$.
  \item \textbf{{Landau fluctuations in $\Delta E$}} (law of total variance):
  \[
    \mathrm{{Var}}(N_{{pe}})=\langle N_{{pe}}\rangle+\!\Bigl(\tfrac{{\langle N_{{pe}}\rangle}}{{\langle\Delta E\rangle}}\Bigr)^{{\!2}}\!\mathrm{{Var}}(\Delta E),
  \]
  $\Rightarrow\;F=1+c\langle N_{{pe}}\rangle$,\enspace
  $c=\mathrm{{Var}}(\Delta E)/\langle\Delta E\rangle^2$
  (collection-efficiency terms absorbed in $c$).
  \item \textbf{{Poisson limit $F\!\to\!1$:}} only at $\langle N_{{pe}}\rangle\!\to\!0$
        (counting noise dominates over energy-loss fluctuations).
  \item \textbf{{Effective relative fluctuation:}} $\sqrt{{c}}\approx{sqrt_c:.3f}$
        is the combined energy-loss + collection jitter,
        consistent with a thin (10 mm) plastic absorber.
\end{{enumerate}}
\column{{0.42\textwidth}}
\begin{{block}}{{Fit result (see next slide)}}
$c = {c:.6f}\pm{c_err:.6f}$\par\medskip
The \emph{{same}} Landau $\Delta E$ fluctuation that shapes
the asymmetric $N_{{pe}}$ spectrum (next section)
also generates the super-Poissonian Fano factor:
\textbf{{one physical mechanism, two observables}}.
\end{{block}}
\end{{columns}}
""",
    )


def moyal_motivation_slide() -> str:
    """T5b: why Landau (not Gauss) and why Moyal as proxy."""
    return frame(
        r"Why Landau (and why Moyal as analytic proxy)",
        r"""
\begin{columns}[T]
\column{0.50\textwidth}
\textbf{Why not Gaussian:}\par\smallskip
\begin{itemize}\small
  \item In 10 mm of plastic the muon undergoes many collisions, but
        the energy-transfer cross-section ($\mathrm{d}\sigma/\mathrm{d}\varepsilon\propto1/\varepsilon^2$,
        Rutherford/M{\o}ller) has a hard tail.
  \item Rare large transfers ($\delta$-rays) dominate the variance;
        the variance per collision diverges and the Central Limit Theorem
        does not apply.
  \item Result: asymmetric distribution with a hard right tail ---
        the Landau distribution.
  \item Gaussian (Vavilov$\to$Gauss) limit only applies to \emph{thick} absorbers.
\end{itemize}
\column{0.48\textwidth}
\textbf{Why Moyal:}\par\smallskip
\begin{itemize}\small
  \item Landau density has no closed form (contour integral).
  \item The Moyal function is a closed-form analytic approximation
        to the Landau core and early tail, enabling stable MINUIT fits.
  \item \textit{Known limitation:} extreme Moyal tail undershoots
        true Landau tail; fit quality degrades at high $N_{pe}$.
\end{itemize}
\medskip
\textbf{Moyal $\otimes$ Gauss:}\par\smallskip
\begin{itemize}\small
  \item Photo-electron counting (Poisson) and collection fluctuations
        broaden the peak symmetrically.
  \item Full model: Moyal $\otimes$ Gaussian.
  \item Connects to Fano slide: same Landau $\Delta E$ fluctuation,
        two observables (spectrum shape + super-Poissonian $F$).
\end{itemize}
\end{columns}
""",
    )


def timing_gate_slide_exec12(timing_gate: pd.DataFrame, tail_all: pd.DataFrame) -> str:
    """T3 (Option B): rename slide, add no-time-window definition box."""
    body = rf"""
\begin{{columns}}[T]\column{{0.48\textwidth}}
\scriptsize\resizebox{{\textwidth}}{{!}}{{\begin{{tabular}}{{rrrrr}}\toprule
$x$ & side & EndTop [ps] & End-only [ps] & ratio \\ \midrule
{timing_table(timing_gate)}
\bottomrule\end{{tabular}}}}
\column{{0.50\textwidth}}
\includegraphics[width=0.88\textwidth]{{figs/exec09_tail_comparison.png}}
\end{{columns}}
\begin{{alertblock}}{{Confirmed mechanism}}
Across $x=0,\pm400$ mm, EndTop lowers $t99$ by
{abs(tail_all.t99_ns_difference.mean()):.2f} ns and suppresses late fractions
at more than {abs(tail_all.frac_gt20_difference_sigma).min():.1f} sigma.
\end{{alertblock}}
\begin{{block}}{{\footnotesize No time acceptance window}}
\scriptsize
In all published analyses (EXEC\_08b--EXEC\_09), End-photon arrival times
are reported without any $t<t_{{max}}$ acceptance cut: all detected End photons
in the ROOT TTree are used
(\texttt{{exec08b\_timing\_gate.py}}, \texttt{{exec09\_timing\_mechanism.py}}).
The EndTop tail narrowing persists over the full arrival-time range
and is therefore physical, not an analysis artefact.
\end{{block}}
"""
    return frame(
        "End arrival-time tails without any time acceptance window: narrower EndTop is physical",
        body,
    )


def top_slide_exec12(summary: pd.DataFrame, tN_summary: pd.DataFrame | None) -> str:
    """T5d: top timing estimate slide with derivation line and optional t_N comparison."""
    table_rows = top_estimate_table(summary)

    comparison_block = ""
    if tN_summary is not None:
        t4_top = tN_summary[
            (tN_summary["N"] == 4) & (tN_summary["group"] == "top_nearest")
        ]
        if not t4_top.empty:
            sim_rows = []
            for x_mm in KEY_POSITIONS:
                row = t4_top[t4_top["x_mm"] == x_mm]
                if row.empty:
                    continue
                est_row = summary[
                    (summary.x_beam_mm == x_mm) & (summary.grupo == "nearest_top")
                ]
                if est_row.empty:
                    continue
                npe = float(est_row.iloc[0].npe_mean)
                sigma_est_ps = float(est_row.iloc[0].sigma_est_ps)
                s4_ns = float(row.iloc[0].sigma_fit_ns)
                s4_ps = s4_ns * 1000.0
                s4_err_ns = float(row.iloc[0].sigma_fit_err_ns)
                s4_err_ps = s4_err_ns * 1000.0 if math.isfinite(s4_err_ns) else math.nan
                npe_sem = float(est_row.iloc[0].npe_rms) / math.sqrt(N_EVENTS)
                sigma_est_err = 0.5 * sigma_est_ps * npe_sem / npe
                s4_str = (
                    rf"${s4_ps:.0f}\pm{s4_err_ps:.0f}$"
                    if math.isfinite(s4_ps) and math.isfinite(s4_err_ps)
                    else r"—"
                )
                sim_rows.append(
                    rf"{x_mm:+d} & ${sigma_est_ps:.0f}\pm{sigma_est_err:.0f}$ & {s4_str} \\"
                )
            if sim_rows:
                comparison_block = rf"""
\vspace{{0.25cm}}
\begin{{block}}{{\footnotesize Analytic estimate vs.\ simulated $\sigma(t_4)$ (T4)}}
\scriptsize\begin{{tabular}}{{rcc}}\toprule
$x$ [mm] & $\sigma_{{t,est}}$ [ps] (analytic) & $\sigma(t_4)$ [ps] (simulated)\\ \midrule
{chr(10).join(sim_rows)}
\bottomrule\end{{tabular}}
\par\scriptsize $\sigma(t_4)$ is the photon-counting lower bound; $\sigma_{{t,est}}$
includes the SPE pulse-shape contribution ($\sqrt{{\tau_r\tau_d}}$) and is
therefore larger --- consistent with the SPE-convolution jitter seen in $\sigma_{{group}}$.
\end{{block}}"""

    return frame(
        "Top readout: analytic timing estimate only",
        rf"""
\centering\small
\begin{{tabular}}{{rrr}}\toprule
$x$ [mm] & nearest Top $N_{{pe}}$ & $\sigma_{{t,est}}$ [ps] \\ \midrule
{table_rows}
\bottomrule\end{{tabular}}
\vspace{{0.35cm}}
\begin{{block}}{{Definition, derivation, and scope}}
$\sigma_{{t,est}}=\sqrt{{\tau_r\cdot\tau_d}}/\sqrt{{N_{{pe}}}}
=1.122\,\mathrm{{ns}}/\sqrt{{N_{{pe}}}}$
\par\scriptsize
Standard single-photodetector timing scaling for a threshold at low $N_{{pe}}$ with
double-exponential emission ($\tau_r\cdot\tau_d$ product, rise $\times$ decay;
see e.g.\ Knoll, \emph{{Radiation Detection and Measurement}}, ch.\ 8).
$\tau_r=0.7$ ns, $\tau_d=1.8$ ns from the verified SSLG4 MPT.
This is an \textbf{{analytic orientation only}}, not a simulated timing resolution.
Top has no test-beam counterpart.
\end{{block}}{comparison_block}
""",
    )


def tN_position_slide(run_idx: int, x_mm: int) -> str:
    """T4: one slide showing the 6-panel t_N figure for a key position."""
    return image_frame(
        rf"Run {run_idx} $|$ x = {x_mm:+d} mm $|$ Time-to-threshold $t_N$",
        f"figs/exec12_tN_{x_mm}mm.png",
        height=0.84,
    )


def tN_synthesis_slide() -> str:
    """T4: sigma_fit(t_N) vs x synthesis slide."""
    return frame(
        r"$\sigma(t_N)$ vs beam position: synthesis",
        r"""
\centering\includegraphics[width=\textwidth,height=0.65\textheight,keepaspectratio]%
{figs/exec12_tN_summary.png}
\begin{block}{\footnotesize Interpretation and connection to $\sigma_{\rm group}$}
\scriptsize
$\sigma(t_N)$ is the single-side timing resolution at an ideal digital threshold of
$N$ detected PE; it connects the cumulative $\langle N(t)\rangle$ curves to the
published $\sigma_{\rm group}$ values.
At high $N_{pe}$ (x=$\pm690,\pm650,\pm400$ near their respective End),
$\sigma(t_4)\sim15$ ps is the pure photon-counting limit;
$\sigma_{\rm group}\sim50$ ps is larger because the analog SUM4 discriminator
adds jitter from the SPE pulse convolution ($\tau_r=0.7$ ns).
At $x=0$ (moderate $N_{pe}$, both Ends equidistant), photon-counting noise dominates
and both converge to $\sim120$ ps.
\end{block}
""",
    )


def glossary_exec12() -> str:
    """T6: extended glossary including exec12 terms."""
    return frame(
        "Backup: glossary",
        r"""
\scriptsize\begin{description}
 \item[HOOK\_WALK] Reserved insertion point for future leading-edge time-walk correction.
 Fixed 4-PE threshold means larger pulses cross earlier; no correction applied in any campaign
 (\texttt{congruent\_sum4\_timing.C:305--315}).
 \item[FPT] First photon time: first detected photon per event for the stated channel/group.
 \item[SUM4] Analog sum of four SiPMs: \{0--3\}, \{4--7\}, \{8--11\}, \{12--15\}
 (\texttt{congruent\_sum4\_timing.C:218--221}).
 \item[Window-track collection effect] Local ($\pm2$ mm) increase of first-pass solid angle when a
 track passes offset from a Top window centre. Empirically $\approx9.8$\% at $\approx12\,\sigma$
 (five-run test); null at centre and midpoint; symmetric. Geometric working hypothesis.
 \item[Timing gate] The comparison script \texttt{exec08b\_timing\_gate.py} tests whether
 EndTop $\sigma_{\rm group}$ is smaller than End-only $\sigma_{\rm group}$ at $x=0,\pm400$ mm.
 It applies no time acceptance window on individual photon arrivals; all ROOT TTree photons are used.
 \item[$t_N$ / time-to-threshold] Arrival time of the $N$-th detected photoelectron per event
 in a given group. $\sigma(t_N)$ is the single-side photon-counting timing resolution at
 an ideal $N$-PE digital threshold. Excludes SPE pulse convolution jitter present in $\sigma_{\rm group}$.
 \item[Fano factor] $F\equiv\mathrm{Var}(N_{pe})/\langle N_{pe}\rangle$.
 Equals 1 for pure Poisson; grows linearly with $\langle N_{pe}\rangle$ when Landau
 $\Delta E$ fluctuations are present: $F=1+c\langle N_{pe}\rangle$.
 \item[Moyal] Closed-form analytic approximation to the Landau density, used for
 stable MINUIT fits. Undershoots the extreme Landau right tail.
 \item[(exact geometry)] Nearest SiPM agrees with maximum, or covered by window-track exception.
 x=300 mm is a statistical tie.
 \item[Intrinsic sigma] Excludes SPTR and electronics jitter; 88 ps hardware value not directly comparable.
\end{description}
""",
    )


def make_tex(output_dir: pathlib.Path, mode: str, with_arrival: bool = False, with_pedagogy: bool = False) -> pathlib.Path:
    if with_pedagogy:
        with_arrival = True  # exec12 always includes arrival slides
    tN_summary: pd.DataFrame | None = None
    if with_pedagogy:
        tN_csv = output_dir / "exec12_tN_summary.csv"
        if tN_csv.is_file():
            tN_summary = pd.read_csv(tN_csv)
    positions = pd.read_csv(output_dir / "per_position_exec07.csv").sort_values("x_beam_mm")
    positions_indexed = positions.set_index("x_beam_mm")
    localization = pd.read_csv(output_dir / "top_localization_gate.csv").set_index("x_beam_mm")
    summary = pd.read_csv(output_dir / "summary_exec07.csv")
    fits = pd.read_csv(output_dir / "fit_results_exec07.csv").set_index("side")
    profiles = pd.read_csv(output_dir / "exec08b_window_dip_profiles.csv")
    timing_gate = pd.read_csv(output_dir / "exec08b_timing_gate.csv")
    tails = pd.read_csv(output_dir / "exec09_tail_comparison.csv")
    fano_fit = pd.read_csv(output_dir / "exec10_fano_fit.csv").iloc[0]
    landau = pd.read_csv(output_dir / "exec10_landau_mpv.csv")
    velocity_fits = pd.read_csv(output_dir / "exec10_velocity_fits.csv")
    selected = positions if mode == "full" else positions[positions.x_beam_mm.isin(KEY_POSITIONS)]
    nearest = summary[summary.grupo == "nearest_top"]
    end_clusters = summary[summary.grupo.isin(["end_left_A_SUM4", "end_right_A_SUM4"])]
    tail_all = tails[tails.group == "all_end"]
    center_top = stats(summary, 0, "nearest_top")
    edge_top = stats(summary, -690, "nearest_top")
    center_ends = stats(summary, 0, "ends_all")
    edge_ends = stats(summary, -690, "ends_all")

    _title_line = (
        r"\title{EXEC\_12: EndTop photon budget --- pedagogical review}"
        if with_pedagogy else
        r"\title{EXEC\_10: EndTop photon budget and intrinsic timing}"
    )
    _subtitle_line = (
        r"\subtitle{t_N distributions, Fano/Landau motivation, window-track mechanism, full scan}"
        if with_pedagogy else
        r"\subtitle{Landau-dominated Npe, clarified timing densities, and apparent-slope diagnosis}"
    )
    preamble = (
        r"""\documentclass[aspectratio=169]{beamer}
\IfFileExists{beamerthememetropolis.sty}{\usetheme{metropolis}}{\usetheme{default}}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{xcolor}
\usepackage{tikz}
\definecolor{shipblue}{RGB}{20,74,130}
\definecolor{warningorange}{RGB}{210,105,30}
\setbeamercolor{structure}{fg=shipblue}
\setbeamercolor{alerted text}{fg=warningorange}
\setbeamertemplate{navigation symbols}{}
"""
        + _title_line + "\n"
        + _subtitle_line + "\n"
        + r"""\author{Rene Rios (ULS)}
\date{June 11, 2026}
\begin{document}
"""
    )
    title = frame(
        "Configuration and provenance",
        r"""
\begin{columns}[T]\column{0.62\textwidth}
\textbf{Branch:} \texttt{feat/endtop-sslg4}\\
\textbf{Commits:} \texttt{0201413}, \texttt{7b4bc80}, \texttt{ad26d78},
\texttt{bcefa24}, \texttt{bdfe955}; this deck is generated from those published sources
\begin{itemize}\small
  \item 31 positions $\times$ 2000 events; 16 End + 70 Top channels; jitter fixed to zero.
  \item EJ-204 via SSLG4 \texttt{OPSC-101}: yield 10400/MeV, $\tau_r=0.7$ ns,
        $\tau_d=1.8$ ns, ABSLENGTH=160 cm (verified effective MPT inputs).
  \item Broadcom AFBR-S4N66P024M: PDE 63\% at 420 nm.
\end{itemize}
\column{0.35\textwidth}
\begin{alertblock}{Scope}
Top 70 SiPM is a simulation coverage study. Real hardware: 32 Top SiPMs.
All reported timing is intrinsic, without SPTR or electronics jitter.
\end{alertblock}
\end{columns}
""",
    )
    geometry = frame(
        "EndTop geometry and channel map",
        r"""
\begin{columns}[T]\column{0.55\textwidth}
\begin{itemize}
 \item EJ-204 bar: $1400\times60\times10$ mm$^3$.
 \item End L IDs 0--7; End R IDs 8--15.
 \item Top L IDs 16--50: $-692,-672,\ldots,-12$ mm.
 \item Top R IDs 51--85: $+12,+32,\ldots,+692$ mm.
 \item Each Top element and window: $6\times6$ mm$^2$.
\end{itemize}
\column{0.42\textwidth}
\begin{block}{Exact constraints}
First centers are 5 mm from each bar end. Pitch is 20 mm except the deliberate
central $-12/+12$ mm pair, whose pitch is 24 mm.
\end{block}
\end{columns}
""",
    )
    faces = frame(
        "Instrumented and wrapped faces",
        r"""
\centering\begin{tabular}{lll}\toprule
Face & Optical state & Readout \\ \midrule
$-X$, $+X$ & open & 8 End SiPMs each \\
$+Y$ & Mylar with 70 exact windows & 70 Top SiPMs \\
$-Y$, $-Z$, $+Z$ & complete Mylar & none \\
\bottomrule\end{tabular}
\vspace{0.5cm}
\begin{block}{External layer}
Black Tedlar is deliberately not modeled; it only blocks ambient light.
The non-reflected Mylar fraction is already absorbed in the optical model.
\end{block}
""",
    )
    optical = frame(
        "Window optics: why timing and collection change",
        r"""
\centering
\begin{tikzpicture}[scale=0.78,>=stealth]
 \draw[fill=cyan!12] (0,0) rectangle (5,4);
 \node at (2.5,3.6) {EJ-204 bar};
 \draw[line width=3pt,gray] (0,0) -- (5,0);
 \node[below] at (2.5,0) {Mylar, fixed $R=0.98$: reflected photons return};
 \draw[->,very thick,blue] (1.0,2.5) -- (2.1,0.15);
 \draw[->,very thick,blue] (2.1,0.15) -- (3.2,2.3);
 \draw[fill=cyan!12] (8,0) rectangle (13,4);
 \node at (10.5,3.6) {EJ-204 bar};
 \draw[line width=3pt,gray] (8,0) -- (9.6,0);
 \draw[line width=3pt,gray] (11.4,0) -- (13,0);
 \draw[fill=yellow!25] (9.6,-0.45) rectangle (11.4,0);
 \node[below] at (10.5,-0.45) {index-matched window};
 \draw[->,very thick,red] (9.4,2.5) -- (10.5,-0.3);
 \node[align=center] at (10.5,1.0) {PDE: detected\\otherwise absorbed\\never returns};
\end{tikzpicture}
\begin{alertblock}{Consequence}
Every window is an optical drain. Long, multiply reflected paths have more chances
to meet one, so late End photons are preferentially removed.
\end{alertblock}
""",
    )
    if with_arrival:
        _parts = []
        for _i, (_, _row) in enumerate(selected.iterrows()):
            _parts.append(position_frame(_row, summary, localization))
            _x = int(_row.x_beam_mm)
            _parts.append(image_frame(
                rf"Run {_i + 1} $|$ x = {_x:+d} mm $|$ Photon arrival $\langle N(t)\rangle$",
                f"figs/exec11_arrival_{_x}mm.png",
            ))
        position_slides = "".join(_parts)
    else:
        position_slides = "".join(position_frame(row, summary, localization) for _, row in selected.iterrows())
    integrated = "".join([
        image_frame("Photon budget versus beam position", "figs/P1_npe_vs_x.png"),
        image_frame("Top localization: full coverage diagonal", "figs/P2_npe_heatmap_top.png"),
        image_frame(
            "Attenuation and the historical mean-time slope", "figs/fits_attenuation_velocity.png",
            rf"$\lambda_{{eff,L}}={fits.loc['left','lambda_eff_cm']:.2f}\pm{fits.loc['left','lambda_eff_error_cm']:.2f}$ cm, "
            rf"$\lambda_{{eff,R}}={fits.loc['right','lambda_eff_cm']:.2f}\pm{fits.loc['right','lambda_eff_error_cm']:.2f}$ cm; "
            r"the mean-time slope is tail-biased and is not a propagation velocity.",
        ),
        image_frame(
            "Landau dominance: Fano factor versus mean Npe", "figs/exec10_fano_vs_mean.png",
            rf"$F=1+c\langle N_{{pe}}\rangle$, $c={fano_fit.c:.6f}\pm{fano_fit.c_error:.6f}$, "
            rf"$\chi^2/ndf={fano_fit.chi2_ndf:.2f}$. Poisson $F=1$ is approached only at low light.",
        ),
        image_frame(
            "High-light example: Moyal-Gaussian Landau proxy", "figs/exec10_landau_fit_example.png",
            r"Moyal approximates Landau. Some channels have poor chi2 or an unidentifiable Gaussian width; all are recorded in \texttt{exec10\_landau\_mpv.csv}.",
        ),
        image_frame(
            "Arrival-time examples: absolute rate and normalized shape", "figs/P5_tdist_examples.png",
            "Left: entries per event per ns. Right: the same distributions normalized to area 1. FPT is the first detected photon time in each event.",
        ),
        image_frame("Mean arrival time and FPT versus x", "figs/P6_tmean_vs_x.png"),
    ])
    landau_slide = frame(
        "Nearest-Top Landau/Moyal MPVs at key positions",
        rf"""
\centering\small
\begin{{tabular}}{{rrrrr}}\toprule
$x$ [mm] & channel & MPV [Npe] & Landau width [Npe] & $\chi^2/ndf$ \\ \midrule
{landau_table(positions_indexed, landau)}
\bottomrule\end{{tabular}}
\begin{{block}}{{Interpretation}}
The right tail originates primarily from Landau-like energy-loss fluctuations; Poisson
counting dominates only in the few-PE limit. MPVs are future inputs for the pending
test-beam ToT(NPE) calibration comparison.
\end{{block}}
""",
    )
    dip = frame(
        "Window-track collection effect: five-run test",
        rf"""
\begin{{columns}}[T]\column{{0.47\textwidth}}
\scriptsize\resizebox{{\textwidth}}{{!}}{{\begin{{tabular}}{{lrr}}\toprule
Run & contrast [Npe] & significance \\ \midrule
{dip_table(profiles)}
\bottomrule\end{{tabular}}}}
\column{{0.51\textwidth}}
\includegraphics[width=\textwidth]{{figs/exec08b_window_dip_profiles.png}}
\end{{columns}}
\begin{{block}}{{Finding}}
Local symmetric $\pm2$ mm effect: about 9.8\% at about 12 sigma; null when centered
and at midpoint. The scan and channel map remain valid.
\end{{block}}
""",
    )
    timing = image_frame(
        "Intrinsic End SUM4 timing", "figs/P7_deltaT_end.png",
        r"$\sigma_{group}=\sigma(\Delta T_{AB})/\sqrt{2}$; no SPTR/jitter. "
        r"The 88 ps hardware context is not directly comparable.",
    )
    timing_gate_slide = frame(
        "Timing gate lifted: narrower EndTop is physical",
        rf"""
\begin{{columns}}[T]\column{{0.48\textwidth}}
\scriptsize\resizebox{{\textwidth}}{{!}}{{\begin{{tabular}}{{rrrrr}}\toprule
$x$ & side & EndTop [ps] & End-only [ps] & ratio \\ \midrule
{timing_table(timing_gate)}
\bottomrule\end{{tabular}}}}
\column{{0.50\textwidth}}
\includegraphics[width=0.88\textwidth]{{figs/exec09_tail_comparison.png}}
\end{{columns}}
\begin{{alertblock}}{{Confirmed mechanism}}
Across $x=0,\pm400$ mm, EndTop lowers $t99$ by
{abs(tail_all.t99_ns_difference.mean()):.2f} ns and suppresses late fractions
at more than {abs(tail_all.frac_gt20_difference_sigma).min():.1f} sigma.
\end{{alertblock}}
\tiny Left panels show absolute End photons/event/ns. Right panels normalize
each distribution to area 1 and isolate shape: the $t>10$ ns tail is suppressed in EndTop
because multiply reflected photons have more opportunities to meet a Top window.
""",
    )
    velocity_slide = frame(
        "Apparent timing slopes: no propagation estimator passes",
        rf"""
\begin{{columns}}[T]\column{{0.62\textwidth}}
\includegraphics[width=\textwidth]{{figs/exec10_velocity_estimators.png}}
\column{{0.36\textwidth}}
\scriptsize\resizebox{{\textwidth}}{{!}}{{\begin{{tabular}}{{lrr}}\toprule
Estimator & slope-derived [cm/ns] & $\chi^2/ndf$ \\ \midrule
{velocity_table(velocity_fits)}
\bottomrule\end{{tabular}}}}
\begin{{alertblock}}{{Refuted predictions}}
Mean FPT is not close to $c/n\simeq19$ cm/ns. The fixed-threshold $f(t>10\,ns)$
increases with distance. Every global line fit is rejected.
\end{{alertblock}}
\end{{columns}}
""",
    )
    top_slide = frame(
        "Top readout: analytic timing estimate only",
        rf"""
\centering\small
\begin{{tabular}}{{rrr}}\toprule
$x$ [mm] & nearest Top $N_{{pe}}$ & $\sigma_{{t,est}}$ [ps] \\ \midrule
{top_estimate_table(summary)}
\bottomrule\end{{tabular}}
\vspace{{0.35cm}}
\begin{{block}}{{Definition and scope}}
$\sigma_{{t,est}}=\sqrt{{0.7\,\mathrm{{ns}}\times1.8\,\mathrm{{ns}}}}/\sqrt{{N_{{pe}}}}
=1.122\,\mathrm{{ns}}/\sqrt{{N_{{pe}}}}$. This is an analytic orientation,
not a simulated timing resolution; Top has no test-beam counterpart.
\end{{block}}
""",
    )
    conclusions = frame(
        "Numerical conclusions",
        rf"""
\small\begin{{itemize}}
 \item At center: End sum {npe_text(center_ends)}, nearest Top {npe_text(center_top)};
       at $x=-690$ mm: End sum {npe_text(edge_ends)}, nearest Top {npe_text(edge_top)}.
 \item $\lambda_{{eff}}$: L ${fits.loc['left','lambda_eff_cm']:.2f}\pm{fits.loc['left','lambda_eff_error_cm']:.2f}$ cm,
       R ${fits.loc['right','lambda_eff_cm']:.2f}\pm{fits.loc['right','lambda_eff_error_cm']:.2f}$ cm;
       no global time estimator yields an acceptable propagation-velocity fit.
 \item Nearest-Top var/mean: median {nearest.var_over_mean.median():.2f}, range
       {nearest.var_over_mean.min():.2f}--{nearest.var_over_mean.max():.2f};
       $F=1+({fano_fit.c:.6f}\pm{fano_fit.c_error:.6f})\langle N_{{pe}}\rangle$ with
       $\chi^2/ndf={fano_fit.chi2_ndf:.2f}$: Landau-like energy-loss fluctuations dominate.
 \item Intrinsic End $\sigma_{{group}}$: mean {end_clusters.sigma_grupo_ps.mean():.1f} ps,
       range {end_clusters.sigma_grupo_ps.min():.1f}--{end_clusters.sigma_grupo_ps.max():.1f} ps.
 \item EndTop/End-only timing ratio: ${timing_gate.observed_sigma_ratio.min():.3f}\pm{timing_gate.observed_sigma_ratio_error.max():.3f}$
       to ${timing_gate.observed_sigma_ratio.max():.3f}\pm{timing_gate.observed_sigma_ratio_error.max():.3f}$.
 \item Three findings: window-track dip characterized; late-tail drain confirmed;
       complete 86-channel coverage and localization gate validated.
\end{{itemize}}
""",
    )
    backup = "".join([
        frame(
            "Backup: validation and gates",
            r"""
\begin{itemize}
 \item 31/31 ROOT files pass: event IDs 0--1999, matching beam x, aggregate IDs 0--85.
 \item Top localization passes all positions: strict nearest maximum except the documented
       $x\equiv10\pmod{20}$ window-track effect (maximum among nearest two, deficit below 15\%).
 \item At x=300 mm the two candidates differ by only 0.08 sigma and pass as a statistical tie.
 \item MPT checks pass: RINDEX 1.58, yield 10400/MeV, decay 1.8 ns, rise 0.7 ns,
       ABSLENGTH 160 cm, emission peak 408.8 nm.
\end{itemize}
""",
        ),
        image_frame("Backup: ID18 impact maps", "figs/exec08b_id18_impact_maps.png"),
        frame(
            "Backup: timing methodology",
            r"""
\begin{itemize}
 \item SUM4 groups: \{0..3\}, \{4..7\}, \{8..11\}, \{12..15\}.
 \item Single-PE response: normalized 0.5/5 ns double exponential.
 \item Absolute leading-edge threshold: 4 PE equivalent.
 \item Estimator: $\sigma(\Delta T_{AB})/\sqrt{2}$.
 \item \texttt{HOOK\_WALK}: pending; no time-walk correction exists or was invented.
 \item EXEC\_09 anti-artifact audit verifies common threshold, pulse, t=0 and estimator.
\end{itemize}
""",
        ),
        frame(
            "Backup: glossary",
            r"""
\scriptsize\begin{description}
 \item[HOOK\_WALK] Reserved insertion point for future leading-edge time-walk correction.
 Fixed 4-PE threshold means larger pulses cross earlier; no correction is applied
 in any campaign (\texttt{congruent\_sum4\_timing.C:305--315}).
 \item[FPT] First photon time: first detected photon in an event for the stated channel/group.
 \item[SUM4] Analog sum of four SiPMs: \{0--3\}, \{4--7\}, \{8--11\}, \{12--15\}.
 \item[(exact geometry)] Nearest geometric SiPM agrees with the maximum, or is covered by
 the documented window-track exception; x=300 mm is explicitly a statistical tie.
 \item[Temporal density] Absolute plots are entries/event/ns; normalized plots integrate to one
 and compare shape only.
 \item[Central 24 mm pair] The sole pitch exception is $-12/+12$ mm because 1384 mm is not
 divisible by 20 mm.
 \item[Intrinsic sigma] Excludes SPTR and electronics jitter; the hardware 88 ps value is not
 directly comparable.
\end{description}
""",
        ),
        image_frame("Backup: EndTop versus End-only tail audit", "figs/exec09_tail_comparison.png"),
    ])
    if with_pedagogy:
        # T5c: statistics section with readingboxes and T5a/T5b motivation slides
        _lambda_l = fits.loc["left", "lambda_eff_cm"]
        _lambda_r = fits.loc["right", "lambda_eff_cm"]
        _lambda_l_e = fits.loc["left", "lambda_eff_error_cm"]
        _lambda_r_e = fits.loc["right", "lambda_eff_error_cm"]
        integrated_exec12 = "".join([
            image_frame_wb(
                "Photon budget versus beam position", "figs/P1_npe_vs_x.png",
                why=(
                    "Quantify $N_{pe}$ collected at each End and nearest-Top channel"
                    " as the muon track scans the bar."
                ),
                axes=(
                    r"x: beam position [mm]; y: mean $N_{pe}$ per event"
                    r" (End sums left/right and nearest Top)."
                ),
                takeaway=(
                    rf"$N_{{pe}}$ decreases exponentially with track-to-End distance"
                    rf" ($\lambda_{{eff}}\approx{0.5*(_lambda_l+_lambda_r):.0f}$ cm"
                    rf" $\ll$ ABSLENGTH $=160$ cm, system effective length);"
                    r" nearest-Top $N_{pe}$ peaks directly under the SiPM window."
                ),
                height=0.62,
            ),
            image_frame_wb(
                "Top localization: full coverage diagonal", "figs/P2_npe_heatmap_top.png",
                why=(
                    "Verify that the highest-$N_{pe}$ Top channel is always"
                    " the geometrically nearest one across all 31 positions."
                ),
                axes=(
                    r"x: beam position [mm]; y: Top channel ID; colour: mean $N_{pe}$."
                ),
                takeaway=(
                    r"Strict diagonal coverage; window-track exception positions"
                    r" ($x\equiv10\pmod{20}$ mm) are labelled and pass the"
                    r" $<15$\% deficit criterion."
                ),
                height=0.62,
            ),
            image_frame_wb(
                "Attenuation and the historical mean-time slope", "figs/fits_attenuation_velocity.png",
                note=(
                    rf"$\lambda_{{eff,L}}={_lambda_l:.2f}\pm{_lambda_l_e:.2f}$ cm, "
                    rf"$\lambda_{{eff,R}}={_lambda_r:.2f}\pm{_lambda_r_e:.2f}$ cm; "
                    r"the mean-time slope is tail-biased and is not a propagation velocity."
                ),
                why=(
                    r"Measure the effective attenuation length $\lambda_{eff}$"
                    r" and test whether the mean-time slope encodes a light-propagation velocity."
                ),
                axes=(
                    r"Left: mean $N_{pe}$ vs track-to-End distance [cm];"
                    r" exponential fit gives $\lambda_{eff}$."
                    r" Right: mean arrival time [ns] vs same distance."
                ),
                takeaway=(
                    rf"$\lambda_{{eff}}\approx{0.5*(_lambda_l+_lambda_r):.0f}$ cm"
                    r" $\ll$ ABSLENGTH $=160$ cm:"
                    r" $\lambda_{eff}$ is a system effective length combining"
                    r" bulk absorption, Mylar reflection losses, and Top window drainage---not"
                    r" the material property. The mean-time slope is tail-biased; no propagation"
                    r" velocity is reliably extracted."
                ),
                height=0.55,
            ),
            fano_motivation_slide(float(fano_fit.c), float(fano_fit.c_error)),
            image_frame_wb(
                "Landau dominance: Fano factor versus mean $N_{pe}$",
                "figs/exec10_fano_vs_mean.png",
                note=(
                    rf"$F=1+c\langle N_{{pe}}\rangle$, "
                    rf"$c={fano_fit.c:.6f}\pm{fano_fit.c_error:.6f}$, "
                    rf"$\chi^2/ndf={fano_fit.chi2_ndf:.2f}$."
                    r" Poisson $F=1$ is approached only at low light."
                ),
                why=(
                    r"Measure $F=\mathrm{Var}(N_{pe})/\langle N_{pe}\rangle$ to distinguish"
                    r" Poisson ($F=1$) from Landau-dominated fluctuations."
                ),
                axes=(
                    r"x: $\langle N_{pe}\rangle$ per event; y: $\mathrm{Var}/\langle N_{pe}\rangle$."
                    r" Dashed: Poisson $F=1$. Line: linear fit $F=1+c\langle N_{pe}\rangle$."
                ),
                takeaway=(
                    r"$F$ grows linearly: Landau $\Delta E$ fluctuations dominate"
                    r" counting noise at all light levels. $\sqrt{c}$ is the effective"
                    r" relative energy-loss fluctuation of this 10 mm absorber."
                ),
                height=0.55,
            ),
            moyal_motivation_slide(),
            image_frame_wb(
                "High-light example: Moyal--Gaussian Landau proxy",
                "figs/exec10_landau_fit_example.png",
                note=(
                    r"Moyal approximates Landau. Some channels have poor $\chi^2$ or"
                    r" an unidentifiable Gaussian width; all recorded in"
                    r" \texttt{exec10\_landau\_mpv.csv}."
                ),
                why=(
                    r"Fit the $N_{pe}$ spectrum at a high-yield position to"
                    r" quantify the Landau MPV and width."
                ),
                axes=(
                    r"x: $N_{pe}$ per event; y: counts."
                    r" Data (step), Moyal$\otimes$Gauss fit (curve)."
                ),
                takeaway=(
                    r"Core described by Moyal proxy; hard right tail from $\delta$-rays."
                    r" MPVs are future inputs for test-beam ToT calibration."
                ),
                height=0.55,
            ),
            image_frame_wb(
                "Arrival-time examples: absolute rate and normalized shape",
                "figs/P5_tdist_examples.png",
                note=(
                    r"Left: entries per event per ns."
                    r" Right: area-normalized distributions. FPT = first detected photon time per event."
                ),
                why=(
                    r"Illustrate the arrival-time structure of individual Top channels"
                    r" and contrast absolute rate vs shape."
                ),
                axes=(
                    r"x: arrival time [ns]; y: left is entries/event/ns, right is area-normalized density [1/ns]."
                ),
                takeaway=(
                    r"FPT distribution captures the prompt scintillation leading edge;"
                    r" the broad tail reflects multiple-reflection photons."
                ),
                height=0.55,
            ),
            image_frame_wb(
                "Mean arrival time and FPT versus x", "figs/P6_tmean_vs_x.png",
                why=(
                    r"Track how the mean photon arrival time and FPT scale with beam position."
                ),
                axes=(
                    r"x: beam position [mm]; y: mean arrival time [ns] and mean FPT [ns]."
                ),
                takeaway=(
                    r"Both estimators increase with track-to-End distance; the mean is"
                    r" more tail-biased than FPT and yields a larger apparent slope."
                ),
                height=0.68,
            ),
        ])

        # T4: key position slides = table + arrival + t_N
        key_positions_df = positions[positions.x_beam_mm.isin(KEY_POSITIONS)]
        key_pos_parts = []
        for _i, (_, _row) in enumerate(key_positions_df.iterrows()):
            _x = int(_row.x_beam_mm)
            key_pos_parts.append(position_frame(_row, summary, localization))
            key_pos_parts.append(image_frame(
                rf"Run {_i + 1} $|$ x = {_x:+d} mm $|$ Photon arrival $\langle N(t)\rangle$",
                f"figs/exec11_arrival_{_x}mm.png",
            ))
            key_pos_parts.append(tN_position_slide(_i + 1, _x))
        key_position_slides = "".join(key_pos_parts)

        # T1: full scan section = table + arrival for all 31 positions
        full_pos_parts = []
        for _j, (_, _row) in enumerate(positions.iterrows()):
            _x = int(_row.x_beam_mm)
            full_pos_parts.append(position_frame(_row, summary, localization))
            full_pos_parts.append(image_frame(
                rf"Run {_j + 1} $|$ x = {_x:+d} mm $|$ Photon arrival $\langle N(t)\rangle$",
                f"figs/exec11_arrival_{_x}mm.png",
            ))
        full_position_slides = "".join(full_pos_parts)

        # T5c: readingbox on timing slides
        timing_wb = frame(
            "Intrinsic End SUM4 timing",
            (
                rf"\centering\includegraphics[width=\textwidth,height=0.60\textheight,"
                rf"keepaspectratio]{{figs/P7_deltaT_end.png}}"
                rf"\par\scriptsize $\sigma_{{group}}=\sigma(\Delta T_{{AB}})/\sqrt{{2}}$;"
                r" no SPTR/jitter. The 88 ps hardware context is not directly comparable."
                "\n"
                + readingbox(
                    r"Measure the intrinsic Same-End SUM4 coincidence timing resolution"
                    r" without SPTR or electronics jitter.",
                    r"x: beam position [mm]; y: $\sigma_{group}$ [ps]."
                    r" Lines for left and right End groups.",
                    r"$\sigma_{group}$ decreases toward the nearer End (higher $N_{pe}$,"
                    r" earlier FPT); consistent with photon-counting scaling $\propto1/\sqrt{N_{pe}}$.",
                )
            ),
        )
        velocity_wb = frame(
            "Apparent timing slopes: no propagation estimator passes",
            rf"""
\begin{{columns}}[T]\column{{0.62\textwidth}}
\includegraphics[width=\textwidth]{{figs/exec10_velocity_estimators.png}}
\column{{0.36\textwidth}}
\scriptsize\resizebox{{\textwidth}}{{!}}{{\begin{{tabular}}{{lrr}}\toprule
Estimator & slope-derived [cm/ns] & $\chi^2/ndf$ \\ \midrule
{velocity_table(velocity_fits)}
\bottomrule\end{{tabular}}}}
\begin{{alertblock}}{{Refuted predictions}}
Mean FPT is not close to $c/n\simeq19$ cm/ns. The fixed-threshold $f(t>10\,ns)$
increases with distance. Every global line fit is rejected.
\end{{alertblock}}
{readingbox(
    r"Test whether any timing estimator reveals a light-propagation velocity.",
    r"x: track-to-End distance [cm]; y: timing estimator [ns]. Fits shown.",
    r"No estimator yields a statistically acceptable propagation-velocity fit;"
    r" the slope is dominated by the long late-tail bias.",
)}
\end{{columns}}
""",
        )

        # T2 + T3 + T5d: window/timing section for exec12
        window_section_exec12 = (
            landau_slide
            + window_track_intro_slide()
            + dip
            + timing_wb
            + timing_gate_slide_exec12(timing_gate, tail_all)
            + velocity_wb
            + top_slide_exec12(summary, tN_summary)
        )

        # T6: conclusions with t_N synthesis
        _tN4_range = ""
        if tN_summary is not None:
            t4s = tN_summary[(tN_summary["N"] == 4) & tN_summary["sigma_fit_ns"].notna()]
            if not t4s.empty:
                t4_near = t4s[t4s["group"] == "end_near"]["sigma_fit_ns"].dropna() * 1000.0
                t4_top = t4s[t4s["group"] == "top_nearest"]["sigma_fit_ns"].dropna() * 1000.0
                if len(t4_near) > 0:
                    _tN4_range = (
                        rf" \item $\sigma(t_4)$ (near-End SUM4):"
                        rf" {t4_near.min():.0f}--{t4_near.max():.0f} ps (high-$N_{{pe}}$ to low);"
                        rf" pure photon-counting lower bound, factor $\sim3$ below"
                        rf" $\sigma_{{group}}$ at high $N_{{pe}}$ due to SPE pulse-convolution"
                        rf" jitter in the analog discriminator."
                    )
                if len(t4_top) > 0:
                    _tN4_range += (
                        rf" Nearest-Top $\sigma(t_4)$:"
                        rf" {t4_top.min():.0f}--{t4_top.max():.0f} ps."
                    )

        conclusions_exec12 = frame(
            "Numerical conclusions",
            rf"""
\small\begin{{itemize}}
 \item At center: End sum {npe_text(center_ends)}, nearest Top {npe_text(center_top)};
       at $x=-690$ mm: End sum {npe_text(edge_ends)}, nearest Top {npe_text(edge_top)}.
 \item $\lambda_{{eff}}$: L ${_lambda_l:.2f}\pm{_lambda_l_e:.2f}$ cm,
       R ${_lambda_r:.2f}\pm{_lambda_r_e:.2f}$ cm
       ($\ll$ ABSLENGTH $=160$ cm; system effective length, not material property).
 \item Nearest-Top $F=1+({fano_fit.c:.6f}\pm{fano_fit.c_error:.6f})\langle N_{{pe}}\rangle$
       ($\chi^2/ndf={fano_fit.chi2_ndf:.2f}$): Landau energy-loss fluctuations dominate.
       Same mechanism: asymmetric $N_{{pe}}$ spectrum (Moyal$\otimes$Gauss) and $F>1$.
 \item Intrinsic End $\sigma_{{group}}$: mean {end_clusters.sigma_grupo_ps.mean():.1f} ps,
       range {end_clusters.sigma_grupo_ps.min():.1f}--{end_clusters.sigma_grupo_ps.max():.1f} ps.
 \item EndTop/End-only timing ratio: ${timing_gate.observed_sigma_ratio.min():.3f}$
       to ${timing_gate.observed_sigma_ratio.max():.3f}$; no time-acceptance window in any script.
{_tN4_range}
 \item Four findings: window-track mechanism defined and characterised; late-tail drain confirmed;
       complete 86-channel coverage and localisation validated; $t_N$ photon-counting timing quantified.
\end{{itemize}}
""",
        )

        # T6: backup with extended glossary
        backup_exec12 = "".join([
            frame(
                "Backup: validation and gates",
                r"""
\begin{itemize}
 \item 31/31 ROOT files pass: event IDs 0--1999, matching beam x, aggregate IDs 0--85.
 \item Top localization passes all positions: strict nearest maximum except the documented
       $x\equiv10\pmod{20}$ window-track effect (maximum among nearest two, deficit below 15\%).
 \item At x=300 mm the two candidates differ by only 0.08 sigma and pass as a statistical tie.
 \item MPT checks pass: RINDEX 1.58, yield 10400/MeV, decay 1.8 ns, rise 0.7 ns,
       ABSLENGTH 160 cm, emission peak 408.8 nm.
\end{itemize}
""",
            ),
            image_frame("Backup: ID18 impact maps", "figs/exec08b_id18_impact_maps.png"),
            frame(
                "Backup: timing methodology",
                r"""
\begin{itemize}
 \item SUM4 groups: \{0..3\}, \{4..7\}, \{8..11\}, \{12..15\}.
 \item Single-PE response: normalized 0.5/5 ns double exponential.
 \item Absolute leading-edge threshold: 4 PE equivalent.
 \item Estimator: $\sigma(\Delta T_{AB})/\sqrt{2}$.
 \item \texttt{HOOK\_WALK}: pending; no time-walk correction exists or was invented.
 \item EXEC\_09 anti-artifact audit verifies common threshold, pulse, t=0 and estimator.
\end{itemize}
""",
            ),
            glossary_exec12(),
            image_frame("Backup: EndTop versus End-only tail audit", "figs/exec09_tail_comparison.png"),
        ])

        content = (
            preamble + title + geometry + faces + optical
            + "\\section{Key positions}\n" + key_position_slides
            + tN_synthesis_slide()
            + "\\section{Photon budget and statistics}\n" + integrated_exec12
            + "\\section{Window effect and timing}\n" + window_section_exec12
            + conclusions_exec12
            + "\\section{Full scan, position by position}\n" + full_position_slides
            + "\\appendix\n" + backup_exec12 + "\\end{document}\n"
        )
        _prefix = "exec12_report"
    else:
        content = (
            preamble + title + geometry + faces + optical
            + "\\section{Position-by-position}\n" + position_slides
            + "\\section{Photon budget and statistics}\n" + integrated
            + "\\section{Window effect and timing}\n" + landau_slide + dip + timing + timing_gate_slide + velocity_slide + top_slide
            + conclusions + "\\appendix\n" + backup + "\\end{document}\n"
        )
        _prefix = "exec11_report" if with_arrival else "exec10_report"
    tex_path = output_dir / f"{_prefix}_{mode}.tex"
    tex_path.write_text(content, encoding="utf-8")
    return tex_path


def compile_tex(tex_path: pathlib.Path) -> pathlib.Path:
    for _ in range(2):
        completed = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tex_path.parent,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if completed.returncode:
            raise RuntimeError(completed.stdout[-5000:])
    return tex_path.with_suffix(".pdf")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    parser.add_argument("--positions", choices=("full", "key", "both"), default="both")
    parser.add_argument(
        "--with-arrival", action="store_true", default=False,
        help=(
            "Interleave EXEC-11 photon-arrival companion slides after each "
            "position slide. Produces exec11_report_*.pdf; exec10_report_*.pdf "
            "is unchanged. Off by default."
        ),
    )
    parser.add_argument(
        "--with-pedagogy", action="store_true", default=False,
        help=(
            "Generate EXEC-12 pedagogical deck: T1 full-scan format, T2 window-track "
            "mechanism slides, T3 renamed timing-gate slide, T4 t_N figures, T5 "
            "Fano/Landau/readingbox layer, T6 extended glossary. Implies --with-arrival. "
            "Produces exec12_report_*.pdf. Off by default."
        ),
    )
    args = parser.parse_args()
    modes = ("full", "key") if args.positions == "both" else (args.positions,)
    for mode in modes:
        pdf = compile_tex(make_tex(
            args.output_dir, mode,
            with_arrival=args.with_arrival,
            with_pedagogy=args.with_pedagogy,
        ))
        print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
