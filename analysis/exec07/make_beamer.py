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


def make_tex(output_dir: pathlib.Path, mode: str, with_arrival: bool = False) -> pathlib.Path:
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

    preamble = r"""\documentclass[aspectratio=169]{beamer}
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
\title{EXEC\_10: EndTop photon budget and intrinsic timing}
\subtitle{Landau-dominated Npe, clarified timing densities, and apparent-slope diagnosis}
\author{Rene Rios (ULS)}
\date{June 11, 2026}
\begin{document}
"""
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
    args = parser.parse_args()
    modes = ("full", "key") if args.positions == "both" else (args.positions,)
    for mode in modes:
        pdf = compile_tex(make_tex(args.output_dir, mode, with_arrival=args.with_arrival))
        print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
