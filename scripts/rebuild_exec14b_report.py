#!/usr/bin/env python3
"""Repair paths and replace empirical Beamer text with generated CSV fragments."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex"
KEY_POSITIONS = {-690, -650, -400, 0, 400, 650, 690}
DISPERSION = {8: -690, 12: -650, 16: -400, 20: 0, 24: 400, 28: 650, 32: 690}


def frames(text: str) -> list[tuple[int, int, int, int]]:
    result = []
    for match in re.finditer(r"\\begin\{frame\}", text):
        end = text.find(r"\end{frame}", match.end())
        if end < 0:
            raise RuntimeError(f"unterminated frame at offset {match.start()}")
        body_start = text.find("\n", match.end()) + 1
        result.append((match.start(), body_start, end, end + len(r"\end{frame}")))
    return result


def replace_body(text: str, frame_number: int, body: str) -> str:
    items = frames(text)
    _, body_start, body_end, _ = items[frame_number - 1]
    return text[:body_start] + "\n" + body.rstrip() + "\n\n" + text[body_end:]


def dispersion_body(x_mm: int) -> str:
    stem = f"figs/exec13_tn_{x_mm}mm"
    return rf"""\begin{{columns}}[T]
\column{{0.60\textwidth}}
\centering
\includegraphics[width=\textwidth,height=0.78\textheight,keepaspectratio]%
{{{stem}_top.png}}
\column{{0.38\textwidth}}
\centering
\includegraphics[width=\textwidth,height=0.37\textheight,keepaspectratio]%
{{{stem}_endL.png}}\\[2pt]
\includegraphics[width=\textwidth,height=0.37\textheight,keepaspectratio]%
{{{stem}_endR.png}}
\end{{columns}}
\par\scriptsize Error bars = RMS (event-to-event). Dashed band: statistical floor
$\sigma_{{stat}}(t_n)\approx\sqrt{{n}}\,\tau_d/\langle N_{{pe}}\rangle$
($\tau_d=1.5$ ns). See dispersion-decomposition frame for physical interpretation."""


def position_inputs(text: str) -> str:
    for frame_number in list(range(5, 30, 4)) + list(range(52, 113, 2)):
        items = frames(text)
        _, start, end, _ = items[frame_number - 1]
        body = text[start:end]
        match = re.search(r"figs/muon_(-?\d+)mm_geometry\.png", body)
        if not match:
            raise RuntimeError(f"frame {frame_number}: cannot determine position")
        x_mm = int(match.group(1))
        filename = f"key_position_{x_mm}.tex" if x_mm in KEY_POSITIONS else f"position_{x_mm}.tex"
        pattern = re.compile(
            r"  \\tiny\s+  \\begin\{tabular\}\{lr\}.*?  \\end\{tabular\}", re.DOTALL
        )
        replacement = rf"  \tiny" + "\n" + rf"  \input{{../tables/{filename}}}"
        body, count = pattern.subn(lambda _: replacement, body)
        if count == 0 and rf"\input{{../tables/{filename}}}" not in body:
            raise RuntimeError(f"frame {frame_number}: position table not found")
        text = text[:start] + body + text[end:]
    return text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    text = args.report.read_text(encoding="utf-8")

    text = text.replace(r"\graphicspath{{../figs/}}", r"\graphicspath{{../}}", 1)
    if r"\input{../tables/exec14_macros.tex}" not in text:
        text = text.replace(
            r"\usepackage{amsmath}",
            "\\usepackage{amsmath}\n\\input{../tables/exec14_macros.tex}",
            1,
        )
    text = re.sub(
        r"\\title\{.*?\}\n\\subtitle\{.*?\}\n",
        lambda _: "\\title{EJ-230 EndTop photon budget and intrinsic timing}\n"
        "\\subtitle{Standalone OPSC-106 analysis report}\n",
        text,
        count=1,
    )
    text = text.replace(r"\date{June 12, 2026}", r"\date{June 13, 2026}", 1)

    text = position_inputs(text)
    text = replace_body(text, 1, r"\input{../tables/provenance.tex}")
    for frame_number, x_mm in DISPERSION.items():
        text = replace_body(text, frame_number, dispersion_body(x_mm))

    frame33 = r"""\centering\includegraphics[width=\textwidth,height=0.65\textheight,keepaspectratio]%
{figs/exec13_tN_summary.png}
\input{../tables/tN_summary.tex}"""
    text = replace_body(text, 33, frame33)
    text = text.replace(r"\lambda_{\rm eff}\approx33", r"\lambda_{\rm eff}\approx\LambdaEffMean")
    text = text.replace(r"\lambda_{eff}\approx33", r"\lambda_{eff}\approx\LambdaEffMean")
    text = replace_body(
        text,
        36,
        r"""\centering\includegraphics[width=\textwidth,height=0.62\textheight,keepaspectratio]%
{figs/P2_npe_heatmap_top.png}
\input{../tables/top_localization_summary.tex}""",
    )

    _, start37, end37, _ = frames(text)[36]
    body37 = text[start37:end37]
    body37 = re.sub(
        r"\\par\\scriptsize \$\\lambda_\{eff,L\}=.*?propagation velocity\.",
        lambda _: r"\par\scriptsize \input{../tables/attenuation_fit.tex}",
        body37,
        count=1,
    )
    body37 = body37.replace(r"$\lambda_{eff}\approx31$ cm", r"$\lambda_{eff}\approx\LambdaEffMean$ cm")
    text = text[:start37] + body37 + text[end37:]

    _, start38, end38, _ = frames(text)[37]
    body38 = text[start38:end38]
    body38 = body38.replace(r"$\sqrt{c}\approx0.249$", r"$\sqrt{c}\approx\FanoSqrtC$")
    body38 = body38.replace(r"$c = 0.062547\pm0.000117$", r"$c = \FanoC\pm\FanoCerr$")
    text = text[:start38] + body38 + text[end38:]

    _, start39, end39, _ = frames(text)[38]
    body39 = text[start39:end39]
    body39 = re.sub(
        r"\\par\\scriptsize \$F=1\+c\\langle N_\{pe\}\\rangle\$.*?low light\.",
        lambda _: r"\par\scriptsize \input{../tables/fano_fit.tex}",
        body39,
        count=1,
    )
    text = text[:start39] + body39 + text[end39:]

    text = replace_body(text, 44, r"\input{../tables/landau_key_positions.tex}")
    _, start45, end45, _ = frames(text)[44]
    body45 = text[start45:end45]
    body45 = re.sub(
        r"  \\item Empirically:.*?symmetric in mirror \(C2\)\.",
        lambda _: r"  \input{../tables/window_dip_mechanism.tex}",
        body45,
        count=1,
        flags=re.DOTALL,
    )
    text = text[:start45] + body45 + text[end45:]

    text = replace_body(
        text,
        46,
        r"""\begin{columns}[T]\column{0.47\textwidth}
\scriptsize\resizebox{\textwidth}{!}{\begin{tabular}{lrr}\toprule
Run & contrast [Npe] & significance \\ \midrule
\input{../tables/window_dip_test.tex}
\bottomrule\end{tabular}}
\column{0.51\textwidth}
\includegraphics[width=\textwidth]{figs/exec08b_window_dip_profiles.png}
\end{columns}
\begin{block}{Finding}
\input{../tables/window_dip_summary.tex}
\end{block}""",
    )
    text = replace_body(
        text,
        47,
        r"""\centering\includegraphics[width=\textwidth,height=0.60\textheight,keepaspectratio]%
{figs/P7_deltaT_end.png}
\par\scriptsize $\sigma_{group}=\sigma(\Delta T_{AB})/\sqrt{2}$; no SPTR or electronics jitter.
\begin{block}{\footnotesize CSV-derived timing summary}
\scriptsize\input{../tables/sum4_timing.tex}
\end{block}""",
    )
    text = replace_body(
        text,
        48,
        r"""\begin{columns}[T]
\column{0.48\textwidth}
\scriptsize\resizebox{\textwidth}{!}{\begin{tabular}{rrrrr}\toprule
$x$ & side & EndTop [ps] & End-only [ps] & ratio \\ \midrule
\input{../tables/endtop_endonly_ratio.tex}
\bottomrule\end{tabular}}
\vspace{2pt}
\begin{alertblock}{\footnotesize Confirmed mechanism}
\scriptsize\input{../tables/endtop_endonly_tails.tex}
\end{alertblock}
\column{0.50\textwidth}
\includegraphics[width=0.88\textwidth,height=0.42\textheight,keepaspectratio]{figs/exec09_tail_comparison.png}
\begin{block}{\footnotesize No time acceptance window}
\scriptsize All End-photon arrivals are used without a $t<t_{max}$ cut.
EndTop tail narrowing persists over the full range: physical, not an artefact.
\end{block}
\end{columns}""",
    )
    text = replace_body(
        text,
        49,
        r"""\begin{columns}[T]\column{0.62\textwidth}
\includegraphics[width=\textwidth]{figs/exec10_velocity_estimators.png}
\column{0.36\textwidth}
\scriptsize\resizebox{\textwidth}{!}{\begin{tabular}{lrr}\toprule
Estimator & slope-derived [cm/ns] & $\chi^2/ndf$ \\ \midrule
\input{../tables/velocity_estimators.tex}
\bottomrule\end{tabular}}
\begin{alertblock}{Refuted predictions}
Every global line fit is rejected; slope-derived values are not accepted propagation velocities.
\end{alertblock}
\end{columns}""",
    )
    text = replace_body(
        text,
        50,
        r"""\begin{columns}[T]
\column{0.42\textwidth}
\scriptsize\setlength{\tabcolsep}{4pt}
\begin{tabular}{rrr}\toprule
$x$ [mm] & $N_{pe}$ & $\sigma_{t,est}$ [ps] \\ \midrule
\input{../tables/top_timing_estimates.tex}
\bottomrule\end{tabular}
\column{0.56\textwidth}
\begin{block}{\footnotesize Definition, derivation, and scope}
\scriptsize $\sigma_{t,est}=\sqrt{\tau_r\tau_d}/\sqrt{N_{pe}}$.
Threshold at low $N_{pe}$, double-exponential emission; analytic orientation only.
\end{block}
\begin{block}{\footnotesize $\sigma_{t,est}$ vs.\ simulated $\sigma(t_4)$}
\tiny\setlength{\tabcolsep}{3pt}\begin{tabular}{rcc}\toprule
$x$ & $\sigma_{t,est}$ [ps] & fitted $\sigma(t_4)$ [ps]\\ \midrule
\input{../tables/top_timing_comparison.tex}
\bottomrule\end{tabular}
\end{block}
\end{columns}""",
    )
    text = replace_body(text, 51, r"\input{../tables/numerical_conclusions.tex}")
    text = replace_body(text, 114, r"\input{../tables/validation_summary.tex}")
    text = replace_body(
        text,
        116,
        r"""\begin{enumerate}\small
 \item \textbf{Noise rejection at hardware parity:} 4 PE is the nominal fixed threshold outside
       the correlated-noise spectrum.
 \item \textbf{Slope optimum:} lower thresholds expose first-photon order statistics; higher
       thresholds reduce slope and increase walk.
 \item \textbf{Far-end efficiency:} \input{../tables/threshold_rationale.tex}
 \item \textbf{Walk:} no time-walk correction is applied in this campaign.
\end{enumerate}""",
    )
    text = replace_body(
        text,
        117,
        r"""\begin{itemize}
 \item SUM4 groups: \{0..3\}, \{4..7\}, \{8..11\}, \{12..15\}.
 \item Single-PE response: normalized 0.5/5 ns double exponential.
 \item Absolute leading-edge threshold: nominal 4 PE equivalent; the CSV-derived
       far-end rationale is reported on the preceding backup slide.
 \item Estimator: $\sigma(\Delta T_{AB})/\sqrt{2}$.
 \item No time-walk correction is applied in this campaign.
 \item The EndTop/End-only anti-artifact audit uses a common threshold, pulse,
       primary time convention and estimator.
\end{itemize}""",
    )
    _, start118, end118, _ = frames(text)[117]
    body118 = text[start118:end118]
    body118 = re.sub(
        r" \\item\[HOOK\\_WALK\].*?\(\\texttt\{sum4\\_timing\.C:305--315\}\)\.",
        lambda _: r" \item[Time-walk] No time-walk correction is applied in the EJ-230 campaign.",
        body118,
        count=1,
        flags=re.DOTALL,
    )
    body118 = re.sub(
        r" \\item\[Intrinsic sigma\].*?directly comparable\.",
        lambda _: r" \item[Intrinsic sigma] Excludes SPTR and electronics jitter.",
        body118,
        count=1,
        flags=re.DOTALL,
    )
    body118 = re.sub(
        r" Excluded from all EXEC\\_12 results.*?electronics jitter\.",
        lambda _: r" Excluded from all EJ-230 results; all reported timing is intrinsic.",
        body118,
        count=1,
        flags=re.DOTALL,
    )
    body118 = re.sub(
        r" \\item\[Window-track effect\].*?Geometric working hypothesis\.",
        lambda _: r" \item[Window-track effect] \input{../tables/window_dip_summary.tex}",
        body118,
        count=1,
        flags=re.DOTALL,
    )
    text = text[:start118] + body118 + text[end118:]

    if re.search(r"\\[1-9]", text):
        raise RuntimeError("residual regex backreference remains after rebuild")
    if len(frames(text)) != 119:
        raise RuntimeError(f"frame count changed: {len(frames(text))}")
    args.report.write_text(text, encoding="utf-8")
    print(f"rebuilt {args.report} with {len(frames(text))} frames")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
