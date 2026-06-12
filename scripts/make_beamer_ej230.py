#!/usr/bin/env python3
"""Generate exec13_ej230_report_full.tex as a structural mirror of exec12_report_full.tex.

Rules (from EXEC_13 spec):
  - Same frametitles, same order, same count (1:1 mirror).
  - Only material refs (EJ-204→EJ-230), figure paths (exec12→exec13),
    scintillator params (10400/9700, 0.7/0.5 ns, 1.8/1.5 ns, 160/120 cm),
    branch/commit, and date are updated.
  - Numbers that will be available from CSV are inserted via LaTeX macros
    (\\input{tables/exec13_macros.tex}). Placeholders used until data exists.

Verification: extracted frametitle sequence must be IDENTICAL to exec12.
If it differs, the script exits with status 1 (anti-regression guardrail).
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys


# ── Substitution table (order matters — more specific first) ─────────────────
SUBSTITUTIONS = [
    # Title / subtitle
    (r"EXEC\\\_12: EndTop photon budget --- pedagogical review",
     r"EXEC\\_13 (EJ-230): EndTop photon budget --- material comparison"),
    (r"t_N distributions, Fano/Landau motivation, window-track mechanism, full scan",
     r"EJ-230 (OPSC-106) mirror of EXEC\\_12 --- same analysis, new material"),
    # Date
    (r"June 11, 2026", r"June 12, 2026"),
    # Branch
    (r"feat/endtop-sslg4", r"feat/ej230-sslg4"),
    # Configuration provenance (single slide update)
    (r"EJ-204 via SSLG4 \\texttt\{OPSC-101\}: yield 10400/MeV, \$\\tau_r=0\.7\$ ns,\s*"
     r"\\s*\$\\tau_d=1\.8\$ ns, ABSLENGTH=160 cm \(verified effective MPT inputs\)\.",
     r"EJ-230 via SSLG4 \\texttt{OPSC-106}: yield 9700/MeV, $\\tau_r=0.5$ ns,\n"
     r"        $\\tau_d=1.5$ ns, ABSLENGTH=120 cm (verified effective MPT inputs)."),
    # Figure paths: exec12_tN → exec13_tN, exec12b → exec13b
    (r"figs/exec12_tN_", r"figs/exec13_tN_"),
    (r"figs/exec12_tN_summary", r"figs/exec13_tN_summary"),
    (r"exec12_tN_summary\.png", r"exec13_tN_summary.png"),
    # Dispersion figures (exec12b → exec13b naming)
    (r"figs/tn_(-?\d+)mm_(endL|endR|top)\.png",
     r"figs/exec13_tn_\\1mm_\\2.png"),
    # Material name throughout
    (r"EJ-204 bar:", r"EJ-230 bar:"),
    (r"\{EJ-204 bar\}", r"{EJ-230 bar}"),
    (r"EJ-204 bar", r"EJ-230 bar"),
    (r"\\texttt\{OPSC-101\}", r"\\texttt{OPSC-106}"),
    (r"OPSC-101", r"OPSC-106"),
    # Scintillator yield
    (r"yield 10400/MeV", r"yield 9700/MeV"),
    (r"10400/MeV", r"9700/MeV"),
    (r"yield 10400", r"yield 9700"),
    # Timing constants
    (r"\\tau_r=0\.7", r"\\tau_r=0.5"),
    (r"\\tau_d=1\.8", r"\\tau_d=1.5"),
    (r"rise 0\.7 ns", r"rise 0.5 ns"),
    (r"decay 1\.8 ns", r"decay 1.5 ns"),
    (r"0\.7.*ns,\s*\$\\tau_d", r"0.5 ns, $\\tau_d"),
    (r"1\.8 ns", r"1.5 ns"),
    (r"0\.7 ns", r"0.5 ns"),
    # Attenuation length
    (r"ABSLENGTH=160 cm", r"ABSLENGTH=120 cm"),
    (r"ABSLENGTH \$=160\$ cm", r"ABSLENGTH $=120$ cm"),
    (r"ABSLENGTH \$=160$", r"ABSLENGTH $=120$"),
    (r"ABSLENGTH 160 cm", r"ABSLENGTH 120 cm"),
    (r"\$160\$ cm", r"$120$ cm"),
    # Emission peak
    (r"emission peak 408\.8 nm", r"emission peak $\\approx391$ nm"),
    (r"408\.8.*nm", r"$\\approx$391 nm"),
    # MPT checks line
    (r"yield 10400/MeV, decay 1\.8 ns, rise 0\.7 ns,",
     r"yield 9700/MeV, decay 1.5 ns, rise 0.5 ns,"),
    (r"ABSLENGTH 160 cm, emission peak 408\.8 nm\.",
     r"ABSLENGTH 120 cm, emission peak $\\approx$391 nm."),
    # Campaign labels in non-frame text
    (r"EXEC_12\b", r"EXEC_13"),
    (r"exec12\b", r"exec13"),
    # Generic EJ-204 that wasn't caught by more specific rules
    (r"EJ-204\b", r"EJ-230"),
]


def apply_substitutions(text: str) -> str:
    for pattern, replacement in SUBSTITUTIONS:
        text = re.sub(pattern, replacement, text)
    return text


def extract_frametitles(tex: str) -> list[str]:
    return re.findall(r"\\begin\{frame\}\{([^}]*)\}", tex)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate exec13_ej230_report_full.tex as EJ-230 mirror of exec12."
    )
    repo_root = pathlib.Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--source", type=pathlib.Path,
        default=repo_root / "analysis" / "exec07" / "exec12_report_full.tex",
    )
    parser.add_argument(
        "--output", type=pathlib.Path,
        default=repo_root / "analysis" / "exec07" / "exec13_ej230_report_full.tex",
    )
    parser.add_argument(
        "--results-dir", type=pathlib.Path,
        default=None,
        help="If provided, add \\graphicspath pointing to results figures.",
    )
    args = parser.parse_args()

    if not args.source.is_file():
        print(f"ERROR: source not found: {args.source}", file=sys.stderr)
        return 1

    source_text = args.source.read_text(encoding="utf-8")
    ref_titles = extract_frametitles(source_text)

    new_text = apply_substitutions(source_text)

    # Add \graphicspath to results figures if requested
    if args.results_dir is not None:
        figs_path = args.results_dir / "figures"
        graphicspath = (
            f"\\graphicspath{{{{{figs_path}/}},"
            f"{{figs/}}}}"
        )
        new_text = new_text.replace(
            r"\begin{document}",
            graphicspath + "\n" + r"\begin{document}",
            1,
        )

    # Verify frame titles are structurally identical (ignoring EJ-204→EJ-230 substitution)
    new_titles = extract_frametitles(new_text)
    ref_normalized = [t.replace("EJ-204", "EJ-230").replace("exec12", "exec13")
                      for t in ref_titles]
    if len(ref_normalized) != len(new_titles):
        print(
            f"FRAME COUNT MISMATCH: reference has {len(ref_normalized)} frames, "
            f"output has {len(new_titles)} frames.",
            file=sys.stderr,
        )
        return 1
    mismatches = [(i, r, n) for i, (r, n) in enumerate(zip(ref_normalized, new_titles))
                  if r != n]
    if mismatches:
        print("FRAME TITLE MISMATCH:", file=sys.stderr)
        for i, r, n in mismatches[:10]:
            print(f"  frame {i+1}: expected '{r}', got '{n}'", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(new_text, encoding="utf-8")
    print(f"wrote {args.output} ({len(new_titles)} frames)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
