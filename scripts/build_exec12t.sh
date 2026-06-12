#!/usr/bin/env bash
set -euo pipefail
out="${1:?results/exec12t_<timestamp> required}"
python3 analysis/exec12t_make_products.py --out "$out"
(cd "$out/report" && lualatex -halt-on-error -interaction=nonstopmode exec12t_timing_position_report.tex >/dev/null)
(cd "$out/report" && lualatex -halt-on-error -interaction=nonstopmode exec12t_timing_position_report.tex >/dev/null)
(cd "$out/beamer" && lualatex -halt-on-error -interaction=nonstopmode exec12t_timing_position_beamer.tex >/dev/null)
(cd "$out/beamer" && lualatex -halt-on-error -interaction=nonstopmode exec12t_timing_position_beamer.tex >/dev/null)
mkdir -p "$out/beamer/rendered"
pdftoppm -png -r 120 "$out/beamer/exec12t_timing_position_beamer.pdf" "$out/beamer/rendered/slide"
python3 scripts/make_contact_sheet.py "$out/beamer/rendered" "$out/beamer/beamer_contact_sheet.png"
pdfinfo "$out/beamer/exec12t_timing_position_beamer.pdf" > "$out/logs/beamer_pdfinfo.txt"
grep -E "Overfull|Underfull|Warning|Error" "$out/beamer/exec12t_timing_position_beamer.log" > "$out/logs/beamer_warnings.txt" || true
