#!/usr/bin/env bash
# Strict compile and machine-readable preflight for the standalone EJ-230 report.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
results="$repo/results_ej230_analysis"
report="$results/report"
logs="$results/logs"
tex="exec13_ej230_report_full.tex"
pdf="$report/exec13_ej230_report_full.pdf"
mkdir -p "$logs"

python3 "$repo/scripts/audit_beamer_assets.py" \
    2>&1 | tee "$logs/exec14b_asset_audit_after.log"
python3 "$repo/scripts/audit_beamer_numbers.py" \
    2>&1 | tee "$logs/exec14b_number_audit.log"
python3 "$repo/scripts/check_exec14b_frame_parity.py" \
    2>&1 | tee "$logs/exec14b_frame_parity.log"
python3 "$repo/scripts/check_exec14b_asset_parity.py" \
    2>&1 | tee "$logs/exec14b_asset_parity.log"
python3 "$repo/scripts/check_exec14b_figure_frame_consistency.py" \
    2>&1 | tee "$logs/exec14b_figure_frame_consistency.log"
python3 "$repo/scripts/audit_exec14b_raster_text.py" \
    2>&1 | tee "$logs/exec14b_raster_text_audit.log"

duplicates="$(
    grep '^\\newcommand' "$results"/tables/*.tex \
        | sed 's/.*{\\\([^}]*\)}.*/\1/' \
        | sort \
        | uniq -d
)"
if [[ -n "$duplicates" ]]; then
    echo "ERROR: duplicate macros:" >&2
    echo "$duplicates" >&2
    exit 1
fi

(
    cd "$report"
    latexmk -pdf \
        -interaction=nonstopmode \
        -halt-on-error \
        -file-line-error \
        "$tex"
) 2>&1 | tee "$logs/exec14b_latexmk.log"

fatal_pattern='File .* not found|using draft setting|Undefined control sequence|Command .* already defined|There were undefined references|LaTeX Warning: Reference .* undefined|\\[1-9]'
if rg -n "$fatal_pattern" "$logs/exec14b_latexmk.log"; then
    echo "ERROR: fatal pattern found in strict LaTeX log" >&2
    exit 1
fi

pdftotext "$pdf" - > "$logs/exec14b_pdftotext.txt"
if grep -F 'figs/' "$logs/exec14b_pdftotext.txt"; then
    echo "ERROR: literal figs/ path found in final PDF" >&2
    exit 1
fi
if grep -E 'exec13_tn_.*\\[12]|File.*not found' "$logs/exec14b_pdftotext.txt"; then
    echo "ERROR: broken path text found in final PDF" >&2
    exit 1
fi

pdfimages -list "$pdf" | tee "$logs/exec14b_pdfimages.txt"
if [[ "$(wc -l < "$logs/exec14b_pdfimages.txt")" -le 2 ]]; then
    echo "ERROR: pdfimages found no embedded images" >&2
    exit 1
fi

python3 "$repo/scripts/preflight_exec14b_pdf.py" \
    2>&1 | tee "$logs/exec14b_pdf_preflight.log"
echo "EXEC_14B strict preflight PASS"
