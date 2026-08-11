#!/usr/bin/env bash
# Compile EJ-230 Beamer report on MSI from downloaded results tree.
# Usage: ./scripts/build_report_msi.sh [path/to/results_ej230/]
set -euo pipefail

RESULTS_DIR="${1:-${RESULTS_DIR:-$HOME/ej230/results_ej230}}"
report_dir="$RESULTS_DIR/report"
tex_file="$report_dir/exec13_ej230_report_full.tex"

if [[ ! -f "$tex_file" ]]; then
    echo "ERROR: $tex_file not found." >&2
    echo "  Download results from t0minidaq first (scp the tarball), then unpack." >&2
    exit 1
fi

if ! command -v pdflatex &>/dev/null; then
    echo "ERROR: pdflatex not found. Install TeX Live on MSI." >&2
    exit 1
fi

echo "Compiling Beamer: $tex_file"
(
    cd "$report_dir"
    pdflatex -interaction=nonstopmode "$(basename "$tex_file")"
    pdflatex -interaction=nonstopmode "$(basename "$tex_file")"
)
pdf="$report_dir/exec13_ej230_report_full.pdf"
echo "Done: $pdf"
