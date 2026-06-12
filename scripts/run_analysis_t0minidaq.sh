#!/usr/bin/env bash
# EXEC_13 EJ-230 — Full analysis pipeline for t0minidaq (headless, no X).
# Runs: ROOT → CSVs → figures (PNG) → LaTeX tables → Beamer .tex
# If pdflatex is available, also compiles the Beamer PDF.
# At the end, generates a tarball of results_ej230/ and prints the scp command.
set -euo pipefail

RESULTS_DIR="${RESULTS_DIR:-$HOME/ej230/results_ej230}"
export RESULTS_DIR

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

data_dir="$RESULTS_DIR/data"
csv_dir="$RESULTS_DIR/csv"
figs_dir="$RESULTS_DIR/figures"
tables_dir="$RESULTS_DIR/tables"
report_dir="$RESULTS_DIR/report"
logs_dir="$RESULTS_DIR/logs"

mkdir -p "$csv_dir" "$figs_dir" "$tables_dir" "$report_dir" "$logs_dir"

echo "======================================================================"
echo " EXEC_13 / EJ-230 analysis pipeline  ($(date '+%Y-%m-%d %H:%M:%S'))"
echo " RESULTS_DIR = $RESULTS_DIR"
echo "======================================================================"

# ── Step 0: Python dependency check / install ─────────────────────────────
echo ""
echo "--- Step 0: Python dependencies ---"
REQUIRED=(uproot numpy pandas matplotlib scipy awkward)
MISSING=()
for pkg in "${REQUIRED[@]}"; do
    python3 -c "import $pkg" 2>/dev/null || MISSING+=("$pkg")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "Installing missing packages: ${MISSING[*]}"
    python3 -m pip install --user "${MISSING[@]}"
else
    echo "All required Python packages present."
fi

# Verify matplotlib runs headless
python3 -c "
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0],[0])
plt.close(fig)
print('matplotlib Agg backend: OK')
"

# ── Step 1: t_N analysis (ROOT → CSV + figures) ───────────────────────────
echo ""
echo "--- Step 1: t_N analysis (ROOT → CSV + PNG) ---"
python3 "$REPO_DIR/analysis/exec13/exec13_tN_analysis.py" \
    --data-dir "$data_dir" \
    --figs-dir "$figs_dir" \
    --csv-dir  "$csv_dir"  \
    2>&1 | tee "$logs_dir/exec13_tN_analysis.log"

# ── Step 2: LaTeX table generation from CSVs ─────────────────────────────
echo ""
echo "--- Step 2: LaTeX tables from CSV ---"
python3 - "$csv_dir" "$tables_dir" <<'PY'
import sys, pathlib, math
import pandas as pd

csv_dir    = pathlib.Path(sys.argv[1])
tables_dir = pathlib.Path(sys.argv[2])
tables_dir.mkdir(parents=True, exist_ok=True)

csv_path = csv_dir / "exec13_tN_summary.csv"
if not csv_path.is_file():
    print(f"WARNING: {csv_path} not found; skipping table generation.")
    sys.exit(0)

df = pd.read_csv(csv_path)

# Table: sigma_fit(t_4) for all key positions and groups
sub4 = df[df["N"] == 4].copy()
table_lines = [
    r"\begin{tabular}{lrrrr}",
    r"\hline",
    r"$x$ [mm] & Group & $\langle t_4 \rangle$ [ps] & RMS [ps] & $\sigma_\mathrm{fit}$ [ps] \\",
    r"\hline",
]
for _, row in sub4.iterrows():
    mean_ps  = f"{row['mean_ns']*1000:.0f}" if math.isfinite(row['mean_ns'])  else "---"
    rms_ps   = f"{row['rms_ns']*1000:.0f}"  if math.isfinite(row['rms_ns'])   else "---"
    sigma_ps = f"{row['sigma_fit_ns']*1000:.0f}" if math.isfinite(row['sigma_fit_ns']) else "---"
    table_lines.append(
        f"  {row['x_mm']:+d} & {row['group']} & {mean_ps} & {rms_ps} & {sigma_ps} \\\\"
    )
table_lines += [r"\hline", r"\end{tabular}"]
out_path = tables_dir / "exec13_tN_t4_table.tex"
out_path.write_text("\n".join(table_lines) + "\n")
print(f"  wrote {out_path}")

# Single-value macros (used by \input{} in Beamer)
macros_path = tables_dir / "exec13_macros.tex"
with macros_path.open("w") as f:
    for _, row in sub4.iterrows():
        safe_x = str(row['x_mm']).replace('-', 'neg').replace('+', '')
        safe_g = row['group'].replace('_', '')
        for col, suffix in [('mean_ns', 'meanps'), ('sigma_fit_ns', 'sigmafit')]:
            val = row[col]
            if math.isfinite(val):
                f.write(rf"\newcommand{{\exec13{safe_g}x{safe_x}{suffix}}}"
                        rf"{{{val*1000:.0f}}}" + "\n")
print(f"  wrote {macros_path}")
PY

# ── Step 3: Generate Beamer .tex ─────────────────────────────────────────
echo ""
echo "--- Step 3: Beamer .tex generation ---"
python3 "$REPO_DIR/scripts/make_beamer_ej230.py" \
    --results-dir "$RESULTS_DIR" \
    --output "$report_dir/exec13_ej230_report_full.tex" \
    2>&1 | tee "$logs_dir/make_beamer_ej230.log" || {
    echo "WARNING: make_beamer_ej230.py failed or not yet available; skipping."
}

# Copy static Beamer .tex from repo if the generator did not produce it yet
if [[ ! -f "$report_dir/exec13_ej230_report_full.tex" ]]; then
    static_tex="$REPO_DIR/analysis/exec07/exec13_ej230_report_full.tex"
    if [[ -f "$static_tex" ]]; then
        cp "$static_tex" "$report_dir/exec13_ej230_report_full.tex"
        echo "  copied static tex to $report_dir/exec13_ej230_report_full.tex"
    fi
fi

# ── Step 4: Compile Beamer if pdflatex available ─────────────────────────
echo ""
echo "--- Step 4: Beamer compilation ---"
TEX_FILE="$report_dir/exec13_ej230_report_full.tex"
if [[ -f "$TEX_FILE" ]] && command -v pdflatex &>/dev/null; then
    echo "pdflatex found — compiling..."
    (
        cd "$report_dir"
        pdflatex -interaction=nonstopmode "$(basename "$TEX_FILE")" \
            2>&1 | tee "$logs_dir/pdflatex1.log"
        pdflatex -interaction=nonstopmode "$(basename "$TEX_FILE")" \
            2>&1 | tee "$logs_dir/pdflatex2.log"
    )
    echo "  PDF: $report_dir/exec13_ej230_report_full.pdf"
elif [[ -f "$TEX_FILE" ]]; then
    echo "BEAMER PENDIENTE: compilar en MSI con scripts/build_report_msi.sh"
    echo "  (pdflatex not found on this machine)"
else
    echo "  No .tex found; generation step incomplete."
fi

# ── Step 5: Tarball for download ─────────────────────────────────────────
echo ""
echo "--- Step 5: Tarball for download ---"
tarball="/home/rrios/ej230/results_ej230_$(date +%Y%m%d_%H%M).tar.gz"
# If not running as rrios, use generic path
if [[ ! -d "/home/rrios/ej230" ]]; then
    tarball="$HOME/ej230/results_ej230_$(date +%Y%m%d_%H%M).tar.gz"
fi
tar -czf "$tarball" -C "$(dirname "$RESULTS_DIR")" "$(basename "$RESULTS_DIR")"
echo ""
echo "======================================================================"
echo " Analysis complete."
echo " Tarball: $tarball"
echo ""
echo " To download from MSI:"
echo "   scp rrios@t0minidaq.cern.ch:$tarball ."
echo "======================================================================"
