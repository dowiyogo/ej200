#!/usr/bin/env bash
# Regenerate EJ-230 photon-budget products affected by inherited material labels/constants.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
results="${RESULTS_DIR:-$repo/results_ej230_analysis}"
data="${MAIN_DATA_DIR:-/home/reriosto/SHiP/t0minidaq/results_ej230/data}"
logs="$results/logs"
figs="$results/figs"
mkdir -p "$logs" "$figs"

python3 "$repo/analysis/exec07_photon_budget.py" \
    --data-dir "$data" \
    --output-dir "$results" \
    --skip-channel-pdfs \
    2>&1 | tee "$logs/exec14b_photon_budget_repair.log"

python3 "$repo/analysis/exec07/exec12b_tn_dispersion.py" \
    --data-dir "$data" \
    --output-dir "$results" \
    --tau-d 1.5 \
    2>&1 | tee "$logs/exec14b_dispersion_repair.log"

for x_mm in -690 -650 -400 0 400 650 690; do
    for panel in top endL endR; do
        cp "$figs/tn_${x_mm}mm_${panel}.png" "$figs/exec13_tn_${x_mm}mm_${panel}.png"
    done
done

echo "EXEC_14B main-analysis repair complete: EJ-230 labels/constants and 21 dispersion panels"
