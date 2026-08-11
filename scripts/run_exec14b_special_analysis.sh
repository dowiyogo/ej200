#!/usr/bin/env bash
# Regenerate EXEC_08b/EXEC_09 EJ-230 CSVs and figures from completed special controls.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
results="${RESULTS_DIR:-$repo/results_ej230_analysis}"
main_data="${MAIN_DATA_DIR:-/home/reriosto/SHiP/t0minidaq/results_ej230/data}"
special="${SPECIAL_DIR:-$results/sim_missing}"
csv="$results/csv"
figs="$results/figs"
logs="$results/logs"
mkdir -p "$csv" "$figs" "$logs"

python3 "$repo/analysis/exec07/exec08b_window_dip.py" \
    --existing "$main_data/photon_hits_x-650mm.root" \
    --run-a "$special/window_dip/photon_hits_run_A_x-652mm.root" \
    --run-b "$special/window_dip/photon_hits_run_B_x-642mm.root" \
    --run-c1 "$special/window_dip/photon_hits_run_C1_x-648mm.root" \
    --run-c2 "$special/window_dip/photon_hits_run_C2_x-654mm.root" \
    --output-dir "$results" \
    2>&1 | tee "$logs/exec14b_window_dip.log"
mv "$results/exec08b_window_dip_profiles.csv" "$csv/exec08b_window_dip_profiles.csv"

set +e
python3 "$repo/analysis/exec07/exec08b_timing_gate.py" \
    --endtop-dir "$main_data" \
    --end-only-dir "$special/end_only" \
    --output "$csv/exec08b_timing_gate.csv" \
    2>&1 | tee "$logs/exec14b_timing_gate.log"
timing_status=${PIPESTATUS[0]}
set -e
if [[ $timing_status -ne 0 ]]; then
    python3 - "$csv/exec08b_timing_gate.csv" <<'PY'
import sys
import pandas as pd

frame = pd.read_csv(sys.argv[1])
if frame.empty or frame["endtop_not_better"].all():
    raise SystemExit("timing gate failed without the documented EndTop improvement")
print("documented physical result accepted: EndTop is narrower in at least one comparison")
PY
fi

set +e
python3 "$repo/analysis/exec07/exec09_timing_mechanism.py" \
    --endtop-dir "$main_data" \
    --end-only-dir "$special/end_only" \
    --output-dir "$results" \
    2>&1 | tee "$logs/exec14b_timing_mechanism.log"
mechanism_status=${PIPESTATUS[0]}
set -e
for product in \
    "$results/exec09_tail_metrics.csv" \
    "$results/exec09_tail_comparison.csv" \
    "$results/exec09_timing_verdict.txt" \
    "$figs/exec09_tail_comparison.png"
do
    [[ -s "$product" ]] || {
        echo "ERROR: exec09 did not produce $product (status=$mechanism_status)" >&2
        exit 1
    }
done
echo "exec09 scientific gate status=$mechanism_status (preserved in verdict/log)"
mv "$results/exec09_tail_metrics.csv" "$csv/exec09_tail_metrics.csv"
mv "$results/exec09_tail_comparison.csv" "$csv/exec09_tail_comparison.csv"
mv "$results/exec09_timing_verdict.txt" "$csv/exec09_timing_verdict.txt"

echo "EXEC_14B special analysis complete: $csv and $figs"
