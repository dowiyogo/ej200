#!/usr/bin/env bash
set -u
set -o pipefail

REPO=/home/rrios/ej204
RUN_DIR=${RUN_DIR:-/home/rrios/ej204/runs/t0minidaq_endtop_scan_20260618_204959}
OUT_DIR=${OUT_DIR:-$RUN_DIR/analysis_corefit}
BIN_PS=${BIN_PS:-10}
N_BOOTSTRAP=${N_BOOTSTRAP:-200}

cd "$REPO" || exit 1

mkdir -p "$OUT_DIR/logs"

{
    echo "============================================================"
    echo "  EndTop core-fit analysis"
    echo "  Start: $(date '+%F %T')"
    echo "  Host:  $(hostname)"
    echo "  Repo:  $REPO"
    echo "  Branch: $(git branch --show-current 2>/dev/null)"
    echo "  Commit: $(git rev-parse --short HEAD 2>/dev/null)"
    echo "  Run dir: $RUN_DIR"
    echo "  Out dir: $OUT_DIR"
    echo "  Bin ps:  $BIN_PS"
    echo "  Bootstrap: $N_BOOTSTRAP"
    echo "  ROOT: $(root-config --version 2>/dev/null || echo unknown)"
    echo "  Python: $(python3 --version 2>&1)"
    echo "============================================================"
} | tee "$OUT_DIR/logs/analysis.log"

python3 scripts/analyze_t0minidaq_endtop_corefit.py \
    --run-dir  "$RUN_DIR" \
    --out-dir  "$OUT_DIR" \
    --bin-ps   "$BIN_PS" \
    --n-bootstrap "$N_BOOTSTRAP" \
    2>&1 | tee -a "$OUT_DIR/logs/analysis.log"

status=${PIPESTATUS[0]}

{
    echo "============================================================"
    echo "  End: $(date '+%F %T')"
    echo "  Exit status: $status"
    echo "  Outputs: $OUT_DIR"
    echo "============================================================"
} | tee -a "$OUT_DIR/logs/analysis.log"

exit "$status"
