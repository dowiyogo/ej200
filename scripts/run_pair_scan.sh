#!/usr/bin/env bash
# scripts/run_pair_scan.sh
# EXEC_08 pair-scan batch runner — for execution on t0minidaq.
#
# Resume-safe: skips positions whose .root output exists and is valid.
# Runs each position sequentially (multi-thread per run: /run/numberOfThreads 16).
# After all positions complete, invokes the ROOT analysis macro.
#
# Usage:
#   cd /home/reriosto/SHiP/ej200
#   nohup ./scripts/run_pair_scan.sh > pairscan.log 2>&1 &
#
# Environment overrides (export before running):
#   BUILD_DIR   — path to build directory (default: ./build)
#   MACRO_DIR   — path to macros (default: ./macros/pairscan)
#   OUTPUT_DIR  — path to results directory (default: ./results/pairscan_2026-06-11)
#   EVENTS      — events per position (default: 3000, set in macros)
#   WORKERS     — not used here (set per-macro via /run/numberOfThreads)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

BUILD_DIR="${BUILD_DIR:-$REPO_DIR/build}"
SIM="$BUILD_DIR/ej200_bar_sim"
MACRO_DIR="${MACRO_DIR:-$REPO_DIR/macros/pairscan}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_DIR/results/pairscan_2026-06-11}"
MASTER_LOG="$REPO_DIR/pairscan.log"

# Pair configuration (must match pair_scan_config.json)
ID_A=28
ID_B=29
EVENTS=3000

# ── Pre-flight checks ─────────────────────────────────────────────────────────
if [[ ! -x "$SIM" ]]; then
    echo "[ERROR] Simulator not found: $SIM" >&2
    echo "        Build with: mkdir build && cd build && cmake .. && make -j 24" >&2
    exit 1
fi

if [[ ! -d "$MACRO_DIR" ]]; then
    echo "[ERROR] Macro directory not found: $MACRO_DIR" >&2
    echo "        Run: python3 scripts/gen_pair_scan_macros.py" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

timestamp() { date '+%Y-%m-%dT%H:%M:%S'; }

log() { echo "[$(timestamp)] $*" | tee -a "$MASTER_LOG"; }

log "EXEC_08 pair scan starting"
log "  Build:   $BUILD_DIR"
log "  Macros:  $MACRO_DIR"
log "  Output:  $OUTPUT_DIR"
log "  Pair:    IDs ($ID_A, $ID_B)"
log "  Events:  $EVENTS"
log "  Commit:  $(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)"

# ── Validation helper ─────────────────────────────────────────────────────────
validate_root() {
    local path="$1"
    local expected_x="$2"
    python3 - "$path" "$expected_x" "$EVENTS" <<'PY'
import sys
import uproot
import numpy as np
path, expected_x, expected_events = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
with uproot.open(path) as rf:
    tree = rf["sipm_hits"]
    event_ids = set()
    positions = set()
    for chunk in tree.iterate(["event_id", "gun_x_mm"], step_size="100 MB", library="np"):
        event_ids.update(np.unique(chunk["event_id"]).tolist())
        positions.update(np.unique(chunk["gun_x_mm"]).tolist())
ok_events = len(event_ids) == expected_events
ok_pos    = abs(list(positions)[0] - expected_x) < 0.1 if len(positions) == 1 else False
if not (ok_events and ok_pos):
    raise SystemExit(1)
print(tree.num_entries)
PY
}

# ── Per-position loop ─────────────────────────────────────────────────────────
n_done=0
n_skip=0
n_error=0
macros=("$MACRO_DIR"/pairscan_x*.mac)
# Sort macros numerically by x value (lexicographic sort works for these filenames)
IFS=$'\n' sorted_macros=($(printf '%s\n' "${macros[@]}" | sort)); unset IFS

for mac in "${sorted_macros[@]}"; do
    base="$(basename "$mac" .mac)"
    final="$OUTPUT_DIR/${base}.root"
    # Extract x_true from filename: pairscan_x-462.0mm.mac -> -462.0
    x_str="${base#pairscan_x}"     # remove "pairscan_x" prefix
    x_str="${x_str%mm}"            # remove "mm" suffix
    x_val="$x_str"

    # Resume: skip if valid output already exists
    if [[ -f "$final" ]] && entries="$(validate_root "$final" "$x_val" 2>/dev/null)"; then
        log "SKIP $base (entries=$entries)"
        (( n_skip++ )) || true
        continue
    fi

    # Workdir under OUTPUT_DIR
    work="$OUTPUT_DIR/.work_${base}"
    rm -rf "$work"
    mkdir -p "$work"
    # Symlink sslg4 data directory if present
    if [[ -d "$BUILD_DIR/sslg4" ]]; then
        ln -sf "$BUILD_DIR/sslg4" "$work/sslg4"
    fi
    # Symlink data directory
    if [[ -d "$BUILD_DIR/data" ]]; then
        ln -sf "$BUILD_DIR/data" "$work/data"
    fi

    pos_log="$OUTPUT_DIR/run_${base}.log"
    log "START $base (x=${x_val} mm)"

    (
        cd "$work"
        "$SIM" -m "$mac" 2>&1 | tee "$pos_log.tmp"
    )
    mv "$pos_log.tmp" "$pos_log"

    candidate="$work/photon_hits_run000.root"
    if entries="$(validate_root "$candidate" "$x_val")"; then
        mv "$candidate" "$final.tmp"
        mv "$final.tmp" "$final"
        rm -rf "$work"
        log "DONE $base (entries=$entries)"
        (( n_done++ )) || true
    else
        log "ERROR $base — invalid output; retained in $work"
        (( n_error++ )) || true
    fi
done

log "Loop complete: done=$n_done skip=$n_skip error=$n_error"

if [[ $n_error -gt 0 ]]; then
    log "ERROR: $n_error positions failed — check logs in $OUTPUT_DIR"
    exit 1
fi

# ── ROOT analysis ─────────────────────────────────────────────────────────────
log "Starting analysis: analyze_pair_scan.C"
cd "$REPO_DIR"
root -l -b -q "analysis/analyze_pair_scan.C(\"$OUTPUT_DIR\")" 2>&1 | tee -a "$MASTER_LOG"
log "Analysis complete."
log "Results: $OUTPUT_DIR/analysis/"
