#!/usr/bin/env bash
# EJ-230 (EXEC_13) position scan — parametrized by --threads N.
# Wrappers with explicit machine names: run_scan_t0minidaq_24t.sh (24t)
#                                       run_scan_msi_16t.sh        (16t)
set -euo pipefail

# ── Configurable parameters ──────────────────────────────────────────────────
THREADS=16
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

repo="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="${BUILD_DIR:-$repo/build}"
sim="$build_dir/ej200_bar_sim"

RESULTS_DIR="${RESULTS_DIR:-$HOME/ej230/results_ej230}"
output_dir="$RESULTS_DIR/data"
log_dir="$RESULTS_DIR/logs"

events="${EVENTS_PER_POINT:-2000}"
positions=(-690 -670 -650 -600 -550 -500 -450 -400 -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350 400 450 500 550 600 650 670 690)

mkdir -p "$output_dir" "$log_dir"

validate_root() {
    python3 - "$1" "$2" "$events" <<'PY'
import sys
import numpy as np
import uproot

path, expected_x, expected_events = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
with uproot.open(path) as root_file:
    tree = root_file["sipm_hits"]
    event_ids = set()
    positions_seen = set()
    for chunk in tree.iterate(["event_id", "gun_x_mm"], step_size="200 MB", library="np"):
        event_ids.update(np.unique(chunk["event_id"]).tolist())
        positions_seen.update(np.unique(chunk["gun_x_mm"]).tolist())
if len(event_ids) != expected_events or positions_seen != {expected_x}:
    raise SystemExit(1)
print(tree.num_entries)
PY
}

# Each resume creates a timestamped subdirectory under data/ so ROOT files are
# never overwritten (rule: corrupt ROOT preserved, never deleted).
session_dir="$output_dir/session_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$session_dir"

echo "=== EJ-230 scan (EXEC_13) | threads=$THREADS events=$events ==="
echo "    sim:     $sim"
echo "    results: $RESULTS_DIR"
echo "    session: $session_dir"

for index in "${!positions[@]}"; do
    x="${positions[$index]}"
    # Final file lives in data/ (not session subdirectory) to allow resume detection
    final="$output_dir/photon_hits_x${x}mm.root"
    work="$session_dir/.work_x${x}mm"
    log="$log_dir/run_x${x}mm.log"

    if [[ -f "$final" ]] && entries="$(validate_root "$final" "$x" 2>/dev/null)"; then
        echo "DONE (cached): x=${x} mm entries=${entries}"
        continue
    fi

    rm -rf "$work"
    mkdir -p "$work"
    ln -sf "$build_dir/sslg4" "$work/sslg4"
    seed1=$((710000 + 1009 * index))
    seed2=$((910000 + 2017 * index))

    cat > "$work/run.mac" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads $THREADS
/det/readout EndTop
/det/scintillator OPSC-106
/sipm/model AFBR-S4N66P024M
/run/initialize
/sipm/jitterSigma 0 ns
/gun/particle mu-
/gun/energy 1 GeV
/muon/angle 0
/random/setSeeds $seed1 $seed2
/muon/gunX $x mm
/run/beamOn $events
EOF

    echo "START: x=${x} mm events=${events} threads=${THREADS}"
    (
        cd "$work"
        "$sim" -m run.mac 2>&1 | tee "$log.tmp"
    )
    mv "$log.tmp" "$log"

    candidate="$work/photon_hits_run000.root"
    if entries="$(validate_root "$candidate" "$x")"; then
        mv "$candidate" "$final.tmp"
        mv "$final.tmp" "$final"
        rm -rf "$work"
        echo "DONE: x=${x} mm entries=${entries}"
    else
        echo "ERROR: invalid ROOT at x=${x} mm; retained in $work" >&2
        exit 1
    fi
done

echo "=== EJ-230 scan complete. Results: $output_dir ==="
