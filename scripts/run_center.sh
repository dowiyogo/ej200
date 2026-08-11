#!/usr/bin/env bash
# EJ-230 single-point validation run at x=0 mm — parametrized by --threads N.
# Wrappers: run_center_t0minidaq_24t.sh / run_center_msi_16t.sh
set -euo pipefail

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
log_dir="$RESULTS_DIR/logs"
mkdir -p "$log_dir"

events="${EVENTS_CENTER:-200}"
x=0

echo "=== EJ-230 center validation | x=${x} mm threads=${THREADS} events=${events} ==="

work="$RESULTS_DIR/data/.work_center_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$work"
ln -sf "$build_dir/sslg4" "$work/sslg4"

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
/random/setSeeds 720001 920001
/muon/gunX $x mm
/run/beamOn $events
EOF

log="$log_dir/center_x${x}mm_$(date +%Y%m%d_%H%M%S).log"
(
    cd "$work"
    "$sim" -m run.mac 2>&1 | tee "$log.tmp"
)
mv "$log.tmp" "$log"

echo "=== Center run done. Log: $log ==="
echo "    ROOT: $work/photon_hits_run000.root"
