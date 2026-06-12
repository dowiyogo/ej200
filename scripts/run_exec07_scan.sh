#!/usr/bin/env bash
# EXEC_07 phase 2 campaign. Do not run without explicit phase-2 approval.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="${BUILD_DIR:-$repo/build-exec07}"
sim="$build_dir/ej200_bar_sim"
output_dir="${OUTPUT_DIR:-$repo/results/exec07_endtop_2000}"
events="${EVENTS_PER_POINT:-2000}"
workers="${WORKERS:-16}"
positions=(-690 -670 -650 -600 -550 -500 -450 -400 -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350 400 450 500 550 600 650 670 690)

mkdir -p "$output_dir"

validate_root() {
    python3 - "$1" "$2" "$events" <<'PY'
import sys
import numpy as np
import uproot

path, expected_x, expected_events = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
with uproot.open(path) as root_file:
    tree = root_file["sipm_hits"]
    event_ids = set()
    positions = set()
    for chunk in tree.iterate(["event_id", "gun_x_mm"], step_size="200 MB", library="np"):
        event_ids.update(np.unique(chunk["event_id"]).tolist())
        positions.update(np.unique(chunk["gun_x_mm"]).tolist())
if len(event_ids) != expected_events or positions != {expected_x}:
    raise SystemExit(1)
print(tree.num_entries)
PY
}

for index in "${!positions[@]}"; do
    x="${positions[$index]}"
    final="$output_dir/photon_hits_x${x}mm.root"
    work="$output_dir/.work_x${x}mm"
    log="$output_dir/run_x${x}mm.log"

    if [[ -f "$final" ]] && entries="$(validate_root "$final" "$x" 2>/dev/null)"; then
        echo "DONE: x=${x} mm entries=${entries}"
        continue
    fi

    rm -rf "$work"
    mkdir -p "$work"
    ln -s "$build_dir/sslg4" "$work/sslg4"
    seed1=$((710000 + 1009 * index))
    seed2=$((910000 + 2017 * index))

    cat > "$work/run.mac" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads $workers
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

    echo "START: x=${x} mm events=${events} workers=${workers}"
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
        echo "ERROR: invalid output at x=${x} mm; retained in $work" >&2
        exit 1
    fi
done

echo "EXEC_07 scan complete in $output_dir"
