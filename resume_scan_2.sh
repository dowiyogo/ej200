#!/usr/bin/env bash
set -u

repo="$(cd "$(dirname "$0")" && pwd)"
output_root="$repo/results/scan_end_wrapped_2026-06-09/resume_missing_2"
sim="$repo/build/ej200_bar_sim"
mkdir -p "$output_root"

# Positions and seeds copied from the original documented 31-position macro.
positions=(350 400 450 500 550 600 650 670 690)
seed1=(332199 333208 334217 335226 336235 337244 338253 339262 340271)
seed2=(732289 733302 734315 735328 736341 737354 738367 739380 740393)

validate_root() {
    python3 - "$1" "$2" <<'PY'
import sys
import numpy as np
import uproot

path = sys.argv[1]
expected_x = float(sys.argv[2])
with uproot.open(path) as root_file:
    tree = root_file["sipm_hits"]
    event_ids = set()
    positions = set()
    for chunk in tree.iterate(["event_id", "gun_x_mm"], step_size="200 MB", library="np"):
        event_ids.update(np.unique(chunk["event_id"]).tolist())
        positions.update(np.unique(chunk["gun_x_mm"]).tolist())
if len(event_ids) != 10000 or positions != {expected_x}:
    raise SystemExit(1)
print(tree.num_entries)
PY
}

for i in "${!positions[@]}"; do
    x="${positions[$i]}"
    run_dir="$output_root/x_${x}mm"
    root_file="$run_dir/photon_hits_run000.root"
    mkdir -p "$run_dir"

    if [[ -f "$root_file" ]] && entries="$(validate_root "$root_file" "$x" 2>/dev/null)"; then
        echo "DONE: ${x}mm  entradas=${entries} (ya existente)"
        continue
    fi

    rm -f "$root_file"
    cat > "$run_dir/run.mac" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads 16
/det/readout End
/det/scintillatorMaterial EJ204
/run/initialize
/sipm/jitterSigma 0 ns
/gun/particle mu-
/gun/energy 1 GeV
/muon/angle 0
/random/setSeeds ${seed1[$i]} ${seed2[$i]}
/muon/gunX ${x} mm
/run/beamOn 10000
EOF

    echo "START: ${x}mm seeds=${seed1[$i]},${seed2[$i]}"
    (
        cd "$run_dir"
        "$sim" -m run.mac > run.log 2>&1
    )
    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "ERROR: ${x}mm"
        continue
    fi
    if entries="$(validate_root "$root_file" "$x" 2>/dev/null)"; then
        echo "DONE: ${x}mm  entradas=${entries}"
    else
        echo "ERROR: ${x}mm"
    fi
done
