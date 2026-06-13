#!/usr/bin/env bash
# Reproduce the EJ-230 EXEC_08b directed runs and EXEC_09 End-only controls.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="${BUILD_DIR:-$repo/build}"
sim="$build_dir/ej200_bar_sim"
output_root="${OUTPUT_ROOT:-$repo/results_ej230_analysis/sim_missing}"
workers="${WORKERS:-16}"

if [[ ! -x "$sim" ]]; then
    echo "ERROR: simulation binary not found: $sim" >&2
    exit 1
fi

validate_root() {
    python3 - "$1" "$2" "$3" "$4" <<'PY'
import sys
import numpy as np
import uproot

path, expected_x, expected_events, expected_max_channel = (
    sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
)
with uproot.open(path) as root_file:
    tree = root_file["sipm_hits"]
    if tree.num_entries <= 0:
        raise SystemExit(f"{path}: empty sipm_hits tree")
    event_ids = set()
    positions = set()
    channels = set()
    for chunk in tree.iterate(
        ["event_id", "global_id", "gun_x_mm"], step_size="200 MB", library="np"
    ):
        event_ids.update(np.unique(chunk["event_id"]).astype(int).tolist())
        positions.update(np.unique(chunk["gun_x_mm"]).astype(float).tolist())
        channels.update(np.unique(chunk["global_id"]).astype(int).tolist())
if event_ids != set(range(expected_events)):
    raise SystemExit(f"{path}: expected {expected_events} complete events, got {len(event_ids)}")
if positions != {expected_x}:
    raise SystemExit(f"{path}: gun_x_mm={positions}, expected {expected_x}")
if not channels or min(channels) != 0 or max(channels) != expected_max_channel:
    raise SystemExit(f"{path}: channel range {min(channels)}..{max(channels)}, expected 0..{expected_max_channel}")
print(f"entries={int(tree.num_entries)} events={len(event_ids)} channels=0..{max(channels)} gun_x_mm={expected_x:g}")
PY
}

run_one() {
    local family="$1"
    local label="$2"
    local x="$3"
    local events="$4"
    local readout="$5"
    local seed1="$6"
    local seed2="$7"
    local max_channel="$8"
    local family_dir="$output_root/$family"
    local final="$family_dir/photon_hits_${label}.root"
    local work="$family_dir/.work_${label}"
    local log="$family_dir/run_${label}.log"

    mkdir -p "$family_dir"
    if [[ -f "$final" ]] && validate_root "$final" "$x" "$events" "$max_channel" >/dev/null 2>&1; then
        echo "DONE existing: $final"
        validate_root "$final" "$x" "$events" "$max_channel"
        return
    fi
    if [[ -f "$final" ]]; then
        echo "ERROR: preserving invalid existing ROOT instead of overwriting it: $final" >&2
        exit 1
    fi

    if [[ -e "$work" ]]; then
        echo "ERROR: preserving incomplete work directory instead of deleting it: $work" >&2
        exit 1
    fi
    mkdir -p "$work"
    ln -s "$build_dir/sslg4" "$work/sslg4"
    {
        printf '%s\n' \
            "/control/verbose 0" \
            "/run/verbose 0" \
            "/event/verbose 0" \
            "/tracking/verbose 0" \
            "/run/numberOfThreads $workers" \
            "/det/readout $readout" \
            "/det/scintillator OPSC-106" \
            "/sipm/model AFBR-S4N66P024M" \
            "/run/initialize" \
            "/sipm/jitterSigma 0 ns" \
            "/gun/particle mu-" \
            "/gun/energy 1 GeV" \
            "/muon/angle 0" \
            "/random/setSeeds $seed1 $seed2" \
            "/muon/gunX $x mm" \
            "/run/beamOn $events"
    } > "$work/run.mac"

    echo "START: family=$family label=$label x=$x events=$events readout=$readout workers=$workers"
    (
        cd "$work"
        "$sim" -m run.mac 2>&1 | tee "$log.tmp"
    )
    mv "$log.tmp" "$log"
    validate_root "$work/photon_hits_run000.root" "$x" "$events" "$max_channel"
    mv "$work/photon_hits_run000.root" "$final.tmp"
    mv "$final.tmp" "$final"
    rm -rf "$work"
    echo "DONE: $final"
}

# EXEC_08b: exact inherited directed-run positions and seed pairs, 2000 events.
run_one window_dip run_A_x-652mm -652 2000 EndTop 808201 808211 85
run_one window_dip run_B_x-642mm -642 2000 EndTop 808301 808307 85
run_one window_dip run_C1_x-648mm -648 2000 EndTop 808401 808421 85
run_one window_dip run_C2_x-654mm -654 2000 EndTop 808501 808519 85

# EXEC_09: selected positions from the inherited 10k-event End-only campaign.
run_one end_only run007 -400 10000 End 317064 717094 15
run_one end_only run015 0 10000 End 325136 725198 15
run_one end_only run023 400 10000 End 333208 733302 15

echo "EXEC_14B special simulations complete: $output_root"
