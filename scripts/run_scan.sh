#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --threads N --host NAME [--events E] [--positions FILE]" >&2
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$script_dir/.." && pwd)"
threads=""
host_name=""
events=2000
positions_file="$script_dir/endonly_positions.txt"
mylar_reflectivity="${MYLAR_REFLECTIVITY:-0.90}"
mylar_specular_lobe="${MYLAR_SPECULAR_LOBE:-1.0}"
mylar_sigma_alpha_deg="${MYLAR_SIGMA_ALPHA_DEG:-0.1}"

while (($#)); do
    case "$1" in
        --threads)
            threads="${2:-}"
            shift 2
            ;;
        --host)
            host_name="${2:-}"
            shift 2
            ;;
        --events)
            events="${2:-}"
            shift 2
            ;;
        --positions)
            positions_file="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 2
            ;;
    esac
done

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer" >&2; exit 2; }
[[ "$events" =~ ^[1-9][0-9]*$ ]] || { echo "--events must be a positive integer" >&2; exit 2; }
[[ -n "$host_name" ]] || { echo "--host is required" >&2; exit 2; }
[[ -f "$positions_file" ]] || { echo "positions file not found: $positions_file" >&2; exit 2; }

host_tag="$(printf '%s' "$host_name" | tr -c '[:alnum:]_.-' '_')"
build_dir="${BUILD_DIR:-$repo/build}"
output_root="${OUTPUT_ROOT:-$repo/output}"
sim="$build_dir/ej200_bar_sim"

if [[ ! -f "$build_dir/CMakeCache.txt" ]]; then
    cmake -S "$repo" -B "$build_dir" -DCMAKE_BUILD_TYPE=RelWithDebInfo
fi
cmake --build "$build_dir" --parallel "$threads"

timestamp="$(date +%Y%m%d_%H%M%S)"
output_dir="$output_root/endonly_mylar_${host_tag}_${timestamp}"
suffix=1
while [[ -e "$output_dir" ]]; do
    printf -v output_dir '%s/endonly_mylar_%s_%s_%02d' "$output_root" "$host_tag" "$timestamp" "$suffix"
    ((suffix += 1))
done
mkdir -p "$output_dir"

positions=()
declare -A seen_positions=()
while IFS= read -r raw || [[ -n "$raw" ]]; do
    position="${raw%%#*}"
    position="${position//[[:space:]]/}"
    [[ -z "$position" ]] && continue
    [[ "$position" =~ ^-?[0-9]+([.][0-9]+)?$ ]] || {
        echo "invalid position '$position' in $positions_file" >&2
        exit 2
    }
    [[ -z "${seen_positions[$position]:-}" ]] || {
        echo "duplicate position '$position' in $positions_file" >&2
        exit 2
    }
    seen_positions["$position"]=1
    positions+=("$position")
done < "$positions_file"
(( ${#positions[@]} > 0 )) || { echo "positions file is empty: $positions_file" >&2; exit 2; }

commit="$(git -C "$repo" rev-parse HEAD)"
start_iso="$(date --iso-8601=seconds)"
{
    echo "commit=$commit"
    echo "git_dirty=$(git -C "$repo" status --porcelain | wc -l)"
    echo "host=$host_name"
    echo "threads=$threads"
    echo "events_per_position=$events"
    echo "positions_file=$(realpath "$positions_file")"
    echo "positions=${positions[*]}"
    echo "top_surface=mylar"
    echo "mylar_reflectivity=$mylar_reflectivity"
    echo "mylar_specular_lobe=$mylar_specular_lobe"
    echo "mylar_sigma_alpha_deg=$mylar_sigma_alpha_deg"
    echo "start=$start_iso"
} > "$output_dir/run_metadata.txt"

echo "Output: $output_dir"
for index in "${!positions[@]}"; do
    x="${positions[$index]}"
    x_tag="${x//./p}"
    work="$output_dir/work_x${x_tag}mm"
    log="$output_dir/run_x${x_tag}mm.log"
    final="$output_dir/photon_hits_x${x_tag}mm.root"
    mkdir -p "$work"
    ln -s "$build_dir/sslg4" "$work/sslg4"

    seed1=$((730000 + 1009 * index))
    seed2=$((930000 + 2017 * index))
    cat > "$work/run.mac" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout EndTop
/ship/geom/topSurface mylar
/ship/geom/mylar/reflectivity $mylar_reflectivity
/ship/geom/mylar/specularLobe $mylar_specular_lobe
/ship/geom/mylar/sigmaAlpha $mylar_sigma_alpha_deg deg
/det/scintillator OPSC-101
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

    echo "START x=$x mm events=$events threads=$threads"
    (
        cd "$work"
        "$sim" -t "$threads" -m run.mac
    ) > "$log" 2>&1

    candidate="$work/photon_hits_run000.root"
    [[ -s "$candidate" ]] || { echo "missing ROOT output for x=$x; retained in $work" >&2; exit 1; }
    root -l -b -q -e \
        "TFile file(\"$candidate\", \"READ\"); gSystem->Exit(file.IsZombie() ? 1 : 0);" \
        >/dev/null
    [[ ! -e "$final" ]] || { echo "refusing to overwrite $final" >&2; exit 1; }
    mv "$candidate" "$final"
    echo "DONE x=$x mm root=$final"
done

python3 "$repo/analysis/endonly_sum4.py" --input-dir "$output_dir" --events "$events"

{
    echo "end=$(date --iso-8601=seconds)"
    echo "status=complete"
} >> "$output_dir/run_metadata.txt"
echo "End-only Mylar scan complete: $output_dir"
