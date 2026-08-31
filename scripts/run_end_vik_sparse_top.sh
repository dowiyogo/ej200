#!/bin/bash
# run_end_vik_sparse_top.sh — END+Vikuiti + N sparse top SiPMs (EndSparseTop mode)
# Scintillators: EJ-230 (OPSC-106) and EJ-204 (OPSC-101)
# N_sparse: 4, 8, 14, 20 uniformly-spaced top SiPMs
# 13 positions: -600 to +600 mm (step 100 mm), 5000 events/position
#
# Usage:
#   ./run_end_vik_sparse_top.sh            # all configs
#   SCINTS="ej230" ./run_end_vik_sparse_top.sh  # only EJ-230
#   NTOP="4 8"     ./run_end_vik_sparse_top.sh  # only N=4 and N=8

set -euo pipefail

NEVENTS=${NEVENTS:-5000}
BUILD=/home/rrios/ej200/build_end_vikuiti
OUTDIR=${OUTDIR:-/home/rrios/ej200/results/scan_end_vik_sparse_top}
SCINTS=${SCINTS:-"ej230 ej204"}
NTOP_LIST=${NTOP:-"4 8 14 20"}
POSITIONS=(-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600)

mkdir -p "$OUTDIR"

declare -A SCINT_CODE=(
    [ej200]="OPSC-100"
    [ej204]="OPSC-101"
    [ej230]="OPSC-106"
)

run_one() {
    local scint="$1" n_top="$2" xpos="$3"
    local code="${SCINT_CODE[$scint]}"
    local tag="${scint}_ntop${n_top}"
    local outfile="${OUTDIR}/photon_hits_${tag}_x${xpos}.root"

    [[ -f "$outfile" && $(stat -c%s "$outfile") -gt 10000 ]] && { echo "  SKIP $outfile"; return; }

    local mac="${BUILD}/_sparse_${tag}_x${xpos}.mac"
    cat > "$mac" <<MAC
/control/verbose 0
/run/verbose     0
/event/verbose   0
/tracking/verbose 0
/det/readout EndSparseTop
/det/nSparseTopSiPMs ${n_top}
/det/scintillator ${code}
/run/initialize
/gun/particle  mu-
/gun/energy    1 GeV
/muon/angle    0
/muon/gunX     ${xpos} mm
/run/beamOn ${NEVENTS}
MAC

    echo -n "  RUN ${tag} x=${xpos} mm ... "
    local t0=$SECONDS
    (cd "$BUILD" && ./ej200_bar_sim "$(basename "$mac")" > "_sparse_${tag}_x${xpos}.log" 2>&1)
    local dt=$(( SECONDS - t0 ))
    rm -f "$mac"

    local tmp="${BUILD}/photon_hits_run000.root"
    if [[ -f "$tmp" && $(stat -c%s "$tmp") -gt 1000 ]]; then
        mv "$tmp" "$outfile"
        echo "OK (${dt}s)"
        rm -f "${BUILD}/_sparse_${tag}_x${xpos}.log"
    else
        echo "FAILED (${dt}s)"
        rm -f "$tmp" "${BUILD}/_sparse_${tag}_x${xpos}.log"
        return 1
    fi
}

# Count total runs
n_total=0
for s in $SCINTS; do for n in $NTOP_LIST; do for x in "${POSITIONS[@]}"; do n_total=$((n_total+1)); done; done; done

echo "========================================================"
echo " END+Vikuiti+SparseTop scan — $(date '+%Y-%m-%d %H:%M')"
echo " Scintillators : $SCINTS"
echo " N sparse top  : $NTOP_LIST"
echo " Positions     : ${POSITIONS[*]} mm"
echo " Events/pt     : $NEVENTS"
echo " Total runs    : $n_total"
echo " Est. time     : ~$(( n_total * 360 / 3600 )) h"
echo " Output dir    : $OUTDIR"
echo "========================================================"

n_done=0
for scint in $SCINTS; do
    for n_top in $NTOP_LIST; do
        echo "── ${scint^^} N_top=${n_top} ──"
        for x in "${POSITIONS[@]}"; do
            run_one "$scint" "$n_top" "$x" || true
            n_done=$((n_done+1))
            echo "    progress: ${n_done}/${n_total}"
        done
    done
done

echo "========================================================"
echo " Scan complete — $(date '+%Y-%m-%d %H:%M')"
echo " ROOT files: $(ls "$OUTDIR"/*.root 2>/dev/null | wc -l) / ${n_total}"
echo "========================================================"
