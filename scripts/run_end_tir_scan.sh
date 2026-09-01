#!/bin/bash
# run_end_tir_scan.sh — End-only readout, TIR only (no Vikuiti panels)
# EJ-204 and EJ-230, 13 positions

set -euo pipefail

NEVENTS=${NEVENTS:-5000}
BUILD=/home/rrios/ej200/build_t0minidaq
OUTDIR=/home/rrios/ej200/results/scan_end_tir

POSITIONS=(-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600)
mkdir -p "$OUTDIR"

run_one() {
    local scint="$1" tag="$2" xpos="$3"
    local outfile="${OUTDIR}/photon_hits_${tag}_x${xpos}.root"
    [[ -f "$outfile" && $(stat -c%s "$outfile") -gt 5000 ]] && { echo "  SKIP $outfile"; return; }

    local mac="${BUILD}/_tir_${tag}_x${xpos}.mac"
    cat > "$mac" <<MAC
/control/verbose 0
/run/verbose     0
/event/verbose   0
/tracking/verbose 0
/det/readout End
/det/scintillator ${scint}
/run/initialize
/gun/particle  mu-
/gun/energy    1 GeV
/muon/angle    0
/muon/gunX     ${xpos} mm
/run/beamOn ${NEVENTS}
MAC
    echo -n "  RUN ${tag} x=${xpos} mm ... "
    local t0=$SECONDS
    (cd "$BUILD" && ./ej200_bar_sim "$(basename "$mac")" > "_tir_${tag}_x${xpos}.log" 2>&1)
    local dt=$(( SECONDS - t0 ))
    rm -f "$mac"

    local tmp="${BUILD}/photon_hits_run000.root"
    if [[ -f "$tmp" && $(stat -c%s "$tmp") -gt 1000 ]]; then
        mv "$tmp" "$outfile"
        echo "OK (${dt}s)"
        rm -f "${BUILD}/_tir_${tag}_x${xpos}.log"
    else
        echo "FAILED (${dt}s)"
        rm -f "$tmp" "${BUILD}/_tir_${tag}_x${xpos}.log"
        return 1
    fi
}

echo "================================================================"
echo " END+TIR scan — $(date '+%Y-%m-%d %H:%M')"
echo " Positions: ${POSITIONS[*]} mm  |  Events/pt: $NEVENTS"
echo "================================================================"

echo "── EJ-204 (OPSC-101) ──"
for x in "${POSITIONS[@]}"; do run_one OPSC-101 ej204 "$x"; done

echo "── EJ-230 (OPSC-106) ──"
for x in "${POSITIONS[@]}"; do run_one OPSC-106 ej230 "$x"; done

echo "================================================================"
echo " Scan complete — $(date '+%Y-%m-%d %H:%M')"
echo " ROOT files: $(ls "$OUTDIR"/*.root 2>/dev/null | wc -l) / 26"
echo "================================================================"
