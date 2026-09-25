#!/bin/bash
# ============================================================
# run_end_vikuiti_scan.sh
# Position scan — End-only readout, Vikuiti 3M air-gap
# 3 scintillators: EJ-200, EJ-204, EJ-230
# ============================================================

set -euo pipefail

NEVENTS=${NEVENTS:-5000}
BUILD=/home/rrios/ej200/build_end_vikuiti
OUTDIR=/home/rrios/ej200/results/scan_end_vikuiti

POSITIONS=(-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600)

mkdir -p "$OUTDIR"

run_one() {
    local scint="$1"   # OPSC-100 | OPSC-101 | OPSC-106
    local tag="$2"     # ej200 | ej204 | ej230
    local xpos="$3"
    local nevents="$4"

    local outfile="${OUTDIR}/photon_hits_${tag}_x${xpos}.root"
    if [[ -f "$outfile" && $(stat -c%s "$outfile") -gt 10000 ]]; then
        echo "  SKIP  $outfile"
        return 0
    fi

    local macfile="${BUILD}/_scan_${tag}_x${xpos}.mac"
    cat > "$macfile" <<MACEOF
/control/verbose 0
/run/verbose     0
/event/verbose   0
/tracking/verbose 0
/det/readout End
/det/scintillator ${scint}
/run/initialize
/sipm/jitterSigma 0.000 ns
/gun/particle  mu-
/gun/energy    1 GeV
/muon/angle    0
/muon/gunX     ${xpos} mm
/run/beamOn ${nevents}
MACEOF

    echo -n "  RUN   ${tag}  x=${xpos} mm  ... "
    local t0=$SECONDS
    (cd "$BUILD" && ./ej200_bar_sim "$(basename "$macfile")" \
        > "${BUILD}/_scan_${tag}_x${xpos}.log" 2>&1)
    local dt=$(( SECONDS - t0 ))
    rm -f "$macfile"

    local tmproot="${BUILD}/photon_hits_run000.root"
    if [[ -f "$tmproot" && $(stat -c%s "$tmproot") -gt 10000 ]]; then
        mv "$tmproot" "$outfile"
        echo "OK (${dt}s)  →  $(basename "$outfile")"
        rm -f "${BUILD}/_scan_${tag}_x${xpos}.log"
    else
        echo "FAILED (${dt}s)"
        cat "${BUILD}/_scan_${tag}_x${xpos}.log" | grep -E "Exception|error" | head -3 || true
        rm -f "$tmproot" "${BUILD}/_scan_${tag}_x${xpos}.log"
        return 1
    fi
}

echo "================================================================"
echo " End-only + Vikuiti scan — $(date '+%Y-%m-%d %H:%M')"
echo " Positions : ${POSITIONS[*]} mm"
echo " Events/pt : $NEVENTS"
echo " Output    : $OUTDIR"
echo "================================================================"

echo ""
echo "── EJ-200  (OPSC-100) ──────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one OPSC-100 ej200 "$x" "$NEVENTS"
done

echo ""
echo "── EJ-204  (OPSC-101) ──────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one OPSC-101 ej204 "$x" "$NEVENTS"
done

echo ""
echo "── EJ-230  (OPSC-106) ──────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one OPSC-106 ej230 "$x" "$NEVENTS"
done

echo ""
echo "================================================================"
echo " Scan complete — $(date '+%Y-%m-%d %H:%M')"
nfiles=$(ls "$OUTDIR"/*.root 2>/dev/null | wc -l)
echo " ROOT files: $nfiles / 39"
echo "================================================================"
