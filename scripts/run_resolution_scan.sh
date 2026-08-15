#!/bin/bash
# ============================================================
# run_resolution_scan.sh
# Full position scan for timing resolution and group velocity
# 4 configurations: EJ-204/EJ-230 × TIR-only/Vikuiti 3M
#
# NOTE: Must run from each build dir so SSLG4 macro path resolves.
# ============================================================

set -euo pipefail

NEVENTS=${NEVENTS:-5000}
BUILD_TIR=/home/rrios/ej200/build_t0minidaq
BUILD_VIK=/home/rrios/ej200/build_vikuiti
OUTDIR=/home/rrios/ej200/results/scan_resolution

# Scan positions in mm (±700mm bar; SiPMs at ±700mm)
POSITIONS=(-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600)

mkdir -p "$OUTDIR"

# run_one <build_dir> <scint> <tag> <xpos> <nevents>
run_one() {
    local bdir="$1"
    local scint="$2"    # OPSC-101 (EJ-204) or OPSC-106 (EJ-230)
    local tag="$3"      # ej204_tir | ej230_tir | ej204_vik | ej230_vik
    local xpos="$4"
    local nevents="$5"

    local outfile="${OUTDIR}/photon_hits_${tag}_x${xpos}.root"
    if [[ -f "$outfile" && $(stat -c%s "$outfile") -gt 10000 ]]; then
        echo "  SKIP  $outfile"
        return 0
    fi

    # Write macro to a temp file in the build dir (so relative sslg4/ paths resolve)
    local macfile="${bdir}/_scan_tmp_${tag}_x${xpos}.mac"
    cat > "$macfile" <<MACEOF
/control/verbose 0
/run/verbose     0
/event/verbose   0
/tracking/verbose 0
/det/readout EndTop
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

    # Run from build dir so sslg4/macros/ relative paths resolve
    (
        cd "$bdir"
        ./ej200_bar_sim "$(basename "$macfile")" > "${bdir}/_scan_tmp_${tag}_x${xpos}.log" 2>&1
    )
    local dt=$(( SECONDS - t0 ))
    rm -f "$macfile"

    local tmproot="${bdir}/photon_hits_run000.root"
    if [[ -f "$tmproot" && $(stat -c%s "$tmproot") -gt 10000 ]]; then
        mv "$tmproot" "$outfile"
        echo "OK (${dt}s)  →  $(basename "$outfile")"
        rm -f "${bdir}/_scan_tmp_${tag}_x${xpos}.log"
    else
        echo "FAILED (${dt}s)"
        cat "${bdir}/_scan_tmp_${tag}_x${xpos}.log" | grep -E "Exception|error|Error" | head -5 || true
        rm -f "$tmproot" "${bdir}/_scan_tmp_${tag}_x${xpos}.log"
        return 1
    fi
}

echo "================================================================"
echo " Resolution scan — $(date '+%Y-%m-%d %H:%M')"
echo " Positions : ${POSITIONS[*]} mm"
echo " Events/pt : $NEVENTS"
echo " Output    : $OUTDIR"
echo "================================================================"

echo ""
echo "── EJ-204  TIR-only ────────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one "$BUILD_TIR" OPSC-101 ej204_tir "$x" "$NEVENTS"
done

echo ""
echo "── EJ-230  TIR-only ────────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one "$BUILD_TIR" OPSC-106 ej230_tir "$x" "$NEVENTS"
done

echo ""
echo "── EJ-204  Vikuiti 3M ──────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one "$BUILD_VIK" OPSC-101 ej204_vik "$x" "$NEVENTS"
done

echo ""
echo "── EJ-230  Vikuiti 3M ──────────────────────────────────────────"
for x in "${POSITIONS[@]}"; do
    run_one "$BUILD_VIK" OPSC-106 ej230_vik "$x" "$NEVENTS"
done

echo ""
echo "================================================================"
echo " Scan complete — $(date '+%Y-%m-%d %H:%M')"
echo " Files in: $OUTDIR"
nfiles=$(ls "$OUTDIR"/*.root 2>/dev/null | wc -l)
echo " ROOT files: $nfiles"
echo "================================================================"
