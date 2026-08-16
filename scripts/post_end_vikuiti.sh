#!/bin/bash
# ============================================================
# post_end_vikuiti.sh
# Run after run_end_vikuiti_scan.sh completes:
#   1. Check all 39 ROOT files are present
#   2. Run bar_end_vikuiti_scan.cxx analysis
#   3. Print compilation command for WSL2
# ============================================================

set -euo pipefail

REPO=/home/rrios/ej200
OUTDIR=$REPO/results/scan_end_vikuiti
BUILDDIR=$REPO/build_end_vikuiti
ANALYSIS=$REPO/analysis/bar_end_vikuiti_scan.cxx

echo "================================================================"
echo " Post-scan analysis — $(date '+%Y-%m-%d %H:%M')"
echo "================================================================"

# ── 1. Check ROOT files ─────────────────────────────────────────────
echo ""
echo "── Checking ROOT files ──────────────────────────────────────────"
TAGS=(ej200 ej204 ej230)
POSITIONS=(-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600)
MISSING=0
for tag in "${TAGS[@]}"; do
    for x in "${POSITIONS[@]}"; do
        f="$OUTDIR/photon_hits_${tag}_x${x}.root"
        if [[ -f "$f" && $(stat -c%s "$f") -gt 10000 ]]; then
            echo "  OK   $f"
        else
            echo "  MISS $f"
            MISSING=$((MISSING+1))
        fi
    done
done

if [[ $MISSING -gt 0 ]]; then
    echo ""
    echo "  WARNING: $MISSING files missing or too small. Analysis may be incomplete."
else
    echo ""
    echo "  All 39 ROOT files present. ✓"
fi

# ── 2. Run analysis ─────────────────────────────────────────────────
echo ""
echo "── Running analysis macro ──────────────────────────────────────"
cd "$REPO"
root -l -b -q "$ANALYSIS" 2>&1 | tee /tmp/end_vikuiti_analysis.log

echo ""
echo "── Output plots ────────────────────────────────────────────────"
ls -lh "$BUILDDIR"/end_vik_*.{png,pdf} 2>/dev/null || echo "  No plots found."

# ── 3. Print WSL2 compilation command ────────────────────────────────
echo ""
echo "================================================================"
echo " To compile the Beamer on WSL2:"
echo ""
echo "   ssh -p 2222 reriosto@128.141.225.8"
echo "   cd /home/reriosto/SHiP/ej200"
echo "   git pull origin feat/bar-end-vikuiti"
echo "   # copy PDFs from lxplus or local (see below)"
echo "   cd presentations/end_vikuiti_2026-08-16"
echo "   pdflatex talk.tex && pdflatex talk.tex"
echo ""
echo " Then scp the PDF back:"
echo "   scp -P 2222 reriosto@128.141.225.8:/home/reriosto/SHiP/ej200/presentations/end_vikuiti_2026-08-16/talk.pdf \\"
echo "       $REPO/presentations/end_vikuiti_2026-08-16/"
echo "================================================================"
