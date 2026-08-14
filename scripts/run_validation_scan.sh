#!/usr/bin/env bash
# ============================================================
# run_validation_scan.sh
# Physics validation scan — feat/endtop-sslg4
# 7 positions x 500 events x 24 threads (~5 min total)
#
# Checks:
#   1. Npe at END SiPMs > 10 PE/event (was ~0.3 with old physics)
#   2. v_eff from timing slope ~ 15-19 cm/ns (physical range)
#   3. No crash at any position
#
# Usage:
#   cd /home/rrios/ej204
#   source /home/tdship/opt/geant4-v11.4.0-install/bin/geant4.sh
#   ./scripts/run_validation_scan.sh
# ============================================================
set -euo pipefail

REPO_DIR="/home/rrios/ej204"
BUILD_DIR="${REPO_DIR}/build_t0minidaq"
GEANT4_SETUP="/home/tdship/opt/geant4-v11.4.0-install/bin/geant4.sh"
MACRO_DIR="${REPO_DIR}/macros/validation_scan"
SIM="${BUILD_DIR}/ej200_bar_sim"
SSLG4_DATA="${BUILD_DIR}/sslg4"
EXPECTED_BRANCH="feat/endtop-sslg4"

# ── Load Geant4 ───────────────────────────────────────────────
if ! command -v geant4-config &>/dev/null; then
    [[ -f "$GEANT4_SETUP" ]] && source "$GEANT4_SETUP" || {
        echo "ERROR: geant4.sh not found at $GEANT4_SETUP" >&2; exit 1; }
fi

cd "$REPO_DIR"

# ── Pre-flight ────────────────────────────────────────────────
current_branch=$(git branch --show-current)
[[ "$current_branch" != "$EXPECTED_BRANCH" ]] && {
    echo "ERROR: branch '$current_branch' != '$EXPECTED_BRANCH'" >&2; exit 1; }
[[ ! -x "$SIM" ]] && { echo "ERROR: $SIM not executable" >&2; exit 1; }
[[ ! -d "$SSLG4_DATA" ]] && { echo "ERROR: $SSLG4_DATA missing" >&2; exit 1; }

COMMIT=$(git rev-parse --short HEAD)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="${REPO_DIR}/runs/validation_${TIMESTAMP}"
mkdir -p "${RUN_DIR}/logs"

echo "============================================================"
echo "  Validation Scan — EndTop EJ-204"
echo "  Branch : $current_branch @ $COMMIT"
echo "  Output : $RUN_DIR"
echo "  Macros : 7 positions x 500 events x 24 threads"
echo "============================================================"

PASS=0; FAIL=0
declare -a ROOT_FILES

for mac in $(ls "$MACRO_DIR"/validation_*.mac | sort); do
    base=$(basename "$mac" .mac)
    pos=$(echo "$base" | grep -o 'x[+-]*[0-9]*mm' | head -1)
    pos_dir="${RUN_DIR}/${pos}"
    mkdir -p "$pos_dir"
    ln -sf "$SSLG4_DATA" "${pos_dir}/sslg4"

    echo -n "  Running $pos ... "
    log="${RUN_DIR}/logs/${base}.log"

    (cd "$pos_dir" && "$SIM" -m "$mac" > "$log" 2>&1) && rc=0 || rc=$?

    root_f=$(find "$pos_dir" -maxdepth 1 -name "*.root" 2>/dev/null | head -1)
    if [[ $rc -eq 0 && -n "$root_f" ]]; then
        echo "OK  -> $(basename $root_f)"
        ROOT_FILES+=("$root_f")
        PASS=$((PASS+1))
    else
        echo "FAIL (rc=$rc, root=$([ -n "$root_f" ] && echo yes || echo no))"
        FAIL=$((FAIL+1))
    fi
done

echo ""
echo "============================================================"
echo "  Positions OK: $PASS / 7    FAIL: $FAIL"
echo "  Run dir: $RUN_DIR"
echo "============================================================"

# Write root file list for analysis
printf '%s\n' "${ROOT_FILES[@]}" > "${RUN_DIR}/root_files.txt"
echo "${RUN_DIR}" > "${REPO_DIR}/runs/last_validation_run.txt"

if [[ $FAIL -gt 0 ]]; then
    echo "  Some positions FAILED — check logs in ${RUN_DIR}/logs/"
    exit 1
fi

echo ""
echo "  Run analysis:"
echo "    python3 scripts/analyze_validation.py ${RUN_DIR}"
