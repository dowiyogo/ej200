#!/usr/bin/env bash
# ============================================================
# run_t0minidaq_endtop_scan_5000.sh
# Longitudinal EndTop scan — EJ-200 / feat/endtop-sslg4
# 31 positions x 5000 events x 24 threads
# Server: t0minidaq / AlmaLinux 9
#
# Usage:
#   ./scripts/run_t0minidaq_endtop_scan_5000.sh [OPTIONS]
#
# Options:
#   --dry-run       List what would be done, no simulation
#   --smoke         Run 1 event at x=-690 mm to validate app (no crash)
#   --first-only    Run only the first position (x=-690 mm, 5000 events)
#   --from N        Start from macro index N (1-based, inclusive)
#   --to   N        Stop  at macro index N (1-based, inclusive)
#
# Environment variables:
#   STOP_ON_ERROR=1 (default) — abort on first failure
#   STOP_ON_ERROR=0           — continue and record failures
# ============================================================
set -u
set -o pipefail

# ── Configurable ─────────────────────────────────────────────
STOP_ON_ERROR=${STOP_ON_ERROR:-1}
REPO_DIR="/home/rrios/ej204"
BUILD_DIR="${REPO_DIR}/build_t0minidaq"
GEANT4_SETUP="/home/tdship/opt/geant4-v11.4.0-install/bin/geant4.sh"
EXPECTED_BRANCH="feat/endtop-sslg4"
EVENTS_PER_POINT=5000
THREADS=24
MACRO_DIR="${REPO_DIR}/macros/t0minidaq_scan_5000"
SIM="${BUILD_DIR}/ej200_bar_sim"
# SSLG4 data dir — must be symlinked into each run dir so SSLG4
# can find its material-definition macros (sslg4/macros/oscnt/...).
SSLG4_DATA="${BUILD_DIR}/sslg4"

# ── Parse arguments ──────────────────────────────────────────
DRY_RUN=0
SMOKE=0
FROM_IDX=1
TO_IDX=31

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)    DRY_RUN=1 ;;
        --smoke)      SMOKE=1 ;;
        --first-only) FROM_IDX=1 ; TO_IDX=1 ;;
        --from)       FROM_IDX="$2" ; shift ;;
        --to)         TO_IDX="$2"   ; shift ;;
        *)
            echo "ERROR: argumento desconocido: $1" >&2
            echo "Uso: $0 [--dry-run] [--smoke] [--first-only] [--from N] [--to N]" >&2
            exit 1 ;;
    esac
    shift
done

# ── Load Geant4 if not already loaded ────────────────────────
if ! command -v geant4-config &>/dev/null; then
    if [[ -f "$GEANT4_SETUP" ]]; then
        # shellcheck source=/dev/null
        source "$GEANT4_SETUP"
    else
        echo "ERROR: no se encontró geant4.sh en $GEANT4_SETUP" >&2
        exit 1
    fi
fi

# ── Enter repo ───────────────────────────────────────────────
cd "$REPO_DIR" || { echo "ERROR: no se puede entrar a $REPO_DIR" >&2; exit 1; }

# ── Pre-flight checks ────────────────────────────────────────
echo "=== PRE-FLIGHT CHECKS ==="

if [[ ! -x "$SIM" ]]; then
    echo "ERROR: ejecutable no encontrado o no ejecutable: $SIM" >&2; exit 1
fi
echo "  Ejecutable:    $SIM  ✓"

if [[ ! -d "$SSLG4_DATA" ]]; then
    echo "ERROR: directorio SSLG4 no encontrado: $SSLG4_DATA" >&2; exit 1
fi
echo "  SSLG4 data:    $SSLG4_DATA  ✓"

current_branch=$(git branch --show-current 2>/dev/null || echo "UNKNOWN")
if [[ "$current_branch" != "$EXPECTED_BRANCH" ]]; then
    echo "ERROR: branch actual '$current_branch' != esperada '$EXPECTED_BRANCH'" >&2; exit 1
fi
echo "  Branch:        $current_branch  ✓"

if [[ ! -d "$MACRO_DIR" ]]; then
    echo "ERROR: directorio de macros no encontrado: $MACRO_DIR" >&2; exit 1
fi

mapfile -t ALL_MACROS < <(ls "$MACRO_DIR"/endtop_5000_[0-9]*.mac 2>/dev/null | sort)
TOTAL_MACROS=${#ALL_MACROS[@]}
if [[ $TOTAL_MACROS -eq 0 ]]; then
    echo "ERROR: no se encontraron macros en $MACRO_DIR" >&2; exit 1
fi
echo "  Macros totales: $TOTAL_MACROS"

if [[ $FROM_IDX -lt 1 || $TO_IDX -gt $TOTAL_MACROS || $FROM_IDX -gt $TO_IDX ]]; then
    echo "ERROR: rango inválido --from $FROM_IDX --to $TO_IDX (total: $TOTAL_MACROS)" >&2
    exit 1
fi

# ── Gather environment info ───────────────────────────────────
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "UNKNOWN")
G4VERSION=$(geant4-config --version 2>/dev/null || echo "UNKNOWN")
G4PREFIX=$(geant4-config --prefix 2>/dev/null || echo "UNKNOWN")
ROOT_VERSION=$(root-config --version 2>/dev/null || echo "UNKNOWN")
HOSTNAME_VAL=$(hostname)
START_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
START_DATETIME=$(date '+%Y-%m-%d %H:%M:%S')

# ── Run directory ─────────────────────────────────────────────
if [[ $SMOKE -eq 1 ]]; then
    MODE_TAG="smoke"
elif [[ $DRY_RUN -eq 1 ]]; then
    MODE_TAG="dryrun"
else
    MODE_TAG="scan"
fi

RUN_DIR="${REPO_DIR}/runs/t0minidaq_endtop_${MODE_TAG}_${START_TIMESTAMP}"

if [[ $DRY_RUN -eq 0 ]]; then
    mkdir -p "${RUN_DIR}/logs"
    mkdir -p "${RUN_DIR}/macros_used"
    mkdir -p "${RUN_DIR}/outputs"
    MASTER_LOG="${RUN_DIR}/logs/master_scan.log"
    SUMMARY_TSV="${RUN_DIR}/summary.tsv"
    echo -e "index\tposition_mm\tmacro\tlog\tstatus\tstart_time\tend_time\tduration_seconds" \
        > "$SUMMARY_TSV"
else
    MASTER_LOG="/dev/null"
    SUMMARY_TSV="/dev/null"
fi

# ── Print header ──────────────────────────────────────────────
print_header() {
    local mode_label="SCAN COMPLETO"
    [[ $DRY_RUN -eq 1 ]] && mode_label="DRY-RUN"
    [[ $SMOKE   -eq 1 ]] && mode_label="SMOKE (1 evento)"
    cat <<HEADER
============================================================
  EndTop Scan 5000 — t0minidaq   [${mode_label}]
============================================================
  Fecha/hora inicio : ${START_DATETIME}
  Hostname          : ${HOSTNAME_VAL}
  Repo              : ${REPO_DIR}
  Branch            : ${current_branch}
  Commit            : ${COMMIT}
  Geant4            : ${G4VERSION}  (${G4PREFIX})
  ROOT              : ${ROOT_VERSION}
  Ejecutable        : ${SIM}
  SSLG4 data        : ${SSLG4_DATA}
  Macros totales    : ${TOTAL_MACROS}
  Rango a correr    : [${FROM_IDX}..${TO_IDX}]
  Eventos/posición  : ${EVENTS_PER_POINT}
  Threads           : ${THREADS}
  Directorio salida : ${RUN_DIR}
  STOP_ON_ERROR     : ${STOP_ON_ERROR}
============================================================
HEADER
}

# ── Core: run one macro file ──────────────────────────────────
# Usage: run_one_macro <label> <pos_mm> <macro_path> <pos_dir> <log_file> [status_out_file]
# Writes "OK" or "FAIL..." to status_out_file (required).
# All sim output goes to: terminal (stdout) + log_file + MASTER_LOG.
run_one_macro() {
    local label="$1"
    local pos_mm="$2"
    local macro_path="$3"
    local pos_dir="$4"
    local log_file="$5"
    local status_file="$6"

    # ── KEY: symlink sslg4 data so SSLG4 can load material macros ──
    # SSLG4 looks for sslg4/macros/oscnt/<model>.mac relative to CWD.
    ln -sf "$SSLG4_DATA" "${pos_dir}/sslg4"

    local t_start t_start_epoch
    t_start=$(date '+%Y-%m-%d %H:%M:%S')
    t_start_epoch=$(date +%s)

    echo "" | tee -a "$MASTER_LOG"
    echo "START ${label}  ${t_start}" | tee -a "$MASTER_LOG"

    # Run sim from pos_dir; pipe output to terminal + log + master log
    local SIM_RC=0
    (
        cd "$pos_dir"
        stdbuf -oL -eL "$SIM" -m "$macro_path" 2>&1
    ) | tee -a "$log_file" | tee -a "$MASTER_LOG" || true
    # Capture sim exit code via PIPESTATUS of the subshell
    SIM_RC=${PIPESTATUS[0]}

    local t_end t_end_epoch duration
    t_end=$(date '+%Y-%m-%d %H:%M:%S')
    t_end_epoch=$(date +%s)
    duration=$(( t_end_epoch - t_start_epoch ))

    # Determine status
    local STATUS
    if [[ $SIM_RC -ne 0 ]]; then
        STATUS="FAIL(exit=${SIM_RC})"
    else
        local root_file
        root_file=$(find "$pos_dir" -maxdepth 1 -name "*.root" 2>/dev/null | head -1)
        if [[ -n "$root_file" ]]; then
            STATUS="OK"
        else
            STATUS="FAIL(no-root)"
        fi
    fi

    # Print result line
    if [[ "$STATUS" == "OK" ]]; then
        echo "OK    ${label}  ${t_end}  duration=${duration}s" | tee -a "$MASTER_LOG"
    else
        echo "FAIL  ${label}  ${STATUS}  ${t_end}  duration=${duration}s" \
            | tee -a "$MASTER_LOG"
    fi

    # Write TSV row
    local macro_basename
    macro_basename=$(basename "$macro_path")
    echo -e "$(basename "$log_file" .log | sed 's/run_//')\t${pos_mm}\t${macro_basename}\t$(basename "$log_file")\t${STATUS}\t${t_start}\t${t_end}\t${duration}" \
        >> "$SUMMARY_TSV"

    # Return status via file (avoids $() stdout-capture problems)
    echo "$STATUS" > "$status_file"
}

# ── DRY-RUN mode ─────────────────────────────────────────────
if [[ $DRY_RUN -eq 1 ]]; then
    print_header
    echo ""
    echo "=== DRY-RUN: macros que se ejecutarían ==="
    for (( i=0; i<TOTAL_MACROS; i++ )); do
        macro_idx=$(( i + 1 ))
        [[ $macro_idx -lt $FROM_IDX || $macro_idx -gt $TO_IDX ]] && continue
        macro_path="${ALL_MACROS[$i]}"
        macro_basename=$(basename "$macro_path")
        pos_raw="${macro_basename##*_x}"
        pos_mm="${pos_raw%mm.mac}"
        pos_mm="${pos_mm#+}"
        label="[$(printf '%03d' "$macro_idx")/${TOTAL_MACROS}] x=${pos_mm}mm"
        echo "  ${label}  ->  ${macro_basename}"
    done
    echo ""
    echo "=== DRY-RUN completado ==="
    echo "  Macros que se ejecutarían: $(( TO_IDX - FROM_IDX + 1 ))"
    echo "  Eventos por posición:      ${EVENTS_PER_POINT}"
    echo "  Threads:                   ${THREADS}"
    echo ""
    echo "  Smoke test (1 evento):     ./scripts/run_t0minidaq_endtop_scan_5000.sh --smoke"
    echo "  Primera posición (5000e):  ./scripts/run_t0minidaq_endtop_scan_5000.sh --first-only"
    echo "  Scan completo:             ./scripts/run_t0minidaq_endtop_scan_5000.sh"
    exit 0
fi

# ── SMOKE mode ────────────────────────────────────────────────
if [[ $SMOKE -eq 1 ]]; then
    print_header | tee "$MASTER_LOG"
    echo ""
    echo "=== SMOKE: 1 evento en x=-690 mm (macro temporal, no modifica las finales) ==="

    smoke_macro_src="${ALL_MACROS[0]}"
    smoke_dir="${RUN_DIR}/outputs/x-690mm"
    mkdir -p "$smoke_dir"

    # Create a 1-event temporary macro (does not touch the 5000-event originals)
    smoke_mac="${smoke_dir}/smoke_1evt.mac"
    sed 's|/run/beamOn [0-9][0-9]*|/run/beamOn 1|' "$smoke_macro_src" > "$smoke_mac"
    cp "$smoke_mac" "${RUN_DIR}/macros_used/smoke_1evt.mac"

    STAT_FILE=$(mktemp)
    run_one_macro \
        "[SMOKE] x=-690mm" \
        "-690" \
        "$smoke_mac" \
        "$smoke_dir" \
        "${RUN_DIR}/logs/run_smoke_x-690mm.log" \
        "$STAT_FILE"

    smoke_result=$(cat "$STAT_FILE")
    rm -f "$STAT_FILE"

    echo ""
    echo "============================================================"
    if [[ "$smoke_result" == "OK" ]]; then
        echo "  SMOKE PASSED — la app corrió 1 evento sin abortar."
        echo "  ROOT output: $(find "$smoke_dir" -name '*.root' | head -1)"
        echo ""
        echo "  Listo para el scan completo:"
        echo "    tmux new -s ej204_scan_5000"
        echo "    cd /home/rrios/ej204"
        echo "    source /home/tdship/opt/geant4-v11.4.0-install/bin/geant4.sh"
        echo "    ./scripts/run_t0minidaq_endtop_scan_5000.sh"
    else
        echo "  SMOKE FAILED — status: ${smoke_result}"
        echo "  Log: ${RUN_DIR}/logs/run_smoke_x-690mm.log"
    fi
    echo "============================================================"
    [[ "$smoke_result" == "OK" ]] && exit 0 || exit 1
fi

# ── SCAN mode ─────────────────────────────────────────────────
print_header | tee "$MASTER_LOG"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for (( i=0; i<TOTAL_MACROS; i++ )); do
    macro_idx=$(( i + 1 ))

    if [[ $macro_idx -lt $FROM_IDX || $macro_idx -gt $TO_IDX ]]; then
        SKIP_COUNT=$(( SKIP_COUNT + 1 ))
        continue
    fi

    macro_path="${ALL_MACROS[$i]}"
    macro_basename=$(basename "$macro_path")
    pos_raw="${macro_basename##*_x}"
    pos_mm="${pos_raw%mm.mac}"
    pos_mm="${pos_mm#+}"
    label="[$(printf '%03d' "$macro_idx")/${TOTAL_MACROS}] x=${pos_mm}mm"
    pos_tag="x${pos_mm}mm"

    pos_dir="${RUN_DIR}/outputs/${pos_tag}"
    mkdir -p "$pos_dir"
    cp "$macro_path" "${RUN_DIR}/macros_used/${macro_basename}"
    log_file="${RUN_DIR}/logs/run_${pos_tag}.log"

    STAT_FILE=$(mktemp)
    run_one_macro \
        "$label" "$pos_mm" "$macro_path" "$pos_dir" "$log_file" "$STAT_FILE"
    pos_result=$(cat "$STAT_FILE")
    rm -f "$STAT_FILE"

    if [[ "$pos_result" == "OK" ]]; then
        PASS_COUNT=$(( PASS_COUNT + 1 ))
    else
        FAIL_COUNT=$(( FAIL_COUNT + 1 ))
        if [[ $STOP_ON_ERROR -eq 1 ]]; then
            {
                echo ""
                echo "ERROR: ${label} falló (${pos_result}). STOP_ON_ERROR=1 → abortando."
                echo "       Para continuar ante fallos: STOP_ON_ERROR=0 $0"
            } | tee -a "$MASTER_LOG"
            break
        fi
    fi
done

# ── Final summary ─────────────────────────────────────────────
END_DATETIME=$(date '+%Y-%m-%d %H:%M:%S')
{
    echo ""
    echo "============================================================"
    echo "  RESUMEN FINAL"
    echo "============================================================"
    echo "  Fin               : ${END_DATETIME}"
    echo "  Directorio salida : ${RUN_DIR}"
    echo "  OK                : ${PASS_COUNT}"
    echo "  FAIL              : ${FAIL_COUNT}"
    echo "  Saltadas          : ${SKIP_COUNT}"
    echo "  Resumen TSV       : ${SUMMARY_TSV}"
    echo "============================================================"
} | tee -a "$MASTER_LOG"

[[ $FAIL_COUNT -gt 0 ]] && exit 1
exit 0
