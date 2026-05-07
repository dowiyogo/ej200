#!/usr/bin/env bash
# run_audit.sh — ejecuta el diagnostico completo de yield de fotones
# Uso: bash diag/run_audit.sh [ruta_al_build]
# Ejemplo: bash diag/run_audit.sh build/

set -euo pipefail

BUILD="${1:-build}"
SIM="${BUILD}/ej200_bar_sim"
ROOT_EXE="${ROOT_EXE:-root}"

last_run_summary() {
    awk '
        /^=== EJ-200 Bar Run Summary ===/ { block=$0; capture=1; next }
        capture {
            block=block ORS $0
            if ($0 ~ /^==============================/) {
                last=block
                capture=0
            }
        }
        END { if (last) print last }
    ' "$1"
}

last_boundary_census() {
    awk '
        /^=== Boundary Census/ { block=$0; capture=1; next }
        capture {
            block=block ORS $0
            if ($0 ~ /^====================================/) {
                last=block
                capture=0
            }
        }
        END { if (last) print last }
    ' "$1"
}

if [[ ! -x "$SIM" ]]; then
    echo "ERROR: ejecutable no encontrado en $BUILD/"
    echo "Compila primero: cmake --build $BUILD -j"
    exit 1
fi

mkdir -p diag

echo "=============================================="
echo " DIAGNOSTICO DE YIELD — EJ-200 SHiP TD"
echo " $(date)"
echo "=============================================="

echo ""
echo "[1/4] Corrida de simulacion (run.mac, x=0)..."
"$SIM" -m macros/run.mac 2>&1 | tee diag/sim_output.log

echo ""
echo "--- Resumen de la simulacion ---"
last_run_summary diag/sim_output.log || true
last_boundary_census diag/sim_output.log || true

echo ""
echo "[2/4] Contando entradas en el ROOT generado..."
ROOTFILE=$(ls photon_hits_run*.root 2>/dev/null | tail -1 || true)
if [[ -z "$ROOTFILE" ]]; then
    echo "  AVISO: no se encontro photon_hits_run*.root en directorio actual"
    echo "  Buscando en $BUILD/..."
    ROOTFILE=$(ls "${BUILD}"/photon_hits_run*.root 2>/dev/null | tail -1 || true)
fi

if [[ -n "$ROOTFILE" ]]; then
    echo "  Archivo ROOT: $ROOTFILE"
    if command -v "$ROOT_EXE" &>/dev/null; then
        N_ENTRIES=$("$ROOT_EXE" -l -b -q \
            "diag/yield_audit.C(\"${ROOTFILE}\")" \
            2>/dev/null | grep "Entradas en TTree" | awk '{print $5}' || echo "ERROR")
        echo "  Entradas en TTree: ${N_ENTRIES}"
    else
        echo "  AVISO: ROOT no encontrado. No se pudo contar el TTree."
    fi
else
    echo "  AVISO: no se encontro archivo ROOT para analizar"
fi

echo ""
echo "[3/4] Ejecutando analisis ROOT..."
if [[ -n "${ROOTFILE:-}" ]] && command -v "$ROOT_EXE" &>/dev/null; then
    "$ROOT_EXE" -l -b -q "diag/yield_audit.C(\"${ROOTFILE}\")" 2>&1 | tee diag/root_audit.log
else
    echo "  AVISO: ROOT no encontrado o no hay archivo ROOT. Saltando analisis." | tee diag/root_audit.log
fi

echo ""
echo "[4/4] Generando informe..."
{
    echo "YIELD AUDIT REPORT"
    echo "=================="
    echo "Fecha: $(date)"
    echo "Branch: $(git branch --show-current 2>/dev/null || echo 'desconocido')"
    echo "Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'desconocido')"
    echo ""
    echo "--- Salida de simulacion ---"
    last_run_summary diag/sim_output.log || echo "(no encontrado)"
    echo ""
    last_boundary_census diag/sim_output.log || echo "(Boundary Census no disponible)"
    echo ""
    echo "--- Analisis ROOT ---"
    cat diag/root_audit.log 2>/dev/null || echo "(no disponible)"
    echo ""
    echo "--- Interpretacion esperada ---"
    echo "Referencia experimental (Betancourt 2020):"
    echo "  ~175 fotones/evento en ambos ends combinados para muon central"
    echo "  sigma_t ~ 85 ps"
    echo ""
    echo "Criterios de diagnostico:"
    echo "  N_pe/evento < 10   -> Problema severo de transporte o generacion"
    echo "  N_pe/evento 10-50  -> Reflector insuficiente o absorcion bulk alta"
    echo "  N_pe/evento 50-150 -> Geometria SiPM o PDE limitando"
    echo "  N_pe/evento > 150  -> Aceptable para comparar con experimento"
    echo ""
    echo "Relacion Bar->Mylar / Mylar->SiPM esperada:"
    echo "  Si Bar->Mylar >> Mylar->SiPM * 100 -> fotones atrapados o absorbidos"
    echo "  Si Mylar->World > Mylar->SiPM * 10 -> reflector insuficiente"
} > diag/yield_audit_report.txt

echo ""
echo "=============================================="
echo " Informe guardado en: diag/yield_audit_report.txt"
echo " Figuras guardadas en: diag/yield_audit.pdf"
echo "=============================================="
cat diag/yield_audit_report.txt
