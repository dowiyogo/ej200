# summary_2026-08-16 — 16 agosto 2026

**Estado:** OBSOLETO (síntesis preliminar; supersedido por v6)
**Fase óptica:** Fase 7 — references `../../build_end_vikuiti` (air gap + Mylar,
`dielectric_dielectric` R=0.95, commit `5576687`, 2026-08-14) y `../../build_t0minidaq`
(geometría EndTop, misma fase óptica). Datos del 2026-08-16.
**Dataset:** `build_end_vikuiti/` (End-only) + `build_t0minidaq/` (EndTop);
EJ-200, EJ-204, EJ-230; configuraciones TIR y Vikuiti.

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

El deck es una síntesis preliminar del 2026-08-16 con σ(T₀)=59 ps reportado para
EJ-230 + END + Vikuiti — valor del análisis rápido inmediato al scan, previo a la
optimización detallada de los parámetros α, β del estimador.

## Por qué se conserva

Síntesis multimaterial (EJ-200/204/230) y multi-readout (TOP vs END, TIR vs Vikuiti)
realizada el mismo día que el scan; registro del primer panorama completo de resultados.
