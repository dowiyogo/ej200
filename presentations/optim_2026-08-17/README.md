# optim_2026-08-17 — 17 agosto 2026

**Estado:** OBSOLETO (supersedido por best_est_2026-08-17 y v6)
**Fase óptica:** Fase 7 — datos de `analysis/optim/` (scan_end_vikuiti 2026-08-16;
Fase 7 commit `5576687`, 2026-08-14).
**Dataset:** `analysis/optim/` — Phase II: combinación END-only + END+TOP sparse
(scan_end_vik_sparse_top_v2); EJ-200, EJ-204, EJ-230.

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

Segunda iteración del análisis de optimización; introduce readout EndSparseTop
combinado (Phase II). Supersedido por el análisis de best_est con ROOT verificado.

## Por qué se conserva

Registro de la segunda iteración del análisis que introduce el readout END+TOP
combinado y el estimador BLUE, antes de la verificación con ROOT en best_est_2026-08-17.
