# best_est_2026-08-17 — 17 agosto 2026

**Estado:** OBSOLETO (supersedido por v6, que incorpora y amplía este análisis)
**Fase óptica:** Fase 7 — datos de `analysis/optim/root_best_est/`
generado con código ROOT verificado (new_analysis_plots.C, m*=7, σ_END=50.494 ps);
Fase 7 (commit `5576687`, 2026-08-14); geometría END-only + END+TOP sparse.
**Dataset:** `analysis/optim/root_best_est/` ← `results/scan_end_vikuiti/` +
`results/scan_end_vik_sparse_top_v2/`; σ_END=53.68 ps (GEN-3 END+TOP),
σ_BLUE=15.21 ps, σ_TOP=15.20 ps.

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Los números canónicos
σ_END=53.68 ps, σ_BLUE=15.21 ps, σ_TOP=15.20 ps corresponden a R=0.95; necesitan
re-simulación con R=0.98 para validarse. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

Los números de este deck son los mismos que aparecen en FINAL_NUMBERS.md de v6 y en
talk_v6.tex (B1/B2 slides).

## Por qué se conserva

Registro de la best estimate canónica con análisis ROOT verificado y valores
explícitos hardcodeados (σ_END=53.68 ps, σ_BLUE=15.21 ps, σ_TOP=15.20 ps).
Fuente primaria de los números presentados en v6.
