# optim_2026-08-16 — 16 agosto 2026

**Estado:** OBSOLETO (primera iteración; supersedido por optim_2026-08-17 y v6)
**Fase óptica:** Fase 7 — datos de `analysis/optim/phase_ab_optimal.csv`
generado desde `results/scan_end_vikuiti/` (scan ejecutado 2026-08-16, dos días
después del commit de Fase 7 `5576687`, 2026-08-14).
**Dataset:** `analysis/optim/phase_ab_optimal.csv` ← `results/scan_end_vikuiti/`;
geometría End-only (sin SiPMs TOP); EJ-200, EJ-204, EJ-230.

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

Primera iteración del análisis de optimización (parámetros α, β del estimador SUM4);
supersedido por optim_2026-08-17 con readout EndSparseTop combinado.

## Por qué se conserva

Registro de la primera iteración del análisis de optimización de parámetros;
punto de partida del flujo de análisis que culmina en v6.
