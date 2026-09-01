# presentation v5 — agosto 2026

**Estado:** OBSOLETO (supersedido por v6)
**Fase óptica:** Fase 7 — air gap 0.10 mm + Mylar 0.05 mm, `dielectric_dielectric` R=0.95,
`GROUPVEL_air = c/n_bar = 18.97 cm/ns` (commit `5576687`, 2026-08-14)
**Dataset:** `analysis/optim/` (END-only + END+TOP sparse; scan_end_vikuiti 2026-08-16
y scan_end_vik_sparse_top_v2); EJ-230, Vikuiti R=0.95

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

El deck es versión intermedia, suprimido por v6.

## Por qué se conserva

Versión intermedia que introduce el estimador BLUE END+TOP y la comparación
de materiales EJ-200/EJ-204/EJ-230 antes de la consolidación en v6.
