# presentation v4 — agosto 2026

**Estado:** OBSOLETO (supersedido por v5 y v6)
**Fase óptica:** Fase 7 — air gap 0.10 mm + Mylar 0.05 mm, `dielectric_dielectric` R=0.95,
`GROUPVEL_air = c/n_bar = 18.97 cm/ns` (commit `5576687`, 2026-08-14)
**Dataset:** `analysis/optim/phase_ab_optimal.csv` ← `results/scan_end_vikuiti/` (2026-08-16);
best-estimate: `analysis/optim/root_best_est/` (fig_end_mscan.pdf, m*=7, σ=50.494 ps)

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

El deck es versión intermedia, suprimido por v5 y v6.

## Por qué se conserva

Versión intermedia que documenta la trayectoria napkin → Geant4 con estimadores
de primer orden (IQR, BLUE) antes de la introducción del estimador de orden-estadístico.
