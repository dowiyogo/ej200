# end_vikuiti_2026-08-16 — 16 agosto 2026

**Estado:** OBSOLETO (template con placeholders; no completado; supersedido por v6)
**Fase óptica:** Fase 7 — referencia `../../build_end_vikuiti` con air gap + Mylar,
`dielectric_dielectric` R=0.95, `GROUPVEL_air = c/n_bar = 18.97 cm/ns`
(commit `5576687`, 2026-08-14). El build_end_vikuiti vigente el 2026-08-16 fue
compilado con Fase 7 (scan_end_vikuiti ejecutado ese mismo día).
**Dataset:** `build_end_vikuiti/` — figuras previstas en el template pero posiblemente
no generadas en el momento del commit ("placeholders ready to fill after scan completes").

## Defectos conocidos

**Reflectividad R = 0.95 (corregida en `c7acb7a`):** Datos producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()` — error de configuración; spec Vikuiti 3M ESR es R≈0.98.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Dataset necesita
re-simulación con R=0.98. Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

El deck es un template: el commit original indica que las figuras son placeholders
preparados antes de que el scan completara.

## Por qué se conserva

Plantilla para la comparación End-only EJ-200/EJ-204/EJ-230 con air gap Vikuiti.
El contenido fue absorbido y completado en v6.
