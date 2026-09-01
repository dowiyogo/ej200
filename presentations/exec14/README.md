# exec14 — EJ-204 vs EJ-230 Beamer Compendium — 15–16 junio 2026

**Estado:** OBSOLETO (datos de rama paralela feat/endonly-mylar, previa al rediseño óptico principal)
**Fase óptica:** Fuera de la cronología principal — ver descripción.

Configuración de `feat/endonly-mylar` @ SHA `3ae135f` (pre-fix skin surface, 2026-06-14):
- Reflector: volúmenes físicos Mylar (`G4LogicalBorderSurface`, `dielectric_dielectric`, R=0.95)
- Geometría: End-only (sin SiPMs TOP), sin cámara de aire separada entre bar y reflector
- Esta configuración no coincide con ninguna de las Fases 3–7 de la cronología de
  `feat/bar-end-vikuiti`; pertenece a la rama paralela `feat/endonly-mylar`.
- El skin surface fix global se aplicó el 2026-06-18 (commit `84e902c` en main);
  los datos EXEC_14 son del 2026-06-14, cuatro días antes.

**Dataset:**
- EndOnly EJ-204: `/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/` (31 posiciones)
- EndOnly EJ-230: `/home/reriosto/SHiP/t0minidaq/endonly_mylar_230/` (31 posiciones)
- EndTop EJ-204: `/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000/` (86 canales)
- EndTop EJ-230: `/home/reriosto/SHiP/t0minidaq/results_ej230/data/` (86 canales, τ_d=1.5 ns)

Los datos ROOT están en t0minidaq (fuera de este repositorio). Los números de análisis
están preservados en `analysis/exec14/outputs/` (CSV y meta.json).

## Defectos conocidos

- Los datos EndOnly usan la geometría de rama paralela previa al rediseño óptico principal.
  Sus resultados σ_t NO son comparables directamente con los de talk_v6 (GEN-3, Fase 7).
- No existe dataset EndTop equivalente para EJ-230 en la misma época que EJ-204 EndOnly.
- Los datos EndTop de EJ-204 (`exec07_endtop_2000`) y EJ-230 (`ej230_endtop`) usan la
  geometría EndTop con SiPMs TOP (86 canales), producida con la rama sslg4/exec07.

## Por qué se conserva

Única comparación directa EJ-204 vs EJ-230 con análisis SUM4 y figuras de order-statistics
(F7: ⟨t_n⟩ vs n por canal Top). Los CSV de análisis en `analysis/exec14/outputs/`
preservan los números sin necesidad de reejecutar el análisis sobre los ROOT files externos.

## Fuente

Importado desde `https://github.com/dowiyogo/ej200_orchestrator` (19 commits, rama main).
Deck canónico: `exec14_deck_aclarado.tex` (615 líneas, commit EXEC_13 close + EXEC_15).
Documentación de QA: `docs/branch_diagnosis/exec14_qa/`.
Framework de análisis: `analysis/exec14/`.
