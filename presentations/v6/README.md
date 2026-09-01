# presentation v6 — 2026-09-01

**Estado:** ACTUAL
**Fase óptica:** Fase 7 — air gap 0.10 mm + Mylar 0.05 mm, `dielectric_dielectric` R=0.95,
`GROUPVEL_air = c/n_bar = 18.97 cm/ns` (commit `5576687`, 2026-08-14).
Material de reflector simulado: G4_MYLAR con border surface `dielectric_dielectric` R=0.95.
El volumen se llama `VikuitiXXXLV`; el nombre es histórico (Vikuiti ESR como candidato físico).
**Dataset:**
- END-only (GEN-3 física): `results/scan_end_vikuiti/` (2026-08-16) — figM1/M2, fig_sigma_t_x, fig_npe_x
- END+TOP sparse (GEN-3): `results/scan_end_vik_sparse_top_v2/` — FINAL_NUMBERS.md

## Defectos conocidos

**Reflectividad R = 0.95 en datos históricos (corregida en `c7acb7a`):** Los datasets
`scan_end_vikuiti` y `scan_end_vik_sparse_top_v2` fueron producidos con `REFLECTIVITY=0.95`
en `CreateBarSkinReflector()`. Los números del deck (σ_END=53.68 ps, σ_BLUE=15.21 ps,
σ_TOP=15.20 ps) corresponden a R=0.95 y necesitan re-simulación con R=0.98 para validarse.
Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01).
Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md`.

Inconsistencias narrativas documentadas en `docs/branch_diagnosis/TALKV6_CONSISTENCY.md`:

- Figuras figM1/M2/fig_sigma_t_x/fig_npe_x etiquetadas como "GEN-2" en CONFIGURATION_AUDIT.md
  pero provienen de `scan_end_vikuiti` con `build_end_vikuiti` Fase 7 (GEN-3 física, END-only).
- Atribución B1 (Δσ_END=3.2 ps a SiPMs TOP) plausible según CONFIGURATION_AUDIT pero
  no demostrada explícitamente en el deck (comparación confunde dos variables simultáneas).
- "Vikuiti ESR" como nombre del reflector: R corregida a 0.98 (commit `c7acb7a`); el deck puede
  mantener "Vikuiti ESR (R=0.98)" sin discrepancia en nuevas simulaciones. `REFLECTOR_LABEL.md`
  superado — ver `REFLECTIVITY_CHANGE.md §6`. Macro `\Rvikuiti{0.95}` y ocurrencias hardcoded
  en el tex necesitan actualización al regenerar el deck con datos R=0.98.
- Línea 860 (`talk_v6.tex:860`): describe la discrepancia R=0.95/0.98 como defecto activo —
  obsoleta una vez el deck se regenere con datos R=0.98.

## Por qué se conserva

Versión actual de la presentación principal del proyecto.
