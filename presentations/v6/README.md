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

Ningún defecto de fase óptica. Inconsistencias narrativas documentadas en
`docs/branch_diagnosis/TALKV6_CONSISTENCY.md` y `docs/branch_diagnosis/REFLECTOR_LABEL.md`:

- Figuras figM1/M2/fig_sigma_t_x/fig_npe_x etiquetadas como "GEN-2" en CONFIGURATION_AUDIT.md
  pero provienen de `scan_end_vikuiti` con `build_end_vikuiti` Fase 7 (GEN-3 física, END-only).
- Atribución B1 (Δσ_END=3.2 ps a SiPMs TOP) plausible según CONFIGURATION_AUDIT pero
  no demostrada explícitamente en el deck (comparación confunde dos variables simultáneas).
- "Vikuiti ESR" como nombre del reflector: la simulación usa G4_MYLAR R=0.95, no R≈0.98 (spec).
  Propuesta de sustitución de texto en REFLECTOR_LABEL.md, pendiente de aprobación.
- Línea 860 (`talk_v6.tex:860`): referencia a un comentario de código ya corregido por commit `71eeb20`.

## Por qué se conserva

Versión actual de la presentación principal del proyecto.
