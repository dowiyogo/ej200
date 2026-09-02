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

**[OPEN] FIX-03 — w_END unresolved at n=5000:**
CP-B w-scan (n=5000, 1000-resample bootstrap) shows σ_IQR(w) is flat within one bootstrap SE
over w ∈ [0, 0.10] (total variation 0.26 ps ≈ SE 0.17 ps). Three covariance conventions give
w_END ∈ {0.013, 0.051, 0.086} with σ_IQR^event within 0.18 ps. Resolution would require
≥50k events at x=0 under the same hybrid geometry; this currently changes no design conclusion.
Sidecar: analysis/optim/root_best_est/blue_wscan_x0.{csv,meta.json,root}

**[OPEN] FIX-01b — figM1/figM2 figure sidecars missing (phase_ab.py):**
`analysis/optim/phase_ab.py` saves `fig1_sigma_vs_k_x0.png`, `fig2_optimal_vs_x.png`,
`fig3_sigma_vs_efficiency.png` but emits no CSV or `.meta.json` sidecars.
figM1/figM2/fig_sigma_t_x/fig_npe_x in `presentations/v6/figs/` lack provenance metadata.
Adding proper sidecar infrastructure to phase_ab.py is estimated >30 min;
deferred pending re-simulation campaign with R=0.98.

**SPTR diapositiva I1 (corregida en `4851c43`→ tex en siguiente commit):** El valor
"~100 ps" no tenía procedencia (DS105 no tabula SPTR). Reemplazado por medida publicada
del 6×6 mm² NUV-MT (Lee et al., IEEE TRPMS 9(4) 2025): 137±4 ps FWHM intrínseco
(σ≈58 ps); 172 ps FWHM con electrónica (σ≈73 ps); Vbias=48 V (~15.5 V OV).
Conversiones FWHM→σ explicitadas en la diapositiva. Fórmula de promedio marcada como
pendiente (estadístico de orden ≠ promedio). Crosstalk óptico 23% añadido a efectos
no incluidos. Ver `docs/branch_diagnosis/SPTR_PROVENANCE.md`.

**SPTR diapositiva I1 — actualización disponible (banco propio, 2026-09-01):** El barrido
de intensidad láser a V_OV = 10 V (punto de operación real) y umbral FastIC+ = 35 da,
por ajuste σ²(N) = SPTR²/N + σ_elec²: **SPTR ≈ 85–92 ps RMS** a 10 V OV; σ_elec ≤ 18 ps.
Coherente con Lee 2025 (58 ps a 15.5 V OV): SPTR empeora al bajar OV, dirección esperada.
El valor original ~100 ps (sin procedencia) resultó razonable como cota conservadora.
La diapositiva I1 puede actualizarse con el rango propio indicando OV y carácter derivado.
**No editar el `.tex` hasta que la medida formal con datos de EOS esté disponible.**
Ver `analysis/timing/TODO.md` y `docs/branch_diagnosis/ELECTRONICS_PARAMETERS.md`.

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
