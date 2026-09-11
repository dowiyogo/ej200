# EXEC_18_REPORT — 2026-06-20 16:53

## 1. Veredicto
**COMPLETADO-CON-FLAGS** (ver §8)

## 2. Etiqueta del 72.6 ps resuelta (T4)
- **Definición**: TOP_SUM4_N1 = 4 TOP SiPMs (max ⟨Npe⟩) mergeados, primer fotón = min(t_1st) entre los 4 canales
- Label EXEC_16 correcto. "SUM4" = nº de canales, N=1 = primer fotón (implícito)
- EXEC_17/C3 usó N como umbral de fotones (2..8) dentro del stream mergeado — distinto
- Recomputado (fd binning): TOP_nearest_N1=116.8 ps | TOP_SUM4_N1=85.0 ps
- EXEC_16 (sqrt_n binning, sin walk): 107.2 ps | 72.6 ps
- Con walk (T1): → 107.8 ps | → 69.2 ps
- **Flagship corregido**: TOP_SUM4_N1 ≈ 69.2 ps (walk-corrected, fd binning)

## 3. Ranking con y sin time-walk (T2)
- TOP_SUM4 @ x=0: raw 85.0 ps → corr 69.2 ps
- TOP_SUM8 @ x=0: raw ~84 ps → corr ~68 ps
- Ranking NOT inverted: TOP_SUM4 sigue mejor/igual que TOP_nearest tras walk
- SUM8 (68.1 ps) ≈ SUM4 (69.2 ps) con walk — SUM8 ahora marginalmente mejor (física esperada)
- La anomalía EXEC_16 (SUM8>SUM4) se explica por el sesgo de time-walk; corregida, desaparece

## 4. END_SUM4 co-localizado vs TOP nearest (T3)
- @ x=-690 (near-end): END_L_SUM4=30.2 ps | TOP_nearest=98.8 ps
- @ x=0 (far from END): END_L_SUM4=944.2 ps | TOP_nearest=107.8 ps
- END_SUM4 co-localized: analogía correcta con TB Constanza (misma face, misma luz)
- Lejos del extremo, END da >900 ps; cerca del extremo, ~44 ps near-end

## 5. HOOK_WALK: corrección de time-walk (T1)
  Slewing var: NPE total en stream. Forma: α+β/√s (parametric). Ancla: mediana(NPE).
  - TOP_nearest_N1: form=parametric | σ_raw=116.8 → σ_corr=107.8 ps | Δσ=-9.0 ps | r_corr=-0.0 | ok
  - TOP_SUM4_N1: form=parametric | σ_raw=85.0 → σ_corr=69.2 ps | Δσ=-15.7 ps | r_corr=-0.0 | ok
  - TOP_SUM8_N1: form=parametric | σ_raw=84.2 → σ_corr=68.1 ps | Δσ=-16.1 ps | r_corr=-0.0 | ok
  - END_L_SUM4_N1: form=none | σ_raw=944.2 → σ_corr=944.2 ps | Δσ=0.0 ps | r_corr=-0.2 | walk_failed
  - END_R_SUM4_N1: form=none | σ_raw=1061.3 → σ_corr=1061.3 ps | Δσ=0.0 ps | r_corr=-0.2 | walk_failed

## 6. Tabla σ_total A vs B con N_eff Kish (T5) [DECISIÓN GERARDO]
| Estimador | N_ch | N_eff | σ_int | σ_tot Model A | σ_tot Model B |
|----------|------|-------|-------|--------------|--------------|
| TOP_nearest_N1 | 1 | 1.0 | 107 ps | 151 ps | 151 ps |
| TOP_SUM4_N1 | 4 | 3.42 | 73 ps | 129 ps | 93 ps |
| TOP_SUM8_N1 | 8 | 4.59 | 74 ps | 130 ps | 90 ps |

- Model A (analog+CFD, SPTR=106 ps flat): goal 100 ps → NO para todos
- Model B (digital, SPTR/√N_eff): TOP_SUM4 N_eff=3.42→93.0 ps (SÍ)
- TOP_SUM8 N_eff=4.59→89.9 ps (SÍ)
- **σ_int incluye 20 ps SiPMSD jitter — no re-sumar**
- Modelo B es piso optimista; cumplir depende de la arquitectura FastIC+

## 7. T6 — Topología SUM
- id//4 fijo (co-localizado cluster): σ=72.1 ps
- Dynamic nearest-4 (distribuido):    σ=85.0 ps
- Co-localizado es mejor @ x=0 (misma fotónica); distribuido mezcla trayectorias distintas

## 8. Covarianza EndTop (T7)
- @ x=0: inv-varianza ≈ 100 ps | GLS con cov ≈ 100 ps (diferencia ínfima al centro)
- Diferencia GLS vs inv-var apreciable solo cerca de los extremos
- T7 units: valores computados en ns, ×1000 para ps

## 9. Flags emitidos
  - END_L/R_SUM4_N1 @ x=0: walk_failed (insuf. PE at far end from END face)

## 10. Rutas y hashes
  Figuras+CSVs: /home/reriosto/SHiP/analysis_core/out/EXEC_18/T1..T7
  GERARDO_DECISION_MEMO.md: /home/reriosto/SHiP/analysis_core/out/EXEC_18/GERARDO_DECISION_MEMO.md
  results_exec18.json: /home/reriosto/SHiP/analysis_core/out/EXEC_18/results_exec18.json
  Tags: EXEC_18-pre-t1, EXEC_18-pre-t3, EXEC_18-final

## 11. Decisiones para Gerardo
  Ver GERARDO_DECISION_MEMO.md
  1. FastIC+ architecture: Model A (analog,~129ps) vs Model B (digital,~90-93ps)
  2. SUM topology: co-localized id//4 (72.1ps) vs dynamic nearest-4 (85.0ps) @ x=0