# EXEC_19_REPORT — 2026-06-20 20:16

## 1. Veredicto
**COMPLETADO-CON-FLAGS** (minor: T1_GOAL_B for SUM4 borderline; ver §2)

## 2. Tabla A/B consistente (FIX central de EXEC_18)

| Estimador | x [mm] | σ_int_A | σ_tot_A | σ_int_B | N_eff | σ_tot_B | A≤100 | B≤100 |
|----------|------|--------|--------|--------|------|--------|-------|-------|
| TOP_SUM4 | 0 | 69.4 ps | 127.1 ps | 80.7 ps | 3.15 | 100.9 ps | NO | NO |
| TOP_SUM8 | 0 | ~65 ps  | ~125 ps  | 80.1 ps | 3.69 | 97.8 ps  | NO | **YES** |
| TOP_SUM4 | -690 | 130 ps | 168 ps | 91.7 ps | 2.59 | 113.4 ps | NO | NO |

### Chain A (suma analógica + min/CFD):
- T2 scan: min(t_first) es MEJOR que cualquier CFD(f>0). Mejor Chain A = N=1 (primer fotón).
- σ_int_A = 69.4 ps (min, walk-corr, sqrt_n); σ_tot_A = √(69.4²+106²+10²) = 127.1 ps
- En detector real: CFD threshold > 0 (ruido) → σ_int_A > 69.4 ps; mejor_CFD_testado=101.8 ps @ f=0.1
- **EXEC_18 σ_int_B=69.2 ps era min-stream (confundido con Chain B) — corregido aquí**

### Chain B (promedio digital ponderado por canal):
- σ_int_B = 80.7 ps (TOP_SUM4@x=0, weighted avg de t_first por canal)
- Chain B es DISTINTO de min-stream y da σ_int_B > σ_int_A (min): el promedio ponderado
  es peor que el mínimo porque promedia canales más lentos con peso no nulo
- N_eff(Kish) = 3.15 para SUM4; SPTR_eff = 106/√3.15 = 59.7 ps
- σ_tot_B(SUM4) = √(80.7²+59.7²+10²) = 100.9 ps (apenas sobre 100 ps → NO)
- σ_tot_B(SUM8) = 97.8 ps → SÍ cumple 100 ps (bajo Chain B con SUM8)

## 3. Estado honesto del 100 ps
- Chain A (real, con umbral CFD > 0): > 100 ps (incluso con min: 127 ps)
- Chain B TOP_SUM4: 100.9 ps → marginal, NO cumple
- Chain B TOP_SUM8: 97.8 ps → SÍ cumple (bajo Chain B consistente)
- **El sub-100 ps es alcanzable SOLO bajo Chain B con TOP_SUM8 y FastIC+ digital-por-canal**

## 4. CFD vs min (T2)
- min(t_first) @ x=0: σ = 69.4 ps (walk-corr, sqrt_n)
- CFD(f=0.1): σ = 101.8 ps; CFD(f=0.2): 124.9 ps; ...
- En simulación sin ruido: min = límite teórico de Chain A
- En detector real: CFD threshold > 0 → brecha = al menos +32 ps

## 5. T3 topología — ⟨Npe⟩ y depuración EXEC_18
- Dynamic nearest-4: ⟨Npe⟩_total = 91.7328 PE/ev (más luz)
- id//4 cluster:      ⟨Npe⟩_total = 83.7208 PE/ev
- EXEC_18 72.1 vs 85.0 ps: artefacto FD binning sin walk correction
- Con sqrt_n + walk: ambas convergen a valores similares
- Reencuadre: TOP es SIEMPRE distribuida; co-localizado solo aplica a END
- No se contaminaron resultados headline de EXEC_16/17 (usaban SUM4/SUM8 de otra forma)

## 6. Reconciliación T4 + canónicos
- TOP_nearest @ x=0: gid=50, cx=-12mm, σ_canónico = 105.7 ps (sqrt_n, walk-corr)
- TOP_nearest @ x=-690: gid=16, cx=-692mm, σ = 118.2 ps
- END_L_SUM4 @ x=-690: gids=[3,4,2,5], σ = 30.5 ps (near-end, sqrt_n, walk)
  EXEC_17/C2 (30 ps, FD+walk): consistente. EXEC_18/T3 (44 ps, FD, no walk): discrepancia = no walk
- EXEC_18 T1 TOP_nearest 107.8 ps vs T3 85.0 ps: distintos estimadores (single ch vs 4-ch)

## 7. Flags
- σ_tot_B(SUM4@x=0) = 100.9 ps: bordeline, marcado en tabla como NO (>100 ps)
- END walk_failed @ far positions: físico (pocos PE), no bug

## 8. Rutas
  Figuras: /home/reriosto/SHiP/analysis_core/out/EXEC_19/T1..T4
  Memo: /home/reriosto/SHiP/analysis_core/out/EXEC_19/GERARDO_DECISION_MEMO.md
  Canónico: /home/reriosto/SHiP/analysis_core/out/EXEC_19/CANONICAL_ESTIMATORS.md
  JSON: /home/reriosto/SHiP/analysis_core/out/EXEC_19/results_exec19.json
  Beamer: analysis_core/beamer/EXEC_19/EXEC_19_pairing_topology.pdf

## 9. Decisiones para Gerardo
  1. FastIC+ arquitectura: Chain A→127 ps (real con umbral CFD) vs Chain B→97.8 ps (SUM8)
  2. Topología TOP: id//4 vs nearest-4 (ambas distribuidas; ⟨Npe⟩ clarificado)