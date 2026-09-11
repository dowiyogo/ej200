# EXEC_20_REPORT — 2026-06-20 22:34

## 1. Veredicto
**COMPLETADO**

## 2. END-only baseline (T1)
- Estimador: t_avg = (t_L + t_R)/2, walk-corr, sqrt_n, all 8 SiPMs per face
- σ_END(x=0)   = 882 ps (0.37+0.38 PE/ev — photon-starved)
- σ_END(x=-690) = 609 ps
- v_eff = -14.70 cm/ns

## 3. TOP contribution (T3)
- Mejora EndTop/END-only @ x=0: 86.2%
- Mejora media (31 pos): 84.5%
- σ_TOP @ x=0 = 105.7 ps; σ_EndTop = 122.2 ps

## 4. Reconciliation Gerardo (T4)
- Nuestro t_avg requiere AMBOS extremos → baja eficiencia en extremos de scan
- Gerardo probablemente usa un estimador diferente (confirmar)
- Near-end σ_END_L_SUM4 ≈ 30.5 ps (EXEC_19 T4)

## 5. Build END-only (T2a/T2b)
- T2a: PASS — OPSC-101, 16 SiPMs [0,15], no TOP hits, reflector skin fix ✓
- Optical diff: R=0.90 (endonly) vs R=0.95 (reference); jitter=0 in both
- T2b scan: Scan output at /home/reriosto/SHiP/ej200_endonly/output/endonly_mylar_msi_20260620_222723 (check run_metadata.txt for status)

## 6. Canonical fix (T5)
- Binning: sqrt_n (FD artifact EXEC_18 resolved)
- dynamic nearest-4: 69.4 ps ≤ id//4: 76.4 ps (as expected)

## 7. Flags
  ninguno

## 8. Rutas
  Figuras: /home/reriosto/SHiP/analysis_core/out/EXEC_20/T1..T5
  Memo Gerardo: /home/reriosto/SHiP/analysis_core/out/EXEC_20/GERARDO_DECISION_MEMO.md
  Canónico: /home/reriosto/SHiP/analysis_core/out/EXEC_20/T5/CANONICAL_ESTIMATORS.md
  JSON: /home/reriosto/SHiP/analysis_core/out/EXEC_20/results_exec20.json

## 9. Decisiones para Gerardo
  1. Topología SUM: dynamic-4 vs id//4 (ambas distribuidas para TOP)
  2. Estimador END-only: t_avg (position-independent) vs single-end (mejor por extremo)
  3. Confirmar estimador exacto de Gerardo para comparación apples-to-apples
  4. Electrónica (SPTR/FastIC): diferida por decisión de René