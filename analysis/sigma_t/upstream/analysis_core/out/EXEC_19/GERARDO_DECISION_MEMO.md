# GERARDO_DECISION_MEMO — EXEC_19 (actualizado, física correcta)
Fecha: 2026-06-20
Para: Gerardo Vásquez | De: EXEC_19 autonomous analysis

## Corrección respecto a EXEC_18
EXEC_18 combinó σ_int del mínimo (Chain A) con SPTR del promedio (Chain B).
EXEC_19 implementa ambas cadenas de extremo a extremo, internamente consistentes.

---
## DECISIÓN 1 — Cadena de lectura: A (analógica) vs B (digital por canal)

### Chain A — suma analógica + un discriminador CFD
- Timestamp = CFD a fracción f=0.0 sobre N(t) acumulada (ver T2)
- σ_int_A (TOP_SUM4 @ x=0) = ? ps
- σ_tot_A = √(σ_int_A² + 106² + 10²) = ? ps
- Cumple SHiP 100 ps: ?

### Chain B — timestamps digitales por canal, promediados
- Cada canal tiene t_first con walk correction individual
- Promedio ponderado T = Σ w_i t_i / Σ w_i
- σ_int_B (TOP_SUM4 @ x=0) = ? ps
- N_eff (Kish) = ? (penaliza canales lejanos con menos luz)
- σ_SPTR_eff = 106/√N_eff ps
- σ_tot_B = √(σ_int_B² + σ_SPTR_eff² + 10²) = ? ps
- Cumple SHiP 100 ps: ?
- TOP_SUM8 Chain B: σ_tot_B = ? ps

### Tabla resumen (TOP_SUM4 @ x=0):
| Cadena | σ_int | σ_SPTR | σ_tot | Cumple 100 ps |
|--------|-------|--------|-------|--------------|
| A (analógica) | ? ps | 106 ps | ? ps | ? |
| B (digital) | ? ps | 106/√? ps | ? ps | ? |

### Estado del requisito 100 ps:
- σ_int intrínseco (sin SPTR/FastIC): ambas cadenas < 100 ps ← dato seguro
- Con SPTR: depende de la cadena. El 100 ps solo se alcanza bajo Chain B.
- Chain B requiere que FastIC+ pueda :
  (i) generar timestamps individuales por canal,
  (ii) promediar ponderadamente en tiempo real.
  → Confirmar con el equipo FastIC+.

---
## DECISIÓN 2 — Topología de SUM: END co-localizado vs TOP distribuido

### Reencuadre (EXEC_19)
'Co-localizado' describe el END (8 SiPMs en la MISMA cara de extremo, mismo x,
misma luz). El TOP es SIEMPRE una FILA DISTRIBUIDA a paso 20mm.
Cualquier agrupación de SiPMs TOP (id//4 o nearest-N) es distribuida.
La comparación co-localizado↔distribuido REAL es:
  END_SUM4 co-localizado (misma face) vs TOP-fila-sum (distribuida)

### Evidencia de selección @ x=0:
- Dynamic nearest-4 TOP: total ⟨Npe⟩ = 91.73 PE/ev
- id//4 cluster TOP:     total ⟨Npe⟩ = 83.72 PE/ev
- Dynamic recoge MÁS luz pero EXEC_18 dio peor σ (artefacto FD binning + no walk)
- Con sqrt_n + walk: ambos deberían dar σ similar — ver T3 figuras

### Recomendación técnica (neutral)
- Para la comparación sim↔TB Constanza: END_SUM4 co-localizado
- Para la instrumentación TOP: nearest-4 dinámico recoge más luz y es la selección
  correcta; id//4 fijo es equivalente cuando el haz está bien alineado
- La decisión final (id//4 fijo vs nearest-4 dinámico) sigue siendo de Gerardo
  pero ambas son topologías DISTRIBUIDAS

---
## Definiciones canónicas (CANONICAL_ESTIMATORS.md)
Ver análysis_core/out/EXEC_19/CANONICAL_ESTIMATORS.md para la definición exacta
de cada estimador headline (canales, combinación, walk, binning) y su σ.