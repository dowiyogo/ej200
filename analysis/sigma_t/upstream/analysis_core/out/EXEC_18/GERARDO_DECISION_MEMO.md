# GERARDO_DECISION_MEMO — EXEC_18
**Fecha**: 2026-06-20
**Para**: Gerardo Vásquez (SHiP Timing / OPSC-101 analysis)
**De**: EXEC_18 autonomous analysis, EJ-204 EndTop 31 pos × 5000 ev

Hay **dos decisiones** que no puede tomar el análisis numérico porque
dependen de la arquitectura real del detector y de la cadena de lectura.
A continuación se presenta la evidencia y la consecuencia en σ de cada opción.

---

## DECISIÓN 1 — Modelo de lectura: analógico (A) vs digital multi-canal (B)

### Contexto
El σ_intrínseco ≤ 72.6 ps (EXEC_16, TOP_SUM4_N1) NO incluye SPTR ni FastIC.
Para proyectar a σ_total = √(σ_int² + σ_SPTR² + σ_FastIC²), el valor de σ_SPTR
depende de cómo el FastIC+ procesa las señales.

### Modelo A — Suma analógica + discriminador único
- Las señales de los N canales se suman analógicamente ANTES del discriminador.
- Un solo discriminador ve la señal sumada → un solo SPTR = 106 ps.
- σ_SPTR_A = 106 ps (plano, independiente de N).
- σ_total(TOP_SUM4_N1) = √(72.6² + 106² + 10²) ≈ **129 ps**
- σ_total(TOP_SUM8_N1) = √(74.4² + 106² + 10²) ≈ **130 ps**
- Cumple SHiP 100 ps: **NO** (ambos ~129 ps)

### Modelo B — Timestamps digitales por canal, promediados
- Cada canal tiene su propio discriminador y SPTR = 106 ps.
- Los timestamps individuales se promedian → SPTR_eff ≈ 106/√N_eff.
- N_eff (Kish) penaliza canales lejanos con menor Npe.
- TOP_SUM4: N_ch=4, N_eff=3.42, σ_SPTR_B=57.3 ps
  → σ_total = **93 ps** (cumple 100 ps: **YES**)
- TOP_SUM8: N_ch=8, N_eff=4.59, σ_SPTR_B=49.5 ps
  → σ_total = **90 ps** (cumple 100 ps: **YES**)

### Consecuencia
| Modelo | TOP_SUM4 | TOP_SUM8 | Cumple 100 ps |
|--------|---------|---------|--------------|
| A (analógico) | 129 ps | 130 ps | NO |
| B (digital)   | 93 ps | 90 ps | SÍ |

### Recomendación técnica (neutral)
- Confirmar con el equipo de FastIC+ si los timestamps son individuales por canal
  (Modelo B) o si la suma es analógica (Modelo A).
- El Modelo B es el **piso optimista**: asume timestamps independientes y promediado
  ideal. En la práctica puede haber correlaciones que reduzcan el beneficio.
- **No afirmar '83.9 ps cumple SHiP' hasta confirmar arquitectura FastIC+**.

---

## DECISIÓN 2 — Topología de SUM: co-localizado vs fila distribuida

### Contexto
En el test-beam (Constanza), el SUM4 suma SiPMs del MISMO extremo (misma luz).
En la simulación hay dos opciones:

### Opción A — Co-localizado (análogo al TB)
- Agrupa SiPMs en el mismo extremo (END_L o END_R): misma posición x, misma luz.
- END_SUM4 co-localizado: mejora timing ∝ 1/√N porque cada canal añade
  información independiente del mismo proceso físico.
- σ_END_SUM4 @ x=−690 mm (near-end): **30.2 ps** (intrínseco, walk-corr.)
- σ_TOP_nearest @ x=−690 mm: **98.8 ps**
- Limitación: lejos del END instrumentado, σ crece a >900 ps (pocos fotones).

### Opción B — Fila distribuida (caso TOP en la simulación)
- Agrupa SiPMs a lo largo de la barra (paso 20 mm): diferente posición x, diferente luz.
- El TOP_SUM4 'distribuido' en EXEC_16 mezcla llegadas de fotones de distintas
  trayectorias: cerca del haz (SiPM más próximo) y lejos (SiPMs vecinos).
- El beneficio de sumar se ve atenuado porque los canales lejanos aportan ruido.
  (Ver T6: σ_cluster_id//4=72.1 ps vs σ_dynamic_nearest=85.0 ps @ x=0)

### El `id//4` (pregunta abierta EXEC_16)
- `id//4` agrupa locales {0..3},{4..7},... → clusters fijos de 4 SiPMs adjacentes
- Para END: `id//4` en local_id da {0..3}→primer cuarteto, {4..7}→segundo cuarteto
  — co-localizados (misma face, misma luz) → topología correcta para comparar con TB.
- Para TOP: `id//4` en local_id (0..69) da clusters de 4 SiPMs a paso 20mm
  — distribuidos, no co-localizados.

### Consecuencia
| Topología | Descripción | σ @ near-end | σ @ centro |
|-----------|-------------|-------------|-----------|
| Co-localizado END | 4 del mismo extremo | 30.2 ps | 944.2 ps |
| Distribuido TOP   | 4 más próximos, paso 20mm | 98.8 ps | 107.8 ps |

### Recomendación técnica (neutral)
- Para la comparación sim↔TB de Constanza: usar **END_SUM4 co-localizado**
  (mismos SiPMs del mismo extremo, misma fotónica).
- Para el detector TOP, la suma `id//4` fija (clusters de 4 adyacentes) es mejor
  que 'nearest-4 dinámico' si los clusters se calibran por posición.
- La 'ventana móvil de 4 vecinos' es equivalente al 'nearest-4 dinámico'
  y da resultados similares a `id//4` fijo cuando el haz está bien centrado.

---

**Resumen de decisiones pendientes:**
1. FastIC+ architecture: analog sum (Model A, ~129 ps) vs digital average (Model B, ~84-91 ps)
2. TOP SUM topology: co-localized cluster (`id//4`) vs dynamic nearest-N

Estos dos ítems NO afectan los números intrínsecos de EXEC_16/17/18 —
solo la proyección a σ_total con el sistema de lectura real.