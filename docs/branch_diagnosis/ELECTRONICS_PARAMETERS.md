# ELECTRONICS_PARAMETERS.md
# Parámetros de electrónica y pulso — hoja de referencia

Dispositivo físico: **AFBR-S4N66P014M / 024M** (Broadcom NUV-MT, 6×6 mm², celda 40 µm)  
Punto de operación del banco: **V_OV = 10 V** (V_bias ≈ 38 V; máximo 16 V por DS105)  
Electrónica de lectura: **FastIC+** — discriminación y TDC de 25 ps dentro del ASIC;
salida SLVS / Aurora 64B/66B a 1.28 Gbps; lectura analógica externa descartada  
Fuente: `TD_34th_SHIP_Gvasquez.pdf`; W4_Weekly_meeting_Timinig_Detector_Status.pdf

Este documento es la fuente de referencia para el capítulo de metodología.
Cada valor va acompañado de su fuente con sección y de su estado de confianza.

Estados admitidos:  
`MEDIDO` · `DATASHEET` · `PUBLICADO (otro punto de operación)` ·  
`PUBLICADO (otro dispositivo)` · `PUBLICADO (otro montaje)` · `DERIVADO (banco propio)` ·
`NO MEDIDO` · `SIN VERIFICAR`

---

## Modelo de pulso SPE

| Parámetro | Símbolo | Valor | Fuente | Estado |
|-----------|---------|-------|--------|--------|
| Constante de subida — montaje Peña-Rodríguez | τ_r | 2.0 ns | arXiv:2411.16710 §4 (circuito con cancelación polo-cero) | PUBLICADO (otro montaje) |
| Constante de caída — montaje Peña-Rodríguez | τ_f | 3.0 ns | arXiv:2411.16710 §4 | PUBLICADO (otro montaje) |
| Constante de subida — SiPM intrínseco (AFBR-S4N66P024M) | τ_r | **No publicado** | — | **NO MEDIDO** |
| Constante de caída — SiPM intrínseco (AFBR-S4N66P024M) | τ_f | 55.0 ns | DS105 tabla "Optical and Electrical Features" | DATASHEET |
| Constante de subida — montaje FastIC+ en este banco | τ_r | **2.0 ns** | Capturas de osciloscopio banco (1 GHz, 6.4 GS/s); ajuste exponencial del flanco; TIA ~20 Ω + discriminador no altera flanco de subida | DERIVADO (banco propio) — confirmar formalmente |
| Constante de caída — montaje FastIC+ en este banco | τ_f | **3.0 ns** | Ídem; red de cancelación polo-cero reduce 55 ns → 3 ns | DERIVADO (banco propio) — confirmar formalmente |

**Nota:** El par (τ_r = 2 ns, τ_f = 55 ns) del código anterior mezclaba la constante de
subida del montaje acortado de Peña-Rodríguez con la de caída intrínseca del DS105.
No describía ningún sistema físico real. Ver `ANALYSIS_CONSOLIDATION.md §6.1`.

La configuración `fastic_measured` en `pulse_models.py` usa ahora (2, 3) ns, derivado de
las capturas del banco. Requiere confirmación formal midiendo τ_r/τ_f con el FastIC+ en
la cadena completa.

---

## SPTR (Single Photon Time Resolution)

| Parámetro | Valor | Conversión explícita | Fuente | Estado |
|-----------|-------|---------------------|--------|--------|
| SPTR intrínseco FWHM (AFBR-S4N66P014M, OV ≈ 15.5 V) | 137 ± 4 ps FWHM | — | Lee et al., IEEE TRPMS 9(4) 406–411 (2025), DOI 10.1109/TRPMS.2024.3518479 | PUBLICADO (otro punto de operación) |
| SPTR intrínseco σ (derivado de Lee) | 137 / 2.355 ≈ **58.2 ps** σ | `FWHM_TO_SIGMA = 1/2.355` (en pulse_models.py) | Ídem | PUBLICADO (otro punto de operación) |
| SPTR detector FWHM (AFBR-S4N66P014M, OV ≈ 15.5 V) | 172 ps FWHM | — | Lee et al. (2025) | PUBLICADO (otro punto de operación) |
| SPTR detector σ (derivado de Lee) | 172 / 2.355 ≈ **73.0 ps** σ | `FWHM_TO_SIGMA = 1/2.355` (en pulse_models.py) | Ídem | PUBLICADO (otro punto de operación) |
| SPTR σ — banco propio, cota inferior (AFBR-S4N66P014M + FastIC+, OV = 10 V) | **85 ps** RMS | Ajuste OLS σ²(N) = SPTR²/N + σ_elec²; pendiente = 7877 ps² → SPTR = 88.8 ps; extremo inferior del rango | DERIVADO (banco propio) |
| SPTR σ — banco propio, cota superior (AFBR-S4N66P014M + FastIC+, OV = 10 V) | **92 ps** RMS | Ídem, extremo superior del rango (N y σ_ToA como rangos, no valores centrales) | DERIVADO (banco propio) |
| SPTR (AFBR-S4N66P024M, OV = 10 V) — medida formal | **Pendiente** | Requiere ROOT de EOS con ToA evento a evento | **NO MEDIDO (formal)** |

**Ajuste del banco:** tres puntos (N=2 σ=63 ps; N=4 σ=43 ps; N=7.5 σ=33.5 ps, midpoints de rangos).
Regresión OLS de σ² vs 1/N: pendiente = 7877 ps² → **SPTR = 88.8 ≈ 89 ps RMS**, intercepto = −6 ps²
(→ σ_elec ≈ 0, cota superior ≤ 18 ps). El jitter de referencia del láser (SYNC ≤ 5 ps) es
despreciable (corrección < 0.5 ps). Coherencia con Lee 2025: 89 ps (10 V) > 58 ps (15.5 V),
dirección correcta (SPTR empeora al bajar OV).

La referencia de colaboración de 106 ps (Hamamatsu a 3–4.5 V OV) cae cerca de este rango —
la cota conservadora del grupo resulta razonable para este sensor a 10 V OV, y mantiene
comparabilidad con el resto de la colaboración.

---

## Electrónica FastIC+

| Parámetro | Valor | Fuente | Estado |
|-----------|-------|--------|--------|
| TDC embebido | 25 ps LSB | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |
| Jitter cuantización TDC (cota inferior) | 25/√12 ≈ **7.1 ps** σ | Derivado de TDC LSB (distribución uniforme) | PUBLICADO (cota inferior) |
| Jitter total FastIC+ (σ_elec) | ≤ **18 ps** σ | Banco propio: ajuste σ²(N), peor caso N=7.5 SPTR=89 ps | DERIVADO (banco propio) — cota superior |
| SUM4 | Suma analógica de 4 canales dentro del ASIC, antes de la discriminación | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |
| Modo de salida | Digital por canal (discriminación + TDC en el ASIC) | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |
| Interfaces de salida | 9 canales SLVS, Aurora 64B/66B a 1.28 Gbps | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |

**Nota sobre el valor de 30 ps:** aparece en W4_Weekly_meeting_Timinig_Detector_Status.pdf
como resolución temporal esperada de un detector completo (teja EJ-228 de 1×1×0.3 cm
con 4 SiPM FBK de 2 mm²), no como jitter del ASIC FastIC+. El banco propio **excluye**
σ_elec = 30 ps: con SPTR = 89 ps y σ_elec = 30 ps el modelo predice σ_ToA = 44 ps a N=7.5,
pero se observa 30–37 ps. El rango derivado es σ_elec ≤ 18 ps (o ≤ 16 ps si SPTR = 92 ps).
Ver `ANALYSIS_CONSOLIDATION.md §6.4`.

---

## Fracción dCFD

| Parámetro | Valor | Fuente | Estado |
|-----------|-------|--------|--------|
| Fracción óptima dCFD (`--fraction`) | 0.14 (14%) | Sin cita | **SIN VERIFICAR** |
| Rango óptimo para SiPMs de 3×3 mm² | 3%–6% | Cattaneo et al. (arXiv:1402.1404 §II-F) | PUBLICADO (otro dispositivo) |
| Salida de `analyze_dCFD_fraction.C` en repositorio | No encontrada | — | **SIN VERIFICAR** |

El 14% es plausiblemente mayor que el óptimo de 3–6% para 3×3 mm² (Derenzo et al.,
PMC nihms596188, predice óptimos mayores con mayor jitter de tránsito), pero no ha sido
verificado para el AFBR-S4N66P024M. Requiere correr `analyze_dCFD_fraction.C` sobre
datos de este detector.

---

## Parámetros del dispositivo (DS105)

| Parámetro | Valor | Fuente | Estado |
|-----------|-------|--------|--------|
| PDE pico | 63% @ 420 nm | DS105 | DATASHEET |
| Diámetro activo | 6×6 mm² | DS105 | DATASHEET |
| Voltaje de ruptura V_BD | 32.5 V | DS105 | DATASHEET |
| Voltaje de ruptura V_br (Peña-Rodríguez) | 45 V | arXiv:2411.16710 Tabla 1 | PUBLICADO — **DISCREPANCIA con DS105** |
| Voltaje de operación banco | V_BD + 10 V ≈ 42.5 V | W4 Weekly meeting | PUBLICADO |
| Crosstalk | 23% | DS105 | DATASHEET |

**Discrepancia V_br:** DS105 especifica V_BD = 32.5 V; Peña-Rodríguez lista V_br = 45 V
para el mismo dispositivo. La diferencia de ≈12.5 V puede deberse a convención de medida
distinta (e.g., medida eléctrica vs. fotónica) o a variación de lote. No se resuelve aquí.

---

## Referencia de la colaboración — 106 ps y ambigüedad de unidades

| Parámetro | Valor | Dispositivo | Fuente | Estado |
|-----------|-------|-------------|--------|--------|
| SPTR σ referencia Geant4 | 106 ps | Hamamatsu S13360-3050CS, 3–4.5 V OV | Estándar colaboración SHiP | PUBLICADO (otro dispositivo) |
| SPTR FWHM colaboración | "100–120 ps FWHM" | Ídem | Notas del grupo | **AMBIGÜEDAD — ver nota** |

**Advertencia de unidades (§1.4 del plan de cierre):** las notas del grupo describen
"100–120 ps FWHM" y "~106 ps RMS" como si fueran equivalentes. **Son incompatibles:**
110 ps FWHM ≡ ≈47 ps σ, no 106 ps. Si el valor inyectado en Geant4 (106 ps) se interpreta
como sigma cuando la medida original era FWHM, el jitter está sobreestimado en un factor
≈2.355. No se resuelve aquí. Requiere consulta con Gerardo para identificar la convención
original. No usar este valor en análisis propios sin aclarar la ambigüedad.

**Coherencia con el banco propio:** 106 ps cae cerca del rango 85–92 ps derivado del banco.
La cota conservadora del grupo resulta razonable para el Broadcom a 10 V OV, y su uso
mantiene comparabilidad con el resto de la colaboración — siempre que la ambigüedad de
unidades se resuelva y el valor sea σ, no FWHM.

---

*Actualizado: 2026-09-01. Ver también `ANALYSIS_CONSOLIDATION.md §6` y `SPTR_PROVENANCE.md`.*
