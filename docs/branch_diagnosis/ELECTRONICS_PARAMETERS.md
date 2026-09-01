# ELECTRONICS_PARAMETERS.md
# Parámetros de electrónica y pulso — hoja de referencia

Dispositivo físico: **AFBR-S4N66P024M** (Broadcom NUV-MT, 6×6 mm²)  
Punto de operación del banco: **V_OV = 10 V** (W4_Weekly_meeting_Timinig_Detector_Status.pdf)  
Electrónica de lectura: **FastIC+** (ASIC desarrollado para SHiP)

Este documento es la fuente de referencia para el capítulo de metodología.
Cada valor va acompañado de su fuente con sección y de su estado de confianza.

Estados admitidos:  
`MEDIDO` · `DATASHEET` · `PUBLICADO (otro punto de operación)` ·  
`PUBLICADO (otro dispositivo)` · `PUBLICADO (otro montaje)` · `NO MEDIDO` · `SIN VERIFICAR`

---

## Modelo de pulso SPE

| Parámetro | Símbolo | Valor | Fuente | Estado |
|-----------|---------|-------|--------|--------|
| Constante de subida — montaje Peña-Rodríguez | τ_r | 2.0 ns | arXiv:2411.16710 §4 (circuito con cancelación polo-cero) | PUBLICADO (otro montaje) |
| Constante de caída — montaje Peña-Rodríguez | τ_f | 3.0 ns | arXiv:2411.16710 §4 | PUBLICADO (otro montaje) |
| Constante de subida — SiPM intrínseco (AFBR-S4N66P024M) | τ_r | **No publicado** | — | **NO MEDIDO** |
| Constante de caída — SiPM intrínseco (AFBR-S4N66P024M) | τ_f | 55.0 ns | DS105 tabla "Optical and Electrical Features" | DATASHEET |
| Constante de subida — montaje FastIC+ en este banco | τ_r | **No medido** | — | **NO MEDIDO** |
| Constante de caída — montaje FastIC+ en este banco | τ_f | **No medido** | — | **NO MEDIDO** |

**Nota:** El par (τ_r = 2 ns, τ_f = 55 ns) del código anterior mezclaba la constante de
subida del montaje acortado de Peña-Rodríguez con la de caída intrínseca del DS105.
No describía ningún sistema físico real. Ver `ANALYSIS_CONSOLIDATION.md §6.1`.

La configuración relevante para este detector es `fastic_measured`, que está **pendiente**
de medida en banco de láser con FastIC+.

---

## SPTR (Single Photon Time Resolution)

| Parámetro | Valor | Conversión explícita | Fuente | Estado |
|-----------|-------|---------------------|--------|--------|
| SPTR intrínseco FWHM (AFBR-S4N66P014M, OV ≈ 15.5 V) | 137 ± 4 ps FWHM | — | Lee et al., IEEE TRPMS 9(4) 406–411 (2025), DOI 10.1109/TRPMS.2024.3518479 | PUBLICADO (otro punto de operación) |
| SPTR intrínseco σ (derivado de Lee) | 137 / 2.355 ≈ **58.2 ps** σ | `FWHM_TO_SIGMA = 1/2.355` (en pulse_models.py) | Ídem | PUBLICADO (otro punto de operación) |
| SPTR detector FWHM (AFBR-S4N66P014M, OV ≈ 15.5 V) | 172 ps FWHM | — | Lee et al. (2025) | PUBLICADO (otro punto de operación) |
| SPTR detector σ (derivado de Lee) | 172 / 2.355 ≈ **73.0 ps** σ | `FWHM_TO_SIGMA = 1/2.355` (en pulse_models.py) | Ídem | PUBLICADO (otro punto de operación) |
| SPTR (AFBR-S4N66P024M, OV = 10 V) | **No medido** | — | — | **NO MEDIDO** |

**Advertencia de punto de operación:** Lee midió a V_OV ≈ 15.5 V; el banco de este proyecto
opera a V_OV = 10 V. El SPTR empeora al bajar el sobrevoltaje (más carga parcial de
microcelda → mayor dispersión del tiempo de disparo). Los 58.2 / 73.0 ps de Lee son
**cotas optimistas** para el punto de operación real de este detector.

---

## Electrónica FastIC+

| Parámetro | Valor | Fuente | Estado |
|-----------|-------|--------|--------|
| TDC embebido | 25 ps LSB | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |
| Jitter cuantización TDC (cota inferior) | 25/√12 ≈ **7.2 ps** σ | Derivado de TDC LSB (distribución uniforme) | PUBLICADO (cota inferior) |
| Jitter total FastIC+ | **No medido** | — | **NO MEDIDO** |
| Interfaces de salida | 9 canales SLVS, Aurora 64B/66B | TD_34th_SHIP_Gvasquez.pdf | PUBLICADO |

**Nota sobre el valor de 30 ps:** aparece en W4_Weekly_meeting_Timinig_Detector_Status.pdf
como resolución temporal esperada de un detector completo (teja EJ-228 de 1×1×0.3 cm
con 4 SiPM FBK de 2 mm²). **No es el jitter del ASIC FastIC+.** Usarlo como
`electronics_sigma_ps` es incorrecto. Ver `ANALYSIS_CONSOLIDATION.md §6.4`.

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

*Actualizado: 2026-09-01. Ver también `ANALYSIS_CONSOLIDATION.md §6` y `SPTR_PROVENANCE.md`.*
