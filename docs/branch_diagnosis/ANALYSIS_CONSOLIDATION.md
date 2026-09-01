# ANALYSIS_CONSOLIDATION.md
# Consolidación de feature/sipm-electronics-response → main

Fecha: 2026-09-01  
Auditoría completa: `docs/branch_diagnosis/ANALYSIS_CONSOLIDATION_AUDIT.md`

---

## 1. Qué se migró

### Archivos nuevos (A) — traídos desde `feature/sipm-electronics-response`

**`analysis/timing/` (creado en esta consolidación):**

| Archivo | Función |
|---------|---------|
| `sipm_waveform_dcfd.py` | Waveform NUV-MT + dCFD 14% (Lv/Peña-Rodríguez) |
| `sipm_waveform_dcfd.cpp` | Mismo algoritmo en C++ |
| `analyze_dCFD.py` | dCFD directo sobre timestamps |
| `analyze_dCFD_fraction.C` | Barrido de fracción dCFD |
| `analyze_dCFD_5thPhoton.C` | Trigger por 5.º fotón |
| `analyze_basic.py` | FPT básico (mínimo de hit times) |
| `resolution_vs_x_dCFD.py` | σ_t vs x con dCFD + jitter FastIC+ |
| `resolution_vs_x_FPT.py` | σ_t vs x con FPT |
| `fpt_vs_n_profile.C` | Perfil ⟨t_n⟩ vs n por SiPM |
| `fpt_vs_n_profile_batch.C` | Versión batch |
| `fpt_vs_n_profile_batch_slides.C` | Pipeline Beamer (ex raíz del repo) |
| `fpt_manifest_to_beamer.py` | TSV → Beamer .tex (ex raíz) |
| `png_to_beamer.py` | PNG → Beamer (ex raíz) |
| `TimeMarkScan.C` | Barrido de modos de time mark |
| `SiPMRankingScan_v2.C` | Ranking SiPMs por FPT (ajuste gaussiano) |
| `SiPMRankingScan_RMS.C` | Ranking SiPMs por RMS |
| `SiPMRankingScan_coreSigma.C` | Ranking SiPMs por σ del núcleo |

**`analysis/` (raíz):**
- `convert_vis_exports.py` — convertidor PDF→PNG para exports de visualización

**`presentations/`:**
- `presentacion_timing_detector_scan_ej-230.tex` — presentación del scan de timing

**`macros/`:**
- `vis_scan_png.mac`, `vis_scan_png_lite.mac`, `vis_scan_png_step.mac` — macros de visualización del scan

### Cambios aplicados a archivos existentes (M neutros)

| Archivo | Cambio |
|---------|--------|
| `analysis/merge_runs.py` | Actualización de nombres de scripts en ayuda |
| `analysis/ResolutionScan_v2.C` | `drawFits` param, diagnóstico de entradas, fallback RMS, `TLatex` |
| `analysis/topreadout_crosstalk.py` | `.to_numpy()` para compat. matplotlib/pandas ≥ 2.0 |
| `analysis/compare_edge_wraps.py` | `.to_numpy()` (5 llamadas) |
| `analysis/edge_resolution.py` | `expected_x`, tolerancia a scan sin ROOT, `.to_numpy()` |
| `analysis/grouped_resolution.py` | `.to_numpy()` en `ax.errorbar` |

---

## 2. Qué NO se migró y por qué

| Elemento | Razón |
|---------|-------|
| `analysis/sipm_waveform_dcfd` (binario, 83 KB) | Compilado para la arquitectura de W2; no ejecutable en general |
| `.claude/settings.json` | Sobreescribiría la configuración de Claude Code en esta máquina |
| `.codex` | Vacío; no tiene valor |
| `audit/runs_file_inventory.txt` | Snapshot de W2 (2026-06-09); no pertenece al árbol de análisis |
| Raíz: `fpt_vs_n_profile.C`, `fpt_vs_n_profile_batch.C` | Duplicados exactos (MD5 idéntico) de `analysis/timing/` — canónico en subdir |

---

## 3. Cambios de geometría bloqueados

Los cuatro cambios siguientes de la rama `feature/sipm-electronics-response` asumen
`kNTopSiPMs = 20` (geometría de esa rama). Main tiene `kNTopSiPMs = 70`.
Aplicarlos descartaría en silencio los canales 36–85 sobre cualquier dataset de 86
canales. **No se aplicaron.**

| Archivo | Cambio bloqueado |
|---------|-----------------|
| `analysis/resolution_vs_x_fixed.py:34` | `N_TOP_SIPMS: 70 → 20` |
| `analysis/grouped_resolution.py:61` | `TOP_IDS: range(16, 86) → range(16, 36)` |
| `analysis/grouped_resolution.C` | Loop dinámico 16–85 → listas estáticas 16–35 |
| `analysis/analyze_hits.C` | `hOccupancy: 86 bins → 60 bins` |

Propuesta de derivación en runtime documentada en `ANALYSIS_CONSOLIDATION_AUDIT.md §2.5b`.

---

## 4. Modelos de pulso — incompatibilidad no resuelta

Los dos marcos de análisis coexisten pero usan parámetros de pulso incompatibles:

| Marco | τ_rise | τ_fall | Dispositivo | Fuente |
|-------|--------|--------|-------------|--------|
| `exec14/engine.py` | 0.5 ns | 5.0 ns | SiPM de ej200_endonly (no identificado) | `common.py` (ej200_endonly exec07) |
| `timing/sipm_waveform_dcfd.py` | 2.0 ns | 55.0 ns | AFBR-S4N66P024M (Broadcom NUV-MT) | τ_fall = 55 ns en DS105; τ_rise sin cita |

El SiPM físico del detector es el **AFBR-S4N66P024M**. Para este dispositivo, DS105
confirma τ_fall = 55 ns. El τ_rise = 0.5 ns de `engine.py` proviene de un modelo de
SiPM diferente. Que los marcos no se solapen funcionalmente no resuelve cuál es el
modelo correcto para este detector.

**Estado:** decisión pendiente del usuario. Los marcos no deben combinarse entre sí
hasta resolver cuál corresponde al dispositivo real.

---

## 5. Datos disponibles vs geometría

Los datasets de Fase 7 en `results/` tienen las siguientes configuraciones:

| Dataset | Config | global_id | Top SiPMs |
|---------|--------|-----------|-----------|
| `scan_end_vikuiti/` | END-only | 0–15 | 0 (sin top) |
| `scan_end_vik_sparse_top_v2/` | Sparse top | 0–19 | 4 (IDs 16–19) |

Los datasets de exec13 (junio 2026, óptica anterior) sí tienen 86 canales con
global_id hasta 85 (70 top SiPMs). Son de una configuración óptica diferente a
Fase 7 y no deben mezclarse con los datasets de Fase 7.

**Consecuencia:** Los análisis que operan sobre top SiPMs (IDs 16+) no tienen datos
de Fase 7 con la geometría completa de 70 canales top. Los únicos datos de Fase 7
disponibles hoy tienen:
- END: resolución completa (16 SiPMs, IDs 0–15)
- TOP: solo 4 SiPMs (IDs 16–19) en `scan_end_vik_sparse_top_v2`

Esto determina qué análisis se pueden correr hoy sobre datos de Fase 7 sin simular:
scripts de END solo (`exec07/`, `exec14/engine.py`, `timing/analyze_dCFD.py` cara 0/1)
y scripts que toleran pocos top SiPMs. Los scripts de ranking de top SiPMs
(`SiPMRankingScan_*.C`, `grouped_resolution.*` con agrupaciones top) requieren una
simulación nueva con la geometría completa de 70 canales.

---

## 6. Parámetros de pulso y electrónica — correcciones aplicadas (2026-09-01)

### 6.1 El error de (2 ns, 55 ns): mezcla de configuraciones

El código original de `sipm_waveform_dcfd.py` usaba `τ_r = 2.0 ns` y `τ_f = 55.0 ns`.
Este par **no describe ninguna configuración física real**. Peña-Rodríguez
(arXiv:2411.16710) caracteriza el AFBR-S4N66P024M y documenta **dos** configuraciones:

| Configuración | τ_r | τ_f | Qué describe |
|---|---|---|---|
| SiPM intrínseco | No medido | 55 ns | Recarga de microcelda (τ_F = R_Q·C_J); DS105 |
| Montaje acortado (§4) | 2 ns | 3 ns | Amplificador RF clase A con cancelación polo-cero |
| Código anterior | 2 ns | 55 ns | **Ninguna de las dos** |

**Origen del error:** el Apéndice A de arXiv:2411.16710 contiene un ejemplo de código con
`Rt = 2e-9` y `Ft = 50e-9`. Alguien tomó el `Rt` del ejemplo (que corresponde al montaje
acortado de §4) y lo combinó con el `τ_f = 55 ns` del datasheet DS105 (valor intrínseco),
sin advertir que los dos nanosegundos están físicamente atados al circuito de acortamiento.

### 6.2 Configuraciones de pulso ahora disponibles

Definidas en `analysis/timing/pulse_models.py` (única fuente de verdad):

| Nombre | τ_r | τ_f | Fuente | Estado |
|--------|-----|-----|--------|--------|
| `penarodriguez_shortened` | 2.0 ns | 3.0 ns | arXiv:2411.16710 §4 | PUBLICADO (otro montaje) |
| `broadcom_intrinsic` | No medido | 55.0 ns | DS105 | DATASHEET (τ_r NO MEDIDO) |
| `fastic_measured` | No medido | No medido | PENDIENTE | NO MEDIDO |

`sipm_waveform_dcfd.py` requiere ahora `--pulse-model` explícito y aborta si no se da.
Los overrides `--tau-r-ns` y `--tau-f-ns` siguen existiendo pero emiten aviso de
salida de configuración validada. Combinar τ de modelos distintos es físicamente incorrecto.

### 6.3 SPTR — parámetros actualizados

| Parámetro CLI | Antes | Ahora | Fuente | Estado |
|--------------|-------|-------|--------|--------|
| `--transit-sigma-ps` | 200 ps (default) | Sin default; requerido | — | — |
| σ recomendado (intrínseco) | — | 137/2.355 ≈ 58.2 ps | Lee et al. IEEE TRPMS 2025 | PUBLICADO (otro OV) |
| σ recomendado (detector) | — | 172/2.355 ≈ 73.0 ps | Lee et al. IEEE TRPMS 2025 | PUBLICADO (otro OV) |

**Aviso de sobrevoltaje (impreso en cada ejecución):** Lee et al. midieron a V_OV ≈ 15.5 V.
Este banco opera a V_OV = 10 V (W4 Weekly meeting). El SPTR empeora al bajar el OV:
los 58.2/73.0 ps son cotas optimistas para el punto de operación real.

La conversión FWHM → σ está explícita en el código: `FWHM_TO_SIGMA = 1/2.355` en
`pulse_models.py`. Ningún σ derivado de un FWHM aparece precalculado.

### 6.4 Jitter electrónico — corrección del valor de 30 ps

`analyze_dCFD.py` usaba `electronics_sigma_ns = 30 ps`. Este valor proviene de
W4_Weekly_meeting_Timinig_Detector_Status.pdf, donde aparece como la **resolución temporal
esperada de un detector completo** (teja de EJ-228 de 1×1×0.3 cm con 4 SiPM FBK de 2 mm²),
no como el jitter del ASIC FastIC+. Es incorrecto usarlo como jitter de electrónica.

Lo único documentado del FastIC+ es el TDC de 25 ps (TD_34th_SHIP_Gvasquez.pdf), cuya
cuantización contribuye ≥ 25/√12 ≈ 7.2 ps. Eso es cota inferior, no jitter total.

`--electronics-sigma-ps` no tiene default en el código actualizado. El usuario debe
especificarlo explícitamente; usar 0 para deshabilitar.

### 6.5 Fracción dCFD

El 14% se conserva como default pero se documenta como ÓPTIMO NO VERIFICADO. No existe
salida de `analyze_dCFD_fraction.C` en este repositorio que lo justifique. Cattaneo et al.
(arXiv:1402.1404) reporta 3–6% para SiPMs de 3×3 mm²; Derenzo et al. (PMC nihms596188)
predice óptimos mayores con mayor jitter (plausible para 6×6 mm², no verificado).

### 6.6 Discrepancia V_br

Peña-Rodríguez (arXiv:2411.16710 Tabla 1) lista V_br = 45 V para el AFBR-S4N66P024M.
DS105 especifica V_BD = 32.5 V. La discrepancia (≈12.5 V) puede deberse a convención de
medida distinta. No se resuelve aquí; se registra. Ver `ELECTRONICS_PARAMETERS.md`.

---

*Fin del documento de consolidación.*
