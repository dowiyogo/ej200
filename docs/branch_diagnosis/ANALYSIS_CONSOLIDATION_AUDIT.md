# §2 AUDITORÍA — feature/sipm-electronics-response → main
# CHECKPOINT OBLIGATORIO — no migrar nada hasta aprobación

Fecha: 2026-09-01  
Rama auditada: `origin/feature/sipm-electronics-response`  
Commit tip: `git show origin/feature/sipm-electronics-response --oneline -1`

---

## Inventario de archivos

### Únicos en la rama (A) — candidatos a migración

**analysis/ (análisis):**

| Archivo | Tipo | TTree compat. |
|---------|------|---------------|
| `analyze_basic.py` | Python / FPT básico | COMPLETA |
| `analyze_dCFD.py` | Python / dCFD | PARCIAL ¹ |
| `analyze_dCFD_5thPhoton.C` | ROOT macro / 5.º fotón | COMPLETA |
| `analyze_dCFD_fraction.C` | ROOT macro / barrido de fracción | COMPLETA |
| `convert_vis_exports.py` | Python / PDF→PNG | N/A (no lee ROOT) |
| `fpt_vs_n_profile.C` | ROOT macro / perfil t_n por SiPM | COMPLETA |
| `fpt_vs_n_profile_batch.C` | ROOT macro / batch | COMPLETA |
| `resolution_vs_x_dCFD.py` | Python / σ_t vs x (dCFD) | COMPLETA |
| `resolution_vs_x_FPT.py` | Python / σ_t vs x (FPT) | COMPLETA |
| `SiPMRankingScan_v2.C` | ROOT macro | COMPLETA |
| `SiPMRankingScan_RMS.C` | ROOT macro | COMPLETA |
| `SiPMRankingScan_coreSigma.C` | ROOT macro | COMPLETA |
| `sipm_waveform_dcfd.py` | Python / reconstrucción waveform | PARCIAL ² |
| `sipm_waveform_dcfd.cpp` | C++ / mismo algoritmo | PARCIAL ² |
| `sipm_waveform_dcfd` | **BINARIO (83 KB)** | **NO MIGRAR** |
| `TimeMarkScan.C` | ROOT macro | COMPLETA |
| `topreadout_crosstalk.py` | Python / ocupación TOP | COMPLETA |

**Raíz del repo (A):**

| Archivo | Disposición propuesta |
|---------|----------------------|
| `fpt_vs_n_profile.C` | DUPLICADO de `analysis/` — no migrar |
| `fpt_vs_n_profile_batch.C` | DUPLICADO de `analysis/` — no migrar |
| `fpt_vs_n_profile_batch_slides.C` | ÚNICO — mover a `analysis/` |
| `fpt_manifest_to_beamer.py` | ÚNICO — mover a `analysis/` |
| `png_to_beamer.py` | ÚNICO — mover a `analysis/` |
| `presentacion_timing_detector_scan_ej-230.tex` | Presentación — mover a `presentations/` |

**Otros (A) — NO migrar:**

| Archivo | Razón |
|---------|-------|
| `.claude/settings.json` | Sobreescribiría config de Claude Code actual |
| `.codex` | Vacío / trivial |
| `audit/runs_file_inventory.txt` | Snapshot informativo de W2 (2026-06-09); no pertenece a analysis/ |
| `macros/vis_scan_png.mac` | Macro de Geant4 — mover a `macros/` directamente |
| `macros/vis_scan_png_lite.mac` | ídem |
| `macros/vis_scan_png_step.mac` | ídem |

### Modificados en la rama (M) — requieren decisión de merge

**Análisis (contenido oportuno):**

| Archivo | Cambio resumido |
|---------|----------------|
| `analysis/merge_runs.py` | 2 líneas: nombre de script actualizado (`analyze.py`→`analyze_dCFD.py`). Trivial; aplicar. |
| `analysis/ResolutionScan_v2.C` | +~90 líneas: `drawFits` param, diagnóstico de entradas, histogramas de ajuste guardados, fallback RMS si el ajuste falla. Mejora sustancial; aplicar. |
| `analysis/topreadout_crosstalk.py` | Cambios pendientes de diff completo. |
| `analysis/edge_resolution.py` | Cambios pendientes de diff completo. **Depende de `resolution_vs_x_fixed.py` (importación explícita).** |
| `analysis/grouped_resolution.C` | Cambios pendientes de diff completo. |
| `analysis/grouped_resolution.py` | Cambios pendientes de diff completo. |
| `analysis/compare_edge_wraps.py` | Cambios pendientes de diff completo. |
| `analysis/analyze_hits.C` | Cambios pendientes de diff completo. |
| `analysis/resolution_vs_x_fixed.py` | Cambios pendientes de diff completo. **Requerido por `edge_resolution.py`.** |

**Código fuente (M) — FUERA DE ALCANCE (no migrar):**

30 archivos `src/`, `include/`, `main.cc`, `CMakeLists.txt`, `macros/*.mac` divergen del main actual. Esto indica que la rama tiene una versión de simulación distinta a Fase 7. Los scripts de análisis fueron desarrollados contra esa versión.

---

## §2.1 Compatibilidad TTree — Fase 7

Ramas de Fase 7 (`sipm_hits`): `event_id`, `global_id`, `local_id`, `face_type`, `time_ns`, `x_mm`, `y_mm`, `z_mm`, `gun_x_mm`, `energy_eV`, `wl_nm`, `pde`.

### Ramas consumidas y veredicto por script

| Script | Ramas leídas del TTree | Ausentes en Fase 7 | Veredicto |
|--------|------------------------|-------------------|-----------|
| `sipm_waveform_dcfd.py/cpp` | `event_id`, `face_type`, `time_raw_ns`→`time_ns` (fallback), `gun_x_mm` | `time_raw_ns` (fallback OK) | PARCIAL — funcional |
| `analyze_dCFD.py` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `time_raw_ns`(opt), `sptr_jitter_ns`(opt→0.0), `energy_eV`, `wl_nm`, `pde`, `x_mm`, `y_mm`, `z_mm`, `gun_x_mm`(opt), `sensor_model_id`(opt), `sensor_sptr_sigma_ns`(opt), `electronics_sigma_ns`(opt) | `time_raw_ns`, `sptr_jitter_ns`, `sensor_*` — todos opcionales con fallback graceful | PARCIAL — funcional; columnas derivadas (`sptr_jitter_ns=0`, `time_raw_ns=time_ns`) |
| `analyze_basic.py` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `energy_eV`, `wl_nm`, `pde`, `x_mm`, `y_mm`, `z_mm`, `gun_x_mm`(opt) | ninguna | COMPLETA |
| `resolution_vs_x_dCFD.py` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `energy_eV`, `wl_nm`, `pde`, `x_mm`, `y_mm`, `z_mm`, `gun_x_mm` | ninguna | COMPLETA |
| `resolution_vs_x_FPT.py` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `energy_eV`, `wl_nm`, `pde`, `x_mm`, `y_mm`, `z_mm`, `gun_x_mm` | ninguna | COMPLETA |
| `analyze_dCFD_fraction.C` | `event_id`, `face_type`, `time_ns` | ninguna | COMPLETA |
| `analyze_dCFD_5thPhoton.C` | `event_id`, `face_type`, `time_ns` | ninguna | COMPLETA |
| `fpt_vs_n_profile.C` (ambas copias) | `event_id`, `global_id`, `local_id`, `face_type`, `time_ns`, `x_mm`, `y_mm`, `gun_x_mm` | ninguna | COMPLETA |
| `fpt_vs_n_profile_batch.C` (ambas) | `gun_x_mm` | ninguna | COMPLETA |
| `fpt_vs_n_profile_batch_slides.C` | `gun_x_mm` | ninguna | COMPLETA |
| `SiPMRankingScan_v2.C` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `SiPMRankingScan_RMS.C` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `SiPMRankingScan_coreSigma.C` | `event_id`, `face_type`, `global_id`, `local_id`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `TimeMarkScan.C` | `event_id`, `face_type`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `edge_resolution.C` | `event_id`, `face_type`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `edge_resolution.py` | vía `resolution_vs_x_fixed.load_event_level_data` (`event_id`, `face_type`, `time_ns`, `gun_x_mm`) | ninguna (hereda compat. de `resolution_vs_x_fixed.py`) | COMPLETA |
| `grouped_resolution.C` | `event_id`, `global_id`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `grouped_resolution.py` | `event_id`, `global_id`, `time_ns`, `gun_x_mm` | ninguna | COMPLETA |
| `topreadout_crosstalk.py` | `event_id`, `face_type`, `global_id` | ninguna | COMPLETA |
| `compare_edge_wraps.py` | lee `edge_summary.csv` (CSV, no ROOT) | N/A | N/A |
| `convert_vis_exports.py` | no lee ROOT | N/A | N/A |

**Nota:** `local_id` no es leído directamente por `topreadout_crosstalk.py`; lo computa como `global_id - TOP_GLOBAL_OFFSET`.

---

## §2.2 Parámetros físicos hardcodeados

### `sipm_waveform_dcfd.py` / `sipm_waveform_dcfd.cpp`

| Parámetro | Valor default | Procedencia en código | Contraste con SPTR_PROVENANCE.md |
|-----------|--------------|----------------------|----------------------------------|
| `--transit-sigma-ps` | **200.0 ps** | `help`: "Sigma del tiempo de transito de **Lv et al.** en ps." — citado explícitamente en el texto de ayuda | Lee 2025: 58 ps (intrínseco), 73 ps (detector). **PROCEDENCIA IDENTIFICADA — valor del SiPM de Lv, no del AFBR-S4N66P024M.** |
| `--fraction` | 0.14 (14%) | NINGUNA — sin cita en código | — (parámetro de dCFD, no de SPTR) |
| `--tau-r-ns` | 2.0 ns | `help`: "Constante de subida Broadcom" — sin cita de sección ni documento | **PROCEDENCIA PARCIAL** — atribuido a Broadcom; no aparece en DS105 como constante de tiempo; puede ser medición externa no citada |
| `--tau-f-ns` | 55.0 ns | `help`: "Constante de bajada Broadcom" — sin cita de sección | **PROCEDENCIA IDENTIFICADA** — DS105 tabla "Optical and Electrical Features" especifica τ_fall = 55 ns (tiempo de recarga de la celda) |
| `--electronics-sigma-ps` | 0.0 ps | sin cita | Contrasta con `analyze_dCFD.py` y `resolution_vs_x_dCFD.py` que usan **30 ps FastIC+** |
| `--pulse-sigma-pe` | 0.1 pe | sin cita | — |
| `--tail-window-factor` | 8.0 | sin cita | — |
| `--seed` | 12345 | convención | — |

**Aclaración sobre `--transit-sigma-ps = 200 ps` (corrección C2):**  
El código nombra explícitamente a Lv et al. (2026) como fuente en el texto de ayuda. Lv et al. modelan el tiempo de tránsito de su SiPM como gaussiana con σ = 200 ps, citando la SPTR medida de su propio dispositivo. El error relevante no es de procedencia desconocida sino de **dispositivo equivocado**: 200 ps es el SPTR del SiPM de Lv, no del Broadcom AFBR-S4N66P024M de 6×6 mm², para el cual Lee et al. 2025 mide 137±4 ps FWHM (σ ≈ 58 ps intrínseco). El script, usado con el default, sobreestima el jitter de tránsito en un factor ≈ 3.4 respecto a Lee 2025. **No se corrige en el código; se documenta.**

### `analyze_dCFD.py`

| Parámetro | Valor | Procedencia |
|-----------|-------|-------------|
| `fraction` | 0.14 | NINGUNA — sin cita |
| `electronics_sigma_ps` | **30.0 ps** | FastIC+ (referenciado en `resolution_vs_x_dCFD.py` docstring, sin cita de paper) |
| `N_TOP_SIPMS` | 20 | "debe coincidir con DetectorConstruction" — no se verifica en tiempo de ejecución |
| `BAR_HALF_X` | 700 mm | ídem |
| `TOP_POSITIONS_MM` | `linspace(-665, 665, 20)` | ídem |

### `analyze_dCFD_fraction.C`

| Parámetro | Valor default | Procedencia |
|-----------|--------------|-------------|
| `fraction` | 0.14 | NINGUNA |
| `electronicsSigmaPs` | 30.0 | NINGUNA (misma que `analyze_dCFD.py`) |

### `resolution_vs_x_dCFD.py`

| Parámetro | Valor default | Procedencia |
|-----------|--------------|-------------|
| `fraction` | 0.14 | NINGUNA |
| `electronics_sigma_ps` | **30.0 ps** | "FastIC+ gaussiano de 30 ps RMS" — docstring, sin cita bibliográfica |

### Scripts ROOT (C++) de análisis directo de TTree

`TimeMarkScan.C`, `SiPMRankingScan_v2/RMS/coreSigma.C`, `grouped_resolution.C`, `edge_resolution.C`, `fpt_vs_n_profile.C`, `fpt_vs_n_profile_batch.C`: sin parámetros físicos hardcodeados relevantes; implementan estimadores estadísticos (FPT, RMS, ajuste gaussiano al núcleo) sin asumir σ_SPTR.

---

## §2.3 Duplicados raíz vs analysis/

| Raíz | analysis/ | ¿Idénticos? | Canónico |
|------|-----------|-------------|----------|
| `fpt_vs_n_profile.C` | `analysis/fpt_vs_n_profile.C` | **SÍ** (MD5 idéntico) | `analysis/` |
| `fpt_vs_n_profile_batch.C` | `analysis/fpt_vs_n_profile_batch.C` | **SÍ** (MD5 idéntico) | `analysis/` |
| `fpt_vs_n_profile_batch_slides.C` | sin equivalente en `analysis/` | — | único en raíz |

Los duplicados de raíz NO deben migrarse. Las copias de `analysis/` son las canónicas.

---

## §2.4 Comparación con `analysis/exec14/engine.py`

**Corrección:** `analysis/exec14/engine.py` **sí existe en main** (merge `87795b4`). La búsqueda anterior se hizo en la rama, no en el destino. Comparación:

### `engine.py` (main, `analysis/exec14/`)

| Aspecto | Descripción |
|---------|-------------|
| Propósito | Motor compartido del pipeline EXEC_14; biblioteca importada por `f1_sidecar.py`…`f7_sidecar.py` |
| Estimador de tiempo | SUM4 leading-edge discriminator — superpone pulsos SPR y detecta cruce de 4 pe |
| Modelo de pulso SPE | `v(t) = (exp(-t/τ_fall) − exp(-t/τ_rise))` con `τ_rise = 0.5 ns`, `τ_fall = 5.0 ns` (de `common.py`, ej200_endonly) |
| Umbral | 4 pe (leading edge) |
| Scope | END SiPMs (`global_id < 16`) en `load_end_hits`; ALL hits disponibles vía `load_all_hits` |
| Dependencias externas | `common.py` en ruta hardcodeada `/home/reriosto/SHiP/ej200_endonly/analysis/exec07` (solo funciona en W1); PyROOT para `TH1::Fit` |
| Jitter | `JITTER_PER_HIT_NS = 20 ps` (ya en `time_ns` desde SiPMSD.cc) + `READOUT_JITTER_QUADRATURE_PS = 20 ps` (en sidecars); riesgo de doble conteo documentado como `JITTER_DOUBLE_COUNT_RISK` |
| Arquitectura | Biblioteca Python; no corre solo |

### `sipm_waveform_dcfd.py` (feature/sipm-electronics-response)

| Aspecto | Descripción |
|---------|-------------|
| Propósito | Script standalone de reconstrucción de waveform + dCFD |
| Estimador de tiempo | dCFD al 14% sobre waveform reconstruida evento por evento |
| Modelo de pulso SPE | `v(t) = A*(1−exp(−t/τ_r))*exp(−t/τ_f)` con `τ_r = 2 ns`, `τ_f = 55 ns` (Broadcom NUV-MT) |
| Umbral | 14% de la amplitud pico (dCFD) |
| Scope | Cara configurable (`--face`); por defecto END Left |
| Dependencias | uproot, numpy, scipy, pandas — sin ROOT |
| Jitter | `transit_sigma_ps = 200 ps` por fotón (Lv et al., dispositivo distinto) |
| Arquitectura | Script autónomo; produce CSV y PNG por sí solo |

### Veredicto de solapamiento

Los dos marcos **NO se reemplazan mutuamente**:

| Diferencia clave | `engine.py` | `sipm_waveform_dcfd.py` |
|-----------------|-------------|------------------------|
| Modelo de pulso | τ_r=0.5 ns, τ_f=5.0 ns (SPR de ej200_endonly) | τ_r=2 ns, τ_f=55 ns (NUV-MT Broadcom) |
| Discriminador | Leading-edge @ 4 pe | dCFD @ 14% |
| Jitter de tránsito | No modelado (ya en `time_ns`) | 200 ps σ por fotón (Lv, dispositivo distinto) |
| Portabilidad | Solo W1 (common.py hardcodeado) | Standalone |
| Rol | Motor de pipeline EXEC_14 | Herramienta exploratoria de waveform |

**Los pulsos son modelos diferentes para dispositivos distintos** (SPR ej200_endonly vs NUV-MT Broadcom). Deben coexistir, documentando explícitamente que usan parámetros de pulso incompatibles entre sí. **No fusionar ni hacer que uno llame al otro.**

---

## §2.5 Diffs de archivos modificados (M) — completados

### `analysis/resolution_vs_x_fixed.py` (M) — 1 línea ⚠ GEOMETRÍA

```diff
-N_TOP_SIPMS = 70
+N_TOP_SIPMS = 20
```

**⚠ NO APLICAR como diff directo — requiere decisión de diseño (ver §2.5b).**

El valor 70 es correcto para main: `DetectorConstruction.hh` (main) define `kNTopSiPMs = 70`, IDs 16–85. exec13 tiene datos con global_id hasta 85 (CSV confirma 70 top SiPMs). El valor 20 corresponde a la geometría de la rama, no de main. Aplicar este cambio descartaría en silencio los canales 36–85 al procesar cualquier dataset de 86 canales.

### `analysis/edge_resolution.py` (M) — ~50 líneas

Cuatro categorías de cambios:
1. **`resolution_vs_x`:** añade `columns=[...]` explícito al `pd.DataFrame(rows, ...)` → evita DataFrame vacío sin columnas cuando hay cero filas. Correctivo neutro.
2. **`dead_fraction_by_x`:** añade parámetro `expected_x`; si se pasa, itera sobre todas las posiciones esperadas aunque no haya datos. Mejora de completitud.
3. **`main`:** añade `--x-min`, `--x-max`, `--x-step`; tolera scan sin archivos ROOT. Mejora de robustez.
4. **Plots:** `.to_numpy()` explícito — compatibilidad pandas ≥ 2.0.

**Veredicto:** Los cambios 1–4 son neutros. No contienen geometría hardcodeada propia. **APLICAR (los cambios no-geométricos).** La importación `from resolution_vs_x_fixed import load_event_level_data` no usa `N_TOP_SIPMS`, por lo que `edge_resolution.py` no se ve afectado por la constante.

**Veredicto:** Mejoras de robustez y compatibilidad. Depende de `resolution_vs_x_fixed.py` con `N_TOP_SIPMS=20`. **APLICAR (junto con `resolution_vs_x_fixed.py`).**

### `analysis/topreadout_crosstalk.py` (M) — 1 línea

```diff
-ax.bar(occ["relative_index"], occ["hits_per_event"], width=0.8,
+ax.bar(occ["relative_index"].to_numpy(), occ["hits_per_event"].to_numpy(), width=0.8,
```

**Veredicto:** `.to_numpy()` para compatibilidad matplotlib/pandas ≥ 2.0. **APLICAR.**

### `analysis/grouped_resolution.C` (M) ⚠ GEOMETRÍA

Elimina el loop dinámico `range(16, 86)` (70 SiPMs) y lo reemplaza por IDs estáticos `16..35` (20 SiPMs). Afecta agrupaciones `single_top`, `sum4_top`, `full_top`.

**⚠ NO APLICAR como diff directo.** En main, `TOP_IDS = range(16, 86)` es correcto (70 SiPMs, confirmado en exec13 CSV con global_id hasta 85). Ver §2.5b.

### `analysis/grouped_resolution.py` (M) ⚠ GEOMETRÍA

`TOP_IDS = range(16, 86)` → `range(16, 36)`; también `.to_numpy()` en `ax.errorbar`.

**⚠ NO APLICAR el cambio de TOP_IDS.** El `.to_numpy()` sí es neutro y puede aplicarse de forma separada. Ver §2.5b.

### `analysis/compare_edge_wraps.py` (M) — 5 líneas

`.to_numpy()` en 5 llamadas a `ax.plot` y `ax.errorbar`.

**Veredicto:** Compatibilidad matplotlib/pandas ≥ 2.0. **APLICAR.**

### `analysis/analyze_hits.C` (M) — 2 líneas ⚠ GEOMETRÍA

```diff
-TH1D *hOccupancy = new TH1D("hOccupancy", ..., 86, 0, 86);
+TH1D *hOccupancy = new TH1D("hOccupancy", ..., 60, 0, 60);
```

(+ eliminación de `\n` final del archivo)

**⚠ NO APLICAR el cambio de 86→60.** En main la geometría completa tiene 86 canales (IDs 0–85). El histograma con 86 bins cubre exactamente todos los canales. El valor 60 no tiene justificación en el código ni en el diff; silenciaría los canales 60–85. El EOF fix (`\n`) sí es neutro pero no justifica aplicar el cambio completo. Ver §2.5b.

### §2.5b — Verificación de geometría y propuesta de derivación en tiempo de ejecución

**Verificación realizada con archivo:línea:**

| Fuente | N_TOP | N_END | Rango global_id top | Verificado en |
|--------|-------|-------|---------------------|---------------|
| `include/DetectorConstruction.hh:38` (main) | **70** | 8/lado | 16–85 | `kNTopSiPMs = 70` |
| `analysis/resolution_vs_x_fixed.py:34` (main) | **70** | 8 | — | `N_TOP_SIPMS = 70` |
| `analysis/grouped_resolution.py:61` (main) | **70** | — | range(16, 86) | `TOP_IDS = tuple(range(16, 86))` |
| `analysis/exec13/exec13_f3_npe_profile.csv` | **70** | — | hasta 85 | CSV confirma global_id=85 en datos reales |
| `results/scan_end_vikuiti/` (Fase 7 END-only) | 0 | 8/lado | 0–15 | uproot: max global_id=15 |
| `results/scan_end_vik_sparse_top_v2/` (Fase 7 sparse) | **4** | 8/lado | 0–19 | uproot: max global_id=19 |
| Feature branch `DetectorConstruction.hh` | **20** | — | 16–35 | modifica src/ — no migrar |

**Conclusión:** Los cambios de geometría (70→20, range(16,86)→range(16,36), 86→60) ajustan los scripts a la geometría de la rama `feature/sipm-electronics-response`, que tiene una configuración diferente a main (kNTopSiPMs=20 vs 70). Aplicarlos en main **silenciaría los canales 36–85 sin error**, produciendo resultados correctos solo sobre un tercio de los datos de top en cualquier dataset de 86 canales.

Nota adicional: los datos de Fase 7 en `results/` son END-only (0 top) o sparse (4 top), así que los scripts con N_TOP=20 funcionarían sobre esos datos sin pérdida observable — pero fallarían silenciosamente sobre cualquier dataset de 86 canales (incluyendo los de exec13 ya en main).

**Propuesta de derivación en tiempo de ejecución (no ejecutar aún — requiere aprobación):**

En lugar de una constante hardcodeada, los scripts deberían derivar la geometría del TTree o de un módulo compartido:

```python
# Opción A — derivar del TTree en load time:
n_top_observed = int(df["global_id"].max()) + 1 - N_END_SIPMS * 2
if n_top_observed != EXPECTED_N_TOP:
    print(f"WARNING: {n_top_observed} top SiPMs observados, esperados {EXPECTED_N_TOP}. "
          f"Revisa la geometría de la simulación.")
TOP_IDS = tuple(range(N_END_SIPMS * 2, N_END_SIPMS * 2 + n_top_observed))

# Opción B — módulo compartido analysis/geometry.py:
# Único punto de verdad; todos los scripts importan de ahí.
# Cuando cambie la geometría, un solo commit actualiza geometry.py y todos los scripts quedan alineados.
```

Para ROOT (C++): leer `tree->GetMaximum("global_id")` al inicio y construir las agrupaciones dinámicamente, como hacía el código original antes del cambio de la rama.

La Opción A (derivación de datos) es más robusta: detecta y avisa si los datos no coinciden con el esperado. La Opción B (módulo compartido) es más simple pero sigue requiriendo un cambio manual cuando cambia la geometría.

**Estado:** propuesta pendiente de decisión del usuario antes de §3.

---

### `analysis/merge_runs.py` (M) — 2 líneas

Actualización de nombres de scripts en mensaje de ayuda: `analyze.py` → `analyze_dCFD.py`, `resolution_vs_x.py` → `resolution_vs_x_dCFD.py`. **APLICAR.**

### `analysis/ResolutionScan_v2.C` (M) — ~90 líneas

Añade `drawFits` parámetro, diagnóstico de entradas, `vector<TH1D*> fitHists`, fallback RMS cuando el ajuste gaussiano falla, `TLatex`. Mejora sustancial de robustez. **APLICAR.**

---

## §2.5c — Decisión pendiente: modelos de pulso τ_r / τ_f

`engine.py` y `sipm_waveform_dcfd.py` coexisten en el mismo árbol pero usan parámetros de pulso incompatibles:

| Parámetro | `engine.py` (EXEC_14) | `sipm_waveform_dcfd.py` |
|-----------|----------------------|------------------------|
| τ_rise | 0.5 ns (SPR ej200_endonly `common.py`) | 2.0 ns ("Broadcom") |
| τ_fall | 5.0 ns (SPR ej200_endonly `common.py`) | 55.0 ns ("Broadcom"; DS105 τ_fall=55ns) |
| Dispositivo | SiPM de ej200_endonly (no identificado) | AFBR-S4N66P024M (Broadcom NUV-MT) |

Que no se solapen funcionalmente no significa que ambos sean correctos para el detector de este proyecto. El SiPM físico es el Broadcom AFBR-S4N66P024M. Si el parámetro relevante es τ_fall del AFBR-S4N66P024M, DS105 confirma 55 ns → `sipm_waveform_dcfd.py` tiene la procedencia correcta y `engine.py` no. Si el parámetro de `engine.py` proviene de un SPR diferente y ese SPR es el que corresponde al detector, es al revés.

**Estado:** decisión pendiente del usuario. No se resuelve aquí; se registra para que quede visible en el árbol antes de §3.

---

## §2.6 Hallazgo adicional — divergencia de simulación

La rama modifica **30 archivos** de `src/`, `include/`, `main.cc`, `CMakeLists.txt` y macros de Geant4 respecto a main. Los scripts de análisis fueron desarrollados contra esa versión de simulación, que puede tener un esquema de TTree diferente al de Fase 7.

Riesgo concreto: si la versión de simulación de la rama emite ramas adicionales (p.ej. `sptr_jitter_ns`, `sensor_model_id`) que la Fase 7 no emite, los scripts que esperan esas ramas producirán columnas de cero en lugar de datos físicos, sin error. Se documenta como riesgo; no bloquea la migración de scripts TTree-compatibles confirmados.

---

## Resumen de veredictos para §3

| Categoría | Acción propuesta |
|-----------|-----------------|
| Scripts `analysis/` únicos con TTree COMPLETA | **MIGRAR** con `git mv` desde la rama |
| `analyze_dCFD.py`, `sipm_waveform_dcfd.py/.cpp` (PARCIAL) | **MIGRAR** con nota de compatibilidad parcial |
| `analysis/sipm_waveform_dcfd` (binario 83 KB) | **NO MIGRAR** |
| Duplicados raíz (`fpt_vs_n_profile.C`, `fpt_vs_n_profile_batch.C`) | **NO MIGRAR** (canónico en `analysis/`) |
| `fpt_vs_n_profile_batch_slides.C`, `png_to_beamer.py`, `fpt_manifest_to_beamer.py` | **MIGRAR** a `analysis/` |
| `presentacion_timing_detector_scan_ej-230.tex` | **MIGRAR** a `presentations/` |
| `macros/vis_scan_png*.mac` | **MIGRAR** a `macros/` |
| M neutros — `.to_numpy()`, `expected_x`, `--x-min/max/step`, tolerancia scan | **APLICAR** |
| M de geometría — `N_TOP: 70→20`, `TOP_IDS: 16-85→16-35`, `hOccupancy: 86→60` | **NO APLICAR** — ajustan geometría de la rama, romperían datos de 86 canales en main (ver §2.5b) |
| `engine.py` vs `sipm_waveform_dcfd.py` | **COEXISTEN** — modelos de pulso distintos (ver §2.4); cuál es correcto para este detector es decisión pendiente (ver §2.5c) |
| Modificados `src/`, `include/`, `main.cc`, `CMakeLists.txt` | **NO MIGRAR** (out of scope per §5 hard rules) |
| `.claude/settings.json`, `.codex` | **NO MIGRAR** |
| `audit/runs_file_inventory.txt` | **NO MIGRAR** (snapshot externo de W2) |

---

## Correcciones aplicadas al informe

| Corrección | Estado |
|-----------|--------|
| C1: §2.4 rehecho — `engine.py` existe en main; comparación completa realizada | ✓ |
| C2: `transit-sigma-ps=200ps` reclasificado como PROCEDENCIA IDENTIFICADA (Lv et al., dispositivo distinto) | ✓ |
| `tau-f-ns=55ns`: PROCEDENCIA IDENTIFICADA — DS105 tabla Optical and Electrical Features | ✓ |
| `tau-r-ns=2.0ns`: PROCEDENCIA PARCIAL — atribuido a Broadcom en código; no aparece en DS105 | ✓ |
| 7 diffs M: leídos y auditados; 4 neutros aprobados; 3 de geometría bloqueados | ✓ |
| Verificación de geometría: `kNTopSiPMs=70` en main confirmado (DetectorConstruction.hh:38, exec13 CSV) | ✓ |
| Propuesta de derivación en tiempo de ejecución documentada en §2.5b | ✓ |
| Decisión pendiente τ_r/τ_f documentada en §2.5c | ✓ |

---

*Fin del checkpoint §2 (completo). Esperando aprobación para proceder con §3.*
