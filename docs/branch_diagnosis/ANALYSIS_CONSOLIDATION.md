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

## 6. SPTR — parámetros de sipm_waveform_dcfd

| Parámetro | Valor default | Procedencia | Estado |
|-----------|--------------|-------------|--------|
| `--transit-sigma-ps` | 200 ps | Lv et al. (2026) — SPTR de su SiPM, no del AFBR-S4N66P024M | **Dispositivo equivocado**; Lee 2025: σ ≈ 58 ps intrínseco para 6×6 mm² NUV-MT |
| `--tau-f-ns` | 55.0 ns | DS105 tabla Optical and Electrical Features | Correcto para AFBR-S4N66P024M |
| `--tau-r-ns` | 2.0 ns | Atribuido a "Broadcom"; no en DS105 | Sin cita formal |
| `--fraction` (dCFD) | 0.14 | Sin cita | Sin origen documentado |
| `--electronics-sigma-ps` | 0 ps | Sin cita | Contrasta con 30 ps (FastIC+) de `analyze_dCFD.py` |

Ver `docs/branch_diagnosis/SPTR_PROVENANCE.md` para la fuente Lee 2025.

---

*Fin del documento de consolidación.*
