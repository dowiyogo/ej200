# analysis/ — Índice de scripts de análisis

Scripts de análisis de la simulación Geant4 de EJ-200 (y variantes EJ-204, EJ-230).
Todos leen el TTree `sipm_hits` producido por la simulación, salvo que se indique.

Geometría de referencia (main): `kNEndSiPMs = 8` por lado (IDs 0–15),
`kNTopSiPMs = 70` (IDs 16–85), total 86 canales. Ver
`include/DetectorConstruction.hh:37-39`.

---

## Subdirectorios

| Directorio | Contenido |
|-----------|-----------|
| `timing/` | Estimadores de timestamp: waveform+dCFD, FPT, k-ésimo fotón, barrido de modos. Ver `timing/README.md`. |
| `exec07/` | Pipeline EXEC_07 — end-only, análisis del scan de junio 2026 (óptica anterior). |
| `exec13/` | Pipeline EXEC_13 — EndTop, 86 canales, datos de junio 2026. |
| `exec14/` | Pipeline EXEC_14 — motor `engine.py` (SUM4 leading-edge, SPR ej200_endonly). Ver `exec14/README.md`. |

---

## Scripts en este nivel

### Utilidades generales

| Script | Función |
|--------|---------|
| `merge_runs.py` | Combina múltiples `photon_hits_runNNN.root` en un solo archivo. |
| `convert_vis_exports.py` | Convierte PDFs de visualización (`trace_x_*.pdf`) a PNG vía ImageMagick. |

### Análisis de posición / borde

| Script | Función |
|--------|---------|
| `edge_resolution.C` | σ_t vs posición x para región de borde. |
| `edge_resolution.py` | Python equivalente; produce `edge_summary.csv`. |
| `compare_edge_wraps.py` | Overlay de σ_t y N_pe para tres configuraciones de wrap (mylar/air/black). |
| `resolution_vs_x_fixed.py` | σ_t vs x con FPT; `N_TOP_SIPMS = 70` (geometría main). |
| `resolution_vs_x.py` | Versión básica (legado). |
| `resolution_vs_x_root.C` | Versión ROOT. |

### Agrupaciones y ranking

| Script | Función |
|--------|---------|
| `grouped_resolution.C` | σ_t vs x para esquemas de agrupación de SiPMs (main: IDs 16–85). |
| `grouped_resolution.py` | Python equivalente. |

### Óptica

| Script | Función |
|--------|---------|
| `topreadout_crosstalk.py` | Crosstalk fotónico entre top SiPMs. |

### Análisis básico

| Script | Función |
|--------|---------|
| `ResolutionScan.C` | Scan de resolución temporal básico. |
| `ResolutionScan_v2.C` | Con `drawFits`, diagnóstico de entradas, fallback RMS. |
| `analyze_hits.C` | Histogramas de ocupancia, tiempo y espectro (rango: 86 canales = 0–85). |
| `analyze.py` | Análisis Python básico (legado; ver `timing/analyze_basic.py` para versión actual). |
| `analyze_dCFD.C` | dCFD en ROOT (standalone, sin integración). |

### Legado / test bench

| Script | Función |
|--------|---------|
| `congruent_sum4_timing.C` | SUM4 congruente (predecesor de `exec14/engine.py`). |
| `bar_comparison_4configs.cxx` | Comparación de 4 configuraciones de barra. |
| `bar_end_vikuiti_scan.cxx` | Scan de barra end Vikuiti. |
| `high_stats_position_scan.C` | Scan de alta estadística. |
| `preliminary_position_scan.C` | Scan preliminar. |
| `tb_mirror_*.py`, `tb_mirror_*.C` | Test bench con espejo. |
| `exec07_photon_budget.py` | Photon budget EXEC_07. |

---

## Notas de compatibilidad

- Los scripts de este nivel fueron desarrollados con `kNTopSiPMs = 70`. Los scripts de
  `timing/` fueron importados de una rama con `kNTopSiPMs = 20`. Antes de correr los
  scripts de `timing/` sobre datos de 86 canales, verificar el rango de `global_id`
  asumido. Ver `docs/branch_diagnosis/ANALYSIS_CONSOLIDATION_AUDIT.md §2.5b`.

- `edge_resolution.py` importa desde `resolution_vs_x_fixed.py` (mismo nivel):
  `from resolution_vs_x_fixed import fit_fpt_distribution, gauss, load_event_level_data`.

- `grouped_resolution.py` importa desde `resolution_vs_x_fixed.py`:
  `from resolution_vs_x_fixed import fit_fpt_distribution`.

- `exec14/engine.py` requiere `common.py` en ruta hardcodeada:
  `/home/reriosto/SHiP/ej200_endonly/analysis/exec07` — solo funciona en W1.
