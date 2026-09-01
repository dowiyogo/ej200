# analysis/exec14 — Sidecar Analysis Framework for EXEC_14

Framework de análisis para la sesión EXEC_14 (EJ-204 vs EJ-230 SHiP Timing).
Importado desde `https://github.com/dowiyogo/ej200_orchestrator` (19 commits).
Deck asociado: `presentations/exec14/exec14_deck.tex`.

## Arquitectura

```
datasets.py          ← registro de datasets (metadata, rutas t0minidaq, provenance)
engine.py            ← motor de análisis PyROOT: lectura hits, estimadores SUM4, fit_core
batch_runner.py      ← lanzador de sidecars en lote
manifest_checker.py  ← verifica integridad de beamer_manifest.csv vs outputs/
render_figures.py    ← renderiza PDFs de figura desde sidecars (escribe en outputs/)
beamer_manifest.csv  ← 22 paneles: figura, material, dataset, SHA-256 de ROOT input
f1_sidecar.py        ← F1: distribución tiempo SUM4 (2 gaussianas, 0–2 ns)
f1_tof_probe.py      ← QA-1c: prueba de falsificabilidad ToF geométrico
f2_sidecar.py        ← F2 (original): σ_t vs posición por extremo L/R
f2_new_sidecar.py    ← F2 (QA-3c): σ_t vs N_pe por posición (31 pts, fit Poisson)
f3_sidecar.py        ← F3: event display desde branches x_mm/y_mm/z_mm
f4_sidecar.py        ← F4: perfiles TOP T4 vs T20 (N_pe medio por posición)
f5_sidecar.py        ← F5: σ_t SUM4 L/R en 3 posiciones (0, 400, 690 mm)
f6_sidecar.py        ← F6: correlación N_pe entre IDs 49/50/51/52 (redundancia Top)
f7_sidecar.py        ← F7: estadísticos de orden ⟨t_n⟩ por canal Top
outputs/             ← CSV y meta.json de análisis (ver abajo)
```

## datasets.py — datasets registrados

| Nombre | Material | τ_d | Readout | Ruta en t0minidaq |
|--------|----------|-----|---------|------------------|
| `exec07_endtop_2000` | EJ-204 | 1.8 ns | ENDTOP | `sslg4/exec07_endtop_2000/` |
| `ej230_endtop` | EJ-230 | 1.5 ns | ENDTOP | `results_ej230/data/` |

Los datasets EndOnly (`endonly_mylar_20260614`, `endonly_mylar_230`) están
referenciados directamente en `beamer_manifest.csv` (no en `datasets.py`).

## outputs/

Contiene los CSV y meta.json de los 22 paneles del deck.
Los archivos ROOT, PNG y PDF están excluidos por `.gitignore`.
Los ROOT files de entrada están en t0minidaq (externo a este repo) —
los CSV preservan los números sin necesidad de reejecutar el análisis.

## Para correr sobre datos de Fase 7

Los sidecars acceden a ROOT files en rutas absolutas bajo
`/home/reriosto/SHiP/t0minidaq/` (ver `beamer_manifest.csv`, columna `root_input`).
Para adaptarlos a datos de Fase 7 (scan_end_vikuiti, scan_end_vik_sparse_top_v2):

1. **datasets.py**: añadir nuevas entradas apuntando a los directorios de Fase 7.
2. **beamer_manifest.csv**: actualizar las rutas `root_input` y los `sha256_input`
   para cada panel.
3. **engine.py**: verificar compatibilidad de branches en los nuevos ROOT files
   (Fase 7 usa la misma estructura de árbol que feat/endonly-mylar si los branches
   x_mm/y_mm/z_mm/face_type/global_id están presentes).
4. **f*_sidecar.py**: ajustar `tau_d_ns` según el material (EJ-230: 1.5 ns vs 1.8 ns
   para EJ-204) vía `datasets.get()`.
5. Los sidecars F4/F6/F7 requieren datos EndTop; verificar que los ROOT files de
   Fase 7 incluyan IDs Top (global_id 16–85) si se quieren comparar directamente.

**No se adaptan ahora**: el framework se conserva como código de referencia
para la comparación EJ-204 vs EJ-230 con los datos originales de EXEC_14.
