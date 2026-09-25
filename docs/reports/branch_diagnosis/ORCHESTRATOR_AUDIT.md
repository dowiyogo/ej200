# ORCHESTRATOR_AUDIT.md — Auditoría de ej200_orchestrator

Fecha: 2026-09-01  
Repositorio auditado: `https://github.com/dowiyogo/ej200_orchestrator`  
Clone de trabajo: scratchpad de sesión (solo lectura — no modificado)  
Commits auditados: 19 (rama única `main`, sin tags)  
Alcance: solo lectura, sin ejecutar código, sin abrir ROOT files  

---

## Veredicto de función

**El repositorio es un orquestador de análisis y compilador de Beamer**, no un simulador.
No contiene código Geant4, no corre simulaciones, no escribe datos ROOT primarios.

Su ciclo de trabajo es:

```
ROOT files en t0minidaq (externos, solo lectura)
    └─► sidecar scripts (f1–f7_sidecar.py, engine.py)
            └─► .root + .csv + .meta.json de análisis
                    └─► render_figures.py → PDFs de figuras
                            └─► exec14_deck*.tex → lualatex → PDFs del deck
```

Todos los artefactos intermedios y finales (PDFs de figuras, ROOT de análisis, PDFs
del deck, archivos auxiliares LaTeX) están comprometidos en git.

---

## 1. Inventario de contenido

### 1.1 Scripts de análisis (`analysis/`)

| Archivo | Función |
|---------|---------|
| `datasets.py` | Registro de metadatos de datasets (ver §2) |
| `engine.py` | Motor de análisis: lectura PyROOT, estimadores SUM4, fit_core |
| `batch_runner.py` | Lanzador de sidecars en lote |
| `manifest_checker.py` | Verifica integridad del `beamer_manifest.csv` vs outputs |
| `render_figures.py` | Renderiza PDFs de figura desde sidecars (escribe en `outputs/`) |
| `f1_sidecar.py` | F1: distribución de tiempo SUM4 (2 gaussianas, rango 0–2 ns) |
| `f2_new_sidecar.py` / `f2_sidecar.py` | F2: σ_t vs N_pe por posición (31 puntos, fit Poisson) |
| `f3_sidecar.py` | F3: event display desde branches x_mm/y_mm/z_mm |
| `f4_sidecar.py` | F4: perfiles TOP T4 vs T20 (solo EJ-204 EndTop) |
| `f5_sidecar.py` | F5: σ_t SUM4 L/R en 3 posiciones (0, 400, 690 mm) |
| `f6_sidecar.py` | F6: correlación N_pe entre IDs 49/50/51/52 (Top redundancia) |
| `f7_sidecar.py` | F7: estadísticos de orden ⟨t_n⟩ por canal Top |
| `f1_tof_probe.py` | QA-1c: prueba de falsificabilidad del ToF geométrico |

### 1.2 Decks LaTeX

| Archivo | Líneas | PDF (bytes) | Estado |
|---------|--------|------------|--------|
| `exec14_deck.tex` | 494 | 924 372 | Base, QA-3b final |
| `exec14_deck_aclarado.tex` | 615 | 1 060 555 | **Canónico** (ver §3) |
| `exec14_deck_new.tex` | **no existe** | 939 937 | Artefacto de build (ver §3) |

**Artefactos auxiliares LaTeX comprometidos (archivos intermedios de lualatex):**
`exec14_deck.{aux,log,nav,out,snm,toc}` +
`exec14_deck_aclarado.{aux,nav,out,snm,toc}` +
`exec14_deck_new.{aux,log,nav,out,toc}` — 15 archivos (~500 KB en total).
Estos son artefactos de build; no son fuente y no deben migrarse.

### 1.3 Figuras comprometidas (`outputs/`)

5.1 MB totales. Incluyen:
- 17 PDFs de figura (fig01a–fig07_ej230)
- 14 PNGs (previews)
- 17 ROOT files de análisis
- 17 CSVs de análisis
- 17 meta.json
- `contact_sheet.png` (hoja de contacto)
- Macros ROOT auxiliares: `fig06_ej230_macros.tex`, `fig07_ej230_macros.tex`

### 1.4 Documentación QA

| Archivo | Contenido |
|---------|----------|
| `PLAN_QA0.md` | Plan inicial: inventario de datos, decisiones de diseño, topología SUM4 |
| `EXEC_14_guia_slide_por_slide.md` | Guía slide-por-slide: descripciones, caveats, checklist QA |
| `beamer_manifest.csv` | Registro de 22 paneles: figura, material, dataset, sha256 de ROOT input, script sidecar |
| `MAIN_REFLECTOR_FIX_STRATEGY.md` | Estrategia del fix de volúmenes reflector → skin surface en `main` de ej200 |
| `REFLECTOR_SKIN_FIX_RECOVERY_REPORT.md` | Informe de recuperación del fix (sesión de reparación) |
| `REFLECTOR_SKIN_FIX_FINAL_REMOTE_CLOSURE.md` | Cierre remoto: validación en 9 ramas, builds, ctest, push results |
| `REFLECTOR_SKIN_FIX_RECOVERY_STATUS.csv` | Tabla de estado por rama (pre/post SHA, push status) |

### 1.5 Diffs y logs del fix de reflector (`recovery_diffs/`, `final_remote_closure_logs/`)

428 KB y 172 KB respectivamente. Son snapshots de las ramas de ej200 en el momento
del fix de reflector (2026-06-18). Útiles como archivo histórico; no son fuente.

---

## 2. Relación con datasets.py

`datasets.py` es un **registro de metadatos puro — solo lectura**. No contiene lógica
de análisis. Mapea nombres de datasets (strings) a dicts con:

| Clave | Tipo | Ejemplo |
|-------|------|---------|
| `material` | str | `"EJ-204"`, `"EJ-230"` |
| `tau_d_ns` | float | `1.8`, `1.5` |
| `tau_r_ns` | float | `0.5` |
| `readout` | str | `"ENDONLY"`, `"ENDTOP"` |
| `data_dir` | Path | `/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000` |
| `n_channels` | int | `86` |
| `opsc_code` | str | `"opsc-102"`, `"opsc-106"` |
| `provenance` | dict | citas a rutas de ej200/SSLG4 para auditoría |

Registra dos datasets:
- `exec07_endtop_2000`: EJ-204, EndTop, τ_d=1.8 ns (opsc-102)
- `ej230_endtop`: EJ-230, EndTop, τ_d=1.5 ns (opsc-106)

**Nota:** Los datasets EndOnly (`endonly_mylar_20260614`, `endonly_mylar_230`) **no están
en datasets.py** — los sidecars F1/F2/F3/F5 los referencian directamente vía rutas hardcoded
en `beamer_manifest.csv`. Solo los datasets EndTop están registrados formalmente en datasets.py.

La función pública es `get(name: str) -> dict`; los sidecars la importan como
`from analysis.datasets import get`.

`datasets.py` no duplica contenido de ej200: cita rutas de ej200 como provenance
(para auditoría), pero no copia código ni datos.

---

## 3. Deck canónico

**Canónico: `exec14_deck_aclarado.tex`** (615 líneas, PDF 1.04 MB).

Cronología de los decks:

| Commit | Fecha | Acción |
|--------|-------|--------|
| `801243a` QA-3 | 2026-06-16 | Primera compilación del deck (14 slides, lualatex) |
| `db5a4b0` QA-3b | 2026-06-16 | Traducción EN + diagrama ToF + 3 fixes (16 slides); `exec14_deck.tex` estabilizado |
| `7531e1c` QA-3d | 2026-06-16 | Edita `exec14_deck.tex`; compila a `exec14_deck_new.pdf` (PDF original bloqueado por viewer) |
| `8bf0995` Orquestador | 2026-06-16 | Commits archivos auxiliares LaTeX (`exec14_deck*.{aux,nav,toc,...}`) |
| `526eec2` EXEC_13 close | 2026-06-16 | Crea `exec14_deck_aclarado.tex` (537 → 615 líneas); 5 fixes in-situ; cierre QA |
| `5d70202` EXEC_15 | 2026-06-16 | Añade 3 frames EJ-230 EndTop (F4/F6/F7 EJ-230) a `exec14_deck_aclarado.tex` |

**Por qué `exec14_deck_new.tex` no existe:**
El commit QA-3d explica: "exec14_deck.pdf write-locked by viewer; no overflow warnings; 16 pages."
El compilador lualatex escribió el resultado en `exec14_deck_new.pdf` como workaround.
`exec14_deck_new.tex` nunca fue creado — es idéntico en contenido a `exec14_deck.tex`
en ese punto de la historia. Los archivos `exec14_deck_new.{aux,log,nav,...}` son
artefactos auxiliares del build de ese workaround.

---

## 4. Contenido único vs ej200

### Contenido ÚNICO en el orquestador (no existe en ej200)

| Contenido | Descripción |
|-----------|-------------|
| Beamer deck EJ-204 vs EJ-230 | Comparación de dos materiales de centelleador para SHiP; geometría EndOnly + Mylar; no es la misma presentación que talk_v6 |
| Sidecar framework (f1–f7) | Scripts de análisis PyROOT con fit_core, SUM4 timing, event display, order-statistics Top |
| datasets.py | Registro de metadatos de datasets de simulación |
| beamer_manifest.csv | Registro con SHA256 de inputs ROOT, 22 paneles |
| `PLAN_QA0.md` | Diseño de la sesión EXEC_14, inventario de datos, topología SUM4 verificada |
| `EXEC_14_guia_slide_por_slide.md` | QA slide-por-slide |
| Documentación del fix de reflector (3 archivos) | Proceso interno de reparación de ej200 el 2026-06-18; pre-análisis de ramas |
| `recovery_diffs/*.diff` | Snapshots de 9 ramas de ej200 al momento del fix |
| `final_remote_closure_logs/*.log` | Logs de ctest y geometría de 9 ramas de ej200 |

### Contenido que SOLAPA con ej200

| Contenido | Relación |
|-----------|---------|
| Referencias en `provenance` de datasets.py | Cita rutas de `ej200/src/external/SSLG4` para auditoría — no copia código |
| `MAIN_REFLECTOR_FIX_STRATEGY.md` | Describe cambios a `ej200:main`; el commit del fix está en ej200 como `84e902c` |
| `REFLECTOR_SKIN_FIX_FINAL_REMOTE_CLOSURE.md` | Lista SHAs y estado post-push de 9 ramas de ej200 — duplica información de git history |

---

## 5. Fase óptica de EXEC_14

**Dataset EndOnly:** `endonly_mylar_20260614` (EJ-204) y `endonly_mylar_230` (EJ-230).

- Fecha de colección: 2026-06-14
- Rama de simulación: `feat/endonly-mylar` de ej200, SHA pre-fix = `3ae135f`
- El fix global de reflector (skin surface) se aplicó el 2026-06-18 → SHA post-fix = `6942ee6`
- **Los datos EXEC_14 EndOnly fueron recolectados ANTES del fix de skin surface**

Física en `feat/endonly-mylar` pre-fix (SHA `3ae135f`):
- Reflector: volúmenes físicos `ReflectorXMinus/XPlus/ZMinus/ZPlus/YMinus` como placas (border surfaces)
- Material de reflector: G4_MYLAR (brach "endonly-**mylar**")
- Superficie: `G4LogicalBorderSurface` bar→reflector, propiedades de `CreateBarSkinReflector()` → `dielectric_dielectric`, R=0.95
- Geometría: EndOnly (sin SiPMs TOP), sin cámara de aire explícita
- NO aplica GEN-2 ni GEN-3 del árbol `feat/bar-end-vikuiti` (son ramas distintas)

**Dataset EndTop:** `exec07_endtop_2000` (EJ-204) y `ej230_endtop` (EJ-230).
- EJ-204: `t0minidaq/sslg4/exec07_endtop_2000/` — 86 canales, Top IDs 16-85
- EJ-230: `t0minidaq/results_ej230/data/` — 86 canales, τ_d=1.5 ns (opsc-106)
- Geometría: EndTop, con SiPMs TOP en Y=+30.25 mm (una sola hilera de 70 SiPMs)

**EXEC_14 no es comparable directamente con los resultados de talk_v6:**
- talk_v6 usa GEN-3 física (Phase 7: air gap 0.10 mm + Mylar 0.05 mm, border surface `dielectric_dielectric`)
- EXEC_14 EndOnly usa física pre-Phase-4 (no skin surface, border surface sobre volumen Mylar, sin air gap separado)
- Las métricas σ_END de EXEC_14 no son equivalentes a las de talk_v6

---

## 6. Artefactos LaTeX comprometidos — resumen de tamaños

| Tipo | Archivos | Tamaño total aprox. |
|------|----------|-------------------|
| PDFs del deck | 3 (.pdf) | ~2.9 MB |
| Fuentes .tex | 2 (.tex; exec14_deck_new.tex ausente) | ~55 KB |
| Auxiliares LaTeX | 15 (.aux,.nav,.toc,.snm,.out,.log) | ~200 KB |
| Figuras outputs/ | 17 PDF + 14 PNG + 17 ROOT + 17 CSV + 17 meta.json + miscelánea | 5.1 MB |
| Diffs y logs | recovery_diffs/ + final_remote_closure_logs/ | 600 KB |
| **Total repositorio** | | **~9 MB** |

Los PDFs del deck y los artefactos auxiliares LaTeX (.aux, .nav, .toc, .snm) son
artefactos de build. Los ROOT y CSV de `outputs/` son datos de análisis intermedios.
Ninguno de estos debería migrarse directamente a ej200 sin filtrado.

---

## Propuesta de destino (P2)

| Contenido del orquestador | Destino en ej200 | Observación |
|--------------------------|-----------------|-------------|
| `analysis/f1–f7_sidecar.py`, `engine.py`, `batch_runner.py`, `manifest_checker.py`, `render_figures.py` | `analysis/exec14/` | Migrar fuente Python; omitir `__pycache__/` |
| `analysis/datasets.py` | `analysis/exec14/datasets.py` | Solo dos datasets registrados; conservar provenance |
| `beamer_manifest.csv` | `analysis/exec14/beamer_manifest.csv` | Registro con SHA256 de inputs |
| `exec14_deck_aclarado.tex` | `talks/exec14/exec14_deck.tex` | Solo la fuente canónica; NO incluir PDFs ni aux |
| `PLAN_QA0.md`, `EXEC_14_guia_slide_por_slide.md` | `docs/branch_diagnosis/exec14_qa/` | Documentación QA interna |
| `MAIN_REFLECTOR_FIX_STRATEGY.md`, `REFLECTOR_SKIN_FIX_*.md`, `REFLECTOR_SKIN_FIX_RECOVERY_STATUS.csv` | `docs/branch_diagnosis/reflector_fix/` | Documentación del fix de reflector |
| `recovery_diffs/`, `final_remote_closure_logs/` | NO migrar | Snapshots históricos de 2026-06-18; la información está en git history de ej200 |
| `outputs/*.pdf`, `*.png`, `*.root`, `*.csv`, `*.meta.json` | NO migrar | Artefactos de análisis (reproducibles desde datos + scripts); NO datos fuente |
| `exec14_deck*.pdf`, `exec14_deck*.aux`, etc. | NO migrar | Artefactos de build LaTeX |
| `README.md` (1 línea) | Omitir | Sin contenido |

**Condición para `exec14_deck_aclarado.tex`:** El deck referencia figuras vía
`\graphicspath{{outputs/}}`. La migración requiere ajustar esa ruta o asegurar
que `analysis/exec14/outputs/` exista y esté poblado (con `.gitignore` apropiado
para excluir los artefactos binarios de análisis).

---

## CHECKPOINT — Espera de aprobación antes de P2

La migración (P2) y la promoción de rama (P3) están en espera.

Acciones que requieren aprobación del usuario:
1. Crear directorio `talks/exec14/` con `exec14_deck_aclarado.tex` como `exec14_deck.tex`
2. Crear `analysis/exec14/` con sidecar scripts y datasets.py
3. Crear `docs/branch_diagnosis/exec14_qa/` y `reflector_fix/` con los .md de documentación
4. Commit(s) en `feat/bar-end-vikuiti`
5. Push de `feat/bar-end-vikuiti` a origin
6. P3: promoción de `feat/bar-end-vikuiti` a `main` y cambio de default en GitHub

*Ningún archivo fue modificado durante esta auditoría. Clone en scratchpad de sesión.*
