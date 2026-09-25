# REPO_TRIAGE — Clones ej200 bajo /home/rrios

**Generado:** 2026-08-11 21:38 (t0minidaq)
**Auditoría:** `/home/rrios/repo_triage_20260811_213804/`
**Modo:** solo lectura — no se modificó nada en los repos

---

## Veredicto breve

No se actualizó ningún repo porque todos tienen cambios locales. La acción correcta es guardar y clasificar los cambios antes de hacer pull. Todos los repos están en `ahead=0 / behind=0` respecto a su upstream, por lo que `pull --ff-only` será safe una vez que los working trees estén limpios.

---

## Resumen ejecutivo de riesgos

| repo | branch | archivos sucios | riesgo pull | prioridad acción |
|---|---|---|---|---|
| `ej204` | `feat/endtop-sslg4` | 3 mod + 10 untracked | BAJO | **PRIORITARIO** (código valioso sin commit) |
| `ej230_end` | `feat/ej230-endonly-mylar` | 27 mod + 2 untracked | ALTO | **URGENTE** (mayor volumen, output 2.4 GB) |
| `ej230` | `feat/ej230-sslg4` | 10 mod + 1 untracked | BAJO-MEDIO | NORMAL |
| `ej200` | `exp/pair-scan-2026-06-11` | 43 mod + 1 untracked | MEDIO | NORMAL |
| `ej200_end` | `feat/endonly-mylar` | 4 mod + 1 untracked | BAJO | BAJO |

---

## /home/rrios/ej200

### Estado git

- **branch:** `exp/pair-scan-2026-06-11`
- **upstream:** `origin/exp/pair-scan-2026-06-11`
- **ahead/behind:** `0 / 0`
- **remote:** `https://github.com/dowiyogo/ej200`
- **dirty:** SÍ — 43 archivos modificados + 1 untracked

### Nota: estado de staging

Todos los archivos modificados están **sin staged** (columna 1 vacía, columna 2 = M). No hay nada en el index.

### Cambios locales clasificados

| path | status | categoría | recomendación |
|---|---|---|---|
| `CMakeLists.txt` | M (unstaged) | CONFIG_BUILD | Commit — refactor + GDML opcional para t0minidaq |
| `scripts/run_pair_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit — script del scan pair activo |
| `macros/pairscan/pairscan_x-422.0mm.mac` … `x-462.0mm.mac` (41 archivos) | M (unstaged) | MACRO | Commit — macros del scan generadas programáticamente; son el registro definitivo del scan |
| `build_t0minidaq/` | untracked | GENERATED_BUILD | Agregar a `.gitignore` — 5.7 MB de binarios compilados, regenerables |

### Análisis del CMakeLists.txt

El diff muestra una reorganización estructural (función `configure_ej200_target`, comentarios de sección). El cambio principal es añadir `GDML opcional` en el `find_package`. Cambio **valioso y funcional**.

### Análisis de las macros pairscan

Cada macro cambia solo 2 líneas (semillas random). Son el parámetro del scan de la campaña actual. Conviene commitearlas como bloque.

### Riesgo si se hace pull ahora

**MEDIO** — si el upstream modificó CMakeLists.txt o las macros, habría conflicto. Con 43 archivos modificados y `ahead=0`, el merge sería complejo si hay divergencia upstream.

### Acción recomendada

```
1. Commit CMakeLists.txt + run_pair_scan.sh + macros/pairscan/
2. Agregar build_t0minidaq/ a .gitignore
3. pull --ff-only
```

---

## /home/rrios/ej200_end

### Estado git

- **branch:** `feat/endonly-mylar`
- **upstream:** `origin/feat/endonly-mylar`
- **ahead/behind:** `0 / 0`
- **remote:** `https://github.com/dowiyogo/ej200`
- **dirty:** SÍ — 4 archivos modificados + 1 untracked

### Nota: estado de staging

Todos sin staged.

### Cambios locales clasificados

| path | status | categoría | recomendación |
|---|---|---|---|
| `CMakeLists.txt` | M (unstaged) | CONFIG_BUILD | Commit — probablemente mismo patch GDML/t0minidaq que en otros repos |
| `resume_scan_2.sh` | M (unstaged) | RUN_SCRIPT | Commit — script de reanudación de scan |
| `scripts/run_exec07_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit — adaptación paths t0minidaq |
| `scripts/run_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit — adaptación paths t0minidaq |
| `build_t0minidaq/` | untracked | GENERATED_BUILD | Agregar a `.gitignore` — 6.6 MB, regenerable |

### Riesgo si se hace pull ahora

**BAJO** — solo 4 archivos, todos scripts/config. Si el upstream no tocó exactamente estas rutas, el pull será limpio.

### Acción recomendada

```
1. Commit CMakeLists.txt + scripts/*.sh + resume_scan_2.sh
2. Agregar build_t0minidaq/ a .gitignore
3. pull --ff-only
```

---

## /home/rrios/ej204

### Estado git

- **branch:** `feat/endtop-sslg4`
- **upstream:** `origin/feat/endtop-sslg4`
- **ahead/behind:** `0 / 0`
- **remote:** `https://github.com/dowiyogo/ej200.git`
- **dirty:** SÍ — 3 modificados (2 staged + 1 unstaged) + 10 directorios/archivos untracked

### IMPORTANTE: estado de staging mixto

```
M  .gitignore       ← STAGED (en el index, listo para commit)
M  CMakeLists.txt   ← STAGED (en el index, listo para commit)
 M src/RunAction.cc ← unstaged (solo en working tree, no staged)
```

Hay un commit parcialmente preparado. Antes de cualquier acción conviene verificar si este stage es intencional.

### Cambios locales clasificados

| path | status | categoría | recomendación |
|---|---|---|---|
| `.gitignore` | M staged | CONFIG_BUILD | Commit — añade exclusiones build/, runs/, output/, *.root |
| `CMakeLists.txt` | M staged | CONFIG_BUILD | Commit — hace GDML opcional (OPTIONAL_COMPONENTS gdml) |
| `src/RunAction.cc` | M unstaged | SOURCE_CODE | **Commit URGENTE** — fix funcional: reporta BarSkin correctamente en lugar de OPEN/UNDEFINED; sin este fix la readout config se imprime incorrecta |
| `macros/t0minidaq_scan_5000/` (31 .mac) | untracked | MACRO | Commit — macros definitivas del scan EndTop, 5000 eventos/pos, 24 threads |
| `scripts/run_t0minidaq_endtop_scan_5000.sh` | untracked | RUN_SCRIPT | Commit — script principal del scan con --dry-run/--smoke/--first-only |
| `scripts/analyze_t0minidaq_endtop_corefit.py` | untracked | ANALYSIS_SCRIPT | Commit — análisis Gaussian core fit + bootstrap 200 iter |
| `scripts/run_analysis_t0minidaq_endtop_corefit.sh` | untracked | RUN_SCRIPT | Commit — wrapper bash del análisis |
| `scripts/analyze_propagation_velocity.py` | untracked | ANALYSIS_SCRIPT | Commit — curva de velocidad efectiva de propagación (~15.6 cm/ns) |
| `scripts/analyze_t0minidaq_endtop_simple_std.py` | untracked | ANALYSIS_SCRIPT | Commit — análisis STD exploratorio END_AVG/END_DIFF/END_SUM |
| `scripts/run_analysis_t0minidaq_endtop_simple_std.sh` | untracked | RUN_SCRIPT | Commit — wrapper bash |
| `runs/` | untracked | RAW_DATA + RESULT_OUTPUT | **NO commitear** — 2.6 GB de datos de simulación ROOT + análisis. Agregar a `.gitignore` |
| `build_t0minidaq/` | untracked | GENERATED_BUILD | Agregar a `.gitignore` — 8.1 MB, regenerable con cmake |

### Estructura de runs/

```
runs/
  t0minidaq_endtop_scan_5000_20260618_203124/   ← primer intento
  t0minidaq_endtop_smoke_20260618_203606/        ← smoke test
  t0minidaq_endtop_smoke_20260618_203829/        ← smoke test
  t0minidaq_endtop_scan_20260618_203915/         ← scan intermedio
  t0minidaq_endtop_scan_20260618_204959/         ← scan final (5000 ev/pos, 31 pos)
```

Total: 2.6 GB. Contiene outputs de simulación (ROOT files de hits) y análisis derivados (PNGs, CSVs).

### Riesgo si se hace pull ahora

**BAJO** — los archivos modificados son correcciones locales. `src/RunAction.cc` es el más crítico. Si el upstream ya tiene una versión actualizada de RunAction.cc, podría haber conflicto, pero es un archivo único y el diff es claro.

### Acción recomendada — PRIORITARIA

Este es el repo más activo de la sesión. Tiene código valioso (fix RunAction.cc, scripts de análisis completos, macros de scan) que **no está en el remoto**. Commitear pronto.

```bash
# Sugerido (no ejecutar ahora):
git -C /home/rrios/ej204 add src/RunAction.cc   # añadir el unstaged
git -C /home/rrios/ej204 commit -m "feat(t0minidaq): add EndTop scan macros, analysis scripts, BarSkin fix"
# luego:
git -C /home/rrios/ej204 add macros/t0minidaq_scan_5000/
git -C /home/rrios/ej204 add scripts/run_t0minidaq_endtop_scan_5000.sh
git -C /home/rrios/ej204 add scripts/analyze_t0minidaq_endtop_corefit.py
git -C /home/rrios/ej204 add scripts/run_analysis_t0minidaq_endtop_corefit.sh
git -C /home/rrios/ej204 add scripts/analyze_propagation_velocity.py
git -C /home/rrios/ej204 add scripts/analyze_t0minidaq_endtop_simple_std.py
git -C /home/rrios/ej204 add scripts/run_analysis_t0minidaq_endtop_simple_std.sh
git -C /home/rrios/ej204 commit -m "feat(t0minidaq): add EndTop scan macros and analysis scripts"
```

---

## /home/rrios/ej230

### Estado git

- **branch:** `feat/ej230-sslg4`
- **upstream:** `origin/feat/ej230-sslg4`
- **ahead/behind:** `0 / 0`
- **remote:** `https://github.com/dowiyogo/ej200`
- **dirty:** SÍ — 10 archivos modificados + 1 untracked

### Nota: estado de staging

Todos sin staged.

### Cambios locales clasificados

| path | status | categoría | recomendación |
|---|---|---|---|
| `scripts/build_report_msi.sh` | M (unstaged) | RUN_SCRIPT | Commit — adaptación t0minidaq/MSI |
| `scripts/make_beamer_ej230.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/run_analysis_t0minidaq.sh` | M (unstaged) | RUN_SCRIPT | Commit — específico t0minidaq |
| `scripts/run_center.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_center_msi_16t.sh` | M (unstaged) | RUN_SCRIPT | Commit — 16 threads para MSI |
| `scripts/run_center_t0minidaq_24t.sh` | M (unstaged) | RUN_SCRIPT | Commit — 24 threads para t0minidaq |
| `scripts/run_exec07_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_scan_msi_16t.sh` | M (unstaged) | RUN_SCRIPT | Commit — 16 threads MSI |
| `scripts/run_scan_t0minidaq_24t.sh` | M (unstaged) | RUN_SCRIPT | Commit — 24 threads t0minidaq |
| `build_t0minidaq/` | untracked | GENERATED_BUILD | Agregar a `.gitignore` — 5.6 MB, regenerable |

### Observación

Los scripts `*_t0minidaq_24t.sh` y `*_msi_16t.sh` son adaptaciones del mismo script base para diferentes máquinas/threads. Conviene commitearlos juntos como un bloque de "soporte multi-servidor".

### Riesgo si se hace pull ahora

**BAJO-MEDIO** — 10 scripts sin staged. Si upstream no tocó estos mismos scripts, el pull será limpio. El riesgo sube si el upstream actualizó `run_scan.sh` o similares.

### Acción recomendada

```
1. Commit de todos los scripts modificados en rama feat/ej230-t0minidaq-scripts
2. Agregar build_t0minidaq/ a .gitignore
3. pull --ff-only
```

---

## /home/rrios/ej230_end

### Estado git

- **branch:** `feat/ej230-endonly-mylar`
- **upstream:** `origin/feat/ej230-endonly-mylar`
- **ahead/behind:** `0 / 0`
- **remote:** `https://github.com/dowiyogo/ej200`
- **dirty:** SÍ — 27 archivos modificados + 2 untracked

### Nota: estado de staging

Todos sin staged.

### Cambios locales clasificados

| path | status | categoría | recomendación |
|---|---|---|---|
| `resume_scan_2.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/audit_beamer_assets.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/audit_beamer_numbers.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/audit_exec14b_raster_text.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/build_report_msi.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/check_exec14b_asset_parity.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/check_exec14b_figure_frame_consistency.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/check_exec14b_frame_parity.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/diag_exec14d_endtop_ratio.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/generate_exec14b_tables.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/generate_exec14d_nominal_parameters.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/make_beamer_ej230.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/preflight_exec14b_pdf.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/rebuild_exec14b_report.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `scripts/run_analysis_t0minidaq.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_center.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_center_msi_16t.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_center_t0minidaq_24t.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_exec07_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_exec14b_main_analysis_repair.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_exec14b_special_analysis.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_exec14b_special_sims.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_exec14b_strict_preflight.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_scan.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_scan_msi_16t.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/run_scan_t0minidaq_24t.sh` | M (unstaged) | RUN_SCRIPT | Commit |
| `scripts/validate_exec14b_roots.py` | M (unstaged) | ANALYSIS_SCRIPT | Commit |
| `build_t0minidaq/` | untracked | GENERATED_BUILD | Agregar a `.gitignore` — 5.5 MB, regenerable |
| `output/` | untracked | RESULT_OUTPUT | **NO commitear** — 2.4 GB, 31 archivos ROOT + 101 archivos totales. Mover a `/home/rrios/results_ej230_end/` o agregar a `.gitignore` |

### Estructura de output/

```
output/
  endonly_mylar_t0minidaq_20260614_173944/  ← run de producción
  ... (101 archivos totales, 31 ROOT files)
```

### Riesgo si se hace pull ahora

**ALTO** — 27 archivos modificados, muchos scripts de análisis EXEC14b/EXEC14d que parecen trabajo activo. Si el upstream también modificó alguno de estos scripts (especialmente `run_scan.sh`, `make_beamer_ej230.py`) habrá conflictos. El volumen de `output/` (2.4 GB) también hace que cualquier operación con `--include-untracked` sea lenta.

### Acción recomendada — REQUIERE ATENCIÓN

```
1. ANTES: mover output/ fuera del repo para evitar problemas:
   mv /home/rrios/ej230_end/output/ /home/rrios/results_ej230_end_backup/
2. Commit de todos los scripts en rama feat/exec14b-analysis-scripts
3. Agregar build_t0minidaq/ a .gitignore
4. pull --ff-only
```

---

## Tabla global

| repo | branch | ahead/behind | archivos sucios | motivo principal | acción siguiente | pull seguro después? |
|---|---|---|---|---|---|---|
| `ej204` | `feat/endtop-sslg4` | 0/0 | 3M+10U | RunAction.cc fix + análisis EndTop | Commit staged+unstaged + add scripts | SÍ |
| `ej230_end` | `feat/ej230-endonly-mylar` | 0/0 | 27M+2U | 27 scripts EXEC14b/14d + output 2.4GB | Mover output/, commit scripts | SÍ — tras limpiar |
| `ej230` | `feat/ej230-sslg4` | 0/0 | 10M+1U | Scripts t0minidaq 24t + MSI 16t | Commit scripts | SÍ |
| `ej200` | `exp/pair-scan-2026-06-11` | 0/0 | 43M+1U | CMakeLists + 41 macros pairscan + run script | Commit CMake+macros+script | SÍ |
| `ej200_end` | `feat/endonly-mylar` | 0/0 | 4M+1U | CMakeLists + 3 scripts | Commit scripts | SÍ |

---

## Orden de trabajo recomendado

### 1. ej204 — PRIORITARIO

Tiene código funcional corregido (`src/RunAction.cc` con fix BarSkin) y 6 scripts/herramientas de análisis completamente nuevos que no están en el remoto. Riesgo de pérdida si el disco se llena o hay un problema.

### 2. ej230_end — URGENTE por volumen

Tiene el mayor volumen de cambios locales (27 scripts) más 2.4 GB de output. La carpeta `output/` debe moverse fuera del repo antes de cualquier operación de git que incluya untracked files.

### 3. ej230 — NORMAL

Scripts de adaptación multi-servidor (t0minidaq 24t, MSI 16t). Trabajo limpio, fácil de commitear en bloque.

### 4. ej200 — NORMAL

Las 41 macros pairscan son el producto de la campaña activa. Junto con CMakeLists.txt y run_pair_scan.sh forman un bloque coherente.

### 5. ej200_end — BAJO

Pocos cambios, bajo riesgo. Se puede dejar para el final.

---

## Archivos que deben quedar fuera de git (`.gitignore`)

Agregar a `.gitignore` en **todos los repos** que no lo tengan:

```gitignore
# Build directories
build_t0minidaq/
build/
build-*/
cmake-build-*/

# Simulation outputs (large binary data)
runs/
output/
out/

# ROOT data files
*.root

# Python cache
__pycache__/
*.pyc
*.pyo
```

---

## Archivos valiosos que conviene commitear

### ej204 (PRIORITARIO)

- `src/RunAction.cc` — fix funcional crítico (BarSkin reporting)
- `.gitignore` — ya staged
- `CMakeLists.txt` — ya staged, GDML opcional
- `macros/t0minidaq_scan_5000/*.mac` — 31 macros del scan definitivo EndTop
- `scripts/run_t0minidaq_endtop_scan_5000.sh`
- `scripts/analyze_t0minidaq_endtop_corefit.py`
- `scripts/run_analysis_t0minidaq_endtop_corefit.sh`
- `scripts/analyze_propagation_velocity.py`
- `scripts/analyze_t0minidaq_endtop_simple_std.py`
- `scripts/run_analysis_t0minidaq_endtop_simple_std.sh`

### ej230_end

- Todos los `scripts/audit_*.py`, `scripts/check_*.py`, `scripts/generate_*.py`
- Todos los `scripts/run_exec14b_*.sh`
- `scripts/make_beamer_ej230.py`, `scripts/validate_exec14b_roots.py`
- `resume_scan_2.sh`

### ej230

- Todos los 10 scripts modificados (`run_scan*.sh`, `run_center*.sh`, etc.)

### ej200

- `CMakeLists.txt` — refactor + GDML opcional
- `scripts/run_pair_scan.sh`
- `macros/pairscan/*.mac` — 41 macros del scan pair

### ej200_end

- `CMakeLists.txt`
- `scripts/run_exec07_scan.sh`, `scripts/run_scan.sh`
- `resume_scan_2.sh`

---

## Archivos que deben quedar fuera de git

| repo | path | tamaño | razón |
|---|---|---|---|
| `ej204` | `runs/` | 2.6 GB | ROOT files de simulación, datos crudos no versionables |
| `ej230_end` | `output/` | 2.4 GB | Resultados de análisis (ROOT + derivados), 101 archivos |
| todos | `build_t0minidaq/` | ~6 MB c/u | Binarios compilados, regenerables con cmake |
| todos | `*.root` | variable | Datos crudos binarios |

---

## Plan de acción seguro (T7) — NO ejecutar, solo referencia

### ej204 — commit completo en dos pasos

```bash
# Paso 1: añadir el unstaged y hacer commit del bloque de fixes
git -C /home/rrios/ej204 add src/RunAction.cc
git -C /home/rrios/ej204 commit -m "fix(t0minidaq): BarSkin reflectivity fix + GDML optional + gitignore"

# Paso 2: añadir y commitear los nuevos scripts/macros
git -C /home/rrios/ej204 add macros/t0minidaq_scan_5000/
git -C /home/rrios/ej204 add scripts/run_t0minidaq_endtop_scan_5000.sh
git -C /home/rrios/ej204 add scripts/analyze_t0minidaq_endtop_corefit.py
git -C /home/rrios/ej204 add scripts/run_analysis_t0minidaq_endtop_corefit.sh
git -C /home/rrios/ej204 add scripts/analyze_propagation_velocity.py
git -C /home/rrios/ej204 add scripts/analyze_t0minidaq_endtop_simple_std.py
git -C /home/rrios/ej204 add scripts/run_analysis_t0minidaq_endtop_simple_std.sh
git -C /home/rrios/ej204 commit -m "feat(t0minidaq): EndTop scan macros + analysis scripts (corefit, velocity, std)"

# Después del commit, actualizar
git -C /home/rrios/ej204 pull --ff-only
```

### ej230_end — mover output primero

```bash
# Paso 0: mover output fuera del repo (HACERLO PRIMERO)
mv /home/rrios/ej230_end/output/ /home/rrios/results_ej230_end_backup/

# Paso 1: commit de todos los scripts
git -C /home/rrios/ej230_end add scripts/
git -C /home/rrios/ej230_end add resume_scan_2.sh
git -C /home/rrios/ej230_end commit -m "feat(exec14b/14d): add EXEC14b/14d analysis and audit scripts"

# Paso 2: actualizar
git -C /home/rrios/ej230_end pull --ff-only
```

### ej230 — scripts de adaptación multi-servidor

```bash
git -C /home/rrios/ej230 add scripts/
git -C /home/rrios/ej230 commit -m "feat(t0minidaq): add t0minidaq-24t and MSI-16t script variants"
git -C /home/rrios/ej230 pull --ff-only
```

### ej200 — scan pair

```bash
git -C /home/rrios/ej200 add CMakeLists.txt scripts/run_pair_scan.sh macros/pairscan/
git -C /home/rrios/ej200 commit -m "feat(pairscan): CMakeLists refactor + pair scan macros and run script"
git -C /home/rrios/ej200 pull --ff-only
```

### ej200_end — cambios mínimos

```bash
# Opción A: commit directo
git -C /home/rrios/ej200_end add CMakeLists.txt resume_scan_2.sh scripts/run_exec07_scan.sh scripts/run_scan.sh
git -C /home/rrios/ej200_end commit -m "fix(t0minidaq): update scan scripts and CMakeLists for t0minidaq"
git -C /home/rrios/ej200_end pull --ff-only

# Opción B: stash si no están listos
git -C /home/rrios/ej200_end stash push -m "WIP: t0minidaq script updates"
git -C /home/rrios/ej200_end pull --ff-only
git -C /home/rrios/ej200_end stash pop
```

### Para agregar .gitignore en repos que no tengan build_t0minidaq/ ignorado

```bash
for repo in ej200 ej200_end ej230 ej230_end; do
  echo "build_t0minidaq/" >> /home/rrios/$repo/.gitignore
  echo "*.root" >> /home/rrios/$repo/.gitignore
done
# Luego revisar y commitear el .gitignore en cada repo
```

---

## Directorio de auditoría

```
/home/rrios/repo_triage_20260811_213804/
  REPO_TRIAGE_REPORT.md     ← este archivo
  triage_summary.csv        ← resumen tabular
  ej200/
    git_status.txt          ← output de git status/remote/log/branch
    changes.txt             ← diff --name-status + untracked
    diffs/
      CMakeLists.txt.diff
      run_pair_scan.sh.diff
      pairscan_sample.diff  ← muestra de un mac modificado
  ej200_end/
    git_status.txt
    changes.txt
    diffs/
      CMakeLists.txt.diff
      run_scan.sh.diff
  ej204/
    git_status.txt
    changes.txt
    diffs/
      src_RunAction.cc.diff
      CMakeLists.txt.diff
      gitignore.diff
  ej230/
    git_status.txt
    changes.txt
    diffs/
      run_scan_t0minidaq_24t.sh.diff
  ej230_end/
    git_status.txt
    changes.txt
```
