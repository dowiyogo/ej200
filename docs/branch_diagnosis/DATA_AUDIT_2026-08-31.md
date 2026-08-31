# DATA_AUDIT_2026-08-31 — Auditoría de datasets de simulación

**Rama de trabajo:** `docs/branch-diagnosis-2026-08-31`  
**Fecha de auditoría:** 2026-08-31  
**Fuente de datos:** `/home/reriosto/SHiP/t0minidaq/` (acceso solo lectura vía SSH)  
**Cronología óptica de referencia:** establecida en `IMPACT_surface_model_2026-08-31.md` — no se recalcula aquí

---

## Propósito

Determinar bajo qué modelo óptico se produjo cada dataset de simulación que existe en `t0minidaq`. La motivación directa es el §4 de `IMPACT_surface_model_2026-08-31.md`: ninguna comparación entre datasets es válida sin saber qué fase óptica usó cada uno. El eje temporal (σ_t) es el observable de tesis; varios datasets fueron producidos bajo modelos que no conservan la reflexión interna total (TIR) o que tienen el bug de group-velocity aliasing.

**Método de asignación.** Se usa la jerarquía de evidencia de PROMPT A:
1. Hash de commit registrado en un log de escaneo estructurado
2. Estructura geométrica visible en el log de Geant4 (volúmenes registrados en el overlap check)
3. Rama inferida del material del escintilador (`opsc-101` = EJ-204, `opsc-106` = EJ-230)
4. Fecha del directorio / nombre del directorio (evidencia más débil)

---

## Cronología de referencia (resumen)

| Fase | Commit | Fecha | Descripción |
|------|--------|-------|-------------|
| 0 | `ed61b9a` | 2026-03-19 | Sin superficie óptica |
| 1–2b | varios | 2026-05-06–2026-06-09 | Near-GEN-1: border/skin, `dielectric_metal` |
| **3** | `1f6aca1` | **2026-06-10** | **Near-GEN-1 SSLG4: border surfaces sobre paneles de 0.5 µm, `dielectric_metal`** |
| **4** | `f39b84c` | **2026-06-18** | **GEN-1: `G4LogicalSkinSurface` sobre `barLV`, `dielectric_metal`, sin TIR** |
| 5 | `fc5d8b0` | 2026-06-21 | SK-DD: skin sobre barLV, `dielectric_dielectric` R=0.95, sin air gap |
| 6 | `4cb2c23` | 2026-08-12 | 2T-MYLAR: air gap explícito + paneles Mylar; **bug group-velocity aliasing** |
| 7 | `5576687` | 2026-08-14 | 2T-MYLAR corregido: `GROUPVEL_air = c/n_bar = 18.97 cm/ns` |

Las Fases 6–7 corresponden a la rama `feat/endtop-sslg4`. La rama `feat/bar-end-vikuiti` (REF_BASE) hereda el modelo 2T-VIK (`dielectric_dielectric` R=0.95) introducido en EXEC_22 via `CreateBarSkinReflector()`, que no aparece en la cronología de fases de `feat/endtop-sslg4` porque es el modelo de referencia de la otra rama.

---

## Inventario de datasets

### Tabla maestra

| ID | Directorio en t0minidaq | Fecha | Rama | Hash commit | Fase / Modelo | Método | Confianza |
|----|------------------------|-------|------|------------|---------------|--------|-----------|
| A | `scan_end_wrapped_2026-06-09/` | 2026-06-09/10 | `feat/endtop-sslg4` (inferida) | **no registrado** | Fase 2b o Fase 3 | fecha + nº EXEC | **AMBIGUO** |
| B | `sslg4/exec07_endtop_2000/` | 2026-06-10 | `feat/endtop-sslg4` | **no registrado** | **Fase 3** | geometría + EXEC | **ALTA** |
| C | `sslg4/exec08b_window_dip/` | 2026-06-10 | `feat/endtop-sslg4` | **no registrado** | **Fase 3** | geometría + fecha | **ALTA** |
| D | `results_ej230/` | 2026-06-12 | `feat/ej230-sslg4` (inferida) | **no registrado** | **Fase 3-equiv. en feat/ej230-sslg4** | geometría + fecha + material | **ALTA** |
| E | `endonly_mylar_20260614/` | 2026-06-14 | `feat/endonly-mylar` | **`3ae135f`** | feat/endonly-mylar @ 3ae135f (border + Mylar R=0.90) | hash en metadata | **ALTA** |
| F | `runs/t0minidaq_endtop_scan_20260618_204959/` | 2026-06-18 | `feat/endtop-sslg4` | **`a0368c4`** | **Fase 4 (GEN-1), sin TIR** | hash en master_scan.log | **ALTA** |

---

### Dataset A — `scan_end_wrapped_2026-06-09/`

**CODEX:** EXEC_03b  
**Archivos:** 31 ROOT (`photon_hits_run000.root` – `run030.root`), 1 macro (`scan_end_wrapped.mac`)  
**Sin log de Geant4.** No hay ficheros `.log` en el directorio.  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | no registrado | — |
| Log Geant4 | no existe | — |
| EXEC número | EXEC_03b (el commit de Phase 3 se llama "feat(exec07)"; EXEC_03 < EXEC_07) | débil |
| Fecha directorio | 2026-06-09 (antes del commit `1f6aca1` de Phase 3 el 2026-06-10) | débil |
| Material | EJ-204 (mac: `/det/scintillatorMaterial EJ204`) | rama inferida |
| Eventos/pos | 10 000 ev/pos × 31 posiciones | — |

**Asignación:** AMBIGUO — Fase 2b o Fase 3.

El directorio lleva fecha 2026-06-09 y los archivos ROOT están fechados 2026-06-10 (la simulación corrió de noche). El commit `1f6aca1` que introduce Phase 3 se realizó el 2026-06-10 con el nombre exacto "feat(exec07)". El número de exec (03b) es anterior al exec07 de Phase 3, lo que sugiere que esta ejecución es pre-Phase 3. Sin hash de commit en ningún log, no se puede determinar con certeza si el binario compilado correspondía a Phase 2b o al commit inicial de Phase 3.

**Presentaciones que usan estos datos:** los commits `22adf23` (2026-06-08, "scan analysis macro + beamer source") y `9546c84` (2026-06-08, "high-stats scan analysis + updated beamer") son anteriores al commit Phase 3. El "high-stats scan" al que se refieren es la versión anterior de este dataset con menor número de posiciones; la versión de 31 posiciones en este directorio puede corresponder a estos commits o a una repetición posterior.

---

### Dataset B — `sslg4/exec07_endtop_2000/`

**CODEX:** EXEC_07  
**Archivos:** ≥17 ROOT en el raíz (`photon_hits_x0mm.root` etc.), ≥11 subdirectories de análisis (`exec13_out_1` – `exec13_out_11`, `exec13_manual_out`)  
**Fecha archivos ROOT:** 2026-06-10  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | no registrado en log | — |
| Log Geant4: geometría | `ReflectorYMinusPV (G4Box)`, `ReflectorZMinusPV (G4Box)`, `ReflectorZPlusPV (G4Box)`, `ReflectorYPlusPV (G4SubtractionSolid)` — volúmenes físicos reflectores, **sin `G4LogicalSkinSurface`** | **fuerte** |
| Material | `opsc-101` = EJ-204 | rama inferida |
| Fecha | 2026-06-10 (= fecha de `1f6aca1`) | corroborador |
| Nombre EXEC | EXEC_07 coincide con "feat(exec07)" del commit `1f6aca1` | **fuerte** |
| Eventos/pos | 2000 ev/pos | — |

**Asignación: Fase 3.** Confianza ALTA.

El nombre EXEC_07 coincide exactamente con el commit `1f6aca1` ("feat(exec07): add SSLG4 EndTop 86-channel geometry") que introduce Phase 3. La presencia de volúmenes `ReflectorYMinusPV` / `ReflectorZMinusPV` / `ReflectorZPlusPV` / `ReflectorYPlusPV` en el log de Geant4 confirma que se usaron paneles físicos con border surfaces — modelo Near-GEN-1 SSLG4, sin `G4LogicalSkinSurface`.

**Uso en análisis:** los subdirectorios `exec13_out_*/` y los macros `exec13_manual_v*.C` dentro de este directorio confirman que los resultados de **EXEC_13 EJ-204** (commits `c030e88`, `a7d22a3`, `b956102`) se derivan de estos ROOT files (Fase 3).

---

### Dataset C — `sslg4/exec08b_window_dip/`

**CODEX:** EXEC_08b  
**Archivos:** 4 ROOT (`photon_hits_run_A_x-652mm.root`, `run_B_x-642mm.root`, `run_C1_x-648mm.root`, `run_C2_x-654mm.root`), 2 logs  
**Fecha archivos:** 2026-06-10  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | no registrado | — |
| Log Geant4: geometría | `ReflectorYMinusPV (G4Box)`, `ReflectorZMinusPV (G4Box)`, `ReflectorZPlusPV (G4Box)`, `ReflectorYPlusPV (G4SubtractionSolid)` — idéntico a Dataset B | **fuerte** |
| Material | `opsc-101` = EJ-204 | — |
| Fecha | 2026-06-10, misma que Dataset B | corroborador |
| Posiciones | 4 posiciones en zona end-side (~x=−642 a −654 mm) | — |

**Asignación: Fase 3.** Confianza ALTA.

Misma geometría de volúmenes que Dataset B, misma fecha. Este dataset estudia el "window dip" en la zona de SiPMs (4 posiciones cercanas al extremo). Commit de análisis: `9613fc9` "analysis(exec08b): diagnose Top-window collection pattern".

---

### Dataset D — `results_ej230/`

**Sin nombre CODEX identificado.**  
**Archivos:** 33 ROOT en `data/`, 31 logs en `logs/`  
**Fecha archivos:** 2026-06-12  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | no registrado en los logs (los logs de Geant4 arrancan directamente sin header de escaneo estructurado) | — |
| Log Geant4: geometría | `ReflectorYMinusPV (G4Box)`, `ReflectorZMinusPV (G4Box)`, `ReflectorZPlusPV (G4Box)`, `ReflectorYPlusPV (G4SubtractionSolid)` — volúmenes físicos, **sin `G4LogicalSkinSurface`** | **fuerte** |
| Material | `opsc-106` = **EJ-230** → rama `feat/ej230-sslg4`, no `feat/endtop-sslg4` | **fuerte (discrimina rama)** |
| Fecha | 2026-06-12 < 2026-06-18 (GEN-1) | corroborador |

**Asignación: Fase 3-equivalente en `feat/ej230-sslg4`.** Confianza ALTA por geometría + fecha; AMBIGUO para commit exacto.

Este dataset es de una **rama diferente** (`feat/ej230-sslg4`, escintilador EJ-230) y no está sujeto directamente a la numeración de las 7 fases de `feat/endtop-sslg4`. La geometría con volúmenes reflectores físicos confirma que **no se usó GEN-1** (`G4LogicalSkinSurface`). La fecha anterior al 2026-06-18 corrobora. Los 33 archivos ROOT (2 más que las 31 posiciones estándar) sugieren posiciones adicionales o archivos de análisis.

**Uso en análisis:** commits `5b93f4c` (EXEC_13-230) y `0b26c5f`, `f22de30` (EXEC_14 EJ-230) usan estos datos.

---

### Dataset E — `endonly_mylar_20260614/`

**Sin nombre CODEX identificado en metadata.**  
**Archivos:** 31 ROOT (posiciones x=-690 a x=+690 mm), logs por posición, `run_metadata.txt`, `analysis/sigma_t_sum4.csv`, `analysis/attenuation_curve.csv`  
**Fecha inicio:** `2026-06-14T00:02:12+02:00` (exacta, del campo `start` en `run_metadata.txt`)  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | **`3ae135f`** (`run_metadata.txt`: `commit=3ae135f...`) | **máximo** |
| git_dirty | 4 archivos modificados sin commit en runtime | riesgo reproducibilidad |
| Rama del hash | `feat/endonly-mylar` (única rama que contiene `3ae135f`) | confirma |
| Modelo óptico en `3ae135f`:DC.cc | border surfaces (`G4LogicalBorderSurface`) en −Y, −Z, +Z; `CreateMylarReflector(R=0.90, lobe=1.0, sigma=0.1)` cuando `top_surface=mylar`; **sin `G4LogicalSkinSurface`** | — |
| Parámetros de superficie | `mylar_reflectivity=0.90`, `mylar_specular_lobe=1.0`, `mylar_sigma_alpha_deg=0.1` | — |
| Eventos/pos | 2000 ev/pos × 31 posiciones, 24 threads | — |

**Asignación: `feat/endonly-mylar` @ `3ae135f`.** Confianza ALTA.

Este dataset es de la **rama `feat/endonly-mylar`**, que es una rama separada de `feat/endtop-sslg4`. El modelo óptico de `3ae135f` usa border surfaces sobre paneles (no skin surface en barLV), con Mylar como reflector del extremo: `CreateMylarReflector` (`dielectric_metal`, R=0.90). Es conceptualmente cercano a Phase 3 (border, dielectric_metal) pero sobre la geometría end-only de esa rama.

**Nota de reproducibilidad:** `git_dirty=4` indica que 4 archivos estaban modificados localmente en el host de simulación cuando corrió. Los resultados pueden no ser exactamente reproducibles desde el commit `3ae135f` limpio.

**Uso en análisis:** commits `489d238` ("feat(presentation): add end-only + Mylar EJ-204 Beamer deck pipeline") y `fb3749d` ("docs(beamer): update endonly-mylar presentation with corrections") usan estos datos. El commit `c120c51` ("fix yield 10000→10400") corrige una cuenta errónea en el análisis de atenuación derivado de este dataset.

---

### Dataset F — `runs/t0minidaq_endtop_scan_20260618_204959/`

**Nombre propio (master_scan.log):** "EndTop Scan 5000"  
**Archivos:** 31 ROOT en `outputs/x*/photon_hits_run000.root`; logs en `logs/` (master_scan.log 735 KB + 31 run logs); análisis en `analysis_corefit/` y `analysis_simple_std/`; `summary.tsv`  
**Fecha inicio:** `2026-06-18 20:49:59` (exacta, de `master_scan.log` y `summary.tsv`)  
**Fecha fin:** `2026-06-18 20:55:30` (todas las 31 posiciones, status=OK)  

**Evidencia de fase:**

| Evidencia | Valor | Peso |
|-----------|-------|------|
| Hash commit | **`a0368c4`** (`master_scan.log`: `Commit : a0368c4`) | **máximo** |
| Rama | `feat/endtop-sslg4` (`master_scan.log`: `Branch : feat/endtop-sslg4`) | confirma |
| Ancestría | `git merge-base --is-ancestor f39b84c a0368c4` → **True** → GEN-1 ya está en el árbol de `a0368c4` | **máximo** |
| Commit msg `a0368c4` | "With G4LogicalSkinSurface, every photon reflection at the bar surface creates a Bar→WorldLV boundary step..." — **menciona explícitamente `G4LogicalSkinSurface`** | confirma |
| Commit msg `f39b84c` (ancestro) | "fix(optics): replace reflector volumes with bar skin surface" — es el commit GEN-1 | — |
| Modelo óptico resultante | `G4LogicalSkinSurface` sobre `barLV`, `dielectric_metal` → **sin TIR** (θ_c = 39.3° para n=1.58 queda eliminado) | — |
| Eventos/pos | 5000 ev/pos × 31 posiciones, 24 threads | — |

**Asignación: Fase 4 (GEN-1).** Confianza ALTA.

**Este es el dataset más crítico de la auditoría.** Fue producido exactamente en la transición Phase 3→4, usando un binario compilado desde `a0368c4`, que es descendiente directo del commit GEN-1 (`f39b84c`). El mensaje del propio commit `a0368c4` confirma que `G4LogicalSkinSurface` estaba activo. El modelo GEN-1 elimina la reflexión interna total (TIR) en la barra, reduciendo drásticamente el Npe en los SiPMs de extremo para posiciones alejadas.

**Análisis producidos (solo en t0minidaq, no comprometidos en git):**

Los subdirectorios `analysis_corefit/` y `analysis_simple_std/` (fechados 2026-06-19) contienen métricas de σ_t para readout END y TOP en las 31 posiciones. El fichero `analysis_corefit/summary/ANALYSIS_CORE_FIT_REPORT.md` (generado 2026-06-19 11:51:51) registra el patrón de resultados: el readout END muestra múltiples posiciones con estatus `LOW_EFF` y `LOW_STATS` — resultado consistente con la eliminación de TIR bajo GEN-1 que reduce la colección de fotones en los SiPMs de extremo para posiciones centrales. El readout TOP no se ve afectado por el mismo defecto.

**Comprometido en git:** NO. Ninguno de los análisis de este dataset aparece en el historial de git. Los commits entre `a0368c4` y `fc5d8b0` (Phase 5) son solo cambios de infraestructura: `c27f5d4` (macros de escaneo), `ad65154` (fix BarSkin). No se añadió ninguna presentación ni resultado de análisis derivado de este dataset al repositorio.

---

## Clasificación por fase

### Datasets bajo Fase 3 (Near-GEN-1 SSLG4) o equivalente

Los siguientes datasets fueron producidos bajo el modelo de Phase 3 (border surfaces sobre paneles físicos, `dielectric_metal`, sin air gap explícito, sin `G4LogicalSkinSurface`):

- **Dataset B** (exec07_endtop_2000, EJ-204) — Confianza ALTA
- **Dataset C** (exec08b_window_dip, EJ-204) — Confianza ALTA
- **Dataset D** (results_ej230, EJ-230, rama `feat/ej230-sslg4`) — Confianza ALTA por geometría + fecha

Adicionalmente, **Dataset E** (endonly_mylar_20260614, rama `feat/endonly-mylar`) usa border surfaces sin skin — modelo geométricamente próximo a Phase 3 pero en una rama diferente con Mylar explícito como parámetro de superficie.

### Datasets bajo Fase 4 (GEN-1, sin TIR)

- **Dataset F** (t0minidaq_endtop_scan_20260618_204959, EJ-204) — Confianza ALTA

### Datasets AMBIGUOS o de fase anterior a Phase 3

- **Dataset A** (scan_end_wrapped_2026-06-09, EJ-204) — AMBIGUO (Fase 2b o Fase 3)

### Datasets bajo Fases 5–7

**Ninguno** de los datasets en `t0minidaq` corresponde a las Fases 5 (SK-DD), 6 (2T-MYLAR con bug) o 7 (2T-MYLAR corregido). Confirmación consistente con `IMPACT_surface_model_2026-08-31.md §2`: no existen resultados comprometidos en git para Fases 4–7, y el único dataset de Fase 4 (Dataset F) tampoco fue comprometido.

---

## Vinculación figura → dataset

| EXEC / Presentación | Commits git que la generan | Dataset usado | Fase óptica | ¿Trazabilidad completa? |
|---------------------|---------------------------|---------------|-------------|------------------------|
| EXEC_03b / high-stats scan | `22adf23`, `9546c84` (2026-06-08) | Dataset A | **AMBIGUO** | **NO** — sin hash |
| EXEC_08 (photon budget) | `3e071c2`, `5d01447` | Dataset B | Fase 3 | Sí (geometría corrobora) |
| EXEC_08b (window dip) | `9613fc9` | Dataset C | Fase 3 | Sí |
| EXEC_09 (fit evidence + QA) | `2a9d57b`, `ea0c6d2` | Dataset B (probablemente) | Fase 3 | Sí (fecha coherente) |
| EXEC_13 EJ-204 | `c030e88`, `a7d22a3`, `b956102` | Dataset B (exec13_out_* viven dentro del dir de Dataset B) | Fase 3 | Sí |
| EXEC_13-230 | `5b93f4c` | Dataset D | Fase 3-equiv. feat/ej230-sslg4 | Parcial (sin hash) |
| EXEC_14 EJ-230 | `0b26c5f`, `f22de30`, `0b8be52`, `0e70266` | Dataset D | Fase 3-equiv. feat/ej230-sslg4 | Parcial (sin hash) |
| End-only Mylar EJ-204 | `489d238`, `fb3749d`, `c120c51` | Dataset E | feat/endonly-mylar @ `3ae135f` (git_dirty=4) | Parcial (hash sí, dirty) |
| "EndTop Scan 5000" (análisis en t0minidaq) | **ninguno comprometido** | Dataset F | **Fase 4 (GEN-1)** | **NO** — no existe en git |

---

## Figuras con trazabilidad rota o incompleta

### Trazabilidad rota (no recuperable desde git)

1. **"EndTop Scan 5000" / Dataset F** — Los análisis de σ_t de 5000 ev/pos × 31 posiciones (EJ-204, Fase 4 GEN-1) existen solo en `t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/analysis_corefit/` y `analysis_simple_std/`. **Ningún resultado derivado de este dataset aparece en un commit de git.** Si alguna figura del proyecto se generó manualmente a partir de estos análisis sin comprometer el origen, esa figura no tiene trazabilidad en el repositorio. Además, los resultados son **incorrectos para el observable de tesis** (σ_t del readout END) porque el modelo GEN-1 elimina TIR.

2. **Dataset A (scan_end_wrapped_2026-06-09)** — Sin hash de commit ni log de Geant4. Las figuras producidas a partir de EXEC_03b (commits `22adf23`, `9546c84`) no pueden asignarse con certeza a Phase 2b o Phase 3. La trazabilidad queda AMBIGUA.

### Trazabilidad incompleta (recuperable parcialmente)

3. **Dataset D (results_ej230)** — Sin hash de commit en logs. La asignación a Fase 3-equivalente es confiable por geometría, pero el commit exacto es desconocido. Los resultados de EXEC_13-230 y EXEC_14 son atribuibles a esta fase con confianza ALTA pero no se puede confirmar reproducibilidad exacta.

4. **Dataset E (endonly_mylar_20260614)** — Hash confirmado (`3ae135f`), pero `git_dirty=4` indica 4 archivos modificados en runtime. Los resultados no son exactamente reproducibles desde el commit limpio.

---

## Qué habría que regenerar

Este apartado enumera qué datasets necesitarían nueva simulación para obtener resultados válidos bajo el modelo óptico correcto (Fase 7: 2T-MYLAR con `GROUPVEL` corregido). No se evalúa aquí si la regeneración es necesaria — eso depende del alcance de la tesis. Solo se registra qué está afectado y por qué.

### Dataset F — Regeneración obligatoria para resultados de tesis

**El Dataset F (31 posiciones × 5000 ev/pos, EJ-204) fue producido bajo Fase 4 (GEN-1, sin TIR).** Este es el único escaneo de alta estadística en `feat/endtop-sslg4` y aparece ser el candidato para figuras de σ_t(x) en la tesis. Sus resultados de End readout son físicamente incorrectos bajo GEN-1: la eliminación de TIR reduce drasticamente la eficiencia de colección en posiciones centrales, produciendo métricas `LOW_EFF`/`LOW_STATS` en el análisis existente. Para ser válido para la tesis, este escaneo necesitaría regenerarse bajo Fase 7 (commit `5576687` o posterior en `feat/endtop-sslg4`).

### Dataset A — Regeneración condicional

**El Dataset A (EXEC_03b) tiene fase ambigua** y estadística limitada para el observable de tesis (los 31 ROOT files existen pero no hay log). Si este dataset fue o será utilizado en figuras de tesis, se necesitaría: (a) determinar la fase exacta mediante compilación y ejecución de prueba, o (b) regenerar bajo Fase 7.

### Datasets B, C — Investigación previa, no presentados como resultado final

Los Datasets B y C (Fase 3, 2000 ev/pos) son la base de EXEC_07–EXEC_09 y los análisis intermedios (exec13_manual_*.C en Dataset B). Si estos no son figuras de tesis sino exploración metodológica, no requieren regeneración. Si alguna figura de EXEC_13 EJ-204 o EXEC_08b se presenta en la tesis como resultado cuantitativo, requiere nota explícita de que corresponde a Fase 3 (Near-GEN-1, sin TIR en la barra).

### Dataset D — No regenerable en esta rama

Dataset D (EJ-230) pertenece a la rama `feat/ej230-sslg4`. Su regeneración es independiente de `feat/endtop-sslg4` y requiere actualizar esa rama a un modelo correcto, lo que está fuera del alcance de este diagnóstico.

### Dataset E — Reproducibilidad degradada

Dataset E puede ser recreado desde `3ae135f` en `feat/endonly-mylar`, pero los 4 archivos dirty introducen incertidumbre. Si los resultados del end-only Mylar se presentan en la tesis, conviene anotar la condición `git_dirty=4` y verificar si los archivos modificados afectaban parámetros de simulación.

---

## Notas de auditoría

- **t0minidaq como buffer temporal:** ninguno de los 6 datasets tiene resultados comprometidos en git bajo Fases 4–7. `t0minidaq` actúa como buffer de simulación; los resultados solo son recuperables mientras no se sobrescriban.
- **Ausencia de header estructurado en algunos logs:** los Datasets B, C y D usan logs de Geant4 que no incluyen el bloque `Branch / Commit` del master_scan.log (presente solo en Dataset F). Para futuros escaneos, se recomienda que el script de escaneo registre siempre el hash.
- **git_dirty como riesgo sistémico:** el `git_dirty=4` del Dataset E ilustra que los escaneos se lanzaron en ocasiones con código localmente modificado. Sin registro de los diff correspondientes, los resultados son parcialmente irreproducibles.

---

*Auditoría producida en `docs/branch-diagnosis-2026-08-31`. Basada en inspección directa de logs, metadata y estructura de directorios en t0minidaq (solo lectura). No se abrieron ficheros ROOT para análisis de contenido físico. No se estimaron magnitudes físicas; las categorías LOW_EFF/LOW_STATS citadas en §Dataset F provienen del fichero `ANALYSIS_CORE_FIT_REPORT.md` ya existente en t0minidaq, no de cálculo nuevo.*
