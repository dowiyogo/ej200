# BRANCH_TRIAGE_2026-09-01.md — Clasificación de ramas remotas contra main

**Fecha:** 2026-09-01  
**main:** `e7e1e35` (2T-VIK R=0.98, Phase 7 + SPTR + presentaciones consolidadas)  
**Referencia anterior:** `docs/branch_diagnosis/DIAGNOSIS_2026-08-31.md` (en rama `docs/branch-diagnosis-2026-08-31`), REF_BASE = `feat/bar-end-vikuiti @ b007750`  
**Ramas remotas totales:** 17 (incluyendo main)  
**Método:** `git rev-list --count`, `git diff --stat origin/main...origin/<B>`, `git cherry`, `git worktree list` en los tres clones. El diff manda; `cherry` falla ante squash y divergencia de contexto.

---

## 1. Tabla de clasificación

| Rama | Hash | Veredicto | Commits únicos vs main | Diff código (src/ include/) | Diff análisis | Diff docs/presentaciones | Worktree |
|------|------|-----------|------------------------|-----------------------------|---------------|--------------------------|----------|
| **main** | `e7e1e35` | REFERENCIA | — | — | — | — | HOST `/home/rrios/ej200`; W1 `ej200_edge_scan` |
| **feat/bar-end-vikuiti** | `1bfb827` | CONTENIDA EN MAIN ⚠ | 0 (21 detrás) | vacío | vacío | vacío | W1 `/mnt/d/SHiP/ej200` (checkout principal) |
| **diag/phase7-delta-2026-08-31** | `a18e863` | REQUIERE DECISIÓN | 6 | vacío | 26 archivos, +8636 líneas | 15 archivos, +4871 líneas | — |
| **docs/branch-diagnosis-2026-08-31** | `2a9645f` | REQUIERE DECISIÓN | 4 | vacío | vacío | 3 archivos, +1054 líneas | — |
| **exp/pair-scan-2026-06-11** | `f19c093` | ARCHIVAR | 30 | 2 archivos (−72/+9) | 14 archivos, +3839 líneas | vacío | — |
| **feat/bar-vikuiti** | `219fbe3` | ARCHIVAR | 36 | 2 archivos (−79/+54) | 17 archivos, +4737 líneas | vacío | — |
| **feat/ej204-bar-tir-only** | `09f8b18` | CONFIGURACIÓN ÚNICA | 32 | 2 archivos (−75/+20) | 15 archivos, +4226 líneas | vacío | — |
| **feat/ej204-event-display-tracks** | `47a9a4f` | FIJADA A WORKTREE | 4 | 7 archivos, +371/−73 | vacío | vacío | W1 `ej200_event_display` |
| **feat/ej228-cylinder** | `66b674c` | CONFIGURACIÓN ÚNICA | 3 | 7 archivos, +281/−575 | 2 archivos, +308 líneas | vacío | — |
| **feat/ej228-tir-only** | `0006919` | CONFIGURACIÓN ÚNICA | 4 | 7 archivos, +264/−575 | 6 archivos, +1514 líneas | vacío | — |
| **feat/ej230-bar-tir-only** | `b281aea` | CONFIGURACIÓN ÚNICA | 34 | 2 archivos (−79/+56) | 15 archivos, +4226 líneas | vacío | — |
| **feat/ej230-endonly-mylar** | `04a8047` | CONFIGURACIÓN ÚNICA | 42 | 13 archivos, +149/−221 | 42 archivos, +6490 líneas | 10 archivos, +692 líneas | W2 `ej230_endonly_mylar` |
| **feat/ej230-sslg4** | `5b93f4c` | FIJADA A WORKTREE | 34 | 2 archivos (−78/+49) | 75 archivos, +9490 líneas | 10 archivos, +692 líneas | W2 `/mnt/d/SHiP/ej230` (checkout principal) |
| **feat/endonly-mylar** | `fb3749d` | FIJADA A WORKTREE | 23 | 5 archivos, +149/−81 | 10 archivos, +2571 líneas | 1 archivo, +243 líneas | W1 `ej200_endonly` |
| **feat/endtop-sslg4** | `5576687` | REQUIERE DECISIÓN | 2 | 1 archivo, +17/−5 | 2 archivos, +308 líneas | vacío | — |
| **feature/sipm-electronics-response** | `bb4cf6a` | CONFIGURACIÓN ÚNICA | 15 | 8 archivos, +566/−110 | 24 archivos, +5304 líneas | vacío | — |
| **wip/host-stash-endtop-junio** | `9710f34` | ARCHIVAR | 1 | vacío | 1 archivo, +3/−3 | vacío | — |
| **wip/host-uncommitted-2026-08-31** | `d2c8d4c` | CONFIGURACIÓN ÚNICA | 2 | 3 archivos, +90/−15 | 1 archivo, +102 líneas | vacío | — |

**Nota sobre el diff de código:** Las columnas con valores negativos altos (ej. −575 en ramas ej228) no indican código eliminado por la rama — indican que la rama tiene una base común anterior a los commits de infra de main y el diff acumula las adiciones de main desde ese punto de divergencia.

---

## 2. Ramas con CONFIGURACIÓN ÚNICA — detalle

### feat/ej204-bar-tir-only @ 09f8b18
**Geometría:** Bar 1400×60×10 mm · **Centelleador:** OPSC-101 (EJ-204) · **Readout:** End + Top + EndTop  
**Superficie:** TIR-ONLY — `G4LogicalSkinSurface` en barLV con `dielectric_dielectric` polished. TIR opera correctamente (θ_c=39.3°); los fotones que fallan TIR se pierden al mundo. Sin reflector secundario.  
**Qué aporta:** 32 commits con scripts de pair-scan y análisis de timing EJ-204 sin reflector. Único para estudiar TIR puro como cota inferior del rendimiento (comparación con main 2T-VIK cuantifica el aporte del reflector).  
**Commits únicos:** par scan macros, scripts/run_pair_scan.sh, scripts/run_resolution_scan.sh, análisis de geometría.

### feat/ej228-cylinder @ 66b674c
**Geometría:** Cilindro ∅25 × 25 mm (G4Tubs) · **Centelleador:** EJ-228 (hardcoded) · **Readout:** 8 SiPMs en las caras cap (4 top + 4 bottom)  
**Superficie:** `dielectric_metal` R=0.98 (wrap anular, gap 0.10 mm). Nota: usa `dielectric_metal` — distinto de main que usa `dielectric_dielectric`. El valor R=0.98 coincide ahora con main tras `c7acb7a`.  
**Qué aporta:** única geometría cilíndrica en el repo. Commits: `feat(ej228)`, GROUPVEL fix (5576687), validation scan.

### feat/ej228-tir-only @ 0006919
**Geometría:** Mismo cilindro ∅25 × 25 mm + EJ-228 · **Superficie:** TIR-ONLY (gap 0.10 mm, sin wrap).  
**Qué aporta:** comparación sistemática TIR-only vs 2T-VIK para el cilindro. 4 commits propios: análisis Beamer comparativo + commit de cilindro heredado de feat/ej228-cylinder.

### feat/ej230-bar-tir-only @ b281aea
**Geometría:** Bar 1400×60×10 mm · **Centelleador:** OPSC-106 (EJ-230) · **Superficie:** TIR-ONLY.  
**Qué aporta:** contraparte EJ-230 del pair-scan TIR-only; 34 commits de análisis. Sistemática para cuantificar contribución del reflector en EJ-230 comparada con main 2T-VIK R=0.98.

### feat/ej230-endonly-mylar @ 04a8047
**Geometría:** Bar 1400×60×10 mm · **Centelleador:** OPSC-106 (EJ-230) · **Readout:** End only · **Superficie:** Near-GEN-1 (`dielectric_metal`, Mylar panels directamente en contacto, sin air gap)  
**Novedades desde DIAGNOSIS:** +1 commit `04a8047` "fix(beamer): manual layout corrections" — edición Beamer, sin cambio físico.  
**Worktree:** W2 `ej230_endonly_mylar` (activo — no borrable sin desmontar).  
**Qué aporta:** rama de EJ-230 End-only más completa: 42 commits únicos, 42 scripts de análisis propios (`exec13`, `exec14b`, `generate_exec14b_tables.py`, `audit_beamer_assets.py`), Beamer EJ-230. El modelo Near-GEN-1 afecta el Npe absoluto pero no las comparaciones relativas internas a la rama.

### feature/sipm-electronics-response @ bb4cf6a
**Geometría:** Bar 1400×60×10 mm · **Centelleador:** EJ-200/EJ-204 (configurable por nombre, sin OPSC codes) · **Superficie:** 3T-MYLAR — bar→Mylar `dielectric_dielectric` (Fresnel), Mylar→air TIR a 37.3°; Mylar wrap 25 µm. Sin air gap explícito. `fEdgeWrapMode` configurable (Mylar/air/black para caps laterales).  
**Modelo de electrónica:** cadena SiPM única en el repo: `fTopSiPMPitch` (configura el layout dinámicamente, distinto al número fijo de main).  
**Qué aporta:** 15 commits con la única implementación de respuesta electrónica SiPM y el único modelo de reflector 3T. No actualizado al sistema OPSC. 8 archivos en src/include distintos de main.

### wip/host-uncommitted-2026-08-31 @ d2c8d4c
**Geometría:** Bar 1400×60×10 mm · **Centelleador:** OPSC-101 (EJ-204) · **Superficie:** 2T-VIK (misma que main)  
**Modo único:** `EndSparseTop` — añade `/det/nSparseTopSiPMs N` (hasta 200 SiPMs uniformemente distribuidos); función `SparseTopSiPMCenterX()` en DetectorConstruction.hh/cc; `SiPMSD.cc` actualizado para clasificarlos.  
**Commits únicos:** `d2c8d4c` (run script), `1f75563` (código DetectorConstruction + SiPMSD).  
**Defectos activos** (de DIAGNOSIS_2026-08-31 §8, sin cambio):
- Defecto #1: gun midpoint usa `TopSiPMCenterX()` en lugar de `SparseTopSiPMCenterX()` en modo EndSparseTop.
- Defecto #2: SiPMs sparse colocados a `kBarHalfY - 2·kTopHalfY` (0.25 mm más profundo que los densos).

---

## 3. Caso especial — feat/bar-end-vikuiti

### Estado
`origin/feat/bar-end-vikuiti` apunta a `1bfb827` (commit anterior al merge commit `87795b4` que adoptó su contenido en main).

```
git rev-list --count origin/main..origin/feat/bar-end-vikuiti  → 0
git rev-list --count origin/feat/bar-end-vikuiti..origin/main  → 21
```

Los 21 commits que main tiene por encima de `feat/bar-end-vikuiti`:

| Hash | Descripción |
|------|-------------|
| `e7e1e35` | fix(talk_v6/I1): SPTR corregido (Lee 2025) |
| `4851c43` | docs: SPTR_PROVENANCE.md |
| `8c308a5` | docs: defecto R=0.95 en 9 READMEs |
| `a791603` | docs: REFLECTIVITY_CHANGE.md; REFLECTOR_LABEL superado |
| `c7acb7a` | fix(optics): R=0.95→0.98 en CreateBarSkinReflector |
| `b17b11e` | docs: MAIN_PROMOTION.md |
| `87795b4` | merge: adopta Phase-7 de feat/bar-end-vikuiti |
| `50cf02e` | docs: presentations/README.md master index |
| `1aa7931` | chore: .gitignore |
| `0149240` | docs: exec14 QA + reflector_fix |
| `a53e1db` | fix: CSV outputs exec14 |
| `4e4502e` | feat: exec14 sidecar framework |
| `14d1699` | feat: exec14 deck |
| `12d6dc3` | docs: 7 READMEs presentations/ |
| `313782f` | feat: napkin_first_principles movido |
| `04e676d` | feat: presentation_v6 movida |
| `0b91983` | feat: presentation_v5 movida |
| `0f64634` | feat: presentation_v4 movida |
| `a24212f` | docs: ORCHESTRATOR_AUDIT.md |
| `1361d83` | docs: SYNC_STATUS.md |
| `84e902c` | fix(optics): GEN-1 (commit exclusivo de old main, sobrescrito por el merge) |

### Implicación operativa

**W1 tiene `feat/bar-end-vikuiti` como checkout principal** (`/mnt/d/SHiP/ej200` HEAD = `1bfb827`). Quien trabaje desde W1 verá:
- Sin la corrección de reflectividad R=0.98 (`c7acb7a`)
- Sin las entradas de defecto R=0.95 en los READMEs
- Sin SPTR_PROVENANCE.md ni la diapositiva I1 corregida
- Sin presentations/ consolidado (las presentaciones siguen en `analysis/presentation_v4/v5/v6/`)

**Riesgo:** iniciar nuevos commits desde W1 sobre `feat/bar-end-vikuiti` crearía divergencia con main. Si la rama recibe nuevos pushes sin ff-forward previo, forzará un merge o un push rechazado.

**Qué implica para el flujo:** W1 debe hacer `git fetch && git merge --ff-only origin/main` (o checkout de main) antes de cualquier nuevo desarrollo. No es automático — requiere acción manual en W1.

**Veredicto:** `feat/bar-end-vikuiti` no aporta contenido que main no tenga. Está detrás 21 commits. La rama en sí misma es CONTENIDA EN MAIN. El riesgo operativo es el checkout de W1, no la existencia de la rama.

---

## 4. Nota sobre tags en origin

La instrucción de triage decía "origin no tiene ningún tag". **Esta premisa es incorrecta:**

```
git ls-remote --tags origin
```

Resultado:
```
checkpoint/pre-endtop-sslg4-2026-06-10    2a9d57b
checkpoint/pre-exec11-2026-06-11          75966975  (dos refs apuntan al mismo objeto)
checkpoint/pre-exec11b-2026-06-11         75966975
checkpoint/pre-exec12-beamer-2026-06-11   6e8d3d6e
checkpoint/pre-exec12b-2026-06-11         f767be87
```

Origin tiene al menos 5 tags `checkpoint/*` del ciclo de sesiones 2026-06-10/11. No son tags de archivo de ramas sino snapshots previos a ejecuciones de análisis.

**Consecuencia para el archivado:** para el archivado de ramas sí se requiere crear tags nuevos del estilo `v-archive/<branch-slug>` (que no existen aún) y pushearlos antes de borrar la rama. La infra de tags existe pero los tags de archivo hay que crearlos expresamente.

---

## 5. Ramas con veredicto REQUIERE DECISIÓN — detalle

### diag/phase7-delta-2026-08-31 @ a18e863
**Commits únicos:** 6 — toda la diagnosis 2026-08-31 + scripts de análisis.  
**Código:** ninguna diferencia en src/include.  
**Análisis:** `analysis/optim/` (make_summary_plots.py, phase_ab/cd/e scripts, root_best_est/), `analysis/presentation_v4/v5/v6/` en su ubicación original (antes de la migración P2 a `presentations/`).  
**Docs:** DIAGNOSIS_2026-08-31.md, DATA_AUDIT_2026-08-31.md, IMPACT_surface_model_2026-08-31.md, PHASE7_DELTA_2026-08-31.md — ninguno de estos está en main.  
**Decisión pendiente:** ¿importar los scripts de análisis de optim/ a main? ¿importar las versiones antiguas de presentation_v4/v5/v6/? Los docs de diagnóstico tienen valor histórico; la diagnosis actual (`BRANCH_TRIAGE_2026-09-01.md`) los actualiza pero no los reemplaza.

### docs/branch-diagnosis-2026-08-31 @ 2a9645f
**Commits únicos:** 4 — subconjunto de diag/phase7-delta (solo los 4 docs, sin los scripts de análisis).  
**Docs:** DIAGNOSIS_2026-08-31.md, DATA_AUDIT_2026-08-31.md, IMPACT_surface_model_2026-08-31.md.  
**Decisión pendiente:** el DIAGNOSIS_2026-08-31.md es el informe de referencia citado por este triage. ¿Importar a main? Si se importan los docs, esta rama queda redundante con diag/phase7-delta.

### feat/endtop-sslg4 @ 5576687
**Commits únicos:** 2.  
- `5576687 fix(optics): eliminate group-velocity aliasing bug` — el fix de GROUPVEL está conceptualmente absorbido en main (main tiene `GROUPVEL_air = c/n_bar` desde Phase 7), pero el hash difiere (contexto distinto).  
- `610b189 feat(validation): physics validation scan — 7 pos × 500 events` — datos de validación únicos de 7 posiciones.  
**Configuración:** 2T-MYLAR (`dielectric_metal` para air→Mylar) — físicamente distinta de main 2T-VIK. El código diff en DetectorConstruction.cc (17 ins/5 del) refleja esta diferencia.  
**Decisión pendiente:** ¿conservar como referencia histórica del modelo 2T-MYLAR? ¿archivar? Los datos de validación (610b189) son el único activo tangible no absorbido.

---

## 6. Ramas FIJADAS A WORKTREE — estado actual

No borrables desde origin mientras el worktree esté montado.

| Rama | Worktree | Checkout HEAD | Advertencia |
|------|----------|---------------|-------------|
| `feat/ej204-event-display-tracks` | W1 `ej200_event_display` | `47a9a4f` | Código único (DisplayTrackingAction); cherry-pick antes de archivar |
| `feat/ej230-sslg4` | W2 `/mnt/d/SHiP/ej230` (principal) | `5b93f4c` | GEN-1; superficie inválida; desmontar W2 antes de borrar |
| `feat/endonly-mylar` | W1 `ej200_endonly` | `fb3749d` | GEN-1; desmontar W1 `ej200_endonly` antes de borrar |

**Nota sobre W1 `feat/bar-end-vikuiti`:** el checkout principal de W1 (`/mnt/d/SHiP/ej200`) está en `feat/bar-end-vikuiti` @ `1bfb827`, no en una worktree vinculada. Este checkout queda 21 commits detrás de main — ver §3.

---

## 7. Correspondencia con DIAGNOSIS_2026-08-31.md

La tabla del §7 de DIAGNOSIS_2026-08-31 usaba REF_BASE = `feat/bar-end-vikuiti`. A continuación los cambios de veredicto con main como nueva referencia.

| Rama | Veredicto en DIAGNOSIS-08-31 | Veredicto ahora | Cambio |
|------|------------------------------|-----------------|--------|
| main | REQUIERE DECISIÓN (GEN-1 bug) | REFERENCIA | ✅ GEN-1 resuelto; main es la nueva base |
| feat/bar-end-vikuiti | REF_BASE | CONTENIDA EN MAIN ⚠ | 21 commits detrás; riesgo en W1 |
| diag/phase7-delta-2026-08-31 | (no existía) | REQUIERE DECISIÓN | Nueva rama desde esa sesión |
| docs/branch-diagnosis-2026-08-31 | (no existía) | REQUIERE DECISIÓN | Nueva rama desde esa sesión |
| feat/bar-vikuiti | ARCHIVAR | ARCHIVAR | Sin cambio |
| feat/ej230-bar-tir-only | REQUIERE DECISIÓN | CONFIGURACIÓN ÚNICA | TIR-only es config. válida y única |
| feat/ej204-bar-tir-only | REQUIERE DECISIÓN | CONFIGURACIÓN ÚNICA | Idem, variante EJ-204 |
| feat/ej228-cylinder | MANTENER VIVA | CONFIGURACIÓN ÚNICA | Geometría cilíndrica única |
| feat/ej228-tir-only | MANTENER VIVA | CONFIGURACIÓN ÚNICA | Contraparte TIR-only del cilindro |
| feat/ej230-endonly-mylar | MANTENER VIVA | CONFIGURACIÓN ÚNICA | EJ-230 End-only + análisis EXEC14b |
| feat/ej230-sslg4 | ARCHIVAR | FIJADA A WORKTREE | W2 checkout; imposible borrar ahora |
| feat/endonly-mylar | ARCHIVAR | FIJADA A WORKTREE | W1 worktree; imposible borrar ahora |
| feat/ej204-event-display-tracks | ARCHIVAR | FIJADA A WORKTREE | W1 worktree; código único pendiente de cherry-pick |
| feat/endtop-sslg4 | REQUIERE DECISIÓN | REQUIERE DECISIÓN | Sin cambio; 2T-MYLAR distinto de main |
| feature/sipm-electronics-response | REQUIERE DECISIÓN | CONFIGURACIÓN ÚNICA | Electrónica + 3T-MYLAR únicos |
| wip/host-uncommitted-2026-08-31 | MANTENER VIVA | CONFIGURACIÓN ÚNICA | EndSparseTop único |
| wip/host-stash-endtop-junio | ARCHIVAR | ARCHIVAR | Sin cambio; 1 commit trivial |
| exp/pair-scan-2026-06-11 | ARCHIVAR | ARCHIVAR | Sin cambio; GEN-1 |
