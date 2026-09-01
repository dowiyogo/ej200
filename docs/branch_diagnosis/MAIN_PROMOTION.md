# MAIN_PROMOTION.md — Promoción de feat/bar-end-vikuiti a main

Fecha: 2026-09-01  
Ejecutado manualmente por el usuario. Este documento verifica y registra el estado.

---

## 1. Gate previo — verificación de ancestría

```bash
git rev-list --count feat/bar-end-vikuiti..main   # → 2
git merge-base --is-ancestor feat/bar-end-vikuiti main && echo "FF POSIBLE"
# → "feat IS ancestor of main"
```

`rev-list --count` devolvió **2** (no 0). Los dos commits en `main` que no estaban en
`feat/bar-end-vikuiti` eran:

| Commit | Descripción |
|--------|-------------|
| `87795b4` | merge: adopt Phase-7 optics and SSLG4 material framework from feat/bar-end-vikuiti |
| `84e902c` | fix(optics): replace reflector volumes with bar skin surface (GEN-1 — INVÁLIDO) |

El commit `84e902c` introduce la óptica GEN-1 (`G4LogicalSkinSurface` + `dielectric_metal`)
que hace inválido `main` previo. El merge adoptó el árbol `src/` e `include/` de
`feat/bar-end-vikuiti` explícitamente, sobrescribiendo los archivos conflictivos.

**Método usado (no PR):** merge directo con `--no-commit --no-ff`, seguido de
`git checkout feat/bar-end-vikuiti -- src/ include/` para resolver conflictos.
Commit de merge: `87795b4`.

Verificación post-merge:
```bash
git diff feat/bar-end-vikuiti main -- src/ include/   # → VACÍO
grep -rn "<<<<<<<\|>>>>>>>\|=======" src/ include/     # → ningún marcador de conflicto
```

---

## 2. Estado actual de main

```
git log --oneline -5 main
87795b4 merge: adopt Phase-7 optics and SSLG4 material framework from feat/bar-end-vikuiti
50cf02e docs(presentations): add master index presentations/README.md
1aa7931 chore(.gitignore): add build_*/ and missing LaTeX artifact patterns
0149240 docs(exec14): import QA docs and reflector-fix documentation from ej200_orchestrator
a53e1db fix(analysis/exec14): stage CSV outputs; add !*.csv exception to outputs/.gitignore
```

`feat/bar-end-vikuiti` apunta a `50cf02e` — es **ancestro directo** de `main`.

---

## 3. Verificación óptica post-merge (con archivo:línea)

### 3.1 Interfaz bar→aire

**`src/Materials.cc:214`** — `CreateBarSurface()`:
```cpp
surf->SetType(dielectric_dielectric);   // Materials.cc:214
surf->SetFinish(polished);              // Materials.cc:216
```
Usado en `src/DetectorConstruction.cc:312`:
```cpp
auto* scintAirSurface = Materials::CreateBarSurface();
```
Comentario de cabecera: `DC.cc:309`:
> `//   bar→air : dielectric_dielectric polished → TIR for θ > 39.3°`

**PASS** ✓ — `dielectric_dielectric`, `polished`; TIR activo.

### 3.2 Interfaz aire→reflector (Vikuiti)

**`src/Materials.cc:347`** — `CreateBarSkinReflector()`:
```cpp
surf->SetType(dielectric_dielectric);   // Materials.cc:347  ← NO dielectric_metal
surf->SetFinish(polished);              // Materials.cc:349
const std::vector<G4double> refl = {0.95, 0.95};  // Materials.cc:354
mpt->AddProperty("REFLECTIVITY", energy, refl);    // Materials.cc:357
```
Usado en `src/DetectorConstruction.cc:313`:
```cpp
auto* airReflectorSurface = Materials::CreateBarSkinReflector();
```
Aplicado a la superficie `AirGapToVikuiti_XXX` (`DC.cc:392`).

Comentario de cabecera: `DC.cc:310`:
> `//   air→Vik : dielectric_dielectric R=0.95 (CreateBarSkinReflector, Materials.cc:354)`

**PASS** ✓ — `dielectric_dielectric` (no `dielectric_metal`); REFLECTIVITY=0.95.

### 3.3 Air gap con RINDEX declarado

**`src/DetectorConstruction.cc:266`**:
```cpp
mpt->AddProperty("RINDEX",   {2.0*eV, 4.0*eV}, {1.0, 1.0});   // DC.cc:266
```
El volumen de aire recibe el `MaterialPropertiesTable` con RINDEX=1.0.

Geometría:
- `kAirGapThickness = 0.10 mm` (`DC.cc:35`)
- `kReflectorThickness = 0.05 mm` (`DC.cc:36`)

**PASS** ✓ — air gap explícito con RINDEX declarado.

### 3.4 GROUPVEL_air = c/n_bar

**`src/DetectorConstruction.cc:264–267`**:
```cpp
const G4double vg = CLHEP::c_light / 1.58;                      // DC.cc:264
mpt->AddProperty("GROUPVEL", {2.0*eV, 4.0*eV}, {vg, vg});       // DC.cc:267
```
`vg = c / 1.58 = 18.97 cm/ns`. Comentario adjunto (DC.cc:260–262):
> "Without GROUPVEL, Geant4 uses c instead of c/n, causing superluminal steps
> that corrupt timing. Set GROUPVEL = c/n_bar to avoid this."

**PASS** ✓ — GROUPVEL_air = c/n_bar = 18.97 cm/ns (Fase 7 correcta).

### Resumen de verificación óptica

| Criterio | Estado | Localización |
|----------|--------|-------------|
| bar→air: `dielectric_dielectric`, `polished` | **PASS** | `Materials.cc:214,216` via `DC.cc:312` |
| air→reflector: `dielectric_dielectric` (no `dielectric_metal`) | **PASS** | `Materials.cc:347` via `DC.cc:313,392` |
| Air gap explícito con RINDEX declarado | **PASS** | `DC.cc:35,266` |
| `GROUPVEL_air = c/n_bar = 18.97 cm/ns` | **PASS** | `DC.cc:264,267` |

**Los cuatro criterios pasan. La óptica GEN-1 de 84e902c fue sobrescrita correctamente.**

---

## 4. Pérdida de API respecto a 84e902c

El commit `84e902c` (GEN-1, óptica inválida) tenía una API de detector diferente para
los SiPMs Top:

| API eliminada | Presente en `84e902c` | Presente en `main` actual |
|---------------|----------------------|--------------------------|
| `ComputeNTopSiPMs()` | Sí (pitch dinámico: default 70 mm → 20 Top SiPMs) | **No** |
| `/det/topSiPMPitch` | Sí (macro command Geant4 UI) | **No** |
| `SetTopSiPMPitch()` | Sí | **No** |
| `GetNTopSiPMs()` | Sí | **No** |

`main` actual usa la API de `feat/bar-end-vikuiti`: número fijo de SiPMs Top
configurado en tiempo de compilación (variable `kNTopSiPMs` en `include/DetectorConstruction.hh`).

**Macros afectadas:** ninguna. Búsqueda en `macros/` de `/det/topSiPMPitch`, `pitch`,
`NTopSiPM` no devuelve ningún resultado. Los macros existentes usan `/det/readout` y
`/det/scintillator`, que sí están presentes en la API actual.

**Impacto:** cualquier script externo (no en este repositorio) que usara
`/det/topSiPMPitch` necesita actualizarse para fijar el número de Top SiPMs
en tiempo de compilación o vía otro mecanismo.

---

## 5. Default branch de GitHub

`main` ya era el default en GitHub antes del merge — no se modificó ningún ajuste
en Settings. Sin cambios necesarios.

---

## 6. Estado de los tres clones post-merge

```bash
git ls-remote origin refs/heads/main    # → 87795b48da7ba967a4914be5ba136374d7a69dd0
```

Verificación en W1 y W2 vía SSH (`ssh -p 9022 reriosto@127.0.0.1`):
```
W1 /home/reriosto/SHiP/ej200  → 87795b48da7ba967a4914be5ba136374d7a69dd0
W2 /home/reriosto/SHiP/ej230  → 87795b48da7ba967a4914be5ba136374d7a69dd0
```

**Los tres clones (HOST, W1, W2) ven `main` en `87795b4`. ✓**

---

## 7. Resumen

| Paso | Estado |
|------|--------|
| Gate previo (ancestría) | `feat` ES ancestro de `main` ✓ |
| Commits exclusivos de old main incorporados | `84e902c` GEN-1 sobrescrito por merge ✓ |
| Método de resolución de conflictos | `git checkout feat/bar-end-vikuiti -- src/ include/` ✓ |
| Merge commit | `87795b4` |
| Verificación óptica (4 criterios) | **PASS** ✓ |
| API loss (`ComputeNTopSiPMs`, `topSiPMPitch`) | Documentada — sin macros afectados ✓ |
| Default branch de GitHub | `main` (sin cambio necesario) ✓ |
| Sincronización de los 3 clones | `87795b4` en HOST, W1 y W2 ✓ |

*Documento de verificación — no re-ejecuta ningún paso del merge.*
