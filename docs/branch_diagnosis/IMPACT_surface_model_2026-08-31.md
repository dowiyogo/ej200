# Impacto del modelo de superficie — feat/endtop-sslg4 — 2026-08-31

**Author:** rrios (Claude Sonnet 4.6 assist)  
**Date:** 2026-08-31  
**Rama analizada:** `origin/feat/endtop-sslg4` (HEAD `5576687`)  
**REF_BASE:** `origin/feat/bar-end-vikuiti` (HEAD `b007750`)  
**Método:** solo lectura de git. No se abrieron archivos de datos ni se recorrió `/home/reriosto/SHiP/t0minidaq/`.

---

## §1 — Cronología del modelo de superficie en feat/endtop-sslg4

Commits que modifican `SetType`, `SetFinish`, `REFLECTIVITY`, el air gap, o la forma de aplicar la superficie (skin vs border) en `src/DetectorConstruction.cc` o `src/Materials.cc`.

### Tabla de fases

| Fase | Rango de fechas | Commit desencadenante | Archivo | Cambio |
|------|----------------|-----------------------|---------|--------|
| **0 — Sin superficie** | 2026-03-19 | `ed61b9a` Primer commit | DC.cc | Sin superficie óptica definida |
| **1 — Al-foil skin** | 2026-05-06 | `51ca4ff` add Al-foil reflective surface to wrap | DC.cc + Mat.cc | Primera superficie: `G4LogicalSkinSurface` en WrapLV, `dielectric_metal` |
| **2a — Skin WrapLV** | 2026-05-07 | `79bed03` replace border surface with LogicalSkinSurface on WrapLV | DC.cc | Skin surface sobre WrapLV |
| **2b — Border panels** | 2026-05-07 | `f07ea21` place BarPV in WorldLV with explicit reflector panels | DC.cc | BarPV en WorldLV; `G4LogicalBorderSurface` en panels (`foilHalfT` no especificado en este commit) |
| **3 — Near-GEN-1 SSLG4** | 2026-06-10 | `1f6aca1` add SSLG4 EndTop 86-channel geometry | DC.cc | Border surfaces en panels de `foilHalfT = 0.5 µm`, `CreateBarSkinReflector()` = `dielectric_metal` |
| **4 — GEN-1** | 2026-06-18 | `f39b84c` replace reflector volumes with bar skin surface | DC.cc | Elimina panels; aplica `G4LogicalSkinSurface("BarSkin", barLV, CreateBarSkinReflector())`; `CreateBarSkinReflector()` = `dielectric_metal` |
| **5 — SK-DD** | 2026-06-21 | `fc5d8b0` port exec22 TIR fix to EndTop build | Mat.cc | `CreateBarSkinReflector()` → `dielectric_dielectric` + REFLECTIVITY=0.95; DC.cc sigue usando skin en barLV; sin air gap explícito |
| **6 — 2T-MYLAR (gv-bug)** | 2026-08-12 | `4cb2c23` port EXEC_23 explicit air-gap + Mylar reflector geometry | DC.cc + Mat.cc | Elimina skin en barLV; añade `kAirGapThickness=0.10 mm` + `kReflectorThickness=0.05 mm`; `CreateBarSurface()` (bar→air) + `CreateMylarReflector()` (air→Mylar); group-velocity bug presente |
| **7 — 2T-MYLAR corrected** | 2026-08-14 | `5576687` eliminate group-velocity aliasing bug | DC.cc | Añade `GROUPVEL_air = c/n_bar` en worldMat; sin cambio en el modelo de superficie |

### Detalle por cambio relevante

**`f39b84c` (2026-06-18) — Fase 4 → GEN-1**

`src/DetectorConstruction.cc` (diff resumido):
- **Elimina** bloques `G4Box` para panels de lámina (`reflYMinusSolid`, `reflZSolid`, `reflXSolid`)
- **Elimina** todos los `G4LogicalBorderSurface(barPV → panelPV)`
- **Añade** `G4LogicalSkinSurface("BarSkin", barLV, Materials::CreateBarSkinReflector())`
- En Mat.cc vigente: `CreateBarSkinReflector()` devuelve `dielectric_metal`, Al-mirror espectral, `SetSigmaAlpha(0.0)`

**`fc5d8b0` (2026-06-21) — Fase 5 → SK-DD**

`src/Materials.cc` (diff resumido):
- `CreateBarSkinReflector()`: cambia `SetType(dielectric_metal)` → `SetType(dielectric_dielectric)`
- Elimina curva espectral Al-mirror; añade REFLECTIVITY=0.95 plana sobre [1.5–6.5] eV
- Mensaje de commit confirma: `EXEC_16-20` fueron simulados con `dielectric_metal` y requieren re-simulación
- DC.cc en este commit: sigue usando `G4LogicalSkinSurface("BarSkin", barLV, CreateBarSkinReflector())`; sin air gap

**`4cb2c23` (2026-08-12) — Fase 6 → 2T-MYLAR**

`src/DetectorConstruction.cc`:
- Añade `static constexpr G4double kAirGapThickness = 0.10 * mm` (línea 33)
- Añade `static constexpr G4double kReflectorThickness = 0.05 * mm` (línea 34)
- **Elimina** `G4LogicalSkinSurface("BarSkin", barLV, ...)`
- **Añade** loop de 3 panels (−Y, +Z, −Z): cada panel = air gap LV + Mylar LV
  - `auto* scintAirSurface = Materials::CreateBarSurface()` → `G4LogicalBorderSurface(barPV, airGapPV)` (línea 285)
  - `auto* airReflectorSurface = Materials::CreateMylarReflector()` → `G4LogicalBorderSurface(airGapPV, mylarPV)` (línea 286)
- `Materials::CreateMylarReflector()`: `dielectric_metal`, `unified`, `polished`, REFLECTIVITY=0.95, SPECULARLOBECONSTANT=0.0, SigmaAlpha=0.0

**`5576687` (2026-08-14) — Fase 7 → 2T-MYLAR corrected**

`src/DetectorConstruction.cc`:
- Añade propiedad `GROUPVEL` al worldMat: `c/1.58 = 18.97 cm/ns` (misma que barra)
- **No cambia** ningún parámetro de la superficie óptica
- Razón: bug de aliasing de velocidad de grupo en G4Transportation al cruzar el boundary bar/air-gap después de TIR; la corrección evita que el paso post-rebote use la velocidad del aire. Efecto en tiempo de llegada: ver `GROUP_VELOCITY_AUDIT.md §12` (referenciado en el commit)

---

## §2 — Artefactos versionados por fase

Todos los artefactos se extrajeron con `git log --diff-filter=A --name-only origin/feat/endtop-sslg4 -- '*.pdf' '*.csv' '*.root'`.

| Fecha de commit | Commit | Artefacto | Fase vigente | Modelo de superficie |
|----------------|--------|-----------|-------------|----------------------|
| 2026-03-19 | `ed61b9a` | `build/*.pdf` (10 PDFs), `build/photon_hits.root` | Fase 0 | Sin superficie óptica |
| 2026-06-09 | `2c8179e` | `presentations/sim_status_hi_2026-06-08/geometry_layout.pdf` | Fase 3 | Near-GEN-1 (`1f6aca1`) |
| 2026-06-10 | `9613fc9` | `analysis/exec07/exec08b_timing_gate.csv`, `exec08b_timing_gate_raw.csv`, `exec08b_window_dip_profiles.csv` | Fase 3 | Near-GEN-1 |
| 2026-06-10 | `2a9d57b` | `results/analysis_sigma_vs_x_2026-06-10/deep_dive/` (5 CSVs) | Fase 3 | Near-GEN-1 |
| 2026-06-10 | `5d01447` | `results/analysis_sigma_vs_x_2026-06-10/sigma_vs_x_End.csv` | Fase 3 | Near-GEN-1 |
| 2026-06-10 | `947b10f` | `scan_resume_2_audit.csv` | Fase 3 | Near-GEN-1 |
| 2026-06-11 | `3e071c2`–`a1439fd` | `analysis/exec07/exec07_photon_budget_report.pdf`, `exec09_*.pdf`, `exec10_*.pdf`, `exec11_*.pdf`, `exec12_report_full.pdf`, y 12 CSVs | Fase 3 | Near-GEN-1 |
| 2026-06-13 | `c030e88`–`b956102` | `analysis/exec13/exec13_report.pdf`, `exec13_*.csv` (6 CSVs), `figs/exec13_*.pdf` (11 PDFs) | Fase 3 | Near-GEN-1 |
| 2026-08-14 | `610b189` | `macros/validation_scan/*.mac` (7 macros), `scripts/analyze_validation.py`, `scripts/run_validation_scan.sh` | Fase 7 | 2T-MYLAR corrected |

**Conclusión crítica:**

> Ningún artefacto PDF/CSV/ROOT fue cometido durante la Fase 4 (GEN-1, 2026-06-18 a 2026-06-20), la Fase 5 (SK-DD, 2026-06-21 a 2026-08-11), ni la Fase 6 (2T-MYLAR con bug de group-velocity, 2026-08-12 a 2026-08-13).
>
> El commit de 2026-08-14 (`610b189`) aporta solo macros y scripts de simulación, sin resultados pre-computados. Los datos de EXEC_23 (Fase 7) no están versionados en git.
>
> **Por tanto, todos los artefactos PDF/CSV versionados en feat/endtop-sslg4 fueron producidos bajo el modelo Near-GEN-1 (Fase 3), salvo los de `build/` de 2026-03-19 (Fase 0).**

---

## §3 — Diferencia exacta entre feat/endtop-sslg4 y feat/bar-end-vikuiti

Comparación al HEAD de cada rama. Fuentes extraídas con `git show origin/<branch>:src/...`.

### Geometría de los panels (idéntica en ambas ramas)

| Parámetro | feat/endtop-sslg4 | feat/bar-end-vikuiti | Fuente |
|-----------|-------------------|---------------------|--------|
| Air gap thickness | 0.10 mm | 0.10 mm | `src/DetectorConstruction.cc:33` / `:35` |
| Reflector thickness | 0.05 mm | 0.05 mm | `src/DetectorConstruction.cc:34` / `:36` |
| Panels laterales | −Y, +Z, −Z (3 panels) | −Y, +Z, −Z (3 panels) | `src/DetectorConstruction.cc:285-350` / `:312-380` |
| Cara +Y | Top SiPMs (sin panel) | Top SiPMs (sin panel) | idem |
| Caras ±X | End SiPMs (sin panel) | End SiPMs (sin panel) | idem |
| Tipo de superficie aplicado | `G4LogicalBorderSurface` | `G4LogicalBorderSurface` | idem |

### Superficie bar→air gap (idéntica en ambas ramas)

Función: `Materials::CreateBarSurface()` — **código byte-a-byte idéntico** en ambas ramas.

| Parámetro | Valor | Fuente (ambas ramas) |
|-----------|-------|----------------------|
| Tipo | `dielectric_dielectric` | `src/Materials.cc:214` |
| Modelo | `unified` | `src/Materials.cc:215` |
| Finish | `polished` | `src/Materials.cc:216` |
| SigmaAlpha | 0.0 | `src/Materials.cc:217` |
| REFLECTIVITY | no establecida (Fresnel puro) | — |
| SPECULARLOBECONSTANT | no establecida | — |
| Aplicada como | border: BarPV ↔ AirGapPV | `src/DetectorConstruction.cc:285` / `:312` |

TIR aplica para fotones con θ > arcsin(1/1.58) = 39.3°. Idéntico en ambas ramas.

### Superficie air gap→reflector (DIFIERE)

| Parámetro | feat/endtop-sslg4 | feat/bar-end-vikuiti | Fuente endtop-sslg4 | Fuente REF_BASE |
|-----------|-------------------|--------------------|----------------------|-----------------|
| Función | `CreateMylarReflector()` | `CreateBarSkinReflector()` | `DC.cc:286` | `DC.cc:313` |
| **Tipo** | **`dielectric_metal`** | **`dielectric_dielectric`** | `Materials.cc:288` | `Materials.cc:347` |
| Modelo | `unified` | `unified` | `Materials.cc:289` | `Materials.cc:348` |
| Finish | `polished` | `polished` | `Materials.cc:290` | `Materials.cc:349` |
| SigmaAlpha | 0.0 (default `Materials.hh:42`) | 0.0 | `Materials.hh:42` | `Materials.cc:350` |
| REFLECTIVITY | 0.95 (default `Materials.hh:40`) | 0.95 (vector [1.5–6.5 eV]) | `Materials.hh:40` | `Materials.cc:357` |
| SPECULARLOBECONSTANT | 0.0 (default `Materials.hh:41`) | no establecida | `Materials.hh:41` | — |
| Material del panel | Mylar, n=1.65, 0.05 mm | Vikuiti 3M ESR, 0.05 mm | `DC.cc:219` (`CreateMylar()`) | `DC.cc:224` (Vikuiti material) |
| Aplicada como | border: AirGapPV ↔ MylarPV | border: AirGapPV ↔ VikuitiPV | `DC.cc:286` | `DC.cc:313` |

**Comportamiento físico de la diferencia:**

- `dielectric_metal`: Geant4 no aplica refracción ni ecuaciones de Fresnel en el boundary air→reflector. El fotón es absorbido o reflejado difusamente/especularmente según REFLECTIVITY y los parámetros unified model. No hay transmisión.
- `dielectric_dielectric`: Geant4 aplica ecuaciones de Fresnel. Un fotón que llega al panel Vikuiti puede reflejarse, transmitirse, o (si el ángulo es mayor que el TIR del panel) quedar retenido por TIR secundario. REFLECTIVITY=0.95 en `dielectric_dielectric` unified rige para los fotones que no hacen TIR secundario.

Esta diferencia afecta **solo a los fotones que no hicieron TIR en la interfaz bar→aire** (θ < 39.3° al llegar al air gap). Los fotones TIR no alcanzan el panel reflector. No se puede separar el efecto del tipo de superficie del efecto del material del reflector (Mylar vs Vikuiti) sin una simulación adicional con los modelos cruzados.

---

## §4 — Validez de comparaciones entre ramas

### Comparaciones INVÁLIDAS (confunden ≥2 variables)

| Comparación | Variables simultáneamente distintas |
|-------------|-------------------------------------|
| feat/endtop-sslg4 EXEC_23 (Fase 7) vs feat/bar-end-vikuiti | Tipo de superficie air→reflector (`dielectric_metal` vs `dielectric_dielectric`) **y** material del reflector (Mylar vs Vikuiti) |
| feat/endtop-sslg4 exec07-exec13 (Fase 3) vs feat/bar-end-vikuiti | Modelo de superficie (Near-GEN-1 vs 2T-VIK) **y** material del reflector **y** presencia/ausencia de air gap |
| feat/endtop-sslg4 exec07-exec13 (Fase 3) vs feat/endtop-sslg4 EXEC_23 (Fase 7, no versionado) | Modelo de superficie distinto dentro de la misma rama (Near-GEN-1 vs 2T-MYLAR corrected) |

### Comparaciones VÁLIDAS dentro de feat/endtop-sslg4

| Comparación | Condición |
|-------------|-----------|
| Entre posiciones dentro del mismo exec (exec07-exec13, Fase 3) | Mismo modelo Near-GEN-1; scintillator y readout fijos; válido para tendencias relativas σ_t(x) |
| Entre posiciones dentro de EXEC_23 (Fase 7, datos no versionados) | Mismo modelo 2T-MYLAR corrected; válido si los datos existen fuera de git |
| exec07 vs exec08b vs exec09 vs ... vs exec13 (comparación temporal de campañas) | Válido porque son todos Fase 3; método de análisis puede variar |

### Comparaciones VÁLIDAS entre ramas (bajo condiciones)

| Comparación | Condición de validez |
|-------------|----------------------|
| feat/endtop-sslg4 (Fase 3) vs otros GEN-1/Near-GEN-1 con EJ-204 (ej. exp/pair-scan-2026-06-11) | Mismo modelo óptico (Near-GEN-1 o GEN-1); admisible si ambos tienen el mismo defecto; diferencia de scintillator nula (ambos EJ-204); diferencia en readout o análisis puede existir |
| feat/endtop-sslg4 (Fase 3) vs feat/ej230-sslg4 (Fase GEN-1) | Scintillator difiere (EJ-204 vs EJ-230); modelo óptico también difiere (Near-GEN-1 vs GEN-1 skin). **No válida** sin declarar ambas diferencias |
| Componente TIR únicamente: feat/endtop-sslg4 (Fase 7) vs feat/bar-end-vikuiti | La interfaz bar→air es **idéntica** en ambas ramas (`dielectric_dielectric, polished`). Observables que dependan exclusivamente de TIR (e.g. tiempo de fotones reflejados por TIR antes del primer rebote en el panel) son comparables. Los observables que incluyen fotones que llegan al panel **no son comparables** |

### Nota sobre magnitudes

Este documento no estima la magnitud de los efectos sobre Npe o σ_t. Esa determinación requiere simulación con los modelos cruzados (Mylar con `dielectric_dielectric`, Vikuiti con `dielectric_metal`) y está fuera del alcance de un análisis de código.

---

## §5 — Resumen ejecutivo

1. **feat/endtop-sslg4 pasó por 7 fases ópticas distintas** entre 2026-03-19 y 2026-08-14. Los modelos Near-GEN-1, GEN-1, SK-DD, y 2T-MYLAR (con y sin bug de group-velocity) se sucedieron.

2. **Todos los artefactos versionados (PDF/CSV/ROOT) en la rama son de Fase 3 (Near-GEN-1)**, salvo los de `build/` (Fase 0 sin superficie). Ningún resultado de las Fases 4–7 está versionado en git.

3. **La diferencia crítica con REF_BASE al HEAD es una sola línea** (`CreateMylarReflector()` vs `CreateBarSkinReflector()`), que implica dos diferencias simultáneas: tipo de superficie (`dielectric_metal` vs `dielectric_dielectric`) y material del reflector (Mylar vs Vikuiti). No es posible atribuir ninguna diferencia en σ_t a un único factor sin simulación adicional.

4. **La interfaz bar→air es byte-a-byte idéntica** en ambas ramas; el TIR a 39.3° funciona igual.

5. **Las figuras de exec07-exec13 (Fase 3)** son internamente consistentes bajo Near-GEN-1 y pueden citarse como resultados con ese modelo, con la advertencia de §4. No son comparables directamente con los resultados de EXEC_23 (Fase 7) ni con feat/bar-end-vikuiti.

---

## §6 — Referencias de evidencia

```bash
# Cronología de cambios ópticos en DC.cc y Mat.cc
git log --follow --format='%h %ad %s' --date=short origin/feat/endtop-sslg4 -- src/DetectorConstruction.cc
git log --follow --format='%h %ad %s' --date=short origin/feat/endtop-sslg4 -- src/Materials.cc

# Patch de cada commit clave
git show f39b84c -- src/DetectorConstruction.cc   # GEN-1
git show fc5d8b0 -- src/Materials.cc               # exec22 SK-DD fix
git show 4cb2c23 -- src/DetectorConstruction.cc src/Materials.cc  # 2T-MYLAR
git show 5576687 -- src/DetectorConstruction.cc   # group-velocity fix

# Artefactos con fecha de commit
git log --format='%h %ad' --date=short --diff-filter=A --name-only origin/feat/endtop-sslg4 -- '*.pdf' '*.csv' '*.root'

# Parámetros de superficie exactos al HEAD
git show origin/feat/endtop-sslg4:src/Materials.cc | grep -n 'SetType\|SetFinish\|REFLECTIVITY\|SigmaAlpha'
git show origin/feat/bar-end-vikuiti:src/Materials.cc | grep -n 'SetType\|SetFinish\|REFLECTIVITY\|SigmaAlpha'
git show origin/feat/endtop-sslg4:include/Materials.hh | grep -A3 'CreateMylarReflector'
```
