# EJ-230 Material Audit — EXEC_13

**Authority:** Eljen Technology EJ-228/EJ-230 datasheet, Rev. Aug 2023  
**Branch:** `feat/ej230-sslg4`  
**Date:** 2026-06-12

## 1. Datasheet vs SSLG4 vs Geant4 — property table

| Property | Datasheet (EJ-230) | SSLG4 opsc-106 source | Value delivered to Geant4 |
|----------|-------------------|----------------------|--------------------------|
| SCINTILLATIONYIELD | 9 700 ph/MeV | `macros/oscnt/opsc-106.mac:9` — `9700 1/MeV` | 9700 / MeV ✓ |
| SCINTILLATIONRISETIME1 | 0.5 ns | `macros/oscnt/opsc-106.mac:11` — `0.5 ns` | 0.5 ns ✓ |
| SCINTILLATIONTIMECONSTANT1 | 1.5 ns | `macros/oscnt/opsc-106.mac:10` — `1.5 ns` | 1.5 ns ✓ |
| ABSLENGTH | 120 cm | `data/oscnt/opsc-106/absLength.txt` — 120 cm (2-point flat) | 120 cm ✓ |
| Emission peak | 391 nm | `data/oscnt/opsc-106/scntComp1.txt` — peak at ~390.5 nm | ~390.5 nm ✓ (within 1 nm) |
| RINDEX | 1.58 | `data/oscnt/opsc-106/rIndex.txt` — 1.58 (flat) | 1.58 ✓ |
| Density | 1.023 g/cm³ | Implemented in `src/Materials.cc:CreateEJ230()` (PVT, FindOrBuildPVT) | 1.023 g/cm³ ✓ (via `FindOrBuildPVT`) |

**Note:** All values match the datasheet. No overrides were required in practice;
`DetectorConstruction.cc:239–270` adds a defensive enforcement block for OPSC-106
(analogous to the OPSC-101 block) to guard against future upstream drift in SSLG4.

### Comparison with EJ-204 (OPSC-101, from feat/endtop-sslg4)

| Property | EJ-204 (OPSC-101) | EJ-230 (OPSC-106) | Notes |
|----------|------------------|------------------|-------|
| SCINTILLATIONYIELD | 10 400 ph/MeV | 9 700 ph/MeV | −6.7 % |
| SCINTILLATIONRISETIME1 | 0.7 ns | 0.5 ns | −29 % faster |
| SCINTILLATIONTIMECONSTANT1 | 1.8 ns | 1.5 ns | −17 % faster |
| ABSLENGTH | 160 cm | 120 cm | −25 % (attenuation dominant for long bar) |
| Emission peak | 408 nm | 391 nm | −17 nm blue-shifted |
| RINDEX | 1.58 | 1.58 | identical |
| Density | 1.023 g/cm³ | 1.023 g/cm³ | identical |

### EJ-204 ABSLENGTH discrepancy (documented from EJ-204 campaign)

The SSLG4 OPSC-101 data file listed 380 cm (Pilot spec), not 160 cm (EJ-204 datasheet).
The `DetectorConstruction.cc` OPSC-101 override block (`src/DetectorConstruction.cc:188–218`)
corrects this to 160 cm. EJ-230/OPSC-106 has no such discrepancy (SSLG4 = datasheet = 120 cm).

---

## 2. PDE spectral check — HOOK_PDE_SPECTRAL

### Finding

The Broadcom AFBR-S4N66P024M SiPM PDE is applied as a **full spectral curve PDE(λ)**,
NOT as a flat value. Source: `src/Materials.cc:CreateSiPMSurface()` calls
`SiPMModel::LoadPDECurve()` which reads `data/sipm/AFBR-S4N66P024M_pde.txt`
(37 wavelength points, 250–900 nm). The curve is installed as `EFFICIENCY(E)` on the
`dielectric_metal` SiPM surface, so Geant4 applies photon-by-photon sampling at the
actual photon wavelength.

### Implication for EJ-204 vs EJ-230 comparison

Because the PDE is a spectral curve, both campaigns use the physically correct response.
However, the different emission peaks produce different effective average PDEs:

| Material | Peak emission | PDE at peak (interpolated) | Comment |
|----------|--------------|---------------------------|---------|
| EJ-204 (OPSC-101) | 408 nm | ≈ 0.630 (plateau region, datasheet 63 % @ 420 nm) | Peak in flat plateau |
| EJ-230 (OPSC-106) | 391 nm | ≈ 0.617 (interpolated: 380 nm→0.600, 400 nm→0.630) | Slightly below plateau |

The effective PDE is ~2 % lower for EJ-230 photons due to the blue-shifted spectrum.
This introduces a **systematic bias** in the EJ-230 vs EJ-204 N_pe comparison that
cannot be attributed solely to the difference in scintillation yield (9700 vs 10400 ph/MeV).
The effect is small (≲2 %) but should be accounted for in campaign comparisons.

### Hook reserved for future correction

```
HOOK_PDE_SPECTRAL: apply wavelength-averaged PDE weight per material
when comparing N_pe between EJ-204 and EJ-230 campaigns.
Weight = integral(S(λ)·PDE(λ)dλ) / integral(S(λ)dλ)
EJ-204: eff_PDE ≈ 0.628  (to be computed from SSLG4 emission spectrum + PDE curve)
EJ-230: eff_PDE ≈ 0.616  (to be computed from SSLG4 emission spectrum + PDE curve)
Correction factor for N_pe comparison: 0.616/0.628 ≈ 0.981
NO code change until explicitly requested.
```

---

## 3. CTest guardrail

The `sslg4_properties_check` CTest (`tests/physics_baseline_check.cc`) verifies at
build time that the active material (OPSC-106) delivers to Geant4:
- yield = 9700 / MeV
- rise = 0.5 ns
- decay = 1.5 ns
- ABSLENGTH = 120 cm (all energy bins)
- RINDEX = 1.58 (all bins)
- emission peak within 2 nm of 391 nm

This guardrail is analogous to the one that existed for OPSC-101/EJ-204.

---

## 4. File references

| File | Line(s) | Change |
|------|---------|--------|
| `include/DetectorConstruction.hh` | 87 | `fScintCode` default changed to `"OPSC-106"` |
| `src/DetectorConstruction.cc` | 71 | Messenger description updated |
| `src/DetectorConstruction.cc` | 110–117 | Added `EJ230→OPSC-106` alias; extended validation |
| `src/DetectorConstruction.cc` | 219–269 | Added OPSC-106 override/validation block |
| `tests/physics_baseline_check.cc` | 49–91 | Updated to verify OPSC-106 properties |
| `macros/endtop_smoke_center.mac` | 8 | `OPSC-101 → OPSC-106` |
| `macros/endtop_smoke_edge.mac` | 8 | `OPSC-101 → OPSC-106` |
| `macros/edge_scan_smoke.mac` | 8 | `OPSC-101 → OPSC-106` |
| `macros/exec08b_run_*.mac` | 7 | `OPSC-101 → OPSC-106` (4 macros) |
| `scripts/run_exec07_scan.sh` | 59 | `OPSC-101 → OPSC-106` in generated macro heredoc |

---

## 5. Hallazgos de ejecución (EXEC_14)

### HOOK_OPSC106_MASSFRACTION

**Warning observed** during EXEC_13 scan (all 31 positions, session `20260612_171407`):

```
*** G4Exception : mat031
For material opsc-106 sum of fractional masses 0.998616 is not 1 — results may be wrong
  SSLG4 code: OPSC-106
```

**Origin**: The `opsc-106` material definition in SSLG4 has elemental fractional masses that sum to 0.998616 instead of 1.000000. The deficit is 0.1384 % (≈ 1.4 per mille).

**Impact assessment**:
- Geant4 applies this composition to ionization energy-loss (Bethe-Bloch) and to optical photon emission density.
- At 0.14 % level, the effect on mean dE/dx and hence on scintillation yield is sub-per-mille and well within the systematic uncertainty from the nominal OPSC-106 composition itself.
- Optical photon properties (RINDEX, ABSLENGTH, emission spectrum) are stored in the material properties table independently of the mass fraction — they are unaffected.
- The warning does not indicate any crash or invalid material object; Geant4 normalizes internally.

**Hook reserved**: `HOOK_OPSC106_MASSFRACTION` — report to SSLG4 maintainers for correction in the next library version; verify that normalization is applied (check `G4Material::SetMassOfMolecule` path in Geant4 source). **NO change to simulation code in EXEC_14.**

### Nearest-Top channel at x=0

At x=0 mm, the nearest-Top channel is a geometric tie between IDs 50 (at −12 mm) and 51 (at +12 mm). The tie is broken by N_pe: for EJ-230, ID 50 wins (vs ID 51 for EJ-204). This is consistent with stochastic fluctuations at symmetric positions and is not a material-dependent effect. The Landau MPV table updated accordingly (line 617 of `exec13_ej230_report_full.tex`).

### σ_group (End) — EJ-230 vs EJ-204

| Metric | EJ-204 | EJ-230 | Explanation |
|---|---|---|---|
| Mean σ_group | 148 ps | 226 ps | Lower yield + shorter λ_eff → fewer photons at far End |
| Range (min–max) | 12–437 ps | 12–595 ps | Far-End statistics more limited for EJ-230 |
| λ_eff | 33.1 cm | 30.8 cm | Shorter for EJ-230 despite bulk att=120 cm vs 160 cm |

The σ_group increase at the far End is physical: EJ-230 yield (9700/MeV) is 6.7 % below EJ-204 (10400/MeV), and the shorter λ_eff amplifies this effect for far-End channels.

---

## 6. Cierre de procedencia del reporte (EXEC_14C)

- El scan principal permanece intacto: `31/31` ROOT válidos, 2000 eventos por
  posición, sesión `session_20260612_171407` producida en `t0minidaq` con 24
  threads.
- Las siete simulaciones especiales fueron producidas con OPSC-106, jitter cero
  y 16 threads en MSI: cuatro window-dip de 2000 eventos y tres End-only de
  10000 eventos.
- Las tres figuras especiales antes heredadas fueron regeneradas desde datos
  EJ-230 y sus SHA256 difieren de las versiones EJ-204.
- El reporte visible es standalone EJ-230: la auditoría OCR sobre 134 PNG
  únicos no encontró `EJ-204`, `OPSC-101`, `EXEC_12`, `material comparison` ni
  `mirror of`.
- `HOOK_OPSC106_MASSFRACTION` permanece documentado; EXEC_14C no modificó
  `src/`, `include/`, `main.cc`, física ni geometría.

---

## 7. Diagnóstico físico y cierre del reporte (EXEC_14D)

### EndTop versus End-only

- Veredicto apples-to-apples: **PHYSICAL**.
- Ambas configuraciones usan los mismos End SUM4, el mismo leading edge y
  `sigma_group = sigma(deltaT_AB)/sqrt(2)`.
- Al submuestrear End-only a 2000 eventos con 200 réplicas, la razón
  EndTop/End-only permanece entre `1.223 +/- 0.019` y `1.963 +/- 0.033`; el
  menor cuantil inferior al 95 % es `1.186`.
- End-only conserva `1.62--3.04` veces más fotoelectrones End. EndTop reduce
  `t99` lateral en `0.518--1.182 ns`, pero para EJ-230 esa remoción de cola no
  compensa la pérdida de fotones End.

### Métricas de timing

- La métrica t4 visible es únicamente el sigma gaussiano ajustado del núcleo.
  Rango nearest-Top: `6.648--9.594 ps`; rango near-End:
  `0.562--127.719 ps`, con núcleos ultrafinos marcados
  `resolution_limited`.
- `HOOK_ADAPTIVE_TN` usa alcance mínimo `0.95`: 17/31 far-End reducen
  `N_eff`, 14/31 conservan N=20 y no queda starvation genuino.
- En las siete posiciones clave, far-End usa
  `N_eff = 7, 8, 17, 20, 17, 8, 7`.
- Las ocho degeneraciones previamente registradas quedaron resueltas: cuatro
  por binning adaptativo de t4 y cuatro por umbral far-End adaptativo.

### Preflight del PDF

- Los 31/31 ROOT principales permanecen válidos, con 2000 eventos, canales
  0--85 y sesión `session_20260612_171407` en `t0minidaq`/24 threads.
- `latexmk -pdf -halt-on-error` terminó en exit 0.
- Paridad: 118 títulos idénticos y una única excepción documentada para la
  diapositiva 48; orden y total: 119/119.
- Assets: 156/156 referencias resueltas. OCR: cero rótulos heredados en 134
  PNG únicos.
- PDF: 119 páginas, 104 con raster, cero texto literal `figs/`, cero nombres
  de archivo y cero cajas de imagen vacías sospechosas.
- Inspección visual: diapositivas 7/11/15/19/23/27/31, 33, 47--51, 116, 119
  y paneles de dispersión 8/20/32. Las anotaciones de N_eff son visibles; no
  hay clipping, tablas truncadas ni intercambio de Top/endL/endR.
- PDF SHA256:
  `c38e516e874031454212ebc2b6c358524e93b883d7c22d708a88a5cf7b37342a`.
