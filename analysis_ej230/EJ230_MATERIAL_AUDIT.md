# EJ-230 Material Audit

Branch: feat/ej230-endonly-mylar  
Date: 2026-06-14  
Audited by: analysis_ej230_endonly_mylar.py, validate_group_velocity_ej230.py  
Source files: `src/Materials.cc`, `src/DetectorConstruction.cc`

---

## 1. Scintillator: EJ-230 (OPSC-106)

| Property | Value | Source file:line | Notes |
|----------|-------|-----------------|-------|
| Material name | EJ-230 | Materials.cc:163 | `CreateEJ230()` |
| Geant4 handle | OPSC-106 | Materials.cc:163 | |
| Density | 1.023 g/cm³ | Materials.cc:24 | `G4Material` constructor |
| Refractive index | 1.58 (flat) | Materials.cc:41 | Flat RINDEX for all photon energies |
| Bulk absorption length | 120 cm | Materials.cc:181; DetectorConstruction.cc:183–186 | ABSLENGTH overridden in DetectorConstruction |
| Scintillation yield | 9 700 /MeV | Materials.cc:182 | `AddScintillatorProperties` |
| Emission peak | 391 nm | Materials.cc:181 | Peak wavelength in emission spectrum |
| Rise time | 0.5 ns | Materials.cc:182; DetectorConstruction.cc:194 | τ_rise overridden; 1.5 ns in EJ-230 datasheet, 0.5 ns in sim |
| Decay time | 1.5 ns | Materials.cc:182 | τ_decay |
| Rayleigh scattering | Active | Materials.cc:45–55 | λ_Rayl ≈ 150 cm at 391 nm, ∝ λ⁴ |
| GROUPVEL property | **Not set** | Materials.cc (entire file) | Geant4 computes v_g = c/n = 189.74 mm/ns |

### ABSLENGTH Override (DetectorConstruction.cc:183–186)
The simulation overrides ABSLENGTH for OPSC-106 to 120 cm. This value represents the bulk absorption length along the actual optical path (not the longitudinal bar distance). A photon traversing 70 cm longitudinally may travel significantly more actual path due to reflections.

### Rise Time Override (DetectorConstruction.cc:194)
The EJ-230 datasheet reports τ_rise ≈ 1.5 ns. The simulation uses τ_rise = 0.5 ns (overridden). This affects the rising edge of the arrival-time distribution and the SUM4 timing estimator.

### Derived Quantities
| Quantity | Value | Derivation |
|----------|-------|-----------|
| v_g = c/n | 189.74 mm/ns | c = 299.792 mm/ns, n = 1.58 |
| PDE at 391 nm | 61.65% | Interpolated from AFBR-S4N66P024M_pde.txt |
| Mylar loss/interaction | 10.0% | 1 − R = 1 − 0.90 |

---

## 2. Mylar Reflector (End-only + Mylar configuration)

| Property | Value | Source file:line |
|----------|-------|-----------------|
| Surface type | dielectric\_metal | Materials.cc:339 |
| Optical model | unified | Materials.cc:340 |
| Surface finish | ground | Materials.cc:341 |
| Reflectivity R | 0.90 (flat) | Materials.cc:347 |
| Specular lobe constant | 1.0 | Materials.cc:348 |
| σ\_α (lobe width) | 0.1° | Materials.cc:344 |
| Applied to faces | ±Y, ±Z (all non-end faces) | DetectorConstruction.cc:236 |

Note: The EJ-230 production run\_metadata.txt does **not** include Mylar optical constants. Values above are taken directly from Materials.cc as the authoritative Geant4 source.

Interpretation: R = 0.90 means 10% loss per Mylar interaction. After N bounces, surviving fraction ≈ 0.90^N. This strongly suppresses long-distance light collection: at d ≈ 70 cm, photons may undergo O(10–100) reflections, resulting in substantial attenuation beyond the bulk ABSLENGTH effect.

---

## 3. BarSkinReflector (EndTop reference configuration)

| Property | Value | Source file:line |
|----------|-------|-----------------|
| Reflectivity R | 0.98 (flat) | Materials.cc:324–325 |
| Applied to | All non-SiPM faces in EndTop geometry | DetectorConstruction |

The EndTop configuration uses BarSkinReflector (R = 0.98) on side faces, giving 2% loss per interaction. The 2.5× higher central-region NPE in EndTop vs End-only+Mylar is consistent with this difference in reflectivity, but is also a configuration-level result (different geometries, different wrapping materials).

---

## 4. SiPM Properties

| Property | Value | Source |
|----------|-------|--------|
| Physical device | AFBR-S4N66P014M | Simulation configuration |
| PDE curve source | AFBR-S4N66P024M | `data/sipm/AFBR-S4N66P024M_pde.txt` |
| PDE at 391 nm | 61.65% | Interpolated from PDE file |
| SiPM IDs (end-left) | 0–7 | ROOT file convention |
| SiPM IDs (end-right) | 8–15 | ROOT file convention |
| SiPM IDs (top, excluded) | 16+ | face\_type = 2, excluded from end-only analysis |

---

## 5. Optical Model Summary

```
EJ-230 (OPSC-106) bar:
  ┌─────────────────────────────────────────────────────┐
  │  Mylar (R=0.90, unified, ground, σα=0.1°)          │
  │  ┌───────────────────────────────────────────────┐  │
  │  │  EJ-230 scintillator                          │  │
  │  │  n=1.58, λ_bulk=120cm, yield=9700/MeV        │  │
  │  │  τ_rise=0.5ns, τ_decay=1.5ns                 │  │
  │  │  Rayleigh: active (λ_Rayl≈150cm at 391nm)    │  │
  │  └───────────────────────────────────────────────┘  │
  │  Mylar (R=0.90, unified, ground, σα=0.1°)          │
  └─────────────────────────────────────────────────────┘
  End-left SiPMs (IDs 0-7)      End-right SiPMs (IDs 8-15)
  PDE = 61.65% at 391nm         PDE = 61.65% at 391nm
```

---

## 6. Bootstrap Metadata

| Parameter | Value |
|-----------|-------|
| Type | Position-level pairs bootstrap |
| Engine | C++/OpenMP (bootstrap\_attenuation\_openmp.cpp) |
| N replicas | 200 |
| Seed (left) | 230123 |
| Seed (right) | 230456 |
| M2 failure rate (left) | 0.0% |
| M2 failure rate (right) | 0.0% |
| CI convention | 68% interval: [P16, P84] |

Bootstrap CIs are distinct from covariance diagonal errors reported by the Levenberg-Marquardt fitter. Use bootstrap CIs for the Summary and model comparison; use covariance errors for per-model parameter tables.
