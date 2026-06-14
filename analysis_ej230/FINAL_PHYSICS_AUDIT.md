# Final Physics Audit — EJ-230 End-only + Mylar Analysis

Branch: feat/ej230-endonly-mylar  
Date: 2026-06-14  
Auditor: EXEC analysis pipeline (automated + manual)

---

## 1. Material Constants — Verified Against Source Code

All constants used in the analysis are traced to `src/Materials.cc` or `src/DetectorConstruction.cc`:

| Constant | Value | Source | Status |
|----------|-------|--------|--------|
| Bulk ABSLENGTH | 120 cm | DetectorConstruction.cc:183–186 (override) | ✓ Verified |
| Scintillation yield | 9 700 /MeV | Materials.cc:182 | ✓ Verified |
| Peak emission | 391 nm | Materials.cc:181 | ✓ Verified |
| Refractive index | 1.58 (flat) | Materials.cc:41 | ✓ Verified |
| Density | 1.023 g/cm³ | Materials.cc:24 | ✓ Verified |
| Rise time | 0.5 ns | DetectorConstruction.cc:194 (override) | ✓ Verified |
| Decay time | 1.5 ns | Materials.cc:182 | ✓ Verified |
| Mylar R | 0.90 | Materials.cc:347 | ✓ Verified |
| Mylar specular lobe | 1.0 | Materials.cc:348 | ✓ Verified |
| Mylar σ_α | 0.1° | Materials.cc:344 | ✓ Verified |
| EndTop R (BarSkin) | 0.98 | Materials.cc:324–325 | ✓ Verified |
| PDE at 391 nm | 61.65% | AFBR-S4N66P024M_pde.txt interpolated | ✓ Verified |
| v_g = c/n | 189.74 mm/ns | Derived: 299.792/1.58 | ✓ Verified |
| GROUPVEL property | Not set | Materials.cc scan (validate_group_velocity_ej230.py) | ✓ PASS |

---

## 2. Group Velocity Validation

**PASS**: No GROUPVEL property found in Materials.cc. Geant4 computes v_g = c/n = 189.74 mm/ns from flat RINDEX = 1.58.

All apparent estimator velocities (FPT: 263.6 mm/ns, t50: 275.5 mm/ns, SUM4: 210.3 mm/ns) exceed v_g = 189.7 mm/ns. These are estimator-dependent apparent longitudinal responses, not violations of photon group velocity.

---

## 3. Data Pipeline Integrity

### End-only + Mylar dataset
- 31 positions confirmed
- SiPM IDs 0–7 (end-left), 8–15 (end-right) only
- face_type 2 (Top) excluded
- Pre-computed CSV `attenuation_curve.csv` and `sigma_t_sum4.csv` used directly

### EndTop dataset
- 31 positions confirmed
- face_type 0 (end-left) and face_type 1 (end-right) only
- face_type 2 (Top SiPMs, ~7.7M hits at x=0) explicitly excluded
- Same end-SiPM pipeline ensures apples-to-apples comparison

### SUM4 global_id masking
- EndTop ROOT files contain global_ids beyond 15 (Top channels)
- Masked with: `end_mask = (gid >= 0) & (gid < 16)` in timing extraction
- This ensures only end channels contribute to SUM4 timing

---

## 4. Attenuation Model Audit

### Model definitions
| Model | Formula | Free params |
|-------|---------|------------|
| M1 | N₀ exp(−d/λ) | 2 |
| M2 | A_s exp(−d/λ_s) + A_l exp(−d/λ_l) | 4 (constrained: λ_l = λ_s + exp(log_δ)) |
| M3 | A exp(−d/λ_p) + C | 3 |
| M4 | N₀ exp(−d/λ_tail) for d > 40 cm | 2 |

### Nominal fit results (left side)

| Model | Key parameter | χ²/ndf | AIC | Interpretation |
|-------|--------------|--------|-----|----------------|
| M1 | λ = 29.25 ± 1.61 cm | 2367 | 68649 | Catastrophic failure |
| M2 | λ_s = 8.95 cm, λ_l = 38.86 cm | 101 | 2734 | Empirical two-scale fit; χ²/ndf≫1 |
| M3 | λ_p = 19.98 cm, C = 13.28 PE | 1070 | 29968 | Fails; offset interpretation unclear |
| M4 | λ_tail = 38.13 cm (d>40cm) | 87.8 | n/a | Best single-exp in tail region; still ≫1 |

**Key finding**: M2 identifies two empirical attenuation scales but does not provide an exact description of the NPE profile (χ²/ndf = 101 ≫ 1). AIC comparison: M2 ≪ M3 ≪ M1. All models are empirical approximations.

M1 N₀ = 782 PE vs observed NPE at d≈0 = 2675 PE: factor ≈3.4 inconsistency. This confirms the two-component structure is physical.

### M2 reparametrization
M2 uses a constrained reparametrization to enforce λ_l > λ_s:
```
λ_l = λ_s + exp(log_δ), δ > 0
```
Covariance errors are propagated through this reparametrization.

---

## 5. Bootstrap Audit

| Parameter | Value |
|-----------|-------|
| Type | Position-level pairs bootstrap |
| Engine | C++/OpenMP (bootstrap_attenuation_openmp.cpp) |
| N replicas | 200 |
| Seed left | 230123 |
| Seed right | 230456 |
| M2 failure rate | 0.0% (both sides) |
| D_TAIL_CM | 40.0 cm (constexpr in C++ code) |
| Input units | cm (attenuation_cm.csv, d_left_cm = distance_left_mm / 10) |

**Unit check**: The original attenuation CSV has distance in mm. The C++ bootstrap uses D_TAIL_CM = 40.0 as a cm threshold. The `attenuation_cm.csv` intermediate file provides distances in cm to ensure correct M4 filtering (d > 40 cm ≡ 400 mm from near end, 21 points).

### Bootstrap confidence intervals (68% CI = P16–P84)

| Side | Parameter | Median | P16 | P84 |
|------|-----------|--------|-----|-----|
| Left | M1 λ (cm) | 29.32 | 27.34 | 31.14 |
| Left | M2 λ_s (cm) | 9.12 | 8.22 | 10.20 |
| Left | M2 λ_l (cm) | 39.21 | 37.90 | 40.84 |
| Left | M4 λ_tail (cm) | 38.16 | 37.22 | 39.21 |
| Right | M1 λ (cm) | 29.26 | 27.30 | 31.46 |
| Right | M2 λ_s (cm) | 8.93 | 8.18 | 10.14 |
| Right | M2 λ_l (cm) | 39.05 | 37.62 | 40.64 |
| Right | M4 λ_tail (cm) | 38.01 | 37.07 | 38.98 |

L/R consistency: M4 λ_tail values agree within bootstrap CI. M2 scales agree within bootstrap CI.

**Bootstrap CIs are distinct from covariance diagonal errors.** Covariance errors are reported in the attenuation fit tables; bootstrap CIs are used in the Summary and Appendix E.

---

## 6. Timing Audit

### σ_eq definition
σ_eq = σ(t_L − t_R) / √2, from SUM4 leading-edge (4 PE threshold).
Core Gaussian fit to the (t_L − t_R) distribution.

**This equals the per-end resolution only under equal, independent-end assumptions.**  
At bar ends (e.g., x = ±690 mm), NPE is strongly asymmetric (left ~2675 PE vs right ~11 PE), so the equal-end assumption breaks down and σ_eq is NOT a meaningful per-end resolution at those positions.

### σ_eq values

| x (mm) | σ_eq (ps) | σ_eq err (ps) | Note |
|---------|----------|--------------|------|
| −690 | 485 | 38 | Near-left: asymmetric NPE, σ_eq not per-end resolution |
| −400 | 249 | 17 | |
| 0 | 140 | 4 | Best timing; near-symmetric NPE |
| +400 | 244 | 13 | |
| +690 | 434 | 31 | Near-right: asymmetric NPE |

Intrinsic only. SPTR (~106 ps per detected photon) and FastIC (~10 ps) are not included.

### Apparent estimator velocities

| Estimator | v_app (mm/ns) | vs c/n = 189.74 mm/ns |
|-----------|--------------|----------------------|
| SUM4 Δt | 210.3 | +10.9% |
| FPT slope | 263.6 | +38.9% |
| t50 slope | 275.5 | +45.2% |

All exceed c/n. These are estimator-dependent apparent longitudinal responses, driven by position-dependent photon-population size, reflection path-length distributions, and threshold bias. **NOT evidence of superluminal propagation.**

---

## 7. Photon Budget Audit

Summed over all 31 positions:

| Fate | Count | Fraction |
|------|-------|---------|
| Surface+bulk absorbed | ~1.03 × 10⁹ | 89.8% |
| Escaped | ~5.32 × 10⁷ | 4.63% |
| Entered SiPM | ~6.36 × 10⁷ | 5.54% |
| Detected (×PDE) | ~3.92 × 10⁷ | 3.41% |

PDE efficiency: ε_PDE,eff = N_det / N_SiPM = 0.62 (consistent with PDE at 391 nm = 61.65%).

**Important caveat**: "Surface+bulk absorbed" is a single census category. The simulation cannot distinguish whether a photon was absorbed at the Mylar surface (after running out of reflectivity budget) or in the bulk (via ABSLENGTH exponential). A separate scoring mechanism would be needed to disentangle these.

Mylar reflector crossings: ~9.1 × 10⁹ total — confirming photons undergo many reflections before fate is determined.

---

## 8. EndTop Comparison Audit

| Observable | End-only+Mylar | EndTop | Ratio |
|-----------|----------------|--------|-------|
| Reflector R | 0.90 | 0.98 | — |
| NPE at x=0 (left) | 60.1 PE | 151.0 PE | 2.51× |
| Central-region NPE (|x|≤100mm) | 62.0 PE | 155.1 PE | 2.50× |
| M4 λ_tail (left) | 38.13 cm | 37.86 cm | ~equal |

**Editorial constraint (ADDENDUM FINAL item 6)**: The 2.5× NPE difference is *consistent with* higher-reflectivity boundaries in EndTop (R=0.98 vs R=0.90), but this is a **configuration-level result**. The two geometries differ in:
- Reflector material (Mylar vs BarSkinReflector)
- Reflector model (dielectric_metal/unified vs BarSkinReflector)
- SiPM layout (end-only vs end+top)
- Possible other geometric differences

The similar M4 tail slopes (38.1 vs 37.9 cm) do **not** prove reflector independence and may reflect model error (χ²/ndf ≫ 1 for both).

---

## 9. Addendum Final — All 14 Corrections Applied

| Item | Description | Status |
|------|-------------|--------|
| 1 | M2 language: "identifies empirical scales; χ²/ndf≫1, not exact" | ✓ Applied (Frames 7, 15) |
| 2 | Distinguish covariance vs bootstrap errors | ✓ Applied (Appendix C vs E) |
| 3 | Document bootstrap type (position-level pairs, C++/OpenMP, seed) | ✓ Applied (Appendix E, G) |
| 4 | FPT: "estimator-dependent apparent longitudinal response" | ✓ Applied (Frames 11, 12, F) |
| 5 | Two NPE metrics (center_point x=0 vs central_region |x|≤100mm) | ✓ Applied (Frame 4, 8) |
| 6 | EndTop: "consistent with higher-R boundaries; configuration-level" | ✓ Applied (Frames 8, 9, 15) |
| 7 | M3 C: "finite-range offset" not "recirculated floor" | ✓ Applied (Frames 5, 6) |
| 8 | Photon budget: surface+bulk not separated | ✓ Applied (Frame 14) |
| 9 | Motivation: include escaped category | ✓ Applied (Frame 2) |
| 10 | σ_eq everywhere (not σ_single in labels) | ✓ Applied (all sigma labels) |
| 11 | M2 full parameters (A_s, A_l, λ_s, λ_l) in model table | ✓ Applied (Frame 7) |
| 12 | Summary with approved template | ✓ Applied (Frame 15) |
| 13 | verify_deck.py extended forbidden-phrase checks | ✓ verify_deck_ej230.py created |
| 14 | Don't interrupt analysis | ✓ Completed without interruption |

---

## 10. Prohibited Actions Check

| Prohibition | Status |
|-------------|--------|
| No EJ-204 values in EJ-230 report | ✓ PASS (verify_deck checks this) |
| No simulation modification | ✓ PASS (source unchanged) |
| No git push | ✓ PASS (commit only) |
| No hardcoded values without source citation | ✓ PASS (all traced to Materials.cc) |

---

## 11. Verification Results

- `verify_deck_ej230.py`: **PASS** (all checks)
- `validate_group_velocity_ej230.py`: **PASS** (GROUPVEL absent)
- LaTeX compilation: **PASS** (22 pages, no fatal errors)
- Remaining overflows: 3 hbox (≤15.8pt, two are Beamer column artifacts), 6 vbox (≤26pt, visual clipping only)
