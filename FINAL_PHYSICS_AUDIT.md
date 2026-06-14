# FINAL PHYSICS AUDIT — End-only + Mylar EJ-204
**Branch:** `feat/endonly-mylar`  
**Date:** 2026-06-14  
**Auditor:** EXEC FINAL review  
**Deck:** `presentation/endonly_mylar/`

---

## 1. Simulation material properties

### 1.1 Source tracing — `src/Materials.cc`

| Property | Value | Source line |
|----------|-------|-------------|
| EJ-204 density | 1.023 g/cm³ | `CreateEJ204()` → `FindOrBuildPVT("EJ-204")` (line 84); PVT base at line 24 |
| RINDEX (flat) | 1.58 (4 energy nodes: 2.0–4.0 eV) | line 41: `rIdx[nOpt] = {1.58, 1.58, 1.58, 1.58}` |
| ABSLENGTH (flat) | 160 cm (passed as `attenuation`) | line 99: `160.0 * cm`; stored via `absL[i] = attenuation` (line 42) |
| RAYLEIGH | λ⁴ scaling from `rayleighRef = 1.5 m` at `peakWavelength` | lines 44–50 |
| Scintillation yield | 10 400 photons/MeV | line 99: `10400.0 / MeV` |
| Rise time τ₁ | 0.7 ns | line 100: `0.7 * ns` |
| Decay time τ₂ | 1.8 ns | line 100: `1.8 * ns` |
| Peak emission λ | 408 nm | line 99: `408.0 * nm` (used as `peakWavelength`) |
| Emission spectrum | 21-point digitization of Eljen datasheet | lines 88–97 |

No `GROUPVEL` property is set. Geant4 11.x computes group velocity as
`v_g = c/n` from the flat `RINDEX`, giving **v_g = 189.742 mm/ns**.

### 1.2 Rayleigh scattering lengths (derived)

The code sets `rayleighRef = 1.5 m` at `peakWavelength = 408 nm` and scales
as `L_R(λ) = rayleighRef × (λ/408 nm)⁴`.

| Energy node | λ (nm) | L_R (cm) |
|-------------|--------|----------|
| 2.0 eV | 619.9 | 799.5 |
| 2.6 eV | 476.9 | 279.9 |
| 3.1 eV | 399.9 | 138.5 |
| 4.0 eV | 310.0 | 50.0 |
| **408 nm (peak)** | 408.0 | **150.0** |

Rayleigh scattering is an additional loss mechanism beyond ABSLENGTH.
It is **not** documented in the EJ-204 datasheet; it is a model choice
to approximate wavelength-dependent scattering in the bulk.

### 1.3 Mylar wrapping — `src/Materials.cc` (CreateMylarReflector)

| Property | Value | Note |
|----------|-------|------|
| Surface model | `dielectric_metal` | specular |
| REFLECTIVITY | 0.90 (uniform over energy range) | ⟹ loss per bounce **1−R = 0.10** |

**Previous error:** The deck stated "1−R = 0.90" (= 90% loss per bounce).
**Correct value:** 1−R = 0.10 (10% loss per bounce). Cumulative loss over
many reflections is what dominates far-end attenuation, not any single interaction.

---

## 2. Attenuation analysis — corrected framework

### 2.1 ABSLENGTH vs fitted λ_eff

`ABSLENGTH = 160 cm` is a **material bulk property** applied to the actual
3-D optical path length `s` of each photon (zig-zag due to reflections).
It is **not** a property of the longitudinal distance `d` from the SiPM end.

The quantity fitted from NPE(d) is `λ_eff` (or `λ_tail`), an **effective
system longitudinal scale** that combines:

- Bulk absorption over actual (zig-zag) path: dominated by ABSLENGTH=160 cm
- Rayleigh scattering (L_R ≈ 150 cm at 408 nm)
- Geometry factor: for a bar of cross-section w×h, the mean path per
  longitudinal step depends on w, h, and the angular acceptance cone
- Repeated Mylar losses (1−R = 0.10 per bounce)
- SiPM active area fraction and Geant4 PDE model
- Trigger / analysis acceptance

Fitted λ_eff < λ_bulk is **expected and correct**; it is not inconsistent
with ABSLENGTH=160 cm.

### 2.2 Four attenuation models

All fits use data from 31 positions over ±690 mm (d ≡ distance from
nearest bar end, range 10–690 mm). NPE uncertainties are SEM from 2 000
events per position. χ²/ndf uses raw χ² as −2 log L.

#### M1 — Single exponential (reference)

```
N(d) = N₀ · exp(−d/λ)
```

| Side | λ (cm) | χ²/ndf | AIC | BIC |
|------|--------|--------|-----|-----|
| Left | 32.16 ± 1.85 | 2320.9 | 67 311 | 67 313 |
| Right | 31.41 ± 1.82 | 2475.3 | — | — |

Poor χ²/ndf indicates a single exponential is inadequate. The near-end
NPE falls faster than the far-end tail, as expected from superposition
of prompt and bulk-absorption components.

#### M3 — Exponential + constant floor (exp+offset)

```
N(d) = A · exp(−d/λ_prompt) + C
```

| Side | λ_prompt (cm) | C (PE) | χ²/ndf | AIC | BIC |
|------|---------------|--------|--------|-----|-----|
| Left | 20.55 ± 1.08 | 19.22 ± 2.00 | 954.1 | 26 722 | 26 726 |

Still poor χ²/ndf; the constant floor is unphysical (NPE continues to
fall beyond d = 500 mm). Bootstrap: λ_prompt = 20.71 [19.11, 22.55] cm.

#### M2 — Constrained double exponential (preferred)

```
N(d) = A_s · exp(−d/λ_s) + A_l · exp(−d/λ_l)
```
Ordering constraint: λ_l = λ_s + exp(log δ) → λ_l > λ_s guaranteed.

| Side | λ_s (cm) | λ_l (cm) | A_s | A_l | χ²/ndf | AIC | BIC |
|------|----------|----------|-----|-----|--------|-----|-----|
| Left | 9.46 ± 0.44 | 42.98 ± 0.92 | 2284 ± 101 | 396 ± 25 | 95.3 | 2 582 | 2 588 |

**ΔAICc (M2 vs M1) ≈ −64 729** → M2 strongly preferred.  
Bootstrap (200 replicas, pairs, seed=42):
- λ_s: median = 9.49, 68% CI = [8.75, 10.94] cm  
- λ_l: median = 43.30, 68% CI = [41.72, 45.43] cm  
- **Failure rate: 0% (stable)**

#### M4 — Tail slope only (d > 40 cm)

```
N(d) = N₀ · exp(−d/λ_tail)    [d > D_min = 40 cm]
```

| Side | λ_tail (cm) | χ²/ndf | n_points |
|------|-------------|--------|----------|
| Left | 41.69 ± 0.90 | 96.5 | 21 |
| Right | 41.95 ± 0.92 | — | 21 |

Bootstrap (200 replicas): λ_tail = 41.83 [40.51, 43.32] cm.

### 2.3 λ_tail range sensitivity

| d_min | λ_tail L (cm) | λ_tail R (cm) | χ²/ndf L |
|-------|--------------|--------------|----------|
| 30 cm | 39.80 ± 1.00 | 39.80 ± 1.04 | 185.1 |
| **40 cm (ref)** | **41.69 ± 0.90** | **41.95 ± 0.92** | **96.5** |
| 50 cm | 43.46 ± 0.82 | 43.78 ± 0.80 | 50.8 |

Variation across the three d_min choices: ~3.7 cm (≈ 9%). The trend
(larger d_min → larger λ_tail) reflects the continued curvature of the
attenuation curve even at large d; χ²/ndf improves as near-end points
are excluded, confirming the double-component structure persists.

---

## 3. Near- and far-end NPE

**Previous error in deck:** At x = −690 mm (source at left end), the deck
showed both ends at ~2920 PE and ~2914 PE — swapping near and far.

**Corrected values** from `deck_values.json`:

| Source position | Near-end NPE | Far-end NPE | Ratio near/far |
|----------------|-------------|------------|----------------|
| x = −690 mm | Left (near) = 2926 PE | Right (far) = 16.83 PE | 174 |
| x = 0 mm (center) | Left = 77.2 PE | Right = 77.0 PE | 1.00 |
| x = +690 mm | Right (near) = 2914 PE | Left (far) = 17.0 PE | 171 |

The far/near asymmetry at the bar ends (~×170) is consistent with the
fitted λ_tail ≈ 42 cm and bar half-length of 700 mm (≈ 16.7 × λ_tail).

---

## 4. Photon budget

Source: Geant4 boundary process census, all 31 positions × 2 000 events.

**Panel A — Terminal fates (sums to 100%):**

| Fate | Count | Fraction |
|------|-------|---------|
| Surface + bulk absorbed | 1.10 × 10⁹ | 89.4% |
| Escaped through bar faces | 5.95 × 10⁷ | 4.84% |
| Entered SiPM active area | 7.10 × 10⁷ | 5.77% |
| **Total generated** | **1.231 × 10⁹** | **100%** |

Surface absorption and bulk absorption are **not separately resolved**
in the Geant4 boundary census. The 89.4% figure combines both.

**Panel B — Conditional SiPM conversion:**

| Step | Fraction of SiPM-entering |
|------|--------------------------|
| Detected (PDE applied) | 3.63% of generated = **63%** of SiPM-entering |
| PDE model in simulation | 63% (flat, from AFBR-S4N66P024M datasheet at 420 nm) |

**Previous errors corrected:**
1. `1−R = 0.10` (loss per Mylar bounce), not 0.90.
2. Panel A now closes to exactly 100% (terminal fates, not conversion chain).

---

## 5. Group velocity and velocity estimators

### 5.1 Group velocity

No `GROUPVEL` property is set in `src/Materials.cc`. Geant4 11.x
computes photon propagation speed from `d(ω)/dk` applied to the
flat `RINDEX = 1.58` array, which gives:

```
v_g = c / n = 299.7925 / 1.58 = 189.742 mm/ns
```

This is the **only velocity scale** from the simulation physics.
It does **not** include dispersive corrections (which would require
a non-flat RINDEX or explicit GROUPVEL property).

### 5.2 Apparent estimator velocities

| Estimator | v_app (mm/ns) | Excess over c/n |
|-----------|--------------|-----------------|
| SUM4 (sum of 4 channels) | 213.72 | +24.0 mm/ns (+12.6%) |
| FPT (first photon time) | 264.09 | +74.3 mm/ns (+39.2%) |
| t50 (50th-percentile arrival) | 275.69 | +85.9 mm/ns (+45.3%) |
| **c/n (physics)** | **189.742** | reference |

All three estimator velocities exceed c/n. This is **not a violation of
special relativity** — it is a statistical bias:

- **SUM4**: the trigger fires when the sum of the 4 channel amplitudes
  crosses threshold; near-end channels fire earlier, biasing the centroid.
- **FPT**: the minimum-order statistic of N photon arrival times. With
  more photons at the near end (higher NPE), the statistical minimum
  is pulled earlier, steepening the slope above c/n.
- **t50**: same mechanism as FPT, more pronounced.

These are **apparent estimator velocities**, not propagation speeds.
The quantity to compare with c/n is the group velocity of photons,
not the slope of the trigger time vs position.

### 5.3 FPT bias direction (corrected)

**Previous error:** "fewer photons at far end → later first arrival → steeper slope."  
**Correct statement:** More photons at the near end shift the statistical
minimum of arrival times **earlier**, making the FPT slope steeper than
what c/n alone would predict. The mechanism is order-statistics bias,
not photon deficit.

---

## 6. Timing resolution

### 6.1 σ_eq definition

```
σ_eq(x) = σ(t_L − t_R) / √2
```

This equals the per-end single-channel resolution **only** under:
1. Equal NPE at both ends (σ_L = σ_R)
2. Statistically independent end channels

Condition (1) breaks down strongly near the bar ends (NPE ratio ~170).
The quantity σ_eq is therefore position-dependent and **not** a
position-independent resolution figure.

**Previous label:** σ_LR/√2 (opaque to the assumption)  
**Corrected label:** σ_eq, with explicit caveat

### 6.2 Numerical values

| Position x (mm) | σ_eq (ps) |
|-----------------|-----------|
| −690 | 406.2 ± 26.6 |
| −400 | 210.6 ± 11.2 |
| 0 | 141.5 ± 4.3 |
| +400 | 209.3 ± 8.9 |
| +690 | 372.1 ± 19.1 |

These are **intrinsic simulation values** (optical propagation + SUM4
estimator only). Not included:
- SPTR ≈ 106 ps (SiPM single-photon time resolution, Hamamatsu data)
- FastIC electronics jitter ≈ 10 ps
- ADC sampling jitter

The total resolution σ_total > σ_eq. For the center position:
σ_total ≈ √(141.5² + 106²/NPE_eff) ps (schematic; actual convolution
depends on NPE distribution and SPTR model at the photoelectron level).

---

## 7. EndTop comparison

| Quantity | End-only + Mylar | EndTop (BarSkinReflector R=0.98) |
|----------|-----------------|----------------------------------|
| λ_tail (cm) | 41.69 ± 0.90 | ~41.4 (from exec07 run) |
| NPE at x=0 | 77.2 PE | 194.4 PE |
| NPE ratio (EndTop/Mylar) | — | ~2.5× |
| Readout | End SiPMs only | End SiPMs only |

**Previous claim:** "Same attenuation for Mylar and BarSkinReflector."  
**Corrected claim:** "Numerically similar λ_tail values do not establish
reflector independence of the attenuation length. The shorter-scale
component (λ_s ≈ 9.5 cm) and overall NPE yield differ significantly
between configurations. A controlled comparison varying only the
reflectivity while keeping geometry, source, and analysis identical
is needed."

---

## 8. Problems found and corrections applied

| # | Problem | Location | Correction |
|---|---------|----------|------------|
| 1 | ABSLENGTH described as longitudinal scale | Frame 5 (old) | Clarified as optical path property; added λ⁴ Rayleigh note |
| 2 | Only M1 (single-exp) fit reported | analysis_endonly_mylar.py | Added M2 (constrained dbl-exp), M3 (exp+offset), M4 (tail) |
| 3 | No model comparison criteria | analysis_endonly_mylar.py | Added AIC/AICc/BIC; M2 preferred by ΔAIC ≈ −64 000 |
| 4 | No uncertainty on M2/M3/M4 parameters | analysis_endonly_mylar.py | Added 200-replica pairs bootstrap; 0% failure for M2 |
| 5 | No λ_tail range sensitivity | analysis_endonly_mylar.py | Added d≥30/40/50 cm; spread = 3.7 cm |
| 6 | Near/far NPE swapped in Frame 4 | build_deck.py | Fixed: near=2926 PE, far=16.83 PE at x=−690 mm |
| 7 | Photon budget Panel A not closing to 100% | build_deck.py + analysis | Restructured as Panel A (fates) + Panel B (conversion) |
| 8 | 1−R stated as 0.90 instead of 0.10 | Frame 8 (old) | Corrected to 1−R = 0.10; caveat on cumulative loss |
| 9 | FPT bias direction wrong | Frame 11 (old) | Corrected: more photons → earlier minimum → steeper slope |
| 10 | v_eff interpreted as propagation velocity | Frame 12 (old) | Renamed to v_app; all v_app > c/n explained as estimator bias |
| 11 | σ_LR/√2 presented as per-end resolution | Frame 13 (old) | Renamed σ_eq; qualified equal-end assumption |
| 12 | EndTop claimed "same attenuation" | Frame 9 (old) | Conclusion softened; differences noted |
| 13 | Summary claimed ≤50 ps performance | Frame 15 (old) | Removed; conservative statement on intrinsic resolution |
| 14 | No GROUPVEL discussion | — | Added Appendix F and v_group key in deck_values.json |
| 15 | \pm in range-sensitivity table outside math mode | build_deck.py | Wrapped in $…$ |

---

## 9. Pending validations

| Validation | How to perform | Status |
|-----------|----------------|--------|
| Optical path s vs longitudinal d ratio | Run single-photon simulation with photon tracking; compare mean s to d | Not done (requires new sim run) |
| GROUPVEL independence | Run monoenergetic photon; verify t_arrival/d = 1/v_g | Validation script written (`validate_group_velocity.py`) |
| Rayleigh contribution to λ_eff | Rerun simulation with RAYLEIGH disabled; compare λ_tail | Not done (requires new sim run) |
| 1-R = 0.10 from boundary census | Count photon-Mylar interactions vs absorbed-at-Mylar | Not done (requires custom stepping action) |
| SPTR convolution | Fold 106 ps Gaussian per photoelectron; recompute σ_eq | Documented in Appendix G; not yet implemented |
| Test-beam comparison | April-2026 TB used top readout + Tyvek (different geometry) | No direct comparison possible |

---

## 10. File provenance

| File | Role | Last modified | Commit |
|------|------|--------------|--------|
| `src/Materials.cc` | Simulation material definitions | not modified in this audit | — |
| `presentation/endonly_mylar/analysis_endonly_mylar.py` | Physics analysis + figure/table generation | 2026-06-14 | `6e128af` |
| `presentation/endonly_mylar/build_deck.py` | Beamer deck generator | 2026-06-14 | `1ef8ddf` |
| `presentation/endonly_mylar/deck_values.json` | Single numerical source of truth for deck | 2026-06-14 | `20653d8` |
| `presentation/endonly_mylar/main.tex` | Auto-generated LaTeX (do not edit by hand) | 2026-06-14 | `20653d8` |
| `presentation/endonly_mylar/verify_deck.py` | Traceability checker | 2026-06-14 | `9a97f16` |
| `audit_backup/` | Pre-audit snapshot + patch | 2026-06-14 | `bf41a16` |

All `verify_deck.py` checks pass (0 failures, 40 checks across 5 suites).
