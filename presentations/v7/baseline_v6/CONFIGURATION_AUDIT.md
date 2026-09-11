# CONFIGURATION_AUDIT.md — SHiP Timing Bar Simulation Generations

Generated 2026-08-21. Every configuration that produced numbers shown in any presentation.

---

## GEN-1: Historical END+TIR (before EXEC_21 fix)

- **Branch / build**: ej200_endonly or early exec builds (before exec21-optfix)
- **Reflector**: `dielectric_metal` skin surface on BarLV → **eliminates TIR**
- **Effect**: ~300× photon deficit at END faces (0.37 PE/event reported)
- **EJ-230 att**: 1200 mm
- **Estimator**: first-photon or per-SiPM mean
- **Status**: **INVALID for physics** — surface bug prevents normal TIR
- **Historical note**: this is the "1-PE anomaly" (see CONTENT_AUDIT.md B6)

---

## GEN-3 END-only: scan_end_vikuiti (build_end_vikuiti Phase 7)

> **EXEC_23 FIX-08 (2026-09-02):** Previously mislabeled "GEN-2". Evidence: phase_ab.py reads
> scan_end_vikuiti; presentations/v6/README.md states scan_end_vikuiti = build_end_vikuiti Phase 7
> (air gap + Mylar + border surface, identical geometry to GEN-3 current, but N_TOP=0).

- **Build**: `build_end_vikuiti` Phase 7 (same codebase as GEN-3 current)
- **Data path**: `results/scan_end_vikuiti/`
- **Geometry**: 1400×60×10 mm³ bar; END-only (N_TOP=0, no TOP SiPMs)
- **Air gap**: 0.10 mm (kAirGapThickness) — explicit physical volume
- **Mylar reflector**: 0.05 mm — explicit physical volume
- **Reflector surface**: `CreateBarSkinReflector()` used as BORDER surface (air→Mylar boundary)
  - Type: `dielectric_dielectric`; Model: `unified`; Finish: `polished`
  - REFLECTIVITY = 0.95 (constant)
- **Scintillator**: EJ-230 (att=1200 mm, τ=1.5 ns, Y=9700/MeV, n=1.58)
- **Also tested**: EJ-200, EJ-204
- **SiPM**: 8 per END face = 16 total; wavelength-dependent PDE (AFBR-S4N66P024M)
- **Events**: 5000 per position, 13 positions (x = -600 to +600 mm step 100)
- **Estimator**: m-average END estimator, optimal m* (phase_ab_optimal.csv)
- **Key result**: σ_END = 50.49 ps (m*=7, EJ-230, x=0); bar mean 51.95 ps
- **Reference files**: `analysis/optim/phase_ab_optimal.csv`, figM1/M2/M3

---

## GEN-3: EndSparseTop V2 — CURRENT (scan_end_vik_sparse_top_v2)

- **Build**: `build_end_vikuiti` (same codebase, different macro config)
- **Branch**: `feat/bar-end-vikuiti` (current)
- **Data path**: `results/scan_end_vik_sparse_top_v2/`
- **Geometry**: 1400×60×10 mm³ bar
- **Air gap**: 0.10 mm (kAirGapThickness) — explicit physical volume
- **Mylar reflector**: 0.05 mm (kReflectorThickness) — explicit physical volume
- **Scintillator surface (bar→air)**:
  - `CreateBarSurface()`: `dielectric_dielectric`, `unified`, `polished`, σ=0
  - Implements natural Fresnel + TIR for θ > arcsin(1/1.58) = 39.3°
- **Reflector surface (air→Mylar/Vikuiti)**:
  - `CreateBarSkinReflector()` used as BORDER surface on AirGap→Mylar boundary
  - Type: `dielectric_dielectric`; Model: `unified`; Finish: `polished`
  - REFLECTIVITY = 0.95 (constant, 2-tier: TIR automatic + R=0.95 for non-TIR)
  - ⚠ Function named "SkinReflector" but used as BORDER surface — naming is legacy
- **SiPM surface**: `dielectric_metal`, wavelength-dependent PDE (AFBR-S4N66P024M)
- **END SiPMs**: 8 per face = 16 total; 6×6 mm² active area each
- **TOP SiPMs**: N_TOP sparse, uniformly spaced along bar length
  - N_TOP ∈ {4, 8, 14, 20}; positions `SparseTopSiPMCenterX(idx, nTotal)`
  - For N=20: x_SiPM ∈ [-664.3, -594.3, …, +594.3, +664.3] mm (step ≈70 mm)
- **Scintillators simulated**: EJ-230 (primary); EJ-204, EJ-200 (comparison)
- **Events**: 5000 per position per configuration
- **Positions**: x ∈ {-600, -500, …, +600} mm (13 positions)
- **Files**: `photon_hits_ej230_ntop{4,8,14,20}_x{-600,...,600}.root`
- **Tree**: `sipm_hits`; branches: event_id, face_type (0=END-L, 1=END-R, 2=TOP), local_id, time_ns, pde, x_mm, gun_x_mm

### Key results (GEN-3):

| Quantity | Value | Notes |
|----------|-------|-------|
| v_eff (END, m=8) | 177.6 mm/ns | from Δt vs x linear fit |
| v_eff (END, first photon) | 184.6 mm/ns | from napkin/legacy analysis |
| v_g (group) | 189.7 mm/ns | c/n = c/1.58 |
| σ_END (m*=7/8, x=0) | 53.68 ps | from m-avg in hybrid geometry |
| σ_TOP (k*=3, x=0) | 15.20 ps | from pooled k-th order stat |
| σ_BLUE (x=0) | 15.21 ps | BLUE with IQR-based weights |
| w_END | 0.013 | IQR variance + empirical cov |
| w_TOP | 0.987 | |
| ρ(END,TOP) | 0.16 | empirical correlation |
| σ_END (END-only geometry) | 50.5 ps | m*=7; no TOP SiPMs |
| Global σ_BLUE (k=3,m=8) | mean 19.5 ps | across x ∈ [-600,+600] mm |
| Oracle σ_BLUE | mean 16.1 ps | local optimal k*(x), m*(x) |

---

## IMPORTANT: Cross-configuration differences

### v_eff: 184.6 vs 177.6 mm/ns

- **184.6 mm/ns**: from napkin/legacy analysis using first-photon or near-first estimator
- **177.6 mm/ns**: from v5 analysis using m=8 average estimator
- **Both are correct** for their respective estimators. The difference arises because the m=8 average uses slightly later photons than the first-photon statistic.
- **σ_x from Δt**:
  - Using v_eff=184.6: σ_x ≈ 10.3 mm
  - Using v_eff=177.6: σ_x ≈ 7.9 mm (with σ_Δt from m=8)

### BLUE weights: 0.013 vs 0.086

- **0.013** (best_est_analysis.py): IQR-based σ² for σ_E² and σ_T²; C_ET from np.cov. In BLUE formula w_E = (σ_T²−C_ET)/(σ_E²+σ_T²−2C_ET): σ_T² in numerator, σ_E² only in denominator.
- **0.086** (v5 analysis_v5.py): np.var() σ² (includes tails) for σ_E² and σ_T²; C_ET from np.cov. Same formula role distribution as above.
- **Recommendation**: use 0.013 as the canonical result (consistent with cross-check σ_BLUE)
- The IQR-based method is physically more meaningful since σ_IQR is the reported resolution metric

### Reflector: R=0.95 vs R=0.98

- **R=0.95**: actual code value in `CreateBarSkinReflector()` (constant)
- **R=0.98**: napkin estimate from manufacturer Vikuiti spec, used for analytical predictions only
- The code says "R=0.98 Vikuiti ESR" in comments but implements R=0.95

### END σ: 50.49 vs 53.68 ps

- **50.49 ps**: GEN-3 END-only (scan_end_vikuiti), m*=7, x=0 (from phase_ab_optimal.csv)
- **53.68 ps**: GEN-3 hybrid with N_TOP=20 (scan_end_vik_sparse_top_v2), m=8, x=0
- _Hypothesis_: TOP SiPMs intercept early photons that would otherwise reach END, accounting for
  the ~3 ps difference. Not directly demonstrated — the comparison conflates N_TOP and estimator
  parameters simultaneously. See EXEC_23 FIX-08 note in talk_v6.tex B1.
