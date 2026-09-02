# FINAL_NUMBERS.md — Authoritative results for talk_v6

Generated 2026-08-21.
Source: GEN-3 dataset `results/scan_end_vik_sparse_top_v2/`
Canonical analysis: `analysis/optim/root_best_est/best_est_analysis.py`
Config: EJ-230, Vikuiti ESR (R=0.95 constant), N_TOP=20 unless noted.

---

## First-principles napkin (from napkin_macros.tex)

| Quantity | Value | Source |
|---------|-------|--------|
| n (EJ-230) | 1.58 | Materials.cc |
| θ_c | 39.27° | arcsin(1/1.58) |
| f_TIR (slab) | 0.5484 | 2cos(θ_c)−1 |
| f_TIR/end | 0.2742 | f_TIR/2 |
| v_g = c/n | 189.7 mm/ns | c/1.58 |
| N_bounce (H=10mm face, L/2) | 85.63 | geometric |
| N_bounce (W=60mm face, L/2) | 14.27 | geometric |
| Λ_refl (H-face, R=0.98) | 404.6 mm | napkin, R=0.98 |
| Λ_refl (H-face, R=0.95) | 159.4 mm | from \NapLambdaReflRlow (EXEC_23 FIX-05) |
| N_pe/end (TIR-only napkin) | 570 | EJ-230, x=0 |
| N_pe/end (Geant4) | 701.3 | GEN-3, x=0 |
| Surplus (reflector-recovered) | +23% | (701.3−570)/570 |

---

## Material properties (from Materials.cc)

| Property | EJ-200 | EJ-204 | EJ-230 |
|---------|--------|--------|--------|
| λ_att [mm] | 3800 | 1600 | **1200** |
| τ_rise [ns] | 0.9 | 0.7 | **0.5** |
| τ_decay [ns] | 2.1 | 1.8 | **1.5** |
| Yield [ph/MeV] | 10000 | 10400 | 9700 |
| n | 1.58 | 1.58 | 1.58 |
| N_pe/end (G4, x=0) | 1311 | 941 | **701.3** |
| σ_END END-only, m*-opt, x=0 [ps] | 54.86 | 53.36 | 50.49 |

---

## Effective velocity (from analysis_v5.py, GEN-3, m=8, EJ-230, N_TOP=20)

| Quantity | Value |
|---------|-------|
| v_eff (m=8 average) | **177.6 mm/ns** |
| v_eff (first-photon, legacy) | 184.6 mm/ns |
| v_g | 189.7 mm/ns |
| fit offset a | 0.44 ± 0.38 ps |
| fit slope b | −11.2587 ± 0.0011 ps/mm |
| χ²/ndf | 1250/11 (p ≈ 0) |
| Residual RMS | 14.5 ps |
| σ_x (END Δt, m=8) | **7.9 mm** |

---

## Point performance at x=0 (GEN-3, EJ-230, N_TOP=20, k=3, m=8)

| Estimator | σ_IQR [ps] | C_ET [ps²] | ρ_BLUE = C/(σ_Eσ_T) |
|-----------|-----------|-----------|-------------------|
| END (m=8) | **53.68** | 196.3 | 0.241 |
| TOP (k=3) | **15.20** | — | — |
| **BLUE** | **15.21** | — | — |

ρ = 0.160 (corrcoef with empirical σ's, NOT BLUE-consistent; see G1) — do not use inside BLUE formula.

**EXEC_23 FIX-03 (2026-09-02): w_END not resolved at n=5000 (CP-B w-scan).**

Three covariance conventions (same 5000 events, x=0):

| Convention | C_ET [ps²] | w_END | σ_IQR^event [ps] |
|------------|-----------|-------|-----------------|
| best_est: IQR σ² + np.cov | 196.3 | 0.013 | 15.21 |
| robust IQR (polarization) | 80.5 | 0.051 | 15.12 |
| empirical: np.var + np.cov | 196.3 | 0.086 | 15.03 |

All three within 0.18 ps ≈ bootstrap SE (≈0.17 ps at n=5000). σ_IQR(w) flat to within one SE over w ∈ [0, 0.10] (total variation 0.26 ps). The choice of covariance convention is not resolved by this dataset; no design conclusion depends on it.

Sidecar: analysis/optim/root_best_est/blue_wscan_x0.{csv,meta.json,root}

---

## Global performance (k=3, m=8, N_TOP=20, EJ-230)

| x [mm] | σ_END [ps] | σ_TOP [ps] | σ_BLUE [ps] | w_END | ρ |
|--------|-----------|-----------|------------|-------|---|
| −600 | 59.3 | 16.3 | 16.7 | 0.115 | 0.088 |
| −500 | 61.0 | 23.2 | 24.6 | 0.208 | 0.095 |
| −400 | 56.7 | 17.0 | 18.1 | 0.178 | 0.098 |
| −300 | 56.7 | 17.3 | 18.6 | 0.183 | 0.090 |
| −200 | 53.6 | 23.2 | 24.5 | 0.237 | 0.086 |
| −100 | 53.0 | 16.4 | 16.6 | 0.166 | 0.070 |
|    0 | 53.7 | 15.2 | **15.0** | 0.086 | 0.160 |
| +100 | 53.8 | 16.5 | 17.2 | 0.166 | 0.072 |
| +200 | 55.0 | 24.3 | 24.7 | 0.229 | 0.135 |
| +300 | 55.3 | 17.2 | 18.1 | 0.185 | 0.087 |
| +400 | 57.6 | 17.1 | 18.0 | 0.165 | 0.083 |
| +500 | 59.7 | 23.5 | 24.5 | 0.214 | 0.089 |
| +600 | 61.0 | 16.8 | 16.9 | 0.118 | 0.082 |

Global mean σ_BLUE = 19.5 ps; min = 15.0 ps; max = 24.7 ps.

---

## Oracle performance (k*(x), m*(x), N_TOP=20, EJ-230)

| x [mm] | k* | m* | σ_END | σ_TOP | σ_BLUE |
|--------|----|----|-------|-------|--------|
| −600 | 2 | 10 | 59.0 | 15.4 | 15.4 |
| −500 | 1 | 11 | 60.1 | 17.8 | 17.6 |
| −400 | 2 |  5 | 55.4 | 15.9 | 15.7 |
| −300 | 2 |  6 | 56.3 | 16.0 | 15.6 |
| −200 | 2 |  8 | 53.6 | 17.6 | 18.3 |
| −100 | 2 |  6 | 52.8 | 14.9 | 14.7 |
|    0 | 3 |  8 | 53.7 | 15.2 | **15.0** |
| +100 | 2 |  5 | 53.1 | 15.6 | 15.2 |
| +200 | 2 | 10 | 54.5 | 18.2 | 18.5 |
| +300 | 2 |  7 | 54.9 | 15.8 | 15.8 |
| +400 | 2 |  5 | 57.4 | 15.8 | 15.4 |
| +500 | 1 |  5 | 58.8 | 17.5 | 17.5 |
| +600 | 2 |  5 | 60.9 | 15.3 | 15.1 |

Oracle mean σ_BLUE = 16.1 ps; Oracle−global gap = 3.4 ps (mean).

---

## N_TOP scan (physical files, x=0, oracle k*(N), m=8, EJ-230)

| N_TOP | N_ch | k* | σ_TOP [ps] | σ_END [ps] | σ_BLUE [ps] | w_END |
|-------|------|----|-----------|-----------|------------|-------|
|  0 (END-only) | 16 | — | — | 50.5 | 50.5 | 1.000 |
|  4 | 20 | 1 | 78.2 | 52.1 | 45.3 | 0.765 |
|  8 | 24 | 1 | 32.0 | 50.7 | 32.2 | 0.438 |
| 14 | 30 | 2 | 18.7 | 51.1 | 18.4 | 0.173 |
| **20** | **36** | **3** | **15.2** | **53.7** | **15.0** | 0.086 |

---

## Position resolution

| Method | σ_x [mm] | Notes |
|--------|---------|-------|
| END Δt (m=8) | **7.9** | v_eff=177.6 mm/ns, no edge bias |
| END Δt (first photon) | 10.3 | v_eff=184.6 mm/ns (legacy) |
| TOP N_pe centroid LOO (linear) | 17.2 | +28.8 mm bias at bar ends |
| TOP N_pe centroid LOO (cubic) | 18.4 | +23.8 mm bias at bar ends |

---

## Design decision

EJ-230 + Vikuiti ESR (R=0.95) + 16 END + 20 TOP = **36 total channels**
σ_t (global, k=3, m=8) = **15.2 ps** at x=0; mean 19.5 ps across bar
σ_x (END Δt) = **7.9 mm**
SHiP target: σ_t < 100 ps → achieved with 6.6× margin.
