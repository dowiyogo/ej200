# EXEC_14 — Fit Failures and Analysis Gaps

Generated: 2026-06-12

## Figures not regenerated for EJ-230

### exec09_tail_comparison.png

- **Origin**: `analysis/exec07/exec09_timing_mechanism.py`
- **Requires**: Paired "EndTop" + "End-only" ROOT files at x = −400, 0, +400 mm. The "End-only" campaign (no Top SiPMs) does not exist for EJ-230.
- **Action**: EJ-204 version used verbatim. Appears on two slides (main deck §timing mechanism, backup). The figure shows a purely geometric tail-suppression effect that is independent of scintillator material (same bar geometry, same Top-window cut).
- **Hook reserved**: `HOOK_EXEC09_EJ230` — regenerate if End-only EJ-230 runs become available.

### exec08b_window_dip_profiles.png / exec08b_id18_impact_maps.png

- **Origin**: `analysis/exec07/exec08b_window_dip.py`
- **Requires**: Window-dip scan ROOT files (`photon_hits_run_A/B/C1/C2_*.root`) in `/home/reriosto/SHiP/t0minidaq/sslg4/exec08b_window_dip/`. These runs target x ≈ −652…−648 mm with fine positioning around the inter-SiPM gap. No equivalent EJ-230 run was scheduled in EXEC_13.
- **Action**: EJ-204 versions used verbatim. The window-dip mechanism is a gap-geometry effect (SiPM spacing 20 mm) independent of scintillator material.
- **Hook reserved**: `HOOK_EXEC08B_EJ230`.

## Fit status (EXEC_14 analysis complete)

### exec10 Landau/Moyal fits

- 257 Moyal-Gaussian fits requested (mean ≥ 20 N_pe); **0 failed** (100 % convergence).
- 35 fits have σ_Gauss at the 0.1 N_pe lower bound (identifiability limit, flagged but not failed).
- 31 high-N_pe fits (mean > 50) with χ²/ndf > 5 — expected for End channels near beam (high counts, small Poisson fluctuations relative to the Landau model).

### exec12b / exec13 t_N consistency check

Gaussian fit on t_4 for End `near` SUM4 group at edge positions:

| x [mm] | group | σ_fit [ns] | Note |
|---:|---|---:|---|
| −690 | end_near | NaN | Low N_pe at far End SUM4 → degenerate fit |
| −650 | end_near | NaN | Same |
| −400 | end_near | 0.0212 | OK (21 ps) |
| 0 | end_near | 0.1279 | OK (128 ps) |
| +400 | end_near | 0.0169 | OK (17 ps) |
| +650 | end_near | NaN | Same as −650 |
| +690 | end_near | NaN | Same as −690 |

NaN fits at ±690 and ±650 are expected: the **far** End SUM4 clusters at these positions see < 30 photons on average, insufficient for a reliable Gaussian fit on t_4. The consistency check criterion (ratio in [0.3, 3.0] for positions with valid fits) passes for x = ±400, 0.
