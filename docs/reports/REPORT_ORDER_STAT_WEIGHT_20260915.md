# Corrected order-statistics weighting analysis

**Date:** 2026-09-15  
**Scope:** revised ad-hoc analysis of the existing `derived_events` tree; no simulation  
**Status:** `COMPLETE` — corrected L2 guard passed 21/21; Part B executed

## Physics result

The single-end penalty is strongly position dependent. At x=0, the robust ratio sigma(tL)/sigma(T0) is 1.378, 1.379 and 1.395 for EJ-200, EJ-204 and EJ-230, respectively, reasonably near the independent/equal-width reference sqrt(2)=1.414. Near the left end, the same ratio falls to 0.130, 0.188 and 0.231 at x=-650 mm: the nearby single end is substantially narrower than the naive two-end average because the latter still includes the far end. The right-end ratios show the reflected pattern.

The negative within-cell correlation with total END light is present at every position, but it does **not** quantitatively reproduce the measured even T0(x) structure. Integrating the measured local TProfile slopes over the measured seven-point Npe_END(x) profile predicts quadratic coefficients of -194.91, -256.98 and -194.58 ps/m2. The measured coefficients are -71.81, -29.93 and +2.65 ps/m2. Thus the signed predicted/measured ratios are 2.71, 8.59 and -73.42. The estimate overshoots EJ-200/204 and gets the sign wrong for EJ-230; it cannot be called a fraction successfully explained.

This failure is informative rather than evidence that rho is spurious. A within-cell regression slope is conditional at fixed position and does not by itself identify the change of the unconditional mean between positions. The changing intercept, path-length and spectral populations remain unresolved by this tree. In addition, 12/21 linear TProfile fits have chi2/ndf > 5, so a single local linear slope is often an inadequate summary.

## Corrected guardrail

The L2 triangle inequality applies to ordinary standard deviations, not to IQR widths or fitted Gaussian-core widths. The hard gate was therefore evaluated with event-level RMS values:

`RMS(T0) <= [RMS(tL)+RMS(tR)]/2 + 3*SE_gap`.

`SE_gap` is a paired influence-function delta-method uncertainty. It retains the covariance among tL, tR and T0 from the same event. All 21 cells pass. The closest cases remain well inside the bound:

| Material | closest x [mm] | T0 minus triangle bound [ps] | 3 SE_gap [ps] | Status |
|---|---:|---:|---:|---|
| EJ-200 | -650 | -8.086 | 0.337 | PASS |
| EJ-204 | -650 | -8.340 | 0.328 | PASS |
| EJ-230 | +650 | -8.555 | 0.367 | PASS |

The tree also satisfies `T0=(tL+tR)/2` exactly: the maximum absolute event-level difference is 0 ns.

## Part A — full single-end ratios

Widths use IQR/1.349 and a ROOT binned Gaussian initialized at the histogram peak. Histograms have 4 ps bins over median +/-0.320 ns; the fit window is median +/-2*q68 with q68=(q84-q16)/2. `*` marks a Gaussian fit with chi2/ndf > 5 or failed status/covariance. No unreliable fit is silently treated as valid.

### EJ-200

| x [mm] | L/T0 IQR | R/T0 IQR | L/T0 Gaussian | R/T0 Gaussian | Gaussian reliability L/R/T0 |
|---:|---:|---:|---:|---:|---|
| -650 | 0.1299 | 1.9555 | 0.0294 | 1.9417 | */OK/OK |
| -500 | 0.7515 | 1.8022 | 0.7648 | 1.7530 | OK/OK/OK |
| -200 | 1.1626 | 1.5949 | 1.1730 | 1.5640 | OK/OK/OK |
| +0 | 1.3784 | 1.3526 | 1.3831 | 1.3697 | OK/OK/OK |
| +200 | 1.5408 | 1.1391 | 1.5232 | 1.1205 | OK/OK/OK |
| +500 | 1.8109 | 0.7395 | 1.7823 | 0.7611 | OK/OK/OK |
| +650 | 1.9652 | 0.1389 | 1.9370 | 0.0281 | OK/*/OK |

### EJ-204

| x [mm] | L/T0 IQR | R/T0 IQR | L/T0 Gaussian | R/T0 Gaussian | Gaussian reliability L/R/T0 |
|---:|---:|---:|---:|---:|---|
| -650 | 0.1884 | 1.9579 | 0.1261 | 1.9466 | */OK/OK |
| -500 | 0.6727 | 1.8570 | 0.6863 | 1.8135 | OK/OK/OK |
| -200 | 1.0881 | 1.5757 | 1.1019 | 1.5709 | OK/OK/OK |
| +0 | 1.3793 | 1.3748 | 1.3776 | 1.3866 | OK/OK/OK |
| +200 | 1.6268 | 1.1383 | 1.6423 | 1.1221 | OK/OK/OK |
| +500 | 1.8616 | 0.6737 | 1.8143 | 0.6801 | OK/OK/OK |
| +650 | 1.9773 | 0.1841 | 1.9940 | 0.1341 | OK/*/OK |

### EJ-230

| x [mm] | L/T0 IQR | R/T0 IQR | L/T0 Gaussian | R/T0 Gaussian | Gaussian reliability L/R/T0 |
|---:|---:|---:|---:|---:|---|
| -650 | 0.2313 | 1.9731 | 0.1481 | 1.9638 | */OK/OK |
| -500 | 0.6328 | 1.8932 | 0.6494 | 1.8151 | OK/OK/OK |
| -200 | 1.0752 | 1.6593 | 1.0645 | 1.5937 | OK/OK/OK |
| +0 | 1.3949 | 1.3869 | 1.3599 | 1.3513 | OK/OK/OK |
| +200 | 1.6152 | 1.0866 | 1.5968 | 1.0761 | OK/OK/OK |
| +500 | 1.8807 | 0.6090 | 1.8575 | 0.6248 | OK/OK/OK |
| +650 | 1.9655 | 0.2249 | 1.9399 | 0.1624 | OK/*/OK |

### Practical near-left boundary

The citable boundary uses the robust IQR ratio and threshold `sigma(tL)/sigma(T0)<0.5`. Only x=-650 mm is below threshold in the sampled grid. The quoted crossing is a linear interpolation between -650 and -500 mm, so the bracket is the primary result and the interpolated x is explicitly approximate.

| Material | ratio at -650 | ratio at -500 | sampled crossing bracket [mm] | interpolated x_cross [mm] | distance from left end [mm] |
|---|---:|---:|---|---:|---:|
| EJ-200 | 0.1299 | 0.7515 | (-650, -500) | -560.7 | 139.3 |
| EJ-204 | 0.1884 | 0.6727 | (-650, -500) | -553.5 | 146.5 |
| EJ-230 | 0.2313 | 0.6328 | (-650, -500) | -549.6 | 150.4 |

The Gaussian interpolations are -554.0, -549.9 and -544.7 mm, but all three near-left Gaussian fits at x=-650 are unreliable; they are retained only in the CSV.

### x=0 symmetry sanity flags

All nine non-aborting checks satisfy sigma(T0)<=sigma(tL) and sigma(T0)<=sigma(tR). Values are in ps.

| Material | Method | sigma(tL) | sigma(tR) | sigma(T0) | Status |
|---|---|---:|---:|---:|---|
| EJ-200 | RMS | 105.972 | 105.335 | 76.912 | PASS |
| EJ-200 | IQR | 108.075 | 106.052 | 78.408 | PASS |
| EJ-200 | Gaussian | 106.618 | 105.582 | 77.084 | PASS |
| EJ-204 | RMS | 100.039 | 99.634 | 72.439 | PASS |
| EJ-204 | IQR | 101.010 | 100.677 | 73.233 | PASS |
| EJ-204 | Gaussian | 99.309 | 99.960 | 72.089 | PASS |
| EJ-230 | RMS | 93.655 | 91.997 | 67.772 | PASS |
| EJ-230 | IQR | 93.238 | 92.702 | 66.842 | PASS |
| EJ-230 | Gaussian | 92.789 | 92.199 | 68.231 | PASS |

### Gaussian fit diagnostics

Each entry is chi2/ndf for tL, tR and T0. `*` denotes the six fits rejected by the preregistered reliability threshold of 5.

| Material | x [mm] | tL chi2/ndf | tR chi2/ndf | T0 chi2/ndf |
|---|---:|---:|---:|---:|
| EJ-200 | -650 | 205.15* | 1.38 | 1.71 |
| EJ-200 | -500 | 1.17 | 1.22 | 1.01 |
| EJ-200 | -200 | 0.91 | 1.36 | 1.24 |
| EJ-200 | +0 | 1.06 | 0.95 | 1.40 |
| EJ-200 | +200 | 0.87 | 1.04 | 1.10 |
| EJ-200 | +500 | 1.17 | 1.09 | 1.33 |
| EJ-200 | +650 | 1.17 | 210.11* | 1.38 |
| EJ-204 | -650 | 460.17* | 1.42 | 1.52 |
| EJ-204 | -500 | 1.28 | 1.71 | 1.31 |
| EJ-204 | -200 | 1.44 | 1.51 | 1.16 |
| EJ-204 | +0 | 0.98 | 1.25 | 1.01 |
| EJ-204 | +200 | 1.04 | 1.35 | 0.70 |
| EJ-204 | +500 | 1.23 | 1.32 | 1.36 |
| EJ-204 | +650 | 1.46 | 503.74* | 2.06 |
| EJ-230 | -650 | 358.93* | 1.67 | 2.38 |
| EJ-230 | -500 | 0.64 | 1.63 | 1.60 |
| EJ-230 | -200 | 1.29 | 1.21 | 1.33 |
| EJ-230 | +0 | 1.05 | 1.09 | 1.11 |
| EJ-230 | +200 | 1.47 | 1.21 | 1.29 |
| EJ-230 | +500 | 1.55 | 1.39 | 1.89 |
| EJ-230 | +650 | 1.78 | 356.66* | 2.68 |

The six rejected Gaussian fits are the direct near-end timestamps: tL at x=-650 and tR at x=+650 for each material. All 21 T0 Gaussian fits are reliable by this threshold.

## Part B — per-cell Npe dependence

Each slope comes from a 40-bin `TProfile` of T0 versus Npe_END and a ROOT `pol1` fit over the observed Npe range. Slopes are shown in ps/PE. `*` marks chi2/ndf > 5; formal fit errors for those rows should not be interpreted as complete model uncertainties.

| Material | x [mm] | mean Npe_END | rho(T0,Npe_END) | d<T0>/dNpe [ps/PE] | chi2/ndf |
|---|---:|---:|---:|---:|---:|
| EJ-200 | -650 | 2589.01 | -0.2767 | -0.03224 +/- 0.00091 | 3.11 |
| EJ-200 | -500 | 1657.55 | -0.3604 | -0.02975 +/- 0.00052 | 33.25* |
| EJ-200 | -200 | 1140.15 | -0.3960 | -0.09231 +/- 0.00171 | 8.52* |
| EJ-200 | +0 | 1056.81 | -0.3926 | -0.07459 +/- 0.00106 | 16.57* |
| EJ-200 | +200 | 1138.34 | -0.3911 | -0.09719 +/- 0.00198 | 7.20* |
| EJ-200 | +500 | 1656.23 | -0.3493 | -0.06871 +/- 0.00134 | 4.13 |
| EJ-200 | +650 | 2602.54 | -0.2967 | -0.03348 +/- 0.00093 | 4.22 |
| EJ-204 | -650 | 2546.45 | -0.2760 | -0.03500 +/- 0.00107 | 4.06 |
| EJ-204 | -500 | 1472.88 | -0.3131 | -0.06364 +/- 0.00156 | 4.10 |
| EJ-204 | -200 | 882.57 | -0.3532 | -0.11086 +/- 0.00265 | 5.13* |
| EJ-204 | +0 | 796.28 | -0.3456 | -0.12254 +/- 0.00280 | 3.83 |
| EJ-204 | +200 | 882.24 | -0.3652 | -0.11133 +/- 0.00236 | 5.60* |
| EJ-204 | +500 | 1474.54 | -0.3270 | -0.05282 +/- 0.00123 | 12.03* |
| EJ-204 | +650 | 2535.60 | -0.2670 | -0.03414 +/- 0.00105 | 3.54 |
| EJ-230 | -650 | 2249.60 | -0.2555 | -0.02933 +/- 0.00048 | 4.38 |
| EJ-230 | -500 | 1237.93 | -0.3019 | -0.06855 +/- 0.00082 | 4.10 |
| EJ-230 | -200 | 684.29 | -0.3258 | -0.09656 +/- 0.00219 | 10.34* |
| EJ-230 | +0 | 608.50 | -0.3511 | -0.09666 +/- 0.00241 | 18.87* |
| EJ-230 | +200 | 684.24 | -0.3373 | -0.09827 +/- 0.00233 | 11.40* |
| EJ-230 | +500 | 1240.70 | -0.2841 | -0.02114 +/- 0.00039 | 22.81* |
| EJ-230 | +650 | 2252.03 | -0.2518 | -0.02204 +/- 0.00070 | 9.18* |

The per-cell rho ranges are EJ-200 [-0.3960,-0.2767], EJ-204 [-0.3652,-0.2670], and EJ-230 [-0.3511,-0.2518]. Correlation strength is largest around x=0 to +/-200 mm and weakens toward the ends.

## Chain-rule estimate

No `physics-learnings.md` file was present under `/home/rrios` at analysis time. To avoid inventing or fitting an attenuation model, the calculation uses the measured seven-point Npe_END means directly. Mirrored means and slopes are averaged to isolate the even component, and `dT=beta dN` is integrated from x=0 with a trapezoidal rule through |x|={0,200,500,650} mm. This preserves all components present in the measured light profile and introduces no single-exponential refit.

| Material | |x| [mm] | observed even T0 shift [ps] | Npe-only predicted shift [ps] |
|---|---:|---:|---:|
| EJ-200 | 0 | +0.000 | +0.000 |
| EJ-200 | 200 | -1.354 | -6.980 |
| EJ-200 | 500 | -6.264 | -44.246 |
| EJ-200 | 650 | -31.922 | -82.784 |
| EJ-204 | 0 | +0.000 | +0.000 |
| EJ-204 | 200 | +0.306 | -10.061 |
| EJ-204 | 500 | +1.936 | -60.121 |
| EJ-204 | 650 | -14.307 | -109.642 |
| EJ-230 | 0 | +0.000 | +0.000 |
| EJ-230 | 200 | +1.679 | -7.352 |
| EJ-230 | 500 | +6.889 | -46.831 |
| EJ-230 | 650 | +0.194 | -82.499 |

The predicted and observed seven-point curves were fitted with the same position weights and even-quadratic basis. The recomputed observed coefficients reproduce the previous report within 1e-8 ns/m2.

| Material | measured a2 [ps/m2] | predicted a2 [ps/m2] | signed predicted/measured | measured minus predicted [ps/m2] | Assessment |
|---|---:|---:|---:|---:|---|
| EJ-200 | -71.81 +/- 1.80 | -194.91 | +2.71 (+271%) | +123.10 | same sign, overprediction |
| EJ-204 | -29.93 +/- 1.79 | -256.98 | +8.59 (+859%) | +227.05 | same sign, overprediction |
| EJ-230 | +2.65 +/- 1.75 | -194.58 | -73.42 (-7342%) | +197.23 | opposite sign |

For EJ-230 the measured quadratic coefficient is only 1.52 standard deviations from zero, so the -73.42 ratio is ill-conditioned as a percentage as well as opposite in sign. It is listed because the requested fraction must be explicit.

An even-quartic cross-check also disagrees in shape:

| Material | observed a2 [ps/m2] | predicted a2 [ps/m2] | observed a4 [ps/m4] | predicted a4 [ps/m4] |
|---|---:|---:|---:|---:|
| EJ-200 | +54.53 | -147.64 | -296.84 | -111.07 |
| EJ-204 | +72.02 | -210.24 | -242.58 | -111.23 |
| EJ-230 | +67.97 | -175.38 | -156.84 | -46.08 |

The position-fit diagnostics are listed below. The same `chi2/ndf > 5`
reliability flag is applied here; `*` marks an unreliable fit. The predicted
quartic curves pass this descriptive threshold, but their coefficients still
do not match the observed quartic coefficients.

| Material | observed quadratic chi2/ndf | predicted quadratic chi2/ndf | observed quartic chi2/ndf | predicted quartic chi2/ndf |
|---|---:|---:|---:|---:|
| EJ-200 | 313.07/5 = 62.61* | 43.14/5 = 8.63* | 11.95/4 = 2.99 | 0.98/4 = 0.24 |
| EJ-204 | 208.76/5 = 41.75* | 44.76/5 = 8.95* | 10.03/4 = 2.51 | 2.99/4 = 0.75 |
| EJ-230 | 96.33/5 = 19.27* | 7.66/5 = 1.53 | 8.81/4 = 2.20 | 0.11/4 = 0.03 |

The chain-rule coefficients and ratios are point estimates. No uncertainty is attached because 12/21 local linear profiles fail the stated fit-quality threshold and because the slopes and Npe means are correlated within each cell. Propagating only the ROOT parameter errors as independent Gaussian errors would give a misleadingly precise result.

## Interpretation and limitations

- Part A quantifies the estimator effect itself: around the center, using one end costs about 38-40% in robust width; within the last sampled 50 mm at the near end it improves the width by 77-87% relative to T0.
- Part B does not identify a causal fraction. It tests whether the observed conditional light-time slopes, transported along the measured attenuation profile, reproduce the mean-position curve. They do not.
- The unresolved remainder includes changing conditional intercepts and path/spectral populations. The current tree has no accumulated path length, boundary count, creation wavelength, or track identifier linking transport history to the detected first hit.
- The linearly interpolated 0.5 crossing is a sparse-grid convenience, not a directly sampled detector boundary.

## Reproduction and artifacts

Input: `/home/rrios/ej200/analysis/tsum_veff_20260914/sources/tsum_veff.root`, opened `READ`, SHA-256 `96fc575b583cad37f523dcab7bc8c1c0873f258e39fdb375bfa5739bf7ac78af`. The hash was unchanged after analysis.

```bash
cd /home/rrios/ej200/analysis/order_stat_weight_20260915
root -l -b -q 'macros/order_stat_weight.C+'
```

ROOT version: `6.40.02`. The command completed with exit code 0.

- Macro: `/home/rrios/ej200/analysis/order_stat_weight_20260915/macros/order_stat_weight.C`
- Complete ROOT objects: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/order_stat_weight.root`
- Part A widths and all fit diagnostics: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_a_widths.csv`
- Corrected guards: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_a_guardrails.csv`
- x=0 flags: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_a_x0_sanity.csv`
- Threshold crossings: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_a_thresholds.csv`
- Per-cell profiles: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_b_profiles.csv`
- Chain-rule points and summary: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_b_chain_rule_points.csv`, `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/part_b_chain_rule_summary.csv`
- Metadata and validation: `/home/rrios/ej200/analysis/order_stat_weight_20260915/sources/order_stat_weight.meta.json`, `/home/rrios/ej200/analysis/order_stat_weight_20260915/logs/validate_outputs.log`

No figures were produced, so no per-figure sidecars were required. The numeric analysis has ROOT, CSV and JSON sidecars.
