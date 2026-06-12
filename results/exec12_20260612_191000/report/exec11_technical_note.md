# EXEC_11 technical note: local Top-pair positioning

All values below are generated from `results/exec11_20260612_182454/analysis`; data commit `f431c01`.
Every resolution is a **simulation prediction — EJ-204 — intrinsic optical timing**.

## Method and geometry
The adjacent Top pair (28,29), centered near -442 mm, is reduced from per-hit
`sipm_hits` rows to fourth-hit times and PE-equivalent counts. The two references
were selected programmatically:

| reference | x_true_mm | criterion |
|---|---|---|
| POS_REF_1 | -441.000 | minimum abs(mu_dt_ps), tie -> minimum abs(x_true_mm) |
| POS_REF_2 | -443.000 | maximum mean_npe_A, tie -> minimum abs(x_true_mm) |

## QA and reconstruction
The original fits at -444, -442 and -439 mm converged to narrow subpeaks. Broad,
seeded, iterative v2 fits expose non-Gaussian tails:

| x_true_mm | sigma_dt_ps_v1 | sigma_dt_ps_v2 | chi2_ndf |
|---|---|---|---|
| -444.000 | 24.428 | 110.569 | 940.401 |
| -442.000 | 40.800 | 109.983 | 500.017 |
| -439.000 | 22.958 | 111.419 | 1081.012 |

| reference | method | bias_mm | sigma_x_mm | sigma_x_error_mm |
|---|---|---|---|---|
| POS_REF_1 | temporal | -0.258 | 12.494 | 0.161 |
| POS_REF_1 | ratio_linear | -1.698 | 5.558 | 0.072 |
| POS_REF_1 | BLUE | -1.502 | 5.265 | 0.068 |
| POS_REF_2 | temporal | 0.145 | 11.864 | 0.153 |
| POS_REF_2 | ratio_linear | 0.342 | 5.410 | 0.070 |
| POS_REF_2 | BLUE | 0.315 | 5.152 | 0.067 |

The count ratio is globally nonlinear; the cubic comparison creates ambiguous roots.
BLUE is covariance-aware. This is a one-pair empirical 1-D adaptation inspired by
pair timing, not Lv et al.'s literal 64-SiPM/2016-pair 2-D circle-intersection method.
It must not be directly compared with their 1.5–3 mm results.

## Limitations
Simulation only; one pair; two held-out test positions; no SPTR/electronics; strong
non-Gaussian timing tails; no experimental validation.
