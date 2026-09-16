# EXEC_46 Step 5 — chain-rule identification gate

Date: 2026-09-16

## Gate result

**HALTED_AT_5_1_MAJORITY_IDENTIFICATION_ARTIFACT.**

No simulation was run. The 210,000-event `derived_events` tree was opened read-only.
The common within-cell slope uses position fixed effects and cell-centered Npe. The
between-position slope is an unweighted OLS regression of the seven cell means; a WLS
sensitivity is retained in the sidecar.

## 5.1 within versus between

| material | beta_within [ps/pe] | beta_between [ps/pe] | gap significance | WLS beta_between [ps/pe] |
|---|---:|---:|---:|---:|
| EJ-200 | -0.05694 +/- 0.00076 | -0.02105 +/- 0.00199 | 16.84 sigma | -0.02095 |
| EJ-204 | -0.05423 +/- 0.00080 | -0.00874 +/- 0.00208 | 20.41 sigma | -0.00840 |
| EJ-230 | -0.05661 +/- 0.00089 | -0.00071 +/- 0.00198 | 25.72 sigma | -0.00024 |

The between-position response is much weaker than the conditional response
inside a cell. A local beta therefore does not identify the change of the
unconditional mean between positions.

## Quadratic summary and mandatory stop

| material | observed a2 | original predicted a2 | between predicted a2 | corrected residual a2 | corrected chi2/ndf | registered removed |
|---|---:|---:|---:|---:|---:|---:|
| EJ-200 | -71.81 | -194.91 | -75.51 | +3.70 +/- 7.37 | 76.33/5=15.27 | 97.00% |
| EJ-204 | -29.93 | -256.98 | -35.26 | +5.33 +/- 8.58 | 105.78/5=21.16 | 97.65% |
| EJ-230 | +2.65 | -194.58 | -2.69 | +5.34 +/- 7.67 | 90.08/5=18.02 | 97.29% |

The quadratic basis remains a poor shape description where its chi2/ndf is
large; the table is retained only to compare with the preregistered +123.10,
+227.05, and +197.23 ps/m^2 remnants. The primary result is pointwise:

| material | |x| [mm] | observed [ps] | original prediction [ps] | registered residual [ps] | between prediction [ps] | corrected residual [ps] | removed |
|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | 200 | -1.354 | -6.980 | +5.626 | -1.735 | +0.382 | +93.2% |
| EJ-200 | 500 | -6.264 | -44.246 | +37.982 | -12.632 | +6.368 | +83.2% |
| EJ-200 | 650 | -31.922 | -82.784 | +50.862 | -32.397 | +0.474 | +99.1% |
| EJ-204 | 200 | +0.306 | -10.061 | +10.367 | -0.753 | +1.059 | +89.8% |
| EJ-204 | 500 | +1.936 | -60.121 | +62.057 | -5.922 | +7.858 | +87.3% |
| EJ-204 | 650 | -14.307 | -109.642 | +95.335 | -15.253 | +0.946 | +99.0% |
| EJ-230 | 200 | +1.679 | -7.352 | +9.030 | -0.054 | +1.733 | +80.8% |
| EJ-230 | 500 | +6.889 | -46.831 | +53.720 | -0.450 | +7.338 | +86.3% |
| EJ-230 | 650 | +0.194 | -82.499 | +82.693 | -1.171 | +1.364 | +98.4% |

All three materials exceed the preregistered 50% majority threshold:
the identification correction removes 97.00--98.26% of the registered remnant across
the OLS primary result and WLS sensitivity. Independently, replacing the common
within slope by the common between slope explains 97.18--97.53% of its
own common-slope remnant. Both definitions cross the majority gate.
The corrected a2 uncertainty conservatively adds the between-slope uncertainty
without a covariance cancellation; it does not affect the gate.

This changes the thesis: the remnant is
predominantly an identification artifact before any attribution to optical
mechanisms. Per the explicit gate, Steps 5.2--5.5 and Step 6 were not run.

## Reproducibility

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step5.py
```

`within_between_slopes` and `identification_residual` each have PDF, CSV, ROOT,
and JSON metadata sidecars. No push, merge, deck edit, or simulation occurred.
