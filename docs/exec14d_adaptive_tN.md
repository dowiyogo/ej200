# EXEC_14D adaptive t_N

## HOOK_ADAPTIVE_TN

- `reach_min = 0.950` (CLI: `--reach-min`)
- `minimum_threshold = 4` PE (CLI: `--minimum-threshold`)
- Nominal thresholds: nearest-Top = 4 PE, near-End = 4 PE, far-End = 20 PE.
- `N_eff` is the largest integer threshold at or below `N_nominal` whose event reach is at least `reach_min`.
- If the minimum threshold does not reach `reach_min`, the result is marked as genuine starvation and no timing value is fabricated.

## Far-End threshold audit

| x [mm] | N nominal | N eff | mean Npe | R(nominal) | R(effective) | status |
|---:|---:|---:|---:|---:|---:|:---|
| -690 | 20 | 7 | 13.77 | 10.8% | 96.3% | reduced |
| -670 | 20 | 7 | 14.49 | 14.0% | 97.2% | reduced |
| -650 | 20 | 8 | 15.01 | 16.4% | 95.7% | reduced |
| -600 | 20 | 9 | 17.00 | 25.1% | 96.7% | reduced |
| -550 | 20 | 11 | 19.18 | 41.9% | 95.5% | reduced |
| -500 | 20 | 12 | 21.33 | 54.4% | 96.5% | reduced |
| -450 | 20 | 14 | 23.69 | 69.2% | 95.7% | reduced |
| -400 | 20 | 17 | 28.00 | 88.3% | 95.8% | reduced |
| -350 | 20 | 19 | 31.34 | 94.0% | 95.9% | reduced |
| -300 | 20 | 20 | 34.24 | 98.0% | 98.0% | nominal |
| -250 | 20 | 20 | 38.73 | 99.6% | 99.6% | nominal |
| -200 | 20 | 20 | 44.37 | 100.0% | 100.0% | nominal |
| -150 | 20 | 20 | 50.21 | 100.0% | 100.0% | nominal |
| -100 | 20 | 20 | 59.57 | 100.0% | 100.0% | nominal |
| -50 | 20 | 20 | 66.61 | 100.0% | 100.0% | nominal |
| 0 | 20 | 20 | 76.05 | 100.0% | 100.0% | nominal |
| 50 | 20 | 20 | 66.17 | 100.0% | 100.0% | nominal |
| 100 | 20 | 20 | 58.99 | 100.0% | 100.0% | nominal |
| 150 | 20 | 20 | 50.29 | 100.0% | 100.0% | nominal |
| 200 | 20 | 20 | 44.41 | 99.9% | 99.9% | nominal |
| 250 | 20 | 20 | 38.77 | 99.7% | 99.7% | nominal |
| 300 | 20 | 20 | 34.04 | 97.9% | 97.9% | nominal |
| 350 | 20 | 20 | 31.51 | 95.5% | 95.5% | nominal |
| 400 | 20 | 17 | 27.77 | 88.5% | 96.8% | reduced |
| 450 | 20 | 14 | 23.64 | 70.5% | 96.0% | reduced |
| 500 | 20 | 12 | 20.99 | 52.6% | 96.3% | reduced |
| 550 | 20 | 11 | 19.19 | 41.6% | 95.4% | reduced |
| 600 | 20 | 9 | 16.74 | 25.6% | 96.0% | reduced |
| 650 | 20 | 8 | 15.16 | 16.9% | 96.5% | reduced |
| 670 | 20 | 7 | 14.15 | 12.6% | 97.3% | reduced |
| 690 | 20 | 7 | 13.85 | 11.9% | 96.4% | reduced |

## Forward note

A future cross-sample study must rerun this same reusable hook on every compared data set. Comparing a fixed-t20 series with an adaptive-threshold series would not be apples-to-apples.
