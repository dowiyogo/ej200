# EXEC_14B missing-analysis audit

Initial audit date: 2026-06-13 CEST

| Frame | Figure | Script | Data required | State | Action |
|---:|---|---|---|---|---|
| 8, 12, 16, 20, 24, 28, 32 | `exec13_tn_<x>mm_{top,endL,endR}.png` | `analysis/exec07/exec12b_tn_dispersion.py` | Main EJ-230 scan | Closed: 21 authentic, distinct EJ-230 figures | Repaired all 21 explicit paths and verified figure/frame consistency |
| 44 | `exec08b_window_dip_profiles.png` | `analysis/exec07/exec08b_window_dip.py` | Existing x=-650 plus directed x=-652,-642,-648,-654 EJ-230 runs | Closed: authentic EJ-230 figure and CSV | Ran four directed simulations and regenerated; SHA differs from EJ-204 |
| 113 | `exec08b_id18_impact_maps.png` | `analysis/exec07/exec08b_window_dip.py` | Same five EJ-230 runs | Closed: authentic EJ-230 figure and CSV | Regenerated from directed EJ-230 simulations; SHA differs from EJ-204 |
| 46, 117 | `exec09_tail_comparison.png` | `analysis/exec07/exec09_timing_mechanism.py` | Main EndTop scan plus End-only x=-400,0,+400 EJ-230 runs | Closed: authentic EJ-230 figure, CSV, and tables | Ran/validated three End-only simulations and regenerated; SHA differs from EJ-204 |

All other referenced raster figures resolve to products generated from the main
EJ-230 scan. Final audits report `156/156` references resolved, 134 unique PNGs,
zero forbidden raster matches, 50 CSV-derived table inputs, and consistent
position/Top/endL/endR assignments.
