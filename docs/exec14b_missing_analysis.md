# EXEC_14B missing-analysis audit

Initial audit date: 2026-06-13 CEST

| Frame | Figure | Script | Data required | State | Action |
|---:|---|---|---|---|---|
| 8, 12, 16, 20, 24, 28, 32 | `exec13_tn_<x>mm_{top,endL,endR}.png` | `analysis/exec07/exec12b_tn_dispersion.py` | Main EJ-230 scan | Authentic figures exist; TeX paths contain literal `\1`/`\2` | Repair all 21 explicit paths |
| 44 | `exec08b_window_dip_profiles.png` | `analysis/exec07/exec08b_window_dip.py` | Existing x=-650 plus directed x=-652,-642,-648,-654 EJ-230 runs | Existing result is byte-identical to forbidden EJ-204 figure | Run directed EJ-230 simulations and regenerate |
| 113 | `exec08b_id18_impact_maps.png` | `analysis/exec07/exec08b_window_dip.py` | Same five EJ-230 runs | Existing result is byte-identical to forbidden EJ-204 figure | Run directed EJ-230 simulations and regenerate |
| 46, 117 | `exec09_tail_comparison.png` | `analysis/exec07/exec09_timing_mechanism.py` | Main EndTop scan plus End-only x=-400,0,+400 EJ-230 runs | Existing result is byte-identical to forbidden EJ-204 figure | Run End-only EJ-230 simulations and regenerate |

All other referenced raster figures resolve to products generated from the main
EJ-230 scan. Their CSV support and frame consistency remain subject to the final
number-provenance and PDF preflight audits.
