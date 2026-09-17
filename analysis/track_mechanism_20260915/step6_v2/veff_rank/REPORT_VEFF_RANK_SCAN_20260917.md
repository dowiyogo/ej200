# EXEC_46 effective velocity rank/CFD tables

## Scope and provenance

| field | value |
|---|---|
| campaign | EXEC_46 v2 |
| production ROOT modification | none |
| simulation launch | none |
| NumPy | 2.0.2 |
| SciPy | 1.13.1; import PASS |
| partial cells | 21 |
| AD2 regenerated-cell control | PASS |
| AD2 maximum absolute field difference | 7.105427e-15 deg for theta fields; 0 for other reported numeric fields |
| CFD window [ns] | 8 |
| CFD bin [ns] | 0.02 |
| SPE pulse form | double exponential: `(1-exp(-t/tau_r))*exp(-t/tau_f)` |
| pulse model sources | `analysis/timing/pulse_models.py` for `fastic_measured` and `penarodriguez_shortened`; `(0.5,2.0) ns` fast contrast |
| MUSIC/SAMPIC model | NOT_AVAILABLE in `pulse_models.py` |

The `fastic_measured` pair `(2,3) ns` is documented as requiring confirmation with the complete FastIC+ chain. None of the three evaluated models represents MUSIC/SAMPIC. The comparison with 155 mm/ns is therefore indicative and is not a reproduction of the Betancourt MUSIC/SAMPIC estimator.

## W6 CFD fits

`mean_tdiff_ns` was obtained from the CFD timestamp at each END. Fits use the seven positions with weights from the per-cell timestamp SEM. `v_eff_two_end_mm_ns = 2/abs(slope)`. `chi2_flag=*` marks chi2/ndf > 5. The experimental values are references only.

| material | pulse_model | cfd_fraction | slope_ns_per_mm | slope_error_ns_per_mm | chi2 | ndf | chi2_ndf | chi2_flag | v_eff_two_end_mm_ns | v_eff_error_mm_ns | experimental_reference_mm_ns | theta_reference_experimental_deg |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| EJ-200 | fast_contrast | 0.14 | 0.011418 | 0.000001 | 6353.033053 | 5 | 1270.606611 | * | 175.163556 | 0.018355 | 155 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 0.011539 | 0.000001 | 9878.226897 | 5 | 1975.645379 | * | 173.317819 | 0.019358 | 155 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 0.011628 | 0.000001 | 14968.328855 | 5 | 2993.665771 | * | 171.999485 | 0.016964 | 155 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 0.011788 | 0.000001 | 16834.786337 | 5 | 3366.957267 | * | 169.666780 | 0.018173 | 155 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 0.011628 | 0.000001 | 14968.328855 | 5 | 2993.665771 | * | 171.999485 | 0.016964 | 155 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 0.011788 | 0.000001 | 16834.786337 | 5 | 3366.957267 | * | 169.666780 | 0.018173 | 155 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 0.011367 | 0.000001 | 4120.908114 | 5 | 824.181623 | * | 175.942718 | 0.019871 | 155 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 0.011469 | 0.000001 | 6116.801008 | 5 | 1223.360202 | * | 174.382679 | 0.021060 | 155 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 0.011548 | 0.000001 | 10701.472757 | 5 | 2140.294551 | * | 173.191084 | 0.018623 | 155 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 0.011679 | 0.000001 | 12921.577819 | 5 | 2584.315564 | * | 171.246298 | 0.019716 | 155 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 0.011548 | 0.000001 | 10701.472757 | 5 | 2140.294551 | * | 173.191084 | 0.018623 | 155 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 0.011679 | 0.000001 | 12921.577819 | 5 | 2584.315564 | * | 171.246298 | 0.019716 | 155 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 0.011258 | 0.000001 | 5268.165563 | 5 | 1053.633113 | * | 177.653283 | 0.021224 | 155 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 0.011367 | 0.000001 | 4406.480826 | 5 | 881.296165 | * | 175.944914 | 0.022810 | 155 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 0.011445 | 0.000001 | 7911.981666 | 5 | 1582.396333 | * | 174.748965 | 0.020630 | 155 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 0.011584 | 0.000001 | 9968.282457 | 5 | 1993.656491 | * | 172.645597 | 0.021506 | 155 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 0.011445 | 0.000001 | 7911.981666 | 5 | 1582.396333 | * | 174.748965 | 0.020630 | 155 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 0.011584 | 0.000001 | 9968.282457 | 5 | 1993.656491 | * | 172.645597 | 0.021506 | 155 | 35.4 |

## Per-cell CFD and transport tables

The checkpoint partial CSVs contain, by cell, model and fraction: event count used/excluded; median boundary encounters; median path length; Cherenkov fraction; and theta_eff q05/median/q95. `theta_eff = arccos(d_axial/path_length_mm)`.

Files: `partial/<cell_id>.csv`.

## Fixed configuration references

| item | value |
|---|---:|
| Betancourt effective velocity reference [mm/ns] | 155 |
| Betancourt theta_eff reference [deg] | 35.4 |
| Blondel center effective velocity reference [cm/ns] | 16.1 |
| Blondel x<20 cm effective velocity reference [cm/ns] | approximately 14 |
| simulated bar length [cm] | 140 |
| test-beam bar length [cm] | 150 |
| configured muon momentum | 1 GeV/c; beta 0.995424 |
| test-beam muon momentum | 2.5 GeV/c |
