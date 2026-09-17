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





## AE1: theta_eff tables

`theta_eff = arccos(d_axial/path_length_mm)` for the threshold-crossing photon. The aggregate table uses an events-weighted mean of per-cell quantiles. The per-position table is listed separately. Reference: 35.4 degrees experimental.

| material | pulse_model | cfd_fraction | extreme | theta_q05_weighted_mean_deg | theta_median_weighted_mean_deg | theta_q95_weighted_mean_deg | reference_experimental_theta_deg |
|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | L | 7.7366021 | 29.026505 | 45.840462 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | R | 7.7846398 | 29.027476 | 45.96484 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | L | 8.4509615 | 31.41153 | 49.289962 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | R | 8.3176478 | 31.451229 | 49.253082 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | L | 8.9392557 | 33.220022 | 50.97341 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | R | 8.8395955 | 33.241699 | 50.974836 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | L | 9.7866055 | 35.797037 | 53.81378 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | R | 9.8738918 | 35.821147 | 53.797892 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | L | 8.9392557 | 33.220022 | 50.97341 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | R | 8.8395955 | 33.241699 | 50.974836 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | L | 9.7866055 | 35.797037 | 53.81378 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | R | 9.8738918 | 35.821147 | 53.797892 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | L | 7.5614904 | 28.643283 | 45.262623 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | R | 7.5302248 | 28.696233 | 45.319612 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | L | 8.1923333 | 30.952606 | 48.134306 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | R | 8.1485476 | 30.896731 | 48.138095 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | L | 8.9690642 | 32.943654 | 50.270948 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | R | 8.8435776 | 32.997526 | 50.353777 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | L | 9.7996451 | 35.380857 | 52.939333 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | R | 9.7634963 | 35.460881 | 53.126754 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | L | 8.9690642 | 32.943654 | 50.270948 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | R | 8.8435776 | 32.997526 | 50.353777 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | L | 9.7996451 | 35.380857 | 52.939333 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | R | 9.7634963 | 35.460881 | 53.126754 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | L | 7.4029968 | 28.009239 | 43.167054 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | R | 7.415269 | 27.962944 | 43.183931 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | L | 8.1943829 | 30.556904 | 46.408515 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | R | 8.2535405 | 30.522227 | 46.445138 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | L | 8.9822524 | 32.439419 | 49.056558 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | R | 8.9334763 | 32.555695 | 49.145487 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | L | 9.9079388 | 35.44104 | 52.028209 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | R | 9.8747689 | 35.57087 | 52.019677 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | L | 8.9822524 | 32.439419 | 49.056558 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | R | 8.9334763 | 32.555695 | 49.145487 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | L | 9.9079388 | 35.44104 | 52.028209 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | R | 9.8747689 | 35.57087 | 52.019677 | 35.4 |

### By position

| material | pulse_model | cfd_fraction | x_mm | theta_q05_weighted_mean_deg | theta_median_weighted_mean_deg | theta_q95_weighted_mean_deg | reference_experimental_theta_deg |
|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | -650 | 8.6082544 | 33.285941 | 51.754756 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | -500 | 7.8241957 | 28.618668 | 44.841476 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | -200 | 7.2126818 | 26.755402 | 42.750268 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 0 | 6.9007699 | 25.901081 | 42.194374 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 200 | 7.1595421 | 26.518639 | 42.799105 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 500 | 7.8165687 | 28.819247 | 45.085362 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 650 | 8.8023338 | 33.289953 | 51.893214 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -650 | 9.6462231 | 35.417027 | 55.633492 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -500 | 8.5457198 | 31.64487 | 49.70876 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -200 | 7.6309519 | 28.715513 | 44.919222 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 0 | 7.4288117 | 28.032217 | 44.234215 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 200 | 7.7136185 | 28.731417 | 44.919497 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 500 | 8.4346154 | 31.805854 | 49.577151 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 650 | 9.2901922 | 35.672761 | 55.90832 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -650 | 9.7495108 | 37.483987 | 57.885265 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -500 | 9.2987585 | 33.144667 | 50.631878 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -200 | 8.3497717 | 30.682522 | 46.911827 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 0 | 8.1725122 | 29.955731 | 45.839363 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 200 | 7.8485452 | 30.71314 | 46.848012 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 500 | 8.6790061 | 33.298753 | 50.862668 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 650 | 10.127875 | 37.337224 | 57.839847 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -650 | 10.813234 | 39.72015 | 60.522069 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -500 | 9.6693485 | 35.46073 | 53.627948 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -200 | 9.0749984 | 33.378409 | 49.453335 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 0 | 9.0518171 | 33.342908 | 49.395318 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 200 | 9.3924102 | 33.619042 | 49.524946 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 500 | 10.031676 | 35.524344 | 53.724444 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 650 | 10.778256 | 39.618063 | 60.392789 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -650 | 9.7495108 | 37.483987 | 57.885265 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -500 | 9.2987585 | 33.144667 | 50.631878 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -200 | 8.3497717 | 30.682522 | 46.911827 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 0 | 8.1725122 | 29.955731 | 45.839363 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 200 | 7.8485452 | 30.71314 | 46.848012 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 500 | 8.6790061 | 33.298753 | 50.862668 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 650 | 10.127875 | 37.337224 | 57.839847 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -650 | 10.813234 | 39.72015 | 60.522069 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -500 | 9.6693485 | 35.46073 | 53.627948 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -200 | 9.0749984 | 33.378409 | 49.453335 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 0 | 9.0518171 | 33.342908 | 49.395318 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 200 | 9.3924102 | 33.619042 | 49.524946 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 500 | 10.031676 | 35.524344 | 53.724444 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 650 | 10.778256 | 39.618063 | 60.392789 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -650 | 8.5581069 | 33.022776 | 51.424399 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -500 | 7.4836333 | 28.165251 | 44.307197 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -200 | 6.8815551 | 26.075494 | 42.100616 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 0 | 6.7969517 | 26.169823 | 41.259768 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 200 | 7.2682072 | 26.049822 | 42.1072 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 500 | 7.3402288 | 28.21591 | 44.286558 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 650 | 8.4923202 | 32.989231 | 51.552086 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -650 | 9.0438848 | 34.928941 | 53.892587 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -500 | 8.2539167 | 31.105921 | 48.310546 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -200 | 7.7165879 | 28.649875 | 44.486917 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 0 | 7.3132339 | 27.664678 | 43.567046 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 200 | 7.4549502 | 28.512651 | 44.408527 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 500 | 8.181749 | 30.835643 | 48.246315 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 650 | 9.2287604 | 34.77497 | 54.041465 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -650 | 9.7776963 | 37.035527 | 57.243299 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -500 | 9.0832115 | 33.131433 | 50.003545 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -200 | 8.2582909 | 30.143069 | 46.14774 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 0 | 7.97589 | 29.645209 | 45.172694 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 200 | 8.2422118 | 30.122233 | 46.24005 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 500 | 9.1802551 | 33.41202 | 50.196963 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 650 | 9.8266906 | 37.304639 | 57.182249 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -650 | 10.642705 | 39.118262 | 59.615057 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -500 | 10.155297 | 35.19176 | 53.101709 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -200 | 9.1062125 | 33.023896 | 48.661922 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 0 | 8.7339051 | 32.725268 | 48.52307 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 200 | 9.4177658 | 33.23936 | 48.538957 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 500 | 9.8534566 | 35.460902 | 53.190748 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 650 | 10.561652 | 39.186637 | 59.599841 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -650 | 9.7776963 | 37.035527 | 57.243299 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -500 | 9.0832115 | 33.131433 | 50.003545 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -200 | 8.2582909 | 30.143069 | 46.14774 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 0 | 7.97589 | 29.645209 | 45.172694 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 200 | 8.2422118 | 30.122233 | 46.24005 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 500 | 9.1802551 | 33.41202 | 50.196963 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 650 | 9.8266906 | 37.304639 | 57.182249 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -650 | 10.642705 | 39.118262 | 59.615057 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -500 | 10.155297 | 35.19176 | 53.101709 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -200 | 9.1062125 | 33.023896 | 48.661922 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 0 | 8.7339051 | 32.725268 | 48.52307 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 200 | 9.4177658 | 33.23936 | 48.538957 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 500 | 9.8534566 | 35.460902 | 53.190748 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 650 | 10.561652 | 39.186637 | 59.599841 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -650 | 8.4784601 | 32.674365 | 50.364501 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -500 | 7.4769597 | 27.420756 | 41.986474 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -200 | 6.5922344 | 24.999925 | 38.772142 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 0 | 6.9386196 | 26.421932 | 40.14821 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 200 | 6.576783 | 24.670214 | 38.664774 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 500 | 7.5825474 | 27.436287 | 41.915497 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 650 | 8.2183261 | 32.279162 | 50.37685 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -650 | 8.9934605 | 34.530712 | 52.706544 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -500 | 8.0979114 | 30.134079 | 45.338325 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -200 | 7.6175891 | 28.269477 | 43.105238 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 0 | 7.5712976 | 27.702761 | 42.430877 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 200 | 7.6282124 | 28.139005 | 43.167798 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 500 | 8.5579501 | 30.545355 | 45.735703 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 650 | 9.1013107 | 34.455565 | 52.503299 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -650 | 9.8631104 | 36.326526 | 54.981897 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -500 | 8.9129306 | 32.278753 | 49.554399 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -200 | 8.4348117 | 30.218862 | 45.248173 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 0 | 7.9765907 | 29.608043 | 44.6177 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 200 | 8.5599765 | 30.298428 | 45.146855 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 500 | 8.982745 | 32.38371 | 49.081651 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 650 | 9.9748858 | 36.368579 | 55.076484 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -650 | 10.542006 | 39.162387 | 58.60896 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -500 | 9.9471708 | 35.394844 | 52.001076 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -200 | 9.475492 | 33.25861 | 47.947212 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 0 | 9.0766177 | 32.564346 | 47.198372 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 200 | 9.3858302 | 33.230158 | 47.800758 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 500 | 9.6938662 | 35.254301 | 51.904221 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 650 | 11.118494 | 39.677043 | 58.706998 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -650 | 9.8631104 | 36.326526 | 54.981897 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -500 | 8.9129306 | 32.278753 | 49.554399 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -200 | 8.4348117 | 30.218862 | 45.248173 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 0 | 7.9765907 | 29.608043 | 44.6177 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 200 | 8.5599765 | 30.298428 | 45.146855 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 500 | 8.982745 | 32.38371 | 49.081651 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 650 | 9.9748858 | 36.368579 | 55.076484 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -650 | 10.542006 | 39.162387 | 58.60896 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -500 | 9.9471708 | 35.394844 | 52.001076 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -200 | 9.475492 | 33.25861 | 47.947212 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 0 | 9.0766177 | 32.564346 | 47.198372 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 200 | 9.3858302 | 33.230158 | 47.800758 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 500 | 9.6938662 | 35.254301 | 51.904221 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 650 | 11.118494 | 39.677043 | 58.706998 | 35.4 |

## AE2: geometry table

The exact event-level `<path>/<d_axial>` was not persisted; it is `NOT_AVAILABLE`. The proxy `1/cos(theta_median)` and implied velocity are tabulated without compatibility calculations.

| material | pulse_model | cfd_fraction | extreme | path_over_daxial_exact | v_group_MPT_mm_ns | v_implied_mm_ns | reference_experimental_theta_deg |
|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | L | 1.2320989 | 166.198 | 134.89014 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | R | 1.2368844 | 166.198 | 134.36826 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | L | 1.2743246 | 166.198 | 130.42046 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | R | 1.2778527 | 166.198 | 130.06037 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | L | 1.3144043 | 166.198 | 126.44359 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | R | 1.3202916 | 166.198 | 125.87977 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | L | 1.3782471 | 166.198 | 120.5865 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | R | 1.3767515 | 166.198 | 120.7175 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | L | 1.3144043 | 166.198 | 126.44359 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | R | 1.3202916 | 166.198 | 125.87977 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | L | 1.3782471 | 166.198 | 120.5865 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | R | 1.3767515 | 166.198 | 120.7175 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | L | 1.2331643 | 171.893 | 139.3918 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | R | 1.229897 | 171.893 | 139.76211 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | L | 1.2628237 | 171.893 | 136.11797 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | R | 1.2704087 | 171.893 | 135.30527 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | L | 1.3082152 | 171.893 | 131.39505 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | R | 1.3083088 | 171.893 | 131.38565 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | L | 1.3712846 | 171.893 | 125.3518 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | R | 1.3763534 | 171.893 | 124.89017 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | L | 1.3082152 | 171.893 | 131.39505 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | R | 1.3083088 | 171.893 | 131.38565 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | L | 1.3712846 | 171.893 | 125.3518 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | R | 1.3763534 | 171.893 | 124.89017 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | L | 1.2136948 | 189.74206 | 156.33424 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | R | 1.2116225 | 189.74206 | 156.60164 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | L | 1.2612387 | 189.74206 | 150.44104 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | R | 1.257572 | 189.74206 | 150.87968 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | L | 1.3045156 | 189.74206 | 145.45021 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | R | 1.3043049 | 189.74206 | 145.4737 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | L | 1.3675921 | 189.74206 | 138.7417 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | R | 1.3688913 | 189.74206 | 138.61003 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | L | 1.3045156 | 189.74206 | 145.45021 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | R | 1.3043049 | 189.74206 | 145.4737 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | L | 1.3675921 | 189.74206 | 138.7417 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | R | 1.3688913 | 189.74206 | 138.61003 | 35.4 |

## AE3-AE5: configuration references

`fastic_measured` and `penarodriguez_shortened` are identical `(2,3 ns)` models. The two distinct pulse forms are `(2,3 ns)` and `(0.5,2 ns)`. The MUSIC/SAMPIC chain has no model in `pulse_models.py`; comparison with 155 mm/ns is indicative only. Geometry references: simulated bar 140 cm versus test-beam 150 cm; configured muon 1 GeV/c versus test-beam 2.5 GeV/c. Blondel references are 16.1 cm/ns at center and approximately 14 cm/ns for x<20 cm.

## W6 updated reference column

`difference_from_experimental_reference_mm_ns` is simulated minus 155 mm/ns; no compatibility calculation is made.

| material | pulse_model | cfd_fraction | slope_ns_per_mm | slope_error_ns_per_mm | chi2 | ndf | chi2_ndf | v_eff_two_end_mm_ns | v_eff_error_mm_ns | experimental_reference_mm_ns | difference_from_experimental_reference_mm_ns |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | 0.0114179 | 1.1964581e-06 | 6353.0331 | 5 | 1270.6066 | 175.16356 | 0.018355026 | 155 | 20.163556 |
| EJ-200 | fast_contrast | 0.24 | 0.011539494 | 1.2888458e-06 | 9878.2269 | 5 | 1975.6454 | 173.31782 | 0.019357862 | 155 | 18.317819 |
| EJ-200 | fastic_measured | 0.14 | 0.011627942 | 1.1468608e-06 | 14968.329 | 5 | 2993.6658 | 171.99949 | 0.016964263 | 155 | 16.999485 |
| EJ-200 | fastic_measured | 0.24 | 0.011787811 | 1.2625943e-06 | 16834.786 | 5 | 3366.9573 | 169.66678 | 0.018173035 | 155 | 14.66678 |
| EJ-200 | penarodriguez_shortened | 0.14 | 0.011627942 | 1.1468608e-06 | 14968.329 | 5 | 2993.6658 | 171.99949 | 0.016964263 | 155 | 16.999485 |
| EJ-200 | penarodriguez_shortened | 0.24 | 0.011787811 | 1.2625943e-06 | 16834.786 | 5 | 3366.9573 | 169.66678 | 0.018173035 | 155 | 14.66678 |
| EJ-204 | fast_contrast | 0.14 | 0.011367336 | 1.2838607e-06 | 4120.9081 | 5 | 824.18162 | 175.94272 | 0.019871494 | 155 | 20.942718 |
| EJ-204 | fast_contrast | 0.24 | 0.011469029 | 1.3851143e-06 | 6116.801 | 5 | 1223.3602 | 174.38268 | 0.021060192 | 155 | 19.382679 |
| EJ-204 | fastic_measured | 0.14 | 0.011547939 | 1.2417523e-06 | 10701.473 | 5 | 2140.2946 | 173.19108 | 0.018623275 | 155 | 18.191084 |
| EJ-204 | fastic_measured | 0.24 | 0.011679085 | 1.3446566e-06 | 12921.578 | 5 | 2584.3156 | 171.2463 | 0.019716226 | 155 | 16.246298 |
| EJ-204 | penarodriguez_shortened | 0.14 | 0.011547939 | 1.2417523e-06 | 10701.473 | 5 | 2140.2946 | 173.19108 | 0.018623275 | 155 | 18.191084 |
| EJ-204 | penarodriguez_shortened | 0.24 | 0.011679085 | 1.3446566e-06 | 12921.578 | 5 | 2584.3156 | 171.2463 | 0.019716226 | 155 | 16.246298 |
| EJ-230 | fast_contrast | 0.14 | 0.011257884 | 1.3449932e-06 | 5268.1656 | 5 | 1053.6331 | 177.65328 | 0.021224455 | 155 | 22.653283 |
| EJ-230 | fast_contrast | 0.24 | 0.011367194 | 1.473694e-06 | 4406.4808 | 5 | 881.29617 | 175.94491 | 0.022810288 | 155 | 20.944914 |
| EJ-230 | fastic_measured | 0.14 | 0.011444989 | 1.3511242e-06 | 7911.9817 | 5 | 1582.3963 | 174.74896 | 0.020629775 | 155 | 19.748965 |
| EJ-230 | fastic_measured | 0.24 | 0.011584425 | 1.443019e-06 | 9968.2825 | 5 | 1993.6565 | 172.6456 | 0.021505674 | 155 | 17.645597 |
| EJ-230 | penarodriguez_shortened | 0.14 | 0.011444989 | 1.3511242e-06 | 7911.9817 | 5 | 1582.3963 | 174.74896 | 0.020629775 | 155 | 19.748965 |
| EJ-230 | penarodriguez_shortened | 0.24 | 0.011584425 | 1.443019e-06 | 9968.2825 | 5 | 1993.6565 | 172.6456 | 0.021505674 | 155 | 17.645597 |
