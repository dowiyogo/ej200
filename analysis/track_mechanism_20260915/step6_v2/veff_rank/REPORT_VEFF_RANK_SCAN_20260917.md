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
| EJ-200 | fast_contrast | 0.14 | L | 6.1655246 | 28.638838 | 45.75316 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | R | 6.3311465 | 28.669548 | 45.853782 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | L | 6.7401017 | 31.130667 | 49.171849 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | R | 6.869441 | 31.180475 | 49.165436 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | L | 7.6261336 | 33.00193 | 50.815196 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | R | 7.400945 | 33.008087 | 50.833967 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | L | 8.4520833 | 35.595778 | 53.692025 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | R | 8.6022718 | 35.612958 | 53.664692 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | L | 7.6261336 | 33.00193 | 50.815196 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | R | 7.400945 | 33.008087 | 50.833967 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | L | 8.4520833 | 35.595778 | 53.692025 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | R | 8.6022718 | 35.612958 | 53.664692 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | L | 6.168464 | 28.287187 | 45.133463 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | R | 5.7780954 | 28.33502 | 45.145945 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | L | 6.7727635 | 30.647433 | 48.032118 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | R | 6.6586851 | 30.58333 | 48.036209 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | L | 7.5806397 | 32.703391 | 50.128076 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | R | 7.4990337 | 32.746649 | 50.208739 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | L | 8.6141645 | 35.173396 | 52.747926 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | R | 8.4088779 | 35.286613 | 52.960871 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | L | 7.5806397 | 32.703391 | 50.128076 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | R | 7.4990337 | 32.746649 | 50.208739 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | L | 8.6141645 | 35.173396 | 52.747926 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | R | 8.4088779 | 35.286613 | 52.960871 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | L | 5.6319097 | 27.637058 | 43.049257 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | R | 5.6884649 | 27.585686 | 43.039756 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | L | 6.761136 | 30.246904 | 46.285634 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | R | 6.9508349 | 30.221695 | 46.337471 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | L | 7.5462878 | 32.16744 | 48.969842 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | R | 7.5306748 | 32.299058 | 49.034092 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | L | 8.56298 | 35.253136 | 51.868158 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | R | 8.5821235 | 35.331542 | 51.844527 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | L | 7.5462878 | 32.16744 | 48.969842 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | R | 7.5306748 | 32.299058 | 49.034092 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | L | 8.56298 | 35.253136 | 51.868158 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | R | 8.5821235 | 35.331542 | 51.844527 | 35.4 |

### By position

| material | pulse_model | cfd_fraction | x_mm | theta_q05_weighted_mean_deg | theta_median_weighted_mean_deg | theta_q95_weighted_mean_deg | reference_experimental_theta_deg |
|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | -650 | 3.2960041 | 31.948839 | 51.404586 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | -500 | 7.3741456 | 28.541001 | 44.85329 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | -200 | 7.4246129 | 26.786489 | 42.767468 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 0 | 7.0731044 | 25.950607 | 42.222615 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 200 | 7.3367675 | 26.56386 | 42.828881 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 500 | 7.5475428 | 28.763892 | 45.101057 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | 650 | 3.6861718 | 32.02466 | 51.4464 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -650 | 3.7445675 | 34.394049 | 55.301559 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -500 | 8.1949107 | 31.622394 | 49.575455 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | -200 | 7.822527 | 28.752986 | 44.95069 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 0 | 7.5795328 | 28.068636 | 44.26377 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 200 | 7.8753195 | 28.779577 | 44.94475 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 500 | 8.1994752 | 31.748781 | 49.569296 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | 650 | 4.2170667 | 34.722573 | 55.57498 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -650 | 5.0940541 | 36.693422 | 57.414854 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -500 | 9.1410941 | 33.13876 | 50.515651 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | -200 | 8.4796943 | 30.703424 | 46.926438 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 0 | 8.3110738 | 29.993696 | 45.858818 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 200 | 8.043641 | 30.747304 | 46.891949 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 500 | 8.391018 | 33.270311 | 50.732213 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | 650 | 5.1341999 | 36.488139 | 57.432147 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -650 | 6.1401664 | 38.986817 | 60.190681 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -500 | 9.4743323 | 35.430778 | 53.511609 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | -200 | 9.2043256 | 33.405893 | 49.475737 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 0 | 9.2052775 | 33.379835 | 49.420743 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 200 | 9.5638154 | 33.675012 | 49.556003 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 500 | 9.9603916 | 35.525298 | 53.578474 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | 650 | 6.1419338 | 38.82694 | 60.015261 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -650 | 5.0940541 | 36.693422 | 57.414854 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -500 | 9.1410941 | 33.13876 | 50.515651 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | -200 | 8.4796943 | 30.703424 | 46.926438 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 0 | 8.3110738 | 29.993696 | 45.858818 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 200 | 8.043641 | 30.747304 | 46.891949 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 500 | 8.391018 | 33.270311 | 50.732213 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | 650 | 5.1341999 | 36.488139 | 57.432147 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -650 | 6.1401664 | 38.986817 | 60.190681 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -500 | 9.4743323 | 35.430778 | 53.511609 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | -200 | 9.2043256 | 33.405893 | 49.475737 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 0 | 9.2052775 | 33.379835 | 49.420743 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 200 | 9.5638154 | 33.675012 | 49.556003 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 500 | 9.9603916 | 35.525298 | 53.578474 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | 650 | 6.1419338 | 38.82694 | 60.015261 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -650 | 3.7617036 | 31.775041 | 50.997175 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -500 | 7.105274 | 28.110665 | 44.236921 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | -200 | 7.1236047 | 26.122903 | 42.129255 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 0 | 6.9740758 | 26.210632 | 41.280015 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 200 | 7.4316572 | 26.083416 | 42.135107 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 500 | 6.9042885 | 28.150071 | 44.188416 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | 650 | 2.5123543 | 31.724996 | 51.011038 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -650 | 4.2093267 | 33.8809 | 53.524394 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -500 | 7.8769299 | 31.030191 | 48.295831 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | -200 | 7.8986574 | 28.700617 | 44.505462 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 0 | 7.428646 | 27.703323 | 43.59973 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 200 | 7.6853215 | 28.546233 | 44.432015 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 500 | 7.8113551 | 30.770766 | 48.247497 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | 650 | 4.0998336 | 33.675639 | 53.634214 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -650 | 4.9949757 | 36.183056 | 56.815962 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -500 | 8.8110641 | 33.102002 | 49.88691 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | -200 | 8.3601506 | 30.17309 | 46.162248 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 0 | 8.098408 | 29.686523 | 45.193945 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 200 | 8.4084049 | 30.157697 | 46.260272 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 500 | 8.8865177 | 33.360598 | 50.112741 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | 650 | 5.2193361 | 36.412177 | 56.746772 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -650 | 6.4441783 | 38.351393 | 59.170954 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -500 | 9.9815523 | 35.161823 | 52.808541 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | -200 | 9.2296529 | 33.074644 | 48.695342 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 0 | 8.9050837 | 32.766219 | 48.573762 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 200 | 9.5386725 | 33.266534 | 48.576249 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 500 | 9.6646488 | 35.443085 | 52.983739 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | 650 | 5.8168597 | 38.546335 | 59.172203 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -650 | 4.9949757 | 36.183056 | 56.815962 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -500 | 8.8110641 | 33.102002 | 49.88691 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | -200 | 8.3601506 | 30.17309 | 46.162248 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 0 | 8.098408 | 29.686523 | 45.193945 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 200 | 8.4084049 | 30.157697 | 46.260272 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 500 | 8.8865177 | 33.360598 | 50.112741 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | 650 | 5.2193361 | 36.412177 | 56.746772 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -650 | 6.4441783 | 38.351393 | 59.170954 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -500 | 9.9815523 | 35.161823 | 52.808541 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | -200 | 9.2296529 | 33.074644 | 48.695342 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 0 | 8.9050837 | 32.766219 | 48.573762 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 200 | 9.5386725 | 33.266534 | 48.576249 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 500 | 9.6646488 | 35.443085 | 52.983739 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | 650 | 5.8168597 | 38.546335 | 59.172203 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -650 | 2.4344979 | 31.346206 | 49.896899 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -500 | 7.0545169 | 27.365818 | 41.96362 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | -200 | 6.7371604 | 25.045644 | 38.809225 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 0 | 7.1493567 | 26.476246 | 40.178506 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 200 | 6.7985215 | 24.723085 | 38.706948 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 500 | 7.1218616 | 27.385541 | 41.893318 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | 650 | 2.3253963 | 30.937064 | 49.863029 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -650 | 4.0384714 | 33.458497 | 52.25598 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -500 | 7.7935295 | 30.060129 | 45.328469 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | -200 | 7.805566 | 28.305404 | 43.144427 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 0 | 7.7081446 | 27.74787 | 42.4346 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 200 | 7.795871 | 28.190788 | 43.187605 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 500 | 8.3923517 | 30.489768 | 45.714692 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | 650 | 4.4579637 | 33.387643 | 52.115094 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -650 | 4.7042545 | 35.355691 | 54.591089 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -500 | 8.7552396 | 32.238791 | 49.593593 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | -200 | 8.6106369 | 30.264909 | 45.293832 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 0 | 8.1426818 | 29.644605 | 44.635492 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 200 | 8.7075949 | 30.333316 | 45.177624 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 500 | 8.7788803 | 32.313968 | 49.091219 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | 650 | 5.0700812 | 35.481463 | 54.630919 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -650 | 5.7778659 | 38.485329 | 58.207845 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -500 | 9.7360785 | 35.371574 | 51.811406 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | -200 | 9.6073189 | 33.298446 | 47.967766 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 0 | 9.3021335 | 32.607763 | 47.240689 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 200 | 9.5356425 | 33.266743 | 47.835483 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 500 | 9.4654735 | 35.248135 | 51.683351 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | 650 | 6.5833496 | 38.768384 | 58.247859 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -650 | 4.7042545 | 35.355691 | 54.591089 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -500 | 8.7552396 | 32.238791 | 49.593593 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | -200 | 8.6106369 | 30.264909 | 45.293832 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 0 | 8.1426818 | 29.644605 | 44.635492 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 200 | 8.7075949 | 30.333316 | 45.177624 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 500 | 8.7788803 | 32.313968 | 49.091219 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | 650 | 5.0700812 | 35.481463 | 54.630919 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -650 | 5.7778659 | 38.485329 | 58.207845 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -500 | 9.7360785 | 35.371574 | 51.811406 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | -200 | 9.6073189 | 33.298446 | 47.967766 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 0 | 9.3021335 | 32.607763 | 47.240689 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 200 | 9.5356425 | 33.266743 | 47.835483 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 500 | 9.4654735 | 35.248135 | 51.683351 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | 650 | 6.5833496 | 38.768384 | 58.247859 | 35.4 |

## AE2: geometry table

The exact event-level `<path>/<d_axial>` was not persisted; it is `NOT_AVAILABLE`. The proxy `1/cos(theta_median)` and implied velocity are tabulated without compatibility calculations.

| material | pulse_model | cfd_fraction | extreme | path_over_daxial_exact | path_over_daxial_from_theta_median | v_group_MPT_mm_ns | v_implied_mm_ns | reference_experimental_theta_deg |
|---|---|---|---|---|---|---|---|---|
| EJ-200 | fast_contrast | 0.14 | L | NOT_AVAILABLE | 1.1393956 | 166.198 | 145.86505 | 35.4 |
| EJ-200 | fast_contrast | 0.14 | R | NOT_AVAILABLE | 1.1397294 | 166.198 | 145.82234 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | L | NOT_AVAILABLE | 1.1682373 | 166.198 | 142.26391 | 35.4 |
| EJ-200 | fast_contrast | 0.24 | R | NOT_AVAILABLE | 1.1688514 | 166.198 | 142.18916 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | L | NOT_AVAILABLE | 1.1923894 | 166.198 | 139.38232 | 35.4 |
| EJ-200 | fastic_measured | 0.14 | R | NOT_AVAILABLE | 1.1924726 | 166.198 | 139.37259 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | L | NOT_AVAILABLE | 1.229795 | 166.198 | 135.14285 | 35.4 |
| EJ-200 | fastic_measured | 0.24 | R | NOT_AVAILABLE | 1.2300591 | 166.198 | 135.11384 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | L | NOT_AVAILABLE | 1.1923894 | 166.198 | 139.38232 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.14 | R | NOT_AVAILABLE | 1.1924726 | 166.198 | 139.37259 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | L | NOT_AVAILABLE | 1.229795 | 166.198 | 135.14285 | 35.4 |
| EJ-200 | penarodriguez_shortened | 0.24 | R | NOT_AVAILABLE | 1.2300591 | 166.198 | 135.11384 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | L | NOT_AVAILABLE | 1.1356108 | 171.893 | 151.36611 | 35.4 |
| EJ-204 | fast_contrast | 0.14 | R | NOT_AVAILABLE | 1.1361217 | 171.893 | 151.29806 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | L | NOT_AVAILABLE | 1.1623578 | 171.893 | 147.88304 | 35.4 |
| EJ-204 | fast_contrast | 0.24 | R | NOT_AVAILABLE | 1.1615885 | 171.893 | 147.98098 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | L | NOT_AVAILABLE | 1.1883841 | 171.893 | 144.64432 | 35.4 |
| EJ-204 | fastic_measured | 0.14 | R | NOT_AVAILABLE | 1.1889608 | 171.893 | 144.57416 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | L | NOT_AVAILABLE | 1.2233726 | 171.893 | 140.50748 | 35.4 |
| EJ-204 | fastic_measured | 0.24 | R | NOT_AVAILABLE | 1.225081 | 171.893 | 140.31154 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | L | NOT_AVAILABLE | 1.1883841 | 171.893 | 144.64432 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.14 | R | NOT_AVAILABLE | 1.1889608 | 171.893 | 144.57416 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | L | NOT_AVAILABLE | 1.2233726 | 171.893 | 140.50748 | 35.4 |
| EJ-204 | penarodriguez_shortened | 0.24 | R | NOT_AVAILABLE | 1.225081 | 171.893 | 140.31154 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | L | NOT_AVAILABLE | 1.1287908 | 189.74206 | 168.0932 | 35.4 |
| EJ-230 | fast_contrast | 0.14 | R | NOT_AVAILABLE | 1.1282616 | 189.74206 | 168.17205 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | L | NOT_AVAILABLE | 1.1575913 | 189.74206 | 163.91109 | 35.4 |
| EJ-230 | fast_contrast | 0.24 | R | NOT_AVAILABLE | 1.1572945 | 189.74206 | 163.95313 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | L | NOT_AVAILABLE | 1.1813407 | 189.74206 | 160.61587 | 35.4 |
| EJ-230 | fastic_measured | 0.14 | R | NOT_AVAILABLE | 1.1830531 | 189.74206 | 160.38339 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | L | NOT_AVAILABLE | 1.2245748 | 189.74206 | 154.94526 | 35.4 |
| EJ-230 | fastic_measured | 0.24 | R | NOT_AVAILABLE | 1.2257616 | 189.74206 | 154.79525 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | L | NOT_AVAILABLE | 1.1813407 | 189.74206 | 160.61587 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.14 | R | NOT_AVAILABLE | 1.1830531 | 189.74206 | 160.38339 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | L | NOT_AVAILABLE | 1.2245748 | 189.74206 | 154.94526 | 35.4 |
| EJ-230 | penarodriguez_shortened | 0.24 | R | NOT_AVAILABLE | 1.2257616 | 189.74206 | 154.79525 | 35.4 |

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
