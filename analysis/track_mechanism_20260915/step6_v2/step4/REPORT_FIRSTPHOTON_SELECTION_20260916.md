# EXEC_46 Step 4 — first-photon selection and source order statistics

Date: 2026-09-16

## Checkpoint verdict

Step 4 is complete. No simulation was run; the 21 production ROOT files were read-only.
The corrected min-versus-min comparison closes the near-END winner accounting: Cherenkov
pays a 32--51 ps transport handicap but gains 54--67 ps in the creation time attached to
the first detected photon, leaving it earlier by 4--34 ps depending on material.

The pure scintillation emission minimum follows N^-1/2 and rejects N^-1. The smaller
effective populations inferred from the detection-selected photon are therefore a
transport-selection effect, not the emission order statistic itself.

This revised Step 4 record includes D1--D3; Step 5 is reported separately.

## Input and estimator definitions

The derived tree has 420,000 event-END rows. First and random selections reproduce the
Step 2 tree exactly. The random control uses the same splitmix64 seed and the same event
and END face. Every mean-difference interval below uses 2,000 paired bootstrap replicas.
The KS statistic compares the paired marginal samples; the mean difference and bootstrap
retain event pairing.

`min_creation_scint` is the actual minimum creation time among detected scintillation
photons. `first_scint` and `first_cherenkov` are minima in detection time within each
source. No scintillation lifetime law is applied to Cherenkov.

## C0a — Cherenkov multiplicity and the cone edge

The cone angle uses each selected photon created wavelength; transport uses its detected wavelength.
Group velocities are numerical Geant4 mesh values; c/n is retained only as a phase speed.
Median [q05, q95]; wavelengths weighted by the selected photons. Cone: created wavelength; transport: detected wavelength. Spectral correction is (1/vg_det - 1/vg_created), in ps/mm. Constant-index models have no measured dispersion domain.

| material | population | n(created) | theta_C(beta) [deg] | critical(created) [deg] | edge(beta) [deg] | vg(det) [mm/ns] | transport edge [mm/ns] | spectral delay [ps/mm] | group-phase delay [ps/mm] | created clamp / outside fraction | detected clamp / outside fraction |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EJ-200 | all_source_type_2, d=50 mm | 1.592 [1.557, 1.650] | 50.870 [49.826, 52.486] | 38.917 [37.313, 39.953] | 39.130 [37.514, 40.174] | 171.188 [162.267, 192.515] | 132.794 [128.444, 147.098] | 0 / 0 / 0 | 0.527359 / 0 / 0.679842 | 22.8800% / 22.8800% | 22.8800% / 22.8800% |
| EJ-200 | all_source_type_2, d=200 mm | 1.629 [1.583, 1.647] | 51.936 [50.616, 52.409] | 37.859 [37.390, 39.169] | 38.064 [37.591, 39.384] | 164.344 [161.760, 173.067] | 129.392 [128.176, 133.765] | 0 / 0 / 0 | 0.64966 / 0.496954 / 0.688808 | 1.3750% / 1.3750% | 1.3750% / 1.3750% |
| EJ-200 | all_source_type_2, d=500 mm | 1.625 [1.575, 1.645] | 51.809 [50.372, 52.350] | 37.985 [37.448, 39.411] | 38.191 [37.650, 39.628] | 165.078 [162.069, 174.999] | 129.744 [128.319, 134.785] | 0 / 0 / 0 | 0.637925 / 0.46033 / 0.684349 | 2.1000% / 2.1000% | 2.1000% / 2.1000% |
| EJ-200 | all_source_type_2, d=700 mm | 1.621 [1.568, 1.643] | 51.698 [50.166, 52.310] | 38.095 [37.487, 39.616] | 38.302 [37.690, 39.834] | 165.737 [162.278, 176.744] | 130.063 [128.416, 135.722] | 0 / 0 / 0 | 0.62715 / 0.426659 / 0.681275 | 2.7250% / 2.7250% | 2.7250% / 2.7250% |
| EJ-200 | all_source_type_2, d=900 mm | 1.617 [1.564, 1.642] | 51.596 [50.041, 52.268] | 38.197 [37.529, 39.739] | 38.404 [37.732, 39.959] | 166.357 [162.505, 177.844] | 130.366 [128.522, 136.319] | 0 / 0 / 0 | 0.616837 / 0.405202 / 0.677921 | 3.3737% / 3.3737% | 3.3737% / 3.3737% |
| EJ-200 | all_source_type_2, d=1200 mm | 1.612 [1.557, 1.639] | 51.456 [49.829, 52.208] | 38.335 [37.589, 39.950] | 38.544 [37.792, 40.171] | 167.227 [162.828, 179.884] | 130.793 [128.673, 137.454] | 0 / 0 / 0 | 0.602115 / 0.364369 / 0.673106 | 4.9446% / 4.9446% | 4.9446% / 4.9446% |
| EJ-200 | all_source_type_2, d=1350 mm | 1.610 [1.557, 1.639] | 51.397 [49.826, 52.197] | 38.394 [37.600, 39.953] | 38.603 [37.803, 40.174] | 167.603 [162.890, 181.133] | 130.980 [128.703, 138.402] | 0 / 0 / 0 | 0.595655 / 0.326401 / 0.672167 | 5.3114% / 5.3114% | 5.3114% / 5.3114% |
| EJ-204 | all_source_type_2, d=50 mm | 1.598 [1.585, 1.633] | 51.046 [50.658, 52.038] | 38.743 [37.758, 39.128] | 38.954 [37.962, 39.342] | 179.358 [169.405, 189.185] | 139.480 [133.441, 146.310] | 0 / 0 / 0 | 0.245344 / 0 / 0.463436 | 32.9350% / 32.9350% | 32.9350% / 32.9350% |
| EJ-204 | all_source_type_2, d=200 mm | 1.633 [1.591, 1.633] | 52.038 [50.852, 52.038] | 37.758 [37.758, 38.935] | 37.962 [37.962, 39.148] | 170.949 [168.590, 182.145] | 134.662 [132.893, 141.263] | 0 / 0 / 0 | 0.410222 / 0.182046 / 0.486101 | 59.2850% / 59.2850% | 59.2850% / 59.2850% |
| EJ-204 | all_source_type_2, d=500 mm | 1.633 [1.597, 1.633] | 52.038 [51.023, 52.038] | 37.758 [37.758, 38.765] | 37.962 [37.962, 38.977] | 170.686 [168.564, 179.634] | 134.432 [132.876, 139.647] | 0 / 0 / 0 | 0.420961 / 0.239491 / 0.486747 | 59.4150% / 59.4150% | 59.4150% / 59.4150% |
| EJ-204 | all_source_type_2, d=700 mm | 1.633 [1.592, 1.633] | 52.038 [50.887, 52.038] | 37.758 [37.758, 38.900] | 37.962 [37.962, 39.113] | 170.948 [168.598, 181.617] | 134.593 [132.902, 140.916] | 0 / 0 / 0 | 0.414896 / 0.194332 / 0.485687 | 54.7027% / 54.7027% | 54.7027% / 54.7027% |
| EJ-204 | all_source_type_2, d=900 mm | 1.633 [1.590, 1.633] | 52.031 [50.823, 52.038] | 37.765 [37.758, 38.964] | 37.969 [37.962, 39.177] | 171.174 [168.590, 182.611] | 134.764 [132.897, 141.559] | 0 / 0 / 0 | 0.407805 / 0.171627 / 0.485658 | 51.6261% / 51.6261% | 51.6261% / 51.6261% |
| EJ-204 | all_source_type_2, d=1200 mm | 1.627 [1.587, 1.633] | 51.862 [50.731, 52.038] | 37.933 [37.758, 39.055] | 38.138 [37.962, 39.269] | 171.585 [168.644, 184.121] | 135.029 [132.932, 142.544] | 0 / 0 / 0 | 0.398406 / 0.137066 / 0.484306 | 43.6073% / 43.6073% | 43.6073% / 43.6073% |
| EJ-204 | all_source_type_2, d=1350 mm | 1.626 [1.587, 1.633] | 51.854 [50.717, 52.038] | 37.940 [37.758, 39.068] | 38.146 [37.962, 39.283] | 171.676 [168.655, 184.359] | 135.098 [132.936, 142.700] | 0 / 0 / 0 | 0.395059 / 0.131606 / 0.484198 | 44.0466% / 44.0466% | 44.0466% / 44.0466% |
| EJ-230 | all_source_type_2, d=50 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 37.0000% | 0.0000% / 37.0000% |
| EJ-230 | all_source_type_2, d=200 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.4400% | 0.0000% / 36.4400% |
| EJ-230 | all_source_type_2, d=500 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.9300% | 0.0000% / 36.9300% |
| EJ-230 | all_source_type_2, d=700 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 37.1574% | 0.0000% / 37.1574% |
| EJ-230 | all_source_type_2, d=900 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.4880% | 0.0000% / 36.4880% |
| EJ-230 | all_source_type_2, d=1200 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.9226% | 0.0000% / 36.9226% |
| EJ-230 | all_source_type_2, d=1350 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.6288% | 0.0000% / 36.6288% |
| EJ-200 | primary_like_proxy, d=50 mm | 1.586 [1.557, 1.640] | 50.708 [49.826, 52.226] | 39.077 [37.571, 39.953] | 39.292 [37.774, 40.174] | 172.365 [163.024, 192.515] | 133.399 [128.787, 147.098] | 0 / 0 / 0 | 0.509383 / 0 / 0.669046 | 19.6069% / 19.6069% | 19.6069% / 19.6069% |
| EJ-200 | primary_like_proxy, d=200 mm | 1.613 [1.564, 1.646] | 51.477 [50.030, 52.381] | 38.314 [37.417, 39.750] | 38.523 [37.619, 39.970] | 167.093 [161.907, 177.943] | 130.728 [128.244, 136.373] | 0 / 0 / 0 | 0.604358 / 0.403261 / 0.686696 | 3.4711% / 3.4711% | 3.4711% / 3.4711% |
| EJ-200 | primary_like_proxy, d=500 mm | 1.611 [1.557, 1.642] | 51.424 [49.826, 52.280] | 38.367 [37.518, 39.953] | 38.576 [37.720, 40.174] | 167.431 [162.442, 184.908] | 130.894 [128.493, 141.285] | 0 / 0 / 0 | 0.598609 / 0.213718 / 0.67886 | 6.5670% / 6.5670% | 6.5670% / 6.5670% |
| EJ-200 | primary_like_proxy, d=700 mm | 1.608 [1.557, 1.641] | 51.345 [49.826, 52.254] | 38.446 [37.543, 39.953] | 38.655 [37.746, 40.174] | 167.937 [162.578, 186.639] | 131.146 [128.556, 142.609] | 0 / 0 / 0 | 0.589878 / 0.163536 / 0.676844 | 7.3511% / 7.3511% | 7.3511% / 7.3511% |
| EJ-200 | primary_like_proxy, d=900 mm | 1.607 [1.557, 1.639] | 51.303 [49.826, 52.195] | 38.487 [37.601, 39.953] | 38.697 [37.805, 40.174] | 168.210 [162.897, 185.811] | 131.281 [128.706, 141.975] | 0 / 0 / 0 | 0.585136 / 0.187435 / 0.672062 | 6.9543% / 6.9543% | 6.9543% / 6.9543% |
| EJ-200 | primary_like_proxy, d=1200 mm | 1.604 [1.557, 1.638] | 51.226 [49.826, 52.159] | 38.564 [37.638, 39.953] | 38.774 [37.841, 40.174] | 168.718 [163.097, 187.357] | 131.536 [128.800, 143.157] | 0 / 0 / 0 | 0.576224 / 0.14301 / 0.669036 | 8.4701% / 8.4701% | 8.4701% / 8.4701% |
| EJ-200 | primary_like_proxy, d=1350 mm | 1.605 [1.557, 1.637] | 51.250 [49.826, 52.154] | 38.540 [37.643, 39.953] | 38.750 [37.846, 40.174] | 168.556 [163.124, 185.464] | 131.454 [128.813, 141.711] | 0 / 0 / 0 | 0.579081 / 0.197479 / 0.668631 | 7.0551% / 7.0551% | 7.0551% / 7.0551% |
| EJ-204 | primary_like_proxy, d=50 mm | 1.597 [1.585, 1.633] | 51.033 [50.658, 52.038] | 38.755 [37.758, 39.128] | 38.967 [37.962, 39.342] | 179.496 [169.461, 189.185] | 139.560 [133.508, 146.310] | 0 / 0 / 0 | 0.24261 / 0 / 0.460906 | 28.2959% / 28.2959% | 28.2959% / 28.2959% |
| EJ-204 | primary_like_proxy, d=200 mm | 1.627 [1.585, 1.633] | 51.876 [50.673, 52.038] | 37.918 [37.758, 39.112] | 38.124 [37.962, 39.327] | 171.970 [168.670, 185.150] | 135.434 [132.950, 143.222] | 0 / 0 / 0 | 0.379189 / 0.113434 / 0.483595 | 47.5986% / 47.5986% | 47.5986% / 47.5986% |
| EJ-204 | primary_like_proxy, d=500 mm | 1.627 [1.587, 1.633] | 51.875 [50.717, 52.038] | 37.919 [37.758, 39.069] | 38.125 [37.962, 39.283] | 171.476 [168.657, 184.363] | 134.966 [132.939, 142.702] | 0 / 0 / 0 | 0.400131 / 0.131527 / 0.484146 | 43.7465% / 43.7465% | 43.7465% / 43.7465% |
| EJ-204 | primary_like_proxy, d=700 mm | 1.622 [1.585, 1.633] | 51.718 [50.671, 52.038] | 38.075 [37.758, 39.115] | 38.282 [37.962, 39.329] | 172.070 [168.718, 185.198] | 135.320 [132.987, 143.253] | 0 / 0 / 0 | 0.387607 / 0.112329 / 0.481985 | 39.2767% / 39.2767% | 39.2767% / 39.2767% |
| EJ-204 | primary_like_proxy, d=900 mm | 1.624 [1.586, 1.633] | 51.798 [50.695, 52.038] | 37.996 [37.758, 39.091] | 38.202 [37.962, 39.305] | 171.933 [168.633, 184.759] | 135.276 [132.926, 142.963] | 0 / 0 / 0 | 0.389072 / 0.122421 / 0.484586 | 43.2807% / 43.2807% | 43.2807% / 43.2807% |
| EJ-204 | primary_like_proxy, d=1200 mm | 1.619 [1.585, 1.633] | 51.658 [50.665, 52.038] | 38.135 [37.758, 39.120] | 38.342 [37.962, 39.335] | 172.360 [168.748, 185.300] | 135.492 [133.002, 143.321] | 0 / 0 / 0 | 0.380736 / 0.109989 / 0.481576 | 37.6108% / 37.6108% | 37.6108% / 37.6108% |
| EJ-204 | primary_like_proxy, d=1350 mm | 1.624 [1.586, 1.633] | 51.793 [50.693, 52.038] | 38.001 [37.758, 39.092] | 38.207 [37.962, 39.307] | 171.877 [168.674, 184.787] | 135.280 [132.959, 142.982] | 0 / 0 / 0 | 0.387857 / 0.12177 / 0.483384 | 42.9294% / 42.9294% | 42.9294% / 42.9294% |
| EJ-230 | primary_like_proxy, d=50 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.7418% | 0.0000% / 36.7418% |
| EJ-230 | primary_like_proxy, d=200 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.2711% | 0.0000% / 36.2711% |
| EJ-230 | primary_like_proxy, d=500 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 37.0175% | 0.0000% / 37.0175% |
| EJ-230 | primary_like_proxy, d=700 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.2107% | 0.0000% / 36.2107% |
| EJ-230 | primary_like_proxy, d=900 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.3526% | 0.0000% / 36.3526% |
| EJ-230 | primary_like_proxy, d=1200 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 37.0884% | 0.0000% / 37.0884% |
| EJ-230 | primary_like_proxy, d=1350 mm | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.5827% | 0.0000% / 36.5827% |

| material | population | phase speed(created) [mm/ns] | group index(det) | phase edge(beta=1) [mm/ns] |
|---|---|---|---|---|
| EJ-200 | all_source_type_2, d=50 mm | 188.329 [181.727, 192.515] | 1.751 [1.557, 1.848] | 146.530 [144.533, 147.576] |
| EJ-200 | all_source_type_2, d=200 mm | 183.989 [182.043, 189.352] | 1.824 [1.732, 1.853] | 145.263 [144.638, 146.802] |
| EJ-200 | all_source_type_2, d=500 mm | 184.508 [182.286, 190.331] | 1.816 [1.713, 1.850] | 145.424 [144.718, 147.052] |
| EJ-200 | all_source_type_2, d=700 mm | 184.962 [182.449, 191.160] | 1.809 [1.696, 1.847] | 145.563 [144.771, 147.257] |
| EJ-200 | all_source_type_2, d=900 mm | 185.380 [182.624, 191.656] | 1.802 [1.686, 1.845] | 145.689 [144.828, 147.376] |
| EJ-200 | all_source_type_2, d=1200 mm | 185.950 [182.871, 192.502] | 1.793 [1.667, 1.841] | 145.858 [144.908, 147.573] |
| EJ-200 | all_source_type_2, d=1350 mm | 186.191 [182.917, 192.515] | 1.789 [1.655, 1.840] | 145.929 [144.923, 147.576] |
| EJ-204 | all_source_type_2, d=50 mm | 187.618 [183.570, 189.185] | 1.671 [1.585, 1.770] | 146.335 [145.132, 146.758] |
| EJ-204 | all_source_type_2, d=200 mm | 183.570 [183.570, 188.400] | 1.754 [1.646, 1.778] | 145.132 [145.132, 146.549] |
| EJ-204 | all_source_type_2, d=500 mm | 183.570 [183.570, 187.709] | 1.756 [1.669, 1.779] | 145.132 [145.132, 146.360] |
| EJ-204 | all_source_type_2, d=700 mm | 183.570 [183.570, 188.260] | 1.754 [1.651, 1.778] | 145.132 [145.132, 146.511] |
| EJ-204 | all_source_type_2, d=900 mm | 183.598 [183.570, 188.520] | 1.751 [1.642, 1.778] | 145.141 [145.132, 146.581] |
| EJ-204 | all_source_type_2, d=1200 mm | 184.292 [183.570, 188.888] | 1.747 [1.628, 1.778] | 145.358 [145.132, 146.680] |
| EJ-204 | all_source_type_2, d=1350 mm | 184.323 [183.570, 188.943] | 1.746 [1.626, 1.778] | 145.367 [145.132, 146.695] |
| EJ-230 | all_source_type_2, d=50 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=200 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=500 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=700 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=900 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=1200 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | all_source_type_2, d=1350 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-200 | primary_like_proxy, d=50 mm | 188.980 [182.796, 192.515] | 1.739 [1.557, 1.839] | 146.704 [144.884, 147.576] |
| EJ-200 | primary_like_proxy, d=200 mm | 185.864 [182.157, 191.699] | 1.794 [1.685, 1.852] | 145.833 [144.676, 147.386] |
| EJ-200 | primary_like_proxy, d=500 mm | 186.082 [182.575, 192.515] | 1.791 [1.621, 1.846] | 145.897 [144.813, 147.576] |
| EJ-200 | primary_like_proxy, d=700 mm | 186.403 [182.680, 192.515] | 1.785 [1.606, 1.844] | 145.990 [144.846, 147.576] |
| EJ-200 | primary_like_proxy, d=900 mm | 186.573 [182.923, 192.515] | 1.782 [1.613, 1.840] | 146.040 [144.925, 147.576] |
| EJ-200 | primary_like_proxy, d=1200 mm | 186.887 [183.074, 192.515] | 1.777 [1.600, 1.838] | 146.129 [144.974, 147.576] |
| EJ-200 | primary_like_proxy, d=1350 mm | 186.788 [183.094, 192.515] | 1.779 [1.616, 1.838] | 146.101 [144.980, 147.576] |
| EJ-204 | primary_like_proxy, d=50 mm | 187.669 [183.570, 189.185] | 1.670 [1.585, 1.769] | 146.349 [145.132, 146.758] |
| EJ-204 | primary_like_proxy, d=200 mm | 184.233 [183.570, 189.122] | 1.743 [1.619, 1.777] | 145.339 [145.132, 146.742] |
| EJ-204 | primary_like_proxy, d=500 mm | 184.237 [183.570, 188.944] | 1.748 [1.626, 1.778] | 145.341 [145.132, 146.695] |
| EJ-204 | primary_like_proxy, d=700 mm | 184.881 [183.570, 189.132] | 1.742 [1.619, 1.777] | 145.538 [145.132, 146.744] |
| EJ-204 | primary_like_proxy, d=900 mm | 184.554 [183.570, 189.035] | 1.744 [1.623, 1.778] | 145.439 [145.132, 146.719] |
| EJ-204 | primary_like_proxy, d=1200 mm | 185.127 [183.570, 189.155] | 1.739 [1.618, 1.777] | 145.613 [145.132, 146.750] |
| EJ-204 | primary_like_proxy, d=1350 mm | 184.574 [183.570, 189.041] | 1.744 [1.622, 1.777] | 145.445 [145.132, 146.720] |
| EJ-230 | primary_like_proxy, d=50 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=200 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=500 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=700 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=900 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=1200 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | primary_like_proxy, d=1350 mm | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |

The beta=1 identity arccos(sin(theta_C(lambda))) = arcsin(1/n(lambda)) is checked at every table node and evaluated photon wavelength; dispersion replicates the identity wavelength by wavelength. Arrival-time edge speeds use vg(lambda_det)*sin(theta_C(lambda_created)); the phase-edge formula c*sqrt(1-1/n^2)/n is retained separately. Multi-medium/reflected paths are not reconstructed by this homogeneous-scintillator reference.
Axial velocity is the measured axial displacement divided by pure propagation time.

| material | d [mm] | proxy N_C low/high quintile | v low/high [mm/ns] | high-edge [mm/ns] |
|---|---:|---:|---:|---:|
| EJ-200 | 50 | 5.05 / 26.49 | 131.715 / 143.597 | +8.661 |
| EJ-200 | 200 | 2.66 / 15.20 | 132.744 / 146.752 | +16.149 |
| EJ-200 | 500 | 1.30 / 8.33 | 127.880 / 146.862 | +15.773 |
| EJ-200 | 700 | 1.00 / 6.12 | 124.992 / 146.580 | +15.244 |
| EJ-200 | 900 | 1.00 / 4.80 | 126.496 / 146.221 | +14.720 |
| EJ-200 | 1200 | 1.00 / 3.78 | 127.985 / 144.950 | +12.863 |
| EJ-200 | 1350 | 1.00 / 3.57 | 131.699 / 144.954 | +12.964 |
| EJ-204 | 50 | 7.40 / 37.42 | 135.555 / 145.218 | +5.424 |
| EJ-204 | 200 | 3.75 / 20.32 | 139.049 / 147.635 | +11.715 |
| EJ-204 | 500 | 1.51 / 9.92 | 133.423 / 148.053 | +12.838 |
| EJ-204 | 700 | 1.00 / 6.58 | 128.461 / 147.604 | +11.824 |
| EJ-204 | 900 | 1.00 / 5.02 | 131.154 / 147.617 | +12.074 |
| EJ-204 | 1200 | 1.00 / 3.66 | 130.845 / 146.076 | +9.861 |
| EJ-204 | 1350 | 1.00 / 3.52 | 136.651 / 146.344 | +10.526 |
| EJ-230 | 50 | 7.84 / 36.98 | 140.358 / 145.722 | -0.728 |
| EJ-230 | 200 | 3.53 / 18.63 | 137.931 / 146.024 | -0.426 |
| EJ-230 | 500 | 1.42 / 8.77 | 133.678 / 146.012 | -0.437 |
| EJ-230 | 700 | 1.00 / 5.29 | 129.455 / 144.706 | -1.743 |
| EJ-230 | 900 | 1.00 / 4.02 | 132.398 / 145.231 | -1.219 |
| EJ-230 | 1200 | 1.00 / 2.79 | 129.666 / 141.875 | -4.575 |
| EJ-230 | 1350 | 1.00 / 2.68 | 135.103 / 143.526 | -2.924 |

For the primary-like proxy, increasing N_C moves the minimum toward the
finite-beta cone edge in the low-occupancy long-distance cells and leaves it on
the edge where occupancy is already high. The undifferentiated source-type-2
sample departs above the edge at high N_C because it includes secondary-particle
Cherenkov; both results are retained in `cherenkov_nc_velocity.csv`.

## C0b — angular window versus distance

The timing formula is tested without a fitted speed:

```text
Delta t(alpha,d,lambda) = d/vg(lambda_det) * [1/cos(alpha) - 1/cos(alpha_edge(lambda_created,beta))]
```

The stored final angle is folded to the axial magnitude because the penalty uses
|cos(alpha)|. The upper window is a selection in arrival time, not a second cone
boundary. At q95(alpha), the prediction averages the individual wavelength/distance
counterfactuals; its median/q05/q95 are retained in the CSV. This is not the quantile
of each photon actual penalty.

| material | d [mm] | q95 alpha proxy [deg] | predicted penalty mean [ps] | predicted penalty median [q05,q95] [ns] | empirical q95 excess [ps] |
|---|---:|---:|---:|---:|---:|
| EJ-200 | 50 | 46.006 | 42.69 | 0.04244 [0.03367, 0.05322] | 16.44 |
| EJ-200 | 200 | 49.271 | 301.35 | 0.30387 [0.25557, 0.33286] | 123.38 |
| EJ-200 | 500 | 52.254 | 1040.41 | 1.05750 [0.87737, 1.13581] | 648.64 |
| EJ-200 | 700 | 54.135 | 1751.92 | 1.77568 [1.49195, 1.90256] | 1780.37 |
| EJ-200 | 900 | 54.065 | 2229.52 | 2.26038 [1.91325, 2.42054] | 4740.43 |
| EJ-200 | 1200 | 53.464 | 2782.13 | 2.82286 [2.37507, 3.04077] | 3176.40 |
| EJ-200 | 1350 | 56.469 | 4176.94 | 4.22808 [3.64967, 4.50017] | 3113.28 |
| EJ-204 | 50 | 44.125 | 30.38 | 0.02954 [0.02621, 0.03643] | 13.88 |
| EJ-204 | 200 | 44.167 | 135.30 | 0.14445 [0.10921, 0.14864] | 27.33 |
| EJ-204 | 500 | 48.821 | 702.77 | 0.72571 [0.61473, 0.74140] | 400.21 |
| EJ-204 | 700 | 51.131 | 1275.40 | 1.30455 [1.13594, 1.34766] | 962.80 |
| EJ-204 | 900 | 51.680 | 1745.56 | 1.79421 [1.55994, 1.83633] | 1856.96 |
| EJ-204 | 1200 | 52.135 | 2425.45 | 2.47002 [2.17676, 2.56352] | 2844.63 |
| EJ-204 | 1350 | 50.905 | 2406.90 | 2.47704 [2.14290, 2.53839] | 1869.62 |
| EJ-230 | 50 | 41.849 | 12.24 | 0.01222 [0.01222, 0.01232] | 12.28 |
| EJ-230 | 200 | 45.059 | 126.25 | 0.12624 [0.12624, 0.12624] | 127.36 |
| EJ-230 | 500 | 49.604 | 651.42 | 0.65139 [0.65139, 0.65175] | 657.56 |
| EJ-230 | 700 | 53.979 | 1492.49 | 1.49241 [1.49240, 1.49316] | 1562.89 |
| EJ-230 | 900 | 53.980 | 1919.49 | 1.91939 [1.91939, 1.92018] | 2274.38 |
| EJ-230 | 1200 | 54.262 | 2632.99 | 2.63288 [2.63288, 2.63371] | 4371.43 |
| EJ-230 | 1350 | 54.557 | 3050.14 | 3.05007 [3.05007, 3.05076] | 3642.78 |

This distance-only comparison is not a valid rejection of angular selection:
the low-N_C quintile is N_C=1 at d >= 700 mm, where no order statistic exists.
D1 below conditions on distance and uses N_C as the control variable.

### D1 — angular q95 at fixed distance versus N_C

| material | d [mm] | low/high mean N_C | low/high q95 [deg] | change [deg] |
|---|---:|---:|---:|---:|
| EJ-200 | 50 | 5.05/26.49 | 61.473/41.359 | -20.113 |
| EJ-200 | 200 | 2.66/15.20 | 63.092/40.578 | -22.515 |
| EJ-200 | 500 | 1.30/8.33 | 64.237/41.149 | -23.088 |
| EJ-204 | 50 | 7.40/37.42 | 55.300/40.694 | -14.606 |
| EJ-204 | 200 | 3.75/20.32 | 55.192/39.623 | -15.568 |
| EJ-204 | 500 | 1.51/9.92 | 61.969/39.899 | -22.070 |
| EJ-230 | 50 | 7.84/36.98 | 52.478/40.643 | -11.836 |
| EJ-230 | 200 | 3.53/18.63 | 55.222/40.199 | -15.023 |
| EJ-230 | 500 | 1.42/8.77 | 60.551/40.370 | -20.182 |

q95 narrows from the lowest to highest N_C quintile in 9/9 fixed-distance groups. The distance-only test is
therefore reclassified as badly conditioned rather than a refutation of the
cone-edge order-statistics mechanism.

## C0c — corrected d=50 mm handicap

| material | mirror | N_S | pure-min creation C-S [ps] | selected creation C-S [ps] | transport C-S [ps] | total C-S [ps] | Cher wins | reorder gap [ps] |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | x=-650/left | 2312.0 | -36.88 | -66.59 | 32.44 | -34.15 | 82.02% | 29.71 |
| EJ-200 | x=+650/right | 2304.0 | -36.73 | -66.36 | 31.96 | -34.41 | 82.26% | 29.64 |
| EJ-204 | x=-650/left | 2362.2 | -30.20 | -57.01 | 36.26 | -20.75 | 72.88% | 26.81 |
| EJ-204 | x=+650/right | 2364.2 | -30.45 | -57.11 | 36.75 | -20.36 | 72.61% | 26.66 |
| EJ-230 | x=-650/left | 2089.8 | -28.95 | -54.44 | 50.89 | -3.55 | 51.77% | 25.49 |
| EJ-230 | x=+650/right | 2092.3 | -29.50 | -54.99 | 51.02 | -3.97 | 52.33% | 25.49 |

The winner fractions reproduce 74.10/73.56%, 62.90/63.42%, and
51.77/52.33% for the two mirrors. Comparing true source-specific creation minima
gives a 29--40 ps Cherenkov advantage, directly matching the scale of the roughly
38 ps napkin discrepancy. Selection by detection time adds another 25--30 ps to the scintillation delay. After the 45--51 ps Cherenkov transport handicap, the measured source-specific
detection minima differ by -26 to -4 ps; this reproduces the winner fractions
without comparing a mean from one source with a minimum from the other.

## C1 and C4 — scintillation order statistic and N_eff

| material | model | exponent | chi2/ndf | fitted effective fraction (sqrt model) |
|---|---|---:|---:|---:|
| EJ-200 | N^-1 | 1.0 | 1841.14/12 = 153.43 | -- |
| EJ-200 | N^-1/2 | 0.5 | 5.28/12 = 0.44 | 0.8015 +/- 0.0059 |
| EJ-200 | free | 0.5035 [0.4925, 0.5145] | 5.19/11 = 0.47 | -- |
| EJ-204 | N^-1 | 1.0 | 3723.98/12 = 310.33 | -- |
| EJ-204 | N^-1/2 | 0.5 | 18.33/12 = 1.53 | 0.7684 +/- 0.0049 |
| EJ-204 | free | 0.5135 [0.5060, 0.5210] | 15.32/11 = 1.39 | -- |
| EJ-230 | N^-1 | 1.0 | 5607.00/12 = 467.25 | -- |
| EJ-230 | N^-1/2 | 0.5 | 39.20/12 = 3.27 | 0.7361 +/- 0.0043 |
| EJ-230 | free | 0.5115 [0.5055, 0.5175] | 35.69/11 = 3.24 | -- |

The N^-1/2 model gives chi2/ndf 0.44--3.27; N^-1 is rejected with chi2/ndf 153.43--467.25. The fitted effective fractions 0.801, 0.768 and 0.736
differ by at least 4.3 formal standard deviations, so N_eff/N_scint is not material-independent within this model.
A universal geometric factor is therefore rejected. The preregistered fractions
needed to force the detection-selected handicap, 0.20/0.13/0.09, are even more
strongly material dependent. The low-boundary fractions span 0.019--0.269, 0.021--0.275, and 0.021--0.280; they are strongly distance dependent and do not reproduce the fitted effective fractions. `scintillation_order_points.csv` gives fitted N_eff, the low-boundary count, and their ratio at each of the fourteen mirror-resolved points per material.

### D2 — correction for primary-muon transit

The gun points along -z and enters the 10-mm bar at z=+5 mm. The corrected
creation coordinate is `t_creation-(5 mm-z_creation)/c`; a common upstream flight
offset is absorbed by the fitted intercept.

| material | raw f_eff | corrected f_eff | corrected exponent | corrected chi2/ndf |
|---|---:|---:|---:|---:|
| EJ-200 | 0.8015 +/- 0.0059 | 0.8396 +/- 0.0061 | 0.5235 | 0.74 |
| EJ-204 | 0.7684 +/- 0.0049 | 0.8156 +/- 0.0051 | 0.5415 | 3.52 |
| EJ-230 | 0.7361 +/- 0.0043 | 0.7934 +/- 0.0047 | 0.5435 | 6.54 |

The raw effective-fraction ordering follows the d=50-mm asymptotic minimum
delays 30.2 > 24.7 > 20.6 ps: the faster material suffers the larger fractional
dilution from the same 33.4-ps traversal. Subtracting the transit moves every
fraction toward one but recovers only 19.2--21.7% of the original
deficit. N_eff does not reach N_scint; primary transit is a real contribution but
does not explain the remaining 15--21% deficit. The corrected free exponents also
remain above 0.5, so the exact i.i.d. common-origin law is not restored.

## C2 — fixed-point Cherenkov/width test

No interpolation is used.

| d [mm] | Pearson r | Spearman rho | role |
|---:|---:|---:|---|
| 50 | -0.9584 | -0.8857 | Cherenkov-dominated operating regime |
| 200 | +0.9307 | +0.9429 | spectator regime |

The six measured values at each fixed point are:

| d [mm] | material | mirror | fCher(first) | sigma_IQR(near)/sigma_IQR(T0) |
|---:|---|---|---:|---:|
| 50 | EJ-200 | x=-650/left | 82.02% | 0.1299 |
| 50 | EJ-200 | x=+650/right | 82.26% | 0.1389 |
| 50 | EJ-204 | x=-650/left | 72.88% | 0.1884 |
| 50 | EJ-204 | x=+650/right | 72.61% | 0.1841 |
| 50 | EJ-230 | x=-650/left | 51.77% | 0.2313 |
| 50 | EJ-230 | x=+650/right | 52.33% | 0.2249 |
| 200 | EJ-200 | x=-500/left | 9.65% | 0.7515 |
| 200 | EJ-200 | x=+500/right | 9.68% | 0.7395 |
| 200 | EJ-204 | x=-500/left | 8.86% | 0.6727 |
| 200 | EJ-204 | x=+500/right | 9.09% | 0.6737 |
| 200 | EJ-230 | x=-500/left | 7.18% | 0.6328 |
| 200 | EJ-230 | x=+500/right | 7.02% | 0.6090 |

At 50 mm, more Cherenkov corresponds to a smaller near-END/T0 robust-width
ratio. At 200 mm, where the Cherenkov fraction is about 7--9%, the correlation
reverses. This sign change bounds the operating range of the mechanism; it is not
evidence that either measured fixed-point correlation is internally inconsistent.

## C3 — physical source mixture at the near END

| material | mirror | fCher | separation [ps] | sigma_G [ps] | chi2/ndf | sigma_mix [ps] | qwidth [ps] | angle-time r |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | x=-650/left | 82.02% | -10.93 | 2.34 | 205.2 | 16.41 | 10.46 | +0.950 |
| EJ-200 | x=+650/right | 82.26% | -11.10 | 2.32 | 210.1 | 16.43 | 10.27 | +0.919 |
| EJ-204 | x=-650/left | 72.88% | -9.73 | 10.68 | 460.2 | 17.44 | 13.84 | +0.968 |
| EJ-204 | x=+650/right | 72.61% | -10.62 | 11.05 | 503.7 | 17.47 | 13.67 | +0.946 |
| EJ-230 | x=-650/left | 51.77% | -11.26 | 12.03 | 358.9 | 20.00 | 18.91 | +0.986 |
| EJ-230 | x=+650/right | 52.33% | -10.92 | 13.28 | 356.7 | 20.18 | 19.18 | +0.965 |

The source labels partition every event, and the source-mixture variance
reconstructs the total RMS exactly. The rejected Gaussian combines two populations
with different means and shapes; the Cherenkov component also carries the axial
caustic. The replacement width for these six cells is therefore `sigma_mixture`,
defined as sqrt(within-source variance + between-source variance), with qwidth as
the robust cross-check. It is an intrinsic, zero-jitter limit of the first-photon
estimator in this geometry, not detector performance. A real detector includes
S13360 SPTR and a readout target below 50 ps, both well above this approximately
19-ps optical limit; instrumentation therefore sets the attainable resolution.
The approximately 12-ps source separation is smaller than sigma_mixture, so
Cherenkov and scintillation cannot be tagged event by event from this timestamp.
Conversely, the angle-time correlations of +0.962 to +0.986 are a direct signature
of the Cherenkov cone. The Gaussian-fit failure is physical, not numerical.

## Original first-versus-random selection result

All requested variables, KS statistics, paired mean shifts, and bootstrap intervals
are in `paired_selection_statistics.csv`. `paired_deltas.pdf` shows the distance
dependence and `paired_distributions.pdf` the paired marginals. Cherenkov enrichment
is strongest at 50 mm and is reported independently of the source-order tests.

| material | d [mm] | fCher(first) | fCher(random) | enrichment | first-random detection [ns] |
|---|---:|---:|---:|---:|---:|
| EJ-200 | 50 | 82.140% | 2.185% | 37.593 | -3.0732 |
| EJ-200 | 200 | 9.665% | 1.925% | 5.021 | -3.4020 |
| EJ-200 | 500 | 5.325% | 2.130% | 2.500 | -3.7861 |
| EJ-200 | 700 | 3.975% | 1.910% | 2.081 | -4.0080 |
| EJ-200 | 900 | 3.655% | 1.805% | 2.025 | -4.2488 |
| EJ-200 | 1200 | 2.790% | 1.830% | 1.525 | -4.5561 |
| EJ-200 | 1350 | 2.515% | 1.885% | 1.334 | -4.5168 |
| EJ-204 | 50 | 72.745% | 2.865% | 25.391 | -2.5680 |
| EJ-204 | 200 | 8.975% | 3.085% | 2.909 | -2.8334 |
| EJ-204 | 500 | 5.605% | 2.780% | 2.016 | -3.1651 |
| EJ-204 | 700 | 4.575% | 2.605% | 1.756 | -3.3992 |
| EJ-204 | 900 | 3.950% | 2.550% | 1.549 | -3.5690 |
| EJ-204 | 1200 | 3.035% | 2.720% | 1.116 | -3.7893 |
| EJ-204 | 1350 | 2.700% | 2.935% | 0.920 | -3.7495 |
| EJ-230 | 50 | 52.050% | 3.175% | 16.394 | -2.1412 |
| EJ-230 | 200 | 7.100% | 3.385% | 2.097 | -2.3912 |
| EJ-230 | 500 | 4.515% | 2.885% | 1.565 | -2.6905 |
| EJ-230 | 700 | 3.760% | 2.905% | 1.294 | -2.8053 |
| EJ-230 | 900 | 3.110% | 2.875% | 1.082 | -2.9486 |
| EJ-230 | 1200 | 1.985% | 2.615% | 0.759 | -3.1595 |
| EJ-230 | 1350 | 1.795% | 2.845% | 0.631 | -3.1833 |

## Step 3 corrections

The Step 3 report now states that the finite-beta edge falls inside a partly filled
bin and the modal bin immediately to its right is not a 0.019-degree discrepancy.
It also identifies isolated 41--43 degree bins at roughly one count per bin as
Poisson noise, and records that the upper angular window is temporal selection.

## Reproducibility

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/build_step4_pairs.py --processes 4
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step4.py
```

Each figure has PDF, ROOT, CSV and metadata sidecars. No push, merge, deck edit, or
simulation was performed.
