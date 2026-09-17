# EXEC_46 Step 3 — pure transport g(d) and Cherenkov guiding

Date: 2026-09-16

## Checkpoint verdict

Step 3 is complete. No simulation was run and all production ROOT files were read-only.
Source-separated fitted velocities are compared below with wavelength-resolved predictions.
The all-photon means are much slower because they include recirculated paths. None of
the low-order global models has acceptable absolute chi-square, so these slopes are
diagnostic summaries rather than complete models of g(d).

This is the required checkpoint. Step 4 has not been started.

## Definitions and input checks

The production macros contain `/muon/angle 0` in all 21 cells. The source maps this
to momentum `(0,0,-1)`, perpendicular to the bar x axis. Each actual cell runtime
supplies its own RINDEX; constant and dispersive models are both accepted.
The unchanged Step 3 builder omits wavelengths. These are joined read-only from
sipm_hits by cell, event, face and selected track ID; timestamps and source labels
must match exactly. Photon choices and all original derived columns are preserved.
The configured muon beta is 0.995423528. There is no scalar material cone angle.

Median [q05, q95]; wavelengths weighted by the selected photons. Cone: created wavelength; transport: detected wavelength. Spectral correction is (1/vg_det - 1/vg_created), in ps/mm. Constant-index models have no measured dispersion domain.

| material | population | n(created) | theta_C(beta) [deg] | critical(created) [deg] | edge(beta) [deg] | vg(det) [mm/ns] | transport edge [mm/ns] | spectral delay [ps/mm] | group-phase delay [ps/mm] | created clamp / outside fraction | detected clamp / outside fraction |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EJ-200 | first_by_source, source 1 | 1.620 [1.603, 1.629] | 51.663 [51.205, 51.935] | 38.130 [37.860, 38.585] | 38.337 [38.065, 38.795] | 165.948 [164.349, 168.858] | 130.166 [129.394, 131.606] | 0 / 0 / 0 | 0.623652 / 0.57375 / 0.649603 | 0.0000% / 0.0000% | 0.0000% / 0.0000% |
| EJ-200 | first_by_source, source 2 | 1.617 [1.557, 1.644] | 51.594 [49.826, 52.333] | 38.199 [37.465, 39.953] | 38.406 [37.667, 40.174] | 166.377 [162.221, 181.816] | 130.390 [128.393, 138.940] | 0 / 0 / 0 | 0.61575 / 0.304858 / 0.68196 | 6.1047% / 6.1047% | 6.1047% / 6.1047% |
| EJ-204 | first_by_source, source 1 | 1.618 [1.607, 1.625] | 51.627 [51.305, 51.820] | 38.165 [37.975, 38.485] | 38.373 [38.180, 38.695] | 172.332 [170.354, 175.977] | 135.107 [133.910, 137.348] | 0 / 0 / 0 | 0.404675 / 0.322436 / 0.449059 | 0.0000% / 0.0000% | 0.0000% / 0.0000% |
| EJ-204 | first_by_source, source 2 | 1.630 [1.585, 1.633] | 51.948 [50.681, 52.038] | 37.847 [37.758, 39.104] | 38.052 [37.962, 39.319] | 171.478 [168.632, 185.001] | 134.998 [132.924, 143.129] | 0 / 0 / 0 | 0.398191 / 0.116542 / 0.484675 | 49.4003% / 49.4003% | 49.4003% / 49.4003% |
| EJ-230 | first_by_source, source 1 | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 3.5086% | 0.0000% / 3.5086% |
| EJ-230 | first_by_source, source 2 | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 36.7963% | 0.0000% / 36.7963% |
| EJ-200 | first_overall, source 0 | 1.619 [1.575, 1.631] | 51.640 [50.381, 51.972] | 38.153 [37.823, 39.402] | 38.360 [38.028, 39.619] | 166.093 [164.166, 174.944] | 130.241 [129.314, 134.759] | 0 / 0 / 0 | 0.621007 / 0.460985 / 0.652179 | 3.1257% / 3.1257% | 3.1257% / 3.1257% |
| EJ-204 | first_overall, source 0 | 1.618 [1.593, 1.627] | 51.615 [50.915, 51.864] | 38.177 [37.931, 38.872] | 38.385 [38.136, 39.085] | 172.465 [170.203, 181.216] | 135.208 [133.845, 140.671] | 0 / 0 / 0 | 0.400538 / 0.202833 / 0.451065 | 4.9736% / 4.9736% | 4.9736% / 4.9736% |
| EJ-230 | first_overall, source 0 | 1.580 [1.580, 1.580] | 50.519 [50.519, 50.519] | 39.265 [39.265, 39.265] | 39.481 [39.481, 39.481] | 189.742 [189.742, 189.742] | 146.450 [146.450, 146.450] | 0 / 0 / 0 | 0 / 0 / 0 | 0.0000% / 7.1807% | 0.0000% / 7.1807% |

| material | population | phase speed(created) [mm/ns] | group index(det) | phase edge(beta=1) [mm/ns] |
|---|---|---|---|---|
| EJ-200 | first_by_source, source 1 | 185.105 [183.992, 186.973] | 1.807 [1.775, 1.824] | 145.607 [145.264, 146.154] |
| EJ-200 | first_by_source, source 2 | 185.389 [182.357, 192.515] | 1.802 [1.649, 1.848] | 145.692 [144.741, 147.576] |
| EJ-204 | first_by_source, source 1 | 185.252 [184.466, 186.563] | 1.740 [1.704, 1.760] | 145.651 [145.411, 146.037] |
| EJ-204 | first_by_source, source 2 | 183.938 [183.570, 189.089] | 1.748 [1.620, 1.778] | 145.248 [145.132, 146.733] |
| EJ-230 | first_by_source, source 1 | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-230 | first_by_source, source 2 | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |
| EJ-200 | first_overall, source 0 | 185.201 [183.839, 190.297] | 1.805 [1.714, 1.826] | 145.635 [145.217, 147.044] |
| EJ-204 | first_overall, source 0 | 185.300 [184.284, 188.146] | 1.738 [1.654, 1.761] | 145.665 [145.355, 146.480] |
| EJ-230 | first_overall, source 0 | 189.742 [189.742, 189.742] | 1.580 [1.580, 1.580] | 146.903 [146.903, 146.903] |

The beta=1 identity arccos(sin(theta_C(lambda))) = arcsin(1/n(lambda)) is checked at every table node and evaluated photon wavelength; dispersion replicates the identity wavelength by wavelength. Arrival-time edge speeds use vg(lambda_det)*sin(theta_C(lambda_created)); the phase-edge formula c*sqrt(1-1/n^2)/n is retained separately. Multi-medium/reflected paths are not reconstructed by this homogeneous-scintillator reference.

For every photon:

```text
tprop = t_detection_ns - t_creation_ns
d_direct = |x_detection - x_creation| in three dimensions
rho_detour = path_length_mm / d_direct
v_apparent = d_direct / tprop
```

The primary g(d) calculation uses microscopic d_direct. Nominal distances are used
only to pair the left and right mirror realizations at 50, 200, 500, 700, 900, 1200
and 1350 mm. At 700 mm the two faces share the same x=0 simulated events, so they are
symmetry realizations but are not statistically independent in the strict sampling
sense; the all-photon summary does not store their covariance.

## B1 — withdrawn sparse-grid boundaries

`BOUNDARY_IDENTITY_NOT_ESTABLISHED` is retained and A2 is withdrawn. The correlation
changes sign between |x|=500 and 650 mm. Both reported linear crossings are artifacts
of interpolating across a regime change. The available grid samples only END distances
{50, 200, 500, 700, 900, 1200, 1350} mm; the transition lies wholly in the 50--200 mm
gap. No new simulation is proposed here.

## Fits to mirror-combined g(d)

The mandated linear fit is through the origin, t=d/v. The quadratic and cubic models
are also through the origin. Errors are event-cluster SEMs for all photons and event
SEMs for first-by-source photons. The shaded bands in `g_d.pdf` span the two mirror
means rather than pretending that a photon-IID error describes an event cluster.

### Linear summaries

| sample | material | source | v [mm/ns] | chi2/ndf | vs first-source transport edge median [q05,q95] | vs first-source vg median [q05,q95] |
|---|---|---|---:|---:|---:|---:|
| all_photons | EJ-200 | scintillation | 133.356 +/- 0.009 | 3839379.5/6 = 639896.6* | +2.45% [+1.33, +3.06] | -19.64% [-21.02, -18.86] |
| all_photons | EJ-200 | Cherenkov | 112.976 +/- 0.056 | 79180.6/6 = 13196.8* | -13.36% [-18.69, -12.01] | -32.10% [-37.86, -30.36] |
| all_photons | EJ-204 | scintillation | 136.571 +/- 0.008 | 5766968.0/6 = 961161.3* | +1.08% [-0.57, +1.99] | -20.75% [-22.39, -19.83] |
| all_photons | EJ-204 | Cherenkov | 120.173 +/- 0.041 | 161545.9/6 = 26924.3* | -10.98% [-16.04, -9.59] | -29.92% [-35.04, -28.74] |
| all_photons | EJ-230 | scintillation | 138.520 +/- 0.008 | 6025468.8/6 = 1004244.8* | -5.41% [-5.41, -5.41] | -27.00% [-27.00, -27.00] |
| all_photons | EJ-230 | Cherenkov | 119.639 +/- 0.040 | 184914.5/6 = 30819.1* | -18.31% [-18.31, -18.31] | -36.95% [-36.95, -36.95] |
| first_by_source | EJ-200 | scintillation | 183.908 +/- 0.010 | 48728.9/6 = 8121.5* | +41.29% [+39.74, +42.13] | +10.82% [+8.91, +11.90] |
| first_by_source | EJ-200 | Cherenkov | 149.026 +/- 0.025 | 3756.9/6 = 626.1* | +14.29% [+7.26, +16.07] | -10.43% [-18.03, -8.13] |
| first_by_source | EJ-204 | scintillation | 184.429 +/- 0.009 | 36566.5/6 = 6094.4* | +36.51% [+34.28, +37.73] | +7.02% [+4.80, +8.26] |
| first_by_source | EJ-204 | Cherenkov | 149.933 +/- 0.024 | 3807.1/6 = 634.5* | +11.06% [+4.75, +12.80] | -12.56% [-18.96, -11.09] |
| first_by_source | EJ-230 | scintillation | 186.155 +/- 0.010 | 5571.4/6 = 928.6* | +27.11% [+27.11, +27.11] | -1.89% [-1.89, -1.89] |
| first_by_source | EJ-230 | Cherenkov | 147.863 +/- 0.028 | 6938.3/6 = 1156.4* | +0.96% [+0.96, +0.96] | -22.07% [-22.07, -22.07] |
| first_overall | EJ-200 | first overall | 183.899 +/- 0.010 | 120060.6/6 = 20010.1* | +41.20% [+36.46, +42.21] | +10.72% [+5.12, +12.02] |
| first_overall | EJ-204 | first overall | 184.389 +/- 0.009 | 78837.0/6 = 13139.5* | +36.37% [+31.08, +37.76] | +6.91% [+1.75, +8.33] |
| first_overall | EJ-230 | first overall | 186.113 +/- 0.010 | 21371.5/6 = 3561.9* | +27.08% [+27.08, +27.08] | -1.91% [-1.91, -1.91] |

`*` marks an inadequate absolute fit (all entries above).

### Curvature tests

| sample | material | source | Delta chi2 linear->quadratic | |b2|/err | quadratic chi2/ndf | cubic chi2/ndf |
|---|---|---|---:|---:|---:|---:|
| all_photons | EJ-200 | Cherenkov | 57260.0 | 239.3 | 21920.6/5 = 4384.1* | 9981.3/4 = 2495.3* |
| all_photons | EJ-200 | scintillation | 2406887.9 | 1551.4 | 1432491.6/5 = 286498.3* | 567262.4/4 = 141815.6* |
| all_photons | EJ-204 | Cherenkov | 108801.6 | 329.9 | 52744.3/5 = 10548.9* | 23195.3/4 = 5798.8* |
| all_photons | EJ-204 | scintillation | 3505432.5 | 1872.3 | 2261535.5/5 = 452307.1* | 837363.2/4 = 209340.8* |
| all_photons | EJ-230 | Cherenkov | 115850.0 | 340.4 | 69064.4/5 = 13812.9* | 26066.7/4 = 6516.7* |
| all_photons | EJ-230 | scintillation | 3603775.4 | 1898.4 | 2421693.4/5 = 484338.7* | 865787.4/4 = 216446.8* |
| first_by_source | EJ-200 | Cherenkov | 835.4 | 28.9 | 2921.5/5 = 584.3* | 938.3/4 = 234.6* |
| first_by_source | EJ-200 | scintillation | 31967.4 | 178.8 | 16761.5/5 = 3352.3* | 4712.2/4 = 1178.0* |
| first_by_source | EJ-204 | Cherenkov | 1602.0 | 40.0 | 2205.1/5 = 441.0* | 389.3/4 = 97.3* |
| first_by_source | EJ-204 | scintillation | 23166.2 | 152.2 | 13400.3/5 = 2680.1* | 3289.9/4 = 822.5* |
| first_by_source | EJ-230 | Cherenkov | 5948.6 | 77.1 | 989.6/5 = 197.9* | 329.5/4 = 82.4* |
| first_by_source | EJ-230 | scintillation | 3457.4 | 58.8 | 2114.0/5 = 422.8* | 158.5/4 = 39.6* |
| first_overall | EJ-200 | first overall | 42014.5 | 205.0 | 78046.1/5 = 15609.2* | 44583.8/4 = 11146.0* |
| first_overall | EJ-204 | first overall | 26994.8 | 164.3 | 51842.2/5 = 10368.4* | 32630.0/4 = 8157.5* |
| first_overall | EJ-230 | first overall | 4725.4 | 68.7 | 16646.1/5 = 3329.2* | 12025.0/4 = 3006.3* |

The residual-guided cubic improves chi-square again but remains rejected.
The opening into boundary-count families explains why a single low-order curve
does not describe the conditional transport mean.

## Comparison with the historical first-PE effective velocity

The historical estimator mixes emission and transport and reported the following
finite-difference ranges. The current first-overall row removes creation time event
by event before fitting, while the source-separated rows expose the mixture.

| material | historical local range [mm/ns] | historical linear | pure first-overall | first scintillation | first Cherenkov |
|---|---:|---:|---:|---:|---:|
| EJ-200 | 172.63--183.57 | 182.216 | 183.899 | 183.908 | 149.026 |
| EJ-204 | 175.04--182.56 | 181.652 | 184.389 | 184.429 | 149.933 |
| EJ-230 | 177.26--182.76 | 181.611 | 186.113 | 186.155 | 147.863 |

The source-specific Cherenkov value is not expected to equal the historical
mixed estimator. Comparisons use photon-weighted group-speed and transport-edge
distributions, not c/n. The all-photon rows use explicitly labelled first-source
reference distributions; aggregate all-photon moments cannot reconstruct spectra.

## B2 — mirror consistency

The full table is `mirror_consistency.csv`: 84 source-separated rows plus 21
first-overall diagnostic rows. Compact maxima are:

| sample | material | source | max |Delta mean| [ps] | max |pull| |
|---|---|---|---:|---:|
| all_photons | EJ-200 | Cherenkov | 32.820 | 1.36 |
| all_photons | EJ-200 | scintillation | 4.140 | 1.78 |
| all_photons | EJ-204 | Cherenkov | 16.212 | 1.06 |
| all_photons | EJ-204 | scintillation | 3.491 | 1.27 |
| all_photons | EJ-230 | Cherenkov | 14.693 | 1.35 |
| all_photons | EJ-230 | scintillation | 1.882 | 1.20 |
| first_by_source | EJ-200 | Cherenkov | 8.165 | 0.62 |
| first_by_source | EJ-200 | scintillation | 1.085 | 1.90 |
| first_by_source | EJ-204 | Cherenkov | 10.683 | 1.46 |
| first_by_source | EJ-204 | scintillation | 0.670 | 0.87 |
| first_by_source | EJ-230 | Cherenkov | 5.477 | 1.39 |
| first_by_source | EJ-230 | scintillation | 2.454 | 1.59 |
| first_overall | EJ-200 | first overall | 0.745 | 0.85 |
| first_overall | EJ-204 | first overall | 0.839 | 1.34 |
| first_overall | EJ-230 | first overall | 2.629 | 1.60 |

The largest pull is 2.36 and the largest absolute mirror difference is
39.654 ps. The two realizations are consistent at the declared precision.

## B3 — direct Cherenkov edge test

The primary-cone geometry applies because the gun is perpendicular to x. The
identity with theta_crit is exact only in the beta=1 limit; the configured finite
beta predicts a wavelength-dependent lower edge, tabulated in cherenkov_edge_summary.csv.
However, `source_type==2` records the creator process and not parent track identity.
It therefore includes Cherenkov photons made by secondary charged particles. The
`primary_like_proxy` requires creation x within 0.001 mm of gun x and creation y
within 0.001 mm of zero; it is a geometric proxy, not a recovered parent label.

| material | selection | N | q0.1% [deg] | median [deg] | fraction below theta_crit | modal 0.02-deg bin |
|---|---|---:|---:|---:|---:|---:|
| EJ-200 | all_source_type_2 | 20000 | 2.546 | 39.364 | 11.0050% | 40.20--40.22 |
| EJ-200 | primary_like_proxy | 7592 | 37.480 | 40.076 | 0.0000% | 40.18--40.20 |
| EJ-204 | all_source_type_2 | 20000 | 2.365 | 39.158 | 11.9550% | 39.34--39.36 |
| EJ-204 | primary_like_proxy | 7541 | 37.926 | 39.435 | 0.0265% | 39.34--39.36 |
| EJ-230 | all_source_type_2 | 20000 | 1.650 | 39.644 | 11.5400% | 39.50--39.52 |
| EJ-230 | primary_like_proxy | 11025 | 39.425 | 39.700 | 0.0181% | 39.50--39.52 |

The finite-beta edge is a distribution for dispersive tables. Its modal bin
must be compared with photon-wise edge predictions, not a single 39.5-degree bin.
The upper selected angular window is temporal selection, not another cone boundary.
The beta=1 identity is exact at each wavelength. The fraction below the critical
angle above uses each photon's created wavelength. A primary-like proxy does not
establish parent track identity; secondary contamination remains a limitation.

## B4 — source-separated transport

The prediction distributions above give phase speed, group transport, and finite-beta
cone-edge transport references separately. The fitted-velocity comparison reports
median [q05,q95] differences; no agreement verdict from the old constant-index
campaign is carried over. Significant fitted curvature remains a limitation of a
single global linear summary.

All detected photon linear summaries are path-population
means, not material group velocities. `g_by_boundary.pdf` separates them into the
mandated 0, 1--2, 3--5, 6--10 and >10 encounter families; `g_by_exit_angle.pdf`
shows the corresponding final-angle stratification.

## Microscopic binning and TH2 coverage

The 10 mm microscopic bins contain 314,688,168 photons. The 0--60 ns TH2
excludes 1,528 upper-overflow photons (0.000486%) only
from the raster; all moments and fits retain them. `tprop_variance.csv` contains
the binned mean and conditional variance by material, face and source.

## Reproducibility

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/build_step3_transport.py --processes 4
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step3.py
```

Primary artifacts:

- `all_photon_cell_face.csv`: event-cluster means for 84 cell/face/source groups.
- `first_by_source.root`: first photon per event, END face and source.
- `fit_summary.csv`, `mirror_consistency.csv`, and `cherenkov_edge_summary.csv`.
- Each figure has `.pdf`, `.root`, `.csv`, and `.meta.json` sidecars.

The build wall time was 226.2 s with four read processes. No push, merge, deck
edit, or Geant4 run was performed.
