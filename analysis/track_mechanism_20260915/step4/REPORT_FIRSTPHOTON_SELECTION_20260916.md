# EXEC_46 Step 4 — first-photon selection and source order statistics

Date: 2026-09-16

## Checkpoint verdict

Step 4 is complete. No simulation was run; the 21 production ROOT files were read-only.
The corrected min-versus-min comparison closes the near-END winner accounting: Cherenkov
pays a 45--51 ps transport handicap but gains 54--71 ps in the creation time attached to
the first detected photon, leaving it earlier by 4--26 ps depending on material.

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

The requested beta=1 edge is 146.903 mm/ns. The configured
1 GeV muon gives beta=0.995424 and 146.450 mm/ns.
Axial velocity is the measured axial displacement divided by pure propagation time.

| material | d [mm] | proxy N_C low/high quintile | v low/high [mm/ns] | high-edge [mm/ns] |
|---|---:|---:|---:|---:|
| EJ-200 | 50 | 8.56 / 39.56 | 140.389 / 145.771 | -0.679 |
| EJ-200 | 200 | 4.60 / 23.66 | 139.124 / 146.096 | -0.354 |
| EJ-200 | 500 | 2.17 / 13.08 | 137.028 / 146.195 | -0.255 |
| EJ-200 | 700 | 1.42 / 8.97 | 128.339 / 145.514 | -0.935 |
| EJ-200 | 900 | 1.15 / 7.30 | 128.350 / 146.082 | -0.368 |
| EJ-200 | 1200 | 1.00 / 5.39 | 124.466 / 144.357 | -2.093 |
| EJ-200 | 1350 | 1.00 / 5.42 | 130.035 / 145.385 | -1.064 |
| EJ-204 | 50 | 8.08 / 37.35 | 140.914 / 145.740 | -0.710 |
| EJ-204 | 200 | 3.88 / 20.21 | 138.351 / 146.069 | -0.380 |
| EJ-204 | 500 | 1.56 / 10.14 | 133.471 / 146.097 | -0.353 |
| EJ-204 | 700 | 1.00 / 6.43 | 126.638 / 145.123 | -1.327 |
| EJ-204 | 900 | 1.00 / 5.01 | 129.143 / 145.640 | -0.809 |
| EJ-204 | 1200 | 1.00 / 3.45 | 126.796 / 142.943 | -3.507 |
| EJ-204 | 1350 | 1.00 / 3.45 | 133.874 / 144.441 | -2.009 |
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
Delta t(alpha,d) = d/(c/n) * [1/cos(alpha) - 1/cos(alpha_edge)]
```

The stored final angle is folded to the axial magnitude because the penalty uses
|cos(alpha)|. The upper window is a selection in arrival time, not a second cone
boundary.

| material | d [mm] | q95 alpha proxy [deg] | predicted penalty [ps] | empirical q95 excess [ps] |
|---|---:|---:|---:|---:|
| EJ-200 | 50 | 41.756 | 11.73 | 11.77 |
| EJ-200 | 200 | 43.354 | 83.77 | 84.05 |
| EJ-200 | 500 | 45.718 | 359.81 | 364.07 |
| EJ-200 | 700 | 50.972 | 1078.20 | 1204.42 |
| EJ-200 | 900 | 51.339 | 1446.60 | 1924.92 |
| EJ-200 | 1200 | 53.949 | 2551.53 | 4366.29 |
| EJ-200 | 1350 | 54.019 | 2891.06 | 2938.37 |
| EJ-204 | 50 | 41.753 | 11.72 | 11.81 |
| EJ-204 | 200 | 43.918 | 97.41 | 97.47 |
| EJ-204 | 500 | 48.906 | 594.40 | 597.60 |
| EJ-204 | 700 | 53.858 | 1474.38 | 1570.25 |
| EJ-204 | 900 | 53.923 | 1908.44 | 2316.16 |
| EJ-204 | 1200 | 54.279 | 2637.51 | 4469.57 |
| EJ-204 | 1350 | 54.469 | 3023.70 | 3665.82 |
| EJ-230 | 50 | 41.849 | 12.24 | 12.28 |
| EJ-230 | 200 | 45.059 | 126.25 | 127.36 |
| EJ-230 | 500 | 49.604 | 651.42 | 657.56 |
| EJ-230 | 700 | 53.979 | 1492.49 | 1562.89 |
| EJ-230 | 900 | 53.980 | 1919.49 | 2274.38 |
| EJ-230 | 1200 | 54.262 | 2632.99 | 4371.43 |
| EJ-230 | 1350 | 54.557 | 3050.14 | 3642.78 |

This distance-only comparison is not a valid rejection of angular selection:
the low-N_C quintile is N_C=1 at d >= 700 mm, where no order statistic exists.
D1 below conditions on distance and uses N_C as the control variable.

### D1 — angular q95 at fixed distance versus N_C

| material | d [mm] | low/high mean N_C | low/high q95 [deg] | change [deg] |
|---|---:|---:|---:|---:|
| EJ-200 | 50 | 8.56/39.56 | 52.492/40.578 | -11.913 |
| EJ-200 | 200 | 4.60/23.66 | 53.461/40.065 | -13.396 |
| EJ-200 | 500 | 2.17/13.08 | 55.186/39.953 | -15.233 |
| EJ-204 | 50 | 8.08/37.35 | 50.242/40.604 | -9.638 |
| EJ-204 | 200 | 3.88/20.21 | 55.524/40.081 | -15.444 |
| EJ-204 | 500 | 1.56/10.14 | 60.488/40.214 | -20.274 |
| EJ-230 | 50 | 7.84/36.98 | 52.478/40.643 | -11.836 |
| EJ-230 | 200 | 3.53/18.63 | 55.222/40.199 | -15.023 |
| EJ-230 | 500 | 1.42/8.77 | 60.551/40.370 | -20.182 |

q95 narrows from the lowest to highest N_C quintile in 9/9 fixed-distance groups. The distance-only test is
therefore reclassified as badly conditioned rather than a refutation of the
cone-edge order-statistics mechanism.

## C0c — corrected d=50 mm handicap

| material | mirror | N_S | pure-min creation C-S [ps] | selected creation C-S [ps] | transport C-S [ps] | total C-S [ps] | Cher wins | reorder gap [ps] |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | x=-650/left | 2271.4 | -40.32 | -70.60 | 45.02 | -25.58 | 74.10% | 30.28 |
| EJ-200 | x=+650/right | 2283.4 | -40.13 | -70.44 | 44.73 | -25.71 | 73.56% | 30.32 |
| EJ-204 | x=-650/left | 2335.5 | -34.08 | -61.07 | 48.04 | -13.03 | 62.90% | 26.98 |
| EJ-204 | x=+650/right | 2324.9 | -33.98 | -61.81 | 48.32 | -13.49 | 63.42% | 27.84 |
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
| EJ-200 | N^-1 | 1.0 | 1511.78/12 = 125.98 | -- |
| EJ-200 | N^-1/2 | 0.5 | 13.73/12 = 1.14 | 0.8125 +/- 0.0062 |
| EJ-200 | free | 0.5285 [0.5170, 0.5400] | 7.93/11 = 0.72 | -- |
| EJ-204 | N^-1 | 1.0 | 3807.89/12 = 317.32 | -- |
| EJ-204 | N^-1/2 | 0.5 | 11.49/12 = 0.96 | 0.7722 +/- 0.0049 |
| EJ-204 | free | 0.5105 [0.5030, 0.5175] | 9.69/11 = 0.88 | -- |
| EJ-230 | N^-1 | 1.0 | 5607.00/12 = 467.25 | -- |
| EJ-230 | N^-1/2 | 0.5 | 39.20/12 = 3.27 | 0.7361 +/- 0.0043 |
| EJ-230 | free | 0.5115 [0.5055, 0.5175] | 35.69/11 = 3.24 | -- |

The N^-1/2 model gives chi2/ndf 0.96--3.27; N^-1 is rejected with chi2/ndf 125.98--467.25. The fitted effective fractions 0.812, 0.772 and 0.736
differ by at least 5.1 formal standard deviations, so N_eff/N_scint is not material-independent within this model.
A universal geometric factor is therefore rejected. The preregistered fractions
needed to force the detection-selected handicap, 0.20/0.13/0.09, are even more
strongly material dependent. The low-boundary fractions span 0.019--0.271, 0.021--0.278, and 0.021--0.280; they are strongly distance dependent and do not reproduce the fitted effective fractions. `scintillation_order_points.csv` gives fitted N_eff, the low-boundary count, and their ratio at each of the fourteen mirror-resolved points per material.

### D2 — correction for primary-muon transit

The gun points along -z and enters the 10-mm bar at z=+5 mm. The corrected
creation coordinate is `t_creation-(5 mm-z_creation)/c`; a common upstream flight
offset is absorbed by the fitted intercept.

| material | raw f_eff | corrected f_eff | corrected exponent | corrected chi2/ndf |
|---|---:|---:|---:|---:|
| EJ-200 | 0.8125 +/- 0.0062 | 0.8522 +/- 0.0064 | 0.5525 | 2.33 |
| EJ-204 | 0.7722 +/- 0.0049 | 0.8179 +/- 0.0052 | 0.5345 | 2.39 |
| EJ-230 | 0.7361 +/- 0.0043 | 0.7934 +/- 0.0047 | 0.5435 | 6.54 |

The raw effective-fraction ordering follows the d=50-mm asymptotic minimum
delays 30.2 > 24.7 > 20.6 ps: the faster material suffers the larger fractional
dilution from the same 33.4-ps traversal. Subtracting the transit moves every
fraction toward one but recovers only 20.1--21.7% of the original
deficit. N_eff does not reach N_scint; primary transit is a real contribution but
does not explain the remaining 15--21% deficit. The corrected free exponents also
remain above 0.5, so the exact i.i.d. common-origin law is not restored.

## C2 — fixed-point Cherenkov/width test

No interpolation is used.

| d [mm] | Pearson r | Spearman rho | role |
|---:|---:|---:|---|
| 50 | -0.9955 | -1.0000 | Cherenkov-dominated operating regime |
| 200 | +0.9484 | +0.9429 | spectator regime |

The six measured values at each fixed point are:

| d [mm] | material | mirror | fCher(first) | sigma_IQR(near)/sigma_IQR(T0) |
|---:|---|---|---:|---:|
| 50 | EJ-200 | x=-650/left | 74.10% | 0.1299 |
| 50 | EJ-200 | x=+650/right | 73.56% | 0.1389 |
| 50 | EJ-204 | x=-650/left | 62.90% | 0.1884 |
| 50 | EJ-204 | x=+650/right | 63.42% | 0.1841 |
| 50 | EJ-230 | x=-650/left | 51.77% | 0.2313 |
| 50 | EJ-230 | x=+650/right | 52.33% | 0.2249 |
| 200 | EJ-200 | x=-500/left | 8.44% | 0.7515 |
| 200 | EJ-200 | x=+500/right | 8.90% | 0.7395 |
| 200 | EJ-204 | x=-500/left | 7.44% | 0.6727 |
| 200 | EJ-204 | x=+500/right | 7.64% | 0.6737 |
| 200 | EJ-230 | x=-500/left | 7.18% | 0.6328 |
| 200 | EJ-230 | x=+500/right | 7.02% | 0.6090 |

At 50 mm, more Cherenkov corresponds to a smaller near-END/T0 robust-width
ratio. At 200 mm, where the Cherenkov fraction is about 7--9%, the correlation
reverses. This sign change bounds the operating range of the mechanism; it is not
evidence that either measured fixed-point correlation is internally inconsistent.

## C3 — physical source mixture at the near END

| material | mirror | fCher | separation [ps] | sigma_G [ps] | chi2/ndf | sigma_mix [ps] | qwidth [ps] | angle-time r |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | x=-650/left | 74.10% | -12.21 | 2.34 | 205.2 | 18.92 | 14.19 | +0.984 |
| EJ-200 | x=+650/right | 73.56% | -12.01 | 2.32 | 210.1 | 19.67 | 15.57 | +0.962 |
| EJ-204 | x=-650/left | 62.90% | -11.95 | 10.68 | 460.2 | 19.56 | 17.34 | +0.986 |
| EJ-204 | x=+650/right | 63.42% | -12.02 | 11.05 | 503.7 | 19.84 | 17.27 | +0.962 |
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
| EJ-200 | 50 | 73.830% | 3.190% | 23.144 | -3.0734 |
| EJ-200 | 200 | 8.670% | 3.240% | 2.676 | -3.3617 |
| EJ-200 | 500 | 5.925% | 3.365% | 1.761 | -3.7803 |
| EJ-200 | 700 | 5.030% | 2.990% | 1.682 | -4.0592 |
| EJ-200 | 900 | 4.610% | 3.025% | 1.524 | -4.2062 |
| EJ-200 | 1200 | 3.360% | 3.030% | 1.109 | -4.4985 |
| EJ-200 | 1350 | 2.985% | 3.290% | 0.907 | -4.5416 |
| EJ-204 | 50 | 63.160% | 3.035% | 20.811 | -2.6009 |
| EJ-204 | 200 | 7.540% | 3.070% | 2.456 | -2.8526 |
| EJ-204 | 500 | 4.935% | 2.970% | 1.662 | -3.1715 |
| EJ-204 | 700 | 4.335% | 2.590% | 1.674 | -3.3417 |
| EJ-204 | 900 | 3.540% | 2.705% | 1.309 | -3.4695 |
| EJ-204 | 1200 | 2.740% | 2.565% | 1.068 | -3.7133 |
| EJ-204 | 1350 | 2.200% | 2.680% | 0.821 | -3.7370 |
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
