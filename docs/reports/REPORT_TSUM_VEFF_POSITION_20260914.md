# Test of position independence and constant effective velocity

**Date:** 2026-09-14  
**Input:** 21 current corrected-transport EndTop cells, 10,000 events per cell  
**Simulation performed:** none  
**Classification:** `V_EFF_IS_DISTANCE_DEPENDENT` — symmetric nonlinear first-arrival response

## Short physics verdict

1. **Is `<T0>` constant with x?** No. A constant is rejected for all three materials. The average center-to-`|x|=650 mm` shift is `-31.92 ps` for EJ-200, `-14.31 ps` for EJ-204 and `+0.19 ps` for EJ-230. EJ-230 nevertheless has a resolved nonmonotonic position dependence with an 8.35 ps peak-to-peak excursion.
2. **Is the odd/linear component compatible with zero?** Yes at the available precision. Full-quadratic odd coefficients have z values `-0.17`, `-0.97`, and `+2.30` for EJ-200/204/230. Mirror and equal-distance LEFT/RIGHT tests show no coherent odd pattern; the largest individual deviations are about 2.0--2.4 sigma.
3. **Is there significant even/quadratic curvature?** EJ-200 and EJ-204 have highly significant quadratic coefficients, `-71.81 +/- 1.80` and `-29.93 +/- 1.79 ps/m^2`. The EJ-230 quadratic coefficient alone is only `+2.65 +/- 1.75 ps/m^2`, because its even dependence is hump-shaped rather than parabolic. An even-quartic diagnostic improves its constant-fit chi-square by 89.83 for two parameters and describes the seven points with `chi2/ndf=8.81/4`.
4. **Is a single `v_eff` adequate?** No. The linear propagation fits have `chi2/ndf=11004/12`, `5492/12`, and `2805/12`. Quadratic fits improve chi-square by 7300, 2875, and 830 respectively but remain unacceptable, so even one quadratic `v_eff(d)` is only a diagnostic approximation.
5. **Do LEFT and RIGHT collapse onto the same distance curve?** Yes. Excluding the correlated `x=0` pair, equal-distance collapse tests give probabilities 0.157, 0.715 and 0.147 for EJ-200/204/230. The nonlinear response is common to both sides.
6. **Does mean-first-five improve position independence?** No. It reduces the event-by-event timing width, but for EJ-230 its mean has a 17.48 ps peak-to-peak excursion versus 8.35 ps for first PE. Its constant fit gives `chi2/ndf=637.6/6`, and the even-quadratic coefficient is `+13.26 +/- 1.62 ps/m^2`.
7. **Can the microscopic cause be determined from the current TTree?** No: `MECHANISM_NOT_DIRECTLY_RESOLVABLE_FROM_CURRENT_TTREE`. Detected wavelength is available and shows distance-dependent selection in the first-arrival subset, but accumulated path, boundary count, creation wavelength and a track identifier linking those quantities to a detected SiPM hit are absent.

The result satisfies **Case B — symmetric nonlinear propagation**. Here `g(d)` is the mean response of the first detected photoelectron. It includes propagation, scintillation-emission order statistics, PDE and distance-dependent light collection; it is not the velocity of one optical ray.

## Data and event reconstruction

The analysis reads `/home/rrios/ej200/presentations/v9/sources/timing_events.root`, a ROOT extraction of the 21 production files in `/home/rrios/exec42_20260913/grid/cells/`. The source has 210,000 events and preserves the first 20 sorted detected PE times at both ENDs plus event-level light yields.

For every event:

- LEFT END is global IDs 0--7 and RIGHT END is global IDs 8--15;
- `tL` and `tR` are the earliest detected `time_ns` on each side;
- `Tsum=tL+tR`, `T0=Tsum/2`, and `DeltaT=tR-tL`;
- `Npe_END=Npe_L+Npe_R`;
- mean-first-five is the arithmetic mean of the first five sorted times independently at each END.

No fitted timing-width summary enters the mean calculation. All means and SEMs come from the 10,000 event values in each cell. The derived event tree is `derived_events` in `sources/tsum_veff.root`.

## Measured mean T0

Values below are direct event means; uncertainties are event-level SEM.

| Material | x=-650 | -500 | -200 | 0 | +200 | +500 | +650 mm |
|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | 4.09550 +/- 0.00081 | 4.12153 +/- 0.00081 | 4.12656 +/- 0.00077 | 4.12757 +/- 0.00077 | 4.12587 +/- 0.00079 | 4.12108 +/- 0.00082 | 4.09579 +/- 0.00081 ns |
| EJ-204 | 4.10270 +/- 0.00085 | 4.11834 +/- 0.00082 | 4.11765 +/- 0.00074 | 4.11649 +/- 0.00072 | 4.11594 +/- 0.00074 | 4.11851 +/- 0.00082 | 4.10167 +/- 0.00085 ns |
| EJ-230 | 4.09949 +/- 0.00085 | 4.10557 +/- 0.00079 | 4.10065 +/- 0.00069 | 4.09982 +/- 0.00068 | 4.10235 +/- 0.00070 | 4.10785 +/- 0.00081 | 4.10054 +/- 0.00085 ns |

This separates the present question from timing resolution. The varying quantity is `E[T0|x]`; no statement about a distribution width is needed to reject position independence.

## Position fits of T0 and Tsum

ROOT `TGraphErrors` and `TF1` fits use x in metres and Minuit2/Migrad with tolerance 0.01. Every fit returned status 0 and covariance quality 3. Complete covariance matrices are stored under `position_fits/` in `tsum_veff.root`.

### T0 fits

Coefficients use ns, ns/m and ns/m2. Blank coefficients are fixed absent by the model.

| Material | model | chi2/ndf | probability | a0 [ns] | a1 [ns/m] | a2 [ns/m2] |
|---|---|---:|---:|---:|---:|---:|
| EJ-200 | constant | 1909.1/6 | underflow | 4.116651 +/- 0.000301 | | |
| EJ-200 | linear | 1909.1/5 | underflow | 4.116650 +/- 0.000301 | -0.000057 +/- 0.000679 | |
| EJ-200 | even quadratic | 313.1/5 | 1.53e-65 | 4.130773 +/- 0.000465 | | -0.071809 +/- 0.001797 |
| EJ-200 | full quadratic | 313.1/4 | 1.64e-66 | 4.130773 +/- 0.000465 | -0.000118 +/- 0.000679 | -0.071810 +/- 0.001797 |
| EJ-204 | constant | 486.9/6 | 5.66e-102 | 4.113573 +/- 0.000298 | | |
| EJ-204 | linear | 486.0/5 | 8.41e-103 | 4.113573 +/- 0.000298 | -0.000651 +/- 0.000697 | |
| EJ-204 | even quadratic | 208.8/5 | 3.75e-43 | 4.119038 +/- 0.000443 | | -0.029929 +/- 0.001795 |
| EJ-204 | full quadratic | 207.8/4 | 7.74e-44 | 4.119039 +/- 0.000443 | -0.000677 +/- 0.000697 | -0.029933 +/- 0.001795 |
| EJ-230 | constant | 98.64/6 | 4.82e-19 | 4.102196 +/- 0.000286 | | |
| EJ-230 | linear | 93.38/5 | 1.31e-18 | 4.102202 +/- 0.000286 | +0.001579 +/- 0.000688 | |
| EJ-230 | even quadratic | 96.34/5 | 3.12e-19 | 4.101737 +/- 0.000417 | | +0.002650 +/- 0.001746 |
| EJ-230 | full quadratic | 91.05/4 | 7.88e-19 | 4.101741 +/- 0.000417 | +0.001583 +/- 0.000688 | +0.002664 +/- 0.001746 |

The poor goodness of fit of even the full quadratic is physical information: the seven-point structure is not described by one parabola. An explicitly even quartic was added as a diagnostic after observing this failure. It reduces the constant chi-square by 1897, 477 and 89.8 for EJ-200/204/230, with final probabilities 0.0177, 0.0398 and 0.0659. This confirms even, non-parabolic structure without introducing an odd detector term.

Because `Tsum=2*T0` event by event, its fitted intercepts and coefficients are exactly twice the T0 values, its errors are twice as large, and its chi-square, ndf and probabilities are identical. The independent Tsum fits and covariances were nevertheless executed and are recorded in `model_fits.csv`.

## Mirror symmetry

Differences are the mean at `+|x|` minus the mean at `-|x|`.

| Material | |x| [mm] | Delta T0 [ps] | z |
|---|---:|---:|---:|
| EJ-200 | 200 | -0.69 +/- 1.10 | -0.63 |
| EJ-200 | 500 | -0.45 +/- 1.15 | -0.39 |
| EJ-200 | 650 | +0.29 +/- 1.15 | +0.26 |
| EJ-204 | 200 | -1.72 +/- 1.05 | -1.63 |
| EJ-204 | 500 | +0.17 +/- 1.16 | +0.15 |
| EJ-204 | 650 | -1.03 +/- 1.20 | -0.86 |
| EJ-230 | 200 | +1.71 +/- 0.99 | +1.73 |
| EJ-230 | 500 | +2.27 +/- 1.13 | +2.01 |
| EJ-230 | 650 | +1.04 +/- 1.21 | +0.87 |

`Delta Tsum` is exactly twice `Delta T0` with twice the uncertainty, so its z values are identical. The isolated EJ-230 2.01-sigma point is not accompanied by a consistent sign and significance across all tests. The fitted odd term reaches 2.30 sigma for EJ-230 and less than one sigma for the other materials. This does not meet the Case C pattern.

## Direct propagation response g(d)

Distances follow the established coordinate convention:

```text
dL = 0.7 m + x
dR = 0.7 m - x
```

The fourteen LEFT/RIGHT measurements per material are stored as `propagation_graphs/g_first_pe_<material>_points`. Equal-distance LEFT and RIGHT values collapse statistically:

| Material | collapse chi2/ndf | probability |
|---|---:|---:|
| EJ-200 | 9.30/6 | 0.157 |
| EJ-204 | 3.71/6 | 0.715 |
| EJ-230 | 9.50/6 | 0.147 |
| EJ-230 mean-first-five | 1.86/6 | 0.932 |

The x=0 comparison is excluded from these aggregate probabilities because its LEFT and RIGHT measurements share events and their covariance is not stored. Individual values remain in `propagation_collapse.csv`.

### Required linear and quadratic fits

| Material/estimator | model | b0 [ns] | b1 [ns/m] | b2 [ps/m2] | chi2/ndf | Delta chi2 |
|---|---|---:|---:|---:|---:|---:|
| EJ-200 first PE | linear | 0.245601 +/- 0.000141 | 5.529105 +/- 0.000481 | | 11004.0/12 | |
| | quadratic | 0.236724 +/- 0.000175 | 5.664321 +/- 0.001654 | -124.84 +/- 1.46 | 3703.5/11 | 7300.5 |
| EJ-204 first PE | linear | 0.242451 +/- 0.000143 | 5.531029 +/- 0.000483 | | 5491.7/12 | |
| | quadratic | 0.236899 +/- 0.000177 | 5.611535 +/- 0.001577 | -76.23 +/- 1.42 | 2616.6/11 | 2875.1 |
| EJ-230 first PE | linear | 0.239652 +/- 0.000144 | 5.519243 +/- 0.000470 | | 2805.5/12 | |
| | quadratic | 0.236677 +/- 0.000178 | 5.559875 +/- 0.001486 | -39.10 +/- 1.36 | 1975.3/11 | 830.2 |
| EJ-230 mean first 5 | linear | 0.243075 +/- 0.000100 | 5.704556 +/- 0.000409 | | 17516.6/12 | |
| | quadratic | 0.237769 +/- 0.000125 | 5.789002 +/- 0.001282 | -83.84 +/- 1.21 | 12682.3/11 | 4834.3 |

The quadratic term is overwhelmingly nonzero in the combined propagation fit, but the quadratic model is itself rejected. AIC also favors the quadratic over the linear model by thousands of units; it does not rescue its absolute fit quality. Residual plots show systematic alternating/S-shaped patterns around both fits.

### Effective-velocity magnitude

The formally requested derivative of the quadratic model,

```text
v_eff(d) = 1000 / (b1 + 2*b2*d)  mm/ns,
```

gives the following endpoint values:

| Material | v_eff(50 mm) | v_eff(700 mm) | v_eff(1350 mm) |
|---|---:|---:|---:|
| EJ-200 | 176.93 | 182.16 | 187.71 mm/ns |
| EJ-204 | 178.45 | 181.66 | 184.99 mm/ns |
| EJ-230 | 179.99 | 181.65 | 183.34 mm/ns |

These are retained only as the requested quadratic diagnostic because their parent fits have `chi2/ndf=180--337`. A more direct finite-difference estimate from adjacent sampled distances shows the actual nonmonotonic structure:

- EJ-200: 172.63--183.57 mm/ns;
- EJ-204: 175.04--182.56 mm/ns;
- EJ-230: 177.26--182.76 mm/ns.

The low-distance interval 50--200 mm gives the smallest velocity in every material. EJ-230 then rises to 182.76 mm/ns around 500--700 mm and falls to 180.11 mm/ns over 1200--1350 mm. A single slope cannot represent these variations.

## Algebraic link to T0

For

```text
g2(d) = b0 + b1*d + b2*d^2
```

and half-length `h=0.7 m`,

```text
T0_pred(x) = 0.5 [g2(h+x) + g2(h-x)]
            = b0 + b1*h + b2*(h^2+x^2)
            = constant + b2*x^2.
```

Thus a valid quadratic propagation fit must produce the same x2 coefficient in T0. It does not: direct T0 coefficients are `-0.0718`, `-0.0299`, and `+0.00265 ns/m2`, whereas the combined propagation fits give `-0.1248`, `-0.0762`, and `-0.03910 ns/m2`. The EJ-230 signs even disagree. The prediction overlays preserve this mismatch. It is further evidence that the global quadratic approximation is incomplete rather than evidence for LEFT/RIGHT asymmetry.

## Mean-first-five

For EJ-230, mean-first-five remains left-right symmetric: all mirrored differences are below 0.69 sigma, and the equal-distance collapse probability is 0.932. Its mean position dependence is larger:

| Estimator | endpoint average minus center | peak-to-peak | constant chi2/ndf | even-quadratic coefficient |
|---|---:|---:|---:|---:|
| first PE | +0.19 ps | 8.35 ps | 98.64/6 | +2.65 +/- 1.75 ps/m2 |
| mean first 5 | +2.35 ps | 17.48 ps | 637.58/6 | +13.26 +/- 1.62 ps/m2 |

Mean-first-five therefore improves resolution but worsens calibration uniformity of the mean. These are separate properties and should not be combined into one claim.

## Light-yield dependence

All 21 event-level profiles are stored under `light_profiles/` in `tsum_veff.root`. Pearson coefficients are descriptive:

- `rho(T0,Npe_END)` ranges from -0.396 to -0.252;
- `rho(tL,Npe_L)` ranges from -0.476 to -0.225;
- `rho(tR,Npe_R)` ranges from -0.470 to -0.228;
- for EJ-230, `rho(T0_mean5,Npe_END)` ranges from -0.540 to -0.385.

Events with more detected light tend to have earlier first-arrival estimates. This is a **photon-statistics / first-arrival bias**, or light-yield-dependent timing bias. It is not electronic time walk: SPTR, pulse formation, threshold and TDC effects are absent.

## Spectral and path information

The original `sipm_hits` tree contains detected `energy_eV`, `wl_nm`, PDE and detector coordinates. A ROOT pass over all 21 production files combined equal-distance LEFT/RIGHT spectra and identified the wavelength of the first END hit in every event.

The integrated detected END spectrum is essentially stable from 50 to 1350 mm:

| Material | mean shift, long-short | median shift |
|---|---:|---:|
| EJ-200 | -0.015 nm | -0.009 nm |
| EJ-204 | +0.012 nm | +0.016 nm |
| EJ-230 | -0.071 nm | -0.024 nm |

The first-arrival subset changes more strongly:

| Material | first-arrival mean shift | median shift |
|---|---:|---:|
| EJ-200 | +10.58 nm | +8.00 nm |
| EJ-204 | -2.82 nm | +1.88 nm |
| EJ-230 | -10.75 nm | -0.95 nm |

This demonstrates wavelength-dependent selection among the first detected photons. It does not establish whether the cause is group-velocity dispersion, attenuation, PDE weighting, emission-order statistics, or their combination.

Accumulated path length, number of boundary encounters, creation wavelength and a detected-hit track identifier are unavailable. The `first_bar_encounters` tree records only the first encounter and cannot be joined to `sipm_hits`, which has no `track_id`. Consequently:

`MECHANISM_NOT_DIRECTLY_RESOLVABLE_FROM_CURRENT_TTREE`.

## ROOT objects and files

Primary numerical artifact: `/home/rrios/ej200/analysis/tsum_veff_20260914/sources/tsum_veff.root`.

- `derived_events`: 210,000 event-level rows;
- `distributions/T0_<cell_id>`: full T0 histogram for each cell;
- `light_profiles/T0_vs_NpeEND_<cell_id>` and side-specific profiles;
- `position_graphs/t0_vs_x_<material>` and `tsum_vs_x_<material>`;
- `position_fits/<graph>_<model>` plus `<fit>_covariance`;
- `propagation_graphs/g_first_pe_<material>` and side-tagged point trees;
- `propagation_fits/<graph>_linear|quadratic` and covariance matrices.

Additional ROOT/CSV artifacts preserve the spectral histograms, even-quartic diagnostic, equal-distance goodness tests, finite-difference velocities, fit tables, mirror tests and decision summary.

Every figure under `/home/rrios/ej200/analysis/tsum_veff_20260914/figures/` has a same-name ROOT macro in `macros/`, a `.root` object file, `.pdf`, and `.meta.json` provenance sidecar. Required figures include all three materials for T0, Tsum, propagation fits, residuals, quadratic velocity, mirror tests and seven-panel T0 distributions, plus the EJ-230 mean-first-five comparison.

## Reproduction

These commands read existing ROOT files and do not invoke Geant4:

```bash
cd /home/rrios/ej200/analysis/tsum_veff_20260914
root -l -b -q 'macros/analyze_tsum_veff.C+'
root -l -b -q 'macros/extended_models.C+'
root -l -b -q 'macros/summarize_decision.C+'
root -l -b -q 'macros/collapse_goodness.C+'
root -l -b -q 'macros/spectral_check.C+'
```

The per-figure wrapper macros can then be executed with `root -l -b -q macros/<figure>.C`.

## Statistics gate

`HIGH_STAT_REQUIRED` is **not** declared. Existing statistical precision is already much smaller than the 8--32 ps position effects and the structured propagation residuals. More events would increase the rejection significance without making the linear or quadratic forms adequate. No new simulation was launched, and v9 was not modified.
