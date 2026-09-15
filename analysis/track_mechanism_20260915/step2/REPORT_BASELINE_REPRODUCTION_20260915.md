# EXEC_46 Step 2 — bit-level baseline reproduction and A1--A3 diagnostics

Date: 2026-09-15

## Verdict

**PASS.** The reconstructed EXEC_46 baseline reproduces the registered END-only first-photon result under the declared tolerances.
No simulation was run. The source ROOT files were opened read-only.

The point tolerance is 0.01 ps, set by the precision of the registered table. The quadratic-coefficient tolerance is five times the registered statistical error.

## Baseline point-by-point reproduction

| Material | |x| [mm] | observed old/new [ps] | Npe prediction old/new [ps] |
|---|---:|---:|---:|
| EJ-200 | 200 | -1.354 / -1.353540 | -6.980 / -6.979838 |
| EJ-200 | 500 | -6.264 / -6.263928 | -44.246 / -44.245989 |
| EJ-200 | 650 | -31.922 / -31.922475 | -82.784 / -82.784413 |
| EJ-204 | 200 | +0.306 / +0.306167 | -10.061 / -10.060705 |
| EJ-204 | 500 | +1.936 / +1.936152 | -60.121 / -60.120943 |
| EJ-204 | 650 | -14.307 / -14.307034 | -109.642 / -109.642155 |
| EJ-230 | 200 | +1.679 / +1.678752 | -7.352 / -7.351532 |
| EJ-230 | 500 | +6.889 / +6.888857 | -46.831 / -46.830878 |
| EJ-230 | 650 | +0.194 / +0.193881 | -82.499 / -82.499382 |

The per-cell 40-bin TProfile `pol1` fits reproduce the previous quality result: 12/21 have chi2/ndf > 5 and are flagged as unreliable. The chain-rule curve is reproduced as registered but is not promoted to a valid physical model by this check.

## Position fits

An asterisk marks chi2/ndf > 5; those fits are descriptive projections, not adequate models.

| Material | observed a2 [ps/m2] | predicted a2 [ps/m2] | observed quadratic chi2/ndf | observed quartic a2/a4 | quartic chi2/ndf |
|---|---:|---:|---:|---:|---:|
| EJ-200 | -71.8094 +/- 1.7976 | -194.9126 | 313.07/5 = 62.61* | +54.53 / -296.84 ps | 11.95/4 = 2.99 |
| EJ-204 | -29.9288 +/- 1.7948 | -256.9820 | 208.76/5 = 41.75* | +72.02 / -242.58 ps | 10.03/4 = 2.51 |
| EJ-230 | +2.6503 +/- 1.7464 | -194.5761 | 96.33/5 = 19.27* | +67.97 / -156.84 ps | 8.81/4 = 2.20 |

## Clock comparison

`T0=(tL+tR)/2` is exact in all 210,000 events for each clock. The maximum event-level differences are:

| Quantity | max abs difference [ns] |
|---|---:|
| tL(time)-tL(detection) | 0 |
| tR(time)-tR(detection) | 0 |
| T0(time)-T0(detection) | 0 |
| T0(time)-formula | 0 |
| T0(detection)-formula | 0 |

## A1 — optical-property limitation

| Material | RINDEX, 200--800 nm | ABSLENGTH, 200--800 nm | spectral dependence |
|---|---:|---:|---|
| EJ-200 | 1.58 | 3800 mm | RINDEX constant; ABSLENGTH constant |
| EJ-204 | 1.58 | 1600 mm | RINDEX constant; ABSLENGTH constant |
| EJ-230 | 1.58 | 1200 mm | RINDEX constant; ABSLENGTH constant |

The group-velocity spectral correction is identically zero by construction. ABSLENGTH also has no wavelength dependence in these tables, so this model cannot produce wavelength-selective bulk attenuation. The wavelength-dependent PDE remains active and is the only modeled spectral selection among these proposed paths.

## A2 — Cherenkov and robust-width boundaries

Both boundaries use linear interpolation inside the sparse -650 to -500 mm bracket: ratio=0.5 and first-left Cherenkov fraction=0.5.

| Material | ratio boundary [mm] | Cherenkov boundary [mm] | Cher - ratio [mm] |
|---|---:|---:|---:|
| EJ-200 | -560.69 | -594.94 | -34.25 |
| EJ-204 | -553.50 | -615.11 | -61.61 |
| EJ-230 | -549.61 | -644.05 | -94.44 |

Across the three interpolated boundaries, Pearson r = -0.9631. Across all 21 left-END material-position points, r(ratio, first-Cherenkov fraction) = -0.7751.
**BOUNDARY_IDENTITY_NOT_ESTABLISHED.** The material ordering of the two interpolated crossings is opposite, so the two boundaries cannot be identified as the same measured transition.
With only three materials and one 150 mm transition interval, the boundary correlation is descriptive. No causal identity is forced from the point correlation alone.

## A3 — Cherenkov guiding diagnostic

For the configured n=1.58 and beta approximately one, theta_C=50.73 deg relative to the primary direction and theta_critical=39.27 deg relative to a boundary normal.

The stored `exit_angle_deg` is relative to the SiPM normal at final detection. It is not the incidence angle at the large bar face and cannot directly prove the proposed TIR-cone inequality. The angle and boundary-count contrasts are therefore indirect tests; a direct test requires per-boundary angle history.

All-END summary (mean with median in parentheses):

| Material | selection | source | exit angle [deg] | boundary encounters | N |
|---|---|---|---:|---:|---:|
| EJ-200 | first | scintillation | 11.24 (10.48) | 10.70 (8) | 119118 |
| EJ-200 | first | Cherenkov | 34.85 (39.52) | 5.34 (4) | 20882 |
| EJ-200 | random | scintillation | 41.80 (38.18) | 46.55 (30) | 135574 |
| EJ-200 | random | Cherenkov | 56.19 (44.90) | 80.40 (74) | 4426 |
| EJ-204 | first | scintillation | 11.17 (10.40) | 10.48 (8) | 122310 |
| EJ-204 | first | Cherenkov | 34.54 (39.52) | 5.19 (4) | 17690 |
| EJ-204 | random | scintillation | 40.47 (36.71) | 42.03 (28) | 136077 |
| EJ-204 | random | Cherenkov | 54.33 (44.42) | 71.39 (66) | 3923 |
| EJ-230 | first | scintillation | 11.02 (10.26) | 10.18 (8) | 125137 |
| EJ-230 | first | Cherenkov | 33.07 (39.50) | 5.02 (4) | 14863 |
| EJ-230 | random | scintillation | 39.76 (35.89) | 40.15 (27) | 135863 |
| EJ-230 | random | Cherenkov | 55.55 (44.48) | 69.18 (65) | 4137 |

Near-END summary at x=-650 left plus x=+650 right:

| Material | selection | source | exit angle [deg] | boundary encounters | N |
|---|---|---|---:|---:|---:|
| EJ-200 | first | scintillation | 17.69 (17.17) | 1.03 (1) | 5234 |
| EJ-200 | first | Cherenkov | 44.36 (39.62) | 3.62 (4) | 14766 |
| EJ-200 | random | scintillation | 55.83 (52.79) | 9.10 (5) | 19362 |
| EJ-200 | random | Cherenkov | 63.95 (57.72) | 11.94 (7) | 638 |
| EJ-204 | first | scintillation | 17.70 (17.40) | 1.01 (1) | 7368 |
| EJ-204 | first | Cherenkov | 43.77 (39.62) | 3.56 (4) | 12632 |
| EJ-204 | random | scintillation | 55.55 (52.28) | 8.31 (5) | 19393 |
| EJ-204 | random | Cherenkov | 62.10 (56.47) | 9.97 (7) | 607 |
| EJ-230 | first | scintillation | 17.73 (17.39) | 1.01 (1) | 9590 |
| EJ-230 | first | Cherenkov | 42.52 (39.60) | 3.45 (4) | 10410 |
| EJ-230 | random | scintillation | 55.26 (52.39) | 7.98 (5) | 19365 |
| EJ-230 | random | Cherenkov | 64.04 (58.36) | 9.92 (7) | 635 |

At the near END, first Cherenkov photons have a narrow median exit angle near 39.6 degrees and a median of four prior boundary encounters in every material. First scintillation photons have median exit angles near 17 degrees and one boundary encounter. The random controls are broader and have medians of seven boundaries for Cherenkov and five for scintillation. This supports a distinct directional, promptly selected Cherenkov population, but the stored final-exit angle does not directly establish TIR at the large faces.

## Reproducibility and artifacts

The derived tree contains 21 cells x 10,000 events and has SHA-256 `816c302c3ed37f7eb29855bcad7e734a1b6a77238f60a30361635b68377e341d`. The paired random control uses `minimum splitmix64 priority per event and END face` with seed `0x46a2c0de5eed1234`.

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/build_step2_derived.py --processes 4
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step2.py
```

- `exec46_derived_events.root` and `.meta.json`: 210,000 event rows and exact production provenance.
  It stores `npe_scint_left/right` separately; Step 4 must use those per-END scintillation counts as N rather than total Npe_END.
- `baseline_reproduction.{pdf,root,csv,meta.json}`.
- `cherenkov_boundary.{pdf,root,csv,meta.json}`.
- `cherenkov_guiding_diagnostics.{pdf,root,csv,meta.json}`.
- `baseline_cells.csv`, `baseline_fits.csv`, `clock_comparison.csv`, and `material_optical_properties.csv`.

No push, simulation, merge, or deck edit was performed. Step 3 requires a new checkpoint approval.
