# EXEC_46 I1-I5: BC-408/BC-404 validation and production launch

Date: 2026-09-16.

## I1: Construction of the EJ-204 MPT

OPSC-101 uses the Huggins et al. BC-404 parameterization from arXiv:2608.13710:

| coefficient | value |
|---|---:|
| A | 1.578 |
| B | 0.818 |
| C [nm^-1] | 0.00729 |

The implemented form is `n(lambda_nm)=A+B*exp(-C*lambda_nm)`. The RINDEX table
has 148 rows: 146 wavelength points from 370 through 660 nm at 2 nm spacing,
plus the 200 nm and 800 nm clamp rows. Thus the actual count is 146 points in
the 370-660 nm interval, not the previously quoted 151 points. The previously
quoted 148-row count applies to the complete RINDEX file; it is not the count
of points in that interval.

The 200 nm row holds the 370 nm index through the lower clamp interval, and the
800 nm row holds the 660 nm index through the upper clamp interval. ABSLENGTH
has 4 rows: 26.58 mm through 372 nm, interpolation to 3800 mm at 439 nm, and
3800 mm through 800 nm. OPSC-101 expresses ABSLENGTH in cm, as OPSC-100 does;
the stored values are therefore 2.658 and 380 cm.

Source: `analysis/track_mechanism_20260915/prepare_campaign.py` and the
OPSC-100/OPSC-101 `rIndex.txt` and `absLength.txt` files used by the campaign.

## I2: EJ204_xm650 validation

The authorized validation used 10,000 events, four workers, EndTop geometry,
zero SPTR, and the unchanged OPSC-101 emission, yield, and timing
configuration. The measured metrics are:

| scenario | first Cherenkov [%] | primary caustic q95-q05 [deg] | first-Cher local velocity [mm/ns] | first-scint local velocity [mm/ns] | MPT v_group(408 nm) [mm/ns] |
|---|---:|---:|---:|---:|---:|
| constant-n baseline | 62.90 +/- 0.48 | 2.435 +/- 0.143 | 150.017 +/- 0.101 | 181.255 +/- 0.111 | 189.742 |
| BC-404, 3800 mm | 72.88 +/- 0.44 | 6.026 +/- 0.205 | 149.356 +/- 0.093 | 173.078 +/- 0.096 | 171.893 |

The corrected-minus-baseline Cherenkov fraction is 9.98 +/- 0.66 percentage
points. The I2 gate is **PASS**. The corrected metrics are read from
`i2_metrics.csv`; the gate and difference are recorded in
`i2_bc404_validation/analysis_summary.json`.

The direction agrees with the EJ-200 result: first Cherenkov fraction changes
from 74.10% to 82.02%. Source: the EXEC46 BC-408 validation record referenced
by the Step 5 input chain.

## I3: Combined production grid

The production directory is
`/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800`.

| material | optical status |
|---|---|
| EJ-200 | `CORRECTED_BC408_3800_VALIDATED_F4` |
| EJ-204 | `CORRECTED_BC404_3800_VALIDATED_I2` |
| EJ-230 | `UNCORRECTED_CONSTANT_RINDEX_NO_MEASURED_ANALOG` |

BC-420 was not measured by Huggins et al.; BC-422 is not a valid substitute.
EJ-230 therefore remains on the constant-RINDEX model.

The four MPT hashes in `campaign.json`, independently checked against the
source files, are:

| material | file | SHA-256 |
|---|---|---|
| EJ-200 | rIndex.txt | `15d1f8cf5a62effd9a0f2f9bd1edaeb5164b94c2035a11cfe6a878d51680f6ed` |
| EJ-200 | absLength.txt | `82191c023ed343d9b0f0c6de09ec3f1be12627e4eed6219d8e048529d7755d3a` |
| EJ-204 | rIndex.txt | `f1c77ee162c767cd23be608e0f8081174263352e40a0e1f57b8ca97ccef8636f` |
| EJ-204 | absLength.txt | `9473735344b5bd883a748df0223129b692edb659a7bd9c18f8d32be450749660` |

The recorded DRY_RUN passed with 21 readable macros, 21 pending cells, zero
outputs at dry-run time, and diagnostics OFF. It recorded no timeout.

`prepare_campaign.py` now requires the I2 summary to exist and to declare
`status == "PASS"` before accepting a corrected SSLG4 source and creating the
combined grid. This gate is in addition to the MPT hash checks.

## I4: Detached launch state

The authorized command is:

```bash
python3 analysis/sigma_t/orchestration/detached_grid.py launch \
  --directory /home/rrios/exec46_20260916/full_grid_bc408_bc404_3800
```

Current state from the files on disk: **LAUNCHED - output ROOT files are
present in the campaign cell directories**. This conclusion uses only the
presence of output files; no `status` command was run for this report. The
DRY_RUN record predates those files and retains `EXECUTED=0` and 21 pending
cells, so it documents preparation state rather than the later disk state.

## I5: Step 5 closeout

`beta_within` is the slope of T0 against Npe measured within a fixed position,
with position fixed effects.

`beta_between` is the slope obtained by regressing the means <T0> against the
means <Npe> between the seven positions.

The `between curve` is the fit of <T0> against <Npe> over those seven points;
the tested specification was linear versus pol2.

The registered remnant is:

| material | remnant [ps/m^2] | beta_within significance [sigma] |
|---|---:|---:|
| EJ-200 | +123.10 | 16.844 |
| EJ-204 | +227.05 | 20.407 |
| EJ-230 | +197.23 | 25.723 |

The 7.183 +/- 0.515 ps residual is identified as an artifact of specifying the
between curve linearly. The nested seven-position tests give:

| material | F(1,4) | p |
|---|---:|---:|
| EJ-200 | 140.204538 | 0.000291 |
| EJ-204 | 113.117874 | 0.000443 |
| EJ-230 | 46.449773 | 0.002423 |

The leave-one-out comparison is:

| material | linear RMSE [ps] | pol2 RMSE [ps] | reduction |
|---|---:|---:|---:|
| EJ-200 | 3.5300427 | 0.7844389 | 77.8% |
| EJ-204 | 4.2141425 | 0.9955193 | 76.4% |
| EJ-230 | 3.7888589 | 1.4929028 | 60.6% |

Together with the prior refutation of `beta_within` as the between-position
response at 16.844, 20.407, and 25.723 sigma, these results support the
conclusion that the nonlinearity of T0(x) is compatible with being entirely an
Npe response once it is modeled with curvature and estimated between
positions.

The +123.10 / +227.05 / +197.23 ps/m^2 remnant is explained by two compound
specification errors: transferring `beta_within` to a between-position
response and forcing a linear between curve. This does not prove absence of a
physical mechanism; it only shows that the seven-position design does not
resolve one.

Sources for I5: `analysis/track_mechanism_20260915/step5/REPORT_CHAINRULE_IDENTIFICATION_20260916.md`,
`h1_even_f_test.csv`, `h1_loo_summary.csv`, and `within_between_summary.csv`.

## Open items

No additional quantitative item was available in the specified sources.
