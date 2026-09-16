# EXEC_46 analysis addendum approved after D1–D4

This addendum records the analysis contract accepted on 2026-09-15. It does
not alter simulation physics, geometry, ROOT inputs, or the earlier failed-gate
record.

## Revised Step 1 decay gate

For every cell, read `SCINTILLATIONTIMECONSTANT1` independently from the active
OPSC macro. Select `source_type == 1` and
`abs(x_creation_mm - gun_x_mm) < 0.1 mm`. The subtraction by `gun_x_mm` makes
the local coordinate meaningful at all seven gun positions; at x=0 it is
identical to the originally diagnosed `abs(x_creation_mm) < 0.1 mm` selection.

Evaluate the untruncated mean excess at cuts 3, 4, and 5 times the configured
decay constant. Each cut passes when its mean excess differs from the configured
decay by at most 3%. The 3% tolerance is an analysis assumption. Require at
least 200 contributing events at every cut. The six-cut result at 3, 4, 5, 8,
12, and 20 decay constants for the complete source-type-1 population remains
an informational diagnostic and is never a gate.

## Step 4 order statistics

Do not use the 3.098019 ns mean of the END emission pool, or any other full-pool
mean, as an exponential benchmark. Predict the first photon from the empirical
`t_creation_ns` distribution.

Compare the measured creation-time order statistic with N^-1, N^-1/2, and a
free power-law exponent. Read rise and decay parameters independently from each
material's active configuration.

An independent check of `G4Scintillation::sample_time` confirms the N^-1/2
exponent and finds that the prefactor written for a standard biexponential is
not the exact prefactor of the Geant4 11.4 rejection sampler. Geant4 samples a
decay exponential and accepts it with probability `1-exp(-t/tau_r)`, yielding
(`G4Scintillation.cc:557-570` in the installed Geant4 11.4.0 source)

```text
f_G4(t) = (tau_r+tau_d)/tau_d^2
          * exp(-t/tau_d) * [1-exp(-t/tau_r)]

f_G4(t) -> [(tau_r+tau_d)/(tau_d^2*tau_r)] * t

E[min_N] -> tau_d * sqrt(pi*tau_r/[2*(tau_r+tau_d)]) * N^(-1/2)
```

The standard biexponential expression supplied in the approval,
`sqrt(2*tau_r*tau_d)*Gamma(3/2)*N^(-1/2)`, must still be shown as a distinct
reference curve. It must not be labeled as the prediction of the actual
Geant4 sampler. The numerical verification and all three material-dependent
prefactors are stored in `exec46_order_stat_asymptotic.{csv,json}`.

Analyze scintillation and Cherenkov order statistics separately. Report the
fraction of events whose first photon is Cherenkov versus Npe and gun position.
Treat mechanisms 4 and 5 as a coupled mixture in Step 8 rather than additive,
independent contributions.

## Step 6 creation-position covariate

Include the variance of creation position together with `path_length_mm` in the
mechanism-1 analysis. Preserve the measured displaced population as physical
spatial dispersion; do not collapse it into a purely temporal correction.

## Additions approved before Step 2

The configured RINDEX tables are constant at 1.58 over 200--800 nm. Therefore
`dn/dlambda` is zero and the proposed group-velocity spectral correction is
identically zero by construction. Remove it from Step 6 and record the missing
dispersion as a model limitation. Investigate spectral selection only through
the configured PDE and ABSLENGTH properties. Step 2 must report the latter
directly from each material configuration. If ABSLENGTH is itself constant in
wavelength, say so; do not infer spectral filtering that the model cannot
produce.

Step 2 also compares the robust single-left-END boundary
`sigma_IQR(tL)/sigma_IQR(T0) = 0.5` with the first-photon Cherenkov transition,
using the three materials without forcing an identification. It contrasts
`exit_angle_deg` and `n_boundary_encounters` for scintillation and Cherenkov,
separately for first and deterministic random photons from the same event and
END face. The recorded exit angle is relative to the SiPM normal at detection;
it is not a per-boundary incidence-angle history.

For Step 4, N in the first-scintillation-photon order statistic is the number
of detected scintillation photons at the relevant END. Use the separately
counted scintillation population; never substitute total `Npe_END`.

## Additions approved before Step 3

The A2 boundary comparison is withdrawn.  The correlation between the first
Cherenkov fraction and `sigma_IQR(tL)/sigma_IQR(T0)` changes sign between
`|x|=500` and 650 mm, so neither linear crossing is a measured boundary.  Both
are sparse-grid interpolation artifacts.  The 21 cells provide only seven
END distances, 50, 200, 500, 700, 900, 1200, and 1350 mm; the observed change
of regime lies wholly inside the unmeasured 50--200 mm interval.  Step 3 must
record that gap and must not propose new simulation.

Reparameterize the Step 3 consistency presentation in physical END distance
`d`: for each of the seven distances, combine the left-END realization with
the right-END realization from the mirror cell, while retaining both as
independent consistency controls.  The primary definition of `g(d)` remains
microscopic `d_direct`, as required by the EXEC_46 prompt; nominal distance is
used only to pair the two realizations and display the seven sampled scales.

Before applying the Cherenkov geometry test, verify the gun direction from the
production macro and source.  For a primary perpendicular to the bar axis,
the Cherenkov cone has an axial-angle lower edge equal to
`asin(1/n)=39.27 deg` and is trapped at the large faces for `n > sqrt(2)`.
Test the fine-binned lower edge and caustic for the first Cherenkov photon at
the near END in all three materials.  `source_type == 2` identifies the
creator process but does not encode the creator track or parentage; therefore
separate the exact primary-cone prediction from any observed population of
secondary-particle Cherenkov photons.

The identity between the axial cone edge and the critical angle is exact in
the beta=1 limit.  The configured 1 GeV kinetic-energy muon has beta=0.995424,
so Step 3 must retain the requested beta=1 prediction and also report the
parameter-free finite-beta correction from the gun configuration.

Build `g(d)` separately for scintillation and Cherenkov.  For the Cherenkov
cone edge, compare the fitted axial velocity with the parameter-free
`c*sqrt(1-1/n^2)/n = 146.9 mm/ns` prediction and with the configured group
velocity, `c/n = 189.742 mm/ns`.  Report disagreement explicitly rather than
retuning either value.

## Additions approved before Step 4

Compare the first Cherenkov photon with the first scintillation photon as two
order statistics.  Count the detected Cherenkov population independently per
event and END, and test whether the source-specific minimum approaches the
cone-edge velocity as that population grows.  Stratify the folded axial exit
angle by the seven measured END distances and test its timing penalty directly.
The upper angular window is a timing selection, not a geometric cone boundary.

Correct the Step 3 caustic wording: the finite-beta edge lies inside a partly
filled bin and the mode is in the following bin; do not quote the lower-bin-edge
difference as the agreement.  Treat isolated bins at 41--43 degrees with about
one count as Poisson fluctuations.  The equality of the beta=1 axial cone edge
and the END critical angle follows exactly from
`acos(sin(theta_C)) = asin(1/n)`; the finite-beta correction remains a separate
configured-gun effect.

For scintillation only, fit the first-emission order statistic against the
detected scintillation population, testing `N^-1`, `N^-1/2`, and a free
exponent.  Compare the inferred effective population with photons having zero
to two boundary encounters and test material independence.  Do not apply an
emission-lifetime order-statistics law to Cherenkov photons.

Recover the Cherenkov/width correlation only at measured fixed points: report
the six mirror values independently at END distances 50 and 200 mm without
interpolating.  Decompose the six rejected near-END Gaussian timestamp fits by
source, quantify the angular caustic, and replace the rejected Gaussian width
there with an explicitly defined source-mixture estimator.

## Corrections approved before Step 5

Recondition the angular-window test on Cherenkov multiplicity at fixed distance
for d <= 500 mm. The former distance-only conclusion is invalid where the low
quintile has N_C=1 and no order statistic can select the cone edge.

For the scintillation order statistic, subtract the parent-muon transit estimate
`(5 mm - z_creation_mm)/c` before selecting the minimum and refit the
`N^-1/2` law. The configured gun points along -z and the +z bar entry is 5 mm.
Retain both raw and corrected results.

Treat `sigma_mixture` as the intrinsic zero-jitter optical limit rather than
detector performance. Compare the source separation with that limit and retain
the angle-time correlation as a direct cone signature.

Step 5 first compares the common within-cell slope with the seven-point
between-position slope. The original majority-artifact gate is superseded by
the following E1--E5 revision (2026-09-16).

## Step 5 E1--E5 revision

The between-position regression fits the response being explained. Its roughly
97% reduction of the former quadratic remnant is reabsorption into a fitted
parameter, not an independent mechanism or explanation. The status is
`CHAIN_RULE_WITHIN_SLOPE_REFUTED`. Step 5.2 is cancelled: cell centering fixes
alpha to the observed cell mean and supplies no independent information.

Finish 5.3--5.5, with the descriptive residual near |x|=500 mm as the target:
retain both mirrors, cell SEM and joint uncertainty; compare pol2 and quantile
profiles, the source-mixture finite-difference identity, and the two-count
chain. Test rather than assume exact material universality of beta_within.
Do not reinterpret the between-fit residual as a causally corrected observable.
The remaining approval gate is before Step 6. No simulation or push.
