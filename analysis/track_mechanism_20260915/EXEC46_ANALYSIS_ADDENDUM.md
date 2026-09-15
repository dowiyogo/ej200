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
