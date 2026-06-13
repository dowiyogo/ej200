# EXEC_14/14D/14E fit failures and quality flags

Updated: 2026-06-13

## Landau/Moyal spectrum fits

Source: `results_ej230_analysis/exec10_landau_mpv.csv`.

- 257 fits were requested and 257 converged (`fit_status=ok`); failed fits: **0**.
- 345 low-light spectra are explicitly `not-requested`, not failed fits.
- 35 converged fits place `sigma_gauss` at or numerically near the 0.1 Npe
  identifiability bound: 33 at the bound within `1e-6` and 2 near it.
- 36 converged fits have `sigma_gauss_error > 10 sigma_gauss`.
- 55 converged fits have chi2/ndf above 3.
- 31 high-light fits (`mean > 50 Npe`) have chi2/ndf above 5, consistent with
  the documented limitation of the Moyal proxy in the hard Landau tail.

These quality flags remain in the CSV and are not converted into fit failures.

## Time-to-threshold Gaussian fits (EXEC_14E final)

Sources:

- `results_ej230_analysis/csv/exec13_tN_summary.csv`
- `results_ej230_analysis/csv/exec14e_fixed_tN_summary.csv`

The adaptive-bin fitted `t4` metric approved in EXEC_14D remains unchanged in
the summary slides:

- Four near-End t4 cores (`x = -690, -650, +650, +690 mm`) now converge with
  the documented adaptive-bin fit. They are explicitly flagged
  `resolution_limited`; the fitted values are not replaced by RMS.
- Final missing fitted widths: **0**. Final genuine-starvation points: **0**.

The nearest-Top fitted t4 range is `6.648--9.594 ps`; the near-End fitted t4
range is `0.562--127.719 ps`. Values narrower than the reporting bin-width hook
remain annotated as resolution-limited together with their mean Npe.

For the restored main per-position display, every panel uses a fixed 25 ps bin
width and explicit `N=4`/`N=20` thresholds. The far-End `t20` histogram is
never hidden because of low reach. All 30 non-central far-End fixed `t20` core
fits converge with finite covariance, including the four edge positions:

- `x=-690`: `R(20)=10.85%`, `sigma_fit=1350.6 +/- 132.2 ps`
- `x=-650`: `R(20)=16.35%`, `sigma_fit=1409.8 +/- 121.6 ps`
- `x=+650`: `R(20)=16.90%`, `sigma_fit=1363.4 +/- 118.9 ps`
- `x=+690`: `R(20)=11.90%`, `sigma_fit=1429.2 +/- 146.2 ps`

Therefore no current far-End panel needs the documented “not fitable; see
Backup” fallback. The low-reach panels still show their real histogram and
excluded fraction. The adaptive `N_eff` analysis is preserved in appended
Backup frames rather than used in the main per-position figures.

## Special-control analyses

EXEC_14B ran authentic EJ-230 window-track and End-only controls. Their CSV and PNG
products replace the three previously copied EJ-204 figures. EXEC_14D verifies
that EndTop is wider for all six side/position comparisons after matching the
estimator and bootstrapping End-only to 2000 events. This is a physical result,
not a fit or execution failure.

The final EndTop/End-only fitted-width ratio spans `1.223218--1.963382`; no
ratio was replaced, clipped, or interpreted as an execution failure.
