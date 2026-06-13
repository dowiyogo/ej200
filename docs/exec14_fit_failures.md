# EXEC_14/14D fit failures and quality flags

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

## Time-to-threshold Gaussian fits (EXEC_14D final)

Source: `results_ej230_analysis/csv/exec13_tN_summary.csv`.

The eight edge-position degeneracies recorded by EXEC_14B are resolved:

- Four near-End t4 cores (`x = -690, -650, +650, +690 mm`) now converge with
  the documented adaptive-bin fit. They are explicitly flagged
  `resolution_limited`; the fitted values are not replaced by RMS.
- Four far-End t20 requests at the same positions now use
  `HOOK_ADAPTIVE_TN`. Their data-selected thresholds are respectively
  `N_eff = 7, 8, 8, 7`, all with event reach above the configured 95% minimum.
- Across all 31 positions, 17 far-End thresholds are reduced, 14 retain the
  nominal `N=20`, and none meet the genuine-starvation condition.
- Final missing fitted widths: **0**. Final genuine-starvation points: **0**.

The nearest-Top fitted t4 range is `6.648--9.594 ps`; the near-End fitted t4
range is `0.562--127.719 ps`. Values narrower than the reporting bin-width hook
remain annotated as resolution-limited together with their mean Npe.

## Special-control analyses

EXEC_14B ran authentic EJ-230 window-track and End-only controls. Their CSV and PNG
products replace the three previously copied EJ-204 figures. The EndTop/End-only
timing gate is allowed to return its documented scientific non-zero status when
EndTop is narrower; this is a physical result, not a fit or execution failure.

The final EndTop/End-only fitted-width ratio spans `1.223218--1.963382`; no
ratio was replaced, clipped, or interpreted as an execution failure.
