# EXEC_14/14B fit failures and quality flags

Updated: 2026-06-13

## Landau/Moyal spectrum fits

Source: `results_ej230_analysis/exec10_landau_mpv.csv`.

- 257 fits were requested and 257 converged (`fit_status=ok`); failed fits: **0**.
- 345 low-light spectra are explicitly `not-requested`, not failed fits.
- 35 converged fits place `sigma_gauss` at the 0.1 Npe identifiability bound.
- 36 converged fits have `sigma_gauss_error > 10 sigma_gauss`.
- 55 converged fits have chi2/ndf above 3.
- 31 high-light fits (`mean > 50 Npe`) have chi2/ndf above 5, consistent with
  the documented limitation of the Moyal proxy in the hard Landau tail.

These quality flags remain in the CSV and are not converted into fit failures.

## Time-to-threshold Gaussian fits

Source: `results_ej230_analysis/csv/exec13_tN_summary.csv`.

Eight requested Gaussian widths are undefined because the selected threshold/group
has insufficient support at the edge positions:

| x [mm] | Group | N | Status |
|---:|---|---:|---|
| -690 | `end_near` | 4 | degenerate / no finite sigma fit |
| -690 | `end_far` | 20 | degenerate / no finite sigma fit |
| -650 | `end_near` | 4 | degenerate / no finite sigma fit |
| -650 | `end_far` | 20 | degenerate / no finite sigma fit |
| +650 | `end_near` | 4 | degenerate / no finite sigma fit |
| +650 | `end_far` | 20 | degenerate / no finite sigma fit |
| +690 | `end_near` | 4 | degenerate / no finite sigma fit |
| +690 | `end_far` | 20 | degenerate / no finite sigma fit |

The report does not replace these missing fits with zero or with a value from another
material. Finite RMS values remain available in the source CSV.

## Special-control analyses

EXEC_14B ran authentic EJ-230 window-track and End-only controls. Their CSV and PNG
products replace the three previously copied EJ-204 figures. The EndTop/End-only
timing gate is allowed to return its documented scientific non-zero status when
EndTop is narrower; this is a physical result, not a fit or execution failure.
