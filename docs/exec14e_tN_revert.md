# EXEC_14E fixed-threshold t_N display restoration

Updated: 2026-06-13

## Scope

EXEC_14E restores the main per-position `t_N` figures to the historical
fixed-threshold display while retaining the EXEC_14D adaptive-threshold result
as appended Backup material.

- Main figures: six panels in two rows, with `t4` (`N=4`) above and `t20`
  (`N=20`) below.
- Columns: nearest Top, End-L SUM4, and End-R SUM4.
- Display bin width: fixed `25 ps` in every panel.
- Summary slide: the approved adaptive-bin fitted `sigma(t4)` metric remains on
  the left; the right panel now shows fitted fixed-threshold far-End `t20`.
- Adaptive outputs: preserved as `exec14d_adaptive_tN_*.png`, referenced only
  by nine frames appended after the original 119-frame deck.

The fixed display and adaptive analysis are intentionally separate products.
The former reproduces the requested presentation; the latter remains the
data-reach study and does not overwrite the main figures.

## Fixed t20 fit rule

`FIXED_DISPLAY_BIN_WIDTH_NS = 0.025` is the documented display hook. Each
fixed-threshold panel always plots the reached-event histogram and excluded
fraction. A Gaussian core is drawn only when the two-pass 25 ps fit has finite
parameters, finite covariance, `sigma_fit > 25 ps`, and
`sigma_fit_error <= sigma_fit`. Otherwise the panel retains mean/RMS and, for
far-End `t20`, displays a note pointing to the adaptive Backup.

All 30 non-central far-End fixed `t20` cores converged under this rule. Thus no
main panel requires the non-fitable note in the current EJ-230 dataset. This is
not a fabricated recovery: the edge panels retain their low reach and explicit
excluded fractions.

## Key-position far-End evidence

Source: `results_ej230_analysis/csv/exec14e_fixed_tN_summary.csv`.

| x [mm] | far side | R(20) | excluded | mean Npe | fitted sigma(t20) [ps] |
|---:|:---|---:|---:|---:|---:|
| -690 | End-R | 10.85% | 89.15% | 13.765 | 1350.6 +/- 132.2 |
| -650 | End-R | 16.35% | 83.65% | 15.006 | 1409.8 +/- 121.6 |
| -400 | End-R | 88.35% | 11.65% | 27.998 | 794.3 +/- 27.8 |
| +400 | End-L | 88.50% | 11.50% | 27.773 | 821.6 +/- 30.3 |
| +650 | End-L | 16.90% | 83.10% | 15.161 | 1363.4 +/- 118.9 |
| +690 | End-L | 11.90% | 88.10% | 13.850 | 1429.2 +/- 146.2 |

At `x=0`, neither End is designated far; both explicit End-L and End-R panels
remain present.

## Parity rule

Frames 1--119 remain the parity-controlled deck, with the existing documented
frame-48 title exception. The three EXEC_14E adaptive frames are allowed only
after frame 119. The parity and asset-parity checkers compare the primary
119-frame prefix and reject any non-adaptive appended title. The adaptive table,
synthesis, and each key-position figure occupy separate frames to avoid clipping
and keep the warning annotations legible.

## Final validation

- Primary parity: `119` frames, `118` exact titles plus the documented frame-48
  exception.
- Appended adaptive Backup: `9` frames; final PDF: `128` pages.
- Asset audit: `164/164` references resolved.
- Numeric-input audit: `128` frames, `53` imported/generated table fragments,
  and zero errors.
- OCR audit: `142` unique PNG files, zero forbidden inherited labels.
- Strict `latexmk -halt-on-error`: exit `0`; no overfull boxes, missing files,
  draft images, undefined controls, duplicate commands, or undefined references.
- PDF preflight: `112` pages with raster images, zero literal `figs/`, zero
  filename text, and zero suspicious empty image boxes.
- Visual render at 130 dpi: frames 7/11/15/19/23/27/31, frame 33, the adaptive
  table/synthesis, and all seven full-width adaptive position frames inspected.
  All main frames show six panels and fixed 25 ps bins; Backup warnings are
  legible and do not overlap the plotted distributions.

Final PDF SHA256:
`129567f8a0a753147d6c04f3b9d5d0ff922d4ff06d72b3c12add41aeb080d9cb`.
