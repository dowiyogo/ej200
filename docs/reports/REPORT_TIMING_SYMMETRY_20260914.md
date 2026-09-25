# Definitive timing-symmetry diagnosis

**Date:** 2026-09-14  
**Dataset:** current corrected optical transport, 21 EndTop cells, 10,000 events per cell  
**Gate 1 verdict:** `ASYMMETRY_SOURCE = FIT_DEFINITION`  
**Fit labels:** `GAUSSIAN_FIT_UNSTABLE`; `FIT_INITIALIZATION_STABLE = YES`  
**New simulation:** no

## Executive verdict

The visible left-right irregularity in the published Gaussian-core timing curve is a fit-definition artifact. It is not supported by the centered event distributions, their global RMS, or their robust central 68% widths. The Gaussian width changes by more than its quoted ROOT covariance error when the fit window, histogram binning, bin origin, or fit objective is changed. The strongest dependency is the robust-reference fit window: the median spread across the 21 cells is 3.40 ps and the maximum is 8.67 ps, while typical Gaussian statistical errors are 0.8--0.9 ps.

The primary timing-resolution metric should therefore be

\[
\sigma_{68}(T_0)=\frac{q_{84}(T_0)-q_{16}(T_0)}{2}.
\]

It is directly defined on events, robust against tails, independent of histogram binning and fit initialization, and has a bootstrap uncertainty. The Gaussian width remains useful as a core-shape diagnostic, but its ROOT fit error is not the full estimator uncertainty.

The current 10,000-event cells close the question at Gate 1. A high-statistics Geant4 campaign would reduce statistical errors but would not resolve the already measured 3--9 ps analysis-policy dependence, so none was launched.

## Inputs and observable reconstruction

The analysis reads `/home/rrios/ej200/presentations/v9/sources/timing_events.root`, which is a ROOT extraction of the 21 validated production ROOT files under `/home/rrios/exec42_20260913/grid/cells/`. The source contains 210,000 events. Every cell has `N=10000`, geometry `EndTop`, 16 END and 70 TOP sensors, `SPTR=0`, no electronics, four Geant4 workers, `eventModulo=1`, seed pair `26092601,8349041`, simulation commit `b35ee84acadef12c93506e0720580c91f901fcbf`, and the corrected transport.

For every event the reconstruction uses the same `time_ns` field and sorting/validity logic on both sides:

- LEFT END: global IDs 0--7;
- RIGHT END: global IDs 8--15;
- `tL` and `tR`: first sorted detected PE on each side;
- `T0 = (tL+tR)/2` and `DeltaT = tR-tL`;
- a valid event requires both END lists to contain the requested order statistic.

The extraction check found zero left/right count mismatches across all 210,000 events. No position-dependent cut or side-dependent operation enters `T0`.

## Distribution evidence

All 21 complete `T0` histograms, baseline Gaussian objects, fit covariance matrices, and fit statuses are stored in `timing_fit_stability.root`. The display policy is common to every panel: 4 ps bins, median-relative range of +/-0.32 ns, unit-area normalization, and a Gaussian fit over `median +/- 2 sigma68`. Every binned baseline fit returned status 0 and covariance quality 3.

Median- and mean-centered mirrored pairs were compared with identical 4 ps bins and `[-0.32,+0.32] ns` range. Across all 18 comparisons, ROOT Kolmogorov probabilities span 0.284--0.979 and Anderson-Darling probabilities span 0.183--0.940. The lowest result is the median-centered EJ-200 `|x|=200 mm` Anderson-Darling probability, 0.183; it does not reject a common centered shape. The overlays and density differences show no repeatable side-dependent tail or core deformation.

These probabilities are supporting diagnostics rather than the decision rule. The decision also uses the estimator comparisons and method variations below.

## Fit sensitivity

All tests use ROOT. The binned fits use `Minuit2/Migrad`, tolerance 0.01. The unbinned cross-check uses `RooGaussian` and RooFit. RooFit completed with fit status 0 and covariance quality 3; ROOT reported that optional `libRooFitMore`/GSL integrators were unavailable, which did not prevent these analytic Gaussian fits.

| Variation over 21 cells | Median effect on sigma_G | Largest effect |
|---|---:|---:|
| 2, 4, 8 ps bins and 0/half-bin origins, full envelope | 1.068 ps | 1.943 ps (EJ-230, -650 mm) |
| half-bin origin, 2 ps bins | 0.147 ps | 0.746 ps |
| half-bin origin, 4 ps bins | 0.226 ps | 0.646 ps |
| half-bin origin, 8 ps bins | 0.370 ps | 1.058 ps |
| fit windows +/-1.5, +/-2.0, +/-2.5 sigma68 | 3.403 ps | 8.670 ps (EJ-230, -650 mm) |
| initialization A/B/C, fixed histogram and range | 0.000 ps | 0.000 ps |

The three initializations were: A = median/sigma68, B = histogram mean/RMS, and C = histogram mode/sigma68. Their maximum sigma difference was below `1e-6 ps`, so the result is deterministic and initialization-stable. A random fit seed is neither used nor relevant.

At the endpoint, the window dependence is explicit:

| Cell | +/-1.5 sigma68 | +/-2.0 sigma68 | +/-2.5 sigma68 |
|---|---:|---:|---:|
| EJ-200, -650 mm | 84.21 ps | 79.71 ps | 75.71 ps |
| EJ-200, +650 mm | 84.50 ps | 82.71 ps | 78.00 ps |
| EJ-230, -650 mm | 86.67 ps | 81.23 ps | 78.00 ps |
| EJ-230, +650 mm | 86.72 ps | 81.79 ps | 78.91 ps |

Changing the objective also moves the result. Relative to the binned chi-square fit, the binned-likelihood sigma is higher by a median 0.723 ps and as much as 2.550 ps. The unbinned RooFit sigma differs by a median +0.661 ps and ranges from -0.093 to +2.563 ps. These method shifts are comparable to, or larger than, the nominal ROOT fit errors.

Using a common physical fit window for the `+x` and `-x` distributions does not turn the Gaussian width into a stable estimator: EJ-200 at `|x|=200 mm` remains at `z=3.18`, and EJ-230 at `|x|=500 mm` reaches `z=2.68`. This is consistent with small sample fluctuations being amplified differently by a non-Gaussian core fit; the raw centered shapes and non-parametric widths remain compatible.

## Gaussian, RMS and robust widths

RMS and `sigma68` uncertainties use 500 deterministic event-bootstrap replicas. The fixed seed series begins at `26091443`. Gaussian errors are taken from the ROOT covariance and describe only the chosen model and fit policy.

The recommended robust central values at `x=0` are:

| Material | sigma68(T0) | scan range |
|---|---:|---:|
| EJ-200 | 77.07 +/- 0.72 ps | 77.07--82.49 ps |
| EJ-204 | 72.15 +/- 0.65 ps | 72.15--85.93 ps |
| EJ-230 | 66.95 +/- 0.70 ps | 66.95--86.29 ps |

### Mirrored differences for the first-PE estimator

All differences are `width(+|x|)-width(-|x|)`.

| Material | |x| [mm] | Delta sigma_G [ps] (z) | Delta RMS [ps] (z) | Delta sigma68 [ps] (z) |
|---|---:|---:|---:|---:|
| EJ-200 | 200 | +4.25 +/- 1.27 (+3.34) | +1.51 +/- 0.82 (+1.85) | +1.77 +/- 1.04 (+1.70) |
| EJ-200 | 500 | -0.16 +/- 1.30 (-0.12) | +0.48 +/- 0.79 (+0.61) | +0.82 +/- 1.02 (+0.80) |
| EJ-200 | 650 | +2.99 +/- 1.27 (+2.35) | +0.61 +/- 0.77 (+0.78) | +0.33 +/- 1.10 (+0.30) |
| EJ-204 | 200 | -0.44 +/- 1.17 (-0.38) | -0.44 +/- 0.74 (-0.59) | +0.29 +/- 0.96 (+0.30) |
| EJ-204 | 500 | +0.91 +/- 1.27 (+0.72) | -0.09 +/- 0.83 (-0.10) | +1.27 +/- 1.15 (+1.11) |
| EJ-204 | 650 | -2.30 +/- 1.27 (-1.82) | +0.25 +/- 0.86 (+0.29) | -0.44 +/- 1.12 (-0.39) |
| EJ-230 | 200 | +0.28 +/- 1.13 (+0.24) | +0.89 +/- 0.70 (+1.28) | +0.90 +/- 0.93 (+0.96) |
| EJ-230 | 500 | +2.30 +/- 1.24 (+1.86) | +1.29 +/- 0.82 (+1.58) | +2.16 +/- 1.07 (+2.01) |
| EJ-230 | 650 | +0.56 +/- 1.23 (+0.46) | +0.45 +/- 0.94 (+0.48) | +1.31 +/- 1.30 (+1.01) |

One robust comparison is marginally above 2 sigma: EJ-230 `sigma68` at 500 mm has `z=2.01`. It is not reproduced by RMS at the same point, by either robust estimator at the other positions, or by a consistent sign/material pattern. Together with compatible shape tests, this is an isolated finite-sample fluctuation, not evidence for a physical asymmetry.

## Mean of the first five photoelectrons

The arithmetic mean of the first five detected PEs on each END was evaluated at all seven EJ-230 positions using the same events and one fixed Gaussian policy. It remains smoother and more symmetric than first PE: its mirrored Gaussian differences have `|z|<0.64`, RMS differences `|z|<1.46`, and `sigma68` differences `|z|<0.65`.

At EJ-230, `x=0`, the fixed-policy Gaussian changes from `68.23 +/- 0.83 ps` for first PE to `58.68 +/- 0.66 ps` for mean-first-five, an improvement of 9.55 ps (14.0%). The earlier adaptive-fit values, approximately `67.1 -> 57.6 ps`, therefore describe a real improvement of the estimator, although those exact Gaussian numbers inherit the unstable fit convention. Under the recommended robust definition the corresponding values are `66.95 +/- 0.70 ps -> 59.92 +/- 0.62 ps`, an improvement of 7.03 ps (10.5%). Mean-first-five remains a valid candidate, with about 60 ps robust intrinsic optical timing at the center.

## Gate decision and presentation correction

The Gate 1 conditions for Case A are met:

1. median- and mean-centered mirrored histograms have compatible shapes;
2. RMS and `sigma68` show no coherent physical left-right asymmetry;
3. `sigma_G` changes appreciably under reasonable binning, bin-origin, fit-window, and objective changes;
4. initialization changes converge to the same numerical minimum, excluding random-seed or starting-point explanations;
5. the observable implementation applies exactly the same construction to both ENDs.

The corrected presentation is `presentations/v9p1/`. Its main timing curve and performance table use `sigma68`; Gaussian widths are labeled as core diagnostics; backup slides contain all 21 distributions, mirrored overlays, sensitivity results, width-method comparisons, and the mean-first-five scan. The original `presentations/v9/` remains unchanged.

## Reproduction commands

No command below launches Geant4.

```bash
cd /home/rrios/ej200/analysis/timing_symmetry_20260914
root -l -b -q 'macros/timing_fit_stability.C+'
root -l -b -q 'macros/mirror_diagnostics.C+'
root -l -b -q 'macros/summarize_diagnostics.C+'
for m in fit_grid_EJ200 fit_grid_EJ204 fit_grid_EJ230 mirror_EJ200 mirror_EJ230 width_methods_EJ200 width_methods_EJ204 width_methods_EJ230 fit_stability_EJ200 mean5_EJ230 sigma68_vs_x; do
  root -l -b -q "macros/${m}.C"
done
```

## Artifacts

- Diagnostic ROOT/CSV: `/home/rrios/ej200/analysis/timing_symmetry_20260914/sources/timing_fit_stability.{root,csv}`
- Mirrored shapes/tests: `/home/rrios/ej200/analysis/timing_symmetry_20260914/sources/mirror_diagnostics.{root,csv}`
- Mirrored width table: `/home/rrios/ej200/analysis/timing_symmetry_20260914/sources/symmetry_metrics.{root,csv}`
- Widths and bootstrap errors: `/home/rrios/ej200/analysis/timing_symmetry_20260914/sources/timing_widths.csv`
- ROOT figures and sidecars: `/home/rrios/ej200/analysis/timing_symmetry_20260914/figures/`

Each requested figure has its ROOT macro, ROOT object file, PDF, and JSON metadata sidecar. The original production datasets were read only and retained.
