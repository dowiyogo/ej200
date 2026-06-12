# EXEC_11 pair (28,29) timing and 1-D position reconstruction

## Provenance and status

- Branch: `exp/pair-scan-2026-06-11`
- Input data: `/home/reriosto/SHiP/t0minidaq/pairscan_2026-06-11`
- Data-producing commit: `f431c01`
- Analysis commit embedded in the final figures: `a103ed6`
- Material/readout: the recorded `pairscan.log` states `OPSC-101`, yield 10400 ph/MeV,
  rise 0.7 ns, decay 1.8 ns, attenuation 160 cm, emission peak 408.8 nm, and `Top`
  readout. Repository documentation maps OPSC-101 to EJ-204.
- Pair identity: `DetectorConstruction::TopSiPMCenterX` and global/local ID mapping put
  IDs 28 and 29 on Top at -452 and -432 mm, respectively.
- ROOT phase 0b: **DEFERRED — ROOT runtime unavailable on MSI**. Static inspection finds
  a plausible crash site: `analysis/analyze_pair_scan.C:281-283` fits with option `N`
  (do not store the TF1), while `analysis/analyze_pair_scan.C:357-360` dereferences
  `g_dt_cal.GetFunction("f_lin_dt")` after the CSV is written. The macro was not changed
  because the crash could not be reproduced without ROOT.

## Input gate

- 41 stable ROOT files, integer positions -462.0 through -422.0 mm, no missing or
  duplicate positions, and no files changed during a 31-second size/mtime check.
- All files opened with uproot. Entries range from `14396486` to
  `14768630`; every file contains exactly 3000 unique event IDs,
  IDs 0..2999, no entirely missing events, and zero non-finite `time_ns` values.
- Every position has 3000 events passing 4 PE in both channels, so efficiency is 1.0.
- `SiPMSD.cc:40-113` records a row only when the Geant4 optical boundary invokes the
  sensitive detector after the selected surface `EFFICIENCY`; PDE is not applied again.
- `count_definition = "number of sipm_hits rows per event and SiPM, using the existing EXEC_08 PE-equivalent convention"`. Axes use SiPM hit count / PE-equivalent
  count where physical photoelectron acceptance cannot be independently re-established.

## QA v1 to v2

The v2 delta-t histogram uses 20 ps bins over [-3000,3000] ps. Its Gaussian is seeded
at the maximum bin, starts from the RMS after clipping the outer 1%, and is fitted twice
inside +/-2.5 sigma. The broad-fit guard prevents convergence to the
known narrow subpeak. All 41 fits terminate successfully, but all have chi2/ndf > 5:
the distributions are visibly non-Gaussian mixtures, and that limitation is retained.

The three collapsed v1 fits become:

| x (mm) | sigma v1 (ps) | sigma v2 (ps) |
|---:|---:|---:|
| -444.0 | 24.428 | 110.569 |
| -442.0 | 40.800 | 109.983 |
| -439.0 | 22.958 | 111.419 |

The all-position temporal slope changes from
`6.5099 +/- 0.0330 ps/mm`
to `6.8384 +/- 0.0363 ps/mm`.

## Automatic references

- `POS_REF_1`: `-441.0 mm`, minimum abs(mu_dt_ps), tie -> minimum abs(x_true_mm).
- `POS_REF_2`: `-443.0 mm`, maximum mean_npe_A, tie -> minimum abs(x_true_mm).

Both selected files contain 3000 valid pair-threshold events. Individual-hit time
figures use every photon-hit row; event-time figures use the fourth-hit estimators.
NPE-shape plots use a **Moyal approximation**, not an exact Landau density.

Event-level timing/count correlations:

- `POS_REF_1`: Pearson(delta_t,R) = `-0.1084 +/- 0.0189`, Spearman = `-0.1171`.
- `POS_REF_2`: Pearson(delta_t,R) = `-0.1387 +/- 0.0186`, Spearman = `-0.1231`.

## Reconstruction

The leave-two-out temporal calibration over the other 39 positions is
`mu_delta_t = (6.83857 +/- 0.03627) x
+ (3023.230 +/- 16.037) ps`, with
chi2/ndf `4.704`.

The linear ratio calibration has chi2/ndf `3969.6` and LOOCV
RMSE `0.1017`. A cubic comparison improves training chi2 and
LOOCV RMSE to `0.0748`, but gives worse event-level resolution
and has ambiguous physical roots for about 3.3% of events. Multiple cubic roots are
resolved by choosing the physical root nearest the predefined linear estimate.

- `POS_REF_1` / `temporal`: mean `-441.258 mm`, bias `-0.258 mm`, sigma `12.494 +/- 0.161 mm`.
- `POS_REF_1` / `ratio_linear`: mean `-442.698 mm`, bias `-1.698 mm`, sigma `5.558 +/- 0.072 mm`.
- `POS_REF_1` / `BLUE`: mean `-442.502 mm`, bias `-1.502 mm`, sigma `5.265 +/- 0.068 mm`.
- `POS_REF_2` / `temporal`: mean `-442.855 mm`, bias `+0.145 mm`, sigma `11.864 +/- 0.153 mm`.
- `POS_REF_2` / `ratio_linear`: mean `-442.658 mm`, bias `+0.342 mm`, sigma `5.410 +/- 0.070 mm`.
- `POS_REF_2` / `BLUE`: mean `-442.685 mm`, bias `+0.315 mm`, sigma `5.152 +/- 0.067 mm`.

BLUE diagnostics:

- `POS_REF_1`: corr(x_t,x_R)=`0.1084`, weights=(`0.1359`, `0.8641`), condition number=`5.143`.
- `POS_REF_2`: corr(x_t,x_R)=`0.1387`, weights=(`0.1338`, `0.8662`), condition number=`4.952`.

The temporal event standard deviations are smaller than `sigma_fit/abs(m_t)` by more
than 10 sigma at both references. This guardrail failure is not tuned away: the broad
single-Gaussian v2 fit describes long tails, whereas event RMS used by reconstruction
is narrower. Both quantities and their disagreement remain in `reconstruction_summary.csv`.

## Interpretation

Lv et al. use 64 timing measurements and all 2016 SiPM pairs to intersect geometric
circles and locate an empty region. EXEC_11 uses one adjacent pair and empirical 1-D
calibrations of delta-t and count ratio. It is an inspired one-pair/one-dimensional
adaptation, not a literal reproduction of their circle algorithm. It must not be
compared directly with their 2-3 mm geometric result or approximately 1.5 mm CNN result.

Every position resolution here is a **simulation prediction — Top readout — intrinsic**.
There is no Top-readout test-beam counterpart, SPTR, electronics jitter, or walk
correction. The strong ratio non-linearity and non-Gaussian delta-t tails are open
limitations, not resolved detector performance.

## Reproduction

```bash
OUT=results/exec11_20260612_182454
pytest -q tests/test_exec11_pair_analysis.py
python3 analysis/exec11_pair_analysis.py derive --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py qa --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py detail --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py reconstruct --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py report --output-dir "$OUT"
```

Figures obtain the analysis SHA at runtime through `git rev-parse --short HEAD`.
They were regenerated after the analysis/reconstruction commit and before the final
report-only commit, so the embedded SHA identifies the exact analysis implementation.
