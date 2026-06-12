# EXEC_11 photon-arrival audit

## Provenance

- Per-photon time variable: `time_ns` from `sipm_hits` TTree (`analysis/exec07_photon_budget.py:342`).
- SUM4 clusters: IDs 0..3/4..7/8..11/12..15 (`analysis/congruent_sum4_timing.C:218-221`).
- Nearest Top channel: `TOP_POSITIONS_MM` (`common.py:23`); tie-break (`exec07_photon_budget.py:359-360`).
- Blocking validation: `exec07_photon_budget.py:89-124` (2000 IDs, gun_x_mm, IDs 0..85).

- Mode: `arrival_cumulative`
- Band: `rms`
- End grouping: `sum4_max`
- Top window: 0--30 ns, dt=0.250 ns
- End window: 0--50 ns, dt=0.250 ns
- LE threshold: 4 PE

## Anti-artifact check 1: monotone non-decreasing cumulative curves

Cumulative N(t) per event is computed as `np.cumsum` of a 2D histogram with non-negative integer entries. By the non-negativity of histogram bins, `np.cumsum` is provably non-decreasing for every event. The mean over events of non-decreasing functions is non-decreasing. No numerical exception is possible.

## Anti-artifact check 2: End windows capture the late tail

Mean N_pe at t_max=50 ns across End groups and positions: 258.7 PE.
The 50 ns End window was chosen to exceed the >30 ns late-tail cutoff established in EXEC_09 (`exec09_timing_mechanism.py`, `audit/exec09_timing_mechanism.md`). N_at_tmax is the cumulative N_pe at the window edge; values well above the FPT level confirm the window captures late-arriving photons.

## Anti-artifact check 3: Top-nearest selection vs. common.py

Top-nearest channel IDs agree with `per_position_exec07.csv` for all 31 analyzed positions. Selection is consistent with the existing analysis.

## Metrics summary (key positions)

| x_mm | group | channel | FPT [ns] | LE fire [ns] | N_at_tmax [PE] |
|------|-------|---------|----------|--------------|----------------|
| -690 | top_nearest | 16 | 0.342 | 0.500 | 318.7 |
| -690 | end_left | 4,5,6,7 | 0.257 | 0.500 | 1554.8 |
| -690 | end_right | 8,9,10,11 | 5.509 | 6.250 | 13.8 |
| -650 | top_nearest | 18 | 0.343 | 0.500 | 318.0 |
| -650 | end_left | 0,1,2,3 | 0.427 | 0.500 | 1071.7 |
| -650 | end_right | 8,9,10,11 | 5.348 | 6.250 | 15.0 |
| -400 | top_nearest | 31 | 0.346 | 0.500 | 367.2 |
| -400 | end_left | 0,1,2,3 | 1.408 | 1.750 | 301.8 |
| -400 | end_right | 8,9,10,11 | 4.372 | 5.000 | 28.0 |
| +0 | top_nearest | 50 | 0.351 | 0.500 | 366.0 |
| +0 | end_left | 0,1,2,3 | 2.875 | 3.250 | 75.9 |
| +0 | end_right | 8,9,10,11 | 2.885 | 3.250 | 76.1 |
| +400 | top_nearest | 70 | 0.346 | 0.500 | 366.3 |
| +400 | end_left | 0,1,2,3 | 4.375 | 5.000 | 27.8 |
| +400 | end_right | 8,9,10,11 | 1.407 | 1.750 | 300.0 |
| +650 | top_nearest | 83 | 0.343 | 0.500 | 319.7 |
| +650 | end_left | 0,1,2,3 | 5.338 | 6.250 | 15.2 |
| +650 | end_right | 8,9,10,11 | 0.427 | 0.500 | 1076.8 |
| +690 | top_nearest | 85 | 0.342 | 0.500 | 320.6 |
| +690 | end_left | 0,1,2,3 | 5.499 | 6.250 | 13.8 |
| +690 | end_right | 12,13,14,15 | 0.257 | 0.500 | 1564.5 |

## Scope notes

- No sigma_group metric recomputed or modified.
- Top panels are simulation-only; no test-beam counterpart.
- End panels are test-beam comparable at intrinsic (no SPTR / no jitter) precision.
- EXEC_09/10 PDFs are untouched historical checkpoints.
- runs/ directory is byte-identical to the checkpoint tag.
