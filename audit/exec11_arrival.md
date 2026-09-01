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

Mean N_pe at t_max=50 ns across End groups and positions: 298.4 PE.
The 50 ns End window was chosen to exceed the >30 ns late-tail cutoff established in EXEC_09 (`exec09_timing_mechanism.py`, `audit/exec09_timing_mechanism.md`). N_at_tmax is the cumulative N_pe at the window edge; values well above the FPT level confirm the window captures late-arriving photons.

## Anti-artifact check 3: Top-nearest selection vs. common.py

Top-nearest channel IDs agree with `per_position_exec07.csv` for all 31 analyzed positions. Selection is consistent with the existing analysis.

## Metrics summary (key positions)

| x_mm | group | channel | FPT [ns] | LE fire [ns] | N_at_tmax [PE] |
|------|-------|---------|----------|--------------|----------------|
| -690 | top_nearest | 16 | 0.343 | 0.500 | 350.9 |
| -690 | end_left | 4,5,6,7 | 0.257 | 0.500 | 1703.1 |
| -690 | end_right | 8,9,10,11 | 5.509 | 6.250 | 20.8 |
| -650 | top_nearest | 18 | 0.343 | 0.500 | 351.0 |
| -650 | end_left | 0,1,2,3 | 0.428 | 0.500 | 1194.5 |
| -650 | end_right | 12,13,14,15 | 5.344 | 6.000 | 22.2 |
| -400 | top_nearest | 31 | 0.346 | 0.500 | 411.3 |
| -400 | end_left | 0,1,2,3 | 1.419 | 1.750 | 358.9 |
| -400 | end_right | 8,9,10,11 | 4.383 | 5.000 | 39.2 |
| +0 | top_nearest | 51 | 0.351 | 0.500 | 407.3 |
| +0 | end_left | 0,1,2,3 | 2.901 | 3.250 | 97.8 |
| +0 | end_right | 8,9,10,11 | 2.900 | 3.250 | 98.2 |
| +400 | top_nearest | 70 | 0.346 | 0.500 | 411.4 |
| +400 | end_left | 0,1,2,3 | 4.380 | 5.000 | 39.4 |
| +400 | end_right | 8,9,10,11 | 1.419 | 1.750 | 359.7 |
| +650 | top_nearest | 83 | 0.343 | 0.500 | 350.9 |
| +650 | end_left | 0,1,2,3 | 5.340 | 6.000 | 22.3 |
| +650 | end_right | 8,9,10,11 | 0.428 | 0.500 | 1189.7 |
| +690 | top_nearest | 85 | 0.342 | 0.500 | 353.0 |
| +690 | end_left | 0,1,2,3 | 5.482 | 6.250 | 20.8 |
| +690 | end_right | 12,13,14,15 | 0.257 | 0.500 | 1715.2 |

## Scope notes

- No sigma_group metric recomputed or modified.
- Top panels are simulation-only; no test-beam counterpart.
- End panels are test-beam comparable at intrinsic (no SPTR / no jitter) precision.
- EXEC_09/10 PDFs are untouched historical checkpoints.
- runs/ directory is byte-identical to the checkpoint tag.
