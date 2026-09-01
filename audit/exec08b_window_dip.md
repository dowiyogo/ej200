# EXEC_08b directed window-dip and timing-gate audit

Date: 2026-06-10

## Decision

- The Top-channel maximum mismatch is a **physical, local and symmetric
  track-window collection pattern**, not an analyzer or channel-map defect.
- It affects the strict `maximum == geometrically nearest` gate around offsets
  of 2 mm, but does not invalidate the EXEC_07 scan, N_pe values or times.
- The revised timing gate **fails in all six comparisons**. EndTop has a
  narrower same-end `DeltaT_AB` distribution than End-only, contrary to the
  proposed `1/sqrt(N_pe)`-only gate.
- Therefore EXEC_08 full regeneration and Beamer production were **not
  resumed**.

## Analyzer consistency gate

The existing x=-650 scan file and all short runs were processed by the same
per-event analyzer, `analysis/exec07/exec08b_window_dip.py`. N_pe errors are
`RMS/sqrt(N_events)`.

For the existing x=-650 file:

| Quantity | Value |
|---|---:|
| ID 18, geometrically nearest | 351.045 +/- 2.014 N_pe/event |
| ID 19, next window on +x side | 385.647 +/- 2.137 N_pe/event |
| ID19 - ID18 | 34.602 +/- 2.936 N_pe/event |
| Significance | 11.78 sigma |
| ID19 / ID18 | 1.09857 |

This exactly reproduces the original EXEC_08 cross-check. The approximately
10% effect is not an analyzer artifact.

The `global_id` to hit-position audit also remains valid: detected hits for
ID 18 have mean x=-651.99 mm in the x=-650 scan, consistent with the nominal
ID 18 center at -652 mm.

## Short-run results

All runs use EndTop, SSLG4 OPSC-101, AFBR-S4N66P024M, jitter 0, 16 workers and
2,000 events.

| Dataset | Track offset from ID18 center | ID17 N_pe | ID18 N_pe | ID19 N_pe | Relevant contrast |
|---|---:|---:|---:|---:|---|
| Existing x=-650 | +2 mm | 315.767 +/- 1.816 | 351.045 +/- 2.014 | 385.647 +/- 2.137 | ID19-ID18 = 34.602 +/- 2.936, 11.78 sigma |
| Run A x=-652 | 0 mm | 352.575 +/- 1.895 | 349.396 +/- 1.884 | 352.251 +/- 1.868 | ID19-ID18 = 2.855 +/- 2.653, 1.08 sigma |
| Run B x=-642 | midpoint ID18/19 | 296.368 +/- 1.706 | 413.197 +/- 2.361 | 412.866 +/- 2.380 | ID19-ID18 = -0.331 +/- 3.352, -0.10 sigma |
| Run C1 x=-648 | +4 mm | 297.511 +/- 1.760 | 377.797 +/- 2.217 | 398.878 +/- 2.304 | ID19-ID18 = 21.081 +/- 3.197, 6.59 sigma |
| Run C2 x=-654 | -2 mm exact mirror | 384.000 +/- 2.147 | 349.619 +/- 2.027 | 313.299 +/- 1.835 | ID17-ID18 = 34.381 +/- 2.953, 11.64 sigma |

The requested C1 coordinate, x=-648 mm, is +4 mm from ID18 at -652 mm; it is
not the exact mirror of x=-650. Run C2 at x=-654 mm was added to perform the
actual -2 mm mirror test.

The +2 mm and -2 mm contrasts agree within 0.05 sigma. The residual deficit
when centered is not significant: ID19-ID18 is 1.08 sigma and ID17-ID18 is
1.19 sigma. At the midpoint, IDs 18 and 19 agree within 0.10 sigma.

## Impact-position and mechanism assessment

The ntuple contains detected-hit `x_mm`, `z_mm`, so ID18 impact maps were
generated. It does not contain photon direction or per-photon optical-boundary
status.

| Dataset | mean local x [mm] | RMS local x [mm] | fraction abs(local x)>2.5 mm |
|---|---:|---:|---:|
| Existing x=-650 | +0.009 | 1.735 | 16.87% |
| Run A x=-652 | -0.003 | 1.733 | 16.79% |
| Run B x=-642 | +0.172 | 1.768 | 18.74% |
| Run C1 x=-648 | +0.129 | 1.750 | 17.85% |
| Run C2 x=-654 | -0.011 | 1.733 | 16.70% |

Mechanisms ruled out by code/geometry inspection:

- Lateral escape through coupling is not represented by this model. The SiPM
  surface is `dielectric_metal`, with PDE as `EFFICIENCY` and zero
  `REFLECTIVITY` (`src/Materials.cc:222-248`); a photon at that interface is
  detected or absorbed there.
- The muon does not traverse coupling/SiPM. It travels at y=0 along -Z, while
  the Top SiPM is placed near y=+30 mm. Track length in coupling/SiPM is 0 mm.
- Channel-map corruption is excluded by physical hit positions.

The remaining interpretation is a local geometrical/angular collection
pattern caused by track-to-window alignment. Its exact optical-boundary
mechanism cannot be separated with the current ntuple because photon
directions and per-channel absorbed-boundary counts are absent.

Figures:

- `analysis/exec07/figs/exec08b_window_dip_profiles.png`
- `analysis/exec07/figs/exec08b_id18_impact_maps.png`

Heavy short-run ROOT files remain local in:
`/home/reriosto/SHiP/t0minidaq/sslg4/exec08b_window_dip`.

## Time-walk premise correction

No time-walk correction exists to port:

- `analysis/congruent_sum4_timing.C:3-5` explicitly says no time-walk
  correction is applied.
- `analysis/congruent_sum4_timing.C:305-315` writes
  `HOOK_WALK=not applied in stage 1`.
- `analysis/tb_mirror_sigma_vs_x.C:3-5` says the leading edge includes
  inherent walk but no correction.
- `analysis/tb_mirror_physical_fits.py:209-214` records walk correction as
  `no`.

`HOOK_WALK` remains a future task. No correction was invented or applied in
EXEC_08b.

## Revised timing gate

The same current leading-edge implementation was used for both campaigns:
same-end SUM4 A/B, normalized 0.5/5 ns double-exponential pulse, absolute
4-PE threshold, no walk correction. The reported value is
`sigma(DeltaT_AB)/sqrt(2)`, using the distribution RMS in both datasets.

| x [mm] | side | N_pe EndTop | N_pe End-only | sigma EndTop [ps] | sigma End-only [ps] | observed ratio | predicted sqrt(Npure/Ntop) | difference |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| -400 | L | 691.286 +/- 3.941 | 601.461 +/- 1.579 | 49.965 +/- 0.790 | 66.135 +/- 0.468 | 0.756 +/- 0.013 | 0.933 +/- 0.003 | -13.22 sigma |
| -400 | R | 74.538 +/- 0.469 | 86.931 +/- 0.245 | 219.938 +/- 3.478 | 248.436 +/- 1.757 | 0.885 +/- 0.015 | 1.080 +/- 0.004 | -12.33 sigma |
| 0 | L | 194.433 +/- 1.061 | 209.048 +/- 0.535 | 107.455 +/- 1.699 | 131.397 +/- 0.929 | 0.818 +/- 0.014 | 1.037 +/- 0.003 | -15.10 sigma |
| 0 | R | 195.167 +/- 1.088 | 209.134 +/- 0.536 | 107.202 +/- 1.695 | 129.495 +/- 0.916 | 0.828 +/- 0.014 | 1.035 +/- 0.003 | -14.11 sigma |
| +400 | L | 74.984 +/- 0.481 | 86.458 +/- 0.235 | 223.776 +/- 3.539 | 247.517 +/- 1.750 | 0.904 +/- 0.016 | 1.074 +/- 0.004 | -10.54 sigma |
| +400 | R | 692.455 +/- 3.998 | 598.820 +/- 1.500 | 51.472 +/- 0.814 | 65.936 +/- 0.466 | 0.781 +/- 0.014 | 0.930 +/- 0.003 | -10.79 sigma |

The required condition `sigma_EndTop >= sigma_End-only` fails in all six
comparisons. A pure `1/sqrt(N_pe)` model is insufficient because the 70 Top
windows alter the photon-arrival-time distribution as well as total light.
Preferential extraction of late photons is a plausible interpretation, but it
requires a dedicated arrival-time comparison to establish.

The 88 ps hardware value is context only: it is a different detector without
70 draining windows and includes SPTR/electronics, so it is not directly
comparable.

Machine-readable tables:

- `analysis/exec07/exec08b_window_dip_profiles.csv`
- `analysis/exec07/exec08b_timing_gate.csv`
- `analysis/exec07/exec08b_timing_gate_raw.csv`

## Gate outcome

Part 1 passes with verdict **physical local collection pattern; scan valid**.
Part 2 fails the requested timing gate. Per EXEC_08b Part 3, full EXEC_08
regeneration, final CSV/report and Beamer decks are stopped pending a revised
timing acceptance criterion.
