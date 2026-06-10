# EXEC_09 timing mechanism audit

Date: 2026-06-10  
Branch: `feat/endtop-sslg4`

## Verdict

**CONFIRMED: the 70 Top windows preferentially drain late optical paths.**

At the three common beam positions, the normalized End-photon arrival-time
distribution has a shorter tail in EndTop than in End-only. The result uses
all End photons and is reproduced separately in each of the four End SUM4
clusters. Errors are delete-one-block jackknife errors over event IDs, which
preserve correlations among photons from the same event.

| x beam [mm] | delta t50 [ns] | delta t99 [ns] | significance t99 | delta fraction t>10 ns | significance | delta fraction t>20 ns | significance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -400 | -0.4331 +/- 0.0028 | -1.1539 +/- 0.0209 | -55.2 sigma | -0.015413 +/- 0.000118 | -130.7 sigma | -0.000133 +/- 0.000013 | -10.3 sigma |
| 0 | -0.4819 +/- 0.0028 | -1.1287 +/- 0.0206 | -54.7 sigma | -0.024442 +/- 0.000194 | -126.1 sigma | -0.000241 +/- 0.000017 | -14.4 sigma |
| +400 | -0.4344 +/- 0.0022 | -1.1509 +/- 0.0182 | -63.1 sigma | -0.015304 +/- 0.000163 | -94.0 sigma | -0.000123 +/- 0.000009 | -13.4 sigma |

The early distribution remains similar in shape, although its median advances
by 0.43--0.48 ns. The late-tail suppression is much larger and unambiguous:
all 12 individual SUM4-cluster comparisons have lower EndTop t99, with
significances from 13.4 to 56.9 sigma.

The previous gate assumption, "less End light cannot yield better intrinsic
timing", is therefore invalid for this geometry. `sigma_EndTop <
sigma_End-only` is an expected physical result because the Top windows alter
the arrival-time distribution, not merely its total photon count. The
photo-statistical `1/sqrt(Npe)` scaling alone does not describe this change.

## Anti-artifact trace

The comparison passed these checks before the mechanism was accepted:

| Quantity | Common implementation |
|---|---|
| End groups and estimator | `analysis/exec07_photon_budget.py:392-409`: same A/B SUM4 groups and `sigma(DeltaT_AB)/sqrt(2)`; End-only mirror uses the same implementation through `analysis/tb_mirror_sigma_vs_x.C:3-10`. |
| LE threshold | `analysis/congruent_sum4_timing.C:48-49`: absolute 4 PE-equivalent threshold. |
| Single-PE pulse | `analysis/congruent_sum4_timing.C:45-47,155-167`: normalized double exponential, 0.5/5 ns. |
| LE algorithm | `analysis/congruent_sum4_timing.C:170-213`; exact Python port in `analysis/exec07_photon_budget.py:146-185`. |
| Time-walk correction | None in either campaign: `analysis/congruent_sum4_timing.C:305-315`, `analysis/tb_mirror_sigma_vs_x.C:3-5`, and `analysis/tb_mirror_physical_fits.py:209-214`. |
| Hit-time convention | `src/SiPMSD.cc:75-80`: Geant4 global time plus configured jitter. EXEC_07 scan fixes jitter to zero at `scripts/run_exec07_scan.sh:57-67`. |
| Primary t=0 | Same particle gun definition, position and direction at `src/PrimaryGeneratorAction.cc:18-23`; `gun_x_mm` only records the vertex x coordinate at `src/EventAction.cc:15-20`. |
| ROOT schema | Both compared files expose identical `event_id`, `global_id`, `time_ns`, and `gun_x_mm` branches and types. |

No difference in threshold, pulse, t=0 convention, estimator, or time-walk
correction explains the timing improvement.

## Inputs and products

- EndTop: `/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000`
- End-only: `/home/reriosto/SHiP/t0minidaq/scan_end_wrapped_2026-06-09`
- Positions: `x = -400, 0, +400 mm`
- Script: `analysis/exec07/exec09_timing_mechanism.py`
- Figure: `analysis/exec07/figs/exec09_tail_comparison.png`
- Tables: `analysis/exec07/exec09_tail_metrics.csv`,
  `analysis/exec07/exec09_tail_comparison.csv`

The gate is lifted as a documented physical finding. The hardware 88 ps
reference remains context only: it is a different detector without 70
draining windows and includes electronics effects absent here.
