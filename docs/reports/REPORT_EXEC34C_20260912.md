# EXEC34C — final analysis report

| Executive status | Result |
| --- | --- |
| EXEC34C | COMPLETE |
| G-P.1 | PASS (validated pilot) |
| G-P.2 | PASS (validated pilot) |
| G-P.3 | PASS (validated pilot) |
| G-P.4 | PASS (validated pilot) |
| G-P.5 | PASS (validated pilot) |
| EJ200_xp0 | COMPLETE |
| EJ200_xp200 | COMPLETE |
| EJ200_xm200 | COMPLETE |
| EJ200_xp500 | COMPLETE |
| EJ200_xm500 | COMPLETE |
| EJ200_xp650 | COMPLETE |
| EJ200_xm650 | COMPLETE |
| EJ204_xp0 | COMPLETE |
| EJ204_xp200 | COMPLETE |
| EJ204_xm200 | COMPLETE |
| EJ204_xp500 | COMPLETE |
| EJ204_xm500 | COMPLETE |
| EJ204_xp650 | COMPLETE |
| EJ204_xm650 | COMPLETE |
| EJ230_xp0 | COMPLETE |
| EJ230_xp200 | COMPLETE |
| EJ230_xm200 | COMPLETE |
| EJ230_xp500 | COMPLETE |
| EJ230_xm500 | COMPLETE |
| EJ230_xp650 | COMPLETE |
| EJ230_xm650 | COMPLETE |

**INTRINSIC timing resolution — electronics not included.** All timing values below are in ps. TOP widths are conditional on at least N selected hits; their efficiencies always use the explicit generated-event denominator.

## Input provenance

Used campaign: `/home/rrios/exec34r_20260912/`, resolving each native ROOT exclusively from the latest `SIMULATION_COMPLETE` row per cell in [manifest.jsonl](/home/rrios/exec34r_20260912/manifest.jsonl) (SHA256 `b7fa1855722f99da9082494957f0bdd13f3242bf3bfc0a86bdd91ac9b16e91eb`). There are 21 selected simulation records. No tree scan was used to select attempts.

Excluded campaign: `/home/rrios/exec34b_20260911/`. The task identifies its 21 timeout failures as incomplete simulations. It was not read, edited or deleted in EXEC34C; no failed ROOT was used. The valid manifest carries an inherited provenance reference to that old campaign; it was not followed.

Every included run records **4 workers, eventModulo=1, N_generated=10000, seeds 26092601 8349041, EndTop with N_TOP=70**, with EJ-200=OPSC-100, EJ-204=OPSC-101 and EJ-230=OPSC-106. Values came from the valid manifest and matching simulation metadata, not macro inference. The pilot used 24 workers. Recorded EXEC33 evidence gives exact zero N_pe/end differences at 1/4/12/24 workers; its scope is preserved in every analysis metadata file. No new simulation was run here.

Unique PDE: `/home/rrios/ej200_exec33_20260911/data/sipm/AFBR-S4N66P024M_pde.txt`, SHA256 `6360bd80eabf1ade77ea0dd80f8e56ed99b6ea8dcbc2ce94bf5269108f2b28e7`, identical to the checkpoint and all included simulations. Every ROOT was read across all branches, SHA256 checked and reconciled with master run summary photon totals and hit-bearing event counts. The hit-only ROOT has no generated-event ledger: the generated denominator is the master `Events run=10000`, independently recorded in metadata; unique hit IDs are never used as this denominator. Input evidence: [verified_inputs.json](/home/rrios/exec34c_20260912/verified_inputs.json).

## Checkpoint and known non-blocking differences

The scoped gate passes: analyze.py, top_split.py, end_bridge.h, gate.py and all checkpoint-listed imported primitives have exactly the checkpoint hashes. G-P.1–G-P.5 are PASS. Canonical parity is even TRAIN, odd EVAL. The following known orchestration-only changes were explicitly allowed by the revised task:

| File | Expected checkpoint SHA256 | Actual SHA256 |
| --- | --- | --- |
| grid.py | `0f070550191c2a1c46b4415e113b4d9bfe55f1beead47537c3aaf8a707055be9` | `fab55aca362ce281a98fbbdca3d3c184786904a76101669dfd9f0de77c9b6169` |
| test_grid.py | `e4f32f280750a70eec30ce868cd7d523b28252608248f4fa03ef9e18228c93e3` | `c7610697a3be5a930410a380b3fd23ab237d59e96a8ed25218b4ab3da078d039` |

These differences come from `498d14f`. Current HEAD `d99f994807a882b12ff63e649888dbbb31f1256c` also contains launcher commit `d99f994`, ahead of checkpoint HEAD `b406143da55595e0e5f2ab81ed3d5b47c9e2cf58`, as expected. No computation script was edited. The initial checkpoint-stop report is preserved at [previous_checkpoint_stop_report.md](/home/rrios/exec34c_20260912/previous_checkpoint_stop_report.md); the revised task explicitly permits continuation.

## Pilot result

| Quantity | Validated result [ps] | Efficiency / qualification |
| --- | --- | --- |
| TOP(c), N=1 | 3.090539 ± 0.057955 | 1.000000; chi2/ndf=13.513377 |
| END FitCore ΔT_LR | 119.904019 ± 1.452208 | 1.000000 |
| END FitCore ΔT_LR/√2 | 84.784945 ± 1.026866 | Derived assumption below |
| END FitPeakSeeded ΔT_LR | 118.205413 ± 0.836748 | 1.000000 |
| END FitPeakSeeded ΔT_LR/√2 | 83.583849 ± 0.591670 | Derived assumption below |

The pilot chooses N=1 with frozen channels [50,51,52,49]. Its poor TOP Gaussian fit quality is retained as a limitation, not used to tune the grid. This is the requested checkpoint summary; no historical timing benchmark is used.

## TOP timing results

(a) is the entire N=1..20 curve in TRAIN and EVAL, tabulated per cell in Appendix A, including errors, efficiencies and fit quality. (b) chooses its argmin on the same full sample: **biased by selection**. (c), primary, learns ranking, four channel IDs, winning N and every data-derived fit seed/window/axis in even TRAIN; odd EVAL applies that frozen construction without ranking or index optimization. Reverse parity results from the unchanged validated pipeline are retained only as diagnostics and never substituted for canonical results. No calibration, walk or jitter is introduced.

| Cell | Frozen channels | N(c) | σ(c) EVAL ± error | Eff(c) | fit status(c) | chi2/ndf(c) | N(b) | σ(b) ± error, biased | Eff(b) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 | [51, 50, 52, 49] | 5 | 3.138761 ± 0.068712 | 1.000000 | 0 | 39.086864 | 1 | 2.927797 ± 0.036704 | 1.000000 |
| EJ200_xp200 | [60, 61, 59, 62] | 1 | 3.542872 ± 0.060825 | 1.000000 | 0 | 18.372800 | 1 | 3.479958 ± 0.038609 | 1.000000 |
| EJ200_xm200 | [41, 40, 42, 39] | 1 | 3.515678 ± 0.063384 | 1.000000 | 0 | 18.360375 | 1 | 3.518026 ± 0.043627 | 1.000000 |
| EJ200_xp500 | [75, 76, 74, 77] | 1 | 3.416735 ± 0.060500 | 1.000000 | 0 | 45.858850 | 1 | 3.433641 ± 0.043491 | 1.000000 |
| EJ200_xm500 | [26, 25, 27, 24] | 1 | 3.448661 ± 0.062089 | 1.000000 | 0 | 15.526107 | 1 | 3.404005 ± 0.043461 | 1.000000 |
| EJ200_xp650 | [83, 82, 84, 81] | 3 | 93.819230 ± 3.404098 | 1.000000 | 4 | 189.996242 | 3 | 1.893196 ± 0.050668 | 1.000000 |
| EJ200_xm650 | [18, 19, 17, 20] | 1 | 2.105133 ± 0.046962 | 1.000000 | 0 | 87.503746 | 1 | 1.920410 ± 0.040916 | 1.000000 |
| EJ204_xp0 | [50, 51, 52, 49] | 1 | 3.090539 ± 0.057955 | 1.000000 | 0 | 13.513377 | 1 | 3.019001 ± 0.038284 | 1.000000 |
| EJ204_xp200 | [60, 61, 59, 62] | 1 | 3.566539 ± 0.061620 | 1.000000 | 0 | 22.546296 | 1 | 3.513703 ± 0.043665 | 1.000000 |
| EJ204_xm200 | [41, 40, 42, 39] | 1 | 3.584238 ± 0.068572 | 1.000000 | 0 | 12.632358 | 1 | 3.584832 ± 0.041647 | 1.000000 |
| EJ204_xp500 | [75, 76, 74, 77] | 1 | 3.457896 ± 0.070299 | 1.000000 | 0 | 42.405924 | 1 | 3.660736 ± 0.039339 | 1.000000 |
| EJ204_xm500 | [26, 25, 27, 24] | 1 | 3.535166 ± 0.058888 | 1.000000 | 0 | 13.490991 | 1 | 3.513720 ± 0.044305 | 1.000000 |
| EJ204_xp650 | [83, 82, 84, 81] | 1 | 1.902613 ± 0.048790 | 1.000000 | 0 | 76.673201 | 1 | 1.897312 ± 0.032054 | 1.000000 |
| EJ204_xm650 | [18, 19, 17, 20] | 1 | 2.160126 ± 0.065924 | 1.000000 | 0 | 247.332830 | 1 | 1.973174 ± 0.039936 | 1.000000 |
| EJ230_xp0 | [50, 51, 49, 52] | 4 | 3.139170 ± 0.079393 | 1.000000 | 0 | 85.116292 | 1 | 3.118541 ± 0.035476 | 1.000000 |
| EJ230_xp200 | [60, 61, 59, 62] | 1 | 3.811590 ± 0.056282 | 1.000000 | 0 | 11.808889 | 1 | 3.736705 ± 0.042834 | 1.000000 |
| EJ230_xm200 | [41, 40, 42, 39] | 3 | 3.786918 ± 0.078087 | 1.000000 | 0 | 6.223879 | 1 | 3.661674 ± 0.043397 | 1.000000 |
| EJ230_xp500 | [75, 76, 74, 77] | 1 | 3.602182 ± 0.058664 | 1.000000 | 0 | 20.336125 | 1 | 3.532111 ± 0.047216 | 1.000000 |
| EJ230_xm500 | [26, 25, 27, 24] | 1 | 3.765049 ± 0.056397 | 1.000000 | 0 | 15.632813 | 1 | 3.666294 ± 0.044594 | 1.000000 |
| EJ230_xp650 | [83, 82, 84, 81] | 1 | 2.023848 ± 0.053453 | 1.000000 | 0 | 141.772934 | 1 | 1.921831 ± 0.037298 | 1.000000 |
| EJ230_xm650 | [18, 19, 17, 20] | 1 | 2.281987 ± 0.055858 | 1.000000 | 0 | 53.234297 | 1 | 2.202075 ± 0.040216 | 1.000000 |

Errors are max(ROOT covariance error, fixed-estimator bootstrap error), retaining both components in sidecars. They are conditional on the trained estimator and omit uncertainty from repeated training selections. Fit convergence does not certify Gaussian adequacy; per-cell quality and fit-window acceptance are retained.

## Selection-bias study

| Cell / material | x [mm] | σ(b) ± error | Eff(b) | σ(c) ± error | Eff(c) | (b) − (c) [ps] |
| --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 / EJ-200 | 0 | 2.927797 ± 0.036704 | 1.000000 | 3.138761 ± 0.068712 | 1.000000 | -0.210964 |
| EJ200_xp200 / EJ-200 | 200 | 3.479958 ± 0.038609 | 1.000000 | 3.542872 ± 0.060825 | 1.000000 | -0.062914 |
| EJ200_xm200 / EJ-200 | -200 | 3.518026 ± 0.043627 | 1.000000 | 3.515678 ± 0.063384 | 1.000000 | 0.002348 |
| EJ200_xp500 / EJ-200 | 500 | 3.433641 ± 0.043491 | 1.000000 | 3.416735 ± 0.060500 | 1.000000 | 0.016906 |
| EJ200_xm500 / EJ-200 | -500 | 3.404005 ± 0.043461 | 1.000000 | 3.448661 ± 0.062089 | 1.000000 | -0.044656 |
| EJ200_xp650 / EJ-200 | 650 | 1.893196 ± 0.050668 | 1.000000 | 93.819230 ± 3.404098 | 1.000000 | -91.926034 |
| EJ200_xm650 / EJ-200 | -650 | 1.920410 ± 0.040916 | 1.000000 | 2.105133 ± 0.046962 | 1.000000 | -0.184723 |
| EJ204_xp0 / EJ-204 | 0 | 3.019001 ± 0.038284 | 1.000000 | 3.090539 ± 0.057955 | 1.000000 | -0.071538 |
| EJ204_xp200 / EJ-204 | 200 | 3.513703 ± 0.043665 | 1.000000 | 3.566539 ± 0.061620 | 1.000000 | -0.052836 |
| EJ204_xm200 / EJ-204 | -200 | 3.584832 ± 0.041647 | 1.000000 | 3.584238 ± 0.068572 | 1.000000 | 0.000594 |
| EJ204_xp500 / EJ-204 | 500 | 3.660736 ± 0.039339 | 1.000000 | 3.457896 ± 0.070299 | 1.000000 | 0.202840 |
| EJ204_xm500 / EJ-204 | -500 | 3.513720 ± 0.044305 | 1.000000 | 3.535166 ± 0.058888 | 1.000000 | -0.021446 |
| EJ204_xp650 / EJ-204 | 650 | 1.897312 ± 0.032054 | 1.000000 | 1.902613 ± 0.048790 | 1.000000 | -0.005301 |
| EJ204_xm650 / EJ-204 | -650 | 1.973174 ± 0.039936 | 1.000000 | 2.160126 ± 0.065924 | 1.000000 | -0.186952 |
| EJ230_xp0 / EJ-230 | 0 | 3.118541 ± 0.035476 | 1.000000 | 3.139170 ± 0.079393 | 1.000000 | -0.020629 |
| EJ230_xp200 / EJ-230 | 200 | 3.736705 ± 0.042834 | 1.000000 | 3.811590 ± 0.056282 | 1.000000 | -0.074885 |
| EJ230_xm200 / EJ-230 | -200 | 3.661674 ± 0.043397 | 1.000000 | 3.786918 ± 0.078087 | 1.000000 | -0.125244 |
| EJ230_xp500 / EJ-230 | 500 | 3.532111 ± 0.047216 | 1.000000 | 3.602182 ± 0.058664 | 1.000000 | -0.070071 |
| EJ230_xm500 / EJ-230 | -500 | 3.666294 ± 0.044594 | 1.000000 | 3.765049 ± 0.056397 | 1.000000 | -0.098756 |
| EJ230_xp650 / EJ-230 | 650 | 1.921831 ± 0.037298 | 1.000000 | 2.023848 ± 0.053453 | 1.000000 | -0.102018 |
| EJ230_xm650 / EJ-230 | -650 | 2.202075 ± 0.040216 | 1.000000 | 2.281987 ± 0.055858 | 1.000000 | -0.079912 |

Median (b)−(c) = **-0.070071 ps**, range -91.926034 to 0.202840 ps; 17 negative, 4 positive, 0 exactly zero. This signed difference is a central methodological result, not an error or an extra gate. The unchanged pipeline also retains a legacy bias_diagnostic comparing this difference with a distinct internal TRAIN-minus-EVAL diagnostic; its historical warning string is not used to reject, tune or reinterpret (b)−(c) in EXEC34C.

Absolute magnitudes span 0.000594 to 91.926034 ps; the signed values and uncertainties above must be retained rather than summarized by a universal correction.

The sign is not consistent across cells. This is a finite-sample, overlapping-sample comparison, not a direct measurement of repeated-selection bias; variation in (b)−(c) alone does not establish why a winning index changed. The unchanged per-cell adjacency diagnostic below records whether the selected index is noise dominated. No sigma target or extra optimization was introduced.

| Cell | N(b) | N(c) | (b)−(c) | Adjacent TRAIN curve flat? | Boundary winner? | Nearest-neighbor Δσ / diagnostic error |
| --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 | 1 | 5 | -0.210964 | False | False | N=4: 0.018178 / 0.099664; N=6: 0.379174 / 0.099462 |
| EJ200_xp200 | 1 | 1 | -0.062914 | False | True | N=2: 0.183360 / 0.099820 |
| EJ200_xm200 | 1 | 1 | 0.002348 | False | True | N=2: 0.177460 / 0.088949 |
| EJ200_xp500 | 1 | 1 | 0.016906 | True | True | N=2: 0.045454 / 0.090987 |
| EJ200_xm500 | 1 | 1 | -0.044656 | False | True | N=2: 0.173561 / 0.096451 |
| EJ200_xp650 | 3 | 3 | -91.926034 | False | False | N=2: 1.238656 / 0.088812; N=4: 16.949624 / 0.476076 |
| EJ200_xm650 | 1 | 1 | -0.184723 | False | True | N=2: 0.948248 / 0.084768 |
| EJ204_xp0 | 1 | 1 | -0.071538 | False | True | N=2: 0.272364 / 0.093944 |
| EJ204_xp200 | 1 | 1 | -0.052836 | False | True | N=2: 0.119470 / 0.091754 |
| EJ204_xm200 | 1 | 1 | 0.000594 | False | True | N=2: 0.252801 / 0.096299 |
| EJ204_xp500 | 1 | 1 | 0.202840 | False | True | N=2: 0.162971 / 0.095898 |
| EJ204_xm500 | 1 | 1 | -0.021446 | False | True | N=2: 0.250160 / 0.091125 |
| EJ204_xp650 | 1 | 1 | -0.005301 | False | True | N=2: 1.119959 / 0.092714 |
| EJ204_xm650 | 1 | 1 | -0.186952 | False | True | N=2: 0.811827 / 0.089134 |
| EJ230_xp0 | 1 | 4 | -0.020629 | False | False | N=3: 0.497316 / 0.577889; N=5: 0.390175 / 0.098649 |
| EJ230_xp200 | 1 | 1 | -0.074885 | False | True | N=2: 0.109882 / 0.088866 |
| EJ230_xm200 | 1 | 3 | -0.125244 | False | False | N=2: 0.059510 / 0.102969; N=4: 8.367992 / 0.291087 |
| EJ230_xp500 | 1 | 1 | -0.070071 | False | True | N=2: 0.182824 / 0.088950 |
| EJ230_xm500 | 1 | 1 | -0.098756 | False | True | N=2: 0.127814 / 0.090889 |
| EJ230_xp650 | 1 | 1 | -0.102018 | False | True | N=2: 1.003717 / 0.080109 |
| EJ230_xm650 | 1 | 1 | -0.079912 | False | True | N=2: 0.969133 / 0.090512 |

1/21 completed cells are flagged flat by the existing rule (adjacent TRAIN width differences within one quadrature conditional error); for those cells the winning N is noise dominated and not interpretable as a physical optimum. This diagnostic uses correlated fits and conditional errors; it is not a new hypothesis test. Boundary winners cannot establish an optimum outside N=1..20. No cross-cell fit was made.

## END timing results

Fixed HOOK_MAP: left {0,1,2,3}/{4,5,6,7}, right {8,9,10,11}/{12,13,14,15}. HOOK_ENDRED takes the first finite SUM4-cluster leading-edge crossing at each end; both ends must be finite. Rise 0.5 ns, fall 5 ns, threshold 4 PE-equivalent summed amplitude. FitCore uses its preserved iterative ±2σ primitive; FitPeakSeeded uses its preserved peak-centered ±2 ns primitive. Only the END primitives were invoked.

**Derived normalization assumes equal and statistically independent timing contributions from the two ends; this assumption is not verified in EXEC_34.** The primary observable is σ(ΔT_LR). The systematic sign is explicitly **FitPeakSeeded − FitCore**; both fits use the same events, so no independent-error interpretation is assigned to that difference.

| Cell | Core σ(ΔT) ± err | Core σ(ΔT)/√2 ± err | Peak σ(ΔT) ± err | Peak σ(ΔT)/√2 ± err | Peak−Core ΔT | Peak−Core /√2 | Eff END | Core fit_used_end | Core χ²/ndf | Peak χ²/ndf |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 | 121.272063 ± 1.434561 | 85.752298 ± 1.014388 | 121.138066 ± 0.879270 | 85.657548 ± 0.621738 | -0.133997 | -0.094750 | 1.000000 | True | 0.739712 | 0.869741 |
| EJ200_xp200 | 125.218882 ± 1.509211 | 88.543121 ± 1.067173 | 125.238235 ± 0.909296 | 88.556805 ± 0.642970 | 0.019352 | 0.013684 | 1.000000 | True | 0.889375 | 0.788948 |
| EJ200_xm200 | 126.561557 ± 1.408155 | 89.492535 ± 0.995716 | 124.325342 ± 0.887556 | 87.911292 ± 0.627597 | -2.236215 | -1.581243 | 1.000000 | True | 1.292117 | 1.439844 |
| EJ200_xp500 | 140.519169 ± 1.700404 | 99.362057 ± 1.202368 | 142.800036 ± 1.084821 | 100.974874 ± 0.767084 | 2.280867 | 1.612817 | 1.000000 | True | 0.977813 | 0.776568 |
| EJ200_xm500 | 142.344115 ± 1.696891 | 100.652489 ± 1.199883 | 144.067587 ± 1.104488 | 101.871168 ± 0.780991 | 1.723472 | 1.218678 | 1.000000 | True | 1.318208 | 1.360935 |
| EJ200_xp650 | 161.669244 ± 1.923730 | 114.317419 ± 1.360283 | 162.913189 ± 1.250818 | 115.197021 ± 0.884462 | 1.243945 | 0.879602 | 1.000000 | True | 1.512931 | 1.609462 |
| EJ200_xm650 | 163.902582 ± 1.968590 | 115.896627 ± 1.392003 | 161.744632 ± 1.166015 | 114.370726 ± 0.824497 | -2.157950 | -1.525901 | 1.000000 | True | 0.653279 | 1.218795 |
| EJ204_xp0 | 119.904019 ± 1.452208 | 84.784945 ± 1.026866 | 118.205413 ± 0.836748 | 83.583849 ± 0.591670 | -1.698606 | -1.201096 | 1.000000 | True | 0.662694 | 0.974045 |
| EJ204_xp200 | 125.291950 ± 1.372663 | 88.594787 ± 0.970620 | 124.909898 ± 0.933639 | 88.324636 ± 0.660182 | -0.382052 | -0.270151 | 1.000000 | True | 1.609626 | 1.258404 |
| EJ204_xm200 | 122.433138 ± 1.404099 | 86.573302 ± 0.992848 | 122.855248 ± 0.914440 | 86.871779 ± 0.646607 | 0.422110 | 0.298477 | 1.000000 | True | 1.252328 | 1.371979 |
| EJ204_xp500 | 156.787109 ± 1.906783 | 110.865228 ± 1.348299 | 156.687739 ± 1.193438 | 110.794963 ± 0.843888 | -0.099369 | -0.070265 | 1.000000 | True | 1.452962 | 1.687147 |
| EJ204_xm500 | 159.622684 ± 1.905874 | 112.870283 ± 1.347657 | 157.783430 ± 1.159118 | 111.569733 ± 0.819620 | -1.839255 | -1.300549 | 1.000000 | True | 1.698193 | 1.669764 |
| EJ204_xp650 | 189.312682 ± 2.277612 | 133.864282 ± 1.610515 | 186.980980 ± 1.353396 | 132.215519 ± 0.956996 | -2.331703 | -1.648763 | 1.000000 | True | 1.881758 | 1.384130 |
| EJ204_xm650 | 186.827952 ± 2.229540 | 132.107312 ± 1.576523 | 183.592361 ± 1.358598 | 129.819404 ± 0.960674 | -3.235591 | -2.287908 | 1.000000 | True | 1.800391 | 1.839983 |
| EJ230_xp0 | 113.058696 ± 1.364083 | 79.944571 ± 0.964552 | 113.060102 ± 0.816839 | 79.945565 ± 0.577593 | 0.001406 | 0.000994 | 1.000000 | True | 1.209049 | 1.179282 |
| EJ230_xp200 | 121.107392 ± 1.472104 | 85.635858 ± 1.040935 | 121.310062 ± 0.898911 | 85.779168 ± 0.635626 | 0.202670 | 0.143309 | 1.000000 | True | 1.461041 | 1.344917 |
| EJ230_xm200 | 123.288376 ± 1.405531 | 87.178047 ± 0.993860 | 122.530982 ± 0.901024 | 86.642489 ± 0.637120 | -0.757394 | -0.535558 | 1.000000 | True | 1.993004 | 1.412885 |
| EJ230_xp500 | 169.089159 ± 1.959859 | 119.564091 ± 1.385830 | 160.924162 ± 1.171273 | 113.790566 ± 0.828215 | -8.164996 | -5.773524 | 1.000000 | True | 2.866872 | 3.593429 |
| EJ230_xm500 | 161.078739 ± 1.882996 | 113.899869 ± 1.331479 | 162.996204 ± 1.299090 | 115.255721 ± 0.918596 | 1.917465 | 1.355853 | 1.000000 | True | 2.400331 | 2.759069 |
| EJ230_xp650 | 205.371200 ± 2.410806 | 145.219368 ± 1.704697 | 197.357934 ± 1.454950 | 139.553133 ± 1.028805 | -8.013267 | -5.666235 | 1.000000 | True | 1.713474 | 2.804722 |
| EJ230_xm650 | 205.174899 ± 2.456881 | 145.080563 ± 1.737277 | 197.767003 ± 1.472365 | 139.842389 ± 1.041119 | -7.407896 | -5.238174 | 1.000000 | True | 2.806636 | 2.797157 |

Electronics are excluded because at least four live versions have unresolved FWHM/sigma ambiguity and SPTR_PROVENANCE.md does not establish sqrt(kN) propagation for order statistics. No SPTR, walk correction, ToT cut, BLUE or TOP/END combination is used.

## N_pe and acceptance

N_pe/end = (total END-left hits + total END-right hits)/(2 × 10000). SEM uses all 10000 generated event contributions, including zeros. TOP acceptance is at least N hits in the frozen group, using 5000 generated EVAL events; END acceptance requires both finite end times using 10000 generated events. Core-window exclusions are distinct from missing-hit/nonfinite-end exclusions.

| Cell | N generated | N_pe/end ± SEM | Zero-hit events | TOP N(c) | TOP accepted / discarded | TOP efficiency | TOP inside / outside fit window | END accepted / discarded | END efficiency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 | 10000 | 528.404700 ± 1.308583 | 0 | 5 | 5000 / 0 | 1.000000 | 3796 / 1204 | 10000 / 0 | 1.000000 |
| EJ200_xp200 | 10000 | 569.170250 ± 1.434063 | 0 | 1 | 5000 / 0 | 1.000000 | 4224 / 776 | 10000 / 0 | 1.000000 |
| EJ200_xm200 | 10000 | 570.073600 ± 1.432695 | 0 | 1 | 5000 / 0 | 1.000000 | 4312 / 688 | 10000 / 0 | 1.000000 |
| EJ200_xp500 | 10000 | 828.117450 ± 2.070033 | 0 | 1 | 5000 / 0 | 1.000000 | 4346 / 654 | 10000 / 0 | 1.000000 |
| EJ200_xm500 | 10000 | 828.777050 ± 2.036998 | 0 | 1 | 5000 / 0 | 1.000000 | 4316 / 684 | 10000 / 0 | 1.000000 |
| EJ200_xp650 | 10000 | 1301.270400 ± 3.255622 | 0 | 3 | 5000 / 0 | 1.000000 | 3032 / 1968 | 10000 / 0 | 1.000000 |
| EJ200_xm650 | 10000 | 1294.507200 ± 3.142038 | 0 | 1 | 5000 / 0 | 1.000000 | 3771 / 1229 | 10000 / 0 | 1.000000 |
| EJ204_xp0 | 10000 | 398.142100 ± 1.006536 | 0 | 1 | 5000 / 0 | 1.000000 | 3970 / 1030 | 10000 / 0 | 1.000000 |
| EJ204_xp200 | 10000 | 441.122150 ± 1.112038 | 0 | 1 | 5000 / 0 | 1.000000 | 4253 / 747 | 10000 / 0 | 1.000000 |
| EJ204_xm200 | 10000 | 441.285050 ± 1.084523 | 0 | 1 | 5000 / 0 | 1.000000 | 4160 / 840 | 10000 / 0 | 1.000000 |
| EJ204_xp500 | 10000 | 737.271100 ± 1.876977 | 0 | 1 | 5000 / 0 | 1.000000 | 4122 / 878 | 10000 / 0 | 1.000000 |
| EJ204_xm500 | 10000 | 736.438150 ± 1.856901 | 0 | 1 | 5000 / 0 | 1.000000 | 4242 / 758 | 10000 / 0 | 1.000000 |
| EJ204_xp650 | 10000 | 1267.800900 ± 3.105364 | 0 | 1 | 5000 / 0 | 1.000000 | 3790 / 1210 | 10000 / 0 | 1.000000 |
| EJ204_xm650 | 10000 | 1273.224300 ± 3.177343 | 0 | 1 | 5000 / 0 | 1.000000 | 3974 / 1026 | 10000 / 0 | 1.000000 |
| EJ230_xp0 | 10000 | 304.252450 ± 0.765931 | 0 | 4 | 5000 / 0 | 1.000000 | 4194 / 806 | 10000 / 0 | 1.000000 |
| EJ230_xp200 | 10000 | 342.119400 ± 0.841106 | 0 | 1 | 5000 / 0 | 1.000000 | 4262 / 738 | 10000 / 0 | 1.000000 |
| EJ230_xm200 | 10000 | 342.146100 ± 0.840498 | 0 | 3 | 5000 / 0 | 1.000000 | 3614 / 1386 | 10000 / 0 | 1.000000 |
| EJ230_xp500 | 10000 | 620.348250 ± 1.532946 | 0 | 1 | 5000 / 0 | 1.000000 | 4152 / 848 | 10000 / 0 | 1.000000 |
| EJ230_xm500 | 10000 | 618.967200 ± 1.533016 | 0 | 1 | 5000 / 0 | 1.000000 | 4208 / 792 | 10000 / 0 | 1.000000 |
| EJ230_xp650 | 10000 | 1126.013600 ± 2.770185 | 0 | 1 | 5000 / 0 | 1.000000 | 3829 / 1171 | 10000 / 0 | 1.000000 |
| EJ230_xm650 | 10000 | 1124.801650 ± 2.745651 | 0 | 1 | 5000 / 0 | 1.000000 | 3833 / 1167 | 10000 / 0 | 1.000000 |

## Correlation and paired design

All 21 cells share seeds **26092601 8349041**. This is a **paired design** intended to reduce variance in position/material comparisons. Consequently the seven σ_t(x) points within a material must **not be treated as statistically independent**. Subsequent σ_t(x) fits, material comparisons or hypothesis tests must account for their correlations; a naive independent-point chi-square can underestimate errors. Covariances were not estimated here and no such fit or comparison was performed. The amount of variance reduction is not measured here.

## Manifest

New append-only analysis manifest: [analysis_manifest.jsonl](/home/rrios/exec34c_20260912/analysis_manifest.jsonl). Each RUNNING and terminal record was flushed/fsynced individually. COMPLETE was appended only after all three sidecars and required values verified. The simulation manifest was not modified. No existing COMPLETE analysis needed recomputation.

| Cell | Status | Start UTC | End UTC | Exit / gate exit | Output | Log |
| --- | --- | --- | --- | --- | --- | --- |
| EJ200_xp0 | COMPLETE | 2026-09-12T17:04:26.715243+00:00 | 2026-09-12T17:08:44.881251+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xp0/9c196a7b4ca34120b0fb2cd984a7e6b7) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xp0/9c196a7b4ca34120b0fb2cd984a7e6b7.log) |
| EJ200_xp200 | COMPLETE | 2026-09-12T17:04:26.715200+00:00 | 2026-09-12T17:08:43.310632+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xp200/9c6cd258e7254050ae2c7aafe56ab206) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xp200/9c6cd258e7254050ae2c7aafe56ab206.log) |
| EJ200_xm200 | COMPLETE | 2026-09-12T17:04:26.715733+00:00 | 2026-09-12T17:08:46.332629+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xm200/d175d3907a45409a952adb6ff3144e65) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xm200/d175d3907a45409a952adb6ff3144e65.log) |
| EJ200_xp500 | COMPLETE | 2026-09-12T17:04:26.715656+00:00 | 2026-09-12T17:08:46.822544+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xp500/f247f0fee8594db894aa41defa46d336) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xp500/f247f0fee8594db894aa41defa46d336.log) |
| EJ200_xm500 | COMPLETE | 2026-09-12T17:08:43.311446+00:00 | 2026-09-12T17:13:00.618268+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xm500/03c12d175fb34c13a5c21b39f26bacf1) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xm500/03c12d175fb34c13a5c21b39f26bacf1.log) |
| EJ200_xp650 | COMPLETE | 2026-09-12T17:08:44.882116+00:00 | 2026-09-12T17:13:15.593140+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xp650/aa07d10a9a2241428b9650cd4fe92eb6) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xp650/aa07d10a9a2241428b9650cd4fe92eb6.log) |
| EJ200_xm650 | COMPLETE | 2026-09-12T17:08:46.333505+00:00 | 2026-09-12T17:13:14.453553+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ200_xm650/963e3d2d48ff42098d68bab3c0c6d81a) | [log](/home/rrios/exec34c_20260912/cells/EJ200_xm650/963e3d2d48ff42098d68bab3c0c6d81a.log) |
| EJ204_xp0 | COMPLETE | 2026-09-12T17:08:46.823505+00:00 | 2026-09-12T17:12:58.177558+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xp0/851a021e42bd4ae38e6557a0b1f8d5ab) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xp0/851a021e42bd4ae38e6557a0b1f8d5ab.log) |
| EJ204_xp200 | COMPLETE | 2026-09-12T17:12:58.178491+00:00 | 2026-09-12T17:17:14.854721+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xp200/9cf18e9bc25040398a748111772d3c8a) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xp200/9cf18e9bc25040398a748111772d3c8a.log) |
| EJ204_xm200 | COMPLETE | 2026-09-12T17:13:00.619230+00:00 | 2026-09-12T17:17:15.006600+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xm200/e394b5a220ba48f1903e2b05af8f40c5) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xm200/e394b5a220ba48f1903e2b05af8f40c5.log) |
| EJ204_xp500 | COMPLETE | 2026-09-12T17:13:14.454466+00:00 | 2026-09-12T17:17:33.422958+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xp500/52f633ddad3e4241a5a715023b79a6f4) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xp500/52f633ddad3e4241a5a715023b79a6f4.log) |
| EJ204_xm500 | COMPLETE | 2026-09-12T17:13:15.594326+00:00 | 2026-09-12T17:17:31.154434+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xm500/9be150164ef64b9abe46904583401643) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xm500/9be150164ef64b9abe46904583401643.log) |
| EJ204_xp650 | COMPLETE | 2026-09-12T17:17:14.856126+00:00 | 2026-09-12T17:21:38.833821+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xp650/c1315a70a7f94d6daca8642a2778642a) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xp650/c1315a70a7f94d6daca8642a2778642a.log) |
| EJ204_xm650 | COMPLETE | 2026-09-12T17:17:15.007574+00:00 | 2026-09-12T17:21:35.614275+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ204_xm650/92e3c08aa5d4441e83541b28295cf5c8) | [log](/home/rrios/exec34c_20260912/cells/EJ204_xm650/92e3c08aa5d4441e83541b28295cf5c8.log) |
| EJ230_xp0 | COMPLETE | 2026-09-12T17:17:31.155359+00:00 | 2026-09-12T17:21:12.081297+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xp0/85aadefe4adf40f6a5333539b7771fa0) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xp0/85aadefe4adf40f6a5333539b7771fa0.log) |
| EJ230_xp200 | COMPLETE | 2026-09-12T17:17:33.423927+00:00 | 2026-09-12T17:21:12.420553+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xp200/99dde029fec1435c8f21f4f4863068a6) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xp200/99dde029fec1435c8f21f4f4863068a6.log) |
| EJ230_xm200 | COMPLETE | 2026-09-12T17:21:12.082466+00:00 | 2026-09-12T17:24:48.611432+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xm200/c6a1ea9a1dcb41cfa7f551d21f24bc54) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xm200/c6a1ea9a1dcb41cfa7f551d21f24bc54.log) |
| EJ230_xp500 | COMPLETE | 2026-09-12T17:21:12.421430+00:00 | 2026-09-12T17:24:56.689039+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xp500/9c54e952170d41ff9d0f84ce5954d10d) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xp500/9c54e952170d41ff9d0f84ce5954d10d.log) |
| EJ230_xm500 | COMPLETE | 2026-09-12T17:21:35.615340+00:00 | 2026-09-12T17:25:16.773076+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xm500/b010d660e4ee42a589c63a320ae78d9b) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xm500/b010d660e4ee42a589c63a320ae78d9b.log) |
| EJ230_xp650 | COMPLETE | 2026-09-12T17:21:38.834787+00:00 | 2026-09-12T17:25:26.682103+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xp650/cdd45b8ca3774737a662ca01c25aff5a) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xp650/cdd45b8ca3774737a662ca01c25aff5a.log) |
| EJ230_xm650 | COMPLETE | 2026-09-12T17:24:48.612422+00:00 | 2026-09-12T17:28:16.839555+00:00 | 0 / 0 | [sidecars](/home/rrios/exec34c_20260912/cells/EJ230_xm650/301cdf120f5e48e08bbed88f5f56ae88) | [log](/home/rrios/exec34c_20260912/cells/EJ230_xm650/301cdf120f5e48e08bbed88f5f56ae88.log) |

The existing G-P gate is also saved per cell for audit. It adds no bias-sign or target-width requirement. COMPLETE denotes successful requested observables and artifact verification; the pilot criteria remain separately reported.

## Appendix A — procedure (a), complete canonical curves

### EJ200_xp0

Frozen TOP_SUM4 channels: [51, 50, 52, 49]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.237401 ± 0.053616 | 1.000000 | 0 | 14.504994 | 3.350821 ± 0.057116 | 1.000000 | 0 | 21.033955 |
| 2 | 3.407891 ± 0.068105 | 1.000000 | 0 | 19.714816 | 3.336702 ± 0.070547 | 1.000000 | 0 | 25.014702 |
| 3 | 3.371720 ± 0.318785 | 1.000000 | 0 | 105.272657 | 3.406693 ± 0.426144 | 1.000000 | 0 | 108.734921 |
| 4 | 3.108140 ± 0.069356 | 1.000000 | 0 | 92.201565 | 3.131045 ± 0.063225 | 1.000000 | 0 | 84.590559 |
| 5 | 3.089962 ± 0.071573 | 1.000000 | 0 | 49.899338 | 3.138761 ± 0.068712 | 1.000000 | 0 | 39.086864 |
| 6 | 3.469136 ± 0.069065 | 1.000000 | 0 | 34.171539 | 3.650551 ± 0.064797 | 1.000000 | 0 | 13.194615 |
| 7 | 3.821451 ± 0.073761 | 1.000000 | 0 | 17.907728 | 3.927395 ± 0.079730 | 1.000000 | 0 | 13.546545 |
| 8 | 4.025403 ± 0.080379 | 1.000000 | 0 | 32.027088 | 4.282463 ± 0.093551 | 1.000000 | 0 | 31.295391 |
| 9 | 4.670123 ± 0.103683 | 1.000000 | 0 | 33.888218 | 4.785716 ± 0.111218 | 1.000000 | 0 | 28.117830 |
| 10 | 4.961177 ± 0.114366 | 1.000000 | 0 | 17.266055 | 4.922325 ± 0.112489 | 1.000000 | 0 | 17.317569 |
| 11 | 5.153603 ± 0.122279 | 1.000000 | 0 | 28.767113 | 5.214623 ± 0.112511 | 1.000000 | 0 | 31.063548 |
| 12 | 19.238559 ± 0.810132 | 1.000000 | 0 | 72.257900 | 19.115674 ± 0.591841 | 1.000000 | 0 | 76.156800 |
| 13 | 19.919611 ± 2.918859 | 1.000000 | 0 | 76.664531 | 21.536895 ± 2.537404 | 1.000000 | 0 | 68.058530 |
| 14 | 19.012051 ± 0.354354 | 1.000000 | 0 | 40.047390 | 17.969063 ± 0.318713 | 1.000000 | 0 | 40.176971 |
| 15 | 18.769272 ± 0.327669 | 1.000000 | 0 | 27.295025 | 18.275991 ± 0.305742 | 1.000000 | 0 | 28.035943 |
| 16 | 18.820077 ± 0.456751 | 1.000000 | 0 | 34.233550 | 19.058333 ± 0.484266 | 1.000000 | 0 | 35.172064 |
| 17 | 18.762933 ± 0.582570 | 1.000000 | 0 | 28.739105 | 19.212913 ± 0.609463 | 1.000000 | 0 | 27.869675 |
| 18 | 21.347656 ± 0.649315 | 1.000000 | 0 | 40.916930 | 20.378765 ± 0.764606 | 1.000000 | 0 | 42.475805 |
| 19 | 32.310474 ± 1.458310 | 1.000000 | 0 | 49.569822 | 20.531343 ± 0.757752 | 1.000000 | 0 | 48.176584 |
| 20 | 46.121198 ± 1.179362 | 1.000000 | 0 | 38.814068 | 46.893759 ± 1.189465 | 1.000000 | 0 | 35.732729 |

### EJ200_xp200

Frozen TOP_SUM4 channels: [60, 61, 59, 62]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.531595 ± 0.064892 | 1.000000 | 0 | 26.475884 | 3.542872 ± 0.060825 | 1.000000 | 0 | 18.372800 |
| 2 | 3.714956 ± 0.075849 | 1.000000 | 0 | 3.038040 | 3.684268 ± 0.068002 | 1.000000 | 0 | 2.575631 |
| 3 | 3.780992 ± 0.083028 | 1.000000 | 0 | 13.641123 | 3.832683 ± 0.085819 | 1.000000 | 0 | 12.107536 |
| 4 | 13.282360 ± 0.328505 | 1.000000 | 0 | 107.778398 | 13.307554 ± 0.338080 | 1.000000 | 0 | 111.751132 |
| 5 | 2.509267 ± 0.104553 | 1.000000 | 2 | 310.428571 | 3.694351 ± 0.102484 | 1.000000 | 0 | 92.250921 |
| 6 | 3.940422 ± 0.078094 | 1.000000 | 0 | 24.011783 | 4.094138 ± 0.085517 | 1.000000 | 0 | 23.914161 |
| 7 | 4.145162 ± 0.078863 | 1.000000 | 0 | 11.903615 | 4.327648 ± 0.076271 | 1.000000 | 0 | 14.866351 |
| 8 | 4.572896 ± 0.081751 | 1.000000 | 0 | 19.563189 | 4.793441 ± 0.081201 | 1.000000 | 0 | 19.492641 |
| 9 | 4.819120 ± 0.093504 | 1.000000 | 0 | 17.128252 | 4.984279 ± 0.102717 | 1.000000 | 0 | 19.820361 |
| 10 | 5.191415 ± 0.124422 | 1.000000 | 0 | 33.138237 | 5.404697 ± 0.128081 | 1.000000 | 0 | 32.051878 |
| 11 | 5.697612 ± 0.134324 | 1.000000 | 0 | 54.102551 | 5.448277 ± 0.150010 | 1.000000 | 0 | 59.450843 |
| 12 | 29.888885 ± 0.689778 | 1.000000 | 0 | 75.020967 | 30.249432 ± 0.876539 | 1.000000 | 0 | 81.655834 |
| 13 | 29.760160 ± 0.701384 | 1.000000 | 0 | 59.983240 | 29.160732 ± 0.662036 | 1.000000 | 0 | 57.467183 |
| 14 | 31.725303 ± 1.387434 | 1.000000 | 0 | 49.128608 | 31.522175 ± 1.231340 | 1.000000 | 0 | 42.619941 |
| 15 | 23.428204 ± 0.950752 | 1.000000 | 0 | 40.809760 | 24.258845 ± 0.917250 | 1.000000 | 0 | 38.125110 |
| 16 | 21.901008 ± 0.769489 | 1.000000 | 0 | 36.572089 | 22.667856 ± 0.680571 | 1.000000 | 0 | 33.239954 |
| 17 | 23.398213 ± 0.472332 | 1.000000 | 0 | 19.604125 | 23.150265 ± 0.483342 | 1.000000 | 0 | 21.093199 |
| 18 | 25.256039 ± 0.627769 | 1.000000 | 0 | 23.950536 | 24.667965 ± 0.644847 | 1.000000 | 0 | 24.440833 |
| 19 | 29.075339 ± 1.115733 | 1.000000 | 0 | 27.868342 | 25.539796 ± 0.982134 | 1.000000 | 0 | 30.741007 |
| 20 | 42.238800 ± 2.866615 | 1.000000 | 0 | 24.715713 | 41.219319 ± 2.870640 | 1.000000 | 0 | 26.575843 |

### EJ200_xm200

Frozen TOP_SUM4 channels: [41, 40, 42, 39]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.488803 ± 0.060598 | 1.000000 | 0 | 22.304728 | 3.515678 ± 0.063384 | 1.000000 | 0 | 18.360375 |
| 2 | 3.666263 ± 0.065114 | 1.000000 | 0 | 3.043682 | 3.622136 ± 0.067801 | 1.000000 | 0 | 6.207129 |
| 3 | 3.888298 ± 0.079470 | 1.000000 | 0 | 8.715660 | 3.756635 ± 0.082376 | 1.000000 | 0 | 8.765047 |
| 4 | 13.440007 ± 0.372710 | 1.000000 | 0 | 98.725354 | 14.660399 ± 0.483727 | 1.000000 | 0 | 100.809601 |
| 5 | 3.932374 ± 0.106929 | 1.000000 | 0 | 90.086381 | 3.695431 ± 0.107045 | 1.000000 | 0 | 94.102961 |
| 6 | 3.892069 ± 0.095665 | 1.000000 | 0 | 42.925976 | 3.943153 ± 0.093505 | 1.000000 | 0 | 36.568888 |
| 7 | 4.364757 ± 0.086847 | 1.000000 | 0 | 23.411271 | 4.328049 ± 0.085958 | 1.000000 | 0 | 26.042611 |
| 8 | 4.531145 ± 0.094948 | 1.000000 | 0 | 16.593659 | 4.469627 ± 0.086689 | 1.000000 | 0 | 14.320073 |
| 9 | 5.029069 ± 0.105223 | 1.000000 | 0 | 31.830031 | 5.002380 ± 0.092377 | 1.000000 | 0 | 19.775425 |
| 10 | 5.146522 ± 0.116666 | 1.000000 | 0 | 29.377941 | 5.141568 ± 0.108714 | 1.000000 | 0 | 25.901942 |
| 11 | 5.618139 ± 0.152397 | 1.000000 | 0 | 51.095760 | 5.605271 ± 0.135441 | 1.000000 | 0 | 47.491991 |
| 12 | 30.784424 ± 0.752009 | 1.000000 | 0 | 69.100549 | 28.611602 ± 0.660061 | 1.000000 | 0 | 78.259387 |
| 13 | 30.019842 ± 0.795453 | 1.000000 | 0 | 66.053695 | 28.007819 ± 0.616738 | 1.000000 | 0 | 61.256043 |
| 14 | 32.367680 ± 1.070396 | 1.000000 | 0 | 52.015000 | 27.666624 ± 1.041711 | 1.000000 | 0 | 52.010997 |
| 15 | 25.651655 ± 0.814090 | 1.000000 | 0 | 41.075306 | 23.249261 ± 0.879465 | 1.000000 | 0 | 45.019647 |
| 16 | 22.794589 ± 0.435293 | 1.000000 | 0 | 15.519939 | 22.071390 ± 0.406709 | 1.000000 | 0 | 17.012541 |
| 17 | 24.789729 ± 0.519941 | 1.000000 | 0 | 25.255015 | 23.273090 ± 0.428298 | 1.000000 | 0 | 22.885699 |
| 18 | 25.633839 ± 0.711542 | 1.000000 | 0 | 33.490163 | 24.434980 ± 0.624429 | 1.000000 | 0 | 30.259038 |
| 19 | 29.209054 ± 1.024047 | 1.000000 | 0 | 29.590432 | 29.985834 ± 1.002127 | 1.000000 | 0 | 31.074751 |
| 20 | 43.305380 ± 2.944128 | 1.000000 | 0 | 30.678651 | 41.216408 ± 2.499402 | 1.000000 | 0 | 28.909905 |

### EJ200_xp500

Frozen TOP_SUM4 channels: [75, 76, 74, 77]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.522174 ± 0.063102 | 1.000000 | 0 | 41.096794 | 3.416735 ± 0.060500 | 1.000000 | 0 | 45.858850 |
| 2 | 3.567628 ± 0.065549 | 1.000000 | 0 | 6.899244 | 3.632374 ± 0.067032 | 1.000000 | 0 | 2.171770 |
| 3 | 3.674577 ± 0.073639 | 1.000000 | 0 | 16.975754 | 3.666300 ± 0.070556 | 1.000000 | 0 | 14.922889 |
| 4 | 13.032878 ± 0.296547 | 1.000000 | 0 | 87.164996 | 12.890996 ± 0.305010 | 1.000000 | 0 | 99.660570 |
| 5 | 3.760196 ± 0.103263 | 1.000000 | 0 | 100.697294 | 3.649275 ± 0.103752 | 1.000000 | 0 | 102.067087 |
| 6 | 4.012975 ± 0.091164 | 1.000000 | 0 | 28.054673 | 3.965200 ± 0.087579 | 1.000000 | 0 | 23.341167 |
| 7 | 4.271280 ± 0.090679 | 1.000000 | 0 | 19.362020 | 4.308679 ± 0.076448 | 1.000000 | 0 | 9.756558 |
| 8 | 4.663733 ± 0.088338 | 1.000000 | 0 | 22.171972 | 4.616301 ± 0.081236 | 1.000000 | 0 | 21.900574 |
| 9 | 5.035938 ± 0.102299 | 1.000000 | 0 | 32.561217 | 5.039229 ± 0.100266 | 1.000000 | 0 | 24.898462 |
| 10 | 5.344209 ± 0.115332 | 1.000000 | 0 | 33.381727 | 5.386222 ± 0.112650 | 1.000000 | 0 | 31.661289 |
| 11 | 5.547127 ± 0.139546 | 1.000000 | 0 | 72.532853 | 5.740258 ± 0.132506 | 1.000000 | 0 | 72.823037 |
| 12 | 28.890344 ± 0.666395 | 1.000000 | 0 | 79.152022 | 28.334552 ± 0.602733 | 1.000000 | 0 | 78.650823 |
| 13 | 28.207678 ± 0.709453 | 1.000000 | 0 | 65.366634 | 27.960141 ± 0.631632 | 1.000000 | 0 | 63.225000 |
| 14 | 30.555071 ± 1.357506 | 1.000000 | 0 | 56.910853 | 27.164183 ± 1.328011 | 1.000000 | 0 | 52.355898 |
| 15 | 23.022743 ± 0.886254 | 1.000000 | 0 | 45.986859 | 23.328711 ± 0.871104 | 1.000000 | 0 | 44.585161 |
| 16 | 22.498732 ± 0.385809 | 1.000000 | 0 | 16.120023 | 21.600913 ± 0.405333 | 1.000000 | 0 | 17.648708 |
| 17 | 23.571260 ± 0.429318 | 1.000000 | 0 | 20.436030 | 23.088500 ± 0.478246 | 1.000000 | 0 | 24.037684 |
| 18 | 24.826236 ± 0.616703 | 1.000000 | 0 | 25.915999 | 25.916803 ± 0.749355 | 1.000000 | 0 | 26.314768 |
| 19 | 28.537277 ± 0.907107 | 1.000000 | 0 | 27.656286 | 28.902852 ± 1.089524 | 1.000000 | 0 | 29.644700 |
| 20 | 42.773459 ± 2.350574 | 1.000000 | 0 | 25.540434 | 42.474509 ± 2.924368 | 1.000000 | 0 | 31.246967 |

### EJ200_xm500

Frozen TOP_SUM4 channels: [26, 25, 27, 24]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.503794 ± 0.067193 | 1.000000 | 0 | 17.595380 | 3.448661 ± 0.062089 | 1.000000 | 0 | 15.526107 |
| 2 | 3.677355 ± 0.069195 | 1.000000 | 0 | 4.811795 | 3.613598 ± 0.064991 | 1.000000 | 0 | 2.547067 |
| 3 | 3.731742 ± 0.078391 | 1.000000 | 0 | 10.023200 | 3.618967 ± 0.072205 | 1.000000 | 0 | 9.712910 |
| 4 | 14.968112 ± 0.507817 | 1.000000 | 0 | 95.365858 | 16.285079 ± 0.640806 | 1.000000 | 0 | 99.350012 |
| 5 | 3.975326 ± 0.098186 | 1.000000 | 0 | 94.209165 | 0.426863 ± 0.117813 | 1.000000 | 2 | 311.785714 |
| 6 | 4.255991 ± 0.087748 | 1.000000 | 0 | 62.208494 | 4.413853 ± 0.096177 | 1.000000 | 0 | 60.390512 |
| 7 | 4.339485 ± 0.076785 | 1.000000 | 0 | 36.473428 | 4.375397 ± 0.084573 | 1.000000 | 0 | 36.935291 |
| 8 | 4.534970 ± 0.082852 | 1.000000 | 0 | 29.403469 | 4.575079 ± 0.090919 | 1.000000 | 0 | 36.536095 |
| 9 | 4.977349 ± 0.090464 | 1.000000 | 0 | 27.707921 | 5.009207 ± 0.108381 | 1.000000 | 0 | 36.191379 |
| 10 | 5.394556 ± 0.111267 | 1.000000 | 0 | 26.279980 | 5.361540 ± 0.122282 | 1.000000 | 0 | 29.458607 |
| 11 | 5.741310 ± 0.143597 | 1.000000 | 0 | 46.209086 | 5.580428 ± 0.146435 | 1.000000 | 0 | 50.608282 |
| 12 | 28.592148 ± 0.649730 | 1.000000 | 0 | 70.194874 | 30.182984 ± 0.703912 | 1.000000 | 0 | 78.303064 |
| 13 | 28.678496 ± 0.562823 | 1.000000 | 0 | 51.839138 | 29.717707 ± 0.578285 | 1.000000 | 0 | 56.875404 |
| 14 | 29.116697 ± 1.209567 | 1.000000 | 0 | 46.927193 | 31.170546 ± 0.952858 | 1.000000 | 0 | 48.492211 |
| 15 | 23.759481 ± 0.816902 | 1.000000 | 0 | 44.537741 | 24.035926 ± 0.740472 | 1.000000 | 0 | 44.141222 |
| 16 | 23.617931 ± 0.746623 | 1.000000 | 0 | 38.196554 | 23.893108 ± 0.608651 | 1.000000 | 0 | 34.890265 |
| 17 | 24.232227 ± 0.505533 | 1.000000 | 0 | 23.306510 | 24.040378 ± 0.457333 | 1.000000 | 0 | 25.264430 |
| 18 | 27.655233 ± 0.799919 | 1.000000 | 0 | 27.984237 | 24.930295 ± 0.646419 | 1.000000 | 0 | 32.839817 |
| 19 | 32.423423 ± 1.214455 | 1.000000 | 0 | 25.836762 | 25.828162 ± 0.908393 | 1.000000 | 0 | 32.396020 |
| 20 | 43.577644 ± 2.852535 | 1.000000 | 0 | 24.290398 | 40.896494 ± 1.870353 | 1.000000 | 0 | 31.225635 |

### EJ200_xp650

Frozen TOP_SUM4 channels: [83, 82, 84, 81]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.018774 ± 0.059550 | 1.000000 | 0 | 79.591512 | 2.080535 ± 0.053478 | 1.000000 | 0 | 59.911472 |
| 2 | 3.098580 ± 0.065916 | 1.000000 | 0 | 42.167668 | 3.105758 ± 0.059698 | 1.000000 | 0 | 42.359331 |
| 3 | 1.859924 ± 0.059520 | 1.000000 | 0 | 122.631590 | 93.819230 ± 3.404098 | 1.000000 | 4 | 189.996242 |
| 4 | 18.809548 ± 0.472340 | 1.000000 | 0 | 116.977532 | 2.857919 ± 0.213571 | 1.000000 | 0 | 87.737326 |
| 5 | 2.268119 ± 0.985056 | 1.000000 | 0 | 128.765362 | 2.268406 ± 0.552682 | 1.000000 | 0 | 135.791657 |
| 6 | 5.854539 ± 0.358127 | 1.000000 | 0 | 127.431167 | 5.707776 ± 0.345678 | 1.000000 | 0 | 129.569188 |
| 7 | 6.233486 ± 0.331112 | 1.000000 | 0 | 101.588048 | 6.036865 ± 0.319687 | 1.000000 | 0 | 104.366902 |
| 8 | 6.869459 ± 0.189660 | 1.000000 | 0 | 34.542197 | 6.754725 ± 0.202735 | 1.000000 | 0 | 38.424717 |
| 9 | 6.894237 ± 0.135816 | 1.000000 | 0 | 18.506992 | 6.921714 ± 0.137037 | 1.000000 | 0 | 15.459544 |
| 10 | 7.091618 ± 0.132650 | 1.000000 | 0 | 3.133421 | 7.041744 ± 0.126593 | 1.000000 | 0 | 3.070805 |
| 11 | 7.752756 ± 0.136843 | 1.000000 | 0 | 7.157636 | 7.606943 ± 0.132173 | 1.000000 | 0 | 6.452820 |
| 12 | 7.939100 ± 0.139892 | 1.000000 | 0 | 8.123371 | 7.873288 ± 0.134603 | 1.000000 | 0 | 6.914505 |
| 13 | 7.600465 ± 0.146547 | 1.000000 | 0 | 6.784042 | 7.774082 ± 0.138359 | 1.000000 | 0 | 4.636659 |
| 14 | 7.873679 ± 0.148881 | 1.000000 | 0 | 14.923471 | 7.982915 ± 0.148379 | 1.000000 | 0 | 11.137751 |
| 15 | 7.953930 ± 0.162512 | 1.000000 | 0 | 22.912875 | 7.966760 ± 0.171662 | 1.000000 | 0 | 22.029609 |
| 16 | 126.117465 ± 5.237439 | 1.000000 | 0 | 87.150624 | 96.096128 ± 6.362375 | 1.000000 | 0 | 82.812448 |
| 17 | 44.473171 ± 0.832160 | 1.000000 | 0 | 53.438464 | 46.024431 ± 0.882884 | 1.000000 | 0 | 55.943005 |
| 18 | 47.909226 ± 0.916107 | 1.000000 | 0 | 45.959663 | 47.726337 ± 0.880198 | 1.000000 | 0 | 46.608028 |
| 19 | 55.034707 ± 1.485486 | 1.000000 | 0 | 35.130254 | 55.141483 ± 1.355459 | 1.000000 | 0 | 34.962989 |
| 20 | 61.734253 ± 1.701114 | 1.000000 | 0 | 28.905925 | 59.895152 ± 1.577727 | 1.000000 | 0 | 30.144194 |

### EJ200_xm650

Frozen TOP_SUM4 channels: [18, 19, 17, 20]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.140435 ± 0.056000 | 1.000000 | 0 | 76.173530 | 2.105133 ± 0.046962 | 1.000000 | 0 | 87.503746 |
| 2 | 3.088683 ± 0.063637 | 1.000000 | 0 | 39.078202 | 3.169375 ± 0.065823 | 1.000000 | 0 | 34.189039 |
| 3 | 2.918187 ± 0.053011 | 1.000000 | 0 | 33.787465 | 2.917117 ± 0.058475 | 1.000000 | 0 | 30.477606 |
| 4 | 18.910043 ± 0.535041 | 1.000000 | 0 | 82.277585 | 17.448935 ± 0.956280 | 1.000000 | 0 | 84.961123 |
| 5 | 2.439957 ± 0.070377 | 1.000000 | 0 | 123.390520 | 2.365354 ± 6.666880 | 1.000000 | 0 | 131.321058 |
| 6 | 5.787829 ± 0.388573 | 1.000000 | 0 | 123.459556 | 5.443006 ± 0.382895 | 1.000000 | 0 | 121.021899 |
| 7 | 6.032925 ± 0.367029 | 1.000000 | 0 | 108.476106 | 5.724448 ± 0.455607 | 1.000000 | 0 | 108.884640 |
| 8 | 6.720795 ± 0.248950 | 1.000000 | 0 | 77.300285 | 6.901559 ± 0.198614 | 1.000000 | 0 | 71.406354 |
| 9 | 6.902956 ± 0.127037 | 1.000000 | 0 | 7.967957 | 6.743190 ± 0.134249 | 1.000000 | 0 | 13.766729 |
| 10 | 7.135220 ± 0.133142 | 1.000000 | 0 | 4.181539 | 7.180076 ± 0.137737 | 1.000000 | 0 | 3.118172 |
| 11 | 7.404269 ± 0.140641 | 1.000000 | 0 | 1.678731 | 7.412182 ± 0.142849 | 1.000000 | 0 | 5.092244 |
| 12 | 7.760213 ± 0.140450 | 1.000000 | 0 | 2.927705 | 7.648738 ± 0.135408 | 1.000000 | 0 | 5.407698 |
| 13 | 7.848607 ± 0.145459 | 1.000000 | 0 | 7.676523 | 7.731838 ± 0.141723 | 1.000000 | 0 | 8.570348 |
| 14 | 7.940796 ± 0.161823 | 1.000000 | 0 | 12.272694 | 7.798519 ± 0.156459 | 1.000000 | 0 | 12.420522 |
| 15 | 8.101995 ± 0.170606 | 1.000000 | 0 | 21.561619 | 8.016383 ± 0.171642 | 1.000000 | 0 | 23.586513 |
| 16 | 105.654886 ± 5.059224 | 1.000000 | 0 | 86.685705 | 95.633318 ± 5.447922 | 1.000000 | 0 | 84.205865 |
| 17 | 44.958186 ± 0.810533 | 1.000000 | 0 | 54.695182 | 46.417575 ± 0.862808 | 1.000000 | 0 | 58.035334 |
| 18 | 44.828229 ± 0.934688 | 1.000000 | 0 | 45.385240 | 47.841205 ± 0.887155 | 1.000000 | 0 | 45.435393 |
| 19 | 44.755479 ± 1.221431 | 1.000000 | 0 | 46.314058 | 48.120410 ± 1.280554 | 1.000000 | 0 | 50.584649 |
| 20 | 53.045342 ± 1.599393 | 1.000000 | 0 | 34.466747 | 58.425355 ± 1.692648 | 1.000000 | 0 | 32.744549 |

### EJ204_xp0

Frozen TOP_SUM4 channels: [50, 51, 52, 49]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.100759 ± 0.054975 | 1.000000 | 0 | 16.947409 | 3.090539 ± 0.057955 | 1.000000 | 0 | 13.513377 |
| 2 | 3.373123 ± 0.076179 | 1.000000 | 0 | 18.647971 | 3.336515 ± 0.063596 | 1.000000 | 0 | 12.695868 |
| 3 | 3.401361 ± 0.587677 | 1.000000 | 0 | 113.908183 | 3.528672 ± 0.621628 | 1.000000 | 0 | 112.561446 |
| 4 | 3.177971 ± 0.075047 | 1.000000 | 0 | 60.690324 | 3.229751 ± 0.066552 | 1.000000 | 0 | 61.792418 |
| 5 | 3.426568 ± 0.068739 | 1.000000 | 0 | 22.233376 | 3.496071 ± 0.059630 | 1.000000 | 0 | 24.071622 |
| 6 | 3.529280 ± 0.071204 | 1.000000 | 0 | 18.497167 | 3.553687 ± 0.073186 | 1.000000 | 0 | 26.029858 |
| 7 | 4.013759 ± 0.083348 | 1.000000 | 0 | 26.037678 | 3.926915 ± 0.083318 | 1.000000 | 0 | 22.261155 |
| 8 | 4.062530 ± 0.084597 | 1.000000 | 0 | 13.885036 | 4.094950 ± 0.091469 | 1.000000 | 0 | 13.052819 |
| 9 | 4.864017 ± 0.123334 | 1.000000 | 0 | 29.515269 | 4.486185 ± 0.108219 | 1.000000 | 0 | 28.536234 |
| 10 | 4.908286 ± 0.124864 | 1.000000 | 0 | 23.296145 | 5.150997 ± 0.114684 | 1.000000 | 0 | 18.810215 |
| 11 | 5.185754 ± 0.126473 | 1.000000 | 0 | 27.444386 | 5.295428 ± 0.116273 | 1.000000 | 0 | 26.037140 |
| 12 | 27.122457 ± 1.117279 | 1.000000 | 0 | 66.248963 | 28.564677 ± 0.940251 | 1.000000 | 0 | 72.009938 |
| 13 | 38.379946 ± 2.778949 | 1.000000 | 0 | 54.293606 | 36.806840 ± 3.598040 | 1.000000 | 0 | 55.993847 |
| 14 | 20.217620 ± 0.388557 | 1.000000 | 0 | 31.819996 | 20.291802 ± 0.396002 | 1.000000 | 0 | 34.535886 |
| 15 | 19.322884 ± 0.315249 | 1.000000 | 0 | 17.461094 | 19.374715 ± 0.306038 | 1.000000 | 0 | 19.855442 |
| 16 | 19.481220 ± 0.394029 | 1.000000 | 0 | 21.909405 | 19.391323 ± 0.384044 | 1.000000 | 0 | 22.257674 |
| 17 | 20.293241 ± 0.512528 | 1.000000 | 0 | 23.538258 | 19.472642 ± 0.524019 | 1.000000 | 0 | 23.928625 |
| 18 | 19.564335 ± 0.652893 | 1.000000 | 0 | 22.053558 | 19.913337 ± 0.631020 | 1.000000 | 0 | 20.336688 |
| 19 | 24.914861 ± 1.049673 | 1.000000 | 0 | 33.171703 | 23.740632 ± 0.890734 | 1.000000 | 0 | 29.234499 |
| 20 | 30.195282 ± 2.289293 | 1.000000 | 0 | 35.121015 | 29.093662 ± 1.801249 | 1.000000 | 0 | 29.796197 |

### EJ204_xp200

Frozen TOP_SUM4 channels: [60, 61, 59, 62]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.540698 ± 0.061501 | 1.000000 | 0 | 24.235721 | 3.566539 ± 0.061620 | 1.000000 | 0 | 22.546296 |
| 2 | 3.660169 ± 0.068091 | 1.000000 | 0 | 2.904264 | 3.605214 ± 0.066719 | 1.000000 | 0 | 2.133098 |
| 3 | 3.861713 ± 0.076927 | 1.000000 | 0 | 8.583433 | 3.663486 ± 0.070994 | 1.000000 | 0 | 11.437656 |
| 4 | 11.064914 ± 0.193974 | 1.000000 | 0 | 102.896477 | 3.791949 ± 0.186149 | 1.000000 | 0 | 93.213536 |
| 5 | 3.909865 ± 0.097900 | 1.000000 | 0 | 88.930158 | 3.927597 ± 0.098139 | 1.000000 | 0 | 92.605712 |
| 6 | 4.081339 ± 0.092483 | 1.000000 | 0 | 33.204443 | 4.213191 ± 0.092597 | 1.000000 | 0 | 31.578165 |
| 7 | 4.270951 ± 0.078146 | 1.000000 | 0 | 15.725004 | 4.248369 ± 0.088187 | 1.000000 | 0 | 15.403446 |
| 8 | 4.528094 ± 0.085667 | 1.000000 | 0 | 17.523025 | 4.406559 ± 0.090796 | 1.000000 | 0 | 18.730550 |
| 9 | 5.206871 ± 0.097296 | 1.000000 | 0 | 18.864829 | 5.017690 ± 0.106004 | 1.000000 | 0 | 26.530348 |
| 10 | 5.255902 ± 0.111503 | 1.000000 | 0 | 27.167876 | 5.337696 ± 0.121466 | 1.000000 | 0 | 24.516215 |
| 11 | 5.445547 ± 0.148740 | 1.000000 | 0 | 37.929771 | 5.655395 ± 0.148387 | 1.000000 | 0 | 36.499376 |
| 12 | 26.092444 ± 0.603707 | 1.000000 | 0 | 59.083732 | 27.447741 ± 0.585903 | 1.000000 | 0 | 57.106748 |
| 13 | 28.195636 ± 0.632195 | 1.000000 | 0 | 53.259207 | 27.709164 ± 0.636624 | 1.000000 | 0 | 45.464800 |
| 14 | 27.283084 ± 0.969254 | 1.000000 | 0 | 39.204698 | 28.756000 ± 1.150593 | 1.000000 | 0 | 38.901453 |
| 15 | 24.904073 ± 0.870657 | 1.000000 | 0 | 33.148600 | 25.333690 ± 0.804835 | 1.000000 | 0 | 31.872636 |
| 16 | 23.250872 ± 0.413107 | 1.000000 | 0 | 11.330279 | 23.588479 ± 0.391865 | 1.000000 | 0 | 9.589875 |
| 17 | 23.848032 ± 0.427036 | 1.000000 | 0 | 13.042720 | 23.244511 ± 0.397203 | 1.000000 | 0 | 14.190368 |
| 18 | 24.808321 ± 0.496641 | 1.000000 | 0 | 12.091731 | 24.439656 ± 0.469951 | 1.000000 | 0 | 13.370247 |
| 19 | 28.037961 ± 0.724057 | 1.000000 | 0 | 16.871114 | 26.451840 ± 0.714672 | 1.000000 | 0 | 18.612624 |
| 20 | 32.284456 ± 1.063268 | 1.000000 | 0 | 18.678465 | 30.913873 ± 0.980982 | 1.000000 | 0 | 18.058389 |

### EJ204_xm200

Frozen TOP_SUM4 channels: [41, 40, 42, 39]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.460003 ± 0.062843 | 1.000000 | 0 | 7.986314 | 3.584238 ± 0.068572 | 1.000000 | 0 | 12.632358 |
| 2 | 3.712804 ± 0.072967 | 1.000000 | 0 | 3.507506 | 3.679096 ± 0.069864 | 1.000000 | 0 | 2.759791 |
| 3 | 3.761003 ± 0.082417 | 1.000000 | 0 | 11.587469 | 3.661915 ± 0.079250 | 1.000000 | 0 | 14.317483 |
| 4 | 14.007189 ± 0.400622 | 1.000000 | 0 | 119.266912 | 14.684723 ± 0.438464 | 1.000000 | 0 | 116.942672 |
| 5 | 3.896186 ± 0.097105 | 1.000000 | 0 | 86.604210 | 4.054359 ± 0.101612 | 1.000000 | 0 | 77.232615 |
| 6 | 4.060093 ± 0.088815 | 1.000000 | 0 | 25.584206 | 3.996689 ± 0.084278 | 1.000000 | 0 | 26.938410 |
| 7 | 4.281697 ± 0.091761 | 1.000000 | 0 | 23.070911 | 4.215797 ± 0.088450 | 1.000000 | 0 | 22.469564 |
| 8 | 4.637986 ± 0.095566 | 1.000000 | 0 | 17.322515 | 4.525778 ± 0.084504 | 1.000000 | 0 | 19.580294 |
| 9 | 5.052624 ± 0.094371 | 1.000000 | 0 | 22.883793 | 4.965482 ± 0.101886 | 1.000000 | 0 | 26.108776 |
| 10 | 5.563043 ± 0.120622 | 1.000000 | 0 | 25.129653 | 5.507345 ± 0.132486 | 1.000000 | 0 | 29.607828 |
| 11 | 5.656828 ± 0.155460 | 1.000000 | 0 | 49.719487 | 5.649417 ± 0.143391 | 1.000000 | 0 | 51.003138 |
| 12 | 26.631171 ± 0.702212 | 1.000000 | 0 | 67.407494 | 27.645042 ± 0.604099 | 1.000000 | 0 | 62.206030 |
| 13 | 26.895970 ± 0.751762 | 1.000000 | 0 | 48.818579 | 29.534240 ± 0.634381 | 1.000000 | 0 | 45.537944 |
| 14 | 27.919046 ± 1.232950 | 1.000000 | 0 | 40.963893 | 29.292066 ± 1.084192 | 1.000000 | 0 | 40.354263 |
| 15 | 25.065557 ± 0.859465 | 1.000000 | 0 | 38.174634 | 23.244047 ± 0.814423 | 1.000000 | 0 | 40.195864 |
| 16 | 23.484037 ± 0.762375 | 1.000000 | 0 | 32.846413 | 23.572630 ± 0.798773 | 1.000000 | 0 | 32.029233 |
| 17 | 23.593282 ± 0.393924 | 1.000000 | 0 | 14.316383 | 23.787248 ± 0.435331 | 1.000000 | 0 | 11.924764 |
| 18 | 24.679003 ± 0.462582 | 1.000000 | 0 | 16.057338 | 25.021653 ± 0.515807 | 1.000000 | 0 | 14.638920 |
| 19 | 27.221413 ± 0.717953 | 1.000000 | 0 | 20.543918 | 26.430399 ± 0.747142 | 1.000000 | 0 | 21.748860 |
| 20 | 30.755690 ± 0.946834 | 1.000000 | 0 | 19.609330 | 31.156800 ± 1.126439 | 1.000000 | 0 | 22.827570 |

### EJ204_xp500

Frozen TOP_SUM4 channels: [75, 76, 74, 77]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.536919 ± 0.068461 | 1.000000 | 0 | 45.231770 | 3.457896 ± 0.070299 | 1.000000 | 0 | 42.405924 |
| 2 | 3.699890 ± 0.067153 | 1.000000 | 0 | 3.604401 | 3.627869 ± 0.068728 | 1.000000 | 0 | 3.315690 |
| 3 | 3.895079 ± 0.073199 | 1.000000 | 0 | 9.048393 | 3.757192 ± 0.078733 | 1.000000 | 0 | 9.994029 |
| 4 | 13.358040 ± 0.341317 | 1.000000 | 0 | 94.886477 | 14.187935 ± 0.411274 | 1.000000 | 0 | 101.715865 |
| 5 | 3.766372 ± 0.105071 | 1.000000 | 0 | 96.012056 | 6.310737 ± 0.101276 | 1.000000 | 2 | 278.312500 |
| 6 | 4.138940 ± 0.095405 | 1.000000 | 0 | 20.311097 | 4.136369 ± 0.095222 | 1.000000 | 0 | 27.242128 |
| 7 | 4.142545 ± 0.084427 | 1.000000 | 0 | 24.172233 | 4.206228 ± 0.091091 | 1.000000 | 0 | 30.401102 |
| 8 | 4.631167 ± 0.089928 | 1.000000 | 0 | 28.421474 | 4.716711 ± 0.087243 | 1.000000 | 0 | 24.878894 |
| 9 | 4.882668 ± 0.102074 | 1.000000 | 0 | 26.693219 | 4.880899 ± 0.096450 | 1.000000 | 0 | 24.585213 |
| 10 | 5.166020 ± 0.114226 | 1.000000 | 0 | 25.755842 | 5.359573 ± 0.127190 | 1.000000 | 0 | 24.415975 |
| 11 | 5.421722 ± 0.139052 | 1.000000 | 0 | 45.348656 | 5.650789 ± 0.144732 | 1.000000 | 0 | 39.901867 |
| 12 | 27.499877 ± 0.682826 | 1.000000 | 0 | 57.230029 | 28.899274 ± 0.675659 | 1.000000 | 0 | 58.499857 |
| 13 | 27.309745 ± 0.654064 | 1.000000 | 0 | 42.968131 | 27.259148 ± 0.621298 | 1.000000 | 0 | 46.450888 |
| 14 | 29.214808 ± 1.102355 | 1.000000 | 0 | 35.949915 | 28.742897 ± 1.051398 | 1.000000 | 0 | 38.034865 |
| 15 | 23.889429 ± 0.931365 | 1.000000 | 0 | 36.141278 | 24.801610 ± 0.814133 | 1.000000 | 0 | 33.537399 |
| 16 | 23.351726 ± 0.393042 | 1.000000 | 0 | 13.976405 | 24.351370 ± 0.407006 | 1.000000 | 0 | 9.494439 |
| 17 | 23.848352 ± 0.440138 | 1.000000 | 0 | 13.600309 | 24.407023 ± 0.435681 | 1.000000 | 0 | 12.850941 |
| 18 | 25.196658 ± 0.520604 | 1.000000 | 0 | 12.071452 | 25.473214 ± 0.474934 | 1.000000 | 0 | 16.592752 |
| 19 | 26.579221 ± 0.781661 | 1.000000 | 0 | 20.038850 | 27.967556 ± 0.668410 | 1.000000 | 0 | 17.820554 |
| 20 | 33.297003 ± 1.114409 | 1.000000 | 0 | 18.795548 | 32.458833 ± 1.120708 | 1.000000 | 0 | 20.226196 |

### EJ204_xm500

Frozen TOP_SUM4 channels: [26, 25, 27, 24]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.506608 ± 0.059631 | 1.000000 | 0 | 13.901456 | 3.535166 ± 0.058888 | 1.000000 | 0 | 13.490991 |
| 2 | 3.756768 ± 0.068906 | 1.000000 | 0 | 4.077538 | 3.721412 ± 0.067496 | 1.000000 | 0 | 2.100980 |
| 3 | 3.690282 ± 0.085393 | 1.000000 | 0 | 12.697348 | 3.801495 ± 0.079733 | 1.000000 | 0 | 13.921357 |
| 4 | 16.041427 ± 0.616232 | 1.000000 | 0 | 123.367483 | 14.463885 ± 0.459295 | 1.000000 | 0 | 115.590161 |
| 5 | 3.973806 ± 0.099047 | 1.000000 | 0 | 96.924035 | 3.666982 ± 0.098798 | 1.000000 | 0 | 100.495914 |
| 6 | 4.041951 ± 0.086880 | 1.000000 | 0 | 39.984724 | 4.157922 ± 0.096399 | 1.000000 | 0 | 29.505749 |
| 7 | 4.305264 ± 0.086771 | 1.000000 | 0 | 28.303442 | 4.240138 ± 0.088199 | 1.000000 | 0 | 24.108533 |
| 8 | 4.712240 ± 0.091322 | 1.000000 | 0 | 22.927232 | 4.591827 ± 0.095080 | 1.000000 | 0 | 29.645498 |
| 9 | 5.039413 ± 0.100557 | 1.000000 | 0 | 24.270106 | 5.048726 ± 0.110856 | 1.000000 | 0 | 24.662860 |
| 10 | 5.443589 ± 0.128402 | 1.000000 | 0 | 23.744814 | 4.908510 ± 0.113689 | 1.000000 | 0 | 24.655676 |
| 11 | 5.705417 ± 0.150262 | 1.000000 | 0 | 48.857292 | 5.580032 ± 0.147083 | 1.000000 | 0 | 45.144793 |
| 12 | 24.604187 ± 0.602677 | 1.000000 | 0 | 65.037229 | 26.826717 ± 0.623854 | 1.000000 | 0 | 69.201167 |
| 13 | 28.847837 ± 0.601615 | 1.000000 | 0 | 48.867198 | 27.408529 ± 0.719840 | 1.000000 | 0 | 53.128678 |
| 14 | 30.315957 ± 0.843488 | 1.000000 | 0 | 41.168788 | 27.450926 ± 1.060433 | 1.000000 | 0 | 38.626013 |
| 15 | 24.424844 ± 0.747828 | 1.000000 | 0 | 37.412722 | 24.939060 ± 0.785833 | 1.000000 | 0 | 32.973361 |
| 16 | 24.017942 ± 0.688863 | 1.000000 | 0 | 30.726125 | 22.840609 ± 0.609898 | 1.000000 | 0 | 28.992293 |
| 17 | 23.772763 ± 0.402960 | 1.000000 | 0 | 11.983820 | 23.251546 ± 0.411357 | 1.000000 | 0 | 14.071597 |
| 18 | 24.874205 ± 0.470762 | 1.000000 | 0 | 14.517626 | 23.941147 ± 0.499397 | 1.000000 | 0 | 15.921768 |
| 19 | 25.926999 ± 0.608910 | 1.000000 | 0 | 17.023902 | 27.206212 ± 0.679859 | 1.000000 | 0 | 15.821578 |
| 20 | 30.264680 ± 0.887436 | 1.000000 | 0 | 19.055188 | 31.189617 ± 0.990007 | 1.000000 | 0 | 19.428683 |

### EJ204_xp650

Frozen TOP_SUM4 channels: [83, 82, 84, 81]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.936607 ± 0.048536 | 1.000000 | 0 | 77.874086 | 1.902613 ± 0.048790 | 1.000000 | 0 | 76.673201 |
| 2 | 3.056566 ± 0.078995 | 1.000000 | 0 | 35.078681 | 3.162192 ± 0.073980 | 1.000000 | 0 | 29.859721 |
| 3 | 3.033198 ± 0.055488 | 1.000000 | 0 | 37.710068 | 3.086033 ± 0.056091 | 1.000000 | 0 | 33.615294 |
| 4 | 18.540056 ± 0.451303 | 1.000000 | 0 | 75.730413 | 18.313498 ± 0.532291 | 1.000000 | 0 | 80.207367 |
| 5 | 5.693848 ± 0.432603 | 1.000000 | 0 | 107.903769 | 5.204772 ± 0.461820 | 1.000000 | 0 | 107.183655 |
| 6 | 6.096040 ± 0.400655 | 1.000000 | 0 | 132.699561 | 5.767493 ± 0.419654 | 1.000000 | 0 | 130.801079 |
| 7 | 6.162438 ± 0.281866 | 1.000000 | 0 | 108.356973 | 6.506342 ± 0.327158 | 1.000000 | 0 | 104.320584 |
| 8 | 6.710118 ± 0.256004 | 1.000000 | 0 | 97.426018 | 7.436383 ± 0.377046 | 1.000000 | 0 | 90.391911 |
| 9 | 6.990320 ± 0.147322 | 1.000000 | 0 | 19.358360 | 7.233248 ± 0.149753 | 1.000000 | 0 | 15.901847 |
| 10 | 7.241143 ± 0.129843 | 1.000000 | 0 | 5.114850 | 7.366757 ± 0.149668 | 1.000000 | 0 | 4.350024 |
| 11 | 7.138399 ± 0.135152 | 1.000000 | 0 | 3.466645 | 7.398641 ± 0.156227 | 1.000000 | 0 | 2.111400 |
| 12 | 7.622534 ± 0.141372 | 1.000000 | 0 | 9.965033 | 7.856301 ± 0.143957 | 1.000000 | 0 | 5.819282 |
| 13 | 7.650181 ± 0.148049 | 1.000000 | 0 | 9.037266 | 7.815292 ± 0.144043 | 1.000000 | 0 | 7.353772 |
| 14 | 7.947904 ± 0.162930 | 1.000000 | 0 | 11.597898 | 7.967262 ± 0.159958 | 1.000000 | 0 | 11.206110 |
| 15 | 8.177936 ± 0.189127 | 1.000000 | 0 | 26.851066 | 8.212354 ± 0.175997 | 1.000000 | 0 | 24.687966 |
| 16 | 8.569815 ± 0.230599 | 1.000000 | 0 | 40.790872 | 8.367281 ± 0.200788 | 1.000000 | 0 | 37.736205 |
| 17 | 38.529673 ± 1.003736 | 1.000000 | 0 | 56.373913 | 41.802381 ± 0.940450 | 1.000000 | 0 | 58.894947 |
| 18 | 39.801053 ± 1.018497 | 1.000000 | 0 | 48.197908 | 44.172292 ± 0.883088 | 1.000000 | 0 | 49.276164 |
| 19 | 45.443498 ± 1.399133 | 1.000000 | 0 | 34.956175 | 48.269461 ± 1.554942 | 1.000000 | 0 | 33.704838 |
| 20 | 58.906422 ± 2.254053 | 1.000000 | 0 | 29.492211 | 58.681882 ± 2.237370 | 1.000000 | 0 | 26.455302 |

### EJ204_xm650

Frozen TOP_SUM4 channels: [18, 19, 17, 20]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.168804 ± 0.061810 | 1.000000 | 0 | 256.195780 | 2.160126 ± 0.065924 | 1.000000 | 0 | 247.332830 |
| 2 | 2.980631 ± 0.064222 | 1.000000 | 0 | 51.461031 | 3.112922 ± 0.058436 | 1.000000 | 0 | 44.086834 |
| 3 | 3.009109 ± 0.057558 | 1.000000 | 0 | 63.695549 | 2.974129 ± 0.057395 | 1.000000 | 0 | 61.027365 |
| 4 | 19.244268 ± 0.493911 | 1.000000 | 0 | 100.291596 | 19.334007 ± 0.467027 | 1.000000 | 0 | 94.088841 |
| 5 | 2.556768 ± 1.634970 | 1.000000 | 0 | 111.792998 | 138.663799 ± 1.337131 | 1.000000 | 4 | 154.756116 |
| 6 | 5.975557 ± 0.404939 | 1.000000 | 0 | 125.209985 | 5.951857 ± 0.448379 | 1.000000 | 0 | 128.095485 |
| 7 | 6.434107 ± 0.317487 | 1.000000 | 0 | 87.437247 | 6.188926 ± 0.306614 | 1.000000 | 0 | 92.078647 |
| 8 | 6.690484 ± 0.733213 | 1.000000 | 0 | 77.887806 | 6.878534 ± 1.244623 | 1.000000 | 0 | 74.553343 |
| 9 | 6.922213 ± 0.140659 | 1.000000 | 0 | 8.046324 | 6.865502 ± 0.126677 | 1.000000 | 0 | 5.718914 |
| 10 | 7.224585 ± 0.130970 | 1.000000 | 0 | 4.090781 | 7.442956 ± 0.129848 | 1.000000 | 0 | 4.113168 |
| 11 | 7.278418 ± 0.131960 | 1.000000 | 0 | 1.397083 | 7.526345 ± 0.130630 | 1.000000 | 0 | 5.317948 |
| 12 | 7.596564 ± 0.140104 | 1.000000 | 0 | 4.582672 | 7.653004 ± 0.138823 | 1.000000 | 0 | 8.045359 |
| 13 | 7.796058 ± 0.161937 | 1.000000 | 0 | 5.337715 | 7.759494 ± 0.152578 | 1.000000 | 0 | 5.961699 |
| 14 | 8.093490 ± 0.170152 | 1.000000 | 0 | 16.529140 | 8.335180 ± 0.184447 | 1.000000 | 0 | 14.896622 |
| 15 | 8.069120 ± 0.161032 | 1.000000 | 0 | 30.857792 | 8.613584 ± 0.191692 | 1.000000 | 0 | 31.655384 |
| 16 | 29.190699 ± 1.125084 | 1.000000 | 0 | 70.143565 | 9.047699 ± 0.233881 | 1.000000 | 0 | 44.992438 |
| 17 | 40.258380 ± 0.974838 | 1.000000 | 0 | 49.935599 | 41.070346 ± 0.853739 | 1.000000 | 0 | 49.416594 |
| 18 | 41.278791 ± 1.028878 | 1.000000 | 0 | 42.909570 | 40.478638 ± 0.837198 | 1.000000 | 0 | 44.231505 |
| 19 | 49.313042 ± 1.385485 | 1.000000 | 0 | 29.135843 | 44.864252 ± 1.440936 | 1.000000 | 0 | 36.175903 |
| 20 | 54.611766 ± 1.852230 | 1.000000 | 0 | 26.349551 | 54.178181 ± 1.985313 | 1.000000 | 0 | 28.884634 |

### EJ230_xp0

Frozen TOP_SUM4 channels: [50, 51, 49, 52]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.137654 ± 0.054774 | 1.000000 | 0 | 14.033328 | 3.148920 ± 0.051215 | 1.000000 | 0 | 12.919262 |
| 2 | 3.307159 ± 0.067151 | 1.000000 | 0 | 17.170097 | 3.288032 ± 0.067753 | 1.000000 | 0 | 9.995751 |
| 3 | 3.594065 ± 0.572761 | 1.000000 | 0 | 88.671853 | 3.562326 ± 0.564790 | 1.000000 | 0 | 86.337312 |
| 4 | 3.096749 ± 0.076818 | 1.000000 | 0 | 85.042490 | 3.139170 ± 0.079393 | 1.000000 | 0 | 85.116292 |
| 5 | 3.486924 ± 0.061893 | 1.000000 | 0 | 27.124239 | 3.535336 ± 0.065947 | 1.000000 | 0 | 25.880646 |
| 6 | 3.667414 ± 0.065966 | 1.000000 | 0 | 22.837165 | 3.660010 ± 0.066129 | 1.000000 | 0 | 14.276835 |
| 7 | 3.933762 ± 0.078074 | 1.000000 | 0 | 25.806209 | 3.829380 ± 0.079269 | 1.000000 | 0 | 30.892409 |
| 8 | 4.477221 ± 0.091964 | 1.000000 | 0 | 25.867832 | 4.379631 ± 0.090445 | 1.000000 | 0 | 33.668843 |
| 9 | 4.826366 ± 0.104026 | 1.000000 | 0 | 31.944474 | 4.850227 ± 0.119034 | 1.000000 | 0 | 30.027755 |
| 10 | 5.075411 ± 0.114534 | 1.000000 | 0 | 19.781310 | 4.996309 ± 0.134977 | 1.000000 | 0 | 23.447531 |
| 11 | 5.479934 ± 0.127616 | 1.000000 | 0 | 26.271314 | 5.241767 ± 0.129338 | 1.000000 | 0 | 30.287766 |
| 12 | 32.408542 ± 0.934636 | 1.000000 | 0 | 57.999182 | 31.917054 ± 0.957572 | 1.000000 | 0 | 56.762122 |
| 13 | 32.155296 ± 2.820015 | 1.000000 | 0 | 49.802714 | 31.847150 ± 2.225803 | 1.000000 | 0 | 46.803046 |
| 14 | 19.285849 ± 1.344554 | 1.000000 | 0 | 46.600664 | 21.081876 ± 1.304693 | 1.000000 | 0 | 44.826208 |
| 15 | 18.688003 ± 0.287201 | 1.000000 | 0 | 15.091171 | 19.558306 ± 0.306427 | 1.000000 | 0 | 12.705314 |
| 16 | 19.951623 ± 0.356869 | 1.000000 | 0 | 18.291609 | 19.708522 ± 0.335557 | 1.000000 | 0 | 17.126976 |
| 17 | 20.041331 ± 0.389533 | 1.000000 | 0 | 16.630998 | 19.992070 ± 0.414031 | 1.000000 | 0 | 20.957453 |
| 18 | 19.957434 ± 0.480583 | 1.000000 | 0 | 20.843662 | 20.107234 ± 0.537947 | 1.000000 | 0 | 22.304199 |
| 19 | 20.292666 ± 0.554560 | 1.000000 | 0 | 20.299743 | 20.525842 ± 0.626087 | 1.000000 | 0 | 19.284421 |
| 20 | 20.769602 ± 0.691366 | 1.000000 | 0 | 18.835699 | 21.777515 ± 0.760553 | 1.000000 | 0 | 22.041435 |

### EJ230_xp200

Frozen TOP_SUM4 channels: [60, 61, 59, 62]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.811166 ± 0.055667 | 1.000000 | 0 | 14.391984 | 3.811590 ± 0.056282 | 1.000000 | 0 | 11.808889 |
| 2 | 3.921048 ± 0.069270 | 1.000000 | 0 | 3.068848 | 3.863897 ± 0.069004 | 1.000000 | 0 | 5.161834 |
| 3 | 3.823641 ± 0.082297 | 1.000000 | 0 | 8.913751 | 3.790806 ± 0.077296 | 1.000000 | 0 | 9.404501 |
| 4 | 10.382799 ± 0.172857 | 1.000000 | 0 | 110.877570 | 10.593465 ± 0.167942 | 1.000000 | 0 | 94.381255 |
| 5 | 4.102368 ± 0.096608 | 1.000000 | 0 | 72.819475 | 4.203869 ± 0.107868 | 1.000000 | 0 | 74.950391 |
| 6 | 4.083304 ± 0.086413 | 1.000000 | 0 | 21.004476 | 4.080293 ± 0.084918 | 1.000000 | 0 | 19.157990 |
| 7 | 4.375054 ± 0.084573 | 1.000000 | 0 | 11.309807 | 4.289210 ± 0.085930 | 1.000000 | 0 | 16.462021 |
| 8 | 4.535663 ± 0.082630 | 1.000000 | 0 | 19.503749 | 4.451381 ± 0.090751 | 1.000000 | 0 | 20.733812 |
| 9 | 5.033547 ± 0.103371 | 1.000000 | 0 | 18.737325 | 4.978168 ± 0.092617 | 1.000000 | 0 | 18.679303 |
| 10 | 5.399909 ± 0.126213 | 1.000000 | 0 | 25.739339 | 5.129294 ± 0.127149 | 1.000000 | 0 | 27.305758 |
| 11 | 5.744348 ± 0.158902 | 1.000000 | 0 | 51.292866 | 5.752150 ± 0.136635 | 1.000000 | 0 | 49.531573 |
| 12 | 26.452921 ± 0.669084 | 1.000000 | 0 | 62.010575 | 26.803988 ± 0.759282 | 1.000000 | 0 | 68.812681 |
| 13 | 26.339235 ± 0.560011 | 1.000000 | 0 | 50.688128 | 25.275252 ± 0.537370 | 1.000000 | 0 | 51.656462 |
| 14 | 27.294666 ± 0.839101 | 1.000000 | 0 | 46.272342 | 26.784888 ± 0.769450 | 1.000000 | 0 | 40.101174 |
| 15 | 23.335521 ± 0.874775 | 1.000000 | 0 | 37.044961 | 25.007481 ± 0.804613 | 1.000000 | 0 | 34.870737 |
| 16 | 23.226770 ± 0.745531 | 1.000000 | 0 | 31.421802 | 25.824799 ± 0.687342 | 1.000000 | 0 | 26.771236 |
| 17 | 23.928620 ± 0.568189 | 1.000000 | 0 | 21.320074 | 23.752386 ± 0.524549 | 1.000000 | 0 | 20.848911 |
| 18 | 24.245982 ± 0.456622 | 1.000000 | 0 | 10.187029 | 23.848011 ± 0.455310 | 1.000000 | 0 | 13.155372 |
| 19 | 24.586669 ± 0.492964 | 1.000000 | 0 | 11.382359 | 24.847035 ± 0.505197 | 1.000000 | 0 | 15.610101 |
| 20 | 26.078020 ± 0.561946 | 1.000000 | 0 | 11.952037 | 25.089268 ± 0.548589 | 1.000000 | 0 | 13.361818 |

### EJ230_xm200

Frozen TOP_SUM4 channels: [41, 40, 42, 39]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.728444 ± 0.061014 | 1.000000 | 0 | 9.633481 | 3.731994 ± 0.056655 | 1.000000 | 0 | 9.120477 |
| 2 | 3.774257 ± 0.067567 | 1.000000 | 0 | 3.804944 | 3.857141 ± 0.068543 | 1.000000 | 0 | 2.500800 |
| 3 | 3.714747 ± 0.077700 | 1.000000 | 0 | 10.955159 | 3.786918 ± 0.078087 | 1.000000 | 0 | 6.223879 |
| 4 | 12.082739 ± 0.280525 | 1.000000 | 0 | 81.481173 | 11.838761 ± 0.298476 | 1.000000 | 0 | 91.017661 |
| 5 | 3.980894 ± 0.095543 | 1.000000 | 0 | 70.099619 | 3.934003 ± 0.113672 | 1.000000 | 0 | 73.339981 |
| 6 | 3.949093 ± 0.092516 | 1.000000 | 0 | 27.880403 | 4.147759 ± 0.085962 | 1.000000 | 0 | 20.734813 |
| 7 | 4.274896 ± 0.088377 | 1.000000 | 0 | 18.652742 | 4.328226 ± 0.082730 | 1.000000 | 0 | 18.750598 |
| 8 | 4.652534 ± 0.089868 | 1.000000 | 0 | 14.967567 | 4.623440 ± 0.090763 | 1.000000 | 0 | 12.665376 |
| 9 | 5.028665 ± 0.104722 | 1.000000 | 0 | 25.298731 | 4.914866 ± 0.096787 | 1.000000 | 0 | 25.190220 |
| 10 | 5.254958 ± 0.120602 | 1.000000 | 0 | 28.343851 | 5.253340 ± 0.124994 | 1.000000 | 0 | 28.538937 |
| 11 | 5.584074 ± 0.151420 | 1.000000 | 0 | 45.773277 | 5.858480 ± 0.160712 | 1.000000 | 0 | 42.366087 |
| 12 | 24.340148 ± 0.635551 | 1.000000 | 0 | 56.388132 | 26.103726 ± 0.791732 | 1.000000 | 0 | 59.201415 |
| 13 | 23.548266 ± 0.693184 | 1.000000 | 0 | 40.256394 | 24.812930 ± 0.553267 | 1.000000 | 0 | 46.306379 |
| 14 | 23.804136 ± 0.701145 | 1.000000 | 0 | 37.803569 | 24.362950 ± 0.670395 | 1.000000 | 0 | 39.862132 |
| 15 | 22.521399 ± 0.882333 | 1.000000 | 0 | 36.383652 | 25.564160 ± 1.081123 | 1.000000 | 0 | 31.910937 |
| 16 | 22.716996 ± 0.660440 | 1.000000 | 0 | 26.091558 | 24.110162 ± 0.740345 | 1.000000 | 0 | 28.620485 |
| 17 | 23.154481 ± 0.369787 | 1.000000 | 0 | 7.623956 | 22.034254 ± 0.389174 | 1.000000 | 0 | 8.460265 |
| 18 | 23.528798 ± 0.424537 | 1.000000 | 0 | 10.042654 | 23.172685 ± 0.420882 | 1.000000 | 0 | 9.952537 |
| 19 | 25.510075 ± 0.466897 | 1.000000 | 0 | 10.048239 | 24.513903 ± 0.486578 | 1.000000 | 0 | 12.039038 |
| 20 | 26.419269 ± 0.596533 | 1.000000 | 0 | 11.330860 | 25.821221 ± 0.593810 | 1.000000 | 0 | 13.618953 |

### EJ230_xp500

Frozen TOP_SUM4 channels: [75, 76, 74, 77]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.665159 ± 0.059999 | 1.000000 | 0 | 11.131890 | 3.602182 ± 0.058664 | 1.000000 | 0 | 20.336125 |
| 2 | 3.847983 ± 0.065668 | 1.000000 | 0 | 1.713290 | 3.788495 ± 0.065683 | 1.000000 | 0 | 4.725823 |
| 3 | 3.939437 ± 0.085281 | 1.000000 | 0 | 10.680096 | 3.857126 ± 0.079411 | 1.000000 | 0 | 13.065721 |
| 4 | 12.371387 ± 0.267512 | 1.000000 | 0 | 108.916706 | 12.978427 ± 0.305027 | 1.000000 | 0 | 116.856235 |
| 5 | 3.974420 ± 0.104966 | 1.000000 | 0 | 102.988838 | 4.061128 ± 1.090136 | 1.000000 | 0 | 105.333550 |
| 6 | 3.960734 ± 0.084004 | 1.000000 | 0 | 37.689170 | 3.968114 ± 0.084598 | 1.000000 | 0 | 32.691710 |
| 7 | 4.136744 ± 0.084054 | 1.000000 | 0 | 10.261949 | 4.230386 ± 0.084606 | 1.000000 | 0 | 10.965829 |
| 8 | 4.554233 ± 0.083341 | 1.000000 | 0 | 21.866042 | 4.613017 ± 0.091173 | 1.000000 | 0 | 17.635921 |
| 9 | 4.986158 ± 0.097564 | 1.000000 | 0 | 17.901837 | 4.825679 ± 0.100646 | 1.000000 | 0 | 19.058319 |
| 10 | 5.315255 ± 0.117242 | 1.000000 | 0 | 23.570651 | 5.114697 ± 0.127607 | 1.000000 | 0 | 27.355091 |
| 11 | 5.567368 ± 0.145116 | 1.000000 | 0 | 37.074101 | 5.431931 ± 0.147953 | 1.000000 | 0 | 34.405739 |
| 12 | 28.675617 ± 1.037125 | 1.000000 | 0 | 70.855601 | 28.055245 ± 1.137678 | 1.000000 | 0 | 69.156008 |
| 13 | 22.321456 ± 0.578619 | 1.000000 | 0 | 46.485082 | 23.950579 ± 0.586903 | 1.000000 | 0 | 49.636766 |
| 14 | 27.108610 ± 0.946171 | 1.000000 | 0 | 40.367340 | 26.078174 ± 0.858944 | 1.000000 | 0 | 40.381191 |
| 15 | 25.004128 ± 0.906304 | 1.000000 | 0 | 33.007618 | 25.324415 ± 0.988750 | 1.000000 | 0 | 32.908462 |
| 16 | 23.680199 ± 0.734488 | 1.000000 | 0 | 24.112398 | 24.559606 ± 0.834531 | 1.000000 | 0 | 24.530223 |
| 17 | 22.926019 ± 0.435813 | 1.000000 | 0 | 8.708970 | 23.434744 ± 0.448572 | 1.000000 | 0 | 7.724020 |
| 18 | 22.856516 ± 0.473079 | 1.000000 | 0 | 12.928610 | 22.223052 ± 0.471762 | 1.000000 | 0 | 15.014021 |
| 19 | 23.692343 ± 0.492249 | 1.000000 | 0 | 11.492781 | 24.147381 ± 0.476562 | 1.000000 | 0 | 10.093828 |
| 20 | 24.595416 ± 0.550641 | 1.000000 | 0 | 13.881665 | 25.169075 ± 0.585167 | 1.000000 | 0 | 14.467635 |

### EJ230_xm500

Frozen TOP_SUM4 channels: [26, 25, 27, 24]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3.687153 ± 0.058240 | 1.000000 | 0 | 14.585192 | 3.765049 ± 0.056397 | 1.000000 | 0 | 15.632813 |
| 2 | 3.814967 ± 0.069778 | 1.000000 | 0 | 5.865080 | 3.886682 ± 0.067937 | 1.000000 | 0 | 2.275171 |
| 3 | 3.869743 ± 0.083601 | 1.000000 | 0 | 8.457575 | 4.013107 ± 0.086724 | 1.000000 | 0 | 7.093911 |
| 4 | 12.732231 ± 0.278823 | 1.000000 | 0 | 90.206799 | 11.428820 ± 0.255377 | 1.000000 | 0 | 89.924419 |
| 5 | 3.923595 ± 0.099728 | 1.000000 | 0 | 85.419697 | 3.991629 ± 0.097997 | 1.000000 | 0 | 82.077018 |
| 6 | 4.183928 ± 0.091373 | 1.000000 | 0 | 17.932719 | 4.089805 ± 0.093366 | 1.000000 | 0 | 17.621792 |
| 7 | 4.276880 ± 0.079960 | 1.000000 | 0 | 19.570603 | 4.316580 ± 0.082874 | 1.000000 | 0 | 18.038739 |
| 8 | 4.575086 ± 0.087877 | 1.000000 | 0 | 14.613238 | 4.618291 ± 0.084219 | 1.000000 | 0 | 13.866305 |
| 9 | 4.855392 ± 0.099811 | 1.000000 | 0 | 24.030245 | 5.066429 ± 0.101864 | 1.000000 | 0 | 18.230046 |
| 10 | 5.464782 ± 0.131705 | 1.000000 | 0 | 26.007786 | 5.514390 ± 0.139517 | 1.000000 | 0 | 25.358843 |
| 11 | 5.756866 ± 0.165371 | 1.000000 | 0 | 38.250659 | 5.675947 ± 0.166945 | 1.000000 | 0 | 39.929186 |
| 12 | 25.843981 ± 0.885962 | 1.000000 | 0 | 55.490291 | 26.816717 ± 0.832108 | 1.000000 | 0 | 57.779975 |
| 13 | 24.838952 ± 0.592019 | 1.000000 | 0 | 38.842958 | 25.140974 ± 0.553831 | 1.000000 | 0 | 40.884818 |
| 14 | 25.362348 ± 0.867910 | 1.000000 | 0 | 39.728620 | 26.877638 ± 0.746717 | 1.000000 | 0 | 34.632252 |
| 15 | 24.308636 ± 0.783700 | 1.000000 | 0 | 32.803948 | 24.832210 ± 0.771934 | 1.000000 | 0 | 30.666083 |
| 16 | 23.708151 ± 0.723031 | 1.000000 | 0 | 24.120322 | 23.789129 ± 0.702977 | 1.000000 | 0 | 25.066494 |
| 17 | 22.084093 ± 0.514208 | 1.000000 | 0 | 21.351256 | 23.260749 ± 0.511642 | 1.000000 | 0 | 17.128767 |
| 18 | 23.411337 ± 0.414649 | 1.000000 | 0 | 8.683033 | 23.327098 ± 0.431214 | 1.000000 | 0 | 8.512654 |
| 19 | 25.113685 ± 0.513628 | 1.000000 | 0 | 9.855161 | 23.831556 ± 0.510378 | 1.000000 | 0 | 13.175834 |
| 20 | 26.420919 ± 0.594447 | 1.000000 | 0 | 12.764033 | 26.060270 ± 0.579879 | 1.000000 | 0 | 13.285360 |

### EJ230_xp650

Frozen TOP_SUM4 channels: [83, 82, 84, 81]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.066738 ± 0.050916 | 1.000000 | 0 | 104.015164 | 2.023848 ± 0.053453 | 1.000000 | 0 | 141.772934 |
| 2 | 3.070456 ± 0.061846 | 1.000000 | 0 | 30.712998 | 3.137308 ± 0.072919 | 1.000000 | 0 | 37.587620 |
| 3 | 3.116219 ± 0.053039 | 1.000000 | 0 | 12.285555 | 3.024217 ± 0.058217 | 1.000000 | 0 | 12.035532 |
| 4 | 17.097022 ± 0.347554 | 1.000000 | 0 | 75.716994 | 19.059421 ± 0.506743 | 1.000000 | 0 | 84.675672 |
| 5 | 2.495902 ± 0.888436 | 1.000000 | 0 | 122.580290 | 2.443378 ± 1.436874 | 1.000000 | 0 | 121.601792 |
| 6 | 6.110890 ± 0.425995 | 1.000000 | 0 | 125.299578 | 6.546569 ± 0.406765 | 1.000000 | 0 | 122.510967 |
| 7 | 6.476313 ± 0.286325 | 1.000000 | 0 | 92.968320 | 6.768531 ± 0.290290 | 1.000000 | 0 | 91.014199 |
| 8 | 7.011816 ± 1.045566 | 1.000000 | 0 | 79.986123 | 7.049059 ± 0.967725 | 1.000000 | 0 | 81.365248 |
| 9 | 7.084653 ± 0.165692 | 1.000000 | 0 | 39.535913 | 7.339507 ± 0.160403 | 1.000000 | 0 | 31.481651 |
| 10 | 7.109246 ± 0.137612 | 1.000000 | 0 | 7.454123 | 7.289049 ± 0.129611 | 1.000000 | 0 | 6.755908 |
| 11 | 7.345018 ± 0.144811 | 1.000000 | 0 | 6.098340 | 7.313326 ± 0.129793 | 1.000000 | 0 | 1.103224 |
| 12 | 7.597496 ± 0.136522 | 1.000000 | 0 | 2.535830 | 7.739522 ± 0.136611 | 1.000000 | 0 | 2.954716 |
| 13 | 7.612241 ± 0.141847 | 1.000000 | 0 | 4.394255 | 7.618744 ± 0.134370 | 1.000000 | 0 | 3.695636 |
| 14 | 8.111499 ± 0.148491 | 1.000000 | 0 | 9.602155 | 8.125235 ± 0.157426 | 1.000000 | 0 | 8.084244 |
| 15 | 8.426492 ± 0.202314 | 1.000000 | 0 | 23.026337 | 8.675544 ± 0.190467 | 1.000000 | 0 | 20.450789 |
| 16 | 8.522852 ± 0.229787 | 1.000000 | 0 | 38.394739 | 9.065843 ± 0.246197 | 1.000000 | 0 | 33.556241 |
| 17 | 41.824769 ± 1.798100 | 1.000000 | 0 | 49.409187 | 44.799568 ± 2.292316 | 1.000000 | 0 | 48.763634 |
| 18 | 38.313149 ± 0.898893 | 1.000000 | 0 | 38.855872 | 38.262400 ± 0.961380 | 1.000000 | 0 | 39.025413 |
| 19 | 42.338926 ± 1.816885 | 1.000000 | 0 | 37.857456 | 38.214056 ± 1.960782 | 1.000000 | 0 | 38.799150 |
| 20 | 50.595616 ± 2.477769 | 1.000000 | 0 | 26.128970 | 51.109355 ± 2.173066 | 1.000000 | 0 | 25.779533 |

### EJ230_xm650

Frozen TOP_SUM4 channels: [18, 19, 17, 20]. Denominator 5000 per half. Width units ps.

| N | σ TRAIN ± error | Eff TRAIN | status TRAIN | χ²/ndf TRAIN | σ EVAL ± error | Eff EVAL | status EVAL | χ²/ndf EVAL |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.290780 ± 0.059529 | 1.000000 | 0 | 52.545257 | 2.281987 ± 0.055858 | 1.000000 | 0 | 53.234297 |
| 2 | 3.259913 ± 0.068182 | 1.000000 | 0 | 34.625962 | 3.240781 ± 0.069049 | 1.000000 | 0 | 44.282392 |
| 3 | 3.098978 ± 0.058725 | 1.000000 | 0 | 48.672272 | 3.114973 ± 0.053360 | 1.000000 | 0 | 50.695617 |
| 4 | 16.716385 ± 0.334513 | 1.000000 | 0 | 99.986302 | 17.178283 ± 0.332097 | 1.000000 | 0 | 96.441291 |
| 5 | 2.626372 ± 4.241126 | 1.000000 | 0 | 124.854059 | 0.354411 ± 5.819236 | 1.000000 | 1 | 205.055470 |
| 6 | 6.557011 ± 0.411303 | 1.000000 | 0 | 127.175015 | 6.564618 ± 0.450936 | 1.000000 | 0 | 128.396263 |
| 7 | 6.255029 ± 0.312837 | 1.000000 | 0 | 91.487712 | 15.063279 ± 0.501926 | 1.000000 | 0 | 86.503196 |
| 8 | 7.189062 ± 1.110181 | 1.000000 | 0 | 78.945876 | 7.133386 ± 0.947947 | 1.000000 | 0 | 76.614506 |
| 9 | 6.899586 ± 0.126051 | 1.000000 | 0 | 24.956457 | 7.005435 ± 0.142577 | 1.000000 | 0 | 22.500887 |
| 10 | 7.093183 ± 0.135616 | 1.000000 | 0 | 5.474098 | 7.148207 ± 0.126900 | 1.000000 | 0 | 2.246248 |
| 11 | 7.362746 ± 0.143221 | 1.000000 | 0 | 2.865697 | 7.283152 ± 0.127331 | 1.000000 | 0 | 1.670654 |
| 12 | 7.520954 ± 0.144310 | 1.000000 | 0 | 4.936312 | 7.512381 ± 0.126778 | 1.000000 | 0 | 4.185925 |
| 13 | 8.022157 ± 0.138575 | 1.000000 | 0 | 7.237769 | 7.695248 ± 0.134044 | 1.000000 | 0 | 6.845668 |
| 14 | 8.037249 ± 0.162423 | 1.000000 | 0 | 9.260388 | 7.904931 ± 0.166592 | 1.000000 | 0 | 6.828934 |
| 15 | 8.165888 ± 0.184865 | 1.000000 | 0 | 22.155084 | 8.426444 ± 0.194728 | 1.000000 | 0 | 22.635805 |
| 16 | 8.446580 ± 0.229648 | 1.000000 | 0 | 26.754773 | 8.556470 ± 0.229124 | 1.000000 | 0 | 29.671000 |
| 17 | 39.959301 ± 2.110790 | 1.000000 | 0 | 49.315074 | 39.949864 ± 2.000319 | 1.000000 | 0 | 43.899979 |
| 18 | 37.019493 ± 0.737892 | 1.000000 | 0 | 34.906059 | 35.971334 ± 0.808587 | 1.000000 | 0 | 35.442392 |
| 19 | 43.246780 ± 1.953295 | 1.000000 | 0 | 29.558452 | 39.576070 ± 1.994211 | 1.000000 | 0 | 33.503567 |
| 20 | 48.014850 ± 2.234947 | 1.000000 | 0 | 23.168502 | 48.485006 ± 2.250594 | 1.000000 | 0 | 24.494955 |

## Appendix B — reproducibility

Updated UTC 2026-09-12T20:53:41.864382+00:00. Worktree `/home/rrios/ej200_exec33_20260911`, branch `diag/exec33-20260911`, HEAD `d99f994807a882b12ff63e649888dbbb31f1256c`, status `clean`; main remains `420addf0fd6029d5b2f0e235f472a8ae47f31fac`. Rollback before edits: `pre-exec34c-resume-20260912` → `d99f994807a882b12ff63e649888dbbb31f1256c`; earlier `pre-exec34c-20260912` and all previous tags retained. No computation source change or operational fix was necessary, so no new Git commit was created. Analysis-only orchestration and report scripts are outside the repositories and hashed below.

Checkpoint pipeline `cab884a9d5fb4393592e7a8ac61243ef4483bf54`, computation orchestration revision `0f3db2d09daaa7731e3b4cc54309f3fae8a50a2b`; actual analysis metadata records current pipeline HEAD and last top_split.py commit. Simulation commit `420addf0fd6029d5b2f0e235f472a8ae47f31fac`, Geant4 11.4.0. Four analysis processes were used with no timeout; simulation workers remain four as recorded, not analysis parallelism.

Commands executed (analysis-only):

```bash
/usr/bin/python /home/rrios/exec34c_20260912/execute_analysis.py verify
/usr/bin/python /home/rrios/exec34c_20260912/execute_analysis.py analyze
/usr/bin/python /home/rrios/exec34c_20260912/make_report.py
```

The `analyze` command resumes by verifying existing COMPLETE hashes/configuration and scheduling only pending or failed analysis attempts. Each attempt gets a new directory; simulation ROOTs remain at their native names, without aliases. Exact per-cell analyze.py and gate.py commands, UTC timestamps and exits are in the manifest. Internal TOP and END stage commands and times are in each metadata sidecar.

RNG: master seeds and eventModulo from the simulation manifest; bootstrap numpy.default_rng PCG64, fixed seed 20260618, 300 replicates. The launcher did not record explicit simulation RNG-state file hashes, so none are invented. Per-orientation trained model SHA256 fingerprints are present in each analysis metadata and below. No EVAL mutation of frozen state occurred.

| Cell | Native ROOT SHA256 | Frozen canonical state SHA256 | Sidecars SHA256 |
| --- | --- | --- | --- |
| EJ200_xp0 | `7cc27f115694f2acdcf4483fbb605386d4cb7a4c99a337a0abbacad624f2b3fd` | `c8727a0475cfc78ce6611cb16646164f0e0264752b7b5ec7416b4767aa7d8ecd` | analysis.root: `3a3107cd86b5c86a090e3e97cbc81fa2fec222c875b934fd3c353021c227235d`; analysis.csv: `5aa97b166e5b2e311977bd72a8632c5b74268133de9cc116f03a4e72bb8968d8`; analysis.meta.json: `1f728d4c2446616c6245b462be979d10278a99c3c3da56db64bec5d5ffe5cdd6` |
| EJ200_xp200 | `50175093bb626a28dec5bb3ce4cf137bdcb880640eef6f78c918a1d8c541bee4` | `6f26965063fbe15ba3daa2facd7d60674292f4f7b7f252d18f7da221c4ccd881` | analysis.root: `1c18fbb731c656ef37cf222dc6047bd7732675f83c025f253dd80c198f9dd6c9`; analysis.csv: `ad5ae4dec4e224a3b4915bae2bf452c1a68ac04319c75ad47dc6528a79c2f532`; analysis.meta.json: `a529a8e5e9e309b010cb4ff8db83f06c054fd19ee6dad2c686fcab267c0ad9af` |
| EJ200_xm200 | `a08779771e0bb6c10672cb263f527433c597a594080f5253b73c0c68d028e3f8` | `9c655cf495f4b3b10661645cf633e2275e744705ab3c2c4a1a1b7f7146e292cc` | analysis.root: `9f7da1e3241c5e08a9afec7569925f462a4dd625c6a90c068f8b5e64d99065d7`; analysis.csv: `022139f302d619ab77ddf8ac72f645af601d981422d8026e9012c9d36213080e`; analysis.meta.json: `b459d8bc1f2574c351bc3ed3444314ff71d3ce1b72181d472ed61c565bc8f428` |
| EJ200_xp500 | `d709e7acb92f67ff02987c3970b95b46d8b9b1086fd4735090e4ba5a8d8b4479` | `94f50d840881e9f65938cc9fda6d8c48cc0d229a8804f4e7e1105c60f2275f1c` | analysis.root: `aba98f0cee4b12e927f22dbe2a631fd70c1549fe8f85863338b4247b35c3784e`; analysis.csv: `6ebffcc85cc5fd4e9d397c64ec88e963e4baadc5231f2829684fb40bbfb25725`; analysis.meta.json: `891f78e29c6c8138789791d179a44b0bb2b4fe8fd17922df878e641a71be78f7` |
| EJ200_xm500 | `27b844e5aa9bf8ddf9981dee26df7796809453d0be5bed8f975798f3f69a654d` | `ea3fb68b288b47273a262d3fc3080d176ef6f256fe55ce0d110dec82e5145de8` | analysis.root: `ec313b21b652bb1fb75f5055fc6874bbcd4ae67c6cc158556cd2614aff5d2780`; analysis.csv: `4812c8a3b4084b399aded691769a2b7886da78b3af9e015f9026d5845aeafe81`; analysis.meta.json: `b8226eaaf1a50926bbabd0e07d98482f8624580d84983587c0936de2745bbade` |
| EJ200_xp650 | `63ed88dd00be500f28a06fe3da2a028e2078c9bc3df1d8fb40b4629028e3b0b2` | `eb6da3dafc135e12cc4d71cba74270e12785672811c642518af89ddf6ce5f7f0` | analysis.root: `fc21bb61682c81ba233e0255a068550049e6a65b165259f9196e52bdd8576305`; analysis.csv: `9712ee553bf053749e99e710b975aa8d82c7efb54fa30f4fe3b9b4133210c54f`; analysis.meta.json: `c46ce12e768f578b6b02a33f721d624d3033fefabbe1f741ba892d1a4401fd97` |
| EJ200_xm650 | `141a12a387af769dfdb46764e8997bb78c9531597c5669166738f07bb447756b` | `2ac53f708d4121b5fd246668545cb6bc889c0948b7c86f6bf2a1703b66734a9c` | analysis.root: `3a41d219038d90154748cba79e15e5f6ff31232350df8baeda708a4b970363c9`; analysis.csv: `46efd890821222ca8a18c2dd42504e6f1bfe417c2d41e16814b9228542426c76`; analysis.meta.json: `3e67181036f3171d6d713e954103ce42f1e85e40e0771bc8bfe1341c97695869` |
| EJ204_xp0 | `82b5048fdb6f32a224102cae7ba1ab5c1de0d7fb454f0a66cb910babc5a334f6` | `e133c4b0357a09ce73b7ae4b05332932562bf00f33e60a07e12b739ab074fc8a` | analysis.root: `1555452920108ff7c421e1b48558bd8ba4c9756b33fe6ddd71ebcfbb774baa41`; analysis.csv: `e02f3655784572676e8ce9743c9c64747daf1f455ad2a1768fab8663ad8312a6`; analysis.meta.json: `c3465ef10168561173554435704441fd269555d9f6ac6cc8ecba5e450128d3b5` |
| EJ204_xp200 | `88a0050cdd3510adb92d6683afcb4aace359739adbbaf10f9941411fc254fa05` | `c260d021f94f4dd7bcef84ff6a3fccaa559e4b94d25f0dcb862aa52e114ab96d` | analysis.root: `2f5e85df31403f6ca879a9bbdc614b56864c150f8ba083ee3d0ce1227f3ceaf6`; analysis.csv: `d245b0d6056ac64a7a61e399af5703b82324c3e62da6fd14b6641db1e53948ac`; analysis.meta.json: `7f1913772da4bdbedb8e79450af1d1d06c376ca96e06ddf61ba8957df0d01a81` |
| EJ204_xm200 | `3f93f533cf961523568b002be1685918f62b8697d3c9e2c8f036e3926da240e1` | `4f5502bd85749bd597827f840b13a3126cb228bd4ed686a636497887ca62d7da` | analysis.root: `f6163984d159a77cade955f8a2e335e7a1308b73da6f05d2fee20798e61c6cf9`; analysis.csv: `1cc8e80024c119c0a91dd5c8c0118b7ec28be3ab4ba14c0d65c88a5b70a75876`; analysis.meta.json: `78f6d9a483f9df050e05ebbc8dbc0cf86b6992d92462637ab91682bf0cad363b` |
| EJ204_xp500 | `6294e3ecbef3c5b01186bc195d9c60378a34db1ca3e8325642571f5fa73ae1ee` | `92a04967946ce844e1847342916077d8b272dd1338ec5056f1956d825605c86d` | analysis.root: `44c1c9da7d30bd81ee00b886205d9209e84763574a8c4d06b524d4718d6ab5a2`; analysis.csv: `52c330c3e41638a4b85601038d72c4234307c5987578e24833bd65434db802a7`; analysis.meta.json: `ee169239b7682d2a565871f149c22aa9d38b5098f2f540cf060aa37ad2d74a85` |
| EJ204_xm500 | `5536f0fa5c66cb3df9801b294bdd088e550b56f995b90a10ebab0903854c238b` | `49868a636ee5db658c80d91844ef6b2de2c1718e784d043bcda41a5e475eeb98` | analysis.root: `0c34539bf72e22cc34d51141add527ec9a473606bac3f087220dca50693318e4`; analysis.csv: `b12b07cb1081a7f531af14f5025da31c5d7bd645ad3d900fa1abb9b3e88a26c6`; analysis.meta.json: `e18a4b3ff588ba3f67a02f0e87cf76ed66c402b45907f9014128a65f84c8bf91` |
| EJ204_xp650 | `4cf0542238e33e05281980cadcb03e3a1e6a276f6137739017ec0dc54bdefd09` | `4c7c5f32e99f8c35c4a48afc55f3dfd76ccc7ba136edf557c472421ac25055cb` | analysis.root: `7c7ab395c9188e7d8cd0abf1cfeeaef388343efd44634a447883f320087415a3`; analysis.csv: `bcc87e19c57e981b1131baa4d9eea43b6b301cc48742c7dc40acba1d99c325b8`; analysis.meta.json: `237e9851138118f745012fdf152a60afcc652b97e8c3384a150b2bddbc9a8888` |
| EJ204_xm650 | `6c861a4d7e7ae0d2d7221b44447623fa81232f766f608ff214dd527797761b6e` | `b37849f96dadd0b7ff2ea1074f59848d83614daa6383d326a06cf938934ad34e` | analysis.root: `b508730275a2457de802f6c667dc68e3b6161f00908e13899c887f2344114902`; analysis.csv: `d3ae338eaee66856583c5d7765cfa1f67a9c89c20d8f89c4aeca1aef213314aa`; analysis.meta.json: `63dfe5a2e5d7497f3451c981e29d8f9764d8a3b82b4f9614e387a8b96c814075` |
| EJ230_xp0 | `d51a38f58a75b0cb176b6613dd96e76f32363793ed222622ef38ae49b6d16569` | `f0c7de32e9908d44407a3fc51d883ab0a252d7122d6a8d1daf37543eb66681f8` | analysis.root: `7ac8e10aad64aa2a04f58e9a0ca75a7124f6995664f78510b25762b7a826910f`; analysis.csv: `904d99d471ef28e9419ee779ff75bef69fc8509cd4c820421c628c9f93c087ca`; analysis.meta.json: `366c5efdc2b815cf6190c01fef0dbbb0e600fb4ce5c11c2b6d762f5b8e5314cf` |
| EJ230_xp200 | `27a2f234c91ac934c6ad43469d9296bd529567d3732b9f79a792819465f37419` | `5e8b156fef1751961fa397cc2d71814bcaca4d53affef2ef9432c2721f27d3f0` | analysis.root: `d64d1a35bb31779a983c203504ba30af83683220721dc40a73ff6e1a3325dc1e`; analysis.csv: `a35a6b9d59f92d6de06cf46343565d0286c414e3776a1822b0e683bbab5cbda6`; analysis.meta.json: `3a744365c989a90dc3ef604a62244d512be86dcbfd76607d2fb8b85e3527c60f` |
| EJ230_xm200 | `c1923b4e5a4695b15dec29720a6936bd04c0217dc0c27cb5dfd9f842729e339e` | `e50be5f2388c0b0bd6c11850b7f1f22c36a82bdc2e24dc5cba14890eeb52682a` | analysis.root: `3bceb8a224c088baebcc69178cc1fd9daa590fa6ff0d070d046588b816593fe7`; analysis.csv: `5915fed075d029e41a4d344d87fc8fe5490220ff7679ca2b85c5eeb9e78172ec`; analysis.meta.json: `b31e03c0fee59939950a20d64f1bda2ecbe19ce642f8a62c89a3e441b7d367c7` |
| EJ230_xp500 | `d0251884ff5c8593422a5b169214b47ee9ab9a1cc3a28c7daab0a1b0dd564f88` | `cbaa135e500713f3a37a67480c5909242f9bb70498ddd1f6f242671332ba1dce` | analysis.root: `767c210493be077184ff53855e22c8c6ddace88628046d19c8442195621e2c3d`; analysis.csv: `b332e53f09b07fa0fe988c21239e533d3a90755e4097a51180c83222ee67952e`; analysis.meta.json: `68c2b418acae2e7f8b8ec6400a9c9a1fd7f8b684869bcde8dec6c4d62296c3a9` |
| EJ230_xm500 | `80b0d6571bfbb2bb7043870fb679af06cd508fe7ef3d08441b118e8844b6c5be` | `01feef3c48de44baa18cd26a91ba5e2c09673bf7193f5a59f5b73eead9b19f04` | analysis.root: `dee59c4e3361b7d6b70f8243e47366c3b9ffd303435dcc8936433ecb020e03af`; analysis.csv: `85973001ffa82e5598d617a4153bb9e406fc104179610e5cd29eb5f0b6f7b0a3`; analysis.meta.json: `bc3019be7ca5e5b094e6cfe31247c7648e66d503221f0d4f08ed1465370dfd02` |
| EJ230_xp650 | `be493a4b7666dc41aa916cc74b552316b77b27b1656aeb8110de2a6e45700ac6` | `330de20731d69db1ddcfd0adceb4bf2c9c78b2d91f6ba825f288a63c3ef89cfa` | analysis.root: `77d102d937b6a17696345044f029dbd8e82e9554dc76fb91189c86dd7247bc45`; analysis.csv: `8a351ec455a311c546ae41646f85a0d113803c99fe8d4c543f1b07c3c23194cd`; analysis.meta.json: `39d8afd36eda73cda3e59a2dff1125f856f4616f284ba96f637a250cfabdeeab` |
| EJ230_xm650 | `ee43025ab5599e1245497c6f50a29ac33f3e9aa6af55fb18f6a4a3ef1f9cd6e7` | `16981f2e42ce03d1d326327f49cb0d2703a44a988aa74d2e382ce1e9269614f6` | analysis.root: `9252e25f7c33c426cf1389ec67ce3904808cfce68846515d5c9ab53f47431b6e`; analysis.csv: `d98f1cbb982ed669c485a97467d8f4c529d4eb9672d5425bc799c6858aa945f1`; analysis.meta.json: `152b8b999e49746f4d32268fff92fda26069bc1d0b201d87b9674feb8cbd0cd1` |

| Artifact | SHA256 |
| --- | --- |
| [EXEC34_HANDOFF_20260911.json](/home/rrios/ej200/docs/execution_logs/EXEC34_HANDOFF_20260911.json) | `ad60420f42472145f124e0cc13c2d29f916f15e0979a6b1aebaed433b4832183` |
| [manifest.jsonl](/home/rrios/exec34r_20260912/manifest.jsonl) | `b7fa1855722f99da9082494957f0bdd13f3242bf3bfc0a86bdd91ac9b16e91eb` |
| [input_index.json](/home/rrios/exec34c_20260912/input_index.json) | `1d17292241ee812db7f182bfd7e73fe50d0bea3bf47c4fc2544204fb72a90620` |
| [verified_inputs.json](/home/rrios/exec34c_20260912/verified_inputs.json) | `fc398172442bda0fd08ba0ef142f23739b3d80883103db80c29ac07ee8378331` |
| [execute_analysis.py](/home/rrios/exec34c_20260912/execute_analysis.py) | `26da537ec091fe74752477fc16775888a3d3d0b1ca4986718489c021e9ac324e` |
| [make_report.py](/home/rrios/exec34c_20260912/make_report.py) | `aa321b8854628c0d3c7882a7ae46db737bead9938b9a72ca7f68c8775fcd0a61` |
| [analysis_manifest.jsonl](/home/rrios/exec34c_20260912/analysis_manifest.jsonl) | `7837fb7496350d4a4a0cb1831e288762402de9c4bd3149b58a98d61b492dc079` |

Final checks: three sidecars and hashes verified per COMPLETE cell; ROOT/CSV row values agree; N_pe denominator 10000; all (b)−(c) values finite; both END fits and normalizations present; canonical parity and frozen state hashes correct; physics input provenance matches the valid manifest; imported primitives unchanged; git diff --check passes; Git clean. No simulations, timeout, BLUE, electronics, historical timing comparison, combined arms, deck changes, merge or push occurred.

**FINAL GATE — STOP.**
