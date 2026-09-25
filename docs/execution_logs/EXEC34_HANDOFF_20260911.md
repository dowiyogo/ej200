# EXEC34A — pilot report and EXEC34B handoff

Updated UTC 2026-09-12T04:44:34.959487+00:00. **G-P: PASS. READY FOR EXEC_34B — GRID NOT STARTED.**

The five specified pilot conditions pass. This does not certify the Gaussian model or close the diagnostic warnings below. EXEC34A stops here; no grid cell has been launched.

## Findings that must travel with the handoff

- Canonical TOP(c): **3.090539 ± 0.057955 ps**, N=1, efficiency=100.0%, 5000 generated EVAL events. Its fit converged (`fit_status=0`) but chi2/ndf=13.513377 (ndf=5); the Gaussian core model fits poorly. No numerical sigma reference or extra chi-square acceptance gate was introduced.
- The internal bias is **+0.010220 ps**, whereas (b)-(c) is **-0.071538 ps**: opposite signs, magnitude ratio 7.000. Under the requested diagnostic this is a possible implementation-defect warning, not a physical result. It remains unresolved; neither orientation nor result was chosen to hide it. Mechanical leakage checks passed, which does not turn this mismatch into agreement.
- Both orientations choose the boundary N=1. N=0 is not a photon order statistic and is outside the prescribed range. The two-sided curvature cannot be reported; the available N=2 neighbor is retained below. The winner is not established as a physical optimum.

## Repositories, provenance and closed decisions

Worktree `/home/rrios/ej200_exec33_20260911`, branch `diag/exec33-20260911`. Initial HEAD `1a9bba5a98ae70623dd4afc4c21b33584a1805b2`; final HEAD `b406143da55595e0e5f2ab81ed3d5b47c9e2cf58`; status **clean**. Rollback tag `pre-exec34a-20260911` points to the initial HEAD and was created before the first edit. New TOP(c) orchestration is its own commit `1b5b1b7174f3a7ca76fd340f4daa49fc33b1fd54`; the analyzed orchestration revision is `0f3db2d09daaa7731e3b4cc54309f3fae8a50a2b`. Analysis execution records pipeline HEAD `cab884a9d5fb4393592e7a8ac61243ef4483bf54`; the final HEAD additionally includes the campaign runner/tests. Imported primitives from 30eba3b/fde610f are unchanged. Main remains 420addf. No push, merge, reset, rebase, branch/tag deletion, deck edit, BLUE, electronics injection or TOP/END combination occurred.

The simulation used **EndTop, 70 TOP, EJ-204/OPSC-101, vertical mu- 1 GeV, x=0 mm, N_generated=10000, seeds 26092601/8349041, diagnostics OFF, one process, 24 workers, eventModulo=1**. Geant4 11.4.0; binary SHA256 `b70b61489d0a8992186f547f9765a7da90dbc2c14889d0ffb7877e00ea91ef25`. The recorded 70 TOP placement IDs are 16..85. Simulation jitter is zero. [geometry_evidence.json](/home/rrios/exec34a_20260911/geometry_evidence.json) preserves this evidence.

**INTRINSIC timing resolution — electronics not included.** There are at least four live electronics versions, unresolved FWHM/sigma ambiguity, and SPTR_PROVENANCE.md does not establish sqrt(kN) propagation for order statistics. No historical timing value is used or compared.

## Literal G-P criterion and automatic result

> El piloto PASA si y solo si se cumplen las cinco condiciones:
> 1. sigma_END finito, con `fit_used_end` verdadero, chi2/ndf registrado
>    y n_eff por encima del mínimo del ajuste.
> 2. sigma_TOP finito con `fit_status` válido en al menos 10 de los 20
>    índices N, en AMBAS mitades de la partición.
> 3. Eficiencia por encima del `EFFICIENCY_FLOOR` configurado para el N
>    ganador de (c).
> 4. Los tres sidecars completos y el comando de extremo a extremo
>    reproducible desde el README.
> 5. La diferencia (b) − (c) es finita y calculable.
> 
> NO hay criterio sobre el VALOR de sigma_t: no existe referencia válida
> contra la cual juzgarlo. Un valor inesperado no es un fallo.

| condition | PASS/FAIL | evidence |
| --- | --- | --- |
| G-P.1 | PASS | {"sigma_END_ps": 119.90401871014716, "fit_used_end": true, "chi2_ndf": 0.6626944280583688, "n_eff": 10000, "minimum": 20} |
| G-P.2 | PASS | {"TRAIN_even": 20, "EVAL_odd": 20} |
| G-P.3 | PASS | {"efficiency": 1.0, "efficiency_floor": 0.05, "winning_N": 1} |
| G-P.4 | PASS | {"analysis.root_hash": true, "analysis.csv_hash": true, "root_csv_row_count": true, "root_csv_values": true, "all_generated_end_events": true, "metadata_complete": true, "simulation_provenance": true, "stage_exit_codes": true, "analysis_scripts_match": true, "readme_end_to_end_command": true, "native_input_name": true, "separate_arm_results": true, "finite_uncertainties": true} |
| G-P.5 | PASS | {"b_minus_c_ps": -0.0715375416082944} |

[gate.json](/home/rrios/exec34a_20260911/pilot/analysis/gate.json) contains the machine evaluation and sidecar hashes. G-P.1 uses FitCore, not its RMS fallback. G-P.2 counts 20/20 valid indices in both canonical halves. G-P.3 uses EVAL efficiency strictly above 0.05. G-P.4 opens ROOT/CSV/JSON, compares numeric rows, verifies hashes and the README command. G-P.5 checks the finite difference, not its sign or magnitude.

## TOP estimators and all requested split diagnostics

TOP_SUM4 ranks aggregate TOP counts in TRAIN only. Canonical TRAIN: even event_id; EVAL: odd event_id, each 5000 generated events including zero-hit events. Canonical channel ranking: `[50, 51, 52, 49]`; reverse ranking: `[51, 50, 52, 49]`. Each ranking is frozen for its EVAL. All N=1..20 curves are in [analysis.csv](/home/rrios/exec34a_20260911/pilot/analysis/analysis.csv) and the ROOT/metadata sidecars. For every N the files record sigma, covariance/bootstrap error, efficiency, fit status, n_eff, window counts and fit quality. Sigma is conditional on having at least N selected hits; fewer-hit events are omitted from the width but retained in the efficiency denominator.

The adapter calls preserved ranking, compute_tN and fit_core_gaussian primitives. It freezes per-N seeds (peak, amplitude, MAD), histogram bounds/binning and fit windows as well as channel IDs and winning N. EVAL fits its measurement parameters with that construction; it does not learn a new window, ranking or N. A scoped dependency injection changes neither imported files nor function bodies. The recorded frozen-model fingerprints are unchanged by EVAL. No walk calibration exists. Bootstrap seed 20260618, 300 replicates, sqrt_n binning, ROOT options R Q S 0, minimum 30 and +/-2 MAD-sigma windows are explicit in the README and metadata.

(a) reports all fixed-N curves separately for TRAIN/EVAL. (b) independently ranks/selects on ALL 10000 events and reports the same-sample argmin, labeled **biased by selection**. (c) uses the TRAIN winner on EVAL and is primary. Reported errors are max(covariance error, fixed-estimator bootstrap error), with both components retained. They are conditional errors, not uncertainty over repeated training selections.

| Quantity | sigma ± conditional error [ps] | N | efficiency | chi2/ndf | n_eff |
| --- | --- | --- | --- | --- | --- |
| TRAIN even at its winner | 3.100759 ± 0.054975 | 1 | 1.0 | 16.94740919081031 | 5000 |
| EVAL odd, canonical (c) | 3.090539 ± 0.057955 | 1 | 1.0 | 13.513377186917902 | 5000 |
| TRAIN odd at its winner | 2.957733 ± 0.054493 | 1 | 1.0 | 21.019297041150622 | 5000 |
| EVAL even, reverse (c) | 3.015649 ± 0.053834 | 1 | 1.0 | 18.40439628434015 | 5000 |
| ALL, (b), biased by selection | 3.019001 ± 0.038284 | 1 | 1.0 | 30.640540229143106 | 10000 |

Canonical EVAL has 0 missing-hit discards; 3970 accepted times lie inside the fixed core window and 1030 outside. `n_eff=5000` denotes eligible timestamps, not the count of histogram entries in the fitted core. The per-N histograms/functions are preserved.

**Partition symmetry:** sigma(c) TRAIN even minus sigma(c) TRAIN odd = +0.074890 ± 0.079101 ps, diagnostic significance 0.947. Compatible under the declared 2-sigma diagnostic convention. No average is formed. EVAL event sets are disjoint, but selection cross-uses the halves; quadrature conditional-error significance does not estimate cross-selection covariance.

**Internal bias:** TRAIN sigma at its winner minus EVAL sigma at the same winner = +0.010220 ± 0.079882 ps. **(b)-(c)** = -0.071538 ± 0.069458 ps (quadrature diagnostic errors; overlapping samples/covariance unmodeled). The sign mismatch is explicit. A factor-ten definition of same order was recorded before analysis; the magnitude ratio 7.000 meets it. These full-sample and split estimates overlap and are not statistically independent. Their small differences relative to their conditional errors do not identify the cause of the mismatch, so it is not declared resolved.

**Local minimum diagnostic:**

| TRAIN parity | N | sigma ± error [ps] | efficiency |
| --- | --- | --- | --- |
| even | 1 | 3.100759 ± 0.054975 | 1.0 |
| even | 2 | 3.373123 ± 0.076179 | 1.0 |
| odd | 1 | 2.957733 ± 0.054493 | 1.0 |
| odd | 2 | 3.376525 ± 0.062501 | 1.0 |

Canonical N=2 minus N=1 = 0.272364 ± 0.093944 ps (2.899 diagnostic sigma). `curve_flat=false` under the declared one-error adjacency rule; `boundary_minimum=true`, two-sided curvature=null. This does not establish an optimum below or beyond the allowed range.

## END, separate from TOP

Left SUM4 clusters {0,1,2,3}/{4,5,6,7}; right {8,9,10,11}/{12,13,14,15}. First finite cluster crossing per end; accept only both finite. Normalized difference-of-exponentials, rise=0.5 ns, fall=5 ns, threshold=4 PE-equivalents of summed amplitude. No SPTR/walk correction/ToT cut. The new bridge calls only the preserved END primitives, never the imported plotting mains.

| Fit | sigma(DeltaT_LR) ± error [ps] | sigma(DeltaT_LR)/sqrt(2) ± error [ps] | chi2/ndf | n_eff |
| --- | --- | --- | --- | --- |
| FitCore | 119.904019 ± 1.452208 | 84.784945 ± 1.026866 | 0.6626944280583688 | 10000 |
| FitPeakSeeded | 118.205413 ± 0.836748 | 83.583849 ± 0.591670 | 0.9740445346929113 | 10000 |

FitCore: four iterative +/-2 sigma fits, `fit_used_end=true`; FitPeakSeeded: peak +/-2 ns window. The imported mirror struct does not expose fit status; its positive width comes from its successful-fit branch. Both accept 10000 events, no nonfinite-end discards. FitCore's eligible count is 10000; its primitive does not return the exact final-window entry count. All END event times and exact input histogram definitions are saved to permit inspection.

Systematic sign **FitPeakSeeded - FitCore** = **-1.698606 ps** for sigma(DeltaT_LR), or -1.201096 ps for the derived normalization. The fit errors are correlated because the events are the same; no independent systematic error is asserted.

The derived normalization assumes **"equal and statistically independent timing contributions from the two ends"**. That assumption has **NOT been verified here**. No sigma_TOP/sigma_END combination is formed.

## Performance, memory and prepared grid

N_pe/end = **398.142100 ± 1.006536 SEM**, using all 10000 generated events. Left total 3981965, right total 3980877. Generated events with zero hits anywhere: 0. No historical timing baseline is used.

Simulation total wall **597.12 s**, startup **0.459826 s**, maximum RSS **407.184 MiB**. GNU time -v directly measures process peak RSS; startup is process launch to the master marker after geometry/physics initialization and before beamOn. The run includes ROOT output. Small local code/test tasks overlapped the simulation; no other optical campaign ran. Analysis wall **213.559 s**, direct getrusage peak **2.684 GiB**.

RAM total 124.406 GiB; available immediately before pilot 110.982 GiB. Reserve 20% of available RAM and budget twice the larger simulation/analysis peak = 5.367 GiB per complete independent one-worker pipeline. This permits **16 processes by memory**, so the prepared grid concurrency is **16 × 1 worker**. Considering only simulation, twice the measured 24-worker process RSS would permit 111 independent processes by memory; the CPU cap is 24. No private per-thread RSS is invented. Using the larger full-pipeline peak protects the concurrent analysis phase. Live available RAM can reduce the cap at launch.

EJ-200=**OPSC-100**, directly from [DetectorConstruction.cc](/home/rrios/ej200_exec33_20260911/src/DetectorConstruction.cc) line 116: `if (code == "EJ200") code = "OPSC-100";`. EJ-204=OPSC-101 at 115; EJ-230=OPSC-106 at 117. The campaign includes all three materials and x=0,±200,±500,±650 mm, 10000 events/cell, one worker, same seeds and even-TRAIN orientation.

Campaign configuration: [campaign.json](/home/rrios/exec34b_20260911/campaign.json). Manifest: [manifest.jsonl](/home/rrios/exec34b_20260911/manifest.jsonl), append-only JSONL. There are exactly 21 PENDING entries and no cell output directory. Runner validation was executed without --execute and launched nothing. Scripts/primitives and binary are hash-pinned. A process lock prevents duplicate runners; PASS cells require sidecar verification to skip. Incomplete/inconsistent stages are preserved and refused. On failure, no further jobs are submitted; already-running jobs finish with recorded evidence. No physics retry is automatic.

Start in EXEC34B (not executed here):

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/grid.py --config /home/rrios/exec34b_20260911/campaign.json --execute
```

Resume:

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/grid.py --config /home/rrios/exec34b_20260911/campaign.json --execute --resume
```

## Reproduction and validation

[README](/home/rrios/ej200_exec33_20260911/analysis/sigma_t/README.md) contains the actual one-command chain and analysis-only replay command. The native ROOT is not renamed or linked to a historical filename. Current replay verification: **PASS**. [replay_verification.json](/home/rrios/exec34a_20260911/replay_verification.json) will contain exact numeric comparison of the same-ROOT replay. ROOT container bytes may differ because of timestamps/UUIDs; equality is checked on the numeric CSV and frozen models.

Seven orchestration tests passed: event-order invariance, poisoned-EVAL isolation, explicit generated denominators, frozen histogram/seeds, fixed END reduction and missing ends, append-only planning/resume, and failed-pilot refusal. The available five CTests passed. git diff --check is clean; imported-source diff against original imports is empty; main unchanged. No external electronics reference was used to tune these results.

## Artifact hashes and commands

Machine-readable checkpoint: [EXEC34_HANDOFF_20260911.json](/home/rrios/ej200/docs/execution_logs/EXEC34_HANDOFF_20260911.json). It includes the complete frozen model, reverse orientation, all diagnostic values, both END normalizations, parameters, campaign script hashes and commands. PDE `/home/rrios/ej200_exec33_20260911/data/sipm/AFBR-S4N66P024M_pde.txt`, SHA256 `6360bd80eabf1ade77ea0dd80f8e56ed99b6ea8dcbc2ce94bf5269108f2b28e7`. Raw ROOT `/home/rrios/exec34a_20260911/pilot/photon_hits_run000.root`, SHA256 `99efbf9b85963cf9117405024272d6247bf6eee067b83cf828f791f3d8a9feb3`.

| Sidecar | SHA256 |
| --- | --- |
| [analysis.root](/home/rrios/exec34a_20260911/pilot/analysis/analysis.root) | fc9ab26e26773a0e6a5511f5feeed9be2fcd40fd4deaa530170718f312220231 |
| [analysis.csv](/home/rrios/exec34a_20260911/pilot/analysis/analysis.csv) | e02f3655784572676e8ce9743c9c64747daf1f455ad2a1768fab8663ad8312a6 |
| [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) | 86bab9ccd42f856720eb0be0a342f6495edd6b5696f6575f2f3044a9bf73f84b |


| Start UTC | End UTC | Exit | Command / internal stage | Evidence |
| --- | --- | --- | --- | --- |
| 2026-09-11T22:29:50.212682+00:00 | 2026-09-11T22:29:50.216175+00:00 | 0 | git add analysis/sigma_t/orchestration/run_simulation.py | [command_001.log](/home/rrios/exec34a_20260911/audit/command_001.log) |
| 2026-09-11T22:29:50.341584+00:00 | 2026-09-11T22:29:50.349293+00:00 | 0 | git commit -m 'analysis(sigma_t): add audited single-cell simulation runner for EXEC34' | [command_002.log](/home/rrios/exec34a_20260911/audit/command_002.log) |
| 2026-09-11T22:29:50.404107+00:00 | 2026-09-11T22:39:47.531059+00:00 | 0 | /usr/bin/time -v -o /home/rrios/exec34a_20260911/pilot/resource_usage.txt /home/rrios/exec33_20260911/build_off/ej200_bar_sim -m /home/rrios/exec34a_20260911/pilot/run.mac | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:31:59.431088+00:00 | 2026-09-11T22:32:01.008246+00:00 | 0 | python -m unittest discover -s analysis/sigma_t/orchestration -p test_top_split.py -v | [command_003.log](/home/rrios/exec34a_20260911/audit/command_003.log) |
| 2026-09-11T22:33:05.445458+00:00 | 2026-09-11T22:33:05.450901+00:00 | 0 | git add analysis/sigma_t/orchestration/top_split.py analysis/sigma_t/orchestration/test_top_split.py | [command_004.log](/home/rrios/exec34a_20260911/audit/command_004.log) |
| 2026-09-11T22:33:05.765705+00:00 | 2026-09-11T22:33:05.782822+00:00 | 0 | git commit -m 'analysis(sigma_t): implement new TOP procedure c orchestration with frozen TRAIN state' -m 'Reuse unchanged channel-ranking, order-statistic and Gaussian-fit primitives. Freeze channels, N, fit seeds, histogram axes and windows; parity split has explicit generated-event denominators. Add leakage, ordering and missing-event regression checks.' | [command_005.log](/home/rrios/exec34a_20260911/audit/command_005.log) |
| 2026-09-11T22:39:38.671279+00:00 | 2026-09-11T22:39:38.682107+00:00 | 0 | git add analysis/sigma_t/orchestration analysis/sigma_t/README.md analysis/sigma_t/decision_pending.json | [command_006.log](/home/rrios/exec34a_20260911/audit/command_006.log) |
| 2026-09-11T22:39:38.973449+00:00 | 2026-09-11T22:39:38.986919+00:00 | 0 | git commit -m 'analysis(sigma_t): integrate native ROOT adapter, preserved END fits and automatic pilot gate' | [command_007.log](/home/rrios/exec34a_20260911/audit/command_007.log) |
| 2026-09-11T22:40:01.417338+00:00 | 2026-09-11T22:40:08.999360+00:00 | 0 | ctest --test-dir /home/rrios/exec33_20260911/build_off --output-on-failure | [command_008.log](/home/rrios/exec34a_20260911/audit/command_008.log) |
| 2026-09-11T22:40:01.852493+00:00 | 2026-09-11T22:43:36.239023+00:00 | 0 | /usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/analyze.py --simulation /home/rrios/exec34a_20260911/pilot --output /home/rrios/exec34a_20260911/pilot/analysis | [cell_invocation_2026-09-11T224001.852493_0000.json](/home/rrios/exec34a_20260911/pilot/cell_invocation_2026-09-11T224001.852493_0000.json) |
| 2026-09-11T22:40:02.529159+00:00 | 2026-09-11T22:43:36.046358+00:00 | 0 | /usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/analyze.py --simulation /home/rrios/exec34a_20260911/pilot --output /home/rrios/exec34a_20260911/pilot/analysis | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:40:04.201114+00:00 | 2026-09-11T22:40:12.761032+00:00 | 0 | analyze.py: uproot sipm_hits iterator; explicit generated IDs 0..9999 | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:40:12.761055+00:00 | 2026-09-11T22:41:15.475155+00:00 | 0 | top_split.learn -> freeze -> top_split.evaluate; no EVAL optimization | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:40:55.808556+00:00 | 2026-09-11T22:40:57.318311+00:00 | 0 | python -m unittest discover -s analysis/sigma_t/orchestration -p 'test_*.py' -v | [command_009.log](/home/rrios/exec34a_20260911/audit/command_009.log) |
| 2026-09-11T22:41:15.475162+00:00 | 2026-09-11T22:42:18.152733+00:00 | 0 | top_split.learn -> freeze -> top_split.evaluate; no EVAL optimization | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:41:23.661189+00:00 | 2026-09-11T22:41:23.664980+00:00 | 0 | git add analysis/sigma_t/orchestration/test_end_bridge.py analysis/sigma_t/orchestration/gate.py | [command_010.log](/home/rrios/exec34a_20260911/audit/command_010.log) |
| 2026-09-11T22:41:23.786431+00:00 | 2026-09-11T22:41:23.793616+00:00 | 0 | git commit -m 'test(sigma_t): verify END acceptance and ROOT CSV gate consistency' | [command_011.log](/home/rrios/exec34a_20260911/audit/command_011.log) |
| 2026-09-11T22:42:00.531811+00:00 | 2026-09-11T22:42:00.535268+00:00 | 0 | git diff --check | [command_012.log](/home/rrios/exec34a_20260911/audit/command_012.log) |
| 2026-09-11T22:43:27.054923+00:00 | 2026-09-11T22:43:35.989905+00:00 | 0 | exec34::analyze_end -> LeadingEdgeTime/Earliest -> FitCore and tbmirror::FitPeakSeeded | [analysis.meta.json](/home/rrios/exec34a_20260911/pilot/analysis/analysis.meta.json) |
| 2026-09-11T22:43:36.239531+00:00 | 2026-09-11T22:43:36.368621+00:00 | 0 | /usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/gate.py --analysis /home/rrios/exec34a_20260911/pilot/analysis | [cell_invocation_2026-09-11T224336.239531_0000.json](/home/rrios/exec34a_20260911/pilot/cell_invocation_2026-09-11T224336.239531_0000.json) |
| 2026-09-11T22:43:47.125596+00:00 | 2026-09-11T22:43:47.129159+00:00 | 0 | git add analysis/sigma_t/orchestration/run_cell.py | [command_013.log](/home/rrios/exec34a_20260911/audit/command_013.log) |
| 2026-09-11T22:43:47.254770+00:00 | 2026-09-11T22:43:47.262407+00:00 | 0 | git commit -m 'fix(sigma_t): verify worker and RNG configuration when resuming cells' | [command_014.log](/home/rrios/exec34a_20260911/audit/command_014.log) |
| 2026-09-11T22:46:02.946765+00:00 | 2026-09-11T22:46:04.444095+00:00 | 0 | python -m unittest discover -s analysis/sigma_t/orchestration -p 'test_*.py' -v | [command_015.log](/home/rrios/exec34a_20260911/audit/command_015.log) |
| 2026-09-11T22:46:34.187682+00:00 | 2026-09-11T22:46:34.191422+00:00 | 0 | git add analysis/sigma_t/orchestration/grid.py analysis/sigma_t/orchestration/test_grid.py analysis/sigma_t/README.md | [command_016.log](/home/rrios/exec34a_20260911/audit/command_016.log) |
| 2026-09-11T22:46:34.312801+00:00 | 2026-09-11T22:46:34.319643+00:00 | 0 | git commit -m 'analysis(sigma_t): prepare memory-bounded EXEC34B runner with verified resume' | [command_017.log](/home/rrios/exec34a_20260911/audit/command_017.log) |
| 2026-09-11T22:46:34.492158+00:00 | 2026-09-11T22:50:10.440443+00:00 | 0 | /usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/analyze.py --simulation /home/rrios/exec34a_20260911/pilot --output /home/rrios/exec34a_20260911/pilot_analysis_replay | [cell_invocation_2026-09-11T224634.492158_0000.json](/home/rrios/exec34a_20260911/pilot/cell_invocation_2026-09-11T224634.492158_0000.json) |
| 2026-09-11T22:47:34.795303+00:00 | 2026-09-11T22:47:34.817325+00:00 | 0 | python /home/rrios/exec34a_20260911/audit/prepare_campaign.py | [command_018.log](/home/rrios/exec34a_20260911/audit/command_018.log) |
| 2026-09-11T22:47:34.865438+00:00 | 2026-09-11T22:47:34.888876+00:00 | 0 | python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/grid.py --config /home/rrios/exec34b_20260911/campaign.json | [command_019.log](/home/rrios/exec34a_20260911/audit/command_019.log) |
| 2026-09-11T22:50:10.440769+00:00 | 2026-09-11T22:50:10.568799+00:00 | 0 | /usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/gate.py --analysis /home/rrios/exec34a_20260911/pilot_analysis_replay | [cell_invocation_2026-09-11T225010.440769_0000.json](/home/rrios/exec34a_20260911/pilot/cell_invocation_2026-09-11T225010.440769_0000.json) |
| 2026-09-11T22:50:39.453671+00:00 | 2026-09-11T22:50:39.496675+00:00 | 0 | python /home/rrios/exec34a_20260911/audit/write_handoff.py | [command_020.log](/home/rrios/exec34a_20260911/audit/command_020.log) |
| 2026-09-12T04:43:44.179349+00:00 | 2026-09-12T04:43:45.881747+00:00 | 1 | python /home/rrios/exec34a_20260911/audit/verify_final.py | [command_021.log](/home/rrios/exec34a_20260911/audit/command_021.log) |
| 2026-09-12T04:43:45.910477+00:00 | 2026-09-12T04:43:45.953889+00:00 | 0 | python /home/rrios/exec34a_20260911/audit/write_handoff.py | [command_022.log](/home/rrios/exec34a_20260911/audit/command_022.log) |
| 2026-09-12T04:43:53.717155+00:00 | 2026-09-12T04:43:55.433153+00:00 | 0 | python /home/rrios/exec34a_20260911/audit/verify_final.py | [command_023.log](/home/rrios/exec34a_20260911/audit/command_023.log) |

Internal analysis stages execute inside the recorded analyze.py command; they are function calls, not separate shell commands. Initial read-only inspection/tag creation is recorded in initial_state.json and the Git tag object; its exact shell transcript is retained by the session. No missing timestamp is invented.

Final commit history:

```text
3ae8214 analysis(sigma_t): add audited single-cell simulation runner for EXEC34
1b5b1b7 analysis(sigma_t): implement new TOP procedure c orchestration with frozen TRAIN state
0f3db2d analysis(sigma_t): integrate native ROOT adapter, preserved END fits and automatic pilot gate
cab884a test(sigma_t): verify END acceptance and ROOT CSV gate consistency
2683713 fix(sigma_t): verify worker and RNG configuration when resuming cells
b406143 analysis(sigma_t): prepare memory-bounded EXEC34B runner with verified resume
```

**READY FOR EXEC_34B — GRID NOT STARTED.**
