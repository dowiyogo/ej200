# EXEC_42 Phase A — production-instrumentation grid preparation

**Status: PREPARED; MANUAL GRID NOT LAUNCHED.**

| Gate | Requirement | Measurement | Status |
|---|---|---|---|
| A1 | Diagnostics-OFF and diagnostics-ON builds; non-transport tests only | Both builds succeeded; `sslg4_properties_check` and `readout_config_check` passed in both | PASS |
| G-P.1 | Reproduce 397.1525 pe/end with exact-zero difference and identical per-event counts | 397.1525 pe/end, difference 0; zero differing events for produced and all three detected counts | PASS |
| G-P.2 | Reproduce V1, V2-H1 and matched V5 within declared precision | V1 and V2-H1 exact; V5 matched exact; expected V5 differs by `1.11e-16` | PASS |
| G-P.3 | Measure actual ROOT size and compare with 987.9 MB estimate and D1 | 1,545,540,402 bytes; 1.565× estimate; 2.183× D1 | MEASURED; EXEC_41 estimate undershot |
| A3 preflight | 21 readable macros, pinned binary/resources, disk and RAM sufficient, fresh outputs | 21 pending; 210.97 GB projected versus 1.426 TB available; 3.916 GB conservative RAM budget versus 118.77 GB available | PASS |

## A1 — source and build verification

The worktree began at `b35ee84acadef12c93506e0720580c91f901fcbf` on `diag/exec40-20260913`. The rollback tag `pre-exec42-20260913` points to that commit. Phase-A preparation is committed as `cbd99efeb158845cf2577b9e43a08356646eada9`. Main remains clean at `420addf0fd6029d5b2f0e235f472a8ae47f31fac`.

The two existing preregistrations remain unchanged:

- `analysis/validation/EXEC40_PREREGISTRATION.md`: SHA256 `cd13d28db1811e5c74d9b17d60d1be7aa0971f02ff8768d48528126602fb5531`.
- `analysis/validation/EXEC41_PREREGISTRATION.md`: SHA256 `933433a2882c4753301c675913df1f37db614e4afc8316ae87df6c122021d418`.

The production schema is implemented in `src/FirstEncounterObservation.cc`: `pre_volume_id` and `post_volume_id` use the frozen integer dictionary; `cos_incidence` is a ROOT `float` populated from `G4float`; and raw `boundary_status`, `normal_valid`, `normal_orientation_valid`, `normal_norm`, and `energy_eV` are inside `EJ200_ENABLE_DIAGNOSTICS`. These are executable code paths, not a storage proposal.

Both binaries were rebuilt from the worktree:

```bash
cmake --build /home/rrios/exec40_20260913/build_off -j8
cmake --build /home/rrios/exec40_20260913/build_on -j8
```

Only the non-transport tests were selected:

```bash
ctest --test-dir /home/rrios/exec40_20260913/build_off --output-on-failure -R '^(sslg4_properties_check|readout_config_check)$'
ctest --test-dir /home/rrios/exec40_20260913/build_on --output-on-failure -R '^(sslg4_properties_check|readout_config_check)$'
```

All four selected test executions passed. No unfiltered CTest command or smoke simulation was run.

## A2 — production-schema verification cell

Before launch, the wall-time estimate was approximately seven minutes, based on the prior EJ-204 EndTop N=2,000 four-worker measurement. This was within the authorized 20-minute single-cell limit.

The executed command was:

```bash
cd /home/rrios/exec42_20260913/a2_validation && /usr/bin/time -v -o /home/rrios/exec42_20260913/a2_validation/resource_usage.txt /home/rrios/exec40_20260913/build_off/ej200_bar_sim -m /home/rrios/exec42_20260913/a2_validation/run.mac
```

The macro is byte-identical to the validated EXEC_40 macro, SHA256 `250bf0c9e2b2e423181a76aa769c1065b64c95d9ec711987f762c1d58d353822`: EndTop, N_TOP=70, EJ-204/OPSC-101, vertical 1 GeV muon at x=0, N=2,000, seeds 26092601 and 8349041, four workers, `eventModulo=1`, zero jitter, diagnostics OFF. The binary SHA256 is `ec0d4498aad912728e5a3134f357198527aeeca7c29da944ca838ece4eecf1b4`.

The process exited zero after **404.63 s (6:44.63)** and used **192,131,072 bytes** maximum RSS. The ROOT SHA256 is `8608285fea6f12a5675f7821ebb5c39a45cd53535ac814df59b69b98036ba59d`.

The exact validation command was:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python /home/rrios/ej200_exec40_20260913/analysis/validation/check_exec42_a2.py --candidate /home/rrios/exec42_20260913/a2_validation/photon_hits_run000.root --reference /home/rrios/exec40_20260913/cell_validated/photon_hits_run000.root --exec41 /home/rrios/exec41_20260913 --d1 /home/rrios/ej200_exec26_20260909/build_exec29_20260910/cells/D1/photon_hits_run000.root --out /home/rrios/exec42_20260913/a2_validation
```

### G-P.1

The candidate totals are 40,039,535 scintillation photons, 42,213,806 optical photons, 794,729 left-END detections, 793,881 right-END detections and 9,510,815 TOP detections. Every per-event value in those five fields is identical to EXEC_40. The resulting END yield is exactly 397.1525 pe/end, with exact-zero difference.

### G-P.2

| Observable | Candidate | EXEC_40/41 reference | Absolute difference | Allowed one-SE precision | Status |
|---|---:|---:|---:|---:|---|
| V1, all deposited energy | 1.000256914479 | 1.000256914479 | 0 | 0.000159481 | PASS |
| V1, non-optical denominator | 1.016318589503 | 1.016318589503 | 0 | 0.000164469 | PASS |
| V2-H1 escape | 0.256898235925 | 0.256898235925 | 0 | 0.000070760 | PASS |
| V5 matched detection/incident | 0.612397869359 | 0.612397869359 | 0 | 0.000113950 | PASS |
| V5 incident-spectrum PDE | 0.612462602316 | 0.612462602316 | `1.11e-16` | 0.000026652 | PASS |

All seven V2 face counts are also identical. The production storage change did not alter the simulated physical results.

### G-P.3

The actual ROOT is **1,545,540,402 bytes**, compared with the EXEC_41 estimate of 987,879,523 bytes. The estimate undershot by 56.45%. The measured tree bytes are:

| Tree | Compressed bytes |
|---|---:|
| `sipm_hits` | 707,601,425 |
| `event_observables` | 109,160 |
| `first_bar_encounters` | 826,487,743 |
| `sipm_event_counts` | 10,072,869 |

The discrepancy comes from using an uproot-rewritten encounter artifact to estimate output produced by the Geant4 ROOT analysis writer. The physical/schema validation passes, but future storage projections must use the measured 1,545,540,402-byte file. The new file is 2.183× the D1 reference size of 707,921,457 bytes.

A2 artifacts are under `/home/rrios/exec42_20260913/a2_validation/`: `a2_results.json`, `a2_validation.{csv,root,meta.json}`, `invocation.meta.json`, the ROOT, macro, log, and resource-use record.

## A3 — 21-cell grid preparation

The prepared campaign directory is `/home/rrios/exec42_20260913/grid`. It contains 21 cells: EJ-200/OPSC-100, EJ-204/OPSC-101 and EJ-230/OPSC-106 at x = 0, ±200, ±500 and ±650 mm. Every cell specifies EndTop, N_TOP=70, a vertical 1 GeV muon, N=10,000, seeds 26092601 and 8349041, diagnostics OFF, four workers and `eventModulo=1`.

The launcher is the detached, resumable EXEC_34R implementation introduced at `d99f994`, imported with its tests and extended only to verify local `sslg4` links and material-resource hashes, scale disk projections between different pilot/target event counts, and require fresh attempt directories at first launch. Twelve synthetic launcher tests passed; none invokes Geant4.

The copied macros differ from the prior source macros only in the authorized thread count. Example:

```diff
--- /home/rrios/exec34b_20260911/cells/EJ204_xp0/run.mac
+++ /home/rrios/exec42_20260913/grid/cells/EJ204_xp0/run.mac
@@
-/run/numberOfThreads 1
+/run/numberOfThreads 4
 /run/eventModulo 1
 /random/setSeeds 26092601 8349041
```

The EJ-200 handoff explicitly names `OPSC-100`; no material code was inferred. All OPSC-100, OPSC-101 and OPSC-106 macro/data resources are present and hashed. The PDE path is `/home/rrios/ej200_exec40_20260913/data/sipm/AFBR-S4N66P024M_pde.txt`, SHA256 `6360bd80eabf1ade77ea0dd80f8e56ed99b6ea8dcbc2ce94bf5269108f2b28e7`. All 21 local `sslg4` links resolve to the build resource tree. No attempt/output directory exists yet.

### Concurrency decision

A2 measured 192,131,072 bytes RSS per four-worker process. The launcher budgets twice that measurement plus 256 MiB for ROOT verification, or 652,697,600 bytes per active cell. Six processes therefore require 1,152,786,432 bytes at measured RSS and 3,916,185,600 bytes under the conservative budget. At preflight, 118,767,988,736 bytes were available. An 80%-of-available-memory calculation would allow 145 such conservative process budgets.

The selected concurrency is **six processes × four Geant4 workers**. Memory is not limiting at six; six is retained because it is the already validated EXEC_34R concurrency. The choice is therefore bounded by the measured memory check rather than derived from the number of cores. EXEC_33 supplies the physics invariance evidence: exact-zero event differences for 1, 4, 12 and 24 workers, with 95.7% efficiency at four workers.

Linear scaling of A2 gives approximately 2,023 s per N=10,000 cell. Four waves are needed for 21 cells at concurrency six, giving a rough campaign estimate of 8,093 s or 2.25 hours.

### Dry preflight output

The dry run returned:

```json
{
  "status": "PASS",
  "disk_projected_bytes": 210966264873,
  "disk_projection_rule": "pilot_ROOT_size * target_N/pilot_N * job_count * 1.3",
  "pilot_ROOT_size": 1545540402,
  "pilot_N_generated": 2000,
  "target_N_generated": 10000,
  "job_count": 21,
  "disk_available_bytes": 1425956577280,
  "expected_six_process_rss_bytes": 1152786432,
  "conservative_six_cell_budget_bytes": 3916185600,
  "memory_available_bytes": 118767988736,
  "macros_readable": 21,
  "fresh_attempt_output_directories": 21,
  "diagnostics": false,
  "binary_sha256": "ec0d4498aad912728e5a3134f357198527aeeca7c29da944ca838ece4eecf1b4",
  "timeout_s": null
}
```

The complete dry-run output is `/home/rrios/exec42_20260913/grid/dry_run.txt`, SHA256 `eff61d4c65d1355ad5015292123224d82e5d5c6785253becd5978b34da48b90d`. Its full joblist is:

```text
EJ200_xp0 EJ200_xp200 EJ200_xm200 EJ200_xp500 EJ200_xm500 EJ200_xp650 EJ200_xm650
EJ204_xp0 EJ204_xp200 EJ204_xm200 EJ204_xp500 EJ204_xm500 EJ204_xp650 EJ204_xm650
EJ230_xp0 EJ230_xp200 EJ230_xm200 EJ230_xp500 EJ230_xm500 EJ230_xp650 EJ230_xm650
JOBLIST 21; complete=0; running=0; pending=21; failed=0; EXECUTED=0
```

Campaign configuration SHA256: `4a06fdd4bcc616b0e0827cbc189276496ca52c5820f500cb544a28d56f3e884e`.

## Commands for René

Launch:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python /home/rrios/ej200_exec40_20260913/analysis/sigma_t/orchestration/detached_grid.py launch --directory /home/rrios/exec42_20260913/grid
```

Status:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python /home/rrios/ej200_exec40_20260913/analysis/sigma_t/orchestration/detached_grid.py status --directory /home/rrios/exec42_20260913/grid
```

Resume — intentionally identical to launch:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python /home/rrios/ej200_exec40_20260913/analysis/sigma_t/orchestration/detached_grid.py launch --directory /home/rrios/exec42_20260913/grid
```

The launch action has not been invoked. Phase B must not begin until René explicitly confirms that all 21 cells completed. No push, merge, deck-source edit, reset, rebase, BLUE combination, or electronics injection was performed.
