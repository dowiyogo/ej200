# EXEC34R detached simulation launcher

Prepared for René; the real campaign has **not been launched**. This is a
simulation-only launcher for 21 cells, six concurrent processes with four
Geant4 workers each, eventModulo 1, 10000 generated events per cell. It neither
runs the sigma_t analysis nor modifies its estimators.

## Commands on t0minidaq

The prepared campaign is `/home/rrios/exec34r_20260912`. From any directory:

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/detached_grid.py dry-run
```

**Launch command, printed for René to execute (not executed in EXEC34R):**

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/detached_grid.py launch
```

The same command resumes: it validates and skips `.DONE` cells and retries
unfinished or failed cells. It returns after starting a detached driver.
Use `--directory /absolute/campaign/path` to select a different prepared
campaign. The default is fixed to the EXEC34R directory above.

**One-line status command:**

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/detached_grid.py status
```

It reports complete, running, pending and failed counts and cell IDs. Status
reads `.DONE`, ROOT size and `state.json`, verifying live PID/start ticks/boot
identity through `/proc`; it does not depend on the driver's memory. Launch
additionally rehashes every DONE ROOT before skipping it. Invalid DONE markers
fail closed rather than silently rerunning changed data.

## Session independence and timeout policy

`launch` calls `Popen(start_new_session=True)`, which performs `setsid`. The
driver has no controlling SSH terminal; stdin is `/dev/null`, stdout/stderr
append to `grid_logs/driver.log`, and the command returns immediately. No
`nohup`, shell background operator or persistent interactive session is needed.
The six simulations inherit that independent session. Required Geant4/library
environment values are captured at preparation; incoming G4 overrides from a
new shell are discarded and the verified PDE directory is pinned explicitly.

There is **no timeout**, no deadline parameter and no automatic per-cell retry
within one invocation. Failed cells do not block the other queued cells. If a
timeout is added later, its baseline must be measured with **four workers and
the same event count**. The 597.12 s, 24-worker pilot is never a valid timeout
baseline for four-worker or one-worker simulations. EXEC33 S2 measured 358.38 s
at four workers for 2000 events and estimated 1790.34 s for 10000; the estimate
is not a new four-worker, 10000-event measurement and is not used as a timeout.

A hung cell occupies one of six slots. A process-wide `flock` prevents duplicate
drivers; its descriptor is inherited by simulation children. If the driver
dies, surviving children retain the lock until exit, preventing a new driver
from duplicating their work. Once all orphan children exit, relaunch retries
unmarked cells in fresh attempt directories. A stale RUNNING record with no
live driver/child is reported as failed. Host reboot stops processes, while
persisted states and markers remain available for resume.

## Macros, evidence and preflight

`prepare` was executed once. It read the existing EXEC34B `campaign.json` and
`cells/*/run.mac` and copied the macros into the new campaign. It changed only
the numeric tokens for `/run/numberOfThreads` and `/run/eventModulo`. All other
bytes, including seeds, material, geometry and N, are unchanged. Repeating
`prepare` refuses an existing configuration; use `dry-run`, `status` or `launch`.
Each launch checks copied macros against their source bytes and SHA256 values.
The original campaign manifest and partial outputs remain intact.

Example actual diff (identical thread-only change in all 21 macros):

```diff
--- /home/rrios/exec34b_20260911/cells/EJ200_xp0/run.mac
+++ /home/rrios/exec34r_20260912/cells/EJ200_xp0/run.mac
@@ -1,4 +1,4 @@
-/run/numberOfThreads 1
+/run/numberOfThreads 4
 /run/eventModulo 1
 /random/setSeeds 26092601 8349041
 /control/verbose 0
```

`EXEC34_HANDOFF_20260911.json` explicitly records `EJ200_OPSC_CODE: OPSC-100`.
The seven EJ-200 cells are included. The tested `NOT_FOUND` path excludes those
seven cells and prints their IDs; an absent/ambiguous code aborts preparation.
Evidence files and their hashes are pinned in the campaign configuration.

Every cell's `simulation.meta.json` records four workers, eventModulo, seeds,
binary/macro/PDE hashes and `parallelism_provenance`: EXEC33 tested 1, 4, 12 and
24 workers with exact zero difference in N_pe/end and event counts relative to
S1. The four-worker efficiency was 0.956966. This records the benchmark's
measured scope; it is not presented as a new invariance test of the 21 cells.

Preflight runs synchronously before detachment and again in the driver:

- Pilot ROOT bytes × 21 × 1.3 against current destination filesystem free bytes.
- Six-process expected RSS from the four-worker EXEC33 S2 measurement
  (129780 KiB/process), plus a conservative gate: six times
  `[2 × max(four-worker RSS, 10000-event pilot simulation RSS) + 256 MiB]`,
  leaving 20% of currently available memory in reserve. Pilot RSS is used only
  as a memory bound, never for a timeout. ROOT validation streams 32 MB chunks.
- All included macros readable and exactly matching authorized transformations.
- Binary executable and SHA256 equal to the verified EXEC34A/34B binary, with
  `EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF` in its CMake cache; runtime data available.
- Handoff, scaling, event reproducibility, source config and PDE hashes intact;
  Python `uproot` dependency importable.

The real dry-run passed: projected disk 96,880,108,279 bytes versus
1,512,596,361,216 available; expected six-process RSS 797,368,320 bytes;
conservative six-cell budget 6,734,168,064 bytes versus 118,947,266,560 available.
These are snapshots; launch repeats the measurements and aborts on any failure.

Full 21-cell joblist, preflight and command templates:
`/home/rrios/exec34r_20260912/dry_run.txt`.
Materials are EJ200/OPSC-100, EJ204/OPSC-101 and EJ230/OPSC-106; each has
x = 0, +200, -200, +500, -500, +650, -650 mm. No cells were excluded.

## Completion, logs and analysis handoff

Each attempt uses `cells/<cell_id>/attempts/<unique-id>/`, containing an exact
copy of the prepared macro, `stdout.log`, `simulation.meta.json` and the native
`photon_hits_run000.root`. Failed attempts are retained. Public logs append to
`grid_logs/<cell_id>.log`, including attempt separators, exact commands and
terminal status; driver messages append to `grid_logs/driver.log`.

A cell gets an atomic, fsynced `cells/<cell_id>/.DONE` only after exit code zero
and successful output validation. The native ROOT has **only a hit tree, not a
generated-event ledger**. The validator therefore requires the master simulator
summary `Events run: 10000`, reads every ROOT branch/basket, verifies every hit
ID belongs to 0..9999, and cross-checks the three face photon totals and the
number of hit-bearing events against the log. Zero-hit events are valid and
are not silently removed from the generated denominator. It checks finite hit
times, records the ROOT SHA256 and size, and verifies the reported PDE path.
The real pilot ROOT passed this reader: 55,640,955 hits, 10000 generated events,
10000 hit-bearing events, SHA256
`99efbf9b85963cf9117405024272d6247bf6eee067b83cf828f791f3d8a9feb3`.
Evidence: `/home/rrios/exec34r_20260912/pilot_read_only_validation.json`.

The append-only, flock-protected and fsynced `manifest.jsonl` uses the existing
fields `utc`, `cell_id`, `status`, `output`. `output` remains the directory that
contains the run artifacts: here it is the unique attempt directory. Added
fields include exact shell command and argv, workers, eventModulo, exit code,
wall_s, hashes, start/end timestamps, cell_directory, log and provenance.
`RUNNING` and `FAILED` retain their usual meaning. `SIMULATION_COMPLETE` is a
new status, explicitly distinct from EXEC34B `PASS` (which included analysis).
No analysis PASS is fabricated. Missing output hashes are JSON null; existing
partial ROOTs are hashed on normal failed termination. An abrupt driver crash
leaves disk artifacts and stale state for the next status/resume operation.

For subsequent analysis use each successful manifest row's `output` as the
`--simulation` directory for the existing `analyze.py`; it contains compatible
simulation metadata and native ROOT naming. The launcher does not invoke that
analysis, imported plotting programs, BLUE or any TOP/END combination.

## Tests and delivery

```bash
cd /home/rrios/ej200_exec33_20260911
/usr/bin/python -m unittest discover -s analysis/sigma_t/orchestration -p test_detached_grid.py -v
```

Tests use synthetic ROOT files and mocked `Popen`, never Geant4. They cover
byte-preserving macro changes, a non-executing and non-writing dry-run,
NOT_FOUND exclusion, disk/RAM/diagnostics/binary/macro preflight failures,
ROOT/log corruption and zero-hit events, DONE hashing and resume, duplicate
locks, detachment arguments, preserved failed attempts and continued scheduling
with six slots. Tests use the archived local handoff/macros as fixtures.
All 12 launcher tests passed; the complete orchestration suite passed 20/20.
The recorded suite output is `/home/rrios/exec34r_20260912/unit_tests.log`.

Rollback tag: `pre-exec34r-20260912` points to the pre-edit commit `498d14f`.
Prepared macros/configuration/evidence live outside the worktree; the launcher,
tests and this README are delivered in their own commit. No campaign was
launched, and no push, merge, reset, rebase, branch/tag removal or main edit
was performed.
