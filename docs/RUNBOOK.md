# Runbook — corrected physics integration

The branch `integrate/physics-fixes-20260911` contains the guard and reflector fixes. Main remains at `8349041` until manual promotion. Physics provenance, measurements and qualifications are recorded in [DEFECTS_AND_CORRECTIONS.md](DEFECTS_AND_CORRECTIONS.md).

## Build and diagnostics

Use separate build directories so the diagnostic configuration is explicit. On the validation host:

```bash
cmake -S /home/rrios/ej200 -B /home/rrios/exec31_20260911/build_off \
  -DCMAKE_BUILD_TYPE=Release \
  -DGeant4_DIR=/home/tdship/opt/geant4-v11.4.0-install/lib64/cmake/Geant4 \
  -DEJ200_ENABLE_DIAGNOSTICS=OFF
cmake --build /home/rrios/exec31_20260911/build_off -j8
ctest --test-dir /home/rrios/exec31_20260911/build_off \
  -R 'sslg4_properties_check|readout_config_check' --output-on-failure
```

`EJ200_ENABLE_DIAGNOSTICS` defaults to **OFF** if omitted. ON compiles the optical-boundary census, terminal-fate TrackingAction, explicit-kill annotations, legacy atomic counters and RNG snapshots. OFF removes those sources/actions and optional hot-loop census updates; ordinary scintillation/event totals and hit ROOT output remain available. EXEC_32 keeps the required `Bar -> SiPM (entering)` boundary-encounter tally and its master-only reset/output active with either flag, so the physics smoke assertion remains observable. This historical tally counts encounters toward a SiPM, not unique transmitted photons. The optical boundary lookup required by the corrected escape guard stays enabled in either build. Build ON in a separate directory using `-DEJ200_ENABLE_DIAGNOSTICS=ON`.

## Macro configuration

Set the readout and material before initialization. `/det/readout End` selects END-only (16 active END SiPMs, no TOP; four lateral wrap panels). `/det/readout EndTop` selects 16 END plus **70 TOP** SiPMs (three lateral wrap panels; +Y carries TOP). These counts and the absence of an experimental +Y panel with TOP are checked by `readout_config_check`.

The D1 regression macro is:

```text
/run/numberOfThreads 1
/random/setSeeds 26092601 8349041
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout EndTop
/det/scintillator OPSC-101
/sipm/model AFBR-S4N66P024M
/run/initialize
/sipm/jitterSigma 0 ns
/gun/particle mu-
/gun/energy 1 GeV
/muon/angle 0
/muon/gunX 0 mm
/run/beamOn 2000
```

`/sipm/model AFBR-S4N66P024M` selects the validated Broadcom model (`Broadcom` is an alias). `$EJ200_DATA_DIR` overrides the compiled default PDE root; it is an environment variable, not a UI command. For example, `EJ200_DATA_DIR=/home/rrios/ej200/data` resolves `sipm/AFBR-S4N66P024M_pde.txt` below that directory. Record the resolved path and SHA256. Do not assume the executable hash identifies the loaded PDE or SSLG4 data.

## Output directory and runtime dependencies

Use a fresh output directory per job and make its `sslg4` path resolve to the configured build's SSLG4 runtime tree. Example for a **new** job name:

```bash
mkdir /home/rrios/exec31_new_job
ln -s /home/rrios/exec31_20260911/build_off/sslg4 /home/rrios/exec31_new_job/sslg4
cp /home/rrios/exec31_20260911/cell_off/run.mac /home/rrios/exec31_new_job/run.mac
cd /home/rrios/exec31_new_job
/home/rrios/exec31_20260911/build_off/ej200_bar_sim -m /home/rrios/exec31_new_job/run.mac
```

The example is not an instruction to rerun a completed dataset. Preserve existing outputs. Missing relative SSLG4 macros already caused a startup failure in EXEC_26 (child SIGABRT/−6, wrapper exit 250; see defect record). The Geant4 OPSC-101 `mat031` mass-fraction warning is also present in the historical reference; EXEC_31 does not change those data.

Archive the source commit, binary hash, macro, seeds, N, worker count, explicit/default eventModulo, CWD, runtime data hashes, UTC start/end, command/exit status and stdout alongside `photon_hits_run000.root`. With ON, retain boundary/terminal CSVs and RNG snapshots. For derived figures retain `.root`, `.csv` and `.meta.json` sidecars. END pe/end means the event average `(left+right)/2`, including zero-hit events; report its event-level SEM.

The legacy wrapping log still prints obsolete BarSkin labels and can show `Reflector R: 0`; use the checked panel geometry and measured census for interpretation, not that label. This limitation is inherited and documented.

## Reproducible scheduling

Prefer **one worker per job**, as used by D1 and both EXEC_31 regression cells. The previously validated alternative is **12 workers with `/run/eventModulo 1`**, set before beamOn. The EXEC_30 G0 gate passed against the one-worker END-only D0 control, including zero differences in per-event END/TOP counts; this is a validation for that configuration, not a guarantee for arbitrary thread counts or Geant4 versions. Evidence: [G0](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/audit/G0.json), [csv](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/event_reproducibility.csv) · [root](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/event_reproducibility.root) · [meta.json](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/event_reproducibility.meta.json).

The recorded elapsed-time ratio was **1.468017×** (D0 3341.036238 s, V1 2275.883586 s; N=2000, seeds 26092601 8349041; [G0 provenance](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/audit/G0.json)). Shared census mutexes are a plausible source of contention from inspection, but this was not a profiler measurement or isolated scaling benchmark. OFF removes that diagnostic work; EXEC_31 does not claim a new general MT speedup.

For position scans, prefer independent one-worker jobs, one fresh output directory per position, rather than MT within each job. Record each job's seeds explicitly. Do not use statistical agreement to replace exact event-level reproducibility when testing a diagnostic toggle.

No timing-resolution campaign or deck-number replacement is part of this integration. New production results must cite the corrected commit and runtime data, rather than inherit `talk_v6` provenance.
