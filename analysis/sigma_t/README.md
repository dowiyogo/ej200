# sigma_t — unchanged upstream pipelines and EXEC_33 decision gate

**Status: inventory/import complete; estimator choice pending. No pilot has been run.** TOP and END are independent arms. This package does not apply BLUE, combine their resolutions, or adopt a historical sigma as a reference.

- [INVENTORY.md](INVENTORY.md): scripts, execution order, dependencies, inputs/outputs and path constraints.
- [ESTIMATORS.md](ESTIMATORS.md): exact definitions, source lines, reporting procedures (a/b/c), EXEC_18/19 conflict and electronics versions.
- [provenance/import_manifest.json](provenance/import_manifest.json): byte-identical MSI source import, commit and hashes.
- [provenance/related_import_manifest.json](provenance/related_import_manifest.json): unchanged END/EXEC_12T/optimization/electronics sources and historical CSV.
- [upstream/analysis_core/README.md](upstream/analysis_core/README.md): original upstream documentation, preserved verbatim.

## Mandatory decisions before any pilot

René must explicitly choose all three:

1. **TOP photon-index range**, and the applicable grouping from the inventory (TOP_SUM4 is the merged stream of the four highest-count TOP channels at each position; it is not an absolute four-PE waveform threshold).
2. **TOP reporting procedure:** (a) one fixed global index, including how/where it is chosen; (b) per-position minimum on the same sample, labeled selection-biased; or (c) index chosen on training events and evaluated on independent events, with split and selection rule. No c-style TOP selector exists in the imported source; bootstrap and leave-one-position-out position calibration are not substitutes.
3. **END parameters:** fixed SUM4 map and first-crossing reduction; pulse shape/rise/fall, amplitude threshold, fit variant (`FitCore` or `tbmirror::FitPeakSeeded`), any NPE/ToT/walk cuts, jitter/electronics choice, and whether the reported END width is sigma(ΔT_LR) or sigma(ΔT_LR)/sqrt(2). The live provisional response is rise=0.5 ns, fall=5 ns, threshold=4 PE-equivalent, no injected SPTR, no walk/ToT correction; it is not automatically approved.

These choices are deliberately **unset** in [decision_pending.json](decision_pending.json). The gate in the EXEC_33 request is the reason the selected two-arm pilot command is not executed or labeled validated here. Implementing a new selector or silently changing upstream constants before the choice would violate the unchanged-import requirement.

## Runtime

On MSI, use `/usr/bin/python3.12` with ROOT/PyROOT 6.36.10. Unqualified `python3` is 3.9.25 and fails against that ROOT build. Python 3.12, ROOT, NumPy, uproot, SciPy, Matplotlib and PyYAML imports and the primary CLI `--help` were checked successfully. Pandas is used by auxiliary scripts. See the inventory for the observed versions.

```bash
ssh -p 9022 -o BatchMode=yes -o ConnectTimeout=5 reriosto@localhost 'echo OK'
```

Use a committed checkout/archive on MSI, not edits through an SSH pipe. The source imports retain upstream absolute paths. Paths below refer to the t0minidaq worktree; when using MSI, set `SIGMA_REPO` to the actual committed checkout location there.

## Exact native commands: ROOT → TOP sigma with uncertainties

This is the **existing fixed-index-curve entry point**, not an approved choice of (a/b/c). It computes each configured N separately; it does not choose the winning N.

```bash
SIGMA_REPO=/home/rrios/ej200_exec33_20260911
SIGMA_PACKAGE="$SIGMA_REPO/analysis/sigma_t"
SIGMA_CONFIG=/absolute/path/to/reviewed_top_config.yaml
SIGMA_TOP_OUT=/absolute/path/to/new/top_output
python3.12 "$SIGMA_PACKAGE/upstream/analysis_core/timing_fit_pipeline.py" \
  --config "$SIGMA_CONFIG" --materials EJ-204 --outdir "$SIGMA_TOP_OUT"
```

`reviewed_top_config.yaml` is a new configuration, not an edit to imported source. It must make every active parameter explicit:

```yaml
EXEC_TAG: EXEC_33
ROLLBACK_TAG: pre-exec33-20260911
MATERIALS:
  EJ-204:
    branch: REVIEWED_SIMULATION_REF
    expected_sha: REVIEWED_SIMULATION_SHORT_SHA
    sslg4_code: OPSC-101
    repo_path: /absolute/path/to/simulation/repository
    subdir: REVIEWED_INPUT_SUBDIRECTORY
    file_prefix: REVIEWED_INPUT_PREFIX
    sha_tag: REVIEWED_INPUT_SHA_TAG
INPUT_ROOT_BASE: /absolute/path/to/input/root/directory
TTREE_NAME: sipm_hits
SCAN_POSITIONS_MM: [0]
GROUPS: [TOP_SUM4]
N_VALUES: null  # mandatory gate choice; no default is silently adopted
REQUIRED_BRANCHES: [event_id, face_type, global_id, time_ns, gun_x_mm]
FIT_WINDOW_SIGMAS: 2.0
FIT_OPTIONS: R Q S 0
MIN_EVENTS_FOR_FIT: 30
CHI2_NDF_WARN: 3.0
CORE_NOT_GAUSSIAN_FRACTION: 0.30
EFFICIENCY_FLOOR: 0.05
N_BOOTSTRAP: 300
RANDOM_SEED: 20260618
BINNING_STRATEGY: sqrt_n
OUTPUT_DIR: /absolute/path/to/new/top_output
REPR_POSITIONS_MM: [0]
QA2_MU_MIN_NS: 0.10
QA2_MU_MAX_NS: 25.00
QA2_SIGMA_MAX_NS: 5.00
HOOK_WALK: false
HOOK_SPTR: false
```

These fit settings show the imported YAML defaults, not a new estimator selection. With approved N_VALUES, each input ROOT must match the exact legacy builder:

```text
INPUT_ROOT_BASE/subdir/file_prefix_xp0mm_300ev_sha_tag.root
```

The `_300ev_` token is hard-coded at `timing_fit_pipeline.py:74`; it does not read N from the filename. **Do not relabel a new N=10000 pilot as a historical 300-event run.** A separately versioned path/metadata adapter is needed after the gate for native `photon_hits_run000.root` naming and an explicit generated-event denominator. This constraint is documented rather than hidden by a misleading symlink/rename. Likewise QA-0 checks a branch tip: the pilot metadata must independently identify the simulation binary/source, not assume that a currently checked-out analysis branch produced an old ROOT.

Outputs: per-group results CSV contains `sigma_fit_ps`, `sigma_fit_err_ps`, `bootstrap_err_ps`, `chi2_ndf`, `fit_status`, event count and efficiency. `sigma_fit_ps` is sigma_TOP of the selected order-statistic timestamp distribution, with the ROOT covariance error and bootstrap error recorded separately. The `.root` sidecar stores histograms, functions and graphs; graph error bars use the maximum of fit and bootstrap error (`lib/sidecar.py:75`). Metadata uses the upstream naming convention `*_metadata.json`; the pilot contract requires a normalized `.meta.json` in a separate adapter. No min across N is taken by this entry point.

## Exact native commands: ROOT → END sigma with uncertainty

Two live END fit variants exist. Their waveform/map constants are compile-time values, **not CLI arguments**. The commands below document their actual interfaces; they were not run on an EXEC_33 sample.

**Congruent SUM4 fit**, accepting the native sequential Geant4 ROOT filename and N:

```bash
SIGMA_REPO=/home/rrios/ej200_exec33_20260911
SIGMA_INPUT=/absolute/path/containing/photon_hits_run000.root
SIGMA_END_OUT=/absolute/path/to/new/end_congruent_output
root -l -b -q \
  "$SIGMA_REPO/analysis/sigma_t/upstream/related/420addf/analysis/congruent_sum4_timing.C(\"$SIGMA_INPUT\",\"$SIGMA_END_OUT\",10000)"
```

`summary.csv` contains `sigma_lr_end_ps ± sigma_lr_end_err_ps`, and their division by sqrt(2) as `sigma_single_end_ps ± sigma_single_end_err_ps`; check `fit_used_end`, `chi2_ndf_end`, `n_eff_end`. The full legacy main also computes a TOP split-waveform diagnostic and plots a historical 88-ps line. **That full main is not an acceptable EXEC_33 pilot entry point:** its TOP definition is different and the pilot forbids historical comparisons. A later adapter must call the unchanged END primitives only, not this plotting main.

**Test-beam mirror fit**, the entry point corresponding to the peak-seeded fit description:

```bash
SIGMA_REPO=/home/rrios/ej200_exec33_20260911
SIGMA_AUDIT=/absolute/path/to/reviewed_end_input.csv
SIGMA_END_OUT=/absolute/path/to/new/end_mirror_output
root -l -b -q \
  "$SIGMA_REPO/analysis/sigma_t/upstream/related/420addf/analysis/tb_mirror_sigma_vs_x.C(\"$SIGMA_AUDIT\",\"$SIGMA_END_OUT\")"
```

The exact positional CSV input schema consumed by `ReadAudit` is:

```csv
x_mm,unused1,root_path,unused3,unused4,unused5,status
0,,/absolute/path/to/photon_hits_run000.root,,,,OK
```

The main has a fixed capacity of 10000 event IDs. `sigma_vs_x_End.csv` contains `sigma_LR_ps`, `err_sigma_LR_ps`, `sigma_single_ps`, `err_sigma_single_ps`, chi-square diagnostics, and alternative NPE-cut widths. This fit uses a peak±2 ns window, whereas `FitCore` iterates a ±2-sigma window. These are different fits. Neither full macro provides the complete pilot `.root/.csv/.meta.json` contract by itself.

## End-to-end pilot contract after the gate

The exact native commands above establish what the imported source can execute. They **do not yet constitute a validated, selected two-arm pilot command**. No such existing command was found: the TOP selector choice is pending, (c) is absent, ROOT naming/denominator metadata need an adapter, and the native END mains include extra historical output or lack complete sidecars. Reporting an already-reproducible pilot command now would be false.

After the explicit decisions, a separate commit must provide only the required orchestration/configuration around the preserved primitives, with any new selection procedure clearly distinguished from imported code. Its command must record the approved choices and the following simulation stage, then produce the two independent arm results:

```text
EndTop; N_TOP=70; EJ-204/OPSC-101; mu- 1 GeV vertical; x=0 mm
N=10000; simulation seeds 26092601 8349041; workers selected from S1–S4
EJ200_ENABLE_DIAGNOSTICS=OFF; simulation SiPM jitter=0 ns unless explicitly revised
simulation command: <verified-binary> -m <archived-pilot.mac>
input: photon_hits_run000.root / sipm_hits
TOP stage: approved grouping + photon-index range + reporting procedure + fit/error model
END stage: approved SUM4/map + pulse + leading-edge threshold + reduction + fit/normalization
outputs: separate sigma_TOP ± error and sigma_END ± error; selected/discarded events,
         fit validity/quality, and .root/.csv/.meta.json for each result
```

All NPE values must include generated zero-hit events; historical `np.unique(event_id)` denominators are explicitly insufficient to establish that denominator. The pilot metadata must carry simulation and analysis commits, Geant4 version, N/seeds, material/code, geometry/TOP layout, position/workers/eventModulo, PDE path/hash and every stage command. Fit error is not a substitute for selection uncertainty, and a zero-valued failure field is not a resolution.

**Stop here.** No selected estimator, pilot sigma, BLUE correction or 21-cell campaign is authorized by this README. René's three decisions are required first; the full grid later requires new approval even after a successful pilot.
