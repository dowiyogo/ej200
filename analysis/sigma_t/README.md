# EXEC34 intrinsic timing pipeline — two independent arms

The imported source under `upstream/` is unchanged. New orchestration lives in
`orchestration/`. EXEC34A closes the estimator choices; it stops after the pilot
and its handoff. It never starts the grid automatically.

## End-to-end command

On t0minidaq the checked runtime is `/usr/bin/python` 3.9.25 with PyROOT 6.40.02,
NumPy 1.23.5, uproot 5.6.9, SciPy 1.13.1, PyYAML and Matplotlib. On MSI use
`python3.12` explicitly, with a compatible PyROOT; do not use unqualified python3.
The EXEC34 pilot runs locally. No dependency installation is required.

From a new output directory, this one command runs simulation, reads the native
ROOT, applies the versioned TOP/END primitives, writes the three sidecars and
evaluates G-P. Exit 0 means PASS; exit 34 means G-P FAIL. Other nonzero codes
mean a failed implementation/runtime stage. No failed stage starts a grid.

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/run_cell.py \
  --binary /home/rrios/exec33_20260911/build_off/ej200_bar_sim \
  --output /home/rrios/exec34a_20260911/pilot \
  --workers 24 --material EJ-204 --opsc OPSC-101 --x 0
```

The pilot simulation was launched first with `run_simulation.py` using exactly
these defaults, while the adapter was implemented. To finish or resume that
existing cell use the same command with `--resume`. Resume verifies the native
ROOT hash, the cell parameters and analysis script hashes, and does not simulate
again when its simulation stage is complete. It refuses to overwrite an existing
incomplete stage. To replay only analysis on the SAME ROOT, use a fresh directory:

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/run_cell.py \
  --binary /home/rrios/exec33_20260911/build_off/ej200_bar_sim \
  --output /home/rrios/exec34a_20260911/pilot \
  --resume --analysis-output /home/rrios/exec34a_20260911/pilot_analysis_replay
```

The binary comes from simulation `420addf`, Release, Geant4 11.4.0,
`EJ200_ENABLE_DIAGNOSTICS=OFF`. Rebuild if necessary:

```bash
cmake -S /home/rrios/ej200_exec33_20260911 -B /absolute/new/build-directory \
  -DCMAKE_BUILD_TYPE=Release -DEJ200_ENABLE_DIAGNOSTICS=OFF \
  -DGeant4_DIR=/home/tdship/opt/geant4-v11.4.0-install/lib64/cmake/Geant4
cmake --build /absolute/new/build-directory -j24
```

The runner fixes EndTop, 70 TOP, vertical mu- 1 GeV, 10000 generated events,
seeds 26092601/8349041, eventModulo=1, zero injected jitter. It accepts native
`photon_hits_run000.root`; no historical `_300ev_` alias is created. All event
IDs must be in [0,10000); generated events with no hits remain in denominators.

## TOP parameters and frozen selection

TOP_SUM4 means four TOP channels with highest aggregate hit counts. All N=1..20
are evaluated. Canonical TRAIN is even event_id, EVAL odd; both denominators
are 5000. A second execution of the split orchestration reverses these parities
on the same ROOT. Neither orientation depends on ROOT hit order.

In TRAIN, preserved `best_gids_by_count`, `compute_tN` and
`fit_core_gaussian` learn the channels, curve and minimum. Lowest N breaks an
exact sigma tie; only finite, positive, status=0 fits can win. EFFICIENCY_FLOOR
is evaluated by G-P, not used to change the optimum. For every N, freeze
TRAIN-derived peak/amplitude/MAD seeds, histogram axes/bin count, fit window
and configuration. Bootstrap seed is fixed, not learned. EVAL receives this
construction through a scoped dependency-injection wrapper; upstream files
and function bodies are unchanged. EVAL fits its Gaussian parameters to
measure its sigma; it does not reconstruct seeds/windows, rank channels or
optimize N. The complete model and a before/after fingerprint are persisted.
No walk or other data-driven calibration exists in this configuration.

Fit configuration: sqrt_n binning, window peak +/-2 MAD-sigma, ROOT options
`R Q S 0`, MIN_EVENTS_FOR_FIT=30, N_BOOTSTRAP=300, RANDOM_SEED=20260618,
EFFICIENCY_FLOOR=0.05. These are the imported defaults. The reported uncertainty
is max(ROOT covariance error, fixed-estimator bootstrap error); both components
are retained. It is conditional on the learned TRAIN estimator, not the
unconditional uncertainty of the training-selection procedure.

(a) All twenty fixed-N results are reported separately for TRAIN/EVAL in each
orientation. EVAL curves are diagnostic only. (b) The full 10000-event sample
has its own ranking and same-sample minimum, explicitly `biased by selection`.
(c) Primary result is EVAL sigma at the TRAIN winner using all frozen state.
`b_minus_c` means full-sample (b) minus canonical EVAL (c). The internal diagnostic
is TRAIN sigma at its winner minus EVAL sigma at the SAME winner. Full-sample
(b) overlaps both halves, so these two bias estimates are correlated. Any sign
or order mismatch is explicitly reported as the requested implementation
warning, not reinterpreted as physics or added to the five gate conditions.

Every sigma row includes efficiency, n_eff, missing-hit discards, fit_status,
chi2/ndf, fit and bootstrap error, window count and histogram overflow. Sigma(N)
is conditional: **events with fewer than N hits are omitted**. The generated
partition denominator is never inferred from np.unique(event_id).

Symmetry diagnostic: canonical minus reversed sigma(c), quadrature conditional
error and significance; <=2 sigma is the stated compatibility convention.
No average is formed. Adjacent-point variation <= one quadrature uncertainty
marks a flat/noise-dominated minimum; a boundary has only its available neighbor,
and no N=0 or N=21 result is invented. Cross-index covariance is not estimated.
These diagnostics are not numerical sigma acceptance gates.

## END parameters

`end_bridge.h` includes preserved sources and calls ONLY END primitives. No
imported plotting main, BLUE calculation or historical comparison is executed.
Fixed left clusters {0,1,2,3}/{4,5,6,7}; right {8,9,10,11}/{12,13,14,15}.
Normalized difference-of-exponential pulse, rise=.5 ns, fall=5 ns; threshold
4 PE-equivalents of summed amplitude. Earliest finite cluster crossing per
end, and both ends must be finite. No SPTR, walk correction or ToT cut.

Primary: `FitCore` from congruent_sum4_timing.C:110, four iterative +/-2 sigma
fits, minimum 20. RMS fallback is preserved but `fit_used_end=false` fails G-P.
Systematic: `tbmirror::FitPeakSeeded`, tb_mirror_sigma_vs_x.C:105, peak +/-2 ns;
its returned struct lacks fit status (zero failure fields are documented).
Systematic sign is **FitPeakSeeded - FitCore** on sigma(DeltaT_LR).

Always report sigma(DeltaT_LR) and sigma(DeltaT_LR)/sqrt(2), both with their
fit errors. The latter assumes "equal and statistically independent timing
contributions from the two ends"; that assumption is NOT verified here.

**INTRINSIC timing resolution — electronics not included.** At least four live
electronics variants and FWHM/sigma ambiguity remain; SPTR_PROVENANCE.md does
not establish sqrt(kN) propagation for order statistics. No TOP/END combination.

## Outputs, provenance, validation

`<output>/simulation.meta.json`, `run.mac`, `stdout.log`, `resource_usage.txt`,
`photon_hits_run000.root` preserve the simulation. `<output>/analysis/` contains:

- `analysis.root`: per-N histograms/functions, numeric results, END times and
  hit counts for ALL generated events, exact END fit histogram definitions.
- `analysis.csv`: all TOP curve rows, primary/biased results, both END fits.
- `analysis.meta.json`: full results and frozen models, uncertainties,
  diagnostics, source/analysis/orchestration commits, script/input/PDE hashes,
  commands, UTC timestamps, exit codes, runtime versions and RNG information.
- `gate.json`: literal G-P criterion, five PASS/FAIL decisions with numeric
  evidence, sidecar hashes. Diagnoses do not silently add gate conditions.

`run_cell.py` writes unique invocation journals with subprocess commands and
actual exit codes. The simulator measures process peak RSS with GNU time -v,
startup through the pre-beam initialization marker and total wall time.
Analysis records its own getrusage peak RSS for the later concurrency budget.

```bash
/usr/bin/python -m unittest discover \
  -s /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration -p 'test_*.py' -v
ctest --test-dir /home/rrios/exec33_20260911/build_off --output-on-failure
```

The EXEC34A handoff outside the repository records actual results and the exact
EXEC34B start/resume commands. G-P PASS prepares that handoff; **it does not
launch the grid**. The eventual grid uses the same even-TRAIN orientation.

## Preserved EXEC33 inventory

[INVENTORY.md](INVENTORY.md), [ESTIMATORS.md](ESTIMATORS.md) and manifests under
`provenance/` document the unchanged imported workflows, historical conflicts
and electronics variants. Their original pending-choice statements describe
EXEC33; the choices above supersede them specifically for EXEC34.

## EXEC34B prepared runner (not executed by EXEC34A)

The handoff pins the complete campaign in
`/home/rrios/exec34b_20260911/campaign.json` and prepares an append-only
`manifest.jsonl`. Cells: EJ-200/OPSC-100, EJ-204/OPSC-101, EJ-230/OPSC-106;
x=0,+200,-200,+500,-500,+650,-650 mm; all 10000 events, one worker,
eventModulo=1, same seeds, frozen even-TRAIN convention and intrinsic parameters.
Every new position learns its own TRAIN model and then freezes it for EVAL.

Read-only validation and plan:

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/grid.py \
  --config /home/rrios/exec34b_20260911/campaign.json
```

EXEC34B start (EXEC34A does not execute this command):

```bash
/usr/bin/python /home/rrios/ej200_exec33_20260911/analysis/sigma_t/orchestration/grid.py \
  --config /home/rrios/exec34b_20260911/campaign.json --execute
```

Resume uses the same command plus `--resume`. PASS cells are skipped only after
verifying their sidecars; completed simulation stages can be reused. Incomplete
or inconsistent outputs are preserved and refused, not overwritten or deleted.
A filesystem lock prevents two concurrent runners. A failed cell stops new
submissions; already active cells finish and are recorded. Each cell retains its
stdout and every stage exit code. No change of physics is attempted on failure.

Concurrency is capped by 24 logical CPUs, pending cells, the recorded 20% RAM
reserve and twice the larger measured simulation/analysis RSS. This bounds the
complete one-worker process chain conservatively; simulation-only memory allows
more processes. Available RAM is checked again at launch and may lower the cap.
OMP/OpenBLAS/MKL thread counts are one inside each independent cell. These are
execution limits, not changes to the estimator. No simultaneous-job throughput
is claimed to have been measured by this pilot.
