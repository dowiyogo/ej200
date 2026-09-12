# EXEC35: existing-ROOT robust-width analysis

Primary and validity thresholds were published before application in
`../validity_exec35.json` at 2026-09-12T21:46:12+00:00 (commit b261096).
The immutable external copies are `/home/rrios/exec35_20260912/preregistration.*`.

Run from this worktree, using the isolated same-version NumPy repair:

```bash
/usr/bin/python -m venv --system-site-packages /home/rrios/exec35_20260912/venv
/home/rrios/exec35_20260912/venv/bin/python -m pip install --ignore-installed --no-deps --only-binary=:all: numpy==1.23.5
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/sigma_t/orchestration/run_exec35.py --jobs 4
```

The global NumPy 1.23.5 installation lacks `numpy._typing`, preventing SciPy
imports. The venv supplies the complete 1.23.5 wheel without changing global
packages or analysis primitives. The first two commands are environment setup;
the last command is the reproducible analysis. `--jobs` is analysis process
concurrency, not Geant4 workers. This script cannot launch Geant4.

Inputs are the latest valid COMPLETE records in the EXEC34C analysis manifest
and SIMULATION_COMPLETE records in the EXEC34R manifest. All 21 matching cells
are required. The runner hashes each raw ROOT and PDE, checks the inherited
analysis sidecars and frozen canonical model digest, and verifies effective
counts and frozen fit-window counts against EXEC34C.

TOP has 60 distributions per cell: TRAIN even and EVAL odd at each fixed N=1..20
using the canonical frozen channels, plus full sample N=1..20 using the frozen
full-sample channels. Procedure (a) is each fixed-N curve, not a newly chosen
global N. Procedures (b)/(c) alias their existing Gaussian-selected winning
indices; there is no new robust argmin. END reuses the exact cached generated-event
DeltaT_LR array from the hashed prior analysis ROOT (same SUM4 leading edge).
END widths refer directly to DeltaT_LR. No arms are combined.

The EXEC12T `robust` function body is loaded through AST, without importing its
module. Its original 16/84 request is explicitly adapted to 15.865/84.135;
the imported source is unchanged. The unchanged `phase_ab.rsigma` function is
also AST-loaded; its module's campaign/plotting statements never execute.

Bootstrap uses 500 generated-event resamples, PCG64 seed 35091201, retaining
NaNs in the resampled event population and omitting them only when measuring.
All quantile and FWHM widths/errors use these resamples. ROOT files preserve
the event timestamps, histograms and replicate vectors. The CSV and JSON retain
raw quantiles, bootstrap SE and percentile intervals. The FWHM convention and
validity rules are fully specified in the preregistration. Gaussian TOP points
and bootstrap errors are inherited (300 replicas, seed 20260618; successful
fraction unarchived). END points are inherited and their unchanged FitCore and
FitPeakSeeded primitives receive the new event bootstrap. Neither imported
macro's main is called. No electronics or jitter is added.

The ratio error holds the previous Gaussian denominator fixed; it is not a
joint uncertainty on the ratio. Width errors are conditional on the frozen
selection, not uncertainty from relearning it. All 21 cells share simulation
seeds 26092601 and 8349041. A sign-test p-value is nominal under independent
signs only; no independent-cell significance or sigma(x) fit is justified here.

Outputs are new `/home/rrios/exec35_20260912/cells/<cell>/widths.{root,csv,meta.json}`.
Completion is journaled after hashing. A rerun skips verified COMPLETE outputs;
an incomplete existing directory causes a visible failure rather than overwrite.
There is no analysis timeout or automatic retry changing any parameter.

Tests:

```bash
/home/rrios/exec35_20260912/venv/bin/python -m unittest discover -s analysis/sigma_t/orchestration -p 'test_*exec35.py' -v
/home/rrios/exec35_20260912/venv/bin/python -m unittest discover -s analysis/sigma_t/orchestration -p test_robust_widths.py -v
```
