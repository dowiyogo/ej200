# EXEC38 physical transport reference

This directory contains the preregistered reanalysis of existing EXEC34R and
END-only D0/D3 data. **It never launches a transport simulation.**

- `EXEC38_PREREGISTRATION.md`: hypotheses, scopes and tolerances committed before
  new calculations (`4e08426`). Existing published results were already known;
  V6 is explicitly retrospective.
- `run_exec38.py`: verifies the archived input/cache chain and calculates physical
  END photon profiles and same-end covariance from unchanged recorded timestamps.
- `build_exec38_report.py`: publishes the English report, scientific figures with
  ROOT/CSV/metadata sidecars, and the machine-readable golden contract.
- `GOLDEN_REFERENCE_20260913.json`: one V1–V7 entry each, exact configuration
  references, predictions, measurements/errors, tolerances and missing streams.
  **`ready_for_acceptance=false`**. Profile hypotheses fail and source/first-hit/
  incident-photon observables are missing; this is not a certified baseline.

## Reproduce the numerical analysis

Use the existing environment with NumPy 1.23.5, SciPy 1.13.1 and uproot 5.6.9:

```bash
cd /home/rrios/ej200_exec33_20260911
OPENBLAS_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/run_exec38.py --self-test
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python -u analysis/validation/run_exec38.py --out /home/rrios/exec38_20260913_reproduction
```

The output directory must not exist. The actual original analysis used
`/home/rrios/exec38_20260913`; all source manifests are read-only. The script
checks all consumed cache/metadata hashes, the 21 native ROOT schemas/entries/
sizes and macro hashes, and rehashes the two END-only native ROOTs. Full grid
ROOT hashes are inherited through the exact END-cache provenance already
verified in EXEC36, not falsely claimed to have been recomputed here.

The archived publisher reads `/home/rrios/exec38_20260913` and writes the report
outside the repository, `/home/rrios/REPORT_EXEC38_20260913.md`:

```bash
OPENBLAS_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/build_exec38_report.py
OPENBLAS_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/verify_exec38.py
```

That publisher regenerates its own report/figures/contract; it does not change
source simulation files or recalculate the numerical analysis. Archived V6
results are transcribed with hashes, not recomputed. The additional legacy V5
table is explicitly a schema audit and never an internal-counter acceptance gate.
The verifier checks all contract artifact hashes, all V7 bootstrap errors and
variance identities, the V3/V4 fit chi-square values reconstructed from saved
full covariance matrices, and readable numeric ROOT sidecars for every figure.

## Statistical definitions and future engine interface

V3 counts all generated events, including zeros, and fits seven distances per
end (`700+x`, `700-x`). V4 uses photon-weighted arrival means, not discriminator
times. Their full covariance and fit bootstrap samples are in the V3 ROOT
sidecar; V4 central/error rows also have their own trio. Rejected model parameter
errors are conditional descriptive errors, not physical model uncertainties.

V7 reuses exact EXEC36/37 timestamps. V2 groups are left `{0,1,2,3}` /
`{4,5,6,7}` and right `{8,9,10,11}` / `{12,13,14,15}`; fixed rise 0.5 ns,
fall 5 ns, threshold 4 PE. Secondary V1 uses left `{0,2}` / `{4,6}`, right
plus eight. Both groups must cross. All generated event indices, including
missing crossings, enter the bootstrap; standard deviations condition on paired
finite timestamps. Ordinary ddof=1 standard deviations are required by Pearson's
variance identity; the EXEC35 robust core estimator is distinct and unchanged.
No new pulse, noise, SPTR, TDC, BLUE or arm combination is introduced.

Every new statistical bootstrap uses 500 replicas, PCG64 seed 38091301. Identical
generated-event resamples preserve pairing within equal-N ensembles; different
N populations are separate. Tolerances were frozen before data analysis.

The JSON contract defines physical counts, destinations and arrival timestamps,
not Geant4 callback/status names. Missing V1/V2/V5 streams are explicit required
fields. An alternative engine must resolve applicability and provide those
physical records; it must not satisfy the contract merely by reproducing a
historical Geant4 counter or one of its rejected fitted profiles. CPU exact
worker reproducibility is a separate same-RNG fixture, not a requirement that
a GPU reproduce CPU random-number sequences. No tolerance is retuned here.

Actual simulation seeds, N, worker configuration, material, source commits,
commands and ROOT hashes live in each configuration entry and input metadata.
Published report links point to the complete artifacts. Rollback tag
`pre-exec38-20260913` preserves the initial worktree at `35d6bcc`.
