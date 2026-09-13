# EXEC37: V2 same-end timing on archived END-only D0/D3

No simulation can be launched by this script. Both archived ROOTs must resolve
and pass hashes before either timing calculation starts. The EXEC29/30 manifest
format is `manifest.meta.json` with `cells.<id>.cwd`; resolve that cwd plus
`photon_hits_run000.root`, checking the original ROOT sidecar hash and the
matching `manifest.csv` row. D0/D3 commits are fixed in `exec37_config.json`.

From this worktree:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/sigma_t/orchestration/run_exec37.py
```

This reuses the established EXEC35 Python environment, EXEC36 `width` and its
unchanged Gaussian FitCore/bootstrap, and the exact original LeadingEdgeTime
wrapper. The new adapter uses **2000 generated IDs** from each manifest, not
the hard-coded 10000-event cache adapter of EXEC36. Missing crossings remain
NaNs until the event bootstrap. Native hit counts recompute N_pe/end and SEM
and must match the same run's manifest. Baseline is V2, same end, left IDs
0..3 versus 4..7 and right 8..11 versus 12..15, rise/fall 0.5/5 ns, 4 PE.

Output: `/home/rrios/exec37_20260913/cells/{D0,D3}/results.{root,csv,meta.json}`.
ROOTs retain all 2000 IDs, T1/T2/difference/acceptance, per-event light yields,
histograms and 500 Gaussian/core bootstrap replicas (PCG64 seed 35091201).
The manifest checkpoints completed triples; verified completed cells may resume,
but incomplete existing directories are not overwritten. No parameter retry.

The preregistered C1 uses normalized **robust** point width <=88.4 ps, provided
EXEC35 robust validity passes; Gaussian validity remains a separate comparator.
C2 has no pass/fail and compares matching-end EXEC36 EndTop widths with each
END-only result. The quoted 1.719 benchmark derives from original EXEC29 D0/D1
light yields, not a new same-N timing control. Cross-run covariance is unknown;
the propagated ratio error is explicitly a zero-covariance approximation.
C3 is a conditional, indirect quadrature residual; a robust quantile width is
not a variance, and no test-beam uncertainty was supplied. Bootstrap replicas
outside the real square-root domain are counted and omitted, not clipped.

This is x=0 only, N=2000, with no END-only position scan. No electronics is
injected and no previous report or EXEC36 sanity flag is rewritten. Sensor/PDE,
reference-time and quantization discrepancies persist; “Mylar” remains legacy
nomenclature for the real Vikuiti 3M ESR reflector. D3's historical dielectric
factory is an intentional code control, not an assertion of different real wrap.

Tests:

```bash
/home/rrios/exec35_20260912/venv/bin/python -m unittest discover -s analysis/sigma_t/orchestration -p test_exec37.py -v
```
