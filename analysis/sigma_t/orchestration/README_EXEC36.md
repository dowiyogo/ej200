# EXEC36 same-end test-beam observable

Existing-ROOT re-analysis only. `exec36_config.json` preregisters physical IDs,
primary end, sanity thresholds and parameter grid before calculation. Rollback
tag `pre-exec36-20260913` precedes edits. The external preregistration copies
are `/home/rrios/exec36_20260913/preregistration.{md,json}`.

From `/home/rrios/ej200_exec33_20260911`, execute the three analysis phases in
this order (none can launch Geant4):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/sigma_t/orchestration/run_exec36.py --phase single --jobs 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/sigma_t/orchestration/run_exec36.py --phase groups --jobs 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/sigma_t/orchestration/run_exec36.py --phase sensitivity --jobs 1
```

The environment is the same isolated NumPy 1.23.5 repair documented for EXEC35.
Each phase checkpoints hashed ROOT/CSV/meta triples and resumes COMPLETE cells
only after verifying hashes and source identity. Existing incomplete output
directories are not overwritten. There is no timeout or parameter retry.

`single` processes exactly 21 latest SIMULATION_COMPLETE manifest records and
fixed physical IDs 0/8 relative to gun t=0. It preserves the exact native END
hit columns in its ROOT, so later phases need not reread TOP hits. It keeps all
10000 generated event IDs, including missing crossings. Both single-SiPM
distributions finish before **any** grouped distribution is calculated: groups
and sensitivity require all 21 single phase outputs. The sanity decision uses
the canonical core width, with Gaussian FitCore comparative. Below 100 or above
500 ps is extreme failure; 150–350 ps inclusive passes; the gaps are explicitly
indeterminate. Any baseline extreme failure labels all group results NOT
COMPARABLE until discriminator calibration. Gaussian and robust validity use
the unmodified EXEC35 rule independently.

`groups` computes CH1 minus CH5 on each end, accepting only two finite crossings.
V2 maps left [0,1,2,3] versus [4,5,6,7]; V1 maps [0,2] versus [4,6]. Right adds 8.
These are low/high spatial clusters: DetectorConstruction.cc:407–419 assigns
local index i to y=(i−3.5)*7.5 mm and places left/right at opposite X. V2 matches
the preserved GroupIndex. V1 retains one physical member per declared input
pair; its unused IDs are recorded. Mapping to hardware input 0+2 -> CH1 and
4+6 -> CH5 is a convention awaiting experimental confirmation, not a discovered
wiring diagram. Neither variant is selected as canonical.

Pulse, LeadingEdgeTime, FitCore and GroupIndex are included directly from the
preserved source; no macro main executes. A decisive control reconstructs the
old opposite-end min-of-two-clusters reduction from the new V2 event times and
requires **exact equality of all 10000 DeltaT_LR values** with EXEC35.

`sensitivity` uses EJ204_xp0 and all six declared rise/fall/threshold choices.
It compiles a verbatim source slice containing SprPeakTime, SprNorm, Pulse and
LeadingEdgeTime into a separate namespace with only three constants exposed.
The baseline still calls the original functions. Tests require bit-exact
agreement at baseline and exact cache/union behavior with missing events. The
2/3 ns case is a parameter variation of the **same normalized difference of
exponentials**, not substitution of the product kernel in sipm_waveform_dcfd.py.
It introduces no SPTR, jitter, TDC quantization or canonical calibration choice.

Each row contains raw Gaussian/core widths and their /sqrt(2) versions, paired
event-bootstrap ratio uncertainty, all EXEC35 robust comparators, accepted
counts, IDs, variant, pulse and threshold. Single-SiPM /sqrt(2) helper fields
are explicitly nonphysical and must not be reported as the single-SiPM
observable. The same-end /sqrt(2) assumes equal independent group contributions
as in the test-beam convention; this is not verified by the division. 500
replicas, PCG64 seed 35091201, resample generated event IDs with fixed grouping.
All 21 cells share simulation seeds 26092601 and 8349041; no independent-cell
significance, sigma(x) fit, BLUE or arm combination is performed.

Test command:

```bash
/home/rrios/exec35_20260912/venv/bin/python -m unittest discover -s analysis/sigma_t/orchestration -p test_exec36.py -v
```
