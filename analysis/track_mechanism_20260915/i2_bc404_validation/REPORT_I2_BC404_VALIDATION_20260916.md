# EXEC_46 I2 — BC-404 validation at EJ204_xm650

Date: 2026-09-16. Existing ROOT analysis plus one authorized 10,000-event validation cell.

The log contains the inherited Geant4 `mat031` warning that the OPSC-101 fractional masses do not sum exactly to one. It is also present with the baseline material data and is not introduced by the RINDEX/ABSLENGTH replacement.

| scenario | first Cherenkov [%] | primary caustic width [deg] | first-Cher local velocity [mm/ns] | first-scint local velocity [mm/ns] | MPT v_group(408 nm) [mm/ns] |
|---|---:|---:|---:|---:|---:|
| constant_n_baseline | 62.90 +/- 0.48 | 2.435 +/- 0.143 | 150.017 +/- 0.101 | 181.255 +/- 0.111 | 189.742 |
| bc404_3800 | 72.88 +/- 0.44 | 6.026 +/- 0.205 | 149.356 +/- 0.093 | 173.078 +/- 0.096 | 171.893 |

The corrected-minus-baseline first-photon Cherenkov difference is 9.98 +/- 0.66 percentage points. The preregistered direction gate is **PASS**: 72.88% is above the measured constant-index baseline of 62.90% (the registered rounded value is 62.90%).

The production campaign may proceed only when this gate is PASS and the four MPT hashes are independently verified during campaign preparation.

Reproduce with:

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_i2_bc404_validation.py
```
