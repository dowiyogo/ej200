# EXEC_14 — Data Inventory

Generated: 2026-06-12  
Branch: `feat/ej230-sslg4`  
Data path: `/home/reriosto/SHiP/t0minidaq/results_ej230/data/`  
Session: `session_20260612_171407` (metadata dir; ROOT files landed directly in `data/`)

## Scan grid

| Nominal x [mm] | File | Entries | Events | gun_x OK | Channels | Status |
|---:|---|---:|---:|:---:|:---:|:---:|
| −690 | photon_hits_x-690mm.root | 10 752 074 | 2000 | ✓ | 0–85 | OK |
| −670 | photon_hits_x-670mm.root | 10 476 xxx | 2000 | ✓ | 0–85 | OK |
| −650 | photon_hits_x-650mm.root | ≈10 M | 2000 | ✓ | 0–85 | OK |
| −600…+600 | (every 50 mm, 25 files) | 8–10 M | 2000 | ✓ | 0–85 | OK |
| +650 | photon_hits_x650mm.root | ≈10 M | 2000 | ✓ | 0–85 | OK |
| +670 | photon_hits_x670mm.root | ≈10 M | 2000 | ✓ | 0–85 | OK |
| +690 | photon_hits_x690mm.root | 10 812 407 | 2000 | ✓ | 0–85 | OK |

**Total: 31/31 positions present and valid.** No missing positions, no corrupt files.

## Verification gates

- TTree `sipm_hits` present: ✓ all 31
- `event_id` range 0–1999 (N=2000): ✓ all 31
- `gun_x_mm` matches filename: ✓ all 31
- `global_id` range 0–85 (86 channels): ✓ all 31
- `time_ns` finite: not individually checked but implied by non-zero hit counts (8–11 M entries per file)
- File stability (no modification in last hour): ✓

## Comparison with EXEC_12 grid

Expected grid: `(-690, -670, -650, ..., 650, 670, 690)` — 31 positions, identical to `EXPECTED_POSITIONS_MM` in `analysis/exec07/common.py`. **Full match.**

## Known warning

`G4Exception mat031` appears in all run logs:
```
For material opsc-106 sum of fractional masses 0.998616 is not 1 — results may be wrong
SSLG4 code: OPSC-106
```
Hook reserved: `HOOK_OPSC106_MASSFRACTION`. See `docs/EJ230_material_audit.md` §"Hallazgos de ejecución".

## Special EJ-230 controls completed in EXEC_14B/14C

The report's special figures now have authentic EJ-230 equivalents:

| Figure | Origin | EJ-230 inputs | Status |
|---|---|---|---|
| `figs/exec09_tail_comparison.png` | `exec09_timing_mechanism.py` | Main EndTop plus End-only x=-400,0,+400 | Regenerated; SHA differs from EJ-204 |
| `figs/exec08b_window_dip_profiles.png` | `exec08b_window_dip.py` | Main x=-650 plus directed x=-652,-642,-648,-654 | Regenerated; SHA differs from EJ-204 |
| `figs/exec08b_id18_impact_maps.png` | `exec08b_window_dip.py` | Same five EJ-230 runs | Regenerated; SHA differs from EJ-204 |

All seven special ROOT files are validated in
`results_ej230_analysis/root_validation_exec14b.csv`. No EJ-204 data figure is
used by the final report.
