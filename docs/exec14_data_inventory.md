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

## Figures with no EJ-230 equivalent

Two figures in the Beamer template require special simulation campaigns not executed for EJ-230:

| Figure | Origin | Reason not regenerated for EJ-230 |
|---|---|---|
| `figs/exec09_tail_comparison.png` | `exec09_timing_mechanism.py` | Requires End-only runs (no Top SiPMs) — not in EXEC_13 scope |
| `figs/exec08b_window_dip_profiles.png` | `exec08b_window_dip.py` | Requires window-dip scan runs — not in EXEC_13 scope |
| `figs/exec08b_id18_impact_maps.png` | `exec08b_window_dip.py` | Same as above |

These three figures are purely geometric effects (same geometry for EJ-204 and EJ-230). The EJ-204 versions are used verbatim with explicit provenance notes in the Beamer captions. Documented in `docs/exec14_fit_failures.md`.
