# EXEC_07 Phase-1 validation report

Date: 2026-06-10

## Scope and branch checkpoint

- Working branch: `feat/endtop-sslg4`.
- Preserved branch: `fix/readout-wrapping-config`.
- Preserved branch local/remote HEAD:
  `2a9d57b454cfb2e659c8551b6c8d12c1b2a34b2b`.
- Pushed backup tag: `checkpoint/pre-endtop-sslg4-2026-06-10`.
- Phase-2 scan status: **not launched**.

## Clean build and automated validation

The Phase-1 implementation was configured and built from an empty
`build-exec07-clean` directory with `RelWithDebInfo`. The build completed
without new warnings.

CTest result: **7/7 passed** in 22.40 s.

| Test | Result | Purpose |
|---|---:|---|
| `smoke_test` | pass | Existing baseline smoke |
| `edge_scan_smoke` | pass | Fast edge configuration smoke |
| `sslg4_properties_check` | pass | Effective OPSC-101 MPT and `DumpTable()` |
| `readout_config_check` | pass | 86 unique IDs, exact Top positions, PDE and wrapping |
| `export_endtop_gdml` | pass | Export EndTop GDML |
| `endtop_gdml_check` | pass | Parse GDML and check all 70 Top positions |
| `endtop_balance_smoke` | pass | Reject center open-face-level leakage |

The exported GDML contains 16 End plus 70 Top SiPM placements. Top IDs are
contiguous from 16 through 85 and occur at the requested positions within
0.01 mm. The central pair is at -12/+12 mm, preserving the deliberate 24 mm
gap.

## SSLG4 OPSC-101 check

The actual bar MPT produced by `OrganicScintillatorFactory` passed:

| Property | Effective value |
|---|---:|
| `RINDEX` | 1.58 |
| `SCINTILLATIONYIELD` | 10,400 photons/MeV |
| `SCINTILLATIONTIMECONSTANT1` | 1.8 ns |
| `SCINTILLATIONRISETIME1` | 0.7 ns |
| `ABSLENGTH` | 160 cm |
| emission maximum | 408.8 nm |

No discrepancy was found between vendored SSLG4 OPSC-101 and the requested
EJ-204 datasheet values. SSLG4 OPSC-100, not OPSC-101, contains the known
380 cm EJ-200 attenuation length.

## Manual 50-event EndTop smokes

| Beam x | IDs with hits | Dominant local Top IDs | Generated photons | Detected photons | Bar-to-world escape |
|---:|---:|---|---:|---:|---:|
| 0 mm | 86/86 | 50, 51 | 977,877 | 239,903 | 29,983 (3.07%) |
| -650 mm | 86/86 | 19, 18 | 1,021,774 | 284,947 | 148,763 (14.56%) |

At the center, the nearest Top pair dominates and wrapped-face escape remains
far below the historic open-face failure. Near -650 mm, the nearest Top
elements dominate and the expected fully open left End face raises the escape
fraction. The automated center balance guardrail passed.

The Phase-1 analysis script was smoke-tested only against the 50-event center
file and produced all expected output families plus 96 summary rows. These
smoke plots are not campaign deliverables.

## Geometry and channel-count change map

| File:line | Change |
|---|---|
| `include/DetectorConstruction.hh:37-39` | Authoritative 16 End, 70 Top, 86 total counts |
| `src/DetectorConstruction.cc:45-48` | Fixed mirrored Top position map |
| `src/DetectorConstruction.cc:185-215` | SSLG4 material factory and EJ-204 runtime guards |
| `src/DetectorConstruction.cc:289-308` | Physical +Y reflector with 70 exact windows |
| `src/DetectorConstruction.cc:322-345` | Preserved EXEC_06 End placement |
| `src/DetectorConstruction.cc:347-364` | Top placement and individual border surfaces |
| `tests/readout_config_check.cc:71-97` | Unique IDs and exact Top position checks |
| `tests/readout_config_check.cc:138-152` | EndTop total and central-gap guardrails |
| `analysis/analyze.py:19-25` | Python channel groups updated to 86 |
| `analysis/analyze_hits.C:51` | ROOT channel count updated to 86 |
| `analysis/grouped_resolution.py:61-78` | Python grouped channel map updated |
| `analysis/grouped_resolution.C:38-56` | ROOT grouped channel map updated |
| `analysis/congruent_sum4_timing.C:216-225` | Top midpoint channels updated |
| `analysis/resolution_vs_x.py:30` | Channel count updated |
| `analysis/resolution_vs_x_fixed.py:34` | Channel count updated |
| `analysis/high_stats_position_scan.C:243` | Top midpoint channel updated |
| `analysis/preliminary_position_scan.C:6,183` | Total count and midpoint updated |
| `macros/event_display_top_*.mac`, `macros/scan_angle.mac` | Top midpoint channels updated |

The final semantic hardcode search over active source, tests, macros, scripts,
analysis, CMake and README returned no residual legacy 36-channel map or
20-Top-SiPM constants:

```bash
rg -n '16-35|16–35|16\.\.35|16…35|16, 36|range\(16, 36\)|N_TOP_SIPMS=20|fNTopSiPMs=20|NSIPM_TOP|36 channels|36-channel' \
  src include tests macros scripts analysis CMakeLists.txt README.md
```

## Prepared but deliberately not run

- `scripts/run_exec07_scan.sh`: resumable 31-position campaign, 2,000
  events/position and 16 workers.
- `analysis/exec07_photon_budget.py`: N_pe/Poisson, arrival-time/FPT, position
  trends, End DeltaT and Top timing-estimator outputs.

The requested campaign plots and `summary_exec07.csv` require Phase 2 and are
therefore intentionally absent at this gate.
