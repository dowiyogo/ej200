# ej200_v2 - SHiP Timing Detector (T0) Optical Simulation & Analysis

Geant4 11.x + SSLG4/OPSim simulation and ROOT/PyROOT timing analysis of the SHiP timing-detector plastic scintillator bar (`EJ-200` / `OPSC-100`, `EJ-204` / `OPSC-101`, `EJ-230` / `OPSC-106`) read out by Broadcom `AFBR-S4N66P024M` SiPMs. The baseline `EXEC_07` configuration combines End and Top readout (86 SiPM channels) and produces `.root` files containing the `sipm_hits` `TTree`.

## Repository Structure (Reorganized 2026-09-25)

- **Critical Coupling Catalog (ROOT `.C` macros & Python/CMake modules):** [`docs/catalogo_macros.md`](docs/catalogo_macros.md)
- **Reorganization Execution Log & SHA-256 Verification Ledger:** [`docs/execution_logs/REORG_20260925.md`](docs/execution_logs/REORG_20260925.md)

| Directory | Contents & Active Entry-Points |
|---|---|
| `src/`, `include/`, `main.cc`, `CMakeLists.txt` | Core Geant4 11.x C++ simulation (`ej200_bar_sim`), detector geometry (`DetectorConstruction`), sensitive detector (`SiPMSD`), and vendored `src/external/{OPSimTool,SSLG4}`. |
| `macros/`, `sslg4/` | Geant4 runtime macros (`run.mac`, `scan_*.mac`, `endtop_smoke_*.mac`) and SSLG4 optical property tables (`OPSC-100`, `OPSC-101`, `OPSC-106`). **Protected by `CMakeLists.txt` `file(COPY ...)`** — do not move without updating `CMakeLists.txt` and CTest. |
| `sim/scripts/` | Active standalone simulation scan runners: `sim/scripts/run_exec07_scan.sh`, `sim/scripts/run_end_tir_scan.sh`, `sim/scripts/run_end_vikuiti_scan.sh`, `sim/scripts/post_end_vikuiti.sh`. |
| `scripts/` | Build/CTest/CI helpers (`scripts/build.sh`, `scripts/verify_reproducibility.sh`), CMake-coupled runner (`scripts/run_scan.sh`), and subdirectories:<br>• `scripts/legacy/`: Archived non-active scripts (`resume_scan_2.sh`).<br>• `scripts/_quarantine_do_not_run/`: Quarantined non-executable (`chmod -x`) one-off/scratch scripts that would overwrite reports or CSVs if executed. |
| `analysis/` | Active ROOT/PyROOT and `uproot` analysis pipelines:<br>• `analysis/track_mechanism_20260915/`: **EXEC_46** mechanism & dispersive-optics pipeline (`prepare_campaign.py`, `analyze_step1.py`, `build_step2_derived.py`, `analyze_step2.py`, `build_step3_transport.py`, `analyze_step3.py`, `build_step4_pairs.py`, `analyze_step4.py`, `analyze_step5.py`, `analyze_step5_revision.py`, `analyze_step6_widths.py`, `analyze_veff_rank_cfd.py`, `aggregate_cfd_report.py`, and LaTeX report in `report/main.tex`). All step scripts enforce explicit `--campaign-dir` and `--output-dir` CLI flags.<br>• `analysis/tsum_veff_20260914/`: Effective velocity & $t_{\text{sum}}$ study (`rebuild.sh`, `macros/*.C`).<br>• `analysis/timing_symmetry_20260914/` & `analysis/order_stat_weight_20260915/`: Timing symmetry and order-statistic weighting macros (`macros/*.C`).<br>• `analysis/sigma_t/orchestration/`: Campaign grid orchestrator (`detached_grid.py`, `prepare_grid.py`).<br>• `analysis/validation/`: Statistical validation suite (`run_exec32_suite.py`, `exec31_stats.py`).<br>• `analysis/timing/`: Waveform + dCFD (`sipm_waveform_dcfd.py`, `sipm_waveform_dcfd.cpp`, `pulse_models.py`) and arrival-time estimators.<br>• `analysis/optim/`, `analysis/exec07/`, `analysis/exec13/`, `analysis/exec14/`: Historical campaign pipelines and estimators (`resolution_vs_x_fixed.py`, `edge_resolution.py`, `grouped_resolution.py`, `ResolutionScan_v2.C`). |
| `presentations/` | Self-contained compilable LaTeX Beamer decks:<br>• `presentations/v9p1/` (`talk_v9p1.tex`, `rebuild_v9p1.sh`)<br>• `presentations/v9/` (`talk_v9.tex`, `rebuild_v9.sh`)<br>• `presentations/v8/` (`talk_v8.tex`, `rebuild_v8.sh`)<br>• `presentations/v7/` (`talk_v7.tex`, `slides/*.tex`)<br>• `presentations/v6/`, `v5/`, `v4/`, `napkin_first_principles/`, `exec14/` |
| `data/` | Tracked SiPM PDE curves (`data/sipm/AFBR-S4N66P024M_pde.txt`) and untracked consolidated campaign storage (`data/ej200_campaigns/{raw,derived,quarantine_corrupt_or_duplicate}`, ignored by Git and excluded from CMake `file(COPY)`):<br>• `raw/OPSC-106_EJ230_endtop_exec13_from_ej230/`<br>• `raw/OPSC-101_EJ204_endonly_mylar_20260614_from_ej200_end/`<br>• `raw/OPSC-106_EJ230_endonly_mylar_20260614_from_ej230_end/`<br>• `raw/OPSC-101_EJ204_endtop_scans_202606_08_from_ej204/`<br>• `raw/EJ228_cylinder_tir_vs_vikuiti_20260815_from_ej204/` |
| `docs/` | Project documentation:<br>• `docs/catalogo_macros.md`: Complete coupling catalog of 115 ROOT macros and Python/CMake modules.<br>• `docs/execution_logs/REORG_20260925.md`: Master reorganization log, SHA-256 verification ledger, and manual cleanup commands.<br>• `docs/reports/`: Consolidated technical audit reports (`GROUP_VELOCITY_AUDIT.md`, `REPORT_branch_status_*.md`).<br>• `docs/branch_diagnosis/` & `docs/literature/`: Branch consolidation records and literature notes. |
| `tests/` | CTest C++ and Python regression tests (`check_endtop_balance.py`, `readout_config_check.cc`, `sslg4_properties_check.cc`). |

## Geometry

| Volume/readout | Dimensions/count | Material/model |
|---|---:|---|
| Scintillator bar | 1400 x 60 x 10 mm | SSLG4 `OPSC-101` (EJ-204) |
| End SiPM elements | 16, IDs 0-15 | 6 x 6 mm2, unchanged EXEC_06 placement |
| Top SiPM elements | 70, IDs 16-85 | 6 x 6 mm2, Broadcom AFBR-S4N66P024M |
| Coupling volumes | one per SiPM | SiO2 approximation, n=1.58 |
| Reflector panels | all faces except End faces and Top windows | Mylar-like `dielectric_metal`, R=0.98 |

The Top centers are fixed and symmetric:

```text
IDs 16-50: x = -692, -672, ..., -12 mm
IDs 51-85: x =  +12,  +32, ..., +692 mm
```

All regular Top pitches are 20 mm. The central pair at -12/+12 mm has a
deliberate 24 mm pitch so the first elements remain 5 mm from both bar ends.
The `+Y` physical reflector panel is a chained `G4SubtractionSolid` with 70
exact 6 x 6 mm2 windows at those positions. Each Bar-to-coupling border surface
therefore has priority at the active interface, while the rest of `+Y` remains
wrapped.

The hardware's exterior black Tedlar light-tight layer is deliberately not
modeled. It does not contribute internal reflection; absorption of the
non-reflected fraction `(1-R)` is already implicit in the Mylar optical surface,
and its material budget is negligible for the muon.

## Runtime selectors

Set selectors before `/run/initialize`:

```text
/det/readout EndTop
/det/scintillator OPSC-101
/sipm/model AFBR-S4N66P024M
```

`OPSC-101` is EJ-204 and is the default. `OPSC-100` selects EJ-200. The
historical `/det/scintillatorMaterial EJ204|EJ200` command remains as an alias.

For OPSC-101, the effective MPT is guarded by CTest:

| Property | Effective value |
|---|---:|
| `RINDEX` | 1.58 |
| `SCINTILLATIONYIELD` | 10,400 photons/MeV |
| `SCINTILLATIONRISETIME1` | 0.7 ns |
| `SCINTILLATIONTIMECONSTANT1` | 1.8 ns |
| `ABSLENGTH` | 160 cm |
| emission maximum | 408.8 nm |

The finite-rise-time sampler is enabled globally. The vendored SSLG4 OPSC-101
data already contains 160 cm and 0.7 ns; runtime guards enforce those values if
a future data update drifts.

## SSLG4 and OPSim

Vendored code and data live in:

```text
src/external/OPSimTool/
src/external/SSLG4/
```

Both upstream projects are GPL-3. Their upstream license/copying files and
README files are included in each vendored directory.

References:

- M. Kandemir et al., "OPSim: A simulation toolkit for optical photon
  processes in Geant4", Computer Physics Communications 292 (2023) 108873.
- M. Kandemir et al., "SSLG4: A novel scintillator simulation library for
  Geant4", Computer Physics Communications 306 (2025) 109385.

SSLG4 resolves `sslg4/macros/...` and `sslg4/data/...` relative to the process
working directory. CMake copies that complete runtime tree into the build
directory, so execute the simulator from the build directory. All supplied
CTest and batch scripts do this. The SiPM data loader uses the source `data/`
directory by default and can be overridden with `EJ200_DATA_DIR`.

## Broadcom PDE

`data/sipm/AFBR-S4N66P024M_pde.txt` is digitized from Broadcom
AFBR-S4N66P024M-DS105, Figure 6, at the datasheet's typical 12 V overvoltage.
It spans 250-900 nm and peaks at 63% at 420 nm. The same file is loaded into the
SiPM optical surface as `EFFICIENCY(lambda)` and read back by the SD for ntuple
metadata. `REFLECTIVITY(lambda)=0`.

Optical crosstalk (~23%) and afterpulsing (<1%) are intentionally not modeled
and are not folded into PDE. The product uses a clear epoxy mold compound; the
datasheet does not specify its refractive index, so the current SiO2 n=1.58
coupling remains an explicit approximation.

The physical AFBR-S4N66P024M package is a 2 x 1 array with 7 mm element pitch.
Each simulated placement represents one 6 x 6 mm2 element, so the package pitch
does not alter the EndTop placement defined above.

## Build and validation

```bash
cmake -S . -B build-exec07 -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-exec07 -j8
ctest --test-dir build-exec07 --output-on-failure
```

Important guardrails:

- `sslg4_properties_check`: dumps and validates the effective OPSC-101 MPT.
- `readout_config_check`: validates 86 unique global IDs, exact Top positions,
  Broadcom PDE, zero reflectivity, and the perforated `+Y` reflector.
- `export_endtop_gdml` + `endtop_gdml_check`: exports and parses EndTop GDML.
- `endtop_balance_smoke`: 50 events at x=0 and rejects open-face-level leakage.

Manual Phase-1 smoke macros:

```bash
(cd build-exec07 && ./ej200_bar_sim -m macros/endtop_smoke_center.mac)
(cd build-exec07 && ./ej200_bar_sim -m macros/endtop_smoke_edge.mac)
```

## Phase 2 artifacts

The re-entrant 31-position, 2000-event-per-point campaign is prepared in
`sim/scripts/run_exec07_scan.sh`. It writes `photon_hits_x{pos}mm.root`, validates
completed positions with uproot, and resumes without accepting interrupted ROOT
files. Do not launch it until Phase 2 is explicitly approved.

After the campaign:

```bash
python analysis/exec07_photon_budget.py results/exec07_endtop_2000 \
  --output-dir results/exec07_analysis
```

The analysis produces per-channel and grouped N_pe/Poisson and arrival-time/FPT
PDFs, central N_pe and time trends versus beam position, End-cluster DeltaT,
and `summary_exec07.csv`.
