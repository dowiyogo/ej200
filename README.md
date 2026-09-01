# ej200_v2 - EXEC_07 EndTop optical simulation

Geant4 simulation of the SHiP timing-detector scintillator bar. The EXEC_07
configuration combines End and Top readout, uses SSLG4/OPSim for the active
scintillator, and uses the Broadcom AFBR-S4N66P024M PDE curve.

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
`scripts/run_exec07_scan.sh`. It writes `photon_hits_x{pos}mm.root`, validates
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
