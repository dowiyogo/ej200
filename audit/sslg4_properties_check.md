# EXEC_07 SSLG4 OPSC-101 effective-property check

Date: 2026-06-10

## Provenance

- SSLG4 upstream: `mkandemirr/SSLG4` commit
  `9b1046fa99e05f7969255be0ad539f5db54d62f5`.
- OPSim upstream: `mkandemirr/OPSim` commit
  `a35804a21eac3be8e74d04781a01080e6187a36d`.
- Vendored runtime macro:
  `src/external/SSLG4/macros/oscnt/opsc-101.mac`.
- Vendored wavelength-dependent data:
  `src/external/SSLG4/data/oscnt/opsc-101/`.

The bar material is obtained from
`OrganicScintillatorFactory::GetInstance()->Get(fScintCode, true)` in
`src/DetectorConstruction.cc:185-186`. `Materials::CreateEJ200()` and
`Materials::CreateEJ204()` remain in the tree for historical/reference use but
are not referenced by `DetectorConstruction`.

## Effective OPSC-101 values

| Property | SSLG4 source | Effective runtime value | EJ-204 target | Result |
|---|---|---:|---:|---|
| `RINDEX` | `rIndex.txt:1-2` | 1.58 | 1.58 | pass |
| `SCINTILLATIONYIELD` | `opsc-101.mac:9` | 10,400 photons/MeV | 10,400 photons/MeV | pass |
| `SCINTILLATIONTIMECONSTANT1` | `opsc-101.mac:11` | 1.8 ns | 1.8 ns | pass |
| `SCINTILLATIONRISETIME1` | `opsc-101.mac:12` | 0.7 ns | 0.7 ns | pass |
| `ABSLENGTH` | `absLength.txt:1-2` | 160 cm | 160 cm | pass |
| `SCINTILLATIONCOMPONENT1` peak | `scntComp1.txt` | 408.8 nm | approximately 408 nm | pass |

`tests/physics_baseline_check.cc:42-89` constructs the actual detector material,
checks every value above, confirms the bar material name is `opsc-101`, confirms
finite rise time is enabled, and calls `DumpTable()`.

## Runtime guards

`src/DetectorConstruction.cc:187-215` checks the effective OPSC-101 MPT after
SSLG4 construction. If a future vendored SSLG4 update changes or omits
`ABSLENGTH`, it is replaced with 160 cm. If it changes or omits
`SCINTILLATIONRISETIME1`, it is replaced with 0.7 ns. The current vendored
files already contain both correct values, so no override message is emitted.

## Discrepancy assessment

No discrepancy was found between the vendored SSLG4 OPSC-101 data and the
requested EJ-204 datasheet values. The known 380 cm trap is present only in
SSLG4 OPSC-100 (EJ-200), not OPSC-101.
