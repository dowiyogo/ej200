# Physics Audit Per Branch

Scope: `main`, `feature/edge-scan-and-readout-grouping`,
`feature/sipm-electronics-response`.

Canonical sources cited here:

- `EJ200.pdf` and `EJ200_EJ204_EJ208_EJ212.pdf` (Eljen Technology, Aug 2023):
  EJ-200 rise 0.9 ns, decay 2.1 ns, attenuation 380 cm, peak 425 nm,
  light yield 10000 ph/MeV, RINDEX 1.58.
- `EJ228_EJ230.pdf` (Eljen Technology, Aug 2023): EJ-230 rise 0.5 ns,
  decay 1.5 ns, attenuation 120 cm, peak 391 nm, yield 9700 ph/MeV.
- `A_prototype_for_the_SHiP_timing_detector.pdf` (Betancourt 2020, NIM A 979)
  and `W4_Weekly_meeting_Timinig_Detector_Status.pdf`: prototype SiPM is
  Hamamatsu S13360-6050PE, 6x6 mm2, 50 um pitch, peak PDE around 40%.
- `Protocolo_de_Experimentación__Validación_del_Timing_Detector_para_el_Experimento_SHiP_CERN_2026.pdf`,
  SOP-3.3: FastIC+ electronics jitter 10 ps sigma.
- User-provided thesis meeting references: 31-Oct-2025 geometry
  1400 x 60 x 10 mm; 13-Apr-2026 post-G4 Gaussian jitter; 17-Apr-2026
  progressive rise time required.

## Branch: `main`

### `src/Materials.cc::CreateEJ200()`

| Property | Canonical | Value in branch | OK? |
|---|---:|---:|---|
| `RINDEX` | 1.58 | `{1.58, 1.58, 1.58, 1.58}` | yes |
| `ABSLENGTH` | 380 cm = 3.8 m | `{3.8*m, 3.8*m, 3.8*m, 3.8*m}` | yes |
| `SCINTILLATIONYIELD` | `10000 / MeV` | `10000.0 / MeV` | yes |
| `SCINTILLATIONTIMECONSTANT1` | `2.1 ns` | `2.1 * ns` | yes |
| `SCINTILLATIONRISETIME1` | `0.9 ns` | absent | **no** |
| `RESOLUTIONSCALE` | 1.0 | `1.0` | yes |
| Emission peak | 425 nm | table has `425 -> 1.000` | yes |

Finding: missing `SCINTILLATIONRISETIME1` is P0. Without it, Geant4 uses
instantaneous scintillation onset, which narrows timing distributions
artificially. Source: Eljen `EJ200.pdf`; user meeting 17-Apr-2026.

### `src/SiPMSD.cc` / `include/SiPMSD.hh`

| Item | Canonical | Value in branch | OK? |
|---|---|---|---|
| PDE table comment | Hamamatsu S13360-6050PE | says `S13360-6025` | **no** |
| PDE numeric table | 33 points, 300-940 nm; peak `0.405 @ 460 nm` | exactly 33 points; peak `0.405 @ 460 nm` | yes |
| Default timing smear | `sqrt(106^2 + 10^2) ps ~= 107 ps` | monolithic `fJitterSigma = 20 ps` | **no** |
| UI command | SPTR + electronics jitter | only `/sipm/jitterSigma` | **no** |

Finding: numerical PDE curve is compatible with the SHiP prototype
S13360-6050PE use case (Betancourt 2020 / W4), but comments incorrectly say
S13360-6025. Timing jitter is a single 20 ps sigma, inconsistent with the
canonical S13360-6050PE SPTR (~106 ps sigma) plus FastIC+ jitter (10 ps sigma).

### `src/DetectorConstruction.cc`

| Item | Canonical/base target | Value in branch | OK? |
|---|---:|---:|---|
| `kBarHalfX` | 700 mm | `700.0 * mm` | yes |
| `kBarHalfY` | 30 mm | `30.0 * mm` | yes |
| `kBarHalfZ` | 5 mm | `5.0 * mm` | yes |
| End SiPM active half-size | 3 mm x 3 mm | `kEndHalfY = kEndHalfZ = 3.0 mm` | yes |
| End SiPM count | 8 per side | `kNEndSiPMs = 8` in header | yes |
| Top SiPM default | 20 at 70 mm pitch | `fNTopSiPMs = 20`, `fTopSiPMPitch = 70.0` | yes |
| Reflector | baseline `R=0.98` explicit panels | `BarPV` direct child; `ReflectorYMinus/XMinus/XPlus/ZMinus/ZPlusPV`; border surfaces `BarPV -> Reflector*PV`; `CreateBarSkinReflector()` R=0.98 | yes for current baseline |
| Mylar active volume | absent in current baseline | absent from active geometry; `CreateMylar()` retained as API reference | yes |
| `G4LogicalSkinSurface` | absent | absent | yes |

Note: the final reflector baseline uses explicit panel volumes including X end
panels. A previous audit checklist requested no X panels, but validation showed
removing them reduced end yield to ~6 ph/event/side; `R=0.98` with X panels gave
accepted end yield. This is physically a high-reflectivity fallback, not a
conservative Al-foil curve.

### `src/PrimaryGeneratorAction.cc`

| Item | Canonical | Value in branch | OK? |
|---|---|---|---|
| Particle | `mu-` | `G4MuonMinus::Definition()` and macros `/gun/particle mu-` | yes |
| Energy | `1 GeV` | `fGun.SetParticleEnergy(1.0 * GeV)` | yes |
| Position | `(0, 0, +60 mm)` | `SetParticlePosition({0,0,60*mm})` | yes |
| Direction | `(0,0,-1)` | `SetParticleMomentumDirection({0,0,-1})`; `/muon/angle` default 0 | yes |
| UI commands | `/muon/angle`, `/muon/midpointSiPMs`, `/muon/gunX` | all present | yes |

### `src/RunAction.cc`

Canonical ntuple columns:

`event_id, face_type, global_id, local_id, time_ns, energy_eV, wl_nm, pde, x_mm, y_mm, z_mm, gun_x_mm`

Branch value: exactly those 12 columns. Additional run summary diagnostics are
present: scintillation photons, total detected photons, detection efficiency,
and boundary census. These are non-ntuple additions and do not break analysis.

## Branch: `feature/edge-scan-and-readout-grouping`

### `src/Materials.cc::CreateEJ200()`

| Property | Canonical | Value in branch | OK? |
|---|---:|---:|---|
| `RINDEX` | 1.58 | `{1.58, 1.58, 1.58, 1.58}` | yes |
| `ABSLENGTH` | 380 cm = 3.8 m | `{3.8*m, 3.8*m, 3.8*m, 3.8*m}` | yes |
| `SCINTILLATIONYIELD` | `10000 / MeV` | `10000.0 / MeV` | yes |
| `SCINTILLATIONTIMECONSTANT1` | `2.1 ns` | `2.1 * ns` | yes |
| `SCINTILLATIONRISETIME1` | `0.9 ns` | absent | **no** |
| `RESOLUTIONSCALE` | 1.0 | `1.0` | yes |
| Emission peak | 425 nm | table has `425 -> 1.000` | yes |

Same P0 missing rise-time issue as `main`.

### `src/SiPMSD.cc` / `include/SiPMSD.hh`

Same as `main` before physics baseline cleanup:

| Item | Canonical | Value in branch | OK? |
|---|---|---|---|
| PDE table comment | S13360-6050PE | says `S13360-6025` | **no** |
| PDE numeric table | 33 points, peak `0.405 @ 460 nm` | same as `main` | yes |
| Default timing smear | `~107 ps` total | monolithic `20 ps` | **no** |
| UI command | SPTR + electronics jitter | only `/sipm/jitterSigma` | **no** |

### `src/DetectorConstruction.cc`

| Item | Canonical/base target | Value in branch | OK? |
|---|---:|---:|---|
| `kBarHalfX/Y/Z` | 700 / 30 / 5 mm | 700 / 30 / 5 mm | yes |
| End SiPM active half-size | 3 mm x 3 mm | 3 mm x 3 mm | yes |
| End SiPM count | 8 per side | 8 per side | yes |
| Top SiPM default | 20 at 70 mm pitch | 20 at 70 mm pitch | yes |
| Reflector | `main` baseline explicit panels R=0.98 | segmented `WrapCenterPV/WrapCap*PV`; Mylar TIR/passive cap model; no Al panel reflector | **no, stale baseline** |
| Mylar active volume | absent in canonical `main` baseline | present as nested wrap volumes | **no, stale baseline** |
| `edgeWrap` feature | feature-specific | configurable `mylar|air|black` edge caps | keep as feature behavior, but must be reconciled with baseline geometry |

Finding: this branch's geometry is intentionally feature-rich for edge scans,
but it is now physically stale relative to `main`. On rebase, the baseline
reflector geometry must win unless the edge-wrap design is intentionally
reimplemented on top of explicit panels.

### `src/PrimaryGeneratorAction.cc`

Same as `main`: `mu-`, `1 GeV`, default position `(0,0,+60 mm)`, direction
`(0,0,-1)`, UI `/muon/angle`, `/muon/midpointSiPMs`, `/muon/gunX`.

### `src/RunAction.cc`

Canonical 12 ntuple columns are present exactly. This branch does not include
the latest `main` run-summary scintillation/boundary diagnostics, because it
has not been rebased after the `main` merge.

## Branch: `feature/sipm-electronics-response`

### `src/Materials.cc::CreateEJ200()`

| Property | Canonical | Value in branch | OK? |
|---|---:|---:|---|
| `RINDEX` | 1.58 | `{1.58, 1.58, 1.58, 1.58}` | yes |
| `ABSLENGTH` | 380 cm = 3.8 m | `{3.8*m, 3.8*m, 3.8*m, 3.8*m}` | yes |
| `SCINTILLATIONYIELD` | `10000 / MeV` | `10000.0 / MeV` | yes |
| `SCINTILLATIONTIMECONSTANT1` | `2.1 ns` | `2.1 * ns` | yes |
| `SCINTILLATIONRISETIME1` | `0.9 ns` | absent for EJ-200 | **no** |
| `RESOLUTIONSCALE` | 1.0 | `1.0` | yes |
| Emission peak | 425 nm | table has `425 -> 1.000` | yes |

Additional material: `CreateEJ230()` exists with rise `0.5 ns`, decay
`1.5 ns`, yield `9700/MeV`, attenuation `120 cm`, peak `391 nm`. These match
`EJ228_EJ230.pdf`. However, this feature changes the default scintillator to
EJ-230 via `fScintillatorName = "EJ230"` and `macros/run.mac` contains
`/det/scintillatorMaterial EJ230`. That is out of scope for the canonical
EJ-200 baseline and should be moved to a later material-selection feature.

### `src/SiPMSD.cc` / `include/SiPMSD.hh`

| Item | Canonical default | Value in branch | OK? |
|---|---|---|---|
| Default sensor | S13360-6050PE prototype | `SensorModel::kBroadcomAFBRS4N66P024M` | **no** |
| `run.mac` sensor | S13360-6050PE prototype | `/sipm/sensorModel Broadcom` | **no** |
| S13360 comment/name | S13360-6050PE | `S13360-6025 (legacy)` | **no** |
| S13360 PDE table | 33 points, peak `0.405 @ 460 nm` | same numeric table | yes |
| SPTR default | 106 ps for S13360-6050PE | `150 ps` for Hamamatsu models | **no** |
| Electronics jitter default | 10 ps FastIC+ | `30 ps` | **no** |
| Quadrature sum | `sqrt(SPTR^2 + electronics^2)` | **not used in timestamp**; only `sptrJitterNs` is added, `fElectronicsSigma` is stored but not applied | **no** |
| UI commands | `/sipm/sptrSigma`, `/sipm/electronicsJitter`, legacy alias | `/sipm/sptrSigma`, `/sipm/electronicsSigma`; no legacy `/sipm/jitterSigma` | partial |
| Ntuple additions | allowed if additive | adds `time_raw_ns`, `sptr_jitter_ns`, `sensor_model_id`, `sensor_sptr_sigma_ns`, `electronics_sigma_ns` | yes, additive |

Additional finding: `src/SiPMSD.cc` appears to contain a duplicated
`G4bool SiPMSD::ProcessHits(...)` signature line in the branch snapshot. If
present in the actual checkout, this is a compile blocker independent of the
physics audit.

### `src/DetectorConstruction.cc`

| Item | Canonical/base target | Value in branch | OK? |
|---|---:|---:|---|
| `kBarHalfX/Y/Z` | 700 / 30 / 5 mm | 700 / 30 / 5 mm | yes |
| End SiPM active half-size | 3 mm x 3 mm | 3 mm x 3 mm | yes |
| End SiPM count | 8 per side | 8 per side | yes |
| Top SiPM default | 20 at 70 mm pitch | 20 at 70 mm pitch | yes |
| Reflector | `main` baseline explicit panels R=0.98 | segmented Mylar wrap, no final panel reflector | **no, stale baseline** |
| Mylar active volume | absent in canonical `main` baseline | present | **no, stale baseline** |
| Default scintillator | EJ-200 | EJ-230 | **no for this audit scope** |
| Material selector | future feature | `/det/scintillatorMaterial EJ230|EJ200` | keep later, but default must not override baseline |

### `src/PrimaryGeneratorAction.cc`

Same as canonical: `mu-`, `1 GeV`, default position `(0,0,+60 mm)`, direction
`(0,0,-1)`, UI `/muon/angle`, `/muon/midpointSiPMs`, `/muon/gunX`.

### `src/RunAction.cc`

The first 4 columns match canonical, but the rest are extended:

`event_id, face_type, global_id, local_id, time_ns, time_raw_ns, sptr_jitter_ns, energy_eV, wl_nm, pde, x_mm, y_mm, z_mm, gun_x_mm, sensor_model_id, sensor_sptr_sigma_ns, electronics_sigma_ns`

This is additive in spirit but changes column indices after `time_ns`, so
existing analysis code that assumes `energy_eV` at column 5 will break unless
it reads branches by name. Requirement: keep canonical branch names and ensure
all analysis uses branch names, not positional column indices.

