# EXEC_24 SiPMSD audit

Source inspected:

- `src/SiPMSD.cc`
- `src/DetectorConstruction.cc`
- `src/Materials.cc`
- `out/EXEC_23/nominal_quick50/run_x0mm.log`

## Findings

| question | answer | evidence |
|---|---|---|
| PDE applied in `SiPMSD.cc` as Bernoulli? | no | `SiPMSD::ProcessHits` computes `pde = GetPDE(energy)` and stores it in ntuple column 7, but has no random accept/reject draw. |
| PDE applied elsewhere? | likely yes, via Geant4 optical surface | `CreateSiPMSurface` sets `EFFICIENCY` to the selected PDE curve and `REFLECTIVITY=0`. The EXEC_23 log has 163638 Bar->SiPM entering crossings and 101763 recorded detections, ratio 0.621879, compatible with mean recorded-row PDE 0.616248. |
| Each TTree row is a PE accepted or an incident photon? | row is treated as an accepted SiPM detection by current Geant4 surface setup | The SD is invoked after the SiPM optical surface; applying another Bernoulli in `SiPMSD` would probably double-apply PDE. |
| Track killed after detection? | yes | `SiPMSD.cc` calls `track->SetTrackStatus(fStopAndKill)` after `AddNtupleRow` and `OpticalAudit::RecordDetection`. |
| One photon can produce multiple entries? | not observed | 50-event audit found `duplicates_per_track_mean = 1.0`, `duplicates_per_track_max = 1.0`, and `n_multihit_event_track_gid_total = 0`. |
| Detection depends on real active area or full END face? | real explicit SiPM volumes, not full face | `DetectorConstruction.cc` places 8 `EndSiPMLV` 6x6 mm volumes per end and attaches border surfaces to each physical volume. |
| `global_id` assignment correct? | yes for END-only | left IDs 0-7 and right IDs 8-15 come from physical-volume copy numbers. |

## Interpretation

The 2035.26 rows/event are not caused by duplicate track counting. The code does not apply PDE inside `SiPMSD`, but the SiPM optical border surface is configured to apply PDE through `EFFICIENCY`. The boundary census supports that interpretation: only 62.19% of Bar->SiPM entering crossings become recorded rows.

Therefore `sum(pde) = 1254.22/event` is retained only as a diagnostic reweighting of already recorded rows. It should not be promoted as the baseline PE count unless a separate validation proves the surface-level `EFFICIENCY` is not gating SD invocation in this geometry.

The remaining high Npe is dominated by optical collection/acceptance: 101763 recorded detections out of 1046886 generated photons, or 9.72% raw detection efficiency over both ends.
