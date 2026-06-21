# EXEC_24 - PE budget and SiPM detection audit

## Veredicto

COMPLETADO-CON-FLAGS.

## Resumen ejecutivo

- Los 2035.26 hits/event de EXEC_23 no vienen de doble conteo: `n_rows_raw == n_unique_tracks` y no hay multihit por `(event_id, track_id, global_id)`.
- `SiPMSD.cc` no aplica una Bernoulli explicita de PDE; guarda la rama `pde`. Pero la superficie `SiPMSurface` si configura `EFFICIENCY=PDE`, y el censo `Bar -> SiPM = 163638` frente a `hits = 101763` da una razon 0.6219, compatible con `mean(pde)=0.6162`.
- Por tanto, `sum(pde)=1254.22/event` se reporta como reponderacion diagnostica, no como PE baseline si la eficiencia de superficie esta activa.
- El exceso restante es fisico/optico en esta geometria: de 1.046886e6 fotones generados se registran 101763 detecciones, una eficiencia raw de 9.72% sobre ambos extremos. Con 48% de cobertura activa por extremo y reflector lateral idealizado, esto sigue siendo demasiado alto para promover a baseline.
- No se aplico correccion de simulacion: otra Bernoulli en `SiPMSD` probablemente doble-contaria PDE. No se toco EndTop, no hubo scans largos, no se tuneo reflectividad.

## T1 PE accounting

| metric | value |
|---|---:|
| events | 50 |
| raw rows/event | 2035.26 |
| unique tracks/event | 2035.26 |
| sum_pde/event | 1254.22 |
| duplicates_per_track mean | 1.000 |
| duplicates_per_track max | 1.000 |
| multihit `(event,track,gid)` | 0 |
| END_L raw/event | 1014.48 |
| END_R raw/event | 1020.78 |
| END_L sum_pde/event | 625.28 |
| END_R sum_pde/event | 628.95 |
| Bar -> SiPM entering/event | 3272.76 |
| recorded rows / Bar -> SiPM entering | 0.6219 |

Gate: no duplicate-detection artifact found. PDE handling is split: surface-level `EFFICIENCY` likely gates detections, while `SiPMSD` only records `pde`.

## T2 SiPMSD audit

- PDE aplicada: yes by intended Geant4 optical surface `EFFICIENCY`; no explicit Bernoulli inside `SiPMSD.cc`.
- Track killed after detection: yes.
- One photon = one hit: yes in the 50-event audit.
- Detector area fisica: yes, explicit 8 x 6x6 mm SiPM volumes per end.
- `global_id`: correct for END-only, copied from physical volume IDs 0-7 and 8-15.

Detailed note: `out/EXEC_24/sipmsd_audit.md`.

## T3 Photon budget

| metric | value |
|---|---:|
| LY used | 10400 ph/MeV |
| optical photons generated | 1046886 |
| generated/event | 20937.72 |
| estimated Edep/event | 2.013 MeV |
| collection efficiency raw | 9.7205% |
| collection efficiency unique | 9.7205% |
| collection efficiency sum_pde diagnostic | 5.9903% |

Gate: raw detection efficiency is high. Even the diagnostic `sum_pde/Ngen` remains high for a real baseline, but should not be treated as a second PDE application without an explicit Geant4 surface validation.

## T4 SiPM geometry

| metric | value |
|---|---:|
| bar end face area | 600 mm2 |
| active area per SiPM | 36 mm2 |
| SiPMs per end | 8 |
| active SiPM area per end | 288 mm2 |
| active coverage per end | 48% |

The geometry is not a full-face detector. It is, however, a very large active coverage per end, and the present optical setup couples many reflected photons back to those active patches.

## T5 Optical lifetime

| metric | value |
|---|---:|
| detected photons | 101763 |
| total optical steps | 90395122 |
| steps/generated photon | 86.35 |
| detected time mean | 6.12 ns |
| detected time p95 | 10.86 ns |
| detected time max | 42.06 ns |
| detected TIR mean | 52.99 |
| detected scint-air reflection mean | 53.22 |

No long-lived photon concern was found for detected photons (`p95 < 50 ns`). The 500-event cost is still high because the total optical step count is large, not because detected photons have a very late time tail.

## Corrections

No simulation correction was applied.

| trial | applied | reason |
|---|---|---|
| diagnostic_only | no | no duplicate hits; PDE appears already applied by Geant4 surface `EFFICIENCY`; another Bernoulli in `SiPMSD` would likely double-apply PDE |

## Resultado final rapido

| x_mm | events | Npe_raw | Npe_unique | Npe_pde_expected | paired_fraction | sigma_END_ps |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 50 | 2035.26 | 2035.26 | 2035.26 | 1.000 | 60.12 |

Diagnostic reweighting: `sum(pde) = 1254.22/event`.

## Flags

- `npe_raw_high_concern`
- `sum_pde_reweight_high_concern`
- `collection_eff_moderate_concern`
- `formal_500_event_center_run_interrupted_after_10min` inherited from EXEC_23

## Git

- Branch: `exec24-pe-budget-audit`
- Pre tag: `EXEC_24-pre-20260621_163654`
- Commit: pending at report creation
- Push printed, not executed

## Proximo paso

Before EXEC_25 scan, validate the Geant4 SiPM surface semantics with a tiny controlled PDE test, for example PDE forced to 0 and 1 on a small x=0 run. If that confirms surface gating, the next physical decision is not PDE counting but the optical model: active coverage, reflector realism, air-gap coupling, and whether a physically motivated time window belongs in the analysis.
