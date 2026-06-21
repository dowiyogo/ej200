# EXEC_23 - explicit air gap + reflector optical fix

## Veredicto
COMPLETADO-CON-FLAGS

## Resumen ejecutivo
- Se implemento geometria explicita: scintillator OPSC-101/EJ-204 -> air gap optico -> paneles externos Mylar.
- EXEC_22b quedo descartado porque `polishedbackpainted` y `groundbackpainted` dejaron `<Npe>_END(0) ~= 0.75 PE`.
- Resultado principal en x=0, muestra diagnostica de 50 eventos:
  - `<Npe>_END` antes EXEC_22b: 0.763 PE
  - `<Npe>_END` despues EXEC_23: 2035.26 PE
  - `sigma_END`: 60.12 ps, con bootstrap inestable por muestra pequena
  - speed gate: PASS
- Flags: `npe_high_concern`, `formal_500_event_center_run_interrupted_after_10min`, `scan11_not_run_due_flags`.

## Cambio geometrico
- `src/DetectorConstruction.cc:31`: `air gap = 0.10 mm`, `reflector = 0.05 mm`.
- `src/DetectorConstruction.cc:323`: modo `mylar` ya no instala `G4LogicalSkinSurface` en `BarLV`.
- `src/DetectorConstruction.cc:327-406`: cuatro paneles laterales `AirGap*PV` y cuatro `Mylar*PV`.
- material air gap: `G4_AIR` con `RINDEX=1.0`.
- material reflector: `G4_MYLAR` con `RINDEX=1.65` y `ABSLENGTH=1 um`.
- surfaces:
  - scintillator-air: `dielectric_dielectric/unified/polished`, sin `REFLECTIVITY`.
  - air-reflector: `dielectric_metal/unified/polished`, `REFLECTIVITY=0.95`, `EFFICIENCY=0`.
- caras END/SiPM: los paneles tienen longitud X igual a la barra y solo cubren ±Y/±Z; no hay material delante de ±X.

## Materiales y RINDEX
| volume/material | role | has RINDEX | RINDEX range | notes |
|---|---|---:|---|---|
| WorldLV / G4_AIR | world | yes | 1.0-1.0 | mother optical air |
| BarLV / opsc-101 | scintillator | yes | 1.58-1.58 | guard agregado en DetectorConstruction |
| AirGap*LV / G4_AIR | air_gap | yes | 1.0-1.0 | gap explicito |
| Mylar*LV / G4_MYLAR | reflector | yes | 1.65-1.65 | reflector externo |
| EndSiPMLV / G4_SILICON_DIOXIDE | sipm | yes | 1.58-1.58 | coupling/detector |

## Geometry check
- overlaps: PASS.
- SiPM END intactos: 16 END, 0 TOP en modo END-only.
- reflector no tapa caras de lectura: PASS por construccion y por hits END_L/END_R.
- dimensiones: `out/EXEC_23/geometry_dimensions.json`.

## Auditoria de velocidad
| material | n_steps | max v/(c/n) | superluminal fraction | verdict |
|---|---:|---:|---:|---|
| opsc-101 | 39427292 | 1.0 | 0.0 | PASS |
| G4_AIR | 44824618 | 1.0 | 0.0 | PASS |
| G4_MYLAR | 6141620 | 1.0 | 0.0 | PASS |
| G4_SILICON_DIOXIDE | 1592 | 1.0 | 0.0 | PASS |

Volumenes sin RINDEX: ninguno.

## Boundary/process audit
Categorias registradas: 78. Extractos relevantes:
- Scintillator/AirGap: millones de `FresnelRefraction`, `FresnelReflection` y `StepTooSmall`; TIR detectada via historia por foton.
- AirGap/Reflector: interacciones presentes; `fraction_visited_reflector_boundary = 0.0375`.
- `NoRINDEX`: no domina; no hay volumen relevante sin RINDEX en speed gate.
- Historia detectada en x=0: 101763 fotones detectados.

## Mini-validacion optica
| x_mm | events | <Npe>_END | <Npe>_L | <Npe>_R | paired fraction | sigma_END_ps | bootstrap_ps |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 50 | 2035.26 | 1014.48 | 1020.78 | 1.0 | 60.12 | 10338.41 |

La corrida formal de 500 eventos en x=0 fue interrumpida tras ~10 min porque no terminaba. El ROOT parcial crecio hasta ~12 MB antes de abortar; no se uso como resultado.

## Antes/despues
- EXEC_21: `<Npe>~0.37`, `sigma~600-880 ps`.
- EXEC_22b: `<Npe>~0.76`, sin `sigma_END` medible.
- EXEC_23: `<Npe>=2035.26`, `sigma_END=60.12 ps` en muestra diagnostica de 50 eventos.

## Sensibilidad reflectividad
No se corrio. Motivo: el nominal ya dispara `npe_high_concern` y la corrida formal de 500 eventos en centro excedio el wall-clock practico.

## Scan11 / Scan31
No se corrieron. Motivo: no promover una configuracion con `Npe > 300 PE` y coste alto como baseline.

## EndTop mini-port
No se porto. Motivo: END-only no queda listo como baseline fisico; requiere decision sobre cobertura/reflectividad/cortes de transporte antes de tocar EndTop.

## Flags
- `npe_high_concern`: `<Npe>_END(0)=2035.26 PE`, muy por encima del umbral de preocupacion de 300 PE.
- `formal_500_event_center_run_interrupted_after_10min`: la validacion nominal de 500 eventos no completo en tiempo razonable.
- `scan11_not_run_due_flags`: no se lanzo scan11/scan31.

## Rutas
- geometria: `out/EXEC_23/geometry_check.txt`
- volumenes: `out/EXEC_23/geometry_volumes.csv`
- dimensiones: `out/EXEC_23/geometry_dimensions.json`
- velocidad: `out/EXEC_23/speed_audit_steps.csv`
- velocidad JSON: `out/EXEC_23/speed_audit_summary.json`
- boundary counts: `out/EXEC_23/optical_boundary_counts.csv`
- photon history: `out/EXEC_23/detected_photon_history.csv`
- historia JSON: `out/EXEC_23/optical_history_summary.json`
- Npe: `out/EXEC_23/end_npe_summary.csv`
- pairing: `out/EXEC_23/end_pairing_summary.csv`
- sigma: `out/EXEC_23/end_sigma_summary.csv`
- resultados: `out/EXEC_23/results_exec23.json`
- figuras: `out/EXEC_23/png/`, `out/EXEC_23/speed_ratio_by_material.png`
- logs: `out/EXEC_23/nominal_quick50/run_x0mm.log`

## Git
- branch: `exec23-explicit-airgap`
- tag pre: `EXEC_23-pre-20260621_132258`
- commit: ver `git rev-parse --short HEAD` en la rama `exec23-explicit-airgap`
- push command printed, not executed: `git push origin exec23-explicit-airgap`

## Proximo paso
No usar aun como baseline EndTop. El siguiente paso honesto es decidir si el exceso de luz se debe a cobertura/reflectividad idealizada o a falta de un limite fisico/operativo de transporte optico, y entonces repetir T5 nominal.
