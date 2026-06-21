# EXEC_26 - scintillator-air surface realism

## Veredicto

COMPLETADO-CON-FLAGS.

## Resumen ejecutivo

- La frontera dominante si era scintillator-air: cambiarla de polished ideal a ground reduce Npe de O(1800) a O(20-200) con reflector apagado.
- La rugosidad scint-air por si sola controla el centro con `R_reflector=0`: `ground sigma_alpha=0.20 rad` da 188.9 PE/event, frente a 2035 PE en EXEC_23 y 1781.95 PE en el minimo de EXEC_25.
- Al reintroducir reflector, la coleccion vuelve a subir: `R=0.85` da 762.35 PE/event y `R=0.95` da 1326.05 PE/event.
- Se necesito un modelo efectivo de perdida/contacto parcial. Con `ground sigma_alpha=0.20 rad`, `R=0.95`, y `surface_loss=0.10`, el centro queda en 107.55 PE/event en 20 eventos y 112.41 PE/event en 100 eventos.
- La validacion de bordes falla como candidato de scan: en x=-690 mm y x=+690 mm la coleccion vuelve a ~2768-2894 PE/event, dominada por el extremo cercano.
- Por tanto, hay candidato central para estudiar el mecanismo, pero no hay candidato listo para EXEC_27 scan11. Antes del scan hace falta revisar respuesta de borde o introducir un modelo de acoplo/contacto dependiente de posicion/geometria.

## T2 scint-air roughness, reflector off

| finish | sigma_alpha rad | R_reflector | Npe | recorded/generated | paired fraction | sigma ps | speed |
|---|---:|---:|---:|---:|---:|---:|---|
| polished | 0.00 | 0.00 | 1837.90 | 0.09415 | 1.00 | 58.38 | PASS |
| ground | 0.05 | 0.00 | 977.30 | 0.04498 | 1.00 | 48.66 | PASS |
| ground | 0.10 | 0.00 | 634.05 | 0.03018 | 1.00 | 30.30 | PASS |
| ground | 0.20 | 0.00 | 188.90 | 0.00998 | 1.00 | 64.68 | PASS |
| ground | 0.30 | 0.00 | 62.85 | 0.00325 | 1.00 | 129.77 | PASS |
| ground | 0.50 | 0.00 | 20.25 | 0.00098 | 1.00 | 246.03 | PASS |

Boundary trend: `ground` largely destroys the ideal TIR guide. For example polished has 12.48M total-internal-reflection boundary statuses in 20 events; `ground sigma=0.20` has only 1234 TIR statuses and many more refractions/absorptions.

## T3 reflector reintroduced

Using `ground sigma_alpha=0.20 rad`:

| finish | sigma_alpha rad | R_reflector | Npe | sigma ps | speed |
|---|---:|---:|---:|---:|---|
| ground | 0.20 | 0.85 | 762.35 | 50.41 | PASS |
| ground | 0.20 | 0.95 | 1326.05 | 46.45 | PASS |

Interpretation: roughness controls the TIR-only case, but realistic reflector recovery still over-collects unless additional contact/loss realism is included.

## T4 surface loss/contact model

Effective model: kill a photon with probability `surface_loss` at a scintillator-air surface encounter. This is explicitly a contact/loss model, not a material retune.

Parameters held fixed: `ground sigma_alpha=0.20 rad`, `R_reflector=0.95`.

| model | loss/contact fraction | Npe | sigma ps | paired fraction | speed |
|---|---:|---:|---:|---:|---|
| surface_loss | 0.0 | 1308.55 | 38.99 | 1.00 | PASS |
| surface_loss | 0.1 | 107.55 | 83.49 | 1.00 | PASS |
| surface_loss | 0.3 | 13.55 | 248.62 | 1.00 | PASS |
| surface_loss | 0.5 | 4.05 | 354.72 | 0.80 | PASS |
| surface_loss | 0.7 | 2.10 | 657.30 | 0.50 | PASS |
| surface_loss | 0.9 | 1.15 | n/a | 0.10 | PASS |

Chosen center configuration: `surface_loss=0.10`, the smallest tested loss that brings the center into the 20-300 PE band.

## Candidate

Center candidate:

- name: `loss_R095_0p1`
- scint-air finish: `ground`
- scint-air sigma alpha: `0.20 rad`
- reflector R: `0.95`
- surface loss: `0.10`
- center Npe: `107.55` in 20 events
- center sigma: `83.49 ps`
- paired fraction: `1.0`
- speed gate: PASS

Why it is physical and not tuning: the dominant correction is on the boundary identified by EXEC_25 as the missing physics, and `surface_loss=0.10` is a small effective contact/loss probability compared with the larger values tested. It is selected as the first point in the physically acceptable band, not to match 40 PE.

Comparison:

| reference | Npe center |
|---|---:|
| EXEC_23 | 2035.26 |
| EXEC_25 lowest | 1781.95 |
| EXEC_26 roughness-only center | 188.90 |
| EXEC_26 center candidate | 107.55 |
| EXEC_26 100-event center validation | 112.41 |

## Candidate validation

| x_mm | events | Npe | Npe_L | Npe_R | paired fraction | sigma ps | runtime s | speed |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 100 | 112.41 | 55.54 | 56.87 | 1.00 | 96.22 | 293.96 | PASS |
| -690 | 50 | 2768.48 | 2763.90 | 4.58 | 1.00 | 327.33 | 101.11 | PASS |
| 690 | 50 | 2893.84 | 5.02 | 2888.82 | 1.00 | 289.26 | 99.33 | PASS |

The center validates. The edge points do not: local SiPM proximity dominates and returns O(2800) PE. This is the main EXEC_26 residual flag.

## Flags

- `surface_loss_model_used`
- `edge_validation_high_collection`

## Rutas

- `out/EXEC_26/scint_air_roughness_R0.csv`
- `out/EXEC_26/scint_air_roughness_R0.json`
- `out/EXEC_26/scint_air_plus_reflector.csv`
- `out/EXEC_26/scint_air_plus_reflector.json`
- `out/EXEC_26/surface_loss_bracket.csv`
- `out/EXEC_26/surface_loss_bracket.json`
- `out/EXEC_26/candidate_speed_audit.csv`
- `out/EXEC_26/candidate_speed_audit.json`
- `out/EXEC_26/candidate_validation.csv`
- `out/EXEC_26/candidate_validation.json`
- `out/EXEC_26/results_exec26.json`
- `out/EXEC_26/png/scint_air_roughness_R0.png`
- `out/EXEC_26/png/scint_air_plus_reflector.png`
- `out/EXEC_26/png/surface_loss_bracket.png`
- `out/EXEC_26/png/candidate_time_pdf_x0.png`

Raw ROOT outputs remain under `out/EXEC_26/raw/` for local inspection and are not staged for commit.

## Git

- Branch: `exec26-scint-air-surface-realism`
- Tag: `EXEC_26-pre-20260621_185850`
- Commit: pending at report creation
- Push printed, not executed

## Proximo paso

Do not launch scan11 yet. Review the edge response first: either introduce explicit position/geometric coupling realism near the END SiPMs, reduce/segment effective active area near the end faces, or add a position-dependent contact/wrapping model. After the edge response is controlled, promote the center candidate into an EXEC_27 END-only scan.
