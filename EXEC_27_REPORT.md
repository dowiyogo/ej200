# EXEC_27 - edge estimator vs SiPM coupling audit

## Veredicto

COMPLETADO-CON-FLAGS

## Resumen ejecutivo

- El fallo de borde observado en EXEC_26 es principalmente del estimador `t_avg`, no de la optica/acoplo SiPM.
- Con la configuracion candidata de EXEC_26 (`ground`, `sigma_alpha=0.20 rad`, `R=0.95`, `surface_loss=0.10`), `t_avg` falla en los extremos porque promedia un lado con miles de PE y otro con pocos PE.
- `wmean_static` y `GLS_residual` rescatan todos los puntos T1: `x=-690,-400,0,+400,+690 mm`.
- No se ejecuto T2 angular ni T3 coupling: el gate T1 cerro el diagnostico y no justifica tocar el acoplo SiPM.
- Hay candidato para scan11, pero con estimadores calibrados weighted/GLS, no con `t_avg` como metrica principal de borde.

## T1 estimator decomposition

| x mm | Npe_L | Npe_R | sigma_L ps | sigma_R ps | sigma_tavg ps | sigma_wmean_static ps | sigma_wmean_event ps | sigma_GLS ps |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -690 | 2763.90 | 4.58 | 2.43 | 654.18 | 327.33 | 2.43 | 5.90 | 3.00 |
| -400 | 303.74 | 14.92 | 33.14 | 323.25 | 166.74 | 35.31 | 46.09 | 40.69 |
| 0 | 55.54 | 56.87 | 148.28 | 145.95 | 96.22 | 91.85 | 101.77 | 115.01 |
| 400 | 13.74 | 305.22 | 340.12 | 30.89 | 169.58 | 30.47 | 40.87 | 45.26 |
| 690 | 5.02 | 2888.82 | 578.45 | 1.66 | 289.26 | 1.65 | 4.76 | 2.52 |

Gate T1:

- `sigma_tavg` fails at `x=-690` and `x=+690`.
- `sigma_wmean_static < 250 ps` at all five positions.
- `sigma_GLS_residual < 250 ps` at all five positions.
- All rows pass the speed gate: `missing_rindex=0`, `max_superluminal_fraction=0`.

Diagnosis: the edge problem is the equal-weight timestamp estimator under extreme END_L/END_R light-sharing asymmetry.

## T2 SiPM angular audit

Not run by design.

Reason: T1 satisfies the EXEC_27 gate that weighted/GLS estimators rescue all edge points. Therefore there is no clear evidence requiring SiPM angular/coupling instrumentation in this execution.

This does not prove the SiPM coupling model is fully realistic; it shows the EXEC_26 edge timing failure is not sufficient evidence for changing it.

## T3 coupling model

No coupling model was implemented or scanned.

The PDE surface semantics remain unchanged. No extra quantum-efficiency factor, LY change, ABSLENGTH change, SPTR, FastIC jitter, EndTop, or position-dependent loss was introduced.

## T4 comparability

The historical 85-90 ps number should be compared to calibrated weighted or GLS-like endpoint combinations, not to unweighted `t_avg` at strongly asymmetric edge positions.

Summary:

- `t_avg`: useful at center, fragile at edges.
- `wmean_static`: closest comparable estimator if historical analysis used calibrated endpoint timing weights.
- `wmean_event`: diagnostic amplitude-weighted estimator; it naturally follows the near-end timing at edges.
- `GLS_residual`: intrinsic timing resolution after removing the per-position mean; good for resolution diagnostics.
- near-end only: excellent local timing diagnostic at edges, not a two-ended ToF estimator.
- far-end only: transport/far-side detectability diagnostic, poor at edges with low PE.

See `out/EXEC_27/comparability_matrix.md`.

## Candidate

Existing EXEC_26 optical configuration remains the candidate:

```text
scint-air finish = ground
sigma_alpha = 0.20 rad
reflector R = 0.95
surface_loss = 0.10
```

Recommended timing estimator for scan11:

```text
wmean_static or GLS_residual with position calibration
```

Validation summary:

- Center `x=0`: `Npe_END=112.41`, `sigma_wmean_static=91.85 ps`, `sigma_GLS=115.01 ps`.
- Edge `x=-690`: near side `2763.90 PE`, far side `4.58 PE`, `sigma_wmean_static=2.43 ps`, `sigma_GLS=3.00 ps`.
- Edge `x=+690`: near side `2888.82 PE`, far side `5.02 PE`, `sigma_wmean_static=1.65 ps`, `sigma_GLS=2.52 ps`.

The huge near-end PE at the edges remains a geometry/light-sharing feature to document in scan11, but it no longer blocks the timing metric once equal weighting is removed.

## Flags

- `tavg_edge_failure`
- `weighted_estimators_rescue_edges`
- `gls_residual_rescues_edges`

## Rutas

- `out/EXEC_27/edge_estimator_summary.csv`
- `out/EXEC_27/edge_estimator_summary.json`
- `out/EXEC_27/candidate_validation.csv`
- `out/EXEC_27/candidate_validation.json`
- `out/EXEC_27/comparability_matrix.md`
- `out/EXEC_27/results_exec27.json`
- `out/EXEC_27/png/sigma_estimators_vs_x.png`
- `out/EXEC_27/png/npe_lr_vs_x.png`
- `out/EXEC_27/png/weights_vs_x.png`

Raw ROOT/logs are under `out/EXEC_27/raw/` and are intentionally not committed.

## Git

- Branch: `exec27-edge-estimator-coupling`
- Tag: `EXEC_27-pre-20260621_212620`
- Commit: pending at report-writing time
- Push printed, not executed:

```bash
git push origin exec27-edge-estimator-coupling
```

## Proximo paso

EXEC_28 should run scan11 END-only with calibrated `wmean_static` and GLS/residual estimators. Do not add SiPM coupling unless a later angular audit shows a separate, estimator-independent pathology.
