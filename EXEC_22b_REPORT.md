# EXEC_22b - polishedbackpainted quick test

## Veredicto
NO FUNCIONA

## Cambio probado
- archivo: `src/Materials.cc:356`
- finish anterior: `polished` (EXEC_22, en `CreateMylarReflector`)
- finish nuevo: `polishedbackpainted`
- `REFLECTIVITY`: `src/Materials.cc:366`, valor de corrida `0.95`
- `RINDEX` de superficie: `src/Materials.cc:367`, valor `1.0`
- world/mother `RINDEX` intacto: si, `src/DetectorConstruction.cc:292`
- caras END: no quedaron bloqueadas; hubo hits END_L y END_R en x=0.

Nota: `RunAction` imprime `Reflector R : 0` porque solo mira el mapa legacy `GetReflectorSurfaces()`, que esta vacio tras la migracion a `G4LogicalSkinSurface`; no inspecciona `GetBarSkinSurface()`.

## Auditoria de velocidad
- posiciones: `x = {0, -300, -690} mm`
- eventos: 300 por posicion
- pasos opticos auditados: 19036646 total, 19032419 en centellador
- max `v/(c/n)` en centellador: 1.0
- fraccion superluminica en centellador: 0.0
- volumenes sin `RINDEX`: ninguno
- PASS/FAIL: PASS

## END-only quick physics
| x_mm | N_events | <Npe>_END | <Npe>_END_L | <Npe>_END_R | sigma_END_ps | sigma_boot_ps |
|------|----------|-----------|-------------|-------------|--------------|---------------|
| -690 | 300 | 908.9467 | 908.9100 | 0.0367 | n/a | n/a |
| -300 | 300 | 1.6633 | 1.4633 | 0.2000 | n/a | n/a |
| 0 | 300 | 0.7633 | 0.3733 | 0.3900 | n/a | n/a |

`sigma_END_ps` es `n/a` porque hubo 0 eventos con ambos extremos disparando el leading-edge SUM4; no hay muestra valida para `t_avg = (t_L + t_R)/2`.

## Antes/despues en x=0
- Antes EXEC_22: `<Npe>_END(0) ~= 0.38 PE`, `sigma_END(0) ~= 985 ps`
- Despues EXEC_22b: `<Npe>_END(0) = 0.7633 PE`, `sigma_END(0) = null`

## Variante unica opcional
Se probo `groundbackpainted` solo en `x=0`, 300 eventos, manteniendo `dielectric_dielectric`, `unified`, `REFLECTIVITY=0.95`, `RINDEX=1.0`.

- velocidad: PASS, max `v/(c/n)=1.0`, fraccion superluminica 0.0
- `<Npe>_END(0) = 0.75 PE`
- `<Npe>_END_L(0) = 0.38 PE`
- `<Npe>_END_R(0) = 0.37 PE`
- `sigma_END(0) = null` por 0 triggers de ambos extremos

## Interpretacion
El modelo backpainted no recupero la luz no-TIR en el centro; siguiente paso es geometria explicita air gap + Mylar.

## Rutas
- CSV Npe: `out/EXEC_22b/end_npe_summary.csv`
- CSV sigma: `out/EXEC_22b/end_sigma_summary.csv`
- CSV auditoria velocidad: `out/EXEC_22b/speed_audit.csv`
- JSON velocidad: `out/EXEC_22b/speed_audit_summary.json`
- JSON resultados: `out/EXEC_22b/results_exec22b.json`
- figura velocidad: `out/EXEC_22b/speed_ratio_hist.png`
- ROOT mini-corrida: `out/EXEC_22b/endonly_mylar_exec22b_20260621_124904/`
- variante ground: `out/EXEC_22b/ground_variant/`
- logs: `out/EXEC_22b/endonly_mylar_exec22b_20260621_124904/run_x0mm.log` y macros en `work_x*mm/`

## Git
- branch: `exec21-optfix`
- tag pre: `EXEC_22b-pre-20260621_124414`
- commit hash base: `9b8361f`
- arbol sucio previo no tocado: `presentation/endonly_mylar/main.tex`, `analysis/bootstrap_attenuation`, ficheros LaTeX auxiliares bajo `presentation/endonly_mylar/`
- comando push NO ejecutado: `git push origin exec21-optfix`

## Verificacion tecnica
- `cmake --build build -j16`: OK
- `./build/endonly_geometry_check`: OK
- `root` no estaba disponible como comando durante `scripts/run_scan.sh`; los ROOT se validaron con `uproot`.
