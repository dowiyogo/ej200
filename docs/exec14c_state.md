# EXEC_14C reanudación e inventario real

Fecha de reanudación: 2026-06-13 CEST  
Rama: `feat/ej230-sslg4`  
HEAD encontrado: `48420794c3f2bbf44d7dd0544937442088ba7b74`  
Commit recuperable WIP: `533b36024420ac580ad08005fce3c59ce441db11`  
Checkpoint: `checkpoint/pre-exec14c-20260613-1050`

## Estado encontrado

- La reconstrucción Beamer de EXEC_14B estaba presente y compilaba con Python.
- El TeX ya tenía la convención canónica `\graphicspath{{../}}`, 119 frames y las
  21 rutas explícitas de dispersión, pero EXEC_09 aún era una figura EJ-204.
- Las cuatro simulaciones window-dip y End-only `x=-400 mm` ya estaban completas.
- Durante la reanudación terminaron End-only `x=0` y `x=+400 mm`.
- No se encontró simulación activa al cerrar el inventario.
- No quedaron restos `.work_run*`, `.tmp` ni `.bad`; no fue necesario borrar,
  renombrar ni sobrescribir ningún ROOT.

## Inventario de simulaciones especiales

Fuente de validación: `results_ej230_analysis/root_validation_exec14b.csv`.
Cada fila exige TTree `sipm_hits`, entradas positivas, IDs de evento completos,
`gun_x_mm` exacto y el conjunto completo de canales esperado.

| Familia | x [mm] | Eventos | Canales | Entradas | ROOT | Estado |
|---|---:|---:|---:|---:|---|---|
| window-dip EndTop | -652 | 2000 | 0--85 | 9,842,681 | `photon_hits_run_A_x-652mm.root` | válido |
| window-dip EndTop | -642 | 2000 | 0--85 | 9,773,895 | `photon_hits_run_B_x-642mm.root` | válido |
| window-dip EndTop | -648 | 2000 | 0--85 | 9,808,613 | `photon_hits_run_C1_x-648mm.root` | válido |
| window-dip EndTop | -654 | 2000 | 0--85 | 9,886,166 | `photon_hits_run_C2_x-654mm.root` | válido |
| End-only | -400 | 10000 | 0--15 | 10,998,326 | `photon_hits_run007.root` | válido |
| End-only | 0 | 10000 | 0--15 | 7,004,583 | `photon_hits_run015.root` | válido |
| End-only | +400 | 10000 | 0--15 | 10,987,350 | `photon_hits_run023.root` | válido |

El mismo validador confirmó `31/31` ROOT del scan principal, cada uno con 2000
eventos, `gun_x_mm` correcto y canales 0--85. Los ROOT principales se leyeron
desde `/home/reriosto/SHiP/t0minidaq/results_ej230/data/` y no se modificaron.

## Propagación de constantes EJ-230

- Estimadores activos: `analysis/exec07/common.py` define `TAU_R_NS=0.5` y
  `TAU_D_NS=1.5`; `analysis/exec13/common13.py` coincide.
- `summary_exec07.csv` fue regenerado con esas constantes; en `x=0`, nearest-Top
  tiene `Npe=365.9995` y `sigma_est=45.267904 ps`.
- Las 31 figuras `muon_*_geometry.png` fueron regeneradas con el rótulo visible
  `EJ-230 OPSC-106`.
- Los 21 paneles `exec13_tn_*` fueron regenerados con `tau_d=1.5 ns`; son 21
  archivos y 21 hashes distintos.
- Se corrigió un rótulo heredado EJ-204 en el generador de reporte de
  `analysis/exec07/exec11_time_arrival.py`.
- Las demás coincidencias EJ-204/EXEC_12 revisadas pertenecen a analizadores
  históricos, fuentes de paridad EXEC_12 o patrones de traducción de
  `scripts/make_beamer_ej230.py`; no alimentan números ni figuras del reporte
  EJ-230 final.

## Pendiente al iniciar EXEC_14C

1. Regenerar EXEC_09 desde los tres End-only válidos.
2. Regenerar tablas CSV-derived, ejecutar auditorías de assets/números y compilar
   estrictamente.
3. Preflight, render visual, documentación final y push.

## Cierre EXEC_14C

- EXEC_09 se regeneró desde los tres ROOT End-only EJ-230 válidos. La figura
  `exec09_tail_comparison.png` tiene SHA256
  `90d4670d4ec8870747a39a767d397cdcaff15686c2990aa0515db8b581e923c7`,
  distinto de la figura heredada EJ-204.
- Las figuras especiales window-dip también son EJ-230 auténticas:
  `exec08b_window_dip_profiles.png` tiene SHA256
  `77d1573267c6ec3342534c5329a71e7f5ab141c5217451474bc0e3de7babc00c`
  y `exec08b_id18_impact_maps.png` tiene SHA256
  `7e9a9d18fd426f4f9c95f373b8508c175e472defcc47001a38bdfef6b1dfd017`.
- El auditor de assets resolvió `156/156` referencias, sin backreferences ni
  copias raster EJ-204. El auditor numérico verificó 119 frames, 50 `\input{}`
  y 50 tablas regenerables desde CSV, sin diferencias.
- La compilación estricta con `latexmk -halt-on-error` terminó con exit code 0.
  El PDF final tiene 119 páginas, 104 páginas con raster, 156 referencias de
  imagen incrustadas y cero páginas con texto literal `figs/` o cajas vacías
  sospechosas.
- Se renderizaron las 119 páginas a 130 dpi. Se inspeccionaron 5--8, 17--20,
  35, 37, 39, 41--49, 113, 115, 117 y 119; la última corrección compactó las
  etiquetas A/B/C del diagrama de mecanismo en la página 45.
