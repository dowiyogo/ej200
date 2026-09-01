# Informe de estado de ramas

Fecha del informe: **2026-08-05**  
Repositorio local: `/home/reriosto/SHiP/ej200`  
Remoto GitHub: `git@github.com:dowiyogo/ej200.git`  
Snapshot remoto: despues de `git fetch --all --prune`  
Rama activa en este worktree: **`exec22-endtop-optfix`** (`76b582c`, 2026-06-21 09:15:26 +02:00)

## Resumen ejecutivo

- Hay **20 ramas logicas** contando ramas locales y remotas de `origin`.
- **13 ramas existen en GitHub**. Todas las ramas que existen tanto localmente como en `origin` estan sincronizadas commit a commit despues del `fetch`.
- Hay **8 ramas locales sin rama remota**: `exec21-optfix`, `exec22-endtop-optfix`, `exec23-explicit-airgap`, `exec24-pe-budget-audit`, `exec25-optical-realism-bracket`, `exec26-scint-air-surface-realism`, `exec27-edge-estimator-coupling` y `exec28-endonly-scan11-weighted`.
- Hay **1 rama remota sin rama local**: `origin/feat/ej230-endonly-mylar`. Esta fue actualizada en GitHub durante el `fetch`, de `7e18117` a `e4110ce`.
- La rama mas reciente por fecha de commit es un empate local: **`exec27-edge-estimator-coupling`** y **`exec28-endonly-scan11-weighted`**, ambas en `bd24951`, hechas el **2026-06-21 21:37:53 +02:00**.
- La rama mas reciente publicada en GitHub es **`origin/feat/ej230-endonly-mylar`**, en `e4110ce`, hecha el **2026-06-18 19:39:31 +02:00**.
- `main` local y `origin/main` coinciden en `84e902c` del **2026-06-18 10:29:27 +02:00**.

Notas:

- En la columna `vs origin`, `0/0` significa local y remoto identicos. `local-only` significa que no existe en GitHub. `remote-only` significa que solo existe en GitHub.
- En la columna `vs main`, el formato es `behind/ahead`: commits que tiene `main` y no la rama / commits que tiene la rama y no `main`.
- Las ramas `+` en `git branch -vv` estan actualmente bloqueadas por otros worktrees; no es un problema de sincronizacion, solo indica que estan checkout en otra ruta.

## Matriz de estado

| Rama | Local | GitHub | Ultimo commit | Fecha | vs origin | vs main | Estado |
|---|---:|---:|---|---:|---:|---:|---|
| `backup/photon-budget-worktree` | `81e242d` | `81e242d` | Backup worktree: capture analysis and run changes (photon budget diagnostics) | 2026-06-09 14:29:32 +02:00 | 0/0 | 1/12 | Publicada y sincronizada |
| `diag/photon-budget` | `0ca6855` | `0ca6855` | fix(optics): replace reflector volumes with bar skin surface | 2026-06-18 10:29:23 +02:00 | 0/0 | 1/14 | Publicada y sincronizada |
| `exec21-optfix` | `9b8361f` | - | fix(optics): add REFLECTIVITY=0.95 to dielectric_dielectric surface for Mylar model | 2026-06-21 09:07:38 +02:00 | local-only | 1/63 | No publicada |
| `exec22-endtop-optfix` | `76b582c` | - | fix(optics): port exec22 TIR fix to EndTop build (CreateBarSkinReflector) | 2026-06-21 09:15:26 +02:00 | local-only | 1/46 | No publicada, rama activa |
| `exec23-explicit-airgap` | `bd78211` | - | EXEC_23: implement explicit air gap reflector geometry | 2026-06-21 13:54:29 +02:00 | local-only | 1/64 | No publicada |
| `exec24-pe-budget-audit` | `e20b47d` | - | EXEC_24: audit PE budget and SiPM detection after explicit air gap fix | 2026-06-21 16:45:49 +02:00 | local-only | 1/65 | No publicada |
| `exec25-optical-realism-bracket` | `7c55326` | - | EXEC_25: validate PDE surface semantics and bracket optical collection realism | 2026-06-21 18:14:25 +02:00 | local-only | 1/66 | No publicada |
| `exec26-scint-air-surface-realism` | `3148da1` | - | EXEC_26: bracket scintillator-air surface realism for END light collection | 2026-06-21 19:41:07 +02:00 | local-only | 1/67 | No publicada |
| `exec27-edge-estimator-coupling` | `bd24951` | - | EXEC_27: diagnose edge response via END estimators and SiPM coupling | 2026-06-21 21:37:53 +02:00 | local-only | 1/68 | No publicada, ultima actividad local |
| `exec28-endonly-scan11-weighted` | `bd24951` | - | EXEC_27: diagnose edge response via END estimators and SiPM coupling | 2026-06-21 21:37:53 +02:00 | local-only | 1/68 | No publicada, mismo commit que `exec27` |
| `exp/pair-scan-2026-06-11` | `14ae395` | `14ae395` | fix(tests): replace escape-fraction guard with sipm-entry guard | 2026-06-18 11:04:37 +02:00 | 0/0 | 1/67 | Publicada y sincronizada |
| `feat/ej204-event-display-tracks` | `0a42656` | `0a42656` | fix(tests): replace escape-fraction guard with sipm-entry guard | 2026-06-18 11:04:38 +02:00 | 0/0 | 1/46 | Publicada y sincronizada, sin upstream configurado localmente |
| `feat/ej230-endonly-mylar` | - | `e4110ce` | DIAGNOSTIC: Add PrintGroupVelDiagnostic() to investigate GROUPVEL propagation velocity anomaly | 2026-06-18 19:39:31 +02:00 | remote-only | 1/78 | Solo en GitHub |
| `feat/ej230-sslg4` | `ca2f1c3` | `ca2f1c3` | fix(tests): replace escape-fraction guard with sipm-entry guard | 2026-06-18 11:04:36 +02:00 | 0/0 | 1/70 | Publicada y sincronizada |
| `feat/endonly-mylar` | `6942ee6` | `6942ee6` | fix(optics): replace reflector volumes with bar skin surface | 2026-06-18 10:00:58 +02:00 | 0/0 | 1/61 | Publicada y sincronizada |
| `feat/endtop-sslg4` | `a0368c4` | `a0368c4` | fix(tests): replace escape-fraction guard with sipm-entry guard | 2026-06-18 11:04:36 +02:00 | 0/0 | 1/45 | Publicada y sincronizada |
| `feature/edge-scan-and-readout-grouping` | `031c131` | `031c131` | feat(analysis): emulate FastIC+ Sum-of-N channel grouping for timing resolution | 2026-04-30 15:49:40 -04:00 | 0/0 | 7/0 | Publicada, sincronizada y ya contenida en `main` |
| `feature/sipm-electronics-response` | `bb4cf6a` | `bb4cf6a` | Document and ignore generated simulation runs | 2026-06-09 14:35:34 +02:00 | 0/0 | 11/15 | Publicada y sincronizada, divergente de `main` |
| `fix/physics-baseline` | `26f4767` | `26f4767` | fix(optics): replace reflector volumes with bar skin surface | 2026-06-18 10:29:18 +02:00 | 0/0 | 1/9 | Publicada y sincronizada |
| `fix/readout-wrapping-config` | `c9e61af` | `c9e61af` | fix(optics): replace reflector volumes with bar skin surface | 2026-06-18 08:42:47 +02:00 | 0/0 | 1/17 | Publicada y sincronizada |
| `main` | `84e902c` | `84e902c` | fix(optics): replace reflector volumes with bar skin surface | 2026-06-18 10:29:27 +02:00 | 0/0 | 0/0 | Rama base publicada y sincronizada |

## Ultima rama activa

Si se interpreta "activa" como **la rama con el commit mas reciente**, la respuesta es:

1. **`exec27-edge-estimator-coupling`** y **`exec28-endonly-scan11-weighted`**
   - Commit: `bd24951`
   - Fecha: **2026-06-21 21:37:53 +02:00**
   - Estado GitHub: **no publicadas**
   - Observacion: ambas apuntan exactamente al mismo commit. `exec28` parece existir para separar un worktree o una linea de ejecucion posterior, pero todavia no tiene commit propio distinto de `exec27`.

Si se interpreta "activa" como **la rama mas reciente que existe en GitHub**, la respuesta es:

1. **`origin/feat/ej230-endonly-mylar`**
   - Commit: `e4110ce`
   - Fecha: **2026-06-18 19:39:31 +02:00**
   - Estado local: **no hay rama local equivalente**
   - Observacion: esta ref remota avanzo durante el `fetch`.

La rama checkout actualmente en este directorio es **`exec22-endtop-optfix`**, pero no es la mas reciente.

## Que contiene cada rama y por que existe

### `main`

Rama base publicada. Contiene el estado comun del repositorio a `84e902c`, con el cambio de optica que reemplaza volumenes reflectores por superficie de piel de barra. Existe como base de integracion estable para comparar el resto de trabajos.

### `fix/physics-baseline`

Rama de correccion fisica de base. Sus commits apuntan a valores de materiales EJ-204/EJ-200/EJ-230, rise time finito, seleccion/logging de centelleador activo, guardrail `physics_baseline_check` y actualizaciones de resultados/presentaciones. Existe para fijar parametros fisicos y pruebas antes de generar o defender resultados.

### `fix/readout-wrapping-config`

Rama para configurar lectura y wrapping. Introduce selector de configuracion `End/Top/EndTop`, logging de la configuracion activa y del wrapping por cara, guardrails de geometria, analisis de sigma vs x y documentacion de auditoria. Existe porque la geometria y el readout necesitaban poder alternarse y verificarse de forma reproducible.

### `diag/photon-budget`

Rama diagnostica del presupuesto de fotones. Agrega contadores de presupuesto optico, cierre de contabilidad terminal, comparaciones de timing/yield y un reporte diagnostico. Existe para explicar perdidas de fotones y comparar PE simulados con referencias/MPV externas.

### `backup/photon-budget-worktree`

Snapshot de seguridad del trabajo diagnostico de photon budget, con scripts/analisis y cambios de ejecucion preservados. Existe como respaldo de un worktree de diagnostico, no como linea principal de desarrollo.

### `feature/edge-scan-and-readout-grouping`

Rama historica de edge scan y agrupacion de lectura. Introduce analisis de escaneo cerca del borde, agrupacion tipo FastIC+ Sum-of-N y timing resolution. Esta rama ya esta contenida en `main` (`0` commits ahead), asi que su razon actual es principalmente historica o de trazabilidad.

### `feature/sipm-electronics-response`

Rama de respuesta electronica/SiPM y analisis asociado. Contiene trabajo de waveform/dCFD, displays de eventos, edge scans, configuracion de sensores y correcciones de analisis para posiciones negras o escasas. Existe para explorar el modelo de electronica del sensor y herramientas de visualizacion/analisis, pero diverge de `main`.

### `feat/endtop-sslg4`

Rama EndTop SSLG4. Agrupa analisis y presentaciones EXEC_12/EXEC_13 sobre observables EndTop, dispersion de tiempos `t_N`, perfiles Npe/timing y correcciones de entrega. Existe para estudiar la configuracion con lectura en extremos y top, especialmente para resultados EJ-204.

### `feat/ej204-event-display-tracks`

Rama de visualizacion de eventos EJ-204. Su commit distintivo agrega almacenamiento selectivo de trayectorias opticas para event display, encima de la linea EndTop. Existe para inspeccionar visualmente tracks/fotones sin guardar todo el trafico optico.

### `exp/pair-scan-2026-06-11`

Rama experimental de pair-scan. Incluye geometria, macros, runner batch, analisis ROOT, estudios de observabilidad `y=0`, combinacion de estimadores EndTop, reconstruccion global X y decks/reports reproducibles. Existe para estudiar pares de posiciones y estimadores temporales/espaciales de forma controlada.

### `feat/ej230-sslg4`

Rama EJ-230 SSLG4. Cambia el foco a EJ-230, incorpora reportes EXEC_14, scripts de ejecucion multi-host, runbook EJ230, fixes opticos y guardrails de tests. Existe para estudiar la variante EJ-230 con la configuracion SSLG4/EndTop.

### `feat/ej230-endonly-mylar`

Rama **solo remota** para EJ-230 End-only + Mylar. Elimina SiPMs TOP y ventanas, deja canales solo END, aplica Mylar en todas las caras salvo extremos, anade runbook/reportes y el diagnostico `PrintGroupVelDiagnostic()` para una anomalia de GROUPVEL. Existe para replicar la estrategia End-only+Mylar en EJ-230 y diagnosticar velocidad de propagacion.

### `feat/endonly-mylar`

Rama EJ-204 End-only + Mylar. Introduce hook `topSurface` (`mylar|sipm`), lectura robusta sin canales top, guardrails contra artefactos, pipeline de atenuacion/sigma_t, decks Beamer y auditoria fisica final. Existe para evaluar una geometria con lectura solo en extremos y Mylar en la superficie superior/lateral.

### `exec21-optfix`

Rama local de correccion optica sobre la linea End-only. Sus commits principales cambian la superficie a `dielectric_dielectric` para habilitar TIR y agregan `REFLECTIVITY=0.95` al modelo Mylar. Existe para probar una correccion fina del modelo optico antes de publicar.

### `exec22-endtop-optfix`

Rama local activa en este worktree. Porta la correccion TIR/skin surface de `exec21` a la construccion EndTop (`CreateBarSkinReflector`). Existe para aplicar el arreglo optico de End-only a la variante EndTop.

### `exec23-explicit-airgap`

Rama local que implementa geometria explicita de air gap para el reflector. Existe para pasar de una aproximacion de superficie a una representacion geometrica mas explicita del espacio aire/reflector.

### `exec24-pe-budget-audit`

Rama local de auditoria posterior al air gap. Agrega `EXEC_24_REPORT.md` y revisa presupuesto de PE y deteccion SiPM tras el fix de geometria. Existe para cuantificar si el air gap arreglo o desplazo el problema de coleccion optica.

### `exec25-optical-realism-bracket`

Rama local para acotar realismo optico. Valida semantica de PDE/superficie y explora el rango plausible de coleccion optica. Existe para separar error de modelado de incertidumbre fisica razonable.

### `exec26-scint-air-surface-realism`

Rama local para acotar la superficie centelleador-aire en coleccion END. Existe para estudiar cuanto depende la respuesta de los detalles de frontera entre scintillator y aire.

### `exec27-edge-estimator-coupling`

Rama local mas reciente. Diagnostica respuesta de borde mediante estimadores END y acoplamiento SiPM. Existe para entender degradaciones o sesgos cerca del borde y relacionarlos con el acoplamiento optico/electronico.

### `exec28-endonly-scan11-weighted`

Rama local que apunta al mismo commit que `exec27`. Por nombre, parece reservada para un escaneo End-only `scan11` ponderado, y esta checkout en `/mnt/d/SHiP/ej200_endonly`. Como no tiene commit propio, su contenido real aun es identico a `exec27`.

## Estado de sincronizacion con GitHub

Ramas locales que coinciden con GitHub despues del `fetch`:

- `backup/photon-budget-worktree`
- `diag/photon-budget`
- `exp/pair-scan-2026-06-11`
- `feat/ej204-event-display-tracks`
- `feat/ej230-sslg4`
- `feat/endonly-mylar`
- `feat/endtop-sslg4`
- `feature/edge-scan-and-readout-grouping`
- `feature/sipm-electronics-response`
- `fix/physics-baseline`
- `fix/readout-wrapping-config`
- `main`

Ramas locales que **no estan publicadas** en GitHub:

- `exec21-optfix`
- `exec22-endtop-optfix`
- `exec23-explicit-airgap`
- `exec24-pe-budget-audit`
- `exec25-optical-realism-bracket`
- `exec26-scint-air-surface-realism`
- `exec27-edge-estimator-coupling`
- `exec28-endonly-scan11-weighted`

Ramas que estan en GitHub pero **no existen localmente**:

- `feat/ej230-endonly-mylar`

## Observaciones operativas

- El worktree actual tiene cambios no trackeados previos al informe: `GROUP_VELOCITY_AUDIT.md` y `runs/`.
- Varias entradas de `git worktree list` bajo `/tmp/...` aparecen como `prunable` porque apuntan a rutas inexistentes. Eso no cambia los commits ni la sincronizacion con `origin`, pero conviene limpiar esos worktrees si molestan en operaciones futuras.
- `feat/ej204-event-display-tracks` tiene una rama remota con el mismo hash, pero no aparece con upstream configurado en `git branch -vv`; funcionalmente esta sincronizada, aunque se podria configurar tracking si se va a trabajar ahi.

## Comandos usados

```bash
git fetch --all --prune
git status --short --branch
git remote -v
git branch -a --verbose --no-abbrev
git branch -vv --no-abbrev
git for-each-ref refs/heads refs/remotes/origin --sort=-committerdate --format='...'
git rev-list --left-right --count <rama>...origin/<rama>
git rev-list --left-right --count main...<rama>
git log --all --graph --decorate --oneline --date=short -n 80
git ls-remote --heads origin
git worktree list --porcelain
```
