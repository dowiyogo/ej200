# CODEX_AUDIT — Estado de desarrollo y coherencia física en todas las ramas

Fecha de auditoría: **2026-06-07**  
Repositorio: `/home/reriosto/SHiP/ej200_edge_scan` (ruta resuelta: `/mnt/d/SHiP/ej200_edge_scan`)  
Modo: **solo lectura**, salvo este reporte. No se cambió de rama ni se compiló.

## 1. Resumen ejecutivo

- Se auditaron **6 refs** (3 ramas locales y sus 3 refs remotas), correspondientes a **3 estados de código únicos**.
- `main` coincide con el commit de referencia **`4e46959`**; no hay stashes ni divergencia local/remota.
- Hallazgos independientes agrupados: **3 P0, 6 P1 y 2 P2**.
- Ninguna rama implementa la baseline temporal requerida para EJ-200: rise time `0.9*ns` y jitter `106 ps ⊕ 10 ps`.
- `main` y `feature/edge-scan-and-readout-grouping` omiten el rise time y usan jitter monolítico de `20 ps`.
- `feature/sipm-electronics-response` separa SPTR/electrónica, pero usa defaults incorrectos y no aplica el término electrónico al tiempo registrado.
- Las dos feature branches conservan geometría Mylar segmentada obsoleta; la rama electrónica además cambia silenciosamente a EJ-230 y Broadcom.
- Todas las ramas tienen más de 30 días sin commits a fecha de auditoría.
- **No está listo para generar datasets de presentación**: los P0 sesgan directamente la resolución temporal y no existe el guardrail `physics_baseline_check`.

## 2. Sanidad del entorno

- Rama activa: `main`, tracking `origin/main`, sin ahead/behind.
- Estado inicial: `?? audit/` (directorio sin trackear preexistente; no se inspeccionó ni modificó).
- Stashes: ninguno.
- `main`: `4e46959`, coincide con el commit de referencia solicitado.
- Remoto: `origin = git@github.com:dowiyogo/ej200.git`.
- No hay ramas remotas sin contraparte local ni ramas locales sin contraparte en `origin`.

## 3. Matriz de ramas

En “vs main”, `behind/ahead` representa commits exclusivos de `main` / exclusivos de la rama. “Edad” se calcula al 2026-06-07.

| Ref | Último commit | Fecha | vs main behind/ahead | Merge-base con main | Edad del tip | ¿en origin? | ¿stale >30 d? |
|---|---|---:|---:|---|---:|---|---|
| `main` | `4e46959` fix(optics): add reflector panels and photon-yield diagnostics | 2026-05-07 | 0 / 0 | `4e46959` (2026-05-07) | 31 d | Sí, idéntica | Sí |
| `feature/edge-scan-and-readout-grouping` | `031c131` feat(analysis): emulate FastIC+ Sum-of-N channel grouping for timing resolution | 2026-04-30 | 6 / 0 | `031c131` (2026-04-30) | 38 d | Sí, idéntica | Sí |
| `feature/sipm-electronics-response` | `d870974` fix(analysis): handle sparse black edge scans | 2026-05-06 | 10 / 11 | `4d7b9a8` (2026-04-09; 59 d) | 32 d | Sí, idéntica | Sí |
| `origin/main` | `4e46959` | 2026-05-07 | 0 / 0 | `4e46959` | 31 d | Remota | Sí |
| `origin/feature/edge-scan-and-readout-grouping` | `031c131` | 2026-04-30 | 6 / 0 | `031c131` | 38 d | Remota | Sí |
| `origin/feature/sipm-electronics-response` | `d870974` | 2026-05-06 | 10 / 11 | `4d7b9a8` | 32 d | Remota | Sí |

Conclusión de topología: la “Fase A: push de ambas feature branches” ya está satisfecha en el estado observado; ambas refs locales coinciden con `origin`.

## 4. Matriz de coherencia física

Las refs `origin/*` son idénticas a sus ramas locales y heredan exactamente los mismos valores.

| Parámetro | `main` | `feature/edge-scan-and-readout-grouping` | `feature/sipm-electronics-response` |
|---|---|---|---|
| Rise time EJ-200 | **Ausente** en `src/Materials.cc:100-105` — **P0** | **Ausente** en `src/Materials.cc:100-105` — **P0** | EJ-200 ausente; EJ-230 `0.5*ns` en `src/Materials.cc:184` y es default — **P0** |
| Decay/yield activo | `10000/MeV`, `2.1*ns`, `src/Materials.cc:102-105` — OK | `10000/MeV`, `2.1*ns`, `src/Materials.cc:102-105` — OK | Default EJ-230: `9700/MeV`, `1.5*ns`, `src/Materials.cc:182-186` — **P1** |
| Etiqueta/tabla PDE | S13360-6025, `src/SiPMSD.cc:17-31` — **P1** | S13360-6025, `src/SiPMSD.cc:17-31` — **P1** | S13360-6025 “legacy”, `src/SiPMSD.cc:105,251,320` — **P1** |
| Jitter temporal | Monolítico `20.0e-3 ns`, `include/SiPMSD.hh:43`; aplicado en `src/SiPMSD.cc:94-95` — **P0** | Igual a `main` — **P0** | Default Broadcom: SPTR `30 ps`, electrónica `30 ps`; para Hamamatsu SPTR `150 ps`; electrónica no se suma al tiempo, `src/SiPMSD.cc:167-168,198-199,308` — **P0** |
| Atenuación | `3.8*m`, `src/Materials.cc:77` — **P1** | `3.8*m`, `src/Materials.cc:77` — **P1** | EJ-200 `3.8*m`, `src/Materials.cc:97`; default EJ-230 `120*cm`, `src/Materials.cc:164` — **P1** |
| Reflectores | Paneles explícitos, `R=0.98`, `src/DetectorConstruction.cc:152-196`, `src/Materials.cc:261` — OK | Wrap Mylar segmentado, `src/DetectorConstruction.cc:175-240` — **P1** | Wrap Mylar segmentado, `src/DetectorConstruction.cc:214-280` — **P1** |
| Centelleador default | EJ-200, `src/DetectorConstruction.cc:121` — OK | EJ-200, `src/DetectorConstruction.cc:148` — OK | `fScintillatorName = "EJ230"`, `include/DetectorConstruction.hh:68,93` — **P1** |
| Sensor default | Hamamatsu S13360 implícito, aunque etiquetado 6025 — P1 por etiqueta | Igual a `main` — P1 por etiqueta | Broadcom AFBR-S4N66P024M, `include/SiPMSD.hh:59`; macros también lo fuerzan — **P1** |
| Barra | `2×(700,30,5) mm = 1400×60×10 mm`, `src/DetectorConstruction.cc:22-24` — OK | Igual, `src/DetectorConstruction.cc:23-25` — OK | Igual, `src/DetectorConstruction.cc:24-26` — OK |
| Mapeo SiPM | left `0-7`, right `8-15`, top `16-35` por default, `src/DetectorConstruction.cc:214,224,247`; `include/DetectorConstruction.hh:74` — OK | Igual, `src/DetectorConstruction.cc:258,268,291`; `include/DetectorConstruction.hh:86` — OK | Igual, `src/DetectorConstruction.cc:298,308,331`; `include/DetectorConstruction.hh:92` — OK |

## 5. Hallazgos detallados

### P0 — Bloqueantes para datasets temporales

1. **Rise time EJ-200 ausente o incorrecto**
   - `main` → `src/Materials.cc:100-105` → no existe `SCINTILLATIONRISETIME1` → esperado `0.9*ns` → **P0**.
   - `feature/edge-scan-and-readout-grouping` → `src/Materials.cc:100-105` → ausente → esperado `0.9*ns` → **P0**.
   - `feature/sipm-electronics-response` → `src/Materials.cc:120-125` → EJ-200 sin rise time; `src/Materials.cc:184` → EJ-230 `0.5*ns` activo por default → esperado EJ-200 `0.9*ns` → **P0**.
   - Corrección recomendada: definir `SCINTILLATIONRISETIME1 = 0.9*ns` para EJ-200 y habilitar finite rise time de forma verificable.

2. **Jitter monolítico de 20 ps en `main` y rama de agrupación**
   - Ambas ramas → `include/SiPMSD.hh:43` → `fJitterSigma = 20.0e-3` ns; `src/SiPMSD.cc:94-95` lo aplica como única gaussiana → esperado SPTR `106 ps σ` y FastIC+ `10 ps σ`, sumados en cuadratura → **P0**.
   - `/sipm/jitterSigma` existe (`src/SiPMSD.cc:47-54`), pero no está marcado como alias deprecado.
   - Corrección recomendada: separar ambos términos y conservar `/sipm/jitterSigma` únicamente como alias deprecado.

3. **Split incompleto e incorrecto en `feature/sipm-electronics-response`**
   - `include/SiPMSD.hh:59-61` → default Broadcom, SPTR `30 ps`, electrónica `30 ps`.
   - `src/SiPMSD.cc:308` → Hamamatsu recibe SPTR `150 ps`, no `106 ps`.
   - `src/SiPMSD.cc:167-168` → el tiempo registrado suma solo `sptrJitterNs`; `fElectronicsSigma` solo se escribe como metadata en `src/SiPMSD.cc:199`.
   - Esperado: S13360-6050PE default, SPTR `106 ps σ` y FastIC+ `10 ps σ`, ambos aplicados en cuadratura → **P0**.

### P1 — Incoherencias físicas y scope drift

1. **Longitud de atenuación incorrecta**
   - `main` y `feature/edge-scan-and-readout-grouping` → `src/Materials.cc:77` → `3.8*m` → esperado `160*cm` → **P1**.
   - `feature/sipm-electronics-response` → EJ-200 `3.8*m` en `src/Materials.cc:97`; default EJ-230 `120*cm` en `src/Materials.cc:164` → esperado EJ-200 `160*cm` → **P1**.

2. **PDE identificada como S13360-6025, no S13360-6050PE**
   - `main` y rama de agrupación → `src/SiPMSD.cc:17`, `include/SiPMSD.hh:12`, `src/Materials.cc:153/149`.
   - Rama electrónica → `src/SiPMSD.cc:105,251,320`, `macros/run.mac:25`, `macros/scan.mac:25`.
   - Valor actual: tabla/identificador 6025 o “legacy”; esperado: comentario/identificador y valores validados para **S13360-6050PE** → **P1**.

3. **Geometría reflectora obsoleta en ambas feature branches**
   - Rama de agrupación → `src/DetectorConstruction.cc:175-240` → envoltura Mylar segmentada.
   - Rama electrónica → `src/DetectorConstruction.cc:214-280` → envoltura Mylar segmentada.
   - Esperado: paneles reflectores explícitos de `main` con `R=0.98` (`main:src/DetectorConstruction.cc:152-196`, `main:src/Materials.cc:261`) → **P1**.

4. **Default EJ-230 y propiedades temporales/yield distintas en rama electrónica**
   - `include/DetectorConstruction.hh:68,93` → default `"EJ230"`.
   - `src/Materials.cc:182-185` → `9700/MeV`, rise `0.5*ns`, decay `1.5*ns`.
   - `macros/run.mac:16` y `macros/scan.mac:16` fuerzan EJ-230.
   - Esperado: EJ-200, `10000/MeV`, rise `0.9*ns`, decay `2.1*ns` → **P1**.

5. **Default Broadcom en rama electrónica**
   - `include/SiPMSD.hh:59` → `kBroadcomAFBRS4N66P024M`.
   - `macros/run.mac:26` y `macros/scan.mac:26` fuerzan `/sipm/sensorModel Broadcom`.
   - Esperado: Hamamatsu S13360-6050PE → **P1**.

6. **Macros de edge-wrap engañosos en `main`**
   - `macros/scan_edge_air.mac:11` y `macros/scan_edge_black.mac:11` invocan `/det/edgeWrap`.
   - `include/DetectorConstruction.hh:32` y `src/DetectorConstruction.cc:90-112` documentan/implementan el comando como no-op tras eliminar WrapLV.
   - Riesgo: datasets nominalmente “air” y “black” ejecutan la misma geometría reflectora de `main` → **P1**.

### P2 — Guardrails/documentación

1. **No existe `tests/physics_baseline_check.cc` ni un test equivalente**
   - Ninguna rama contiene `tests/`.
   - Solo existen `smoke_test` y `edge_scan_smoke` en `CMakeLists.txt:20,25` (`24,29` en la rama electrónica) → **P2**.

2. **Comentario de reflectividad inconsistente en `main`**
   - `include/Materials.hh:22` afirma `R(λ) = 0.99`, mientras `src/Materials.cc:261` implementa `0.98`.
   - La implementación coincide con el valor esperado; corregir documentación → **P2**.

## 6. Divergencias entre ramas

- `feature/edge-scan-and-readout-grouping` está totalmente contenida en `main`: está **6 commits behind y 0 ahead**. Su física está stale respecto a los paneles explícitos y diagnósticos de yield incorporados en `main`.
- `feature/sipm-electronics-response` divergió desde `4d7b9a8`: está **10 behind y 11 ahead**. Contiene trabajo útil de análisis/electrónica, pero mezcla scope drift físico: EJ-230, Broadcom, múltiples sensores, wrap Mylar segmentado y macros de visualización.
- `main` tiene el modelo reflector esperado, pero aún conserva baseline temporal y atenuación incorrectas.
- No hay divergencia local/remota: ambas feature branches ya están publicadas en `origin`.

## 7. Geometría, macros, mensajería y tests

- `main`: 22 macros; añade `diag_optical.mac` y `diag_scan5.mac`. `/det/edgeWrap` es un no-op legado.
- `feature/edge-scan-and-readout-grouping`: 20 macros; elimina ambos macros de diagnóstico de `main`.
- `feature/sipm-electronics-response`: 23 macros; añade `vis_scan_png*.mac`, fuerza EJ-230/Broadcom en `run.mac` y `scan.mac`, y cambia `scan_step.mac` de 200 a 20000 eventos.
- `/sipm/jitterSigma`: existe solo en `main` y rama de agrupación; **no** está marcado como alias deprecado.
- Rama electrónica: reemplaza el comando por `/sipm/sptrSigma` y `/sipm/electronicsSigma`; no ofrece alias `/sipm/jitterSigma`.
- `CMakeLists.txt`: rama de agrupación idéntica a `main`; rama electrónica añade copia post-build de macros y cambia regex de los smoke tests.
- Ninguna rama contiene `tests/physics_baseline_check.cc`.
- No se compiló ni ejecutó CTest, conforme a la restricción de Fase 4.

Macros por estado único:

- `main`: `crosscheck_thickness.mac`, `diag_optical.mac`, `diag_scan5.mac`, `event_display_top_3d.mac`, `event_display_top_batch.mac`, `event_display_top_lateral.mac`, `event_display_top_midpoint.mac`, `print_optical_processes.mac`, `run.mac`, `run_1000.mac`, `scan.mac`, `scan_angle.mac`, `scan_angle_step.mac`, `scan_edge.mac`, `scan_edge_air.mac`, `scan_edge_black.mac`, `scan_edge_neg.mac`, `scan_edge_step.mac`, `scan_step.mac`, `test.mac`, `vis.mac`, `vis_photon_traces_xz.mac`.
- Rama de agrupación: lista de `main` sin `diag_optical.mac` ni `diag_scan5.mac`.
- Rama electrónica: lista de rama de agrupación más `vis_scan_png.mac`, `vis_scan_png_lite.mac`, `vis_scan_png_step.mac`.

## 8. Yields de referencia y reproducibilidad

| Magnitud de referencia | Objetivo | Estado de reproducibilidad |
|---|---:|---|
| End-left | ~28 ph/evt | No verificado; no hay resultado versionado y no se ejecutó simulación |
| End-right | ~28 ph/evt | No verificado; no hay resultado versionado y no se ejecutó simulación |
| Top | ~334 ph/evt | No verificado; no hay resultado versionado y no se ejecutó simulación |
| Eficiencia de detección | ~2.06% | No verificado; `main` puede imprimirla, pero no hay baseline test/resultados versionados |

`main` conserva `diag/run_audit.sh`, `diag/yield_audit.C`, `diag_optical.mac` y `diag_scan5.mac`, por lo que ofrece la mejor base para reproducir los yields tras corregir la baseline física. Ambas feature branches eliminan esos diagnósticos. Además, `main:src/Materials.cc:247-249` documenta que el reflector explícito `R=0.98` se eligió para mantener el yield de ends por encima del umbral, pero no contiene los valores de referencia solicitados como artefacto verificable.

## 9. Próximos pasos recomendados

### Fase A — `CODEX_EXEC_01`

- Confirmar que no hace falta push: las dos feature branches ya coinciden con sus refs en `origin`.
- Crear, tras aprobación, el tag de seguridad `checkpoint/pre-physics-baseline-2026-05-08` apuntando a `main` (`4e46959`) y publicarlo.
- No modificar código en esta fase.

### Fase B — `CODEX_EXEC_02`

Crear `fix/physics-baseline` desde `main` con cuatro commits separados:

1. Añadir rise time EJ-200 `SCINTILLATIONRISETIME1 = 0.9*ns` y habilitar finite rise time.
2. Corregir/validar la tabla y comentarios PDE para S13360-6050PE.
3. Implementar SPTR `106 ps σ` y FastIC+ `10 ps σ`, aplicados en cuadratura; mantener `/sipm/jitterSigma` como alias deprecado.
4. Añadir `tests/physics_baseline_check.cc` y registrarlo en CTest, incluyendo checks de rise/decay/yield/attenuation, defaults de sensor y jitter.

Incluir también la corrección de `ABSLENGTH` a `160*cm`, resolver los macros edge-wrap no-op y corregir el comentario `R=0.99`.

### Fase C — `CODEX_EXEC_03`

Pendiente de especificar. Recomendación: definirla solo después de que Fase B pase CTest y se reproduzcan los yields de referencia con semillas/configuración versionadas.

## 10. Anexo — Comandos usados para reproducibilidad

Todos se ejecutaron desde `/home/reriosto/SHiP/ej200_edge_scan`; no se usó `checkout`, `switch`, compilación ni comandos destructivos.

```bash
pwd
git status --porcelain=v1 --branch
git stash list
git rev-parse --short main
git remote -v
git --no-pager branch -avv
git --no-pager log --graph --decorate --oneline --all --date-order -30

git --no-pager log -1 --format='%h|%ci|%an|%s' <ref>
git rev-list --left-right --count main...<ref>
git merge-base main <ref>
git --no-pager log -1 --format='%h|%ci|%an|%s' <merge-base>

git --no-pager diff --stat main <rama>
git --no-pager diff --name-status main <rama>
git --no-pager diff main <rama> -- CMakeLists.txt

git --no-pager grep -nEi \
  'SCINTILLATIONRISETIME1|RISETIME|RISE_TIME|SCINTILLATIONYIELD|SCINTILLATIONTIMECONSTANT|DECAYTIME|FASTTIMECONSTANT|ABSLENGTH|attenuation' \
  <rama> -- 'src/*' 'include/*' '*.mac' CMakeLists.txt

git --no-pager grep -nEi \
  'S13360|6025|6050|PDE|jitterSigma|jitter|SPTR|FastIC|106|10\*ps|20\*ps' \
  <rama> -- 'src/*' 'include/*' '*.mac' CMakeLists.txt

git --no-pager grep -nEi \
  'Reflectivity|0\.98|Reflector|WrapLV|Mylar|G4LogicalSkinSurface|G4LogicalBorderSurface|EJ-?200|EJ-?230|Broadcom|AFBR|1400|barLength|barWidth|barThickness|copyNo|GetCopyNumber|top|end' \
  <rama> -- 'src/*' 'include/*' '*.mac' CMakeLists.txt

git --no-pager grep -nE '^[/#]*\s*/(sipm|det|run)/' <rama> -- '*.mac' 'src/*' 'include/*'
git --no-pager grep -nE 'add_test|CTest|physics_baseline_check' <rama> -- CMakeLists.txt 'tests/*'
git ls-tree -r --name-only <rama>
git show <rama>:<ruta>
```
