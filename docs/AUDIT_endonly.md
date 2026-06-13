# Auditoría previa: rama End-only + Mylar

Fecha: 2026-06-13  
Rama auditada: `feat/endtop-sslg4`  
HEAD base real: `c6e7843b2486fb3d073f4d56f7a1f0767e14419e`

El HEAD esperado en el encargo era `0201413`, pero la rama local real estaba cinco
commits por delante de `origin/feat/endtop-sslg4` (`5783a0d`). El remoto no estaba
por delante, por lo que no se hizo `pull`. Esta auditoría describe el código real
en `c6e7843`; no modifica ni corrige parámetros fuera del alcance.

## Hallazgo de coordenadas que condiciona la implementación

El encargo denomina extremos a `±Z` y pide envolver `±X/+Y/-Y`. En el código
base, sin embargo, la dimensión de 1.4 m y el barrido longitudinal son **X**:
la barra tiene semilongitudes `(700, 30, 5) mm`
(`src/DetectorConstruction.cc:28-31`), `/muon/gunX` mueve el haz sobre X
(`src/PrimaryGeneratorAction.cc:51-61,91-95`) y los END se colocan en `-X/+X`
(`src/DetectorConstruction.cc:329-343`). La incidencia vertical es sobre Z
(`src/PrimaryGeneratorAction.cc:20-23,121-131`).

Por tanto, preservar los END byte-idénticos y envolver a la vez `±X` es
geométricamente incompatible. La implementación End-only preservará los END
reales en `±X` y uniformará con Mylar las caras no-extremo reales `±Y/±Z`.
Este es el único modo de satisfacer la restricción prioritaria de no tocar los
END ni su acoplamiento.

## Geometría y sensores

La jerarquía declarada contiene 8 SiPM END izquierdos, 8 derechos y 70 TOP
(`include/DetectorConstruction.hh:16-25,37-39`). Cada volumen END representa un
activo individual de `6 x 6 mm²`, con espesor simulado de `0.5 mm`, pitch END
de `7.5 mm` y centros Y `(i-3.5)*pitch`
(`src/DetectorConstruction.cc:33-37,323-343`). Cada volumen físico tiene un
único `copyNo`, que es también el canal global: END izquierdo `0..7`, END
derecho `8..15` (`src/DetectorConstruction.cc:51-60,329-343`). En consecuencia,
los 16 canales END corresponden a 16 SiPM físicos single-element, uno por canal.

La agrupación analógica existente usa cuatro SUM4:

- END izquierdo: A=`{0,1,2,3}`, B=`{4,5,6,7}`.
- END derecho: A=`{8,9,10,11}`, B=`{12,13,14,15}`.

Hay 4 SiPM por suma y 2 sumas por extremo
(`analysis/exec07/common.py:12-22`,
`analysis/congruent_sum4_timing.C:216-225`). El tiempo de cada extremo es el
primer SUM4 que cruza umbral (`analysis/congruent_sum4_timing.C:228-233`).

Los TOP son 70 activos individuales `6 x 6 mm²`, IDs globales `16..85`, sobre
la cara `+Y`; sus centros recorren X con un hueco central de 24 mm
(`src/DetectorConstruction.cc:39-48,347-364`,
`tests/readout_config_check.cc:71-97`).

## Verificación AFBR-S4N66P014M solicitada

El código base **no nombra P014M**. Nombra y solo acepta
`AFBR-S4N66P024M` (`include/DetectorConstruction.hh:87-90`,
`src/SiPMModel.cc:28-44`). No se cambia este hecho porque el encargo exige
documentar discrepancias sin corregir el SiPM.

| Parámetro solicitado | Estado real en la base |
|---|---|
| Single-element, 1 canal | Coincide funcionalmente: cada placement tiene un único `copyNo`/canal (`src/DetectorConstruction.cc:329-343`; `src/SiPMSD.cc:55-68`). |
| Área sensible 6 x 6 mm² | Coincide para END y TOP (`src/DetectorConstruction.cc:33-42`). |
| Package P014M 6.5 x 7.1 mm² | No se modela; solo se modela el activo/acoplamiento 6 x 6 x 0.5 mm³. |
| Ventana epoxy EMC UV | Discrepa: el volumen usa `G4_SILICON_DIOXIDE` con `RINDEX=1.58`, descrito como acoplamiento perfecto (`src/Materials.cc:186-200`). No hay grosor de epoxy DS104. |
| Pitch microcelda 40 µm / ~22428 celdas | No se modela ni aparece en el código. El pitch geométrico END de placements es 7.5 mm, concepto distinto (`src/DetectorConstruction.cc:37`). |
| Rango espectral 250-900 nm | La tabla PDE sí cubre 250-900 nm (`data/sipm/AFBR-S4N66P024M_pde.txt:4,7-42`). El stepping action mata fotones fuera de 300-900 nm, por lo que 250-299 nm no llega al sensor (`src/SteppingAction.cc:97-110`). |
| Pico 420 nm / PDE 63% | Coincide: tabla PDE 0.63 a 420 nm (`data/sipm/AFBR-S4N66P024M_pde.txt:17-18`) y guardrail explícito (`tests/readout_config_check.cc:41-57`). La tabla también da 0.63 a 400 nm. |
| Crosstalk 23% / afterpulsing <1% | Solo documentados; no se simulan (`data/sipm/AFBR-S4N66P024M_pde.txt:5`). |
| DCR, gain, breakdown, overvoltage | No se modelan. |

La superficie de detección es `dielectric_metal`, `unified`, `polished`,
`sigmaAlpha=0`, con la PDE como `EFFICIENCY` y reflectividad cero
(`src/Materials.cc:222-249`). El SD registra únicamente fotones aceptados por el
proceso de frontera (`include/SiPMSD.hh:12-25`).

## Tratamiento óptico por cara en EndTop

No hay `G4LogicalSkinSurface`. Toda la detección y envoltura activa se hace con
`G4LogicalBorderSurface` desde `BarPV` hacia SiPMs o paneles hermanos
(`include/DetectorConstruction.hh:33,47-51`;
`src/DetectorConstruction.cc:254-287,335-343,362-363`).

La superficie común de los paneles es `BarSkinReflector`:
`dielectric_metal`, modelo `unified`, finish `polished`, `sigmaAlpha=0`, y
`REFLECTIVITY=0.98` uniforme entre 350 y 800 nm
(`src/Materials.cc:293-332`).

| Cara real | Estado en configuración `EndTop` |
|---|---|
| `-X` | 8 END activos. No hay panel reflector X cuando END está instrumentado (`src/DetectorConstruction.cc:309-320,322-345`). Las zonas de la cara no cubiertas por los activos no tienen superficie explícita. |
| `+X` | Igual que `-X`, con IDs `8..15` (`src/DetectorConstruction.cc:338-343`). |
| `-Y` | Panel sólido con `BarPV -> ReflectorYMinusPV` (`src/DetectorConstruction.cc:258-274,282-283`). |
| `+Y` | Panel reflector perforado mediante 70 sustracciones exactas y 70 TOP en las ventanas (`src/DetectorConstruction.cc:289-308,347-364`). **No está abierta salvo en las ventanas TOP**; el resto está cerrado por reflector. |
| `-Z` | Panel reflector sólido (`src/DetectorConstruction.cc:260-280,284-285`). |
| `+Z` | Panel reflector sólido (`src/DetectorConstruction.cc:260-280,286-287`). |

El cambio End-only uniformará `-Y/+Y/-Z/+Z` bajo una única superficie Mylar
configurable. En la base ya comparten una superficie común R=0.98; el cambio
del default a R=0.90 y de rugosidad será intencional y macro-configurable. Se
mantendrá el patrón `G4LogicalBorderSurface` con paneles explícitos y el modelo
`unified/dielectric_metal`, porque evita que una skin surface tape las caras
activas y es consistente con la geometría base
(`src/DetectorConstruction.cc:254-256`; `src/Materials.cc:293-315`).

## Sensitive detector y mapa de canales

`ConstructSDandField()` crea o reutiliza un solo `SiPMSD`, le asigna el modelo
PDE y lo vincula a los logical volumes END/TOP que existan
(`src/DetectorConstruction.cc:370-383`). La enumeración se deriva del `copyNo`
y `FaceType`: `0..7 -> face_type 0`, `8..15 -> face_type 1`, el resto
`-> face_type 2` (`src/DetectorConstruction.cc:51-60`;
`src/SiPMSD.cc:55-68,87-110`). EndTop tiene 86 canales; End tiene 16
(`tests/readout_config_check.cc:71-74,100-114,138-152`).

El TTree `sipm_hits` conserva un esquema común por hit con `event_id`,
`face_type`, `global_id`, `local_id`, `time_ns`, energía, longitud de onda,
PDE, posición y `gun_x_mm` (`src/RunAction.cc:115-133`).

## Barrido longitudinal y threads

El barrido de producción EndTop existente usa 31 posiciones
`-690,-670,-650..650(step 50),670,690`, 2000 eventos por posición por default,
y fija X con `/muon/gunX` (`scripts/run_exec07_scan.sh:9-11,35-39,63-68`).
La macro genérica `scan.mac` es distinta: `-650..650` con paso 65 y 200 eventos
(`macros/scan.mac:26-27`; `macros/scan_step.mac:4-5`).

`main()` crea `G4RunManagerType::Default` en
`main.cc:18-19`. No hay parser CLI ni variable de entorno propios para threads.
El script EndTop escribe `/run/numberOfThreads $workers` antes de inicializar,
con `WORKERS=16` por default (`scripts/run_exec07_scan.sh:10,52-62`).

## Pipeline END/SUM4 y paridad anti-artefacto

El pipeline congruente usa:

- pulso SPE doble exponencial normalizado `0.5/5.0 ns`
  (`analysis/congruent_sum4_timing.C:45-49,155-167`);
- umbral leading-edge absoluto `4 PE` y algoritmo de cruce por bisección
  (`analysis/congruent_sum4_timing.C:48-49,170-213`);
- cuatro grupos END SUM4 y primer grupo disparado por extremo
  (`analysis/congruent_sum4_timing.C:216-233`);
- estimador `sigma(deltaT)/sqrt(2)`
  (`analysis/congruent_sum4_timing.C:7-10,243-254`);
- `time_ns = Geant4 global time + jitter`; la campaña fija jitter a cero
  (`src/SiPMSD.cc:75-80`; `scripts/run_exec07_scan.sh:61-68`);
- t=0 común definido por el inicio del evento/particle gun
  (`src/PrimaryGeneratorAction.cc:18-23`; `src/EventAction.cc:15-20`).

La auditoría anti-artefacto existente ya concluye que umbral, pulso, t=0,
estimador y corrección de walk son comunes
(`audit/exec09_timing_mechanism.md:33-49`). Estos valores se conservarán y se
convertirán en guardrail automatizado.

`analysis/congruent_sum4_timing.C` tolera de hecho grupos TOP vacíos porque solo
añade `topDelta` si ambos disparan (`analysis/congruent_sum4_timing.C:339-375`),
pero aún crea columnas/gráficas TOP. En cambio,
`analysis/exec07_photon_budget.py` exige exactamente los 86 canales
(`analysis/exec07_photon_budget.py:89-123`) y usa matrices/perfiles TOP de forma
incondicional (`analysis/exec07_photon_budget.py:335-383`). Ese pipeline no es
robusto a END-only y debe adaptarse o complementarse.

## Guardrails existentes antes del cambio

CTest registra smoke general, smoke de edge, propiedades SSLG4, configuraciones
de readout, export/check GDML EndTop y presupuesto/balance de fotones EndTop
(`CMakeLists.txt:28-65`).

- `sslg4_properties_check` fija para OPSC-101/EJ-204: yield `10400/MeV`,
  rise `0.7 ns`, decay `1.8 ns`, ABSLENGTH `160 cm`, RINDEX `1.58` y pico
  `408.8 nm` (`tests/physics_baseline_check.cc:42-85`). La fuente SSLG4 declara
  los mismos yield/times (`src/external/SSLG4/macros/oscnt/opsc-101.mac:5-13`)
  y datos de ABSLENGTH/RINDEX
  (`src/external/SSLG4/data/oscnt/opsc-101/absLength.txt:1-2`,
  `src/external/SSLG4/data/oscnt/opsc-101/rIndex.txt:1-2`). No se modificarán.
- `readout_config_check` verifica conteos, placements, ventanas y PDE
  (`tests/readout_config_check.cc:100-160`).
- `endtop_balance_smoke` exige luz generada/detectada, fuga menor a 35% y más
  cruces de panel que escapes (`tests/check_endtop_balance.py:19-49`).

No existe antes del cambio un CTest que afirme explícitamente: cero TOP en
modo Mylar, +Y cerrada sin detección TOP, orden resolución intrínseca frente a
total medida, o paridad anti-artefacto. La referencia de 88 ps es solo una línea
cualitativa en la gráfica (`analysis/congruent_sum4_timing.C:267-275`), no un
guardrail automatizado.

## Problemas fuera de alcance documentados, no corregidos

- Nombre/modelo P024M en vez de P014M, ventana/acoplamiento simplificado y
  ausencia de package/microceldas/saturación.
- Filtro de tracking 300-900 nm frente a tabla PDE 250-900 nm.
- Las regiones no cubiertas de las caras END `±X` carecen de panel explícito.
- Comentarios antiguos alternan EJ-200/EJ-204 y algunos describen Mylar físico,
  aunque la geometría efectiva usa paneles de aire con superficie óptica.
- La base local contiene cinco commits no publicados sobre
  `origin/feat/endtop-sslg4`; la nueva rama parte deliberadamente del HEAD local
  solicitado como estado real.

## Addendum posterior a la implementación

El modo nuevo y sus defaults se exponen en
`include/DetectorConstruction.hh:32-35,80-87,104-107` y se registran como
comandos UI en `src/DetectorConstruction.cc:98-123`. `IsTopInstrumented()` exige
ahora explícitamente `topSurface=sipm`
(`src/DetectorConstruction.cc:198-200`), por lo que `mylar` no crea ventanas ni
placements TOP (`src/DetectorConstruction.cc:361-380,419-437`) y el SD solo se
asocia al logical volume END existente (`src/DetectorConstruction.cc:443-455`).

Las caras no-END reales comparten una sola instancia de reflector seleccionada
en `src/DetectorConstruction.cc:323-328` y las cuatro border surfaces se colocan
en `src/DetectorConstruction.cc:344-380`. La superficie Mylar implementada es
`unified/dielectric_metal/ground`, con reflectividad, lobe y sigmaAlpha
configurables; spike y backscatter son cero
(`src/Materials.cc:337-357`). El bloque END quedó sin cambios funcionales
(`src/DetectorConstruction.cc:394-417`) y el guardrail compara ese bloque byte a
byte contra `feat/endtop-sslg4`.

Los threads se resuelven con prioridad CLI > entorno > default y se aplican al
run manager en `main.cc:33-76`. El barrido portable, directorios timestamp,
metadata, macro Mylar, validación ROOT y análisis automático están en
`scripts/run_scan.sh:8-16,47-69,90-105,108-162`.

El análisis END-only conserva pulso `0.5/5 ns`, umbral `4 PE`, primer SUM4 por
extremo y `sigma(deltaT_LR)/sqrt(2)`
(`analysis/endonly_sum4.py:16-21,45-119,202-250`). Emite la curva y fits de
atenuación y la curva temporal en
`analysis/endonly_sum4.py:317-400`; TOP se filtra condicionalmente y se cuenta
como diagnóstico (`analysis/endonly_sum4.py:202-206,233-250,331-365`).

Los cuatro CTest nuevos se registran en `CMakeLists.txt:67-84`: mapa sin TOP,
presupuesto de fotones sin detección TOP, orden físico y paridad anti-artefacto.
El guardrail de orden físico usa el único término total definido en el código:
la suma en cuadratura de la resolución intrínseca y el jitter default de lectura
de 20 ps (`analysis/endonly_sum4.py:20,110-119,345-362`). No compara contra los
88 ps históricos, porque la propia base los identifica como referencia
cualitativa de otro detector (`analysis/congruent_sum4_timing.C:267-275`).

Smoke local ejecutado tras los guardrails: una posición (`x=0`), 200 eventos,
16 threads, ROOT válido, 31 488 hits END, cero hits TOP, 100% de eventos con
ambos extremos disparados y CSVs producidos. El ajuste corto dio
`sigma_single = 128.4 ± 17.9 ps`; este valor no se hardcodea ni se usa como
guardrail.
