# Diagnóstico del bajo rendimiento de luz End

Fecha: 2026-06-09  
Rama: `diag/photon-budget`  
Configuración de las corridas chicas: EJ-204, mu- 1 GeV, rise time finito,
jitter 0, 1000 eventos por punto, `x = 0, 350, 650 mm`.

## Parte A: auditoría óptica read-only

- Barra: `1400 x 60 x 10 mm`; SiPM End: `6 x 6 mm` y 8 por extremo
  (`src/DetectorConstruction.cc:24-32`). La cobertura geométrica es
  `8*36/600 = 48%`.
- Hay paneles reflectores y `G4LogicalBorderSurface` en `-Y`, `-X`, `+X`,
  `-Z` y `+Z` (`src/DetectorConstruction.cc:207-232`). No hay panel `+Y`:
  el código declara explícitamente esa cara abierta para Top
  (`src/DetectorConstruction.cc:188-190`). Esto es una fuga obvia para una
  comparación con el test beam End-only.
- Las zonas no cubiertas de las caras End encuentran los paneles completos
  `+/-X`; no son huecos absorbentes. Los SiPM End están desnudos, embebidos y
  flush con la barra, sin guía de luz (`src/DetectorConstruction.cc:234-267`).
- El reflector es `dielectric_metal`, modelo `unified`, acabado `polished`,
  `R=0.98` constante (`src/Materials.cc:324-344`).
- La barra y el volumen de acoplamiento SiPM tienen `RINDEX=1.58`; no hay gap
  de aire ni TIR en barra-SiPM (`src/Materials.cc:37-52,185-199`).
- La superficie SiPM es `dielectric_dielectric/unified/polished`, sin
  `EFFICIENCY` ni `REFLECTIVITY`; transmite al volumen y `SiPMSD` aplica
  manualmente una PDE Bernoulli, después mata el fotón
  (`src/Materials.cc:221-262`, `src/SiPMSD.cc:100-139`).
- La PDE usada está rotulada Hamamatsu S13360-6025 y alcanza 0.405 a 460 nm
  (`src/SiPMSD.cc:17-33`). La identidad S13360 frente a S14160 sigue abierta.

**Veredicto A:** no hay gap de acoplamiento ni fuga en la fracción no cubierta
de End, pero la cara `+Y` casi completa está abierta. La geometría actual
incluye Top y no representa un setup End-only completamente envuelto.

## Parte B: balance de fotones

Los destinos terminales cierran al 100%. `Direct Bar -> World` cuenta solo la
salida directa desde la barra; las pérdidas que atraviesan el volumen-panel ya
quedan clasificadas como absorción superficial.

| x [mm] | llegan por primera vez a cara End | entran End SiPM | PE End | PE Top | bulk abs. | surface abs. | Bar -> World directo |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.056% | 0.557% | 0.195% | 1.479% | 9.309% | 17.021% | 68.843% |
| 350 | 2.385% | 1.282% | 0.449% | 1.474% | 9.189% | 16.891% | 68.393% |
| 650 | 21.915% | 13.040% | 4.517% | 1.394% | 8.071% | 15.050% | 59.814% |

En el scan existente de 10000 eventos, el centro produce `39.549 PE/evento`
en End (`19.719` izquierda, `19.830` derecha), equivalente a `0.1999%` de los
fotones de centelleo. En `x=-690/+690 mm`, End total sube a
`1926.37/1915.71 PE/evento`, pero el extremo lejano queda en solo
`4.97/4.95 PE/evento`.

**Sumidero dominante:** la salida directa por la cara abierta domina por mucho
al bulk de EJ-204. La atenuación física existe, pero no explica por sí sola la
hambruna central. La absorción superficial acumulada por muchos rebotes es el
segundo sumidero.

Tablas generadas:

- `results/diag_2026-06-09/photon_budget.csv`
- `results/diag_2026-06-09/cheap_efficiency_31_positions.csv`

## Parte C: HOOK_ENDRED

Resultado en el centro, con el mismo SPR y umbral leading-edge provisionales:

| Reducción por extremo | N_eff | sigma_single End intrínseco |
|---|---:|---:|
| primer cluster SUM4 que cruza | 9732 | 467.46 +/- 5.80 ps |
| un pulso SUM8 por extremo | 9986 | 325.05 +/- 3.91 ps |
| promedio de dos tiempos SUM4 | 6825 | 522.44 +/- 8.20 ps |

SUM8 reduce sigma central en 30.5% respecto al hook actual, pero todavía queda
en `3.69 x 88 ps`. El promedio de clusters no ayuda con el umbral provisional.
En los extremos de la barra, SUM8 da `731 ps` a -690 mm y `689 ps` a +690 mm:
la gran luz del extremo cercano no compensa que el extremo lejano tenga menos
de 5 PE/evento. Para `DeltaT=t_L-t_R`, ambos extremos deben disparar.

`analysis/endred_photon_budget.C` acepta `HOOK_SUM_WIDTH=4|8`; la comparación
final con 88 ps requiere conocer la posición y el ancho SUM de la medición de
Gerardo. Figuras y tabla:

- `results/diag_2026-06-09/endred_summary.csv`
- `results/diag_2026-06-09/endred_comparison_vs_position.pdf`
- `results/diag_2026-06-09/best_endred_vs_light.pdf`

## Parte D: Landau MPV

No se encontraron las slides 12/22 de Constanza en el repositorio, bajo
`/home/reriosto/SHiP`, ni mediante la búsqueda disponible en Dropbox. Por eso
no es posible realizar una comparación numérica honesta contra el MPV.

El hook pendiente es:

`MPV_LSB_predicho = N_PE_sim * ganancia_SiPM * e / Q_LSB_energía`

Faltan la ganancia usada, el LSB de carga/energía y el MPV por ASIC/canal. El
LSB temporal de aproximadamente 24 ps no sirve para esta conversión de carga.

## Veredicto integrado

1. **Artefacto/sub-recolección del modelo actual: dominante.** La cara `+Y`
   abierta produce una fuga directa de 68.8% en el centro. Es especialmente
   incongruente con el test beam End-only.
2. **Hambruna física real: secundaria pero relevante.** El bulk de EJ-204
   elimina 9.3% en el centro y el extremo lejano cae a unos 5 PE cerca de
   `|x|=690 mm`.
3. **Hook ENDRED: contribuye, pero no resuelve.** SUM8 baja 467 ps a 325 ps;
   queda muy por encima de 88 ps incluso antes de jitter.

Orden recomendado: primero reproducir la envoltura y geometría End-only real
(en particular decidir/cerrar `+Y`) y validar la carga con el MPV de Constanza;
después fijar posición/ancho SUM y hooks de electrónica. Ajustar solo ENDRED no
puede corregir una pérdida óptica dominante.
