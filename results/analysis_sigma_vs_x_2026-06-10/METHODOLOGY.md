# Metodología espejo del análisis SUM4 de test beam

Referencia de solo lectura:
`/home/reriosto/SHiP/test-beam-ship-april-2026-master/process_sum4.C`.
El macro de Constanza no se ejecuta sobre los ROOT de simulación porque espera
el esquema de datos FastIC+ (`asic`, `ch`, `type`, timestamps y pulse width).

## Receta extraída de `process_sum4.C`

- La lectura SUM4 usa los canales FastIC+ 1 y 5
  (`process_sum4.C:33-35`). `build_master_fastic_events` construye eventos
  alrededor del trigger `ch==8,type==0`, conserva pulsos type 1 antes/después
  del trigger, exige un type 2 en el mismo canal y rechaza eventos sin canal
  válido (`process_sum4.C:245-250, 326-393, 414-485`).
- La distribución temporal se llena con diferencias entre el trigger y el
  tiempo type 1 del canal válido (`process_sum4.C:1152-1166`). Para el fit
  gaussiano se busca el bin de máximo dentro de la ventana de fit, se usa como
  semilla de la media, se inicializa `sigma=0.5 ns` y se imponen límites a
  amplitud, media y sigma (`process_sum4.C:1294-1321`). Los histogramas con
  corte repiten la siembra en el pico y el fit restringido
  (`process_sum4.C:1931-1978`).
- El FWHM se calcula tomando el primer y último bin por encima de la mitad del
  máximo y restando sus centros (`process_sum4.C:57-89`).
- El observable de energía es pulse width type 2. El fit Landau busca primero
  el máximo físico, construye un rango alrededor del pico, inicializa
  amplitud/MPV/sigma y limita los tres parámetros antes del fit
  (`process_sum4.C:553-563, 587-648`).
- Los cortes de energía se aplican exigiendo
  `pulse_width_lsb >= min_energy_lsb` al type 2; el pipeline produce cortes de
  400, 500 y 550 LSB (`process_sum4.C:1736-1753, 1869-1885, 2639-2676`).
- El mapa de correlación asocia pulse width type 2 con cada diferencia
  temporal válida (`process_sum4.C:2304-2332, 2356-2374, 2462-2488`).

## Adaptación operativa a la simulación

- Cada `event_id` ya identifica un evento Geant4; no se puede replicar el
  event builder FastIC+ porque los ROOT simulados no contienen pulsos
  type 0/1/2.
- Se forman SUM4 provisionales `{0..3}`, `{4..7}`, `{8..11}`, `{12..15}`. Para
  cada cluster se suma la respuesta 1-PE doble exponencial existente y se toma
  el primer cruce leading-edge a 4 PE. El tiempo de cada extremo es el primer
  SUM4 que cruza; `DeltaT_LR=t_L-t_R`.
- La cadena existente **incluye el walk inherente al leading-edge, pero no una
  corrección de time-walk**. Esto está declarado en
  `analysis/congruent_sum4_timing.C:1-10, 45-51, 170-213, 305-319`.
- Para espejar el fit temporal se busca el máximo del histograma y se usa como
  semilla. Como `DeltaT_LR` cambia con la posición, la ventana fija de 2-6 ns
  del dato se traslada a una ventana local de 4 ns centrada en el pico.
- NPE total por evento (`NPE_L+NPE_R`) sustituye pulse width como proxy de
  energía. Los cortes comparativos son los percentiles locales p25 y p50; no
  equivalen a 400/500/550 LSB sin calibración `ToT(NPE)`.
- Se reporta `sigma_single=sigma(DeltaT_LR)/sqrt(2)`. Es una resolución
  intrínseca sin SPTR ni jitter electrónico.

Implementación nueva: `analysis/tb_mirror_sigma_vs_x.C`. No modifica ni
ejecuta los macros de Constanza ni los análisis existentes del repositorio.
