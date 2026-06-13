# Rama `feat/endonly-mylar`

Esta rama implementa el caso controlado **END-only + Mylar** para EJ-204:
conserva los 16 SiPM END y su acoplamiento de `feat/endtop-sslg4`, elimina los
70 SiPM TOP en el modo por defecto y cierra uniformemente las caras no-END con
Mylar aluminizado configurable.

## Nota de coordenadas

En el código real la barra de 1.4 m está orientada sobre **X** y los END están
en `-X/+X`. Por eso las caras no-END envueltas son `-Y/+Y/-Z/+Z`, aunque el
encargo original llamaba extremos a `±Z`. Envolver `±X` habría modificado los
END, lo cual estaba explícitamente prohibido. La discrepancia y la geometría
base completa están documentadas en `docs/AUDIT_endonly.md`.

## Modos de geometría

El hook nuevo es:

```text
/ship/geom/topSurface mylar|sipm
```

- `mylar` es el default de esta rama: no coloca TOP, no registra canales TOP y
  usa un panel `+Y` sólido.
- `sipm` recupera el comportamiento EndTop de la rama base, incluidas las 70
  ventanas y placements TOP. Para EndTop completo se usa junto con
  `/det/readout EndTop`.

Hooks Mylar y defaults:

```text
/ship/geom/mylar/reflectivity 0.90
/ship/geom/mylar/specularLobe 1.0
/ship/geom/mylar/sigmaAlpha 0.1 deg
```

La superficie común usa modelo `unified`, tipo `dielectric_metal`, finish
`ground`, `SPECULARLOBECONSTANT=specularLobe`,
`SPECULARSPIKECONSTANT=0`, `BACKSCATTERCONSTANT=0`; la fracción restante es
Lambertiana.

## Threads

El ejecutable acepta `-t N` y respeta esta prioridad:

```text
CLI -t N > G4FORCENUMBEROFTHREADS > default 1
```

Ejemplo:

```bash
./build/ej200_bar_sim -t 16 -m macros/endonly_guard.mac
```

## Barrido

Los scripts configuran/compilan si hace falta, crean siempre un directorio
timestamp nuevo bajo `output/`, validan cada ROOT sin sobrescribir archivos,
guardan logs y `run_metadata.txt`, y ejecutan el análisis END-only al terminar.

MSI, 16 threads:

```bash
scripts/run_msi.sh
```

t0minidaq, 24 threads:

```bash
git fetch && git checkout feat/endonly-mylar && scripts/run_t0minidaq.sh
```

Opciones comunes:

```bash
scripts/run_scan.sh --threads 16 --host msi --events 2000 \
  --positions scripts/endonly_positions.txt
```

Los defaults Mylar del script pueden variarse sin editar fuente:

```bash
MYLAR_REFLECTIVITY=0.92 MYLAR_SPECULAR_LOBE=0.8 \
MYLAR_SIGMA_ALPHA_DEG=0.2 scripts/run_msi.sh
```

Smoke de una posición:

```bash
printf '0\n' >/tmp/endonly_center.txt
scripts/run_scan.sh --threads 16 --host smoke --events 200 \
  --positions /tmp/endonly_center.txt
```

## Salidas de análisis

Cada directorio de barrido contiene `analysis/` con:

- `attenuation_curve.csv`: Npe END izquierdo/derecho frente a X y distancia.
- `attenuation_fits.csv`: ajustes de `lambda_eff` izquierdo, derecho y combinado.
- `sigma_t_sum4.csv`: `sigma(deltaT_LR)/sqrt(2)` frente a X.
- `analysis_metadata.txt`: mapa SUM4, pulso, umbral, estimador y definición de t=0.

El análisis usa exclusivamente IDs END `0..15`; si aparecen hits TOP los ignora
para las métricas END y los contabiliza en `n_top_hits`. Requiere Python con
NumPy, SciPy y uproot.

## Guardrails

```bash
cmake -S . -B build
cmake --build build --parallel 16
ctest --test-dir build --output-on-failure
```

Los cuatro guardrails nuevos verifican:

1. exactamente 16 canales END y cero TOP en `mylar`, con `+Y/-Y/+Z/-Z`
   cerradas por una única superficie Mylar;
2. presupuesto de fotones corto con cero detecciones TOP;
3. orden `sigma_intrínseca < sigma_total_modelada`, donde el total añade en
   cuadratura el jitter de lectura default de 20 ps;
4. paridad anti-artefacto y paridad byte a byte con la base para bloque END,
   SiPMSD, EJ-204, PDE, t=0 y análisis temporal.

Los 88 ps históricos pertenecen a otro detector/hardware y siguen siendo una
comparación cualitativa; no se usan como “total medido” de esta simulación.
