# EJ-230 END-only + Mylar

Replica para EJ-230 de la geometria `feat/endonly-mylar`, partiendo de
`feat/ej230-sslg4`.

- 16 SiPM END, IDs globales `0..15`, con la colocacion y el acoplamiento de la
  rama base.
- EJ-230/SSLG4 `OPSC-106` como centellador por defecto.
- Mylar fijo sobre las caras `-Y`, `+Y`, `-Z` y `+Z`.
- Las caras de extremo `-X` y `+X` quedan exclusivamente para los SiPM END.
- No existen SiPM TOP, ventanas, recortes ni modos alternativos TOP.

La superficie Mylar copia los defaults de `feat/endonly-mylar`: modelo
`unified`, tipo `dielectric_metal`, finish `ground`, reflectividad uniforme
`0.90`, `SPECULARLOBECONSTANT=1.0`, spike/backscatter `0` y
`sigmaAlpha=0.1 deg`. Estos valores no son configurables en esta rama.

## Build sin GDML

GDML es opcional y esta desactivado por defecto:

```bash
cmake -S . -B build -DUSE_GDML=OFF
cmake --build build --parallel 16
ctest --test-dir build --output-on-failure
```

## Barrido

MSI usa 16 threads:

```bash
scripts/run_msi.sh
```

`t0minidaq` usa 24 threads:

```bash
scripts/run_t0minidaq.sh
```

Ambos runners configuran y compilan con `USE_GDML=OFF`, usan las posiciones
longitudinales de `scripts/endonly_positions.txt`, crean un subdirectorio
timestamp nuevo bajo `output/` y ejecutan el analisis al terminar.

Para un smoke corto:

```bash
printf '0\n' >/tmp/ej230_endonly_center.txt
scripts/run_msi.sh --events 10 --positions /tmp/ej230_endonly_center.txt
```

Cada sesion produce ROOT por posicion y, bajo `analysis/`:

- `attenuation_curve.csv` y `attenuation_curve.png`;
- `attenuation_fits.csv`;
- `sigma_t_sum4.csv` y `sigma_t_sum4.png`;
- `analysis_metadata.txt`.

El timing usa SUM4 `{0..3}`, `{4..7}`, `{8..11}`, `{12..15}`, leading edge,
`t_L-t_R` y `sigma_single = sigma_LR / sqrt(2)`.
