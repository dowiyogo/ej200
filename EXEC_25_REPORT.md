# EXEC_25 - PDE semantics and optical realism bracket

## Veredicto

COMPLETADO-CON-FLAGS.

## Resumen ejecutivo

- PDE surface `EFFICIENCY` quedo validada: `zero` produjo 0 hits, `one` produjo recorded/entering = 1.018, y `nominal` produjo recorded/entering = 0.626 frente a mean PDE branch = 0.616.
- Queda prohibido aplicar otra Bernoulli PDE en `SiPMSD`: eso duplicaria PDE.
- El exceso de 2035 PE no lo explica el reflector lateral ideal por si solo. Incluso con `REFLECTIVITY=0`, el run dio 1967.8 PE/event.
- Roughness/finish tampoco corrige: `ground` a R=0.95 dio 1915.15 PE/event.
- Cambiar el air gap entre 0.02 y 0.20 mm mantiene Npe en 1781.95-2269.45 PE/event.
- La ventana temporal no rescata el resultado: en 5 ns quedan 1456-1561 PE/event segun configuracion.
- No hay configuracion candidata para EXEC_26. El siguiente paso fisico es modelar acoplo/SiPM/cobertura o wrapping/contacto parcial no idealizado antes de cualquier scan.

## T1 PDE semantics

| mode | entering Bar->SiPM | recorded hits | recorded/entering | expected |
|---|---:|---:|---:|---|
| nominal | 31045 | 19429 | 0.6258 | close to mean PDE branch 0.6164 |
| zero | 36403 | 0 | 0.0000 | zero hits |
| one | 33174 | 33786 | 1.0184 | approximately all crossings recorded |

Veredicto: PASS.

The `one` ratio is slightly above 1 because the boundary-census denominator is not an exact per-SD invocation counter, but it is close enough for the hard semantic test together with the exact zero-hit result at `EFFICIENCY=0`.

## T3 Reflectivity bracket

| R | Npe | recorded/generated | recorded/entering | sigma ps | runtime s |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 1967.80 | 0.09396 | 0.62558 | 14.88 | 169.9 |
| 0.50 | 1782.40 | 0.09389 | 0.62360 | 51.93 | 171.6 |
| 0.70 | 1882.50 | 0.09400 | 0.61992 | 49.74 | 186.5 |
| 0.85 | 1840.35 | 0.09372 | 0.61847 | 43.94 | 205.3 |
| 0.90 | 1927.30 | 0.09498 | 0.62286 | 53.85 | 247.4 |
| 0.95 | 1931.25 | 0.09806 | 0.62297 | 34.47 | 271.0 |
| 0.98 | 2129.40 | 0.10292 | 0.62309 | 52.21 | 349.7 |

Interpretation: reflectivity has a weak effect on Npe in this setup, but a strong effect on cost at high R. Because R=0 remains near 2000 PE/event, the excess is dominated by prompt/direct capture and effective SiPM/coupling acceptance, not by late recovery from the lateral reflector.

## T4 Roughness / finish bracket

| finish | sigma_alpha deg | Npe | recorded/generated | sigma ps | runtime s |
|---|---:|---:|---:|---:|---:|
| polished | 0.0 | 1881.15 | 0.09788 | 57.91 | 274.2 |
| polished | 0.1 | 1799.60 | 0.09721 | 20.99 | 251.7 |
| polished | 0.3 | 1888.45 | 0.09721 | 28.77 | 287.1 |
| ground | 0.1 | 1915.15 | 0.09702 | 54.53 | 294.8 |

Interpretation: changing finish/roughness at R=0.95 does not move Npe toward a physical 20-300 PE range.

## T5 Gap bracket

| air_gap mm | Npe | recorded/generated | sigma ps | runtime s |
|---:|---:|---:|---:|---:|
| 0.02 | 2061.00 | 0.09851 | 65.81 | 309.4 |
| 0.05 | 1985.00 | 0.09760 | 51.47 | 285.7 |
| 0.10 | 2269.45 | 0.09794 | 44.38 | 312.8 |
| 0.20 | 1781.95 | 0.09596 | 48.02 | 254.5 |

Interpretation: simple air-gap thickness is not the governing parameter. The issue is more likely the continuous idealized coupling/acceptance model than the exact 0.10 mm value.

## T6 Time-window diagnostic

Reference: `t0 = first END PE per event`.

| config | total Npe | Npe 5ns | Npe 10ns | Npe 20ns | Npe 50ns |
|---|---:|---:|---:|---:|---:|
| R050_polished | 1782.40 | 1455.85 | 1741.45 | 1781.10 | 1782.40 |
| R095_polished | 1931.25 | 1561.40 | 1885.65 | 1930.20 | 1931.25 |
| rough_polished_s0p1 | 1799.60 | 1457.95 | 1757.10 | 1798.60 | 1799.60 |

Interpretation: the excess is prompt. A 5-10 ns diagnostic window still leaves O(1500-1900) PE/event.

## Candidate selection

No candidate selected.

Reason: all bracket configurations remain far above the acceptable 20-300 PE range. The lowest observed point is 1781.95 PE/event (`gap_0p2mm`), still high by a factor of about 6 relative to the upper flagged bound of 300 PE.

This is not a tuning failure of R/finish/gap. The bracket says those parameters are not the dominant missing realism.

## Candidate mini-validation

Not run because no candidate satisfied T7.

## Flags

- `no_candidate_in_20_to_300_pe_range`
- `high_collection_configs_remain`
- `prompt_collection_excess`
- `reflectivity_not_primary_driver`
- `gap_thickness_not_primary_driver`

## Rutas

- `out/EXEC_25/pde_semantics_test.csv`
- `out/EXEC_25/pde_semantics_summary.json`
- `out/EXEC_25/bracket_results.csv`
- `out/EXEC_25/bracket_results.json`
- `out/EXEC_25/time_window_diagnostic.csv`
- `out/EXEC_25/results_exec25.json`
- `out/EXEC_25/png/reflectivity_bracket.png`
- `out/EXEC_25/png/roughness_bracket.png`
- `out/EXEC_25/png/gap_bracket.png`
- `out/EXEC_25/png/npe_vs_time_window.png`

Raw ROOT outputs were retained under `out/EXEC_25/raw/` for local inspection and are not staged for commit.

## Git

- Branch: `exec25-optical-realism-bracket`
- Tag: `EXEC_25-pre-20260621_165559`
- Commit: pending at report creation
- Push printed, not executed

## Proximo paso

Implement partial-contact or non-ideal wrapping/coupling before any END-only scan. The most likely next physics handle is not Mylar reflectivity, finish, or gap thickness, but the effective prompt acceptance from 48% active end coverage plus perfect bar-SiPM optical coupling.
