# EXEC_27b + EXEC_28 — chained campaign

| Gate | Literal criterion | Measured | Result |
|---|---|---|---|
| G1 | El residual debe cerrar por debajo del 2 %. Se espera que OpAbsorption en los volúmenes Vikuiti sea el componente dominante, y OpAbsorption en el bulk del centellador el segundo. Si no cierra al 2 % con telemetría completa, existe un tercer defecto no identificado. | residual=0.00000000%; starts=10100121; terminals=10100121; duplicates=0; unknown=0 | PASS |
| G2.1 | 1. **Control positivo, gratis y obligatorio.** En el brazo B, el censo de    fronteras del par AirGap*PV→Vikuiti*PV debe mostrar estados de    reflexión ≥ 95 %. Si no los muestra, la fábrica dielectric_metal no    está reflejando y el A/B no mide lo que dice medir: ABORTA la campaña    y reporta. | YMinus=97.996734%; ZPlus=97.999104%; ZMinus=98.000064%; ALL=97.999521% | PASS |
| G2.2 | 2. **Ganancia.** N_pe/end del brazo B dividido por el del brazo A debe    caer entre 2 y 5. Fuera de esa banda se reporta como tal, sin ajustar. | B/A=1.01154746 | FAIL |
| G2.3 | 3. **Asimetría de mecanismo.** La ganancia en N_pe/end debe superar la    ganancia en N_pe TOP total con razón mayor que 1,5. Si ambos suben por    igual, el mecanismo no es el rescate de fotones no-TIR. | END gain/TOP gain=0.84343045; TOP gain=1.19932528 | FAIL |
| G2.4 | 4. **Convergencia con el deck.** El valor de referencia es    \NapNpesimejB = 941 pe/end (EJ-204, x=0). Si el brazo B cae entre 800 y    1100, los dos defectos explican la totalidad del déficit histórico.    Fuera de ese intervalo, queda un efecto residual por identificar. | B=397.152500 ± 2.409927 pe/end; reference 941 | FAIL |

**Measured verdict.** All 15 cells completed, totaling 16,500 generated events. G1 closes at zero residual, and B0 passes the forward-boundary reflection control at 97.999521%. However, at x=0, END gain is only **1.011547 ± 0.008303**, while TOP gain is **1.199325**. G2.2, G2.3 and G2.4 fail. The two implemented corrections do **not** explain the full historical deficit under the preregistered 800–1100 pe/end test: B0 measures **397.1525 ± 2.4099 pe/end**, versus the stated 941 reference. The N=1000 position pairs do show END gain increasing with |x|, from about 1.02 at ±200 mm to about 1.28 at ±650 mm; this spatial dependence does not reverse the failed center tests. Evidence: C0 (N=500), A0/B0 (N=2000 each), and the twelve position cells (N=1000 each), all with seeds `26092601 8349041`; exact commands and provenance appear below. Gate and gain values use the supplied point-estimate criteria without post-hoc tolerances.

Status: **COMPLETE**, updated 2026-09-10T08:13:51.265544+00:00. Gate table sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/gates.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/gates.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/gates.meta.json). All run claims below refer to cell IDs in the manifest and exact commands in the appendix; every cell uses seeds `26092601 8349041`, one worker, EJ-204/OPSC-101, EndTop 70 TOP, vertical 1 GeV mu−. Only x and requested N vary.

G1 failure qualifies mechanism attribution but does not stop the campaign. Only G2.1 failure aborts. No criteria are tuned after seeing results.

## Block 1: disjoint terminal ledger

C0: N=500, commit `a7a909045653719087b6f755cb9bfa39e8721a25`. Generated scintillation=10,100,121; terminal scintillation=10,100,121; signed residual=0 (0.00000000%). Independent tracking starts and terminal identities both equal 10,100,121 for scintillation and 10,646,879 for all optical photons. Duplicate terminal callbacks, nonterminal callbacks, terminal-without-start, unknown process and unknown volume counts are all zero. The 546,758 non-scintillation optical photons are excluded from the scintillation ledger.

| Last-step defining process | Track logical volume | Explicit stepping kill flag | Count | Fraction of generated |
|---|---|---|---:|---:|
| Transportation | BarLV | none | 3,662,780 | 36.264714% |
| OpAbsorption | BarLV | none | 1,836,800 | 18.185921% |
| Transportation | BarLV | world_guard | 1,310,891 | 12.978963% |
| OpAbsorption | VikuitiZPlusLV | none | 1,104,103 | 10.931582% |
| OpAbsorption | VikuitiZMinusLV | none | 1,102,777 | 10.918453% |
| OpAbsorption | VikuitiYMinusLV | none | 967,815 | 9.582212% |
| Transportation | AirGapZPlusLV | none | 27,072 | 0.268036% |
| Transportation | AirGapZMinusLV | none | 26,514 | 0.262512% |
| Transportation | TopSiPMLV | none | 25,399 | 0.251472% |
| Transportation | AirGapYMinusLV | none | 22,550 | 0.223265% |
| Transportation | TopSiPMLV | world_guard | 4,969 | 0.049197% |
| Transportation | EndSiPMLV | none | 4,862 | 0.048138% |
| Transportation | AirGapZPlusLV | world_guard | 1,251 | 0.012386% |
| Transportation | AirGapZMinusLV | world_guard | 1,209 | 0.011970% |
| Transportation | EndSiPMLV | world_guard | 1,128 | 0.011168% |
| Transportation | AirGapYMinusLV | world_guard | 1 | 0.000010% |

Ledger sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/terminal_fates_run0_scintillation.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/terminal_fates_run0_scintillation.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/terminal_fates_run0_scintillation.meta.json). The terminal map counts unique (run,event,track) identities, selects creator `Scintillation`, and excludes suspended/nonterminal callbacks until termination. All-optical and identity-audit tables are separate. `none` means no explicit SteppingAction kill flag; it is not inferred from the last-step process. The process field is the **step-defining process**, which need not be the process that absorbed/detected a photon at a forced optical boundary. The volume is `track->GetVolume()` at the terminal callback; it is not a guessed destination. Flags identify issued SteppingAction kills, not independently proven priority over earlier kill decisions. Unknown categories are retained.

Observed OpAbsorption breakdown in C0: BarLV=1,836,800; VikuitiZPlusLV=1,104,103; VikuitiZMinusLV=1,102,777; VikuitiYMinusLV=967,815. The combined Vikuiti OpAbsorption count is 3,174,695 (31.432247% of generated scintillation), followed by BarLV OpAbsorption at 1,836,800 (18.185921%) within the bulk-absorption groups. The largest raw key is instead Transportation/BarLV/none at 3,662,780; it pools terminal effects not distinguished by the step-defining process. Thus G1 PASS refers to exact closure, while the absorption ordering is confirmed within OpAbsorption, not asserted as a ranking over heterogeneous terminal groups. No missing sink was imputed.

Instrumentation neutrality against EXEC_27 R2: {
  "rng_all_begin_end_identical": true,
  "event_yields_identical": true,
  "reference": "/home/rrios/ej200_exec26_20260909/build_exec27_20260910/run500",
  "reference_commit": "218241a571489db3c2ebed04dee15eb8cc98a663"
}

## Block 2: center A/B, N=2000 per arm

| x (mm) | N per arm | A N_pe/end ± SEM | B N_pe/end ± SEM | B/A ± approximate SEM | A TOP total | B TOP total | TOP B/A |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2000 | 392.6187 ± 2.1701 | 397.1525 ± 2.4099 | 1.01155 ± 0.00830 | 7,930,138 | 9,510,815 | 1.19933 |

This table and the scan table are display views of [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/ab_comparison.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/ab_comparison.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/ab_comparison.meta.json). A0 uses source `a7a909045653719087b6f755cb9bfa39e8721a25`; B0 uses `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. Both use seeds `26092601 8349041` and the commands below. The factory control is measured, not inferred from its name: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/reflection_panels.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/reflection_panels.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/reflection_panels.meta.json). Forward panel denominators exclude Undefined, NotAtBoundary, SameMaterial and StepTooSmall, all retained in the raw census; the reverse direction is not mixed. Every B0 panel exceeds 95%.

The disjoint terminal census directly shows a strong change in losses despite the small center END gain:

| Scintillation terminal group | A0 count | A0 fraction | B0 count | B0 fraction |
|---|---:|---:|---:|---:|
| OpAbsorption in BarLV | 7,286,978 | 18.188804% | 9,099,837 | 22.727130% |
| OpAbsorption in combined Vikuiti volumes | 12,594,570 | 31.436924% | 197,707 | 0.493779% |
| Other terminals; last-step process Transportation; no stepping flag | 14,950,059 | 37.316388% | 22,188,395 | 55.416215% |
| Explicit SteppingAction world_guard | 5,231,378 | 13.057884% | 8,553,596 | 21.362875% |

Sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/terminal_group_comparison.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/terminal_group_comparison.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/terminal_group_comparison.meta.json). Denominators are the measured scintillation totals, A0=40,062,985 and B0=40,039,535. Each column closes exactly. These are population comparisons, not matched counterfactual photon trajectories. In particular, the Transportation remainder is not renamed to a guessed physical absorber.

## Block 3: complete position scan, N=1000 per arm

| x (mm) | N per arm | A N_pe/end ± SEM | B N_pe/end ± SEM | B/A ± approximate SEM | A TOP total | B TOP total | TOP B/A |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 1000 | 672.4640 ± 5.0485 | 745.2045 ± 6.4346 | 1.10817 ± 0.01268 | 3,541,330 | 4,410,423 | 1.24541 |
| -500 | 1000 | 681.3290 ± 5.2717 | 733.4575 ± 5.6449 | 1.07651 ± 0.01175 | 3,576,651 | 4,339,618 | 1.21332 |
| 200 | 1000 | 432.4315 ± 3.3514 | 440.9770 ± 3.3935 | 1.01976 ± 0.01114 | 3,927,817 | 4,749,281 | 1.20914 |
| -200 | 1000 | 430.2955 ± 3.2128 | 440.3965 ± 3.1994 | 1.02347 ± 0.01066 | 3,915,409 | 4,740,688 | 1.21078 |
| 650 | 1000 | 993.1925 ± 7.8011 | 1267.3105 ± 9.4466 | 1.27600 ± 0.01382 | 3,116,704 | 3,807,703 | 1.22171 |
| -650 | 1000 | 1003.2155 ± 8.2454 | 1286.0155 ± 10.4494 | 1.28189 ± 0.01482 | 3,146,727 | 3,859,947 | 1.22665 |

The gain is not flat in these measurements: both ±200 mm values are near 1.02, both ±500 mm values are higher (1.07651 and 1.10817), and both ±650 mm values are about 1.28. The preregistered weighted fit over the six N=1000 points gives slope **0.000484544 ± 0.000026569 per mm**, which exceeds two approximate SEM. The constant fit gives χ²/dof = **418.393/5**. A straight-line fit also fits poorly, χ²/dof = **85.793/4**; the observed rise is curved, so the linear slope is a trend summary, not a validated response model. Sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/position_fit.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/position_fit.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/position_fit.meta.json).

The observed growth with |x| matches the requested spatial signature. It alone does not identify which individual trajectories were rescued or establish the proposed selective END mechanism: G2.3 fails at the center, and the scan END/TOP gain ratios remain approximately 0.84–1.05. Changing |x| also changes distances to the two ends; |x| is not itself a unique optical path length. No causal attribution beyond the measured comparisons is claimed.

Event means use `(L_i+R_i)/2` over all generated events, including zeros; uncertainties are sample SEM. B/A error propagation and the weighted fits neglect unmeasured cross-cell covariance despite equal seeds, so their uncertainty/significance statements are approximate. Boundary entries are encounters, not independent photons.

Time budget: the first N=2000 simulation took 965.657 s. Its initial projection was 1.609 h for Block 3 and 2.214 h for all planned simulations including C0. After B0, the estimated remaining Block 3 time was 2.210 h. No positions were cut and N was never reduced. Actual summed simulation wall time was **10139.768 s (2.817 h)**; first C0 start to final cell sidecar completion was **10627.868 s (2.952 h)**. Plan: [time_plan.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/time_plan.json).

## Incremental manifest

| cell_id | arm | x_mm | N | commit | seeds | start_utc | end_utc | exit_code | npe_end_mean | npe_end_sem | npe_top_total | wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C0 | A | 0 | 500 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T00:46:51.994661+00:00 | 2026-09-10T00:50:57.191844+00:00 | 0 | 395.394 | 4.509400727129456 | 1999208 | 245.19695864105597 |
| A0 | A | 0 | 2000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T00:52:16.858888+00:00 | 2026-09-10T01:08:22.516399+00:00 | 0 | 392.61875 | 2.1700957454663135 | 7930138 | 965.6572835349943 |
| B0 | B | 0 | 2000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T01:09:07.399450+00:00 | 2026-09-10T01:37:13.807945+00:00 | 0 | 397.1525 | 2.40992722846526 | 9510815 | 1686.4081590180285 |
| A_p500 | A | 500 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T01:38:05.931718+00:00 | 2026-09-10T01:45:23.658783+00:00 | 0 | 672.464 | 5.0484659761865 | 3541330 | 437.7268402610207 |
| B_p500 | B | 500 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T01:45:46.536363+00:00 | 2026-09-10T01:58:44.992521+00:00 | 0 | 745.2045 | 6.434595943391062 | 4410423 | 778.4559326800518 |
| A_m500 | A | -500 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T01:59:12.790780+00:00 | 2026-09-10T02:06:34.361975+00:00 | 0 | 681.329 | 5.271671867945053 | 3576651 | 441.57087370799854 |
| B_m500 | B | -500 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T02:06:58.138411+00:00 | 2026-09-10T02:19:53.882721+00:00 | 0 | 733.4575 | 5.644879459107861 | 4339618 | 775.7440879539354 |
| A_p200 | A | 200 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T02:20:21.279030+00:00 | 2026-09-10T02:28:19.865870+00:00 | 0 | 432.4315 | 3.351429043614736 | 3927817 | 478.5866089130286 |
| B_p200 | B | 200 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T02:28:42.744833+00:00 | 2026-09-10T02:42:38.813795+00:00 | 0 | 440.977 | 3.393521247532931 | 4749281 | 836.068720121053 |
| A_m200 | A | -200 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T02:43:05.856581+00:00 | 2026-09-10T02:51:15.289992+00:00 | 0 | 430.2955 | 3.21276263883959 | 3915409 | 489.4331713520223 |
| B_m200 | B | -200 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T02:51:38.275365+00:00 | 2026-09-10T03:05:27.918555+00:00 | 0 | 440.3965 | 3.1993672012193657 | 4740688 | 829.6429276859853 |
| A_p650 | A | 650 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T03:05:54.675673+00:00 | 2026-09-10T03:12:37.869827+00:00 | 0 | 993.1925 | 7.801116594043693 | 3116704 | 403.1939286349807 |
| B_p650 | B | 650 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T03:13:02.341888+00:00 | 2026-09-10T03:24:25.933532+00:00 | 0 | 1267.3105 | 9.446601961995562 | 3807703 | 683.5913884219481 |
| A_m650 | A | -650 | 1000 | a7a909045653719087b6f755cb9bfa39e8721a25 | 26092601 8349041 | 2026-09-10T03:24:56.308142+00:00 | 2026-09-10T03:31:37.272980+00:00 | 0 | 1003.2155 | 8.245408453819637 | 3146727 | 400.9646073059412 |
| B_m650 | B | -650 | 1000 | f4d90a6fb03eab4f9a81adc81cfeab09824180d0 | 26092601 8349041 | 2026-09-10T03:32:02.045533+00:00 | 2026-09-10T03:43:29.572169+00:00 | 0 | 1286.0155 | 10.449423163706417 | 3859947 | 687.5263987870421 |

Manifest sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/manifest.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/manifest.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/manifest.meta.json). A completion row is appended only after that cell’s analysis sidecars are written and checked. Failed invocations retain their exit code and blank yield fields. `active_cell.json` records an in-flight cell.

## Appendix: exact commands, timing, RNG and provenance

**A0**, arm A, x=0 mm, N=2000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T00:52:16.858888+00:00 → 2026-09-10T01:08:22.516399+00:00; exit 0; simulation wall 965.657284 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "6cce39575e2bc45faba50b5c2fc3a94d6482bb96fef28d6fd53e9be13f214220", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "7ee82cad693586741975198b72034d54252b25e1f4069e8b98148f7ebdf1af04"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A0/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_m200**, arm A, x=-200 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T02:43:05.856581+00:00 → 2026-09-10T02:51:15.289992+00:00; exit 0; simulation wall 489.433171 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "4f8d7f7e8e43fbc44d0d3b4257029e72617510611c80b8f9dfcad84bca7881b2"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m200/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_m500**, arm A, x=-500 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T01:59:12.790780+00:00 → 2026-09-10T02:06:34.361975+00:00; exit 0; simulation wall 441.570874 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "6367f97d82cdc777dfd1fb15ceb8eb7f56324f46747f2648f5cd0efc56cdd3fc"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m500/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_m650**, arm A, x=-650 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T03:24:56.308142+00:00 → 2026-09-10T03:31:37.272980+00:00; exit 0; simulation wall 400.964607 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "44f1f786573a9a9389762a43014494309c177ee6751629c5a0f908ab55106106"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_m650/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_p200**, arm A, x=200 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T02:20:21.279030+00:00 → 2026-09-10T02:28:19.865870+00:00; exit 0; simulation wall 478.586609 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "ee3fe147471fcdb4622759e37e1a158b9246f45689eb5e51ec4a251d42083b51"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p200/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_p500**, arm A, x=500 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T01:38:05.931718+00:00 → 2026-09-10T01:45:23.658783+00:00; exit 0; simulation wall 437.726840 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "1728c346b22f0f7d86e59d2a60fd1f170ebfeb9dd1f55510dbd273f91866910f"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p500/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**A_p650**, arm A, x=650 mm, N=1000; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T03:05:54.675673+00:00 → 2026-09-10T03:12:37.869827+00:00; exit 0; simulation wall 403.193929 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "00cb1981383803149bd88997fb9757d5263d5c8f07b66ad09d13b231ea19a110"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/A_p650/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B0**, arm B, x=0 mm, N=2000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T01:09:07.399450+00:00 → 2026-09-10T01:37:13.807945+00:00; exit 0; simulation wall 1686.408159 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "6cce39575e2bc45faba50b5c2fc3a94d6482bb96fef28d6fd53e9be13f214220", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "d772975da08468fc072fd9165ba6323f087ab018477d6e591788c386cdad4b33"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B0/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_m200**, arm B, x=-200 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T02:51:38.275365+00:00 → 2026-09-10T03:05:27.918555+00:00; exit 0; simulation wall 829.642928 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "ea77e214142f068f65e9ae11a8ec1d25786786f9f3688ddf1ba6f89ecd09c677"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m200/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_m500**, arm B, x=-500 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T02:06:58.138411+00:00 → 2026-09-10T02:19:53.882721+00:00; exit 0; simulation wall 775.744088 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "cc0bd9b090a7ff4d412928e12680ab4f22954d5be076fd92b2e9772451291aba"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m500/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_m650**, arm B, x=-650 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T03:32:02.045533+00:00 → 2026-09-10T03:43:29.572169+00:00; exit 0; simulation wall 687.526399 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "041cbc9f656aaa8eae27c554e6be3b0a76483e8659a487d347bd69716d6354e0"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_m650/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_p200**, arm B, x=200 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T02:28:42.744833+00:00 → 2026-09-10T02:42:38.813795+00:00; exit 0; simulation wall 836.068720 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "3cfaf47e493a62ef15d9c2ee30604bceea2ca934a94f5ca19f5038fd1f6e81ae"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p200/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_p500**, arm B, x=500 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T01:45:46.536363+00:00 → 2026-09-10T01:58:44.992521+00:00; exit 0; simulation wall 778.455933 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "aa13a26531a52748dfdf65692e0d9a4a885a1deff8ced9eabe4ca3f644837dc3"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p500/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**B_p650**, arm B, x=650 mm, N=1000; commit `f4d90a6fb03eab4f9a81adc81cfeab09824180d0`. UTC 2026-09-10T03:13:02.341888+00:00 → 2026-09-10T03:24:25.933532+00:00; exit 0; simulation wall 683.591388 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "9a5446c337137f257669864c32189522d6fc3217db4d4247c28ae2cb369c90f9", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "5aa86095b5876e571309162464482bbd40976557cdc152afb455d58c0f7791ab"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/B_p650/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

**C0**, arm A, x=0 mm, N=500; commit `a7a909045653719087b6f755cb9bfa39e8721a25`. UTC 2026-09-10T00:46:51.994661+00:00 → 2026-09-10T00:50:57.191844+00:00; exit 0; simulation wall 245.196959 s. cwd `/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0`.

```bash
/home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/run.mac
```

RNG hashes: `{"rng_run0_thread-1_begin.rndm": "5adc16913fd73fbeeaceb380e8381a7a31b63758a90f41990636980eac860a95", "rng_run0_thread-1_end.rndm": "2ffb028328522836ecce2e99ca8020627cfdb199731163ec364dcdb477a3a745", "rng_run0_thread0_begin.rndm": "b6c1b5a824b79769d5f128ec867d79593b8661cec0f465378a7416ffb89976f0", "rng_run0_thread0_end.rndm": "5a272d98f4892d783359c3959baf30ed11a864c773a53c7aef77604fe97ced50"}`.

Cell metric sidecars: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/cell.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/cell.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/cell.meta.json); [metadata](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/cell_summary.json); [log](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/cells/C0/run.log). All tables in the cell have ROOT/CSV/JSON counterparts; original hit data are preserved separately.

Build/edit audit: [commands.jsonl](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/commands.jsonl); preregistration: [preregister.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/preregister.json); reproducible executor: [campaign.py](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/campaign.py). Geant4 identity is `geant4-11-04 [MT] (5-December-2025)`, version 11.4.0, installation `/home/tdship/opt/geant4-v11.4.0-install`. Exact arm source/binary hashes are stored before use in [arms.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/arms.json).

C001: 2026-09-10T00:40:01.378619+00:00 → 2026-09-10T00:40:01.381117+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_001.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 tag pre-exec27b-nightly-20260910 218241a
```

C002: 2026-09-10T00:41:13.118852+00:00 → 2026-09-10T00:41:13.131323+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_002.log).

```bash
python /home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/instrument.py
```

C003: 2026-09-10T00:41:31.634629+00:00 → 2026-09-10T00:41:32.583885+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_003.log).

```bash
cmake -S /home/rrios/ej200_exec26_20260909 -B /home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A -DCMAKE_BUILD_TYPE=Release
```

C004: 2026-09-10T00:41:32.607261+00:00 → 2026-09-10T00:41:35.828924+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_004.log).

```bash
cmake --build /home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_A --target ej200_bar_sim -j 8
```

C005: 2026-09-10T00:41:57.048048+00:00 → 2026-09-10T00:41:57.052341+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_005.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 diff --check
```

C006: 2026-09-10T00:41:57.067186+00:00 → 2026-09-10T00:41:57.070805+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_006.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 diff --stat
```

C007: 2026-09-10T00:41:57.085382+00:00 → 2026-09-10T00:41:57.088859+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_007.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 add include/TrackingAction.hh src/TrackingAction.cc src/ActionInitialization.cc src/SteppingAction.cc src/RunAction.cc
```

C008: 2026-09-10T00:41:57.103717+00:00 → 2026-09-10T00:41:57.110685+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_008.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 commit -m 'feat(diag): terminal-fate census by process and volume (EXEC_27b)'
```

C009: 2026-09-10T00:42:38.349588+00:00 → 2026-09-10T00:42:38.369915+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_009.log).

```bash
python /home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/setup_campaign.py
```

C010: 2026-09-10T00:46:51.804180+00:00 → 2026-09-10T00:46:51.836245+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_010.log).

```bash
python -
```

C011: 2026-09-10T00:46:51.851065+00:00 → 2026-09-10T00:51:08.931499+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_011.log).

```bash
python /home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/campaign.py c0
```

C012: 2026-09-10T00:51:16.000534+00:00 → 2026-09-10T00:51:16.002167+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_012.log).

```bash
grep -n 'auto\* airReflectorSurface' /home/rrios/ej200_exec26_20260909/src/DetectorConstruction.cc
```

C013: 2026-09-10T00:51:46.022689+00:00 → 2026-09-10T00:51:46.034035+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_013.log).

```bash
python -
```

C014: 2026-09-10T00:51:46.050952+00:00 → 2026-09-10T00:51:46.053922+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_014.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 diff -- src/DetectorConstruction.cc
```

C015: 2026-09-10T00:51:46.070038+00:00 → 2026-09-10T00:51:46.072995+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_015.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 add src/DetectorConstruction.cc
```

C016: 2026-09-10T00:51:46.087199+00:00 → 2026-09-10T00:51:46.094234+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_016.log).

```bash
git -C /home/rrios/ej200_exec26_20260909 commit -m 'fix(optics): bind dielectric_metal reflector at air-wrap border — DetectorConstruction.cc:313 before→after'
```

C017: 2026-09-10T00:51:46.110578+00:00 → 2026-09-10T00:51:47.062861+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_017.log).

```bash
cmake -S /home/rrios/ej200_exec26_20260909 -B /home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B -DCMAKE_BUILD_TYPE=Release
```

C018: 2026-09-10T00:51:47.086946+00:00 → 2026-09-10T00:51:50.331473+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_018.log).

```bash
cmake --build /home/rrios/ej200_exec26_20260909/build_nightly_20260910/arm_B --target ej200_bar_sim -j 8
```

C019: 2026-09-10T00:52:16.635003+00:00 → 2026-09-10T00:52:16.664392+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_019.log).

```bash
python -
```

C020: 2026-09-10T00:52:16.679985+00:00 → 2026-09-10T03:44:00.058765+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_020.log).

```bash
python /home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/campaign.py full
```

C021: 2026-09-10T08:10:40.733881+00:00 → 2026-09-10T08:10:49.826816+00:00; exit 0; cwd `/home/rrios/ej200`; [output](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/command_021.log).

```bash
python /home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/validate_campaign.py
```

Final gate: no deck/presentations changes, no merge to main and no push are authorized or performed. Main is retained at 8349041. No further correction follows this campaign without explicit approval.

## Final repository and artifact checks

Exactly two commits were created after `218241a`: the telemetry commit `a7a9090`, then the one-line factory binding commit `f4d90a6`. Rollback tag `pre-exec27b-nightly-20260910` remains at `218241a`. Both arms were built cleanly before measurement; A's frozen executable was preserved while the worktree advanced to B. C0 reproduces all 500 event yields and all four RNG files from EXEC_27 R2. No simulation reruns, parameter tuning, geometry/material/gun-filter changes beyond requested x positions, or deck/presentations changes were made. The original material mass-fraction warning and legacy reflector banner are retained; the latter is not used as a physical measurement. Main remains clean at `8349041`; the diagnostic worktree is clean at `f4d90a6`.

All 120 cell table/ntuple triples and the aggregate tables have hashes, required provenance and checked ROOT entry counts; analytical table values match CSV exactly. The complete raw-hit CSV exports are preserved, with hashes and entry counts. Sidecar completion precedes the next simulation in every cell: [csv](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/sidecar_completion_audit.csv) / [root](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/sidecar_completion_audit.root) / [meta.json](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/sidecar_completion_audit.meta.json). [Validation results](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/validation.json) and [validator](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/validate_campaign.py) are retained. Report finalization: [finalize_report.py](/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit/finalize_report.py). No figures were produced. The campaign is complete and stops at the requested final gate.
