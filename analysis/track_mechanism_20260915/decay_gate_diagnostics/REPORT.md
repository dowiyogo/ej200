# EXEC_46 — D1–D4 investigation of the decay-tail gate

Only `EJ200_xp0` was analyzed. The existing gate remains unchanged and the full grid remains stopped.

## Finding

The global creation-time tail is composite: its mean excess increases from 2.4274 ns to 10.2882 ns as the cut increases from 6.3 ns to 42 ns. The runtime MPT has exactly one scintillation component with decay constant 2.1 ns. Photons created within 0.1 mm of the gun axis recover approximately 2.1 ns with the same estimator; spatially displaced photons strongly enrich the global late tail. Wavelength-shifting reemission is not supported by the recorded creator labels or the effective MPT. A unique 2.60 ns decay constant for the emission pool is therefore not justified.

The earlier explanation based on a constant time-origin offset is withdrawn. Exponential memorylessness makes that explanation insufficient. A varying, spatially correlated production-time distribution is a distinct hypothesis; the current tree does not store parent IDs, parent species, or the parent deposition time needed to identify that distribution directly.

## Scope and provenance

- Source ROOT: `/home/rrios/exec46_20260915/full_grid/cells/EJ200_xp0/attempts/c73247d595ca49f4aa18793200edab77/photon_hits_run000.root`
- SHA-256 recomputed from the source: `8e20d372d616e7d58d92c040d542b1691186367b7f667226de43fbac4b71f88e`.
- Simulation: `/home/rrios/exec46_20260915/build_baseline/ej200_bar_sim -m /home/rrios/exec46_20260915/full_grid/cells/EJ200_xp0/attempts/c73247d595ca49f4aa18793200edab77/run.mac`.
- Seeds: `[26092601, 8349041]`; N = 10,000 events.
- Detected photons: 62,596,038; scintillation: 60,368,164; Cherenkov: 2,227,874; primary/other: 0/0.
- No other cell was read. No event was generated. Geometry, physics, the existing gate, and the original Step 1 evidence were not edited.
- Rollback tag: `pre-exec46-d1-d4-20260915` at `7f9ff31`.

Commands (working directory `/home/rrios/ej200`):

```bash
python3 analysis/track_mechanism_20260915/run_d1_dump.py
python3 analysis/track_mechanism_20260915/diagnose_decay.py
python3 analysis/track_mechanism_20260915/decay_followup_summary.py
python3 analysis/track_mechanism_20260915/measure_emission_pool.py
```

Each invocation takes substantially less than 20 minutes. The final D2–D4 collection and figure generation took 31.60 s. The sufficient-statistics cache contains only this cell; no scan of the full grid was invoked.

## D1 — Configuration and effective MPT

The helper reuses the existing production `.o` files and SSLG4/OPSim static libraries, replacing only the executable entry point. It applies the actual attempt macro through initialization, stops before `/run/beamOn`, and dumps the MPT in memory. Generated events = 0. The production executable SHA matches the `.DONE` value. Repository physics sources have no diff relative to simulation commit `4967ec8`.

This is a runtime reconstruction with the preserved binary objects and current configuration files, not a memory snapshot archived during the historical run. `d1_provenance.json` records object/configuration hashes and exact compiler, linker, and invocation arguments. The historical records do not provide a separate original hash of every material table, so byte-for-byte historical table identity cannot be independently proved from those records alone.

| Property | Effective value / state |
|---|---|
| SCINTILLATIONCOMPONENT1 | Present; full spectrum in effective_mpt.json |
| SCINTILLATIONCOMPONENT2 / 3 | Absent / absent |
| SCINTILLATIONTIMECONSTANT1 / 2 / 3 | 2.1 ns / absent / absent |
| SCINTILLATIONYIELD1 / 2 / 3 | 1 / absent / absent |
| YIELDRATIO | Absent |
| SCINTILLATIONRISETIME1 | 0.9 ns; finite rise enabled |
| SCINTILLATIONYIELD | 10000 photons/MeV |
| ABSLENGTH | 3800 mm, at both tabulated energies |
| WLSABSLENGTH / WLSTIMECONSTANT / WLSCOMPONENT | All absent |
| WLSABSLENGTH2 / WLSTIMECONSTANT2 / WLSCOMPONENT2 | All absent |
| RINDEX | 1.58 at both tabulated energies |
| GROUPVEL | 189.74206202531644 mm/ns |

The active scintillation model has one positive decay component, with a finite rise. Its density is proportional to `exp(-t/tau_d)*(1-exp(-t/tau_r))`. The second, subtractive exponential associated with the rise is not an additional slow scintillation component; the mean emission delay of that isolated model is `tau_d + tau_d*tau_r/(tau_d+tau_r) = 2.73 ns`. This sampling form is implemented in local `G4Scintillation.cc:557`.

Full, unabridged macro:

```text
#EJ-200/Pilot F/BC-408 mpt
/control/alias scnt opsc-100
/control/alias pathToDataDir sslg4/data/oscnt/{scnt}

/mpt/{scnt}/addProperty SCINTILLATIONCOMPONENT1 {pathToDataDir}/scntComp1.txt nm unitless
/mpt/{scnt}/addProperty RINDEX {pathToDataDir}/rIndex.txt nm unitless
/mpt/{scnt}/addProperty ABSLENGTH {pathToDataDir}/absLength.txt nm cm

/mpt/{scnt}/addConstProperty SCINTILLATIONYIELD 10000 1/MeV
/mpt/{scnt}/addConstProperty SCINTILLATIONYIELD1 1
/mpt/{scnt}/addConstProperty SCINTILLATIONTIMECONSTANT1 2.1 ns
/mpt/{scnt}/addConstProperty SCINTILLATIONRISETIME1 0.9 ns
/mpt/{scnt}/addConstProperty RESOLUTIONSCALE 1.0
```

Full effective MPT: [machine-readable dump](/home/rrios/exec46_20260915/decay_diagnostics/effective_mpt.json) and [Geant4 DumpTable output](/home/rrios/exec46_20260915/decay_diagnostics/effective_mpt.log). The machine-readable dump includes every loaded vector entry and every defined constant. Geant4 internal units here are MeV, mm, ns.

## D2 — Mean excess versus cut

The estimator is exactly `mean(t_creation_ns - cut | t_creation_ns >= cut, source_type == 1)`, with no upper truncation. Errors below use event-cluster influence functions: photons from the same event are not treated as independent. The companion CSV also retains the photon-i.i.d. errors for comparison.

| Cut / nominal tau | Cut [ns] | Tail photons | Events with tail | tau_fit [ns] | Event-cluster SE [ns] |
|---:|---:|---:|---:|---:|---:|
| 3 | 6.3 | 5,503,662 | 10,000 | 2.427407 | 0.001103 |
| 4 | 8.4 | 2,275,956 | 10,000 | 2.504295 | 0.001765 |
| 5 | 10.5 | 961,890 | 10,000 | 2.601323 | 0.002957 |
| 8 | 16.8 | 85,247 | 9,993 | 3.207598 | 0.013979 |
| 12 | 25.2 | 6,389 | 4,653 | 5.268989 | 0.086063 |
| 20 | 42.0 | 319 | 310 | 10.288152 | 0.717588 |

All five adjacent changes are positive. Accounting for the shared events and nested tails, their differences divided by paired standard errors range from 7.32 to 59.43. This decisively selects the composite-distribution outcome of D2. The maximum recorded scintillation creation time is 111.275089 ns.

Figure: [cut stability](/home/rrios/exec46_20260915/decay_diagnostics/d2_cut_stability.pdf).

## D3 — Spatial locality and creator process

Pearson correlations use exact, unbinned full-range moments:

| Population | rho(abs(x_creation), t_creation) | Mean abs(x_creation) [mm] |
|---|---:|---:|
| all_hits | 0.318850 | 6.321163 |
| scintillation | 0.321988 | 6.552675 |
| Cherenkov | 0.526734 | 0.047946 |

The conditional means below are exact and untruncated; bins are half-open in absolute creation x. Their errors retain event clustering.

| abs(x_creation) [mm] | Photons | Mean t_creation [ns] | SE [ns] | Fraction t >= 10.5 ns |
|---|---:|---:|---:|---:|
| [0, 0.1) | 53,908,590 | 2.931015 | 0.000300 | 1.062% |
| [0.1, 0.5) | 2,191,848 | 2.934241 | 0.001493 | 1.066% |
| [0.5, 1) | 1,299,021 | 2.941280 | 0.001973 | 1.069% |
| [1, 2) | 941,115 | 2.949804 | 0.002342 | 1.101% |
| [2, 5) | 579,354 | 2.979654 | 0.003580 | 1.180% |
| [5, 10) | 121,487 | 3.970545 | 0.140095 | 4.541% |
| [10, 20) | 124,098 | 5.933241 | 0.200453 | 9.898% |
| [20, 50) | 135,539 | 6.166754 | 0.122688 | 10.666% |
| [50, 100) | 170,226 | 6.778097 | 0.064412 | 13.726% |
| [100, 200) | 201,708 | 7.542725 | 0.009565 | 17.732% |
| [200, 400) | 262,309 | 8.651661 | 0.007965 | 25.012% |
| [400, 700) | 432,869 | 10.505832 | 0.042198 | 41.113% |

Only 2.399% of detected scintillation photons originate at abs(x_creation) >= 5 mm, but they comprise 34.824% of the tail above 10.5 ns. Creation positions extend to abs(x) = 699.999978 mm. The spatially displaced late population is therefore measured directly.

Locality cross-check, using the same unbinned estimator on abs(x_creation) < 0.1 mm:

| Cut [ns] | Photons | tau_fit [ns] | Event-cluster SE [ns] |
|---:|---:|---:|---:|
| 6.3 | 4,218,664 | 2.102404 | 0.001028 |
| 8.4 | 1,553,310 | 2.103800 | 0.001676 |
| 10.5 | 572,488 | 2.105881 | 0.002778 |
| 16.8 | 28,815 | 2.110331 | 0.012358 |
| 25.2 | 547 | 2.018451 | 0.086837 |
| 42.0 | 1 | 0.996402 | undefined: one event |

The 6.3–16.8 ns cuts remain close to 2.1 ns for this local proxy. This spatial selection was introduced as a diagnostic in this task; it is not a truth-level selection of photons emitted by the primary muon.

**Deslocalization does not establish optical reemission in this dataset.** `SiPMSD.cc:205` stores `PhysicalObservation::Source(track)`. In `include/PhysicalObservation.hh:22`, that helper inspects the immediate creator process: no creator -> 0, `Scintillation` -> 1, `Cerenkov` -> 2, everything else -> 3. A conventional photon created by `OpWLS` or `OpWLS2` would thus be 3, irrespective of the original photon ancestor. The entire cell has zero source_type=3 hits and no loaded WLS absorption/emission/time properties. The Geant4 WLS implementation returns an effectively infinite interaction length when WLSABSLENGTH is absent (`G4OpWLS.cc:281`).

A photon created by scintillation from a secondary depositing particle remains source_type=1. Such a population can be displaced and late without being optical WLS. This is consistent with D1–D3, but identifying the parent species and deposition delay requires information absent from this tree. No definite parent-process attribution is claimed.

Figure: [creation locality and conditional profile](/home/rrios/exec46_20260915/decay_diagnostics/d3_creation_locality.pdf).

## D4 — Double-exponential tail fits

Binned Poisson fits use exact bin integrals with 0.1 ns bins:

```text
mu_bin = A * { f1 [exp(-(low-cut)/tau1) - exp(-(high-cut)/tau1)]
             + f2 [exp(-(low-cut)/tau2) - exp(-(high-cut)/tau2)] }
f2 = 1-f1; 0 < tau1 < tau2
```

Here f1 and f2 are integrated fractions **conditional on being above the lower fit cut**, extrapolated to infinity. They are not source fractions at emission. Parameter errors use an event-cluster sandwich covariance. Pearson chi2/ndf is reported as requested; sparse bins and intra-event correlations limit a literal asymptotic p-value interpretation.

| Window [ns] | tau1 [ns] | f1 | tau2 [ns] | f2 | chi2/ndf |
|---|---:|---:|---:|---:|---:|
| 6.3–42.0 | 2.296377 ± 0.003430 | 0.943636 ± 0.002698 | 4.591869 ± 0.059282 | 0.056364 ± 0.002698 | 566.637/353 = 1.605 |
| 6.3–25.2 | 2.219979 ± 0.011604 | 0.848976 ± 0.016503 | 3.543869 ± 0.084528 | 0.151024 ± 0.016503 | 240.085/185 = 1.298 |
| 10.5–42.0 | 2.372300 ± 0.007760 | 0.926362 ± 0.004913 | 5.408216 ± 0.123764 | 0.073638 ± 0.004913 | 341.742/311 = 1.099 |

The fast component does **not** recover 2.1 ns. Both lifetimes and weights depend on the chosen fit window. For the primary 6.3–42 ns fit, fixing tau1=2.1 ns increases the Poisson deviance by 1077.13 and gives Pearson chi2/ndf = 2173.37/354 = 6.140 (*unreliable*, >5). The free fit has chi2/ndf=1.605, below the inherited >5 flag threshold, but it is still an imperfect descriptive approximation with structured residuals.

The primary fit extrapolates to 131.29 photons above 42 ns; 319 are observed. The fit therefore does not close the extreme tail either. No component of this fit is promoted to a new material constant.

Figure: [double exponential and residuals](/home/rrios/exec46_20260915/decay_diagnostics/d4_double_exponential.pdf).

## Effective emission-pool scale and proposed gate — NOT APPLIED

The data do not support a single cut-independent exponential lifetime for the whole pool. Two quantities can be reported without conflating their meanings:

1. Conditional late-tail scale at 10.5 ns: **2.601323 ± 0.002957 ns**, a mean excess restricted to that tail.
2. Whole detected-scintillation-pool moment scale, with the recorded global origin t_ref=0: **tau_pool,moment = E[t_creation] = 3.053070 ± 0.000326 ns**. This uses all 60,368,164 detected scintillation photons. It includes finite rise, event production times, spatial composition and detection selection. It is not the material decay constant or a measurement of pure emission delay relative to each parent deposit.

The END pool was additionally measured without selecting first photons. With equal weight for each event and each of its two END faces, the measured moment scale is:

| Pool | Mean global creation time [ns] | Event SE [ns] |
|---|---:|---:|
| Left, equal event weight | 3.098377 | 0.001150 |
| Right, equal event weight | 3.097662 | 0.001146 |
| END, equal event and face weight | **3.098019** | **0.000814** |

The END sample contains 10,247,023 detected scintillation photons. Photon-weighted pooling instead gives 3.098263 ± 0.000802 ns. The two weighting definitions are retained in `emission_pool_by_face.csv`.

For the requested exponential napkin benchmark in Step 4, I propose the **measured END moment scale 3.098019 ± 0.000814 ns**, with the recorded global origin t_ref=0 and equal event/face weighting fixed explicitly. `tau_pool,moment/N` would be labeled an **exponential surrogate**, not a validated prediction of the actual first-photon distribution. It must not be called an intrinsic material lifetime or a pure emission delay: it includes the production-time origin and finite rise. This value is specific to EJ-200 at x=0; no transfer to another material or position is justified here.

A defensible quantitative benchmark for the actual composite pool is the empirical-CDF expression `E[min(t1,...,tN)] = integral_0^infinity [1-F_pool(t)]^N dt`, under an explicitly tested i.i.d. assumption. It preserves the finite rise and mixture. If a single number is desired at each N, define `tau_eff(N) = N * integral [1-F_pool(t)]^N dt`; in general it depends on N. No such Step 4 calculation was run here. There is no justified universal effective decay constant that both summarizes D2 and predicts the first photon.

Proposed replacement for the current global-tail compatibility clause, subject to approval:

- Keep all existing schema, clock, causality, path, speed and uniqueness gates.
- Require the effective MPT component inventory and finite-rise setting to match the intended material configuration. Never infer one component only from a file name.
- Use source_type=1 and abs(x_creation-gun_x)<0.1 mm as a declared local diagnostic proxy. Compare mean excess at cuts 3, 4 and 5 configured tau_d to the configured decay, with a declared 3% relative allowance plus three event-cluster SE. Require at least 200 contributing events at each cut; insufficient support is INCONCLUSIVE. This selection and rule are proposed now, informed by D3; they are not claimed to have been preregistered.
- Preserve and report the full-pool six-cut scan, spatial profile and measured pool-moment scale. A changing global mean excess is a population-composition result, not by itself a material-load failure. Do not tune cuts or replace the mixture with a fitted material lifetime.
- Report the exponential-surrogate benchmark alongside the empirical-pool benchmark in Step 4, with END/face/event weighting explicit. Use independently evaluated events for any predictive validation.

The proposed gate is not implemented. The original failed gate and its artifacts remain intact; resumption requires René's explicit approval.

## Artifacts and reproducibility

All output artifacts are under `/home/rrios/exec46_20260915/decay_diagnostics/`.

- `effective_mpt.json`, `effective_mpt.log`, `opsc-100.mac.txt`, `d1_provenance.json`: full D1 evidence.
- `d2_cut_stability.{pdf,root,csv,meta.json}` and `d2_paired_cut_increases.csv`.
- `d3_creation_locality.{pdf,root,csv,meta.json}`, `d3_creation_profile.csv`, `d3_correlations.csv`, `d3_local_tail_stability.csv`.
- `d4_double_exponential.{pdf,root,csv,meta.json}` and `d4_fixed_fast_component.csv`.
- `decay_sufficient_statistics.npz` with its source/config/code-linked cache metadata; `d1_d4_summary.json`; `emission_pool_by_face.csv` and its metadata.
- The ROOT sidecars contain histograms, profiles and fit bin predictions/parameters. No source ROOT was opened for writing.

No pushes, merges, deck edits, simulations, or full-grid resumption were performed.
