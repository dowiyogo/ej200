# EXEC_26 Phase 1b — observational optical-boundary census
2026-09-09 · English report · base main `8349041140958226a0ac1cb3bb3e30aff2303435` · instrumented worktree `diag/reflector-ab-20260909` at `0336ba9dbb3fd16525eb4baec0928e3cf0955632`.
## Verdict
**The measured air→wrap outcomes match the preregistered “reflector inerte” signature in all three forward panels.** In run **R1** (500 generated events; seeds **26092601, 8349041**; EJ-204, EndTop, 70 TOP, one worker), the forward interface recorded **13.371145% FresnelReflection, 84.620018% FresnelRefraction, 2.008837% Absorption, and 0% TotalInternalReflection** over **3,623,938 eligible encounters**. No forward reflection variant exceeded 90%; the alternative/refutation scenario was not observed. The measured boundary behavior does not support interpreting this configuration as a 98% reflecting interface in this run. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.meta.json).
“Inert” here is the label of the user’s preregistered status-pattern scenario. The measured reflector does return some photons; this is **not** a claim of zero reflection or a measurement of less-than-5% endpoint yield sensitivity. No R=0/R=1 scan or factory A/B was performed, and no conclusion is extended to an unrun geometry or material. Interface refraction was counted; subsequent bulk absorption inside the wrap was not instrumented. All reflector-behavior statements in this report refer to R1 and its recorded metrics.
The same R1 data give **82.304 ± 1.023 detected pe/end per generated event** (mean ± event-level standard error), compared with the deck’s 514.9. `sipm_hits` contains global IDs from **0 through 85**, with hits in **70 distinct TOP channels**. The authorized run is EJ-204/70 TOP, while the deck’s design frame names EJ-230/20 TOP; this is a reported mismatch, not a controlled reproduction with matching settings. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.meta.json).
**Phase 1b is stopped at the required gate. Phase 2 was not started.** No materials, optical surfaces, geometry, gun, transport, existing kill guard or existing random draw was changed. Main remains clean at the original SHA.
## R1: exact run provenance
- Source commit: `0336ba9dbb3fd16525eb4baec0928e3cf0955632`; original physics base: `8349041140958226a0ac1cb3bb3e30aff2303435`.
- Geant4 runtime: **geant4-11-04 [MT], 5-December-2025**, version **11.4.0**; install `/home/tdship/opt/geant4-v11.4.0-install`. Linked libraries and binary hash are preserved.
- Seeds supplied: **26092601 8349041**; engine **MixMaxRng**; `/run/numberOfThreads 1`, observed worker `G4WT0`. Master and worker RNG states are saved before/after the run.
- N: **500 generated events**, all 500 have at least one recorded detection. No zero-hit event was omitted from the mean denominator.
- Material: **OPSC-101 / EJ-204**, readout **EndTop**, 16 END SiPMs and 70 TOP. Gun: mu−, 1 GeV, vertical, x=0 mm; jitter macro value 0 ns, existing jitter implementation unchanged.
- TOP positions: for local i=0…34, x=−692+20i mm; for i=35…69, x=12+20(i−35) mm; y=29.75 mm, z=0. Global/copy IDs 16…85. Coordinates are in every metadata sidecar.
- Start UTC: **2026-09-09T16:40:17.179612+00:00**. End UTC: **2026-09-09T16:41:50.539700+00:00**. Exit code: **0**. Elapsed wall time from ledger: **93.360 s**.
Executed command [C025](#c025), with working directory `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready`:

```bash
/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/run500.mac
```

Full successful log: [command_025.log](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_025.log). Raw hits: [photon_hits_run000.root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/photon_hits_run000.root). Execution macro: [run500.mac](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/run500.mac).
RNG states: [rng_run0_thread-1_begin.rndm](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/rng_run0_thread-1_begin.rndm), [rng_run0_thread-1_end.rndm](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/rng_run0_thread-1_end.rndm), [rng_run0_thread0_begin.rndm](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/rng_run0_thread0_begin.rndm), [rng_run0_thread0_end.rndm](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/rng_run0_thread0_end.rndm). Their hashes are embedded in all derived-table metadata. Recording engine state serializes the engine; it adds no random draw.
## Preregistered prediction and denominator
The following user prediction is quoted literally, before filling its measurement column. The preregistration and the runtime-directory amendment both predate R1: [preregister.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/preregister.json); [preregister_ready.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/preregister_ready.json). No numerical tolerance was invented after the data for the approximate signs.

```text
| Estado | Predicción (reflector inerte) | Medido |
|---|---|---|
| FresnelReflection | 6–20 % | |
| FresnelRefraction | 75–90 % | |
| Absorption | ≈ 2 % | |
| TotalInternalReflection | ≈ 0 % | |

Escenario alternativo: si FresnelReflection o alguna variante de
reflexión supera el 90 %, la hipótesis de reflector inerte queda
refutada y así debe declararse.
```

**Eligible encounter denominator, fixed before analysis:** recorded optical-photon steps with post-step status `fGeomBoundary`, matching **AirGap<stem>PV → Vikuiti<same stem>PV**, for stems YMinus, ZPlus and ZMinus. Directions and panel copy numbers are preserved. Exclude `Undefined`, `NotAtBoundary`, `SameMaterial`, and `StepTooSmall`; retain their counts separately. Counts are interface encounters, not a deduplicated track sample. The conditional hook does not count non-geometrical steps; zeros in NotAtBoundary do not establish absence of such steps in transport.
For R1 forward direction: **3,623,951 raw records − 13 StepTooSmall = 3,623,938 eligible encounters**. Other excluded status counts are zero. Percent = 100 × status_count / 3,623,938. No independent-photon error model is applied to these correlated encounter counts. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.meta.json).
| State | Preregistered prediction | Measured count | Measured percentage |
|---|---|---:|---:|
| FresnelReflection | 6–20 % | 484,562 | 13.371145% |
| FresnelRefraction | 75–90 % | 3,066,577 | 84.620018% |
| Absorption | ≈ 2 % | 72,799 | 2.008837% |
| TotalInternalReflection | ≈ 0 % | 0 | 0.000000% |
R1 table sidecars: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.meta.json). Absorption differs from exactly 2% by **+0.008837 percentage points**; measured TIR is exactly zero. Total forward reflection across all variants equals FresnelReflection, **484,562 encounters (13.371145%)**; all other reflection variants are zero. These observations match the stated ranges without broadening them. [C025](#c025).
### Panel separation and reverse direction
| Direction | Panel | Raw | Eligible | StepTooSmall | FresnelReflection % | FresnelRefraction % | Absorption % | TIR % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| air_to_wrap | YMinus | 1,008,022 | 1,008,022 | 0 | 10.053253 | 87.954926 | 1.991822 | 0.000000 |
| air_to_wrap | ZPlus | 1,305,480 | 1,305,467 | 13 | 14.670306 | 83.308042 | 2.021652 | 0.000000 |
| air_to_wrap | ZMinus | 1,310,449 | 1,310,449 | 0 | 14.629108 | 83.361733 | 2.009159 | 0.000000 |
| wrap_to_air | YMinus | 101,339 | 0 | 101,339 | N/A | N/A | N/A | N/A |
| wrap_to_air | ZPlus | 191,530 | 14 | 191,516 | 0.000000 | 0.000000 | 7.142857 | 92.857143 |
| wrap_to_air | ZMinus | 191,707 | 0 | 191,707 | N/A | N/A | N/A | N/A |
R1 panel table sidecars: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.meta.json). Each forward panel independently falls within 6–20% reflection and 75–90% refraction, with absorption approximately 2% and TIR zero. The aggregate therefore does not hide a failing forward panel. [C025](#c025).
Reverse wrap→air is a different population: **484,576 raw records, of which 484,562 are StepTooSmall**, leaving **14 eligible encounters: 13 TIR and 1 Absorption**. Its 92.857143% TIR fraction is based on those 14 **reverse** encounters and is **not** the forward >90% alternative criterion. YMinus and ZMinus have no eligible reverse observations; N/A is retained rather than displaying a fictitious zero percentage. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.meta.json).
All 43 status names, including Undefined, NotAtBoundary, SameMaterial and StepTooSmall, are exported for every observed PV/copy pair, even at zero count. The complete raw counter also has its own trio: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.meta.json). StepTooSmall is never counted as physical reflection.
## Same-run deck checks
The following table uses the same R1 events, without changing material, layout, gun position or analysis selection to approach a deck value. [C025](#c025).
| Quantity | R1 measurement | Deck reference | Qualification |
|---|---|---|---|
| Observed global_id range | 0–85 | Not supplied | Actual sipm_hits rows, not declared geometry |
| TOP channels with at least one hit | 70 | 20 TOP in design | Unchanged base geometry has 70 TOP |
| Mean detected N_pe/end/generated event | 82.304 ± 1.023 (event SEM) | 514.9 | 82,304 END hits / (2 × 500); EJ-204, 70 TOP versus deck EJ-230, 20 TOP |
| TOP Absorption+Detection encounters / generated scintillation photons | 13.342719 % | 26.6 % | 1,336,894 / 10,019,652; encounter/scintillation normalization, not total-optical unique photon fraction; deck 26.6% is END depletion versus END-only |
| Detected TOP hits / generated scintillation photons | 8.301935 % | Not equivalent to 26.6 % | 831,825 / 10,019,652; distinct from absorption+Detection encounters |
| Requested unique TOP-intercepted / all generated optical photons | UNVERIFIED | Not the definition in the deck | Census has no track/creator identifiers; existing generation counter counts scintillation only |
Table sidecars: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.meta.json). The channel-level hit distribution is preserved in [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.meta.json); all 500 event counts and yields in [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.meta.json); underlying numerical totals in [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.meta.json).
**Mean END calculation:** `(40,916 + 41,388) / (2 × 500) = 82.304`. The 1.022890 pe/end standard error is the sample standard deviation of the 500 per-event END means divided by sqrt(500), not a Poisson assumption about summed hits. The difference from 514.9 is **−432.596 pe/end (−84.015537%)**. No event selection, efficiency rescaling or material correction was made. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.meta.json).
### What the TOP fraction does and does not measure
R1’s existing generation counter records **10,019,652 scintillation photons**. The new counter records **1,336,894 eligible TOP interface encounters**, all from BarPV→TopSiPMPV, comprising **517,662 Absorption + 819,232 Detection**. Their ratio to that generation counter is **13.342719%**. The separate TTree TOP-hit count is **831,825**, giving **8.301935%** relative to the same scintillation counter. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.meta.json).
**The exact fraction of unique TOP-intercepted photons over all generated optical photons is UNVERIFIED from the authorized schema.** The existing counter is scintillation-only; the PV-pair/status map does not store track or creator IDs. It cannot certify an all-optical production denominator or identify/deduplicate individual photons. The two reported normalized quantities above are explicitly identified observational proxies, not silently renamed as the unavailable all-optical fraction. No additional instrumentation or run was introduced to fill that gap.
A further measurable mismatch remains: TTree TOP hits exceed eligible TOP Detection statuses by **12,593** (**831,825 − 819,232**). The run-summary TOP count agrees exactly with the TTree, but not with that status count. All nonzero post-TopSiPMPV counter rows were checked; they are from BarPV and contain only Absorption or Detection, so selecting a different recorded incoming volume does not close the difference. Its causal explanation is UNVERIFIED; neither count was adjusted. [C025](#c025); [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.meta.json).
The deck itself labels **−26.6% “vs END-only”** beside **701.3 END-only and 514.9 hybrid** in `talk_v6.tex:1021–1023`. Its design frame at lines 766–787 explicitly names EJ-230 and N_TOP=20. Thus 26.6% describes relative END-yield depletion, rather than TOP encounters or detected TOP hits divided by generated photons. This single EndTop run has no matched END-only arm, so it does not measure that depletion ratio. The requested numerical comparison is reported; equality would not be a valid reproduction test across different definitions, material and TOP layout. No settings were changed in response. [C028](#c028).
## Instrumentation and preservation checks
The existing rollback tag was verified at the base SHA and was not recreated. All source edits were confined to `/home/rrios/ej200_exec26_20260909`. Exactly the two requested commits were created, in order:
1. `3f0808d130d681a8546878a165a097aaf059ef67` — `feat(diag): BoundaryCensus for optical boundary status by PV pair (EXEC_26)`. Adds the singleton/header, all-status CSV export, and RunAction reset/write/RNG-snapshot lifecycle.
2. `0336ba9dbb3fd16525eb4baec0928e3cf0955632` — `feat(diag): wire BoundaryCensus into SteppingAction — SteppingAction.cc:57 before→after`. Adds only includes and the counter hook before every existing filter. Original function line 51 becomes line 56; the inserted hook comment is line 57.
The existing project already uses namespace `BoundaryCensus` for legacy accessor functions. The new singleton is `BoundaryCensus::Census`, accessed exactly as `BoundaryCensus::Instance().Record(...)`; this preserves the legacy namespace/API. Its `std::map<BoundaryKey,G4long>` is protected by a mutex. Reset and CSV export occur once on the master. The hook finds OpBoundary through the process list once per thread and never uses GetProcessDefinedStep as its sole lookup.
Validation before R1: the simulator built successfully in a fresh Release build directory, with no compiler warnings/errors. The standalone counter test verified 10,000 concurrent increments, reset, CSV quoting, separate directions/copy numbers, all 43 enum names plus unknown-status retention. Its first standalone link omitted libG4clhep and failed; adding the correct library to the **test link command** resolved it without changing simulation source. A subsequent preparation check initially failed because that test output did not yet exist; it then passed after the corrected test ran. All journaled failures are retained.
The protected-source check proved that **no pre-existing line was deleted or modified in RunAction.cc or SteppingAction.cc**: these files contain additions only. DetectorConstruction.cc, Materials.cc, the detector header, PrimaryGeneratorAction.cc, SiPMSD.cc, main.cc and the tracked macro are unchanged from the base. Final diff is exactly four instrumentation files, **201 inserted lines**. Both worktrees are clean in porcelain status; generated run/build evidence is ignored by existing ignore rules. Main remains `8349041`. No Phase 2 factory edit, push, branch/tag deletion, reset, rebase or stash drop was executed. [C030](#c030).
## Operational exceptions and limitations
- **Macro reconciliation:** the existing `macros/endtop_smoke_center.mac` requests 50 events. It was kept byte-for-byte unchanged. The new run-only copy prepends exactly the two approved thread/seed lines and changes only beamOn 50→500 to satisfy the requested N. Source and execution-copy hashes are in metadata. No tracked gun or run macro was edited.
- **One failed startup, zero generated events:** the initial new output directory did not contain the `sslg4/...` relative runtime data expected by the unchanged factory; initialization reported a missing macro and then a fatal missing SCINTILLATIONYIELD property, exit 250. It produced no event sample. All failed-start artifacts were preserved. A second, distinct output directory links to CMake’s existing byte-identical SSLG4 runtime-data tree. The macro, binary, seeds and source were unchanged. R1 is the sole completed 500-event run.
- **Existing material warning:** R1 emits `mat031`, mass fractions summing to 0.996667 rather than 1 for opsc-101. This was not corrected; the result is conditional on the unchanged approved model. R1 subsequently reports EJ-204 yield 10400 ph/MeV, rise 0.7 ns, decay 1.8 ns and attenuation 160 cm.
- **Legacy banner:** RunAction’s old skin-surface summary prints reflector R=0 / OPEN/UNDEFINED. The protected source still binds the original border factory. That banner is not used as the measured interface reflectivity; the new directional status census is the measurement.
- Counts are encounters. The requested aggregate counter has no event/track/creator labels, so no event-bootstrap boundary-fraction error bars or total-optical unique interception estimate can be recovered retrospectively.
- No figure was needed. Every derived result table, including the full raw exported census table, has all three sidecars. ROOT and CSV values were round-trip compared, including zeros and NaN for empty denominators; metadata fields and hashes passed validation.
- The earlier <5% R-endpoint criterion and the one-line factory A/B remain unmeasured. They require the next approval. No post-result source edits or additional simulations were performed.
## Sidecar inventory
Ten validated table sets; each link group below includes `.root`, `.csv`, and `.meta.json`:
- `directional_states`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/directional_states.meta.json)
- `panel_states`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/panel_states.meta.json)
- `prediction_comparison`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.meta.json)
- `event_yields`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/event_yields.meta.json)
- `observed_channels`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/observed_channels.meta.json)
- `deck_metrics`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_metrics.meta.json)
- `top_boundary_states`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/top_boundary_states.meta.json)
- `deck_reproduction`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/deck_reproduction.meta.json)
- `hit_census_reconciliation`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/hit_census_reconciliation.meta.json)
- `boundary_census_run0`: [root](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.root) · [csv](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.csv) · [meta.json](ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.meta.json)
Every metadata file includes source/base commits, Geant4 identity/version/prefix, binary hash, seed pair, N, material, worker count, exact TOP layout, exact run command and cwd, run timestamps/exit code, RNG-state hashes, input checksums and analysis provenance. ROOT tables are TTrees named `table`.
## Appendix — exact commands, timestamps and exit codes
Machine-readable ledger: [commands.jsonl](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/commands.jsonl). Each command below was executed through the preserved [run_command.py](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/run_command.py) wrapper, which uses `GIT_OPTIONAL_LOCKS=0`; stdout and stderr are combined in its linked complete log. Working directories are explicit. Collection/edit/analysis scripts are retained in the same audit directory. R1 is [C025](#c025).
<a id="c001"></a>

**C001** — start `2026-09-09T16:34:52.738166+00:00`; end `2026-09-09T16:34:52.740337+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_001.log).
```bash
git -C /home/rrios/ej200 worktree list --porcelain
```

<a id="c002"></a>

**C002** — start `2026-09-09T16:34:52.755922+00:00`; end `2026-09-09T16:34:52.757999+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_002.log).
```bash
git -C /home/rrios/ej200 rev-parse 'pre-exec26-reflector-ab-20260909^{commit}'
```

<a id="c003"></a>

**C003** — start `2026-09-09T16:34:52.773534+00:00`; end `2026-09-09T16:34:52.778354+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_003.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 status --porcelain=v2 -b
```

<a id="c004"></a>

**C004** — start `2026-09-09T16:34:52.794163+00:00`; end `2026-09-09T16:34:52.798514+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_004.log).
```bash
rg -n 'namespace BoundaryCensus|UserSteppingAction|BeginOfRunAction|EndOfRunAction|#include' /home/rrios/ej200_exec26_20260909/src/SteppingAction.cc /home/rrios/ej200_exec26_20260909/src/RunAction.cc
```

<a id="c005"></a>

**C005** — start `2026-09-09T16:36:05.810913+00:00`; end `2026-09-09T16:36:05.822625+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_005.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/instrument_counter.py
```

<a id="c006"></a>

**C006** — start `2026-09-09T16:36:05.837368+00:00`; end `2026-09-09T16:36:05.841295+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_006.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 diff --check
```

<a id="c007"></a>

**C007** — start `2026-09-09T16:36:05.855776+00:00`; end `2026-09-09T16:36:05.859389+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_007.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 add include/BoundaryCensus.hh src/BoundaryCensus.cc src/RunAction.cc
```

<a id="c008"></a>

**C008** — start `2026-09-09T16:36:05.874182+00:00`; end `2026-09-09T16:36:05.881792+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_008.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 commit -m 'feat(diag): BoundaryCensus for optical boundary status by PV pair (EXEC_26)'
```

<a id="c009"></a>

**C009** — start `2026-09-09T16:36:44.913467+00:00`; end `2026-09-09T16:36:44.915820+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_009.log).
```bash
rg -n 'UserSteppingAction|#include|GetDefinition|fGeomBoundary' /home/rrios/ej200_exec26_20260909/src/SteppingAction.cc
```

<a id="c010"></a>

**C010** — start `2026-09-09T16:36:44.932108+00:00`; end `2026-09-09T16:36:44.943525+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_010.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/wire_counter.py
```

<a id="c011"></a>

**C011** — start `2026-09-09T16:36:44.959095+00:00`; end `2026-09-09T16:36:44.963464+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_011.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 diff --check
```

<a id="c012"></a>

**C012** — start `2026-09-09T16:36:44.978317+00:00`; end `2026-09-09T16:36:44.981399+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_012.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 add src/SteppingAction.cc
```

<a id="c013"></a>

**C013** — start `2026-09-09T16:36:44.995476+00:00`; end `2026-09-09T16:36:45.002693+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_013.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 commit -m 'feat(diag): wire BoundaryCensus into SteppingAction — SteppingAction.cc:57 before→after'
```

<a id="c014"></a>

**C014** — start `2026-09-09T16:36:45.018493+00:00`; end `2026-09-09T16:36:46.045862+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_014.log).
```bash
cmake -S /home/rrios/ej200_exec26_20260909 -B /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909 -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

<a id="c015"></a>

**C015** — start `2026-09-09T16:36:46.061446+00:00`; end `2026-09-09T16:36:48.975353+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_015.log).
```bash
cmake --build /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909 --target ej200_bar_sim -j 8
```

<a id="c016"></a>

**C016** — start `2026-09-09T16:38:03.180933+00:00`; end `2026-09-09T16:38:04.323599+00:00`; exit **1**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_016.log).
```bash
c++ -std=c++17 -pthread -I/home/rrios/ej200_exec26_20260909/include -I/home/tdship/opt/geant4-v11.4.0-install/include/Geant4 /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/test_counter.cc /home/rrios/ej200_exec26_20260909/src/BoundaryCensus.cc -o /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/test_counter
```

<a id="c017"></a>

**C017** — start `2026-09-09T16:38:40.753742+00:00`; end `2026-09-09T16:38:40.771000+00:00`; exit **1**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_017.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/verify_prepare.py
```

<a id="c018"></a>

**C018** — start `2026-09-09T16:39:03.154066+00:00`; end `2026-09-09T16:39:04.283107+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_018.log).
```bash
c++ -std=c++17 -pthread -I/home/rrios/ej200_exec26_20260909/include -I/home/tdship/opt/geant4-v11.4.0-install/include/Geant4 /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/test_counter.cc /home/rrios/ej200_exec26_20260909/src/BoundaryCensus.cc -L/home/tdship/opt/geant4-v11.4.0-install/lib64 -Wl,-rpath,/home/tdship/opt/geant4-v11.4.0-install/lib64 -lG4clhep -o /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/test_counter
```

<a id="c019"></a>

**C019** — start `2026-09-09T16:39:22.132237+00:00`; end `2026-09-09T16:39:22.143559+00:00`; exit **0**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_019.log).
```bash
/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/test_counter
```

<a id="c020"></a>

**C020** — start `2026-09-09T16:39:22.158120+00:00`; end `2026-09-09T16:39:22.204000+00:00`; exit **0**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_020.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/verify_prepare.py
```

<a id="c021"></a>

**C021** — start `2026-09-09T16:39:22.220594+00:00`; end `2026-09-09T16:39:22.222490+00:00`; exit **0**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_021.log).
```bash
sha256sum /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/ej200_bar_sim /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/preregister.json
```

<a id="c022"></a>

**C022** — start `2026-09-09T16:39:22.238420+00:00`; end `2026-09-09T16:39:22.244759+00:00`; exit **0**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_022.log).
```bash
ldd /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/ej200_bar_sim
```

<a id="c023"></a>

**C023** — start `2026-09-09T16:39:33.436596+00:00`; end `2026-09-09T16:39:34.107310+00:00`; exit **-6**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_023.log).
```bash
/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500/run500.mac
```

<a id="c024"></a>

**C024** — start `2026-09-09T16:40:16.971096+00:00`; end `2026-09-09T16:40:17.006536+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_024.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/stage_runtime.py
```

<a id="c025"></a>

**C025** — start `2026-09-09T16:40:17.179612+00:00`; end `2026-09-09T16:41:50.539700+00:00`; exit **0**. Cwd: `/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_025.log).
```bash
/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/ej200_bar_sim -m /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/run500.mac
```

<a id="c026"></a>

**C026** — start `2026-09-09T16:43:37.700019+00:00`; end `2026-09-09T16:43:41.274323+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_026.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/analyze_phase1b.py
```

<a id="c027"></a>

**C027** — start `2026-09-09T19:37:10.351469+00:00`; end `2026-09-09T19:37:10.354302+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_027.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 show 8349041:src/SiPMSD.cc
```

<a id="c028"></a>

**C028** — start `2026-09-09T19:37:10.370065+00:00`; end `2026-09-09T19:37:10.373834+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_028.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 grep -n -E '514\.9|26\.6|TOP intercept|701\.3|Design Decision|EJ-230|EJ-204' 8349041 -- presentations/v6/talk_v6.tex
```

<a id="c029"></a>

**C029** — start `2026-09-09T19:39:11.500898+00:00`; end `2026-09-09T19:39:11.744886+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_029.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/complete_tables.py
```

<a id="c030"></a>

**C030** — start `2026-09-09T19:39:11.761635+00:00`; end `2026-09-09T19:39:11.766271+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_030.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 diff 8349041140958226a0ac1cb3bb3e30aff2303435 --stat
```

<a id="c031"></a>

**C031** — start `2026-09-09T19:39:11.782508+00:00`; end `2026-09-09T19:39:11.787369+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_031.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 status --porcelain=v2 -b
```

<a id="c032"></a>

**C032** — start `2026-09-09T19:39:11.804598+00:00`; end `2026-09-09T19:39:11.809659+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_032.log).
```bash
git -C /home/rrios/ej200 status --porcelain=v2 -b
```

<a id="c033"></a>

**C033** — start `2026-09-09T19:39:11.825611+00:00`; end `2026-09-09T19:39:11.828279+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_033.log).
```bash
git -C /home/rrios/ej200_exec26_20260909 log --format=fuller 8349041140958226a0ac1cb3bb3e30aff2303435..HEAD
```

<a id="c034"></a>

**C034** — start `2026-09-09T19:43:00.399818+00:00`; end `2026-09-09T19:43:00.416325+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_034.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/write_report.py
```

<a id="c035"></a>

**C035** — start `2026-09-09T19:43:48.078440+00:00`; end `2026-09-09T19:43:48.115093+00:00`; exit **0**. Cwd: `/home/rrios/ej200`. [Complete output](ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/command_035.log).
```bash
python /home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/audit/validate_final.py
```

One helper-launch failure was printed by the execution tool before the wrapper could append a ledger record: the attempted standalone `audit/test_counter` executable did not exist after its initial failed link. The wrapper exited 1 with FileNotFoundError. Its exact UTC timestamp is **UNVERIFIED**; it occurred between the failed standalone-link record and the failed preparation-check record. It did not execute a simulation. All consequential source edits, build, completed tests, startup failure, R1, analysis and final repository checks are present in the ledger.
## Stop condition
The forward data clearly match the preregistered first scenario, with no hidden failing panel. This report ends Phase 1b. **Await a new explicit approval before Phase 2.** The scope remains observational; no proposed destructive action is needed or executed.
