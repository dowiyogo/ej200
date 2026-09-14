# Production optical observations (EXEC40)

Worktree starts from main 420addf; the three new ntuples are present with
`EJ200_ENABLE_DIAGNOSTICS=OFF` and `ON`. The pre-existing `sipm_hits` columns,
detection code, random draws and escape guard are preserved. This claim is
tested against every generated event of archived EXEC29 D1, not inferred from
a code review. There is no authorization to rerun the 21-cell grid.

## Output contract

All tables join on zero-based `event_id`. Event-local track IDs identify optical
photons; ID sets are cleared at the beginning of each event and are private to
each worker. No track user-information ownership is taken.

* `event_observables`: one row per generated event, including zero-hit events.
  `edep_total_MeV` sums `GetTotalEnergyDeposit()` on all steps whose pre-step
  logical volume is `BarLV`, excluding sensor daughters. The separate
  `edep_nonionizing_MeV` is `GetNonIonizingEnergyDeposit()`;
  `edep_ionizing_MeV=total-nonionizing`. `edep_optical_MeV` identifies optical
  deposition included in the total. `produced_scint` and `produced_optical`
  count secondaries at creation on those steps, not detected photons or only
  tracks which survived to a boundary. L/R/T detected counts include zeroes.
* `first_bar_encounters`: one row per optical track's first physical encounter
  involving BarLV. Undefined, NotAtBoundary and StepTooSmall bookkeeping steps
  are excluded. A photon absorbed in bulk before reaching a surface has no row.
  Subsequent reflections never add rows. Pre/post physical volume names and
  copy numbers identify the ordered pair. `exiting_bar` specifies direction.
  `source`: 0 primary, 1 Scintillation, 2 Cerenkov, 3 other creator.
  Raw `boundary_status` uses this Geant4 build's enum. Portable `outcome` is
  0 unknown, 1 reflected, 2 transmitted, 3 absorbed, 4 detected; the raw enum
  retains the distinction between absorption and a missing refractive index.
  `cos_incidence` is pre-step momentum dotted with the oriented global normal.
  The normal points from the pre volume to the post volume: on exit use the
  pre solid's outward normal, on entry into a daughter use the negative of
  the daughter's outward normal. Transform using the corresponding touchable
  history. Do not use the bar's outer box normal on an embedded sensor face.
  The solid must classify the endpoint as `kSurface`; ±1e-5 mm probes verify
  inward/outward orientation. Unit length and signed cosine are checked in the
  ROOT regression. No absolute value or data-dependent sign correction occurs.
  `GetGlobalExitNormal()` is deliberately not called because it writes cached
  navigator state in Geant4 11.4; this uses the permitted equivalent solid normal.
* `sipm_event_counts`: one row per active sensor/event, including zeroes.
  EndTop70 has left global IDs 0–7, right 8–15, TOP 16–85. Sensor IDs derive
  from the configured border map, not a fixed active-channel count in C++.
  `incident` independently counts unique (track, sensor) accepted encounters:
  post-step is fGeomBoundary, post physical volume is a sensor, pre != post,
  and status is Detection, Absorption, or a supported transmission/refraction
  status including SameMaterial. Reflected photons are explicitly excluded.
  At the R=0 dielectric_metal sensitive surface an absorbed/detected photon
  terminates at the surface; requiring subsequent volume transport would
  incorrectly exclude those accepted photons. Transmission through other
  directed interfaces is counted and separately labeled.
  `detected` observes the existing SD callback independently. Unique detections
  are matched to incident identities at event end because SD invocation may
  occur before stepping. Unmatched or duplicate detections fail the observation
  regression; they never manufacture an incident denominator.
  Surface Detection, Absorption and transmission counters audit detection paths.
  `expected_surface_pde_sum` sums the actual directed border's EFFICIENCY at
  each accepted incident energy (Geant4 energy interpolation); missing EFFICIENCY
  increments `unknown_surface_pde`. This is separate from the historical
  emission-weighted reference and the wavelength-interpolated `sipm_hits.pde`.
  `surface_pde_variance_sum=sum(p*(1-p))` and the wavelength sum are diagnostics
  of the independently measured incident population, without additional RNG.

The detailed old all-bounce and terminal-fate censuses retain their diagnostic
flag. First-encounter records scale with produced photons; channel and energy
ledgers are small. Measure cost, and seek a decision before moving any fields.

## Reproduce validation without further transport

Use the existing analysis environment (NumPy 1.23.5, SciPy 1.13.1, uproot 5.6.9):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python \
  analysis/validation/check_exec40.py \
  --cell /home/rrios/exec40_20260913/cell_validated \
  --d1 /home/rrios/ej200_exec26_20260909/build_exec29_20260910/cells/D1 \
  --out /home/rrios/exec40_20260913/analysis
```

This test reconstructs D1's event counts from its ROOT and cross-checks its
archived CSV. It compares L/R/T counts exactly for all 2000 events and checks
the new event/channel ledgers independently against `sipm_hits`. Nonzero G-I.1
differences exit before V1/V2/V5. It checks all first-photon keys, finite signed
angles, normals, per-event energy data and every active channel, then evaluates
the fixed hypotheses. Exit 2 means invariance failed; exit 3 means observation
validation failed. Physical-hypothesis FAIL with valid observations exits 0:
that is a scientific finding, not a software test failure.

Both diagnostic builds must compile. Existing non-transport checks can run with:

```bash
ctest --test-dir /home/rrios/exec40_20260913/build_off \
  -R '^(sslg4_properties_check|readout_config_check)$' --output-on-failure
```

Do not run unfiltered CTest here: its smoke tests launch extra simulations.
The sole cell's exact command, binary and resource hashes, source patch, seeds,
worker count, wall time and output hash are in `cell_validated/invocation.meta.json`.
`run_exec40_cell.py` refuses an existing output directory and checks material
resources against D1 before creating the run. The original failed initialization
is preserved in `cell/`: missing local `sslg4` alias, fatal property lookup,
no event/track generated; it is not a second completed transport cell.

See EXEC40_PREREGISTRATION.md for thresholds, population definitions, resampling
and the first-encounter isotropic-flux applicability caveat. No result changes
parameters or selects a more favorable denominator. Acceptance remains false
until a separately authorized grid confirms the required tests.

## Validated cell outcome

G-I.1 PASS: exact event-by-event L/R/T equality; 397.1525 ± 2.4099272 pe/end.
V1 PASS: 1.00025691 ± 0.00015530. V2 FAIL against the registered flux
hypothesis: escape 0.25689824 ± 0.00007193, with all 41,695,454 first-encounter
normal and identity checks passing. G-I.2 fails incident/detection closure:
169,139 SD detections have no independently accepted incident identity.
Thus V5 is not validated (its numerical emitted-PDE comparison also fails).
The tool intentionally returns 3 on these data; do not weaken the check to
turn that into success. The matched-boundary population agrees with the
incident-spectrum PDE within its uncertainty, a separately labeled post-hoc
diagnostic that does not replace the primary numerator or criterion.

ROOT 2,380,261,978 bytes, 3.362× D1; wall 420.57 s. The report proposes
optional storage changes without applying them. No new cell can be used to
repair or validate missing histories under this one-cell authorization.
