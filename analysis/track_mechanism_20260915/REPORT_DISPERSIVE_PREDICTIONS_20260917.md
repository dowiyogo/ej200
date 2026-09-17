# EXEC_46 — wavelength-resolved optical predictions

2026-09-17. Implementation and focused validation only. No T2–T6 campaign
analysis, simulation, production ROOT write, MPT edit, or push was performed.
The preceding physics-analysis task was not resumed.

`analyze_step2.py`, `analyze_step3.py`, and `analyze_step4.py` now use
`dispersive_optics.py`. The default optical campaign is
`/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2`, overridable with
`EXEC46_CAMPAIGN_DIR`. Each cell resolves its own `sslg4` symlink and reads
its own RINDEX, absorption table, and gun macro. There are no ODR coefficients
or material-specific refractive-index constants in the implementation.
Existing analysis input/output routing is otherwise unchanged: derived inputs
must belong to the selected campaign. This patch does not regenerate them.

## Optical definitions and wavelength choices

- Cone and its axial edge: `wl_nm_created`, since emission fixes the cone.
  Beta comes from the cell's configured primary kinetic energy. These are
  primary-muon predictions; a source-type label alone does not establish
  charged-parent identity.
- Phase speed: `c/n(lambda_created)`. The requested beta=1 phase-edge expression
  `c*sqrt(1-1/n²)/n` is retained separately from the transport prediction.
- Group speed: numerical Geant4 `CalculateGROUPVEL` on the ascending-energy
  runtime mesh, interpolated in energy at `wl_nm` for transport. The existing
  transcription in `exec46_schema.py` was reused, including logarithmic
  differences, midpoint energies and the anomalous-dispersion guard. Effective
  group index is `c/GROUPVEL`; no analytic ODR derivative is substituted.
- Arrival-time edge speed: `vg(lambda_detected)*sin(theta_C(lambda_created,beta))`.
  With dispersion, substituting phase speed here would not reproduce F4.
- Step 2 final-angle threshold uses the **detected** wavelength. Step 3's
  cone/critical-identity diagnostic uses the **created** wavelength; both
  critical angles are retained as separate prediction fields.
- Step 4 penalty is `d/vg(lambda_detected) * [sec(alpha)-sec(alpha_edge)]`.
  Its q95-angle counterfactual is evaluated for each photon's distance and
  wavelengths, then summarized, rather than inserting a single material index.
  `q95_minus_edge_deg` is now the quantile of paired angle-minus-edge values.
- These are homogeneous-scintillator references, not reconstructions of every
  segment through other media or every reflection. Created/detected wavelength
  differences are retained. The spectral delay per length is
  `1/vg_detected - 1/vg_created`; the separate group-minus-phase delay per length
  is `1/vg_detected - n_detected/c`. Neither is declared identically zero for
  dispersive tables. A photon with unchanged wavelength can still have zero
  spectral difference and nonzero group-minus-phase delay.

`arccos(sin(theta_C)) = arcsin(1/n)` is asserted at beta=1 at every table node
and evaluated photon wavelength (absolute tolerance 2e-14 rad). Dispersion
replicates the identity at each wavelength. Finite-beta edges remain distinct.

Reports and CSV/ROOT sidecars carry median, q05 and q95 predictions, with
explicit population labels. Step 2 uses first/random photons by source; Step 3
uses first-by-source and first-overall photons; Step 4 retains its distance,
count-quantile and near-END selections. Step 3 all-hit transport means have no
retained wavelength distribution; their comparison is explicitly labelled as
using the **first-source reference**, not an invented all-hit spectrum.

Creation and detection clamp fractions are separate, with below-370 and
above-660 fractions in the numerical summaries. Outside-domain predictions are
labelled unmeasured. EJ-230 is labelled a constant-index model, not a measured
dispersion table with an interior trusted domain. No photon is removed because
of wavelength. Runtime paths, table hashes and constant/dispersive status are
recorded. Angular display coverage in Step 4 starts at 30 degrees instead of
38 to retain the dispersive edge; the 0.02-degree bin width is unchanged.

## Map of the supplied audit

Line numbers below are the **original audit locations**; named functions are
the stable locations in the revised scripts.

| Script / original lines | Change |
|---|---|
| Step 2: 689–693 | `main`: replace constant/identical assertions with per-cell informative runtime records. |
| Step 2: 694–696 | `main`: replace scalar angles with first/random, source-labelled photon distributions. |
| Step 2: 622–625, 627 | `build_report`: actual RINDEX/ABSLENGTH ranges and type; numerical group-speed and delay distributions replace the zero-correction assertion. |
| Step 2: 643 | `build_report`: wavelength-dependent angle distributions; `make_guiding_diagnostics` applies photon-specific final-angle thresholds. |
| Step 3: 49, 63–74 | Remove scalar index/angle/speed globals; shared runtime evaluator supplies all predictions. |
| Step 3: 332–336 | `verify_gun_and_index`: preserve gun checks; accept and record mixed constant/dispersive tables. |
| Step 3: 529, 537–545 | `render_report`: report photon-weighted angles and phase/group/edge speeds instead of unique scalar predictions. |
| Step 3: 577–578, 656–666 | `render_report`: fitted-speed comparisons use median/q05/q95 prediction distributions with explicit population labels. |
| Step 3: 736–738 | `main`: runtime provenance and prediction distributions replace scalar summary metadata. |
| Step 4: 38–51 | Remove OPSC-100 scalar index and derived globals; shared per-cell tables handle both constant and dispersive cases. |
| Step 4: 255, 258 | `cherenkov_nc_scan`: measured-minus-predicted velocities are paired photon by photon before averaging; prediction distributions are retained. |
| Step 4: 292, 296, 301 | `cherenkov_angle_window`: group-speed/edge-dependent penalties and q95-angle counterfactual distributions. |
| Step 4: 579 | `near_end_mixture`: photon-wise group-speed penalty; mixture weights and width estimators are unchanged. |
| Step 4: 642–643 | `save_bundle`: per-cell optical provenance replaces constant-index metadata. |
| Step 4: 721, 727–744 | `make_figures`: medians with q05–q95 bands replace horizontal/vertical scalar references; sidecars retain those values. |
| Step 4: 918–939 | `render_report`: wavelength-resolved edge/group predictions and updated penalty equation. |

The unchanged Step 3 builder omits wavelengths. `recover_selected_wavelengths`
therefore joins the already-selected photons to read-only `sipm_hits` on cell,
event, face and track ID. It verifies source and timestamps exactly and rejects
missing/duplicate matches. It preserves row order and all original columns;
it never recomputes the first photon or writes the derived ROOT.

Both builders have the same SHA-256 as at task entry. An AST comparison with
the saved pre-edit scripts confirms unchanged selection, transport fits,
clock/width calculations, order-statistics fits, bootstrap and integrity-gate
functions. In modified functions, changes are optical predictions, references,
provenance and their presentation. Pre-existing uncommitted work is excluded
from this commit.

## Validation

V1 passes with the actual EJ-230 table, reproducing the old values exactly at
the requested displayed precision:

| Quantity | New value |
|---|---:|
| theta_C(beta=1), degrees | 50.73475 |
| theta_critical, degrees | 39.26525 |
| c/n = vg, mm/ns | 189.742062 |
| beta=1 phase/transport edge, mm/ns | 146.902903 |
| finite-beta transport edge, mm/ns | 146.449825 |

V2 passes with the actual cell tables; these are **reference wavelengths**, not
photon-population averages:

| Material / reference wavelength | n | numerical group index | vg [mm/ns] | finite-beta transport edge [mm/ns] |
|---|---:|---:|---:|---:|
| EJ-200 / 420 nm | 1.626295995 | 1.818692353 | 164.839565946 | 129.629369302 |
| EJ-204 / 408 nm | 1.619785205 | 1.744067020 | 171.892739541 | 134.839616626 |

EJ-200 finite-beta edge: 37.458195763 degrees at 370 nm and 40.174267051 at
660 nm, span 2.716071288 degrees. The configured beta is 0.995423527580.
The small group-index differences from the supplied analytic/rounded F4
references are retained: tolerances are 1e-6 for n, 4e-6 for group index,
0.001 mm/ns for speed and 0.001 degrees for edge angle. These tolerate mesh
differentiation and reference rounding; they do not change the runtime table.

V3 exercises the actual Step 3 diagnostic with two photons at the same measured
angle but different wavelengths straddling their own critical angles. The old
scalar code returns fraction-below-critical **1.0**, failing the expected
**0.5**. The revised code returns **0.5** and nonzero angle-quantile width.
This failure was reproduced against the saved pre-edit Step 3 script, not
inferred from inspection.

Nine focused tests pass, covering V1–V3, creation/detection wavelength roles,
per-cell runtime resolution within one material, identity/clamps, Step 2
thresholds and sidecars, Step 3 join integrity, and Step 4 temporal penalties.
Only small synthetic photon fixtures and runtime configuration files were read;
the full campaign was not analyzed. Syntax compilation and `git diff --check`
also pass.
The same nine tests also pass against an isolated export of the Git index,
independent of the pre-existing uncommitted changes.

```bash
PYTHONPATH=analysis/track_mechanism_20260915 python3 -m unittest discover \
  -s analysis/track_mechanism_20260915 -p test_dispersive_optics.py -v
```

Verification snapshots and numerical runtime records are preserved in
`/home/rrios/exec46_dispersion_edit_20260917/`. The three runtime RINDEX SHA-256
values (identical within each material's seven cells) are:

```text
EJ-200 15d1f8cf5a62effd9a0f2f9bd1edaeb5164b94c2035a11cfe6a878d51680f6ed
EJ-204 f1c77ee162c767cd23be608e0f8081174263352e40a0e1f57b8ca97ccef8636f
EJ-230 d3042cd86cab09ba7e8f08a97b94d2435aa2858c99aed7341d012362b4d555f6
```

Push command, printed only and **not executed**:

```bash
git push origin diag/exec46-track-mechanism-20260915
```
