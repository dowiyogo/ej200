# Drift Classification

Command basis:

```bash
for branch in feature/edge-scan-and-readout-grouping feature/sipm-electronics-response; do
  for file in src/Materials.cc src/SiPMSD.cc src/DetectorConstruction.cc \
              src/PrimaryGeneratorAction.cc include/Materials.hh include/SiPMSD.hh \
              include/DetectorConstruction.hh; do
    diff_lines=$(git diff main..$branch -- "$file" | wc -l)
    if [ "$diff_lines" -gt 0 ]; then
      echo "=== $branch :: $file ($diff_lines lines diff) ==="
      git diff main..$branch -- "$file"
    fi
  done
done
```

Summary of diff sizes:

| Branch | File | Diff size | Classification |
|---|---|---:|---|
| edge | `src/Materials.cc` | 63 lines | A: stale baseline loss (`CreateBarSkinReflector` absent), not a feature physics change |
| edge | `src/DetectorConstruction.cc` | 283 lines | B/C: intentional edge-wrap feature, but conflicts with current baseline reflector |
| edge | `include/Materials.hh` | 15 lines | A: stale baseline loss (`CreateBarSkinReflector` declaration absent) |
| edge | `include/DetectorConstruction.hh` | 75 lines | B/C: intentional edge-wrap state, but stale hierarchy comments/members |
| electronics | `src/Materials.cc` | 213 lines | B plus A: EJ-230 implementation useful, but default material drift must be reverted |
| electronics | `src/SiPMSD.cc` | 315 lines | B plus A: electronics/sensor feature useful, but defaults and timing formula drift |
| electronics | `src/DetectorConstruction.cc` | 329 lines | B plus A: material selector useful, but stale wrap geometry and EJ230 default |
| electronics | `include/Materials.hh` | 38 lines | B plus A: EJ-230 API useful later, missing current baseline reflector naming |
| electronics | `include/SiPMSD.hh` | 78 lines | B plus A: split fields useful, defaults/canonical names wrong |
| electronics | `include/DetectorConstruction.hh` | 82 lines | B plus A: material selector useful, default and wrap hierarchy stale |

## Detailed Classification

### `feature/edge-scan-and-readout-grouping`

#### `src/Materials.cc` / `include/Materials.hh`

Classification: **A) accidental/stale physical drift**.

This branch lacks the current `main` reflector surface factory
`CreateBarSkinReflector()` and its R=0.98 explicit-panel surface. That is not
an edge-scan feature; it is simply because the branch predates the reflector
merge. On rebase, `main` should win for baseline optical material/surface
definitions.

Keep from feature: `CreateMylar()` and `CreateBlackTape()` remain useful if
edge-wrap studies are reimplemented as optional local perturbations.

#### `src/DetectorConstruction.cc` / `include/DetectorConstruction.hh`

Classification: **B/C mixed**.

Intentional feature behavior:

- `/det/edgeWrap mylar|air|black`
- segmented 50 mm edge caps
- top-readout placement logic and display support

Stale baseline drift:

- active `WrapCenterPV`, `WrapCapLeftPV`, `WrapCapRightPV`
- SiPM border surfaces registered from wrap volumes to SiPMs
- bar sub-volumes shrunk by `kMylarThick` in earlier geometry lineage
- no final explicit reflector panel baseline

Resolution rule: after rebase onto `fix/physics-baseline`, baseline geometry
must win by default. If edge-wrap variation is still needed, it should be
reintroduced as an optional study that does not alter the canonical default
EJ-200 reflector model.

#### `src/PrimaryGeneratorAction.cc`

No physics drift detected.

#### `src/RunAction.cc`

No ntuple drift detected relative to canonical 12 branches. It is merely stale
with respect to `main` diagnostics.

### `feature/sipm-electronics-response`

#### `src/Materials.cc` / `include/Materials.hh`

Classification: **B plus A**.

Useful/intended:

- `CreateEJ230()` values match `EJ228_EJ230.pdf`: rise 0.5 ns, decay 1.5 ns,
  yield 9700 ph/MeV, attenuation 120 cm, peak 391 nm.
- A future material selector is scientifically useful.

Accidental/out-of-scope drift:

- EJ-200 still lacks rise time 0.9 ns.
- The branch sets default scintillator to EJ-230 in detector state and macros.
  Current baseline must remain EJ-200.
- `CreateEJ230()` belongs in a future `feature/scintillator-materials-selection`
  or should be kept dormant with default EJ-200.

#### `src/SiPMSD.cc` / `include/SiPMSD.hh`

Classification: **B plus A**.

Useful/intended:

- Adds sensor-model selection machinery.
- Adds alternative PDE curves for Broadcom AFBR-S4N66P024M and Hamamatsu
  S14160 variants.
- Adds extra ntuple branches (`time_raw_ns`, `sptr_jitter_ns`,
  sensor metadata), which are acceptable as additive branch-name extensions.

Accidental/canonical drift:

- Default sensor is Broadcom, not S13360-6050PE prototype.
- S13360 is named `S13360-6025 (legacy)` rather than S13360-6050PE.
- Hamamatsu SPTR default is 150 ps, not 106 ps.
- Electronics jitter default is 30 ps, not 10 ps.
- `fElectronicsSigma` is stored in the ntuple but is not used in the timestamp
  smearing formula; only `fSPTRSigma` is applied.
- UI name `/sipm/electronicsSigma` differs from requested canonical
  `/sipm/electronicsJitter`.
- No legacy `/sipm/jitterSigma` alias is present.
- Snapshot shows a duplicated `G4bool SiPMSD::ProcessHits(...)` signature line;
  verify during rebase because it may be a compile blocker.

Resolution rule: keep pulse/electronics infrastructure, but make the default
S13360-6050PE + FastIC+ and compute timestamp smearing as:

```cpp
sigmaTotal = sqrt(fSiPMSPTR*fSiPMSPTR + fElectronicsJitter*fElectronicsJitter)
```

Then pass `t_global + Gauss(0, sigmaTotal)` into pulse generation / CFD.

#### `src/DetectorConstruction.cc` / `include/DetectorConstruction.hh`

Classification: **B plus A**.

Useful/intended:

- `/det/scintillatorMaterial EJ230|EJ200` may be useful later.

Accidental/out-of-scope drift:

- Default is EJ230. This contradicts the current EJ-200 baseline.
- Active Mylar wrap geometry is stale relative to `main`.
- SiPMs remain outside wrap volumes, not the final `BarLV` daughter topology.

Resolution rule: default must be EJ-200 after rebase. If selector is retained,
it must be off by default and not alter canonical `run.mac`.

#### `src/RunAction.cc`

Classification: **C with compatibility risk**.

The branch adds ntuple columns. This is acceptable if analysis reads branches by
name. It is risky if older ROOT macros use positional column indices. During
rebase, preserve canonical branch names and update analysis readers to use
branch names.

#### `macros/run.mac` / `macros/scan.mac`

Classification: **A) accidental scope drift**.

The branch sets:

```text
/det/scintillatorMaterial EJ230
/sipm/sensorModel Broadcom
```

These are not valid defaults for the SHiP 2020 EJ-200 + S13360-6050PE baseline.
They should be removed from default production macros and moved to explicit
future-study macros.

