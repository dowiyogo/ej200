# Remediation Plan

No commits, rebases, or merges were performed during this audit. This is the
recommended execution plan for a later, approved iteration.

## 5.1 Create Canonical Physics Baseline Branch

Start from updated `main`:

```bash
git checkout main
git pull --ff-only origin main
git checkout -b fix/physics-baseline
```

### Commit 1: `fix(physics): add EJ-200 rise time from datasheet`

File: `src/Materials.cc`, in `Materials::CreateEJ200()`.

Add:

```cpp
// Rise time from EJ-200 datasheet (Eljen Technology, EJ200.pdf in project)
// Source: EJ-200 Plastic Scintillator datasheet, "Rise Time, ns ... 0.9"
mpt->AddConstProperty("SCINTILLATIONRISETIME1", 0.9 * ns);
```

Justification: `EJ200.pdf` and `EJ200_EJ204_EJ208_EJ212.pdf` specify EJ-200
rise time = 0.9 ns. Without this property, Geant4 uses instantaneous emission
onset, narrowing timing distributions artificially. This also implements the
user/director requirement from the 17-Apr-2026 meeting: progressive rise time
is mandatory for publishable timing.

### Commit 2: `fix(sipm): correct PDE table comment to S13360-6050PE`

Files: `src/SiPMSD.cc`, `include/SiPMSD.hh`.

Change comments only. Do not change PDE values.

Recommended comment:

```cpp
// PDE table — Hamamatsu S13360-6050PE (pixel pitch 50 um, 6x6 mm2)
// Model used in the SHiP 2020 prototype (Betancourt et al., NIM A 979 164398).
// 33 points in 300-940 nm. Peak 40.5% at 460 nm.
// Numeric values are consistent with the S13360-6050PE prototype curve; the
// previous comment ("S13360-6025") was incorrect.
```

Justification: Betancourt 2020 and W4 identify the prototype sensor as
S13360-6050PE. The current numeric PDE curve is acceptable for that prototype
model, but the label is wrong and would be embarrassing in a thesis.

### Commit 3: `feat(sipm): split jitter into SPTR and electronics components`

Files: `include/SiPMSD.hh`, `src/SiPMSD.cc`, and any macros/comments that
mention `/sipm/jitterSigma`.

Canonical fields:

```cpp
// Jitter has two independent physical contributions:
//   1) SPTR (Single Photon Time Resolution) intrinsic to the SiPM
//   2) FastIC+ electronics jitter (SOP-3.3, SHiP protocol)
// Total smearing is sqrt(SPTR^2 + electronics^2).
//
// Defaults for SHiP 2020 prototype:
//   SPTR = 106 ps sigma (Hamamatsu S13360-6050PE, ~250 ps FWHM, literature)
//   Electronics = 10 ps sigma (FastIC+, Protocolo SHiP SOP-3.3)
G4double fSiPMSPTR          = 106.0e-3;  // ns
G4double fElectronicsJitter = 10.0e-3;   // ns

void     SetSiPMSPTR(G4double s)          { fSiPMSPTR = s; }
void     SetElectronicsJitter(G4double s) { fElectronicsJitter = s; }
G4double GetSiPMSPTR() const              { return fSiPMSPTR; }
G4double GetElectronicsJitter() const     { return fElectronicsJitter; }
```

In `ProcessHits()`:

```cpp
const G4double sigmaTotal = std::sqrt(fSiPMSPTR*fSiPMSPTR
                                      + fElectronicsJitter*fElectronicsJitter);
const G4double jitter  = G4RandGauss::shoot(0.0, sigmaTotal);
const G4double time_ns = (track->GetGlobalTime() + jitter) / ns;
```

UI commands:

```text
/sipm/sptrSigma <value> ns
/sipm/electronicsJitter <value> ns
```

Keep legacy compatibility:

```text
/sipm/jitterSigma <value> ns
```

Alias behavior: set SPTR to the provided value and electronics jitter to 0, or
emit a clear warning explaining that the command is deprecated. Choose one
behavior and document it in `README.md`.

Do not introduce Broadcom or S14160 defaults in this baseline branch. Those are
future selectable models.

## 5.2 Rebase Plan for Feature Branches

### `feature/edge-scan-and-readout-grouping`

```bash
git checkout feature/edge-scan-and-readout-grouping
git rebase fix/physics-baseline
```

Conflict policy:

- Physical baseline wins:
  - EJ-200 rise time = 0.9 ns
  - PDE label = S13360-6050PE
  - SPTR/electronics split and defaults
  - `main` reflector topology: `BarPV` direct in `WorldLV`, explicit panels
    with validated R=0.98 fallback
- Keep feature behavior:
  - edge-scan macros and analysis
  - top-readout displays
  - Sum-of-N / FastIC+ grouping analysis
- Rework `/det/edgeWrap`:
  - It must not silently replace the canonical default reflector.
  - Either make it a clearly named non-baseline study mode, or move it to
    dedicated macros/analysis for comparison only.

Expected conflict files:

- `src/DetectorConstruction.cc`
- `include/DetectorConstruction.hh`
- `src/Materials.cc`
- `include/Materials.hh`
- `README.md`

### `feature/sipm-electronics-response`

```bash
git checkout feature/sipm-electronics-response
git rebase fix/physics-baseline
```

Conflict policy:

- Physical baseline wins:
  - default scintillator = EJ-200
  - default sensor = S13360-6050PE
  - default SPTR = 106 ps, electronics jitter = 10 ps
  - timestamp smear uses quadrature total
  - no default Broadcom/EJ230 in `run.mac` or `scan.mac`
  - current `main` reflector topology
- Keep feature behavior:
  - pulse-shape / waveform / CFD analysis
  - additive ntuple branches, if analysis reads by branch name
  - optional sensor-model selector, but default must be S13360-6050PE

Expected conflict files:

- `src/SiPMSD.cc`
- `include/SiPMSD.hh`
- `src/RunAction.cc`
- `src/DetectorConstruction.cc`
- `include/DetectorConstruction.hh`
- `src/Materials.cc`
- `include/Materials.hh`
- `macros/run.mac`
- `macros/scan.mac`

Pipeline resolution for electronics:

```text
t_global
  -> + Gauss(0, sqrt(SPTR^2 + electronics^2))
  -> pulse generation
  -> CFD / dCFD / time-walk correction
  -> t_hit stored in ntuple
```

Do not apply electronics jitter twice. If the pulse/CFD simulation already
models electronics noise explicitly, document whether `fElectronicsJitter`
represents residual timing jitter or is disabled.

## 5.3 Future Coherence Policy

Add a new section to `README.md` or create `docs/PHYSICS_BASELINE.md`:

1. Any branch that changes a physical parameter must cite the source in the
   commit message and code comment.
2. Any physical parameter visible to users must be reflected in README tables.
3. Feature branches may add study modes, but default macros must remain on the
   canonical baseline unless the branch is explicitly a physics-baseline change.
4. Analysis must read ROOT branches by name, not by fixed ntuple column index.

Suggested CTest: `tests/physics_baseline_check.cc`.

Checks:

- `Materials::CreateEJ200()` MPT has:
  - `SCINTILLATIONRISETIME1 == 0.9 ns`
  - `SCINTILLATIONTIMECONSTANT1 == 2.1 ns`
  - `SCINTILLATIONYIELD == 10000 / MeV`
  - `RESOLUTIONSCALE == 1.0`
  - `RINDEX == 1.58`
  - `ABSLENGTH == 3.8 m`
- `SiPMSD` defaults:
  - SPTR = 106 ps
  - electronics jitter = 10 ps
  - total sigma = sqrt(106^2 + 10^2) ps
- Optional smoke check that `run.mac` does not select EJ230 or Broadcom by
  default.

## 5.4 EJ-230 and EJ-204 Future Feature

Do not fold multi-material selection into the immediate physics-baseline fix.
Create a later branch:

```bash
git checkout -b feature/scintillator-materials-selection
```

Canonical values for that branch:

```cpp
// EJ-204 — Eljen Technology, EJ200_EJ204_EJ208_EJ212.pdf (Aug 2023)
//   SCINTILLATIONRISETIME1     = 0.7 * ns
//   SCINTILLATIONTIMECONSTANT1 = 1.8 * ns
//   SCINTILLATIONYIELD         = 10400 / MeV
//   ABSLENGTH                  = 1.6 * m       // 160 cm (datasheet)
//   RINDEX                     = 1.58
//   Peak emission              = 408 nm
//
// Note: W4 and the Protocolo document state 380 cm for EJ-204, but the
// official Eljen datasheet says 160 cm. Use 160 cm in code and leave a note
// for thesis-director discussion.

// EJ-230 — Eljen Technology, EJ228_EJ230.pdf (Aug 2023)
//   SCINTILLATIONRISETIME1     = 0.5 * ns
//   SCINTILLATIONTIMECONSTANT1 = 1.5 * ns
//   SCINTILLATIONYIELD         = 9700 / MeV
//   ABSLENGTH                  = 1.2 * m       // 120 cm
//   RINDEX                     = 1.58
//   Peak emission              = 391 nm
```

For the future sensor selection branch:

- S14160-6010PS: PDE 18% at 460 nm, datasheet
  `s141601310ps_etc_kapd1070e.pdf`; do not use it as default.
- Broadcom AFBR-S4N66P024M: PDE 63% at 420 nm, datasheet
  `AFBRS4N66P024MDS105.pdf`; do not call it FBK in code.

