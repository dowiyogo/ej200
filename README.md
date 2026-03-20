# ej200_v2 — EJ-200 Scintillating Bar + SiPM Simulation (v2)

Improved version of `scintillator_geant4` (Repo1), incorporating the best
practices from `simple_g4sim_scint_sipm` (Repo2).

---

## Geometry

| Volume | Dimensions | Material |
|---|---|---|
| EJ-200 bar | 1400 × 60 × 10 mm | `G4_PLASTIC_SC_VINYLTOLUENE` + EJ-200 optical props |
| End SiPMs (×16) | 6 × 6 × 0.5 mm | Optical coupling (n=1.58) |
| Top SiPMs (×20) | 6 × 6 × 0.5 mm | Optical coupling (n=1.58) |

### SiPM layout

```
Global IDs   Face         Count   Description
  0 –  7     End left     8       8×1 array on −X face
  8 – 15     End right    8       8×1 array on +X face
 16 – 35     Top          20      uniform, −665 → +665 mm (step 70 mm)
```

---

## Key improvements vs Repo1

| Feature | Repo1 | **v2** |
|---|---|---|
| SiPM surface | `LogicalSkinSurface` (global) | `LogicalBorderSurface` per SiPM |
| Optical coupling | Air (n=1.0) → TIR at ~39° | Coupling material n=1.58 (no TIR) |
| PDE | Manual C++ interpolation in SD | `DETECTIONEFFICIENCY` in surface MPT |
| Wrapping | None / default | Tyvek-like, R=0.95, `dielectric_metal` |
| Top SiPMs | Only 2nd half of bar (bug) | Full bar length (−665 to +665 mm) |
| EventAction | Bug (AddEndLeftHit unconditional) | SiPMSD → EventAction → RunAction |
| Materials class | Inline in DetectorConstruction | Separate `Materials` namespace |

### Why the old top SiPMs didn't collect photons

In Repo1, the `LogicalSkinSurface` was applied to the **bar** logical volume,
meaning **all** surfaces of the bar (including those facing the top SiPMs) were
treated as the same surface type. The SiPM volumes were in air (n=1.0), causing
total internal reflection (TIR) at the bar–SiPM interface for photons hitting at
angles > arcsin(1.0/1.58) ≈ 39°. Additionally, the `SkinSurface` effectively
competed with the `BorderSurface` used for end SiPMs, with unpredictable
priority results.

In v2:
1. **Per-SiPM `BorderSurface`** (bar→sipmPhys) takes priority over the wrapping
   surface (bar→world) by G4 convention.
2. **Coupling material** (n=1.58) eliminates TIR: the critical angle is 90°,
   so all photons that reach the SiPM face can enter.
3. **20 top SiPMs** span the full bar length (not just one half).

---

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## Run

```bash
# Standard production (1000 events, mu- at centre)
./ej200_bar_sim -m macros/run.mac

# Full longitudinal scan (21 x positions × 200 events each)
./ej200_bar_sim -m macros/scan.mac

# Interactive visualisation
./ej200_bar_sim
```

---

## Analysis

```bash
pip install uproot numpy matplotlib pandas
python ../analysis/analyze.py photon_hits.root
```

Produces:
- `photons_per_sipm.pdf` — total hits per SiPM
- `top_sipm_profile.pdf` — **key plot**: top SiPM hits vs x-position (validates full coverage)
- `arrival_time.pdf` — timing spectrum per face type
- `wavelength.pdf` — detected photon wavelength spectrum
- `end_asymmetry.pdf` — (L−R)/(L+R) for position reconstruction

---

## PDE table

Hamamatsu S13360-6025 (33 points, 300–940 nm). Peak PDE ≈ 40.5% at 460 nm.
Stored in the `SiPMSurface` MPT (`DETECTIONEFFICIENCY`) and interpolated by
`SiPMSD::GetPDE()` using `G4PhysicsVector::Value()`.

---

## Extending

- **Change SiPM count/positions**: edit `kNTopSiPMs` and `TopSiPMCenterX()` in
  `DetectorConstruction.hh/.cc`.
- **Change particle/energy**: edit `macros/run.mac` or use interactive `/gun/` commands.
- **Add 2D SiPM array on top**: extend the top SiPM placement loop to use both X and Z.
