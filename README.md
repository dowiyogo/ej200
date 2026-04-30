# ej200_v2 — EJ-200/EJ-230 Scintillating Bar + SiPM Simulation (v2)

Improved version of `scintillator_geant4` (Repo1), incorporating the best
practices from `simple_g4sim_scint_sipm` (Repo2).

---

## Geometry

| Volume | Dimensions | Material |
|---|---|---|
| EJ-230 bar (default) | 1400 × 60 × 10 mm | Polyvinyltoluene + EJ-230 optical props |
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
| PDE | Not applied | Manual Bernoulli trial in `SiPMSD::ProcessHits()` via `GetPDE()` |
| Wrapping | None / default | `dielectric_dielectric polished` (TIR-based, bar–air interface) |
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

The muon (1 GeV μ⁻) starts 60 mm above the bar centre on the wide face (+Z),
travels in the **−Z direction**, and traverses the full 10 mm thickness of the bar.
Scans step along the **X axis** (longitudinal direction).

```bash
# Standard production (1000 events, mu- at bar centre x=0)
./ej200_bar_sim -m macros/run.mac

# Full longitudinal scan (21 x positions × 200 events each, −650 → +650 mm)
./ej200_bar_sim -m macros/scan.mac

# Edge scan: 11 positions × 500 events from x=+600 to +700 mm
./ej200_bar_sim -m macros/scan_edge.mac
python analysis/edge_resolution.py

# Compare edge wrapping modes in the last 5 cm of the bar
./ej200_bar_sim -m macros/scan_edge_air.mac
./ej200_bar_sim -m macros/scan_edge_black.mac
python analysis/compare_edge_wraps.py mylar air black

# Interactive visualisation
./ej200_bar_sim
```

For longitudinal scans, prefer `/muon/gunX` over `/gun/position` in macros so
the event-level `gun_x_mm` written to the ntuple tracks the intended scan
coordinate unambiguously.

The SiPM model is selected in the macro after `/run/initialize`, for example:

```tcl
/sipm/sensorModel Broadcom
#/sipm/sensorModel S14160-6010PS
#/sipm/sensorModel S14160-6015PS
```

The scintillator material is selected before `/run/initialize`:

```tcl
/det/scintillatorMaterial EJ230
#/det/scintillatorMaterial EJ200
```

---

## Analysis

```bash
pip install uproot numpy matplotlib pandas scipy

# General plots (single-position or scan data, includes dCFD timing)
python analysis/analyze_dCFD.py photon_hits.root

# Temporal resolution vs longitudinal position (requires scan data)
python analysis/resolution_vs_x_dCFD.py photon_hits.root
# FPT-only legacy variant:
# python analysis/resolution_vs_x_FPT.py photon_hits.root

# Edge timing/light-yield study
python analysis/edge_resolution.py photon_hits_run*.root --label mylar

# Edge-wrap comparison from three directories containing edge_summary.csv
python analysis/compare_edge_wraps.py mylar air black
```

### `analyze_dCFD.py` outputs

- `photons_per_sipm.pdf` — total hits per SiPM
- `top_sipm_profile.pdf` — top SiPM hits vs x-position (validates full coverage)
- `arrival_time.pdf` — timing spectrum per face type
- `wavelength.pdf` — detected photon wavelength spectrum
- `end_asymmetry.pdf` — (L−R)/(L+R) asymmetry distribution

### `resolution_vs_x_dCFD.py` outputs (scan data required)

- `resolution_vs_x.pdf` — **key thesis plot**: σ_t [ps] vs x for end and top SiPMs
- `fpt_dist_end.pdf` — first-photon-time distributions at 3 representative positions (end SiPMs)
- `fpt_dist_top.pdf` — first-photon-time distributions at 3 representative positions (top SiPMs)
- `asymmetry_vs_x.pdf` — end-SiPM charge asymmetry ⟨(N_L−N_R)/(N_L+N_R)⟩ vs x (position reconstruction)
- `n_photons_vs_x.pdf` — mean detected photons per event vs x, by face type

---

## SiPM/MPPC characteristics

`SiPMSD` applies PDE manually with linear interpolation over the selected sensor
curve and stores `sensor_model_id`, `sensor_sptr_sigma_ns`, and
`electronics_sigma_ns` in the ntuple.

Supported sensor curves.  SPTR values are simulation defaults; the supplied
datasheets provide the PDE/electrical curves but do not quote a numeric SPTR
sigma for every model.

- `Broadcom` / `AFBR-S4N66P024M`: Broadcom AFBR-S4N66P024M datasheet curve, 250–900 nm, peak PDE 63% at 420 nm, default SPTR sigma 30 ps.
- `S14160-6010PS` / `S14160-3010PS` / `S14160-1310PS`: Hamamatsu S14160 10 um-pixel MPPC curve, 290–900 nm, peak PDE 18% at 460 nm, default SPTR sigma 150 ps.
- `S14160-6015PS` / `S14160-3015PS` / `S14160-1315PS`: Hamamatsu S14160 15 um-pixel MPPC curve, 290–900 nm, peak PDE 32% at 460 nm, default SPTR sigma 150 ps.
- `Hamamatsu` / `S13360-6025`: legacy Hamamatsu S13360-6025 curve, 300–940 nm, peak PDE about 40.5% at 460 nm.

---

## Extending

- **Change SiPM count/positions**: edit `kNTopSiPMs` and `TopSiPMCenterX()` in
  `DetectorConstruction.hh/.cc`.
- **Change particle/energy**: edit `macros/run.mac` or use interactive `/gun/` commands.
- **Add 2D SiPM array on top**: extend the top SiPM placement loop to use both X and Z.
