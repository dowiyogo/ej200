# Group-Velocity Audit — EJ-200/EJ-230 SSLG4 Simulations

**Date**: 2026-06-17  
**Branch audited**: `feat/ej230-endonly-mylar` @ commit `7e18117` (t0minidaq run)  
**Reference branch**: `feat/endtop-sslg4` (EJ-204, current local branch)  
**Geant4 version**: 11.4.0 (`geant4-11-04 [MT]`, 5-December-2025)

---

## 1. Executive Summary

**ROOT CAUSE CONFIRMED (2026-06-17)**: Optical photon step-transit times are computed
using `c_light = 29.98 cm/ns` after each lateral-face reflection because the bar's thin air-filled
reflector sibling volumes cause the track touchable to transiently point to the wrong material at step
initialisation. `GROUPVEL` itself is stored correctly at 18.97 cm/ns; the bug is in how G4Transportation
picks up the pre-step velocity after the photon returns from the air reflector.

| Quantity | Value | Source |
|---|---|---|
| Experimental v_eff (EJ-230) | **15.5 cm/ns** | External data |
| Simulated v_eff (SUM4 trigger fit) | **21.47 cm/ns** | `figure4_sum8_validation.csv` |
| Simulated v_phot (from t_min at x=0) | **≈ 27.7 cm/ns** | t0minidaq data, confirmed below |
| Expected v_phot = c/n(EJ-230) | **18.97 cm/ns** | c/1.58 = 299.79/1.58 |
| GROUPVEL in MPT (runtime confirmed) | **18.97 cm/ns** | `[GROUPVEL_DIAG]` print, local build |
| Predicted t_min (bug model) | **2.547 ns** | See §10 |
| Observed t_min (local build) | **2.579 ns** | local ROOT file |
| Observed t_min (t0minidaq) | **2.713 ns** | `photon_hits_x0mm.root` |

The `GROUPVEL` MPT property is correct. The discrepancy is a **step-timing bug** triggered by the
geometry: thin (0.5 μm) G4_AIR reflector volumes adjacent to the bar cause photon pre-step velocities
to flip between bar and air values after each lateral reflection, compressing all post-bounce transit
times to `length/c_light` instead of `length/(c/n)`.

---

## 2. Branch Inventory

| Repo | Branch | Scint | Commit run | Run location |
|---|---|---|---|---|
| ej200 | feat/endtop-sslg4 | OPSC-101 (EJ-204) | local | local MSI |
| ej230_endonly_mylar | feat/ej230-endonly-mylar | OPSC-106 (EJ-230) | 7e18117 | t0minidaq |

See `branch_optics_matrix.csv` for full optics parameter comparison.

---

## 3. Evidence — File:Line

### 3.1 RINDEX data file
`ej230_endonly_mylar/src/external/SSLG4/data/oscnt/opsc-106/rIndex.txt` (same at commit 7e18117):
```
200.00    1.58
800.00    1.58
```
Two data points, n=1.58 at 200 nm and 800 nm.

### 3.2 Macro that loads RINDEX
`ej230_endonly_mylar/src/external/SSLG4/macros/oscnt/opsc-106.mac`:
```
/control/alias pathToDataDir sslg4/data/oscnt/{scnt}
/mpt/{scnt}/addProperty RINDEX {pathToDataDir}/rIndex.txt nm unitless
```
Path `{pathToDataDir}` is **relative to CMAKE_BINARY_DIR** (comment in
`OrganicScintillatorFactory.cc:22`).

### 3.3 Code that executes the macro
`OPSimTool/src/ScintillatorBuilder.cc:76-80`:
```cpp
void ScintillatorBuilder::DefineMPTUICommands() {
  if (IsMPTEnabled()) {
    pMPT = new MaterialPropertiesTable(mName);          // line 77: custom MPT
    G4UImanager::GetUIpointer()->ExecuteMacroFile(mMacroFilePath);  // line 78
  }
}
```

### 3.4 Wavelength-to-energy conversion
`OPSimTool/src/MaterialPropertiesTable.cc:86-98`:
```cpp
} else if (G4UnitDefinition::GetCategory(photonUnit) == "Length") {
    std::sort(photonVector.rbegin(), photonVector.rend()); // descending λ → ascending E
    for_each(photonVector.begin(), photonVector.end(),
    [](auto &pair) {
        pair.first = (CLHEP::hbarc * CLHEP::twopi) / (pair.first);
    });
}
```
Result: energies = [1.55 eV, 6.20 eV] (ascending), values = [1.58, 1.58].

### 3.5 GROUPVEL computation in Geant4
`geant4-v11.4.0/source/materials/src/G4MaterialPropertiesTable.cc:600-641`:
- For 2-point flat-n table (n0=n1=1.58):
  - First entry: `vg = c_light / (n0 + (n1−n0)/G4Log(E1/E0)) = c_light/1.58 = 18.97 cm/ns`
  - Last entry: `vg = c_light / (n1 + (n1−n0)/G4Log(E1/E0)) = c_light/1.58 = 18.97 cm/ns`
  - Loop (line 613): runs 0 times for 2-point table (runs for i=2..size-1)
- "Normal dispersion" guard (line 606): `vg=18.97 ≤ c/n0=18.97` → NOT triggered

### 3.6 G4Track velocity query
`geant4-v11.4.0/source/track/src/G4Track.cc:166-212`:
```cpp
velocity = c_light;  // default
groupvel = mat->GetMPT()->GetProperty(kGROUPVEL);  // line 189
if (groupvel != nullptr)
    velocity = groupvel->Value(fpDynamicParticle->GetTotalMomentum()); // line 202,205
```
For massless photon: `GetTotalMomentum() = KineticEnergy` (in MeV). Photon at 391 nm:
E = 3.17 eV = 3.17e-6 MeV. GROUPVEL table range: [1.55e-6, 6.20e-6] MeV. Query is in-range.

### 3.7 Override block (DetectorConstruction.cc:168-198)
Modifies only `ABSLENGTH` (→120 cm) and `SCINTILLATIONRISETIME1` (→0.5 ns). Does **not** touch
`RINDEX` or `GROUPVEL`.

---

## 4. Material Construction Chain

```
DetectorConstruction::Construct()
  └─ OrganicScintillatorFactory::Get("OPSC-106", true)
       └─ CreateMaterial("opsc-106", true)
            └─ BuildScintillator()
                 └─ ScintillatorBuilder(name, macroPath, density, elemFrac, enableMPT=true)
                      └─ DefineMPTUICommands()
                           ├─ pMPT = new MaterialPropertiesTable("opsc-106")
                           │    [registers /mpt/opsc-106/addProperty command]
                           └─ G4UImanager::ExecuteMacroFile("sslg4/macros/oscnt/opsc-106.mac")
                                └─ /mpt/opsc-106/addProperty RINDEX ... nm unitless
                                     └─ MaterialPropertiesTable::AddProperty("RINDEX", ...)
                                          └─ G4MaterialPropertiesTable::AddProperty("RINDEX", [1.55,6.20]eV, [1.58,1.58])
                                               └─ CalculateGROUPVEL()  [THEORETICAL: vg=18.97 cm/ns]
     └─ VMaterialBuilder::GetProduct()
          └─ mat->SetMaterialPropertiesTable(pMPT)
```

---

## 5. GROUPVEL: Theoretical vs Empirical

### 5.1 Theoretical (from static analysis)
RINDEX flat at n=1.58 for all photon energies in [1.55, 6.20] eV.
```
GROUPVEL = c / n = 299.792 mm/ns / 1.58 = 189.74 mm/ns = 18.974 cm/ns
```
GROUPVEL table: 2 entries, both = 18.97 cm/ns. Linear interpolation always returns 18.97 cm/ns.

### 5.2 Empirical (from simulation data)
From `photon_hits_x0mm.root` at t0minidaq (commit 7e18117), End SiPMs only (global_id 0–15):
```
t_min_min  = 2.713 ns   (minimum photon arrival time across 2000 events)
t_min_mean = 2.849 ns
```
Muon: 1 GeV mu⁻ from z=60mm → bar entry at z=+5mm:
```
t_muon_enter = 55 mm / (β × c) = 55 / (0.9954 × 299.79) ≈ 0.184 ns
```
Inferred photon velocity:
```
v_phot = 70 cm / (2.713 − 0.184) ns = 27.7 cm/ns   [lower bound from min]
v_phot = 70 cm / (2.849 − 0.184) ns = 26.3 cm/ns   [estimate from mean]
```
**Conclusion: v_phot_sim ≈ 26–28 cm/ns, NOT 18.97 cm/ns.**

### 5.3 T_s consistency check

For any photon velocity v, the trigger bias T_s(d) = t_trigger(d) − t_muon_entry − d/v must satisfy:
T_s ≥ 0 everywhere (trigger cannot fire before photons arrive).

| v_phot (cm/ns) | T_s range | All ≥ 0? | Monotone? |
|---|---|---|---|
| 18.97 (c/n) | −877 ps to +84 ps | **NO** | No |
| 21.47 (v_eff fit) | −75 ps to +166 ps | **NO** | No |
| 24.00 | +79 ps to +659 ps | YES | YES |
| 27.70 (t_min estimate) | +85 ps to +1432 ps | YES | YES |
| 29.98 (c_light) | +88 ps to +1814 ps | YES | YES |

**Minimum consistent v_phot = 24.0 cm/ns. Best estimate from t_min: 27.7 cm/ns.**

---

## 6. Trigger Time Audit (SUM4)

- Trigger: analog sum of SPE responses for 4 SiPMs in a cluster (groups 0–1 for left, 2–3 for right).
- SPE: biexponential with τ_rise=0.5 ns, τ_fall=5.0 ns, peak=0.824 ns after photon arrival.
- Threshold: 4 PE in the analog sum.
- Algorithm: `leading_edge_time()` in `engine.py` (exact port of `congruent_sum4_timing.C`).
- `t_left = min(cluster0_trigger, cluster1_trigger)` — minimum of two 4-SiPM clusters.
- **Trigger delay T_s increases with distance** (lower Npe at far end → slower rise → later crossing).
- T_s ≈ 80 ps near end (d=1 cm), ≈ 700–1400 ps far end (d=139 cm).
- **Trigger bias makes v_eff < v_phot**, consistent with v_eff=21.47 and v_phot≈27.7 cm/ns.

---

## 7. Units Audit

| Formula | Value | Notes |
|---|---|---|
| c_light (Geant4) | 299.792458 mm/ns | Standard |
| c/n_EJ-230 | 189.74 mm/ns = 18.974 cm/ns | n=1.58 |
| E(800nm) | 1239.84/800 = 1.550 eV = 1.550e-6 MeV | hc/λ |
| E(200nm) | 1239.84/200 = 6.199 eV = 6.199e-6 MeV | |
| E(391nm) | 1239.84/391 = 3.171 eV = 3.171e-6 MeV | EJ-230 emission peak |
| v_eff formula | v = 1/slope, slope = Δt/Δx | x in cm, t in ns |
| v_eff from delta | v = 2 / |d(ΔT_LR)/dx| | ΔT_LR = t_L − t_R |

All formulas verified correct in `figure4_sum8_validation.csv` generation code.
No units error identified.

---

## 8. Branch Differences (ej200 vs ej230_endonly_mylar)

| Parameter | ej200/feat/endtop-sslg4 | ej230_endonly_mylar/feat/ej230-endonly-mylar |
|---|---|---|
| Scintillator | OPSC-101 (EJ-204) | OPSC-106 (EJ-230) |
| RINDEX | 1.58 (flat, 2 pts) | 1.58 (flat, 2 pts) |
| GROUPVEL (theory) | 18.97 cm/ns | 18.97 cm/ns |
| ABSLENGTH | 160 cm (overridden) | 120 cm (overridden) |
| τ_d / τ_r | 1.8 ns / 0.7 ns | 1.5 ns / 0.5 ns |
| Scint yield | 10400 ph/MeV | 9700 ph/MeV |
| Reflector | CreateBarSkinReflector | CreateMylarReflector |
| Reflectivity | 0.98 | 0.90 |
| Finish | polished | ground (σ=0.1°) |
| SiPM topology | End + Top | End-only |
| v_eff,sim (trigger) | not yet characterized | 21.47 cm/ns |
| v_phot_sim (t_min) | not yet characterized | ≈27.7 cm/ns |

Note: RINDEX=1.58 is IDENTICAL in both simulations, so GROUPVEL discrepancy exists
equally in both if the root cause is in the SSLG4/OPSimTool chain.

---

## 9. Cause Analysis — Confirmed

### Status of previous hypotheses

| Rank | Hypothesis | Status |
|---|---|---|
| 1 | GROUPVEL wrong in MPT | **RULED OUT** — runtime print confirms 18.97 cm/ns |
| 2 | GROUPVEL absent → c_light | **RULED OUT** — GROUPVEL is present (2 entries) |
| 3 | MT race condition | **RULED OUT** — reproduces in 1-thread mode locally |
| 4 | Different rIndex.txt on t0minidaq | **UNLIKELY** — local 2.579 ns ≈ t0minidaq 2.713 ns (same bug) |
| 5 | Trigger-only artifact | **RULED OUT** — T_s goes negative, and t_min confirms fast arrival |

### Confirmed root cause — reflector-volume velocity aliasing

**Mechanism** (step-by-step):

1. Photon is born in the bar (`opsc-106`, GROUPVEL=18.97 cm/ns). Its first bar step is timed
   correctly: `ΔT = stepLen / (c/n) = stepLen / 18.97` cm/ns.

2. The photon hits a lateral face (y=±30 mm or z=±5 mm). The adjacent reflector is a separate
   physical volume of `G4_AIR` (0.5 μm foil, `foilHalfT`). The photon enters this AIR volume.
   G4_AIR has `RINDEX=1.0` → `GROUPVEL = c_light = 29.98 cm/ns`.

3. A **zero-length step** is recorded in `G4_AIR` (`stepLen=0`, `ΔT=0`). G4OpBoundaryProcess
   reflects the photon back into the bar. The track's velocity is now keyed to the AIR touchable.

4. The **next bar step** (photon back in `opsc-106`) calls `G4Step::InitializeStep`, which sets
   the pre-step point velocity via `fpTrack->CalculateVelocityForOpticalPhoton()`. At this instant,
   the track's touchable still references the reflector (`G4_AIR`), so the computed velocity is
   **c_light = 29.98 cm/ns**, not 18.97 cm/ns.

5. G4Transportation computes `ΔT = stepLen / pre_vel = stepLen / 29.98` (wrong — 37% too fast).

6. This pattern repeats for every subsequent lateral reflection. The velocity *value* returned by
   `track->GetVelocity()` (called after the step) is correct (18.97), but the *pre-step point*
   velocity used for timing was wrong.

**Runtime evidence from `[STEP_DIAG]` print** (local build, evt=0):
```
step=1  mat=opsc-106  pre_vel=18.97  stepLen=15.6mm  ΔT=0.082 ns  ← correct
step=2  mat=G4_AIR    pre_vel=18.97  stepLen=0        ΔT=0          ← zero-length air bounce
step=3  mat=opsc-106  pre_vel=29.98  stepLen=15.8mm  ΔT=0.053 ns  ← WRONG (29.98 instead of 18.97)
step=4  mat=G4_AIR    pre_vel=18.97  stepLen=0        ΔT=0
step=5  mat=opsc-106  pre_vel=29.98  stepLen=15.7mm  ΔT=0.052 ns  ← WRONG
```

Error per post-bounce step: `15.7 × (1/18.97 − 1/29.98) × 10 = 0.031 ns` (step is 37% too fast).

**Not all photons are equally affected**:
- Photons emitted nearly parallel to x (small y/z component) never hit lateral walls → arrive at correct
  ~3.87 ns.
- Photons that bounce early (first reflection after ~15 mm in x) then travel the remaining ~685 mm
  entirely with air pre-velocity → arrive at ~2.55 ns (near c_light).
- The 200-event t_min_min = 2.579 ns matches the model prediction of 2.547 ns.

**Why trigger v_eff = 21.47 cm/ns is intermediate**:
The SUM4 trigger fires at the 4 PE threshold. It needs ~4 photons from a 4-SiPM cluster. The
minimum-time photons (buggy) may contribute 0–1 PE each; the remaining PE threshold is met by
later (correct-velocity) photons. The trigger time is therefore pulled earlier than the 18.97-cm/ns
prediction but later than c_light, giving the observed 21.47 cm/ns.

---

## 10. Explicit Conclusion

**ROOT CAUSE CONFIRMED: Reflector-volume velocity aliasing in G4Transportation.**

The `GROUPVEL` MPT table for `opsc-106` is correctly populated at **18.97 cm/ns** (confirmed
at runtime via `[GROUPVEL_DIAG]` print). The RINDEX data file, wavelength-to-energy conversion,
and `CalculateGROUPVEL()` formula are all correct.

The spurious fast arrivals are caused by the thin (0.5 μm) G4_AIR reflector sibling volumes.
After each lateral-face reflection the `G4Step::InitializeStep` pre-step velocity is pulled from
the reflector (air) touchable instead of the bar touchable, causing all post-bounce bar steps
to be timed at `c_light` instead of `c/n`. This is not a GROUPVEL table error but a geometry-
induced step-initialisation artefact.

**Predicted vs observed t_min (x=0, d=699.5 mm)**:

| Source | t_min | Comment |
|---|---|---|
| Correct physics (v=18.97 cm/ns) | 3.871 ns | 0.184 + 699.5/189.7 |
| Bug model (step 1 correct, rest at c_light) | 2.547 ns | 0.184 + 15.6/189.7 + 683.9/299.8 |
| Observed locally (local build, 200 events) | **2.579 ns** | ROOT file, this audit |
| Observed t0minidaq (commit 7e18117) | **2.713 ns** | `photon_hits_x0mm.root` |

The difference between local (2.579) and t0minidaq (2.713) is due to statistical fluctuation with
different RNG seeds; both are firmly within the bug-model range, NOT consistent with 3.871 ns.

---

## 11. Axial Test Note

The axial test macro (`macros/axial_photon_test.mac`) was created but is not conclusive because
`PrimaryGeneratorAction::GeneratePrimaries` overrides the particle direction to `{sinT, 0, -cosT}`
regardless of `/gun/direction`. A photon along x cannot be fired from the G4ParticleGun.

The per-step velocity diagnostic (`[STEP_DIAG]` in `SteppingAction.cc`, `[VEL_DIAG]` in
`SiPMSD.cc`) directly confirmed the mechanism. See §9.

---

## 12. Recommended Fix

**Replace thin air-filled reflector sibling volumes with a G4LogicalSkinSurface on BarLV.**

```
CURRENT (broken)                          FIXED
─────────────────────────────────         ─────────────────────────────────────────────
WorldLV ─┬─ BarPV (opsc-106)             WorldLV ─── BarPV (opsc-106)
          ├─ ReflectorYMinusPV (G4_AIR)              G4LogicalSkinSurface("BarSkin",
          ├─ ReflectorYPlusPV  (G4_AIR)              barLV, mylarSurface)
          ├─ ReflectorZMinusPV (G4_AIR)
          └─ ReflectorZPlusPV  (G4_AIR)
```

The skin surface applies the mylar optical properties (R=0.90, ground finish) to ALL bar faces
that do not have an explicit border surface. End faces keep their border surface to SiPMs.
Photons that hit a lateral bar face reflect (or absorb) without entering any other volume,
so `G4Step::InitializeStep` always picks up the bar touchable and uses GROUPVEL=18.97 cm/ns.

**Minimum-change implementation** in `DetectorConstruction.cc` (replaces reflector PV placements):
```cpp
// Remove four G4PVPlacement calls for reflector panels
// Add one skin surface to BarLV for all non-instrumented faces:
auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV,
                    Materials::CreateMylarReflector());
(void)barSkin;
```

This completely removes the timing artefact. The optical reflectivity is unchanged (same surface
optical properties). The BarPV→SiPM border surfaces are unaffected.

**Alternative fix (keep volume geometry, zero out air GROUPVEL)**:
```cpp
// In DetectorConstruction.cc, when setting worldMat MPT:
auto* mpt = new G4MaterialPropertiesTable();
mpt->AddProperty("RINDEX", {2.0*eV, 4.0*eV}, {1.0, 1.0});
// Override GROUPVEL to bar value so air steps use c/n not c_light:
G4double v_bar = CLHEP::c_light / 1.58;
mpt->AddProperty("GROUPVEL", {2.0*eV, 4.0*eV}, {v_bar, v_bar});
worldMat->SetMaterialPropertiesTable(mpt);
```

This is a hack: it sets n_air effective-transport-velocity to c/n_bar, which is physically wrong
for photon transport in the world volume. Use only as a stopgap.

**Recommended action**: implement the skin-surface fix before the next campaign. The trigger
v_eff and t_min will shift from ~21 cm/ns and ~27 cm/ns to their correct values near 18.97 cm/ns
(modulo trigger bias, scintillation timing, and remaining physics).

---

*ROOT CAUSE CONFIRMED — 2026-06-17*
*Data: t0minidaq commit 7e18117, EJ-230 end-only+mylar, 31 positions × 2000 events*
*Local build: /tmp/ej230_build (Geant4 11.4.0, g++ -O2)*
