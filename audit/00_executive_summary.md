# Executive Summary — Cross-Branch Physics Audit

Date: 2026-05-07  
Repository worktree audited: `/home/reriosto/SHiP/ej200_edge_scan`  
Current branch at end of audit: `main`  
Important: no commits, rebases, or merges were performed. Only audit files were
created under `audit/`.

## Branches Audited

| Branch | Purpose | Physics files touched / drift observed |
|---|---|---|
| `main` | Current baseline with explicit reflector panels | `Materials.cc`, `SiPMSD.cc`, `DetectorConstruction.cc`, `RunAction.cc`, `SteppingAction.cc` |
| `feature/edge-scan-and-readout-grouping` | Edge scan, wrapping modes, displays, Sum-of-N grouping | Stale `Materials` and `DetectorConstruction` versus `main`; still uses active Mylar wrap |
| `feature/sipm-electronics-response` | SiPM pulse/electronics response, sensors, CFD | Changes `Materials`, `SiPMSD`, `DetectorConstruction`, `RunAction`; adds EJ-230/Broadcom defaults inconsistent with baseline |

## Compact Coherence Table

| Parameter | Canonical | `main` | edge-scan | electronics |
|---|---|---|---|---|
| EJ-200 rise time | 0.9 ns | **missing** | **missing** | **missing for EJ-200** |
| EJ-200 decay | 2.1 ns | OK | OK | OK |
| EJ-200 yield | 10000 ph/MeV | OK | OK | OK |
| EJ-200 attenuation | 3.8 m | OK | OK | OK |
| EJ-200 RINDEX | 1.58 | OK | OK | OK |
| Reflector | explicit panels R=0.98 baseline | OK | **stale Mylar wrap** | **stale Mylar wrap** |
| SiPM PDE label | S13360-6050PE | **says S13360-6025** | **says S13360-6025** | **says S13360-6025 legacy** |
| PDE peak | 0.405 @ 460 nm | OK | OK | OK for S13360 curve |
| Jitter | SPTR 106 ps + FastIC+ 10 ps | **20 ps monolithic** | **20 ps monolithic** | **split exists but wrong defaults and electronics not applied** |
| Bar geometry | 1400 x 60 x 10 mm | OK | OK | OK |
| Top SiPM default | 20 at 70 mm pitch | OK | OK | OK |
| Default scintillator | EJ-200 | OK | OK | **EJ-230** |
| Default sensor | S13360-6050PE | OK numerically, mislabeled | OK numerically, mislabeled | **Broadcom** |

## Prioritized Discrepancies

### P0 — Blocks Publication

1. **EJ-200 rise time is absent in all branches.**  
   Canonical value: `SCINTILLATIONRISETIME1 = 0.9 ns`. Source:
   `EJ200.pdf`, `EJ200_EJ204_EJ208_EJ212.pdf`, and user/director meeting
   17-Apr-2026. Impact: Geant4 instantaneous onset artificially narrows timing.

2. **Prototype SiPM PDE table is mislabeled.**  
   Code comments say S13360-6025; SHiP prototype is S13360-6050PE
   (Betancourt 2020 / W4). Numeric table can remain unchanged.

### P1 — Affects Interpretation

3. **Jitter model is not canonical.**  
   Canonical default should be SPTR = 106 ps sigma for S13360-6050PE plus
   FastIC+ electronics jitter = 10 ps sigma, summed in quadrature.

4. **Feature branches are geometrically stale.**  
   Both features still use active Mylar wrap volumes while `main` now uses
   `BarPV` directly in `WorldLV` plus explicit reflector panels.

5. **`feature/sipm-electronics-response` changes default material/sensor.**  
   It defaults to EJ-230 + Broadcom, but the canonical baseline must remain
   EJ-200 + S13360-6050PE. Those options should become explicit study modes.

### P2 — Cosmetic / Maintenance

6. Boundary-census accessor names still mention `Mylar` historically, although
   user-facing labels in `main` were updated to `Bar -> reflector/World/SiPM`.

7. `CreateBarSkinReflector()` name is historical; it now returns a border
   surface used on explicit reflector panels, not a skin surface.

## Recommended Action Plan

1. Create `fix/physics-baseline` from `main`.
2. Commit 1: add EJ-200 rise time `0.9 ns` in `Materials::CreateEJ200()`.
3. Commit 2: correct S13360 PDE comments to S13360-6050PE without changing
   numeric PDE values.
4. Commit 3: split SiPM timing into SPTR + electronics jitter with defaults
   `106 ps` and `10 ps`, summed in quadrature; keep `/sipm/jitterSigma` as a
   deprecated compatibility alias.
5. Rebase `feature/edge-scan-and-readout-grouping` onto `fix/physics-baseline`.
   Baseline physics wins; keep edge-scan/display/grouping behavior.
6. Rebase `feature/sipm-electronics-response` onto `fix/physics-baseline`.
   Keep pulse/CFD infrastructure, but restore EJ-200 + S13360-6050PE defaults.
7. Add a future `tests/physics_baseline_check.cc` CTest to lock canonical
   values.
8. Move EJ-204/EJ-230 material selection and Broadcom/S14160 sensor selection
   to explicit future study branches/macros, not default baseline.

## Audit Outputs

- `audit/01_branch_inventory.txt`
- `audit/02_physics_per_branch.md`
- `audit/03_cross_branch_coherence.md`
- `audit/04_drift_classification.md`
- `audit/05_remediation_plan.md`

