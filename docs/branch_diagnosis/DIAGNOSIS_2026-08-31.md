# Branch Diagnosis — ej200 repository — 2026-08-31

**Author:** rrios (Claude Sonnet 4.6 assist)  
**Date:** 2026-08-31  
**REF_BASE:** `origin/feat/bar-end-vikuiti` @ `b007750`  
**Total branches diagnosed:** 16 (1 main + 15 feature/wip/exp)  
**Clones covered:** HOST (`/home/rrios/ej200`), W1 (`/mnt/d/SHiP/ej200`), W2 (`/mnt/d/SHiP/ej230`)

---

## §1 — Scope and method

All data is extracted from code (`git show <branch>:src/...`), `git diff --stat`, and `git cherry`. Branch names are not used as evidence. `INDETERMINADO` appears when the code is ambiguous or not read.

Security constraints respected throughout: no push (CP2 already resolved), no commit in user branches, no reset/rebase/merge, no branch deletions, no file deletions.

---

## §2 — Safety record (CP0–CP2)

### CP0 — Bundles

Three bundles created 2026-08-31, verified with `git bundle verify` (all PASS):

| Bundle | Path | Refs | SHA-256 | Device |
|--------|------|------|---------|--------|
| HOST | `/tmp/rrios/ej200_host_2026-08-31.bundle` | 42 | `b8f571483ca3aa855ee0f06d8f67dd9972d72af564c8a9827ed3690904f8c116` | HOST fs |
| W1 | `/mnt/c/SHiP_backups/ej200_w1_2026-08-31.bundle` | 61 | `1aa4c0bd2f63f7f2e4d6453c651a8c38cb45a5e45b9644d8fe022b1863d4714d` | C: drive (dev 70) |
| W2 | `/mnt/c/SHiP_backups/ej200_w2_2026-08-31.bundle` | 39 | `f2cc78ed7e0c6ceac9080692054f660d3d6b36f40f7d1e56438e959033d0d33a` | C: drive (dev 70) |

W1/W2 `.git` dirs live on D: drive (device 71); bundles on C: drive (device 70) — physically isolated from the WSL2 containers that hold the repos.

### CP1 — Fetch and branch table

All three clones fetched with `git fetch --all --prune --tags`. Branch table rebuilt from remote refs after prune.

### CP2 — ej228 rescue

**Not needed.** Both `feat/ej228-cylinder` (66b674c) and `feat/ej228-tir-only` (0006919) were already present on `origin` before this session. Confirmed via `git ls-remote --heads origin 'feat/ej228-*'`. No push was made.

### Stashes

Zero stashes in all three clones at start of CP3 (verified after user committed WIP branches).

---

## §3 — Physical configuration matrix

Each cell is extracted from code; see §8 for evidence references.

**Surface model codes:**
- **GEN-1**: `G4LogicalSkinSurface` on `barLV` with `dielectric_metal` + REFLECTIVITY=0.98 → TIR eliminated (×1500 photon deficit vs CORRECT)
- **NEAR-GEN-1**: `G4LogicalBorderSurface` on bar→panel boundary, `dielectric_metal`, no air gap → TIR still eliminated
- **TIR-ONLY**: `G4LogicalSkinSurface` on `barLV` with `dielectric_dielectric` polished → TIR works, no secondary reflector (photons that fail TIR lost to world)
- **2T-MYLAR**: two-tier; bar→air `dielectric_dielectric` + air→Mylar panel `dielectric_metal` R=0.95; air gap 0.10 mm, Mylar 0.05 mm
- **2T-VIK**: two-tier; bar→air `dielectric_dielectric` + air→Vikuiti panel `dielectric_metal` R=0.98; air gap 0.10 mm, Vikuiti 0.05 mm ← CURRENT BASELINE
- **3T-MYLAR**: three-tier; bar→Mylar wrap `dielectric_dielectric` (Fresnel) + Mylar→air TIR at 37.3°; no air gap; Mylar thickness 0.025 mm

| Branch | Geometry | Scint (default) | Readout modes | Surface | Air gap | Reflector | Gun default |
|--------|----------|-----------------|---------------|---------|---------|-----------|-------------|
| **feat/bar-end-vikuiti** (REF) | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | **2T-VIK** | 0.10 mm | Vikuiti 3M ESR R=0.98 | Through 10 mm (Z) |
| main | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | **GEN-1** | — | Vikuiti skin on barLV | Through 10 mm (Z) |
| feat/bar-vikuiti | Bar 1400×60×10 mm | OPSC-106 (EJ-230) | End / Top / EndTop | **GEN-1** | — | Vikuiti skin on barLV | Through 10 mm (Z) |
| feat/ej230-bar-tir-only | Bar 1400×60×10 mm | OPSC-106 (EJ-230) | End / Top / EndTop | TIR-ONLY | — | None | Through 10 mm (Z) |
| feat/ej204-bar-tir-only | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | TIR-ONLY | — | None | Through 10 mm (Z) |
| feat/ej228-cylinder | Cylinder ∅25 × 25 mm | EJ-228 (hardcoded) | End (8 cap SiPMs) | **2T-VIK** | 0.10 mm annular | Vikuiti wrap R=0.98 | Through 25 mm (Z caps) |
| feat/ej228-tir-only | Cylinder ∅25 × 25 mm | EJ-228 (hardcoded) | End (8 cap SiPMs) | TIR-ONLY | 0.10 mm annular | None | Through 25 mm (Z caps) |
| feat/ej230-endonly-mylar | Bar 1400×60×10 mm | OPSC-106 (EJ-230) | End only | NEAR-GEN-1 | — | Mylar panels (direct contact) | Through 10 mm (Z) |
| feat/ej230-sslg4 | Bar 1400×60×10 mm | OPSC-106 (EJ-230) | End / Top / EndTop | **GEN-1** | — | Vikuiti skin on barLV | Through 10 mm (Z) |
| feat/endonly-mylar | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End (+ "sipm" top) | **GEN-1** | — | Mylar/Vikuiti skin on barLV | Through 10 mm (Z) |
| feat/endtop-sslg4 | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | **2T-MYLAR** | 0.10 mm | Mylar R=0.95 | Through 10 mm (Z) |
| exp/pair-scan-2026-06-11 | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | **GEN-1** | — | Vikuiti skin on barLV | Through 10 mm (Z) |
| feat/ej204-event-display-tracks | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | **GEN-1** | — | Vikuiti skin on barLV | Through 10 mm (Z) |
| feature/sipm-electronics-response | Bar 1400×60×10 mm | EJ-200/204 (cfg) | INDETERMINADO (pitch-cfg) | **3T-MYLAR** | — | Mylar wrap 25 µm | Through 10 mm (Z) |
| wip/host-uncommitted-2026-08-31 | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop / **EndSparseTop** | **2T-VIK** | 0.10 mm | Vikuiti R=0.98 | Through 10 mm (Z) |
| wip/host-stash-endtop-junio | Bar 1400×60×10 mm | OPSC-101 (EJ-204) | End / Top / EndTop | NEAR-GEN-1 | — | Vikuiti panels (direct contact) | Through 10 mm (Z) |

---

## §4 — Per-branch detail

### REF_BASE — feat/bar-end-vikuiti @ b007750

**Unique commits vs REF_BASE:** 0 (this IS the baseline)  
**Code changes vs main:** +17 commits; introduces 2T-VIK geometry replacing GEN-1 skin surface  
**Key scripts:** `scripts/post_end_vikuiti.sh`, `scripts/run_end_vikuiti_scan.sh`, `scripts/run_end_tir_scan.sh`, `scripts/run_end_vik_sparse_top.sh`  
**Unique artifacts:** `analysis/presentation_v6/` (4-config Beamer comparison)  
**Defects:** Defects #1, #2, #3 (see §8)

---

### main @ 84e902c

**⚠ WARNING:** main has the **GEN-1 bug** (`G4LogicalSkinSurface` on `barLV`, `dielectric_metal` R=0.98).  
**Unique commits vs REF_BASE:** 0 commits unique to main (REF_BASE is 17 commits ahead of main)  
**Worktree-pinned in:** W1 worktree `/mnt/d/SHiP/ej200_edge_scan`  
**Destination:** REQUIERE DECISIÓN — pinned worktree; physics results from main are affected by GEN-1 bug

---

### feat/bar-vikuiti @ 219fbe3

**Unique commits vs REF_BASE:** 34  
**Surface:** GEN-1 bug  
**Default scint:** OPSC-106 (EJ-230) — differs from REF_BASE default  
**Key scripts:** `scripts/run_resolution_scan.sh`, `scripts/run_pair_scan.sh`, `scripts/analyze_geometry.py`  
**Unique artifacts:** 4-config bar comparison Beamer, 41 pairscan macros (`macros/pairscan/`)  
**Destination:** ARCHIVAR — GEN-1 bug; analysis scripts can be cherry-picked to REF_BASE if needed

---

### feat/ej230-bar-tir-only @ b281aea

**Unique commits vs REF_BASE:** 32  
**Surface:** TIR-ONLY (no reflector; photons failing TIR lost to world)  
**Default scint:** OPSC-106 (EJ-230)  
**Key scripts:** Same pair-scan + resolution-scan set as feat/bar-vikuiti  
**Unique artifacts:** EJ-204 vs EJ-230 TIR-only bar timing comparison Beamer  
**Destination:** REQUIERE DECISIÓN — TIR-only is a physically valid "lower bound" configuration for reflector comparisons; analysis results valid within TIR-only interpretation

---

### feat/ej204-bar-tir-only @ 09f8b18

**Unique commits vs REF_BASE:** 30  
**Surface:** TIR-ONLY  
**Default scint:** OPSC-101 (EJ-204)  
**Key scripts:** `scripts/run_pair_scan.sh`, bar timing resolution analysis  
**Destination:** REQUIERE DECISIÓN — same physics interpretation caveat as EJ-230 variant

---

### feat/ej228-cylinder @ 66b674c

**Unique commits vs REF_BASE:** 3  
**Geometry:** G4Tubs cylinder, radius 12.5 mm, half-height 12.5 mm (∅25 × 25 mm)  
**Material:** EJ-228 (hardcoded, `Materials::CreateEJ228()`)  
**Surface:** 2T-VIK — correct (annular air gap 0.10 mm + Vikuiti annular wrap 0.10 mm)  
**Readout:** 4 top-cap SiPMs (z = +12.5 + ε) + 4 bottom-cap SiPMs (z = −12.5 − ε)  
**Key scripts:** `scripts/run_t0minidaq_endtop_scan_5000.sh`, `scripts/analyze_t0minidaq_endtop_corefit.py`, `scripts/run_validation_scan.sh`, 31 t0minidaq scan macros  
**Destination:** MANTENER VIVA — unique geometry (cylinder), correct physics

---

### feat/ej228-tir-only @ 0006919

**Unique commits vs REF_BASE:** 4  
**Geometry:** Same cylinder as feat/ej228-cylinder  
**Surface:** TIR-ONLY (air gap 0.10 mm, no wrap)  
**Key scripts:** Same t0minidaq scan + analysis set as feat/ej228-cylinder  
**Unique artifacts:** Beamer comparing cylinder TIR-only vs cylinder+Vikuiti  
**Destination:** MANTENER VIVA — unique geometry + systematic comparison with feat/ej228-cylinder

---

### feat/ej230-endonly-mylar @ 98af271

**Unique commits vs REF_BASE:** 41 (most scripts unique to this branch)  
**Surface:** NEAR-GEN-1 (Mylar panels directly on bar, no air gap, `dielectric_metal`)  
**Default scint:** OPSC-106 (EJ-230)  
**Readout:** End only (top-view macros deleted vs main)  
**Key scripts:** `scripts/run_center.sh`, `scripts/run_scan.sh`, `scripts/run_t0minidaq.sh`, `scripts/run_analysis_t0minidaq.sh`, `scripts/generate_exec14b_tables.py`, `scripts/rebuild_exec14b_report.py`, `scripts/audit_beamer_assets.py` (10+ unique analysis scripts)  
**Unique artifacts:** EXEC14b/14d analysis, preflight scripts, EJ-230 Beamer  
**Last commit before W2 session:** User committed `main.tex` changes; W2 currently on `feat/ej230-sslg4` (not this branch)  
**Destination:** MANTENER VIVA — richest EJ-230 End-only analysis; Near-GEN-1 bug affects photon yield but not relative comparisons within branch

---

### feat/ej230-sslg4 @ 5b93f4c

**Unique commits vs REF_BASE:** 33  
**Surface:** GEN-1 bug  
**Default scint:** OPSC-106 (EJ-230)  
**Readout:** End/Top/EndTop  
**Key scripts:** `scripts/run_exec14b_*`, `scripts/run_scan.sh`, `scripts/run_analysis_t0minidaq.sh` (same set as feat/ej230-endonly-mylar but with Top readout)  
**Worktree-pinned in:** W2 (current checkout at `/mnt/d/SHiP/ej230`)  
**Destination:** ARCHIVAR — GEN-1 bug, physically superseded by feat/ej230-endonly-mylar for EJ-230 analysis; worktree-pinned (cannot delete without changing W2 HEAD)

---

### feat/endonly-mylar @ fb3749d

**Unique commits vs REF_BASE:** 23  
**Surface:** GEN-1 bug (skin surface, `fTopSurface = "mylar"` default → `CreateMylarReflector()` as skin on barLV)  
**Default scint:** OPSC-101 (EJ-204)  
**Readout:** End only (+ experimental "sipm" top mode)  
**Key scripts:** `scripts/run_scan.sh`, `scripts/run_msi.sh`, `scripts/run_t0minidaq.sh`  
**Worktree-pinned in:** W1 worktree `/mnt/d/SHiP/ej200_endonly`  
**Destination:** ARCHIVAR — GEN-1 bug; worktree-pinned

---

### feat/endtop-sslg4 @ 5576687

**Unique commits vs REF_BASE:** 2  
Commits:  
- `5576687 fix(optics): eliminate group-velocity aliasing bug in EXEC_23 air-gap geometry`  
- `610b189 feat(validation): physics validation scan — 7 pos × 500 events`  
**Surface:** 2T-MYLAR (correct; uses `CreateMylarReflector()` R=0.95 on air→Mylar boundary; air gap 0.10 mm)  
**Default scint:** OPSC-101 (EJ-204)  
**Key scripts:** `scripts/run_t0minidaq_endtop_scan_5000.sh`, 31 t0minidaq macros, `scripts/analyze_t0minidaq_endtop_*.py`  
**Unique artifacts:** Group-velocity fix commit (historical reference), 7-position validation scan  
**Destination:** MANTENER VIVA (as archive reference) or ARCHIVAR — predecessor to REF_BASE; two unique commits that were NOT cherry-picked into REF_BASE (group-velocity fix and validation scan). Check if group-velocity fix is already absorbed; if not, cherry-pick candidate.

---

### exp/pair-scan-2026-06-11 @ f19c093

**Unique commits vs REF_BASE:** 28  
**Surface:** GEN-1 bug  
**Default scint:** OPSC-101 (EJ-204)  
**Key scripts:** `scripts/gen_pair_scan_macros.py`, `scripts/run_pair_scan.sh`, 41 pairscan macros, `scripts/analyze_geometry.py`  
**Destination:** ARCHIVAR — GEN-1 bug; pair-scan scripts exist in feat/bar-vikuiti as well

---

### feat/ej204-event-display-tracks @ 47a9a4f

**Unique commits vs REF_BASE:** 2  
Commits:  
- `47a9a4f fix(tests): update check_endtop_balance smoke test`  
- (one more — smoke test related)  
**Surface:** GEN-1 bug  
**Key macro:** `macros/ej204_x0_first100_tracks.mac` (unique event display macro)  
**Worktree-pinned in:** W1 worktree `/mnt/d/SHiP/ej200_event_display`  
**Destination:** ARCHIVAR — GEN-1 bug; only unique contribution is event display macro + smoke test fix (cherry-pick candidate); worktree-pinned

---

### feature/sipm-electronics-response @ bb4cf6a

**Unique commits vs REF_BASE:** 15  
**Surface:** 3T-MYLAR (correct; Mylar wrap volumes, bar→Mylar `dielectric_dielectric`, Mylar→air TIR at 37.3°)  
**Default scint:** EJ-200/EJ-204 (configurable, uses `fScintillatorName` not OPSC codes)  
**Readout:** Configurable `fTopSiPMPitch` — different SiPM layout model from all other branches  
**Key feature:** Electronic response simulation (SiPM electronics chain); `fEdgeWrapMode` (Mylar/air/black configurable edge caps)  
**Destination:** REQUIERE DECISIÓN — oldest branch, unique electronics model; 3T-MYLAR surface model may be physically interesting; has not been updated to OPSC code system

---

### wip/host-uncommitted-2026-08-31 @ d2c8d4c

**Unique commits vs REF_BASE:** 2 (cherry counts)  
Commits:  
- `d2c8d4c wip: add EndSparseTop run script`  
- `1f75563 wip: uncommitted DetectorConstruction and SiPMSD changes from host`  
**Surface:** 2T-VIK (same as REF_BASE)  
**Default scint:** OPSC-101 (EJ-204)  
**New readout mode:** `EndSparseTop` — adds `/det/nSparseTopSiPMs N` command, up to 200 uniformly-spaced top SiPMs; `SparseTopSiPMCenterX()` function  
**Defect present:** Defect #1 (gun midpoint uses dense TopSiPMCenterX instead of SparseTopSiPMCenterX); Defect #2 (SiPM placement depth 0.25 mm off in EndSparseTop mode)  
**Destination:** MANTENER VIVA — unique EndSparseTop topology; merge candidate into REF_BASE after defect #1 and #2 are resolved

---

### wip/host-stash-endtop-junio @ 9710f34

**Unique commits vs REF_BASE:** 1  
Commit: `9710f34 wip: recover June stash (CMakeLists + exec07 scan script)`  
**Base commit:** `1f6aca1 feat(exec07): add SSLG4 EndTop 86-channel geometry` (on old `feat/endtop-sslg4` from June 2026, NOT on REF_BASE's ancestry)  
**Surface:** NEAR-GEN-1 — `CreateBarSkinReflector()` (dielectric_metal R=0.98) applied as border surfaces on explicit reflector panels with NO air gap between bar and panel  
**Destination:** ARCHIVAR — Near-GEN-1 bug, single WIP commit on stale base; content (CMakeLists + run_exec07_scan.sh) already exists in REF_BASE

---

## §5 — Redundancy grouping

### Group A — Bar, EJ-204, 2T-VIK (CORRECT baseline)

Most advanced: **REF_BASE** (`feat/bar-end-vikuiti`)  
Others: `wip/host-uncommitted-2026-08-31` (extends with EndSparseTop)  
Action: merge EndSparseTop into REF_BASE after defect fixes

### Group B — Bar, EJ-204, GEN-1 (buggy)

Members: `main`, `exp/pair-scan-2026-06-11`, `feat/ej204-event-display-tracks`  
Most advanced analysis-wise: `exp/pair-scan-2026-06-11`  
All results from this group are affected by GEN-1 (Npe ÷1500 vs CORRECT)

### Group C — Bar, EJ-204, TIR-only (no reflector)

Members: `feat/ej204-bar-tir-only`  
No GEN-1 bug; TIR-only interpretation

### Group D — Bar, EJ-230, 2T-MYLAR

Members: `feat/endtop-sslg4` (Mylar R=0.95 air-gap model)  
2 unique commits vs REF_BASE; predecessor to REF_BASE

### Group E — Bar, EJ-230, GEN-1/Near-GEN-1 (buggy)

Members: `feat/bar-vikuiti` (GEN-1 skin), `feat/ej230-sslg4` (GEN-1 skin), `feat/ej230-endonly-mylar` (Near-GEN-1 panels), `wip/host-stash-endtop-junio` (Near-GEN-1 panels)  
Most analysis content: `feat/ej230-endonly-mylar` (41 unique commits, End-only)

### Group F — Bar, EJ-230, TIR-only (no reflector)

Members: `feat/ej230-bar-tir-only`  
TIR-only interpretation; systematic comparison counterpart to feat/ej230-endonly-mylar

### Group G — Bar, EJ-204, Mylar wrap (3-tier)

Members: `feature/sipm-electronics-response`  
Unique electronics model, different SiPM layout scheme

### Group H — Cylinder, EJ-228

Members: `feat/ej228-cylinder` (2T-VIK correct), `feat/ej228-tir-only` (TIR-only)  
Unique geometry; no overlap with bar branches

### Group I — Bar, EJ-204, EndOnly + Mylar skin (GEN-1)

Members: `feat/endonly-mylar`  
Worktree-pinned; GEN-1 bug

---

## §6 — Worktree-pinned branches

These branches are checked out in an active git worktree or clone HEAD. Deleting them from origin while the worktree is mounted causes `git fetch --prune` to leave dangling heads. **Do not delete without first changing the worktree HEAD.**

| Branch | Pinned by | Location |
|--------|-----------|----------|
| `feat/bar-end-vikuiti` | HOST main checkout + W1 main checkout | `/home/rrios/ej200`, `/mnt/d/SHiP/ej200` |
| `main` | W1 linked worktree | `/mnt/d/SHiP/ej200_edge_scan` |
| `feat/endonly-mylar` | W1 linked worktree | `/mnt/d/SHiP/ej200_endonly` |
| `feat/ej204-event-display-tracks` | W1 linked worktree | `/mnt/d/SHiP/ej200_event_display` |
| `feat/ej230-sslg4` | W2 main checkout | `/mnt/d/SHiP/ej230` |

W1 worktrees discovered via `git worktree list` (SSH, 2026-08-31): four worktrees at `/mnt/d/SHiP/ej200`, `ej200_edge_scan`, `ej200_endonly`, `ej200_event_display`.

---

## §7 — Destination proposals

| Branch | Proposal | Rationale |
|--------|----------|-----------|
| **feat/bar-end-vikuiti** | MANTENER VIVA | Active development branch; REF_BASE; 2T-VIK correct physics |
| **main** | REQUIERE DECISIÓN | GEN-1 bug; worktree-pinned in W1 ej200_edge_scan; must not be the merge target until bug is fixed or acknowledged |
| feat/bar-vikuiti | ARCHIVAR | GEN-1 bug; analysis scripts overlap with REF_BASE; cherry-pick resolution scan if needed |
| feat/ej230-bar-tir-only | REQUIERE DECISIÓN | TIR-only, physically valid for systematic comparison; no GEN-1 bug; 32 unique analysis commits |
| feat/ej204-bar-tir-only | REQUIERE DECISIÓN | Same as above, EJ-204 variant; 30 unique commits |
| feat/ej228-cylinder | MANTENER VIVA | Unique cylinder geometry, correct 2T-VIK physics |
| feat/ej228-tir-only | MANTENER VIVA | Systematic TIR-only counterpart to feat/ej228-cylinder |
| feat/ej230-endonly-mylar | MANTENER VIVA | Richest EJ-230 End-only analysis (41 unique commits); Near-GEN-1 bug affects absolute Npe, not relative trends |
| feat/ej230-sslg4 | ARCHIVAR | GEN-1 bug; worktree-pinned in W2; unmount W2 before deleting |
| feat/endonly-mylar | ARCHIVAR | GEN-1 bug; worktree-pinned W1 ej200_endonly; unmount before deleting |
| feat/endtop-sslg4 | REQUIERE DECISIÓN | 2T-MYLAR correct; 2 unique commits not in REF_BASE (group-velocity fix, validation scan); cherry-pick candidate |
| exp/pair-scan-2026-06-11 | ARCHIVAR | GEN-1 bug; pair-scan scripts also in feat/bar-vikuiti |
| feat/ej204-event-display-tracks | ARCHIVAR | GEN-1 bug; unique content = 1 event display macro + smoke test (cherry-pick before archiving); worktree-pinned |
| feature/sipm-electronics-response | REQUIERE DECISIÓN | Unique electronics model; 3T-MYLAR correct surface; not updated to OPSC codes; 15 unique commits |
| wip/host-uncommitted-2026-08-31 | MANTENER VIVA | EndSparseTop topology unique; merge target after defect resolution |
| wip/host-stash-endtop-junio | ARCHIVAR | Near-GEN-1, stale base, 1 WIP commit; content already in REF_BASE |

**Archive command template** (run after user approval; do NOT execute automatically):
```bash
# Tag before deleting (example):
git tag v-archive-<branch-slug> origin/<branch>
git push origin v-archive-<branch-slug>
git push origin --delete <branch>
```

---

## §8 — Code defects

All defects documented only; none fixed. Evidence from `git show HEAD:src/...` and `git show origin/<wip>:src/...`.

### Defect #1 — Wrong SiPM position function in midpoint gun (EndSparseTop mode)

**File:** `src/PrimaryGeneratorAction.cc:108–109`  
**Branches affected:** `wip/host-uncommitted-2026-08-31`; all branches that add EndSparseTop  
**Description:** The `/muon/midpointSiPMs i j` command guard checks `GetNTopSiPMs()` (dense grid count, 70) then calls `TopSiPMCenterX()` (dense 20 mm pitch). In EndSparseTop mode it should call `SparseTopSiPMCenterX()` instead. Effect: gun X position in EndSparseTop midpoint mode does not correspond to any actual SiPM position. In End-only mode the guard always fails (0 top SiPMs), so the bug is inoperative and the fallback x=0 has no physical impact (coincides with SiPMs 34/35 midpoint).

### Defect #2 — EndSparseTop SiPM placement depth off by 0.25 mm

**File:** `src/DetectorConstruction.cc:319` (wip/host-uncommitted-2026-08-31)  
**Branches affected:** `wip/host-uncommitted-2026-08-31`  
**Description:** Sparse top SiPMs placed at `kBarHalfY - 2.0*kTopHalfY` (= 30.0 − 0.5 = 29.5 mm from centre) versus `kBarHalfY - kTopHalfY` (29.75 mm) for dense top. The 0.25 mm difference shifts the SiPM into the bar by one half-thickness extra, affecting TIR acceptance angle for photons heading toward top SiPMs. Pending decision on whether this placement was intentional (recessed SiPM coupling) or a bug.

### Defect #3 — Analysis scripts hardcode SiPM globalId range TOP = 16–85

**Files:** `analysis/congruent_sum4_timing.C:224,308`, `analysis/grouped_resolution.C:46`, `analysis/preliminary_position_scan.C:6,183`, `analysis/high_stats_position_scan.C:243`, `analysis/grouped_resolution.py:61`, `analysis/exec07_photon_budget.py:250,419,570,579`, `analysis/topreadout_crosstalk.py:26`  
**Branches affected:** All branches sharing this analysis/ directory (all bar branches)  
**Description:** `kNTotalSiPMs = 86` is a named constant in `DetectorConstruction.hh` but is not used at runtime; `SetNSparseTopSiPMs` accepts up to 200 SiPMs; in EndSparseTop runs globalIds can exceed 85. Scripts that gate on `globalId >= 16 && globalId <= 85` will silently misclassify sparse-top SiPMs with id > 85 as neither End nor Top. Affects only EndSparseTop runs.

### Defect #4 — GEN-1 bug on main branch (physics validity)

**File:** `src/DetectorConstruction.cc:153–155` (on `origin/main`)  
**Branches affected:** `main`, `feat/bar-vikuiti`, `feat/ej230-sslg4`, `exp/pair-scan-2026-06-11`, `feat/ej204-event-display-tracks`, `feat/endonly-mylar`  
**Description:** `G4LogicalSkinSurface("BarSkin", barLV, dielectric_metal R=0.98)` applied directly to the bar logical volume. This prevents Geant4 from applying Fresnel equations at the bar-air interface, eliminating TIR (θ_c = 39.3° for n=1.58). Photon yield is reduced by factor ~1500 vs CORRECT two-tier model. Results from these branches are physically invalid for absolute photon-yield comparisons. Relative timing-resolution comparisons within a single branch remain self-consistent if consistently re-scaled.

---

## §9 — Evidence references

All commands are read-only; no files modified, no branches altered.

```bash
# Physical matrix: source of truth
git show origin/<branch>:src/DetectorConstruction.cc | grep -En '...'
git show origin/<branch>:src/Materials.cc | grep -A20 'CreateBarSurface\|CreateMylarReflector\|CreateBarSkinReflector'

# Cherry counts
git cherry -v origin/feat/bar-end-vikuiti origin/<branch> | grep '^+' | wc -l

# Diff stats
git diff --stat "origin/main...origin/<branch>" -- src/ include/ CMakeLists.txt macros/ scripts/

# SiPM ID range in analysis
grep -rn '16\b.*85\b\|globalId.*85\|range(16.*86)' analysis/

# W1 worktrees
ssh -p 9022 reriosto@127.0.0.1 bash --noprofile --norc -c 'git -C /home/reriosto/SHiP/ej200 worktree list'
```

---

## §10 — Checklist status

| Criterion | Status |
|-----------|--------|
| 3 bundles + HOST, `bundle verify` green, SHA-256 recorded, different devices | ✓ DONE |
| Zero stashes in all 3 clones | ✓ DONE |
| `feat/ej228-*` confirmed on origin | ✓ DONE (66b674c, 0006919) |
| Every matrix cell has value or INDETERMINADO | ✓ DONE |
| Every destination proposal cites evidence | ✓ DONE |
| Zero tags created, zero branches deleted, zero merges | ✓ DONE (verifiable) |
| Code defects documented (not fixed) | ✓ DONE |
