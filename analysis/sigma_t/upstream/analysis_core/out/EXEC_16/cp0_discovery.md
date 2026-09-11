# CP0 — Discovery Report · EXEC_16 · 2026-06-20

## 1. Git status

| Item | Result |
|------|--------|
| Branch | `feat/endtop-sslg4` ✓ |
| Modified tracked files | **None** — working tree clean ✓ |
| Untracked files | `GROUP_VELOCITY_AUDIT.md`, `runs/` — not analysis artifacts, no action needed |

---

## 2. Path verification

| Constant | Declared path | Status |
|---------|--------------|--------|
| `REPO_DIR` | `/home/reriosto/SHiP/ej200` | EXISTS ✓ |
| `ORCH_DIR` | `/home/reriosto/SHiP/orchestrator` | EXISTS ✓ |
| `DATA_DIR` | `/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959` | EXISTS ✓ |

**⚠ Path note:** `datasets.py` lives at `{ORCH_DIR}/analysis/datasets.py`, not at `{ORCH_DIR}/datasets.py` as declared in the constant comment. Actual path: `/home/reriosto/SHiP/orchestrator/analysis/datasets.py`.

### ROOT files inventory

- Simulation ROOT files: **31** in `outputs/x{pos}mm/photon_hits_run000.root`
- Pre-existing analysis files: 2 (in `analysis_corefit/root/` and `analysis_simple_std/root/`)
- Total disk: 1.9 GB
- File naming pattern: `x{N}mm` where N ∈ signed integers (e.g. `x-100mm`, `x0mm`, `x100mm`)

---

## 3. datasets.py — registration check

Dataset `t0minidaq_endtop_scan_20260618_204959` is **NOT registered** in `DATASETS`.

Registered datasets found:
```
"exec07_endtop_2000":  material=EJ-204, readout=ENDTOP, tau_d=1.8ns, tau_r=0.5ns,
                        n_channels=86, opsc_code="opsc-102",
                        data_dir=/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000
"ej230_endtop":        material=EJ-230, readout=ENDTOP, tau_d=1.5ns, tau_r=0.5ns,
                        n_channels=86, opsc_code="opsc-106"
```

This is the **known 2×2 gap** (EJ-204 EndTop scan from 2026-06-18 not yet registered).
Analysis will read `DATA_DIR` directly; this gap is recorded in provenance.

**Ancillary note:** `exec07_endtop_2000` declares `opsc_code="opsc-102"` but
`DetectorConstruction.cc` line 71 maps EJ-204 → OPSC-101 (not OPSC-102). This is
an inconsistency in the old registry entry, not in the current dataset. Recorded here
for future audit; does not affect EXEC_16.

---

## 4. TTree schema — FULL MATCH

File probed: `outputs/x0mm/photon_hits_run000.root`

| Field | Expected (`TREE_NAME`/`BR`) | Found | Match |
|-------|--------------------------|-------|-------|
| Tree name | `sipm_hits` | `sipm_hits` | ✓ |
| `event_id` | int | int32 | ✓ |
| `face_type` | int | int32 | ✓ |
| `global_id` | int | int32 | ✓ |
| `local_id` | int | int32 | ✓ |
| `time_ns` | float | float64 | ✓ |
| `energy_eV` | float | float64 | ✓ |
| `wl_nm` | float | float64 | ✓ |
| `pde` | float | float64 | ✓ |
| `x_mm` | float | float64 | ✓ |
| `y_mm` | float | float64 | ✓ |
| `z_mm` | float | float64 | ✓ |
| `gun_x_mm` | float | float64 | ✓ |

All 12 branches present and types match. Entries at x=0mm: **588,443 hits**.

---

## 5. Material — OPSC-101 / EJ-204 (mac file + DetectorConstruction overrides)

Source: `/home/reriosto/SHiP/ej200/src/external/SSLG4/macros/oscnt/opsc-101.mac`
and `/home/reriosto/SHiP/ej200/src/DetectorConstruction.cc` (overrides at lines 200–215).

| Property | Read from repo | DS_EJ204 (cross-check) | Gate |
|---------|---------------|----------------------|------|
| SCINTILLATIONYIELD | 10400 /MeV | 10400 /MeV | ✓ OK |
| SCINTILLATIONTIMECONSTANT1 (τ_d) | **1.8 ns** | 1.8 ns | ✓ OK |
| SCINTILLATIONRISETIME1 (τ_r) | **0.7 ns** (overridden by DetConstr) | 0.7 ns | ✓ OK |
| ABSLENGTH | **160 cm** (overridden by DetConstr) | 160 cm | ✓ OK |
| Emission peak (SCINTCOMP1 data file) | not extracted here; file referenced | 408 nm | — |

**τ_d = 1.8 ns read from mac, confirmed by DetectorConstruction.cc — not inherited from EJ-230.**

### Wrapping note

`WRAPPING = "Mylar"` declared in constants. However, `DetectorConstruction.cc` line 89:
> "Legacy no-op. The Mylar wrap volume was removed; reflection is handled by a
> reflector skin surface on BarLV."

The physical Mylar wrap *volume* does not exist in the current geometry. Reflection at
bar edges is handled by the skin surface on BarLV (result of the reflector skin fix,
OPSC-106/branch `feat/endtop-sslg4`). The WRAPPING="Mylar" constant is a label
describing the intended physics (Mylar-like reflectivity), not a physical volume.
**Recommendation:** clarify label to "reflector_skin_surface (Mylar-equivalent)" in
captions. Not an abort — recorded for provenance.

---

## 6. Position mapping (gun_x_mm per file)

31 positions, all with **5000 events**. Mapping is unambiguous (one gun_x_mm value per file).

| Positions [mm] | Count |
|----------------|-------|
| −690, −670, −650, −600, −550, −500, −450, −400, −350, −300, −250, −200, −150, −100, −50 | 15 |
| 0 | 1 |
| +50, +100, +150, +200, +250, +300, +350, +400, +450, +500, +550, +600, +650, +670, +690 | 15 |

Range: [−690, +690] mm ✓  
Total events in scan: **31 × 5000 = 155,000 events**.  
Stated scan: "31 positions, x ∈ [−690,+690], paso ~46 mm" — average step = 1380/30 = **46 mm** ✓  
(Steps are 50 mm uniform except densification at ±650→±670→±690.)

---

## 7. ⛔ ABORT CONDITION — face_type / global_id mismatch

Probed at x=0mm (representative):

| `face_type` | global_id range | N unique gid | N hits (x=0mm) |
|------------|----------------|--------------|----------------|
| 0 | 0..7 | 8 | 1870 |
| 1 | 8..15 | 8 | 1890 |
| 2 | 16..85 | **70** | 584,683 |

**Expected per §0 constants vs found:**

| Constant | Declared | Found | Status |
|---------|---------|-------|--------|
| `END_IDS` | `range(0, 16)` — 16 channels | 16 END channels (gid 0..15, split across face_type=0 and =1) | Correct range, but split into 2 face_types |
| `TOP_IDS` | `range(16, 36)` — **20 channels** | gid 16..85 — **70 channels** | ❌ MISMATCH |
| face_type values | implied: {END, TOP} (2 values) | {0, 1, 2} (3 values) | ❌ MISMATCH |
| Total channels | 36 | **86** | ❌ MISMATCH |

**Interpretation (consistent with exec07_endtop_2000 → n_channels=86, ej230_endtop → max_gid=85):**
- face_type=0 → END_LEFT (gid 0..7, 8 SiPMs)
- face_type=1 → END_RIGHT (gid 8..15, 8 SiPMs)
- face_type=2 → TOP (gid 16..85, 70 SiPMs)
- Total: 86 channels (consistent with existing registry)

`TOP_IDS = range(16, 36)` in the §0 constants block is **wrong for this geometry**.
The correct range is `range(16, 86)`.

**⛔ This is an ABORT per CP0 rules.** The `END_IDS`/`TOP_IDS`/face_type logic in §0 does
not match the actual data. Proceeding with wrong `TOP_IDS` would silently discard 50 of 70
TOP channels in every analysis. I am stopping here and requesting your decision on:

1. Correct `TOP_IDS` to `range(16, 86)` (70 TOP channels)?
2. How to handle `face_type`: use integer values {0, 1, 2} with
   `END_FACE_TYPES = {0, 1}` and `TOP_FACE_TYPE = 2`?
3. Does face_type=0 vs face_type=1 distinguish END_LEFT from END_RIGHT, or
   some other physical meaning? (Needed for the Δt = t_endL − t_endR spatial resolution.)
4. "TOP nearest" (max `<Npe>`) over 70 channels — confirm this is the intended default.

---

## 8. Hardware note (open question, no decision)

Clustering: `id//4` fixed vs moving window of 4 neighbors — decision of Gerardo.
For this analysis: nearest = channel with max `<Npe>` at the beam position.
Noted in provenance.

---

## Summary: STOP conditions for CP0

| # | Item | Status |
|---|------|--------|
| 1 | Git clean + correct branch | ✓ |
| 2 | DATA_DIR exists, 31 ROOT files | ✓ |
| 3 | datasets.py — scan NOT registered (known gap) | ✓ noted |
| 4 | TTree schema — all branches match | ✓ |
| 5 | Material OPSC-101/EJ-204: τ_d=1.8ns, τ_r=0.7ns, λ=160cm | ✓ |
| 6 | Position mapping: 31 positions, 5000 evt/pos, gun_x constant per file | ✓ |
| 7 | face_type {3 values} + TOP_IDS range [16,36]≠[16,85] | **⛔ ABORT — awaiting decision** |
