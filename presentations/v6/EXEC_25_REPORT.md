# EXEC_25_REPORT.md — Internal Consistency Campaign for talk_v6.tex

Generated: 2026-09-02
Authoritative clone: t0minidaq:/home/rrios/ej200, branch main
Rollback tag: pre-exec25-260902 = 1c5a4c5
HEAD before campaign: 1c5a4c5 (38 pages, EXEC_24 result)
HEAD after Phase 1: ba00e0e

---

## §1 — PDF Provenance

| Artifact | sha256 | Bytes | Pages | Built from |
|----------|--------|-------|-------|------------|
| Reviewer PDF (EXEC_24 build) | e47e2014a62a22782455e3a1a5619207108da7504d9c10e6eea56daaf275512b | 472 996 | 38 | 1c5a4c5 |
| EXEC_25 Phase 1 build | 15e2e0beeccead87935fe51214bb5861f1d2d748fe5e179324fcf378573c8550 | 472 856 | 38 | ba00e0e |
| EXEC_25 Phase 1b build (after FIX-23b) | ce20abef20cd3b316f8f0e9eccf4740d055ab85a10f8be8d9a09590d0ba62d73 | 472 856 | 38 | 9f8a3a9 |
| EXEC_25 Phase 2+3 build (FIX-23c, FIX-18b, Phase 2 figs) | 6d6d9fa8257872fef4482e7ce5cc2c991d417ef88eeceac1ed864128b824cbb0 | 472 656 | 38 | 6488075 |
| **EXEC_25 Phase 4 final build (FIX-26)** | **3f5f43c0f640b5a7351babf648942bdeb425aed73eed2d67a2dd138226b113cc** | **472 658** | **38** | **7f68624** |

PDF built on MSI exec25-deckfix branch, latexmk -lualatex.
Phase 4 PDF committed at HOST main c93c833.

---

## §2 — CHECKPOINT 0 Actions (Phase 0)

| Item | Action | Outcome |
|------|--------|---------|
| Stash@{0} on MSI (exec24-working-tree-before-sync) | Dropped | Byte-identical to 1c5a4c5 (empty diff verified); sha fea60d20 |
| Branch exec24-deckfix on MSI | Deleted | Fully contained in exec24-sync at 1c5a4c5 |
| W1 (/home/reriosto/SHiP/ej200) tex sha256 | Verified | 064e6df2... ✓ matches 1c5a4c5 |
| W2 search (find / -maxdepth 4 -type d -name 'ej200*') | No additional clones found | Only 4 known /mnt/d/SHiP/ej200* worktrees |
| Untracked .pptx: SHiP_timing_first_principles.pptx | Identified; NOT added to .gitignore | Manually created PowerPoint (575 KB, 2026-08-19); not a generated prototype; excluded from all commits |
| Rollback tag pre-exec25-260902 | Created | = 1c5a4c5 |
| MSI branch topology (Option A) | exec25-deckfix created from exec24-sync FETCH_HEAD | MSI main stays at 5e90025 in ej200_edge_scan; untouched |

---

## §3 — FIX-16..25 Status

| FIX | Status | Commit | Files | Before → After |
|-----|--------|--------|-------|----------------|
| FIX-16 | DONE | d6e8794 | talk_v6.tex | C1 table: deleted "Napkin/G4 ratio N_pe" row (0.807/0.718, 0.651/0.535, napkin EJ-230 only) and its preceding \midrule. Mixed-estimator diagnostic row not used in any design conclusion. |
| FIX-17 | DONE | 832ca55 | talk_v6.tex | G2 qualitative-BLUE-check: `$w_E \propto \sigma_T^2$ in BLUE` → `$w_E$ monotone increasing in $\sigma_T^2$ in BLUE`. w_E is monotone increasing in σ_T² (not proportional; denominator also varies). |
| FIX-18 | DONE | f885c97 | talk_v6.tex | S17 oracle-gap text line 680: "arise because k*=2 is better there but global uses k=3" → "x=±200 mm (k*=2) and x=±500 mm (k*=1): both below global k=3". Fixes wrong claim that k*=2 at x=±500 (oracle shows k*=1 there). |
| FIX-19 | DONE | 92aab49 | talk_v6.tex | S20 summary \alert{}: "Achieved: 15.2 ps — 6.6× better than target" → "Achieved: 15.2 ps intrinsic (see App.~I2) — 6.6× better than target". Clarifies that 15.2 ps is simulation-only (no SPTR convolution). |
| FIX-20 | DONE | d16a584 | talk_v6.tex | S18 Pareto caption: "20 TOP channels for last 3.4 ps gain" → "20 TOP channels for last 3.4 ps gain (unrelated to the S17 oracle gap)". Disambiguates from the S17 oracle-vs-global 3.4 ps gap (same number, different physics). |
| FIX-21 | DONE | 47a42e4 | talk_v6.tex | S12/S13: three ~13 mm → ~10 mm. (a) Line 499: removed "(path ≈ 13 mm)" parenthetical (redundant with ~10 mm already stated). (b) Line 524 table: ~13 mm → ~10 mm. (c) Line 552 tikzpicture label: ~13 mm → ~10 mm. 10 mm = bar height H. |
| FIX-22 | DONE | c029726 | talk_v6.tex | \date{August 2026} → \date{September 2026}. Page 1 verified. |
| FIX-23 | DONE | 660bbc9 | talk_v6.tex | 10 occurrences of literal 0.95 → \Rvikuiti: napkin photon-loss item, optical-surface listing, design-summary bullet, sim-model R item, A2 heading, B1 GEN config table (×2), exec21-optfix note, J2 dataset table. Skipped: line 277 computed exponent ($= 0.95^{85.6}$); macro definition line 33. |
| FIX-23b | DONE | 9f8a3a9 | talk_v6.tex | Missed literal caught by grep check: line 277 `$= 0.95^{85.6}$` → `$= \Rvikuiti^{85.6}$`. After this, grep returns exactly one hit (line 33 macro definition). |
| FIX-23c | DONE | 84b4dde | talk_v6.tex | Three occurrences that were originally text-mode were incorrectly wrapped in `$...$` by FIX-23, adding math-italic R and thin spaces around `=`, causing pdftotext reflow on pages 21, 23, 25. Reverted to text-mode `R=\Rvikuiti` (lines 821, 888, 942). Pages 21, 23, 25 now EMPTY in gate. |
| FIX-24 | DONE | 47d0a2f | FINAL_NUMBERS.md, README.md | FINAL_NUMBERS.md: (1) Λ_refl(R=0.95) ~240 mm → 159.4 mm (EXEC_23 FIX-05 source); (2) σ_END END-only row 47.6/52.1/49.6 → 54.86/53.36/50.49 ps (EXEC_23 FIX-01 source, m*-opt x=0 from results_macros.tex). README.md: add "MSI branch topology" paragraph (Option A decision). |
| FIX-25 | DONE | ba00e0e | EXEC_23_REPORT.md, CONFIGURATION_AUDIT.md | EXEC_23_REPORT.md FIX-03 row: note G1/G2 "combined frame" was split into G1+G2 by EXEC_24 FIX-11. CONFIGURATION_AUDIT.md BLUE weights §: clarify formula role — in w_E=(σ_T²−C_ET)/(σ_E²+σ_T²−2C_ET), σ_T² appears in numerator and denominator; σ_E² appears in denominator only. |
| FIX-18b | DONE | 8073eec | talk_v6.tex | S17 oracle-gap line 680: "both below global k=3" → "global k=3 is suboptimal there". "Both below" misread as k ordering statement; new wording correctly says k=3 is suboptimal at ±200 and ±500 mm. |
| FIX-26 | DONE | a09b7bc | talk_v6.tex | S13 "Why TOP beats END": `~100 pooled photons` → `~1300 pooled photons`. Source: top_npe_diag.csv x=0: N_pe^TOP mean=1341, median=1242. F1 d_near sentence also added but reverted (see FIX-26b). |
| FIX-26b | DONE | 7f68624 | talk_v6.tex | Revert F1 appendix d_near sentence: frame overflow 12.3 pt (Overfull \\vbox at line 1103). Per CHECKPOINT 2/3: "ONLY IF it fits without overflow, otherwise skip." |

---

## §4 — Commit Log (pre-exec25-260902..HEAD)

```
c93c833 EXEC_25 Phase 4: commit final PDF — talk_v6.pdf at 7f68624
7f68624 EXEC_25 FIX-26b: revert F1 d_near sentence — frame overflows by 12.3pt
41c94c5 EXEC_25: track SHiP_timing_first_principles.pptx in napkin_first_principles/
8765cd0 EXEC_25 Phase 2: add CSV + meta.json sidecars for v5_veff_fit and v5_top_position_loo
3b4cc31 EXEC_25 Phase 3 addendum: top_npe_diag — add hybrid vs END-only N_pe^END comparison
a09b7bc EXEC_25 FIX-26: S13 ~100→~1300 pooled photons; F1 add d_near explanation sentence
c449fdc EXEC_25: update EXEC_25_REPORT.md — Phase 1b, 2, 3 status and gate results
6488075 EXEC_25 Phase 3: top_npe_diag.py — N_pe^TOP and N_pe^END diagnostic (CSV + meta.json)
26d7355 EXEC_25 Phase 2: fix v5 figure encoding bugs (mojibake title; v_eff_err precision)
8073eec EXEC_25 FIX-18b: S17 oracle-gap — replace "both below global k=3" with "global k=3 is suboptimal there"
84b4dde EXEC_25 FIX-23c: keep originally text-mode R=\Rvikuiti in text mode (S21, A2, B1)
9f8a3a9 EXEC_25 FIX-23b: replace missed literal 0.95 in survival formula (line 277)
249869f EXEC_25: add EXEC_25_REPORT.md — Phase 1 campaign report
ba00e0e EXEC_25 FIX-25: fix formula descriptions in EXEC_23_REPORT and CONFIGURATION_AUDIT
47d0a2f EXEC_25 FIX-24: update stale values in FINAL_NUMBERS.md; add MSI topology note to README
660bbc9 EXEC_25 FIX-23: replace 10 occurrences of literal R=0.95 with \Rvikuiti macro
47a42e4 EXEC_25 FIX-21: S12/S13 — TOP path label ~13 mm → ~10 mm (bar height)
d16a584 EXEC_25 FIX-20: S18 Pareto — note 3.4 ps channel gain is unrelated to S17 oracle gap
92aab49 EXEC_25 FIX-19: S20 summary slide — label 15.2 ps as intrinsic, add App.~I2 ref
f885c97 EXEC_25 FIX-18: S17 oracle-gap text — specify k*=2 at x=±200 mm and k*=1 at x=±500 mm
d6e8794 EXEC_25 FIX-16: C1 — remove Napkin/G4 ratio N_pe row from material table
832ca55 EXEC_25 FIX-17: G2 qualitative-BLUE-check — \propto → monotone-increasing-in
c029726 EXEC_25 FIX-22: update date August → September 2026
```

---

## §5 — Overflow Inventory

### Baseline (pre-exec25-260902 = 1c5a4c5, 38 pages)

| Log line | Type | Size | Status |
|----------|------|------|--------|
| 953 | Underfull \hbox | badness 10000 | Pre-existing; B1 footnote paragraph end |

### EXEC_25 Phase 1 (ba00e0e, 38 pages)

| Log line | Type | Badness | Cause |
|----------|------|---------|-------|
| 953 (×3) | Underfull \hbox | 10000, 1270, 1270 | Pre-existing; B1 footnote paragraph — same location, paragraph now reflows into 3 lines each slightly underfull. All mild. |

No Overfull \vbox or \hbox warnings. **CLEAN.**

### EXEC_25 Phase 2+3 build (6488075, 38 pages)

No Overfull \vbox or \hbox warnings. **CLEAN.** (Underfull at B1 footnote expected to persist.)

### EXEC_25 Phase 4 final build (7f68624, 38 pages)

No Overfull \vbox or \hbox warnings. **CLEAN.**
(FIX-26b reverted the F1 sentence that caused 12.3 pt overflow at line 1103 in the intermediate build a09b7bc.)

---

## §6 — MUST-NOT-CHANGE Gate

Named constants verified present (as literals or via macros) in talk_v6.tex after all FIX commits.
FIX-23 replaced literal `0.95` with `\Rvikuiti` (expands to `0.95`; rendered value unchanged).
FIX-21 changed `13` → `10` (13 mm was not a named constant; 10 mm = bar height, physically correct).
FIX-26 changed `~100` → `~1300` (`~100` was not a named constant; `~1300` is derived from N_pe^TOP mean=1341 at x=0, top_npe_diag.csv).
No named constant values were altered.

---

## §7 — Visual Inspection (exec25-deckfix build, 80 dpi)

| Page | Frame | FIX verified | Verdict |
|------|-------|-------------|---------|
| 1 | Title | FIX-22 | "September 2026" ✓ |
| 12 | TOP Estimator | FIX-21 | "~10 mm", τ_END ∈ [3.689, 5.829 ns] ✓ |
| 13 | Near-Field vs Far-Field | FIX-21 | table "~10 mm", diagram label "~10 mm" ✓ |
| 17 | Oracle vs Global | FIX-18 | "x=±200 mm (k*=2) and x=±500 mm (k*=1): both below global k=3" ✓ |
| 18 | N_TOP Scan / Pareto | FIX-20 | "3.4 ps gain (unrelated to the S17 oracle gap)" in caption ✓ |
| 20 | Design Decision | FIX-19, FIX-23 | "15.2 ps intrinsic (see App. I2)" ✓; R=\Rvikuiti renders 0.95 ✓ |
| 27 | C1 Material Table | FIX-16, FIX-23 | Napkin/G4 ratio row absent; Napkin bulk survival = last row ✓ |
| 33 | G2 Covariance | FIX-17 | "w_E monotone increasing in σ_T²" ✓; all caveat text intact ✓ |

---

## §8 — MSI Worktree Status (Phase 4)

```
/mnt/d/SHiP/ej200                → exec25-deckfix @ 7f68624 (Phase 4 final, FIX-26b)
/mnt/d/SHiP/ej200_edge_scan      → main @ 5e90025 (frozen pre-EXEC_23; do not touch)
/mnt/d/SHiP/ej200_endonly        → feat/endonly-mylar @ fb3749d
/mnt/d/SHiP/ej200_event_display  → feat/ej204-event-display-tracks @ 47a9a4f
/tmp/exec25_old                  → detached HEAD @ 1c5a4c5 (cleanup pending)
W1 /home/reriosto/SHiP/ej200     → exec25-deckfix @ 7f68624 (main repo on MSI; ff'd to Phase 4)
```

tex SHA-256 verified identical on HOST and MSI at 7f68624: `7ff0ce972fdb3c02b992442019d440542661a95f3d35332f0a4bf7962e80f808`

---

## §8b — CHECKPOINT 1 Gate (pdftotext, pre-exec25-260902 vs 8073eec)

Pages tested: 1, 4, 5, 6, 9, 12, 13, 17, 18, 20, 21, 23, 24, 25, 26, 27, 33, 37.

| Page | Result | Content |
|------|--------|---------|
| 1 | DIFFERS — expected | FIX-22: `August 2026` → `September 2026` |
| 4,5,6,9 | EMPTY ✓ | |
| 12 | DIFFERS — expected | FIX-21: removed `(path ≈ 13 mm)` + column reflow |
| 13 | DIFFERS — expected | FIX-21: `∼ 13 mm` → `∼ 10 mm` (×2 in table) |
| 17 | DIFFERS — expected | FIX-18+18b: sentence replaced on oracle-gap peaks |
| 18 | DIFFERS — expected | FIX-20: `(unrelated to the S17 oracle gap)` added |
| 20 | DIFFERS — expected | FIX-19: `intrinsic (see App. I2)` + column reflow |
| 21,23,24,25,26 | EMPTY ✓ | FIX-23c corrected text-mode spacing |
| 27 | DIFFERS — expected | FIX-16: `Napkin/G4 ratio N_pe` row absent |
| 33 | DIFFERS — expected | FIX-17: `monotone increasing in σ_T²` |
| 37 | EMPTY ✓ | |

Gate PASSES. Check A and B PASS.

---

## §8c — Phase 2: Figure Encoding Fixes

| Figure | Bug | Fix | Gate |
|--------|-----|-----|------|
| v5_top_position_loo.pdf | UTF-8 em dash `—` in ROOT TLatex → mojibake `â€"` in title | Replace `—` with `--` (ASCII) in DrawLatex call | pdftotext diff: title line only + cosmetic column-position shifts ✓ |
| v5_veff_fit.pdf | `v_eff_err` formatted `:.0f` → `± 0 mm/ns` | Change to `:.2f` → `± 0.02 mm/ns` | pdftotext diff: uncertainty line only + cosmetic column shifts ✓ |

Source: `presentations/v5/scripts/analysis_v5.py` (Sections 1 and 4, run as scratchpad/fix_phase2_figs.py).
Files replaced in-place: `presentations/v6/figs/v5_top_position_loo.pdf`, `presentations/v6/figs/v5_veff_fit.pdf`.

**Phase 2 sidecars (commit 8765cd0):**
CSV + meta.json sidecars added for both figures, documenting binning, source ROOT path, SHA-256, and generating command.
- `v5_veff_fit.csv`: x, dt_median_ps, dt_err_ps, residual_ps, n_ev (13 rows)
- `v5_veff_fit_meta.json`: fit params (a=0.44 ps, b=−11.259 ps/mm, v_eff=177.6 mm/ns), SHA-256 of all 13 ROOT files
- `v5_top_position_loo.csv`: x, sigma_lin_mm, bias_lin_mm, sigma_cub_mm, bias_cub_mm (13 rows; mean σ_x linear=17.2 mm)
- `v5_top_position_loo_meta.json`: TOP SiPM positions, calibration method, SHA-256 of all 13 ROOT files
Generated by `gen_v5_fig_sidecars.py` (scratchpad); uproot+numpy replication of analysis_v5.py §1 and §4.

---

## §8d — Phase 3: top_npe_diag Diagnostic

Script: `analysis/top_npe_diag/top_npe_diag.py`
Outputs: `analysis/top_npe_diag/top_npe_diag.csv`, `analysis/top_npe_diag/top_npe_diag_meta.json`
Gate: N_pe^END/face(x=0) = 701.28 ∈ 701.3 ± 1 ✓

**Data sources:**
- N_pe^TOP (pooled, all 20 SiPMs): `scan_end_vik_sparse_top_v2`, EJ-230, ntop=20, face_type=2, raw hit count
- N_pe^END (per-face mean): `scan_end_vikuiti`, EJ-230 END-only, face_type 0+1, raw count / 2

**Phase 3 table (x, N_pe^END/face, N_pe^TOP mean/median/IQR-σ, d_nearest, d_second):**

| x [mm] | N_pe^END/face | N_pe^TOP mean | N_pe^TOP median | N_pe^TOP IQR-σ | d_near [mm] | d_sec [mm] |
|--------|--------------|---------------|-----------------|----------------|-------------|------------|
| -600 | 947.23 | 1101.00 | 1021.50 | 160.12 | 5.7 | 64.3 |
| -500 | 859.15 | 1182.20 | 1093.00 | 173.46 | 24.3 | 45.7 |
| -400 | 799.43 | 1257.98 | 1169.00 | 181.62 | 15.6 | 54.3 |
| -300 | 751.66 | 1320.84 | 1223.00 | 186.81 | 14.4 | 55.6 |
| -200 | 725.37 | 1344.32 | 1251.00 | 197.92 | 25.6 | 44.4 |
| -100 | 711.09 | 1419.35 | 1317.00 | 203.11 | 4.5 | 65.5 |
|  **0** | **701.28** | 1341.05 | 1242.00 | 192.92 | 34.5 | 34.5 |
| +100 | 705.63 | 1412.99 | 1314.00 | 200.89 | 4.5 | 65.5 |
| +200 | 725.07 | 1345.91 | 1247.00 | 199.41 | 25.6 | 44.4 |
| +300 | 751.04 | 1319.78 | 1225.00 | 190.51 | 14.4 | 55.6 |
| +400 | 795.17 | 1253.48 | 1166.00 | 177.91 | 15.6 | 54.3 |
| +500 | 863.40 | 1174.58 | 1086.00 | 167.53 | 24.3 | 45.7 |
| +600 | 951.15 | 1097.09 | 1020.00 | 163.83 | 5.7 | 64.3 |

N_pe^END/face is highest near the END faces (x=±600) because photons have a shorter path to END. N_pe^TOP peaks near x=±100 mm where the four nearest TOP SiPMs (d_near=4.5 mm) strongly intercept — but dips at x=0 (d_near=34.5 mm). N_pe^END at x=0 (701.28) reproduces the FINAL_NUMBERS canonical value 701.3 ✓.

**Phase 3 addendum (commit 3b4cc31):** Extended top_npe_diag.py and CSV to include hybrid vs END-only comparison.

New columns: npe_end_hybrid (END hits/face from hybrid scan, N_TOP=20), npe_end_endonly (alias for npe_end_per_face), delta_pe (hybrid − endonly per face per event), delta_pct (%).

**x=0 addendum row:**

| Quantity | Value |
|----------|-------|
| npe_end_endonly (scan_end_vikuiti) | 701.28 hits/face/event |
| npe_end_hybrid (scan_end_vik_sparse_top_v2) | 514.94 hits/face/event |
| delta_pe | −186.34 hits/face/event |
| delta_pct | −26.6 % |

The 20 TOP SiPMs absorb ~26.6 % of photons that would otherwise reach the END faces at x=0. Photon-count half of OPEN-06.

---

## §9 — Open Items Carried Forward

All OPEN-01..06 items from EXEC_23_REPORT.md remain open. No new open items from EXEC_25 Phase 1.

**OPEN-07:** SHiP_timing_first_principles.pptx — **RESOLVED** (commit 41c94c5). Moved to `presentations/napkin_first_principles/` and tracked. README.md updated with one-line entry: "Manual PowerPoint prototype (2026-08-19), superseded by talk_v3.tex / talk_v6.tex."

---

## §10 — Synchronization State (Phase 4)

| Location | Branch | HEAD | Ahead of origin |
|----------|--------|------|----------------|
| HOST main | main | c93c833 | 37 commits |
| MSI exec25-deckfix (/home/reriosto/SHiP/ej200) | exec25-deckfix | 7f68624 | — |
| MSI W1 (/mnt/d/SHiP/ej200) | exec25-deckfix | 7f68624 | — |
| GitHub origin/main | main | 5e90025 | — |

tex SHA-256 HOST = MSI = `7ff0ce972fdb3c02b992442019d440542661a95f3d35332f0a4bf7962e80f808` (at 7f68624).

**Phase 4 push commands (PRINT ONLY — do not execute):**
```bash
# 1. MSI cleanup
git -C /home/reriosto/SHiP/ej200 stash drop
git -C /home/reriosto/SHiP/ej200 branch -d exec24-sync

# 2. W1 fast-forward (after push)
git -C /mnt/d/SHiP/ej200 fetch origin && git -C /mnt/d/SHiP/ej200 merge --ff-only origin/main

# 3. Push HOST main to GitHub
git push origin main

# 4. Push rollback tag
git push origin pre-exec25-260902
```

No push executed. No merge. No rebase. No branch deleted.
