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

tex sha256 (HOST and MSI identical): `448917693bd6c92ec219e32be88731a2105c0dbbca197b8949ac748466f604d4`

PDF built on MSI exec25-deckfix branch, latexmk -lualatex, mtime 18:10:12.

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
| FIX-24 | DONE | 47d0a2f | FINAL_NUMBERS.md, README.md | FINAL_NUMBERS.md: (1) Λ_refl(R=0.95) ~240 mm → 159.4 mm (EXEC_23 FIX-05 source); (2) σ_END END-only row 47.6/52.1/49.6 → 54.86/53.36/50.49 ps (EXEC_23 FIX-01 source, m*-opt x=0 from results_macros.tex). README.md: add "MSI branch topology" paragraph (Option A decision). |
| FIX-25 | DONE | ba00e0e | EXEC_23_REPORT.md, CONFIGURATION_AUDIT.md | EXEC_23_REPORT.md FIX-03 row: note G1/G2 "combined frame" was split into G1+G2 by EXEC_24 FIX-11. CONFIGURATION_AUDIT.md BLUE weights §: clarify formula role — in w_E=(σ_T²−C_ET)/(σ_E²+σ_T²−2C_ET), σ_T² appears in numerator and denominator; σ_E² appears in denominator only. |

---

## §4 — Commit Log (pre-exec25-260902..HEAD)

```
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

---

## §6 — MUST-NOT-CHANGE Gate

Named constants verified present (as literals or via macros) in talk_v6.tex after all FIX commits.
FIX-23 replaced literal `0.95` with `\Rvikuiti` (expands to `0.95`; rendered value unchanged).
FIX-21 changed `13` → `10` (13 mm was not a named constant; 10 mm = bar height, physically correct).
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

## §8 — MSI Worktree Status

```
/mnt/d/SHiP/ej200                → exec25-deckfix @ ba00e0e (EXEC_25 build)
/mnt/d/SHiP/ej200_edge_scan      → main @ 5e90025 (frozen pre-EXEC_23; do not touch)
/mnt/d/SHiP/ej200_endonly        → feat/endonly-mylar @ fb3749d
/mnt/d/SHiP/ej200_event_display  → feat/ej204-event-display-tracks @ 47a9a4f
W1 /home/reriosto/SHiP/ej200     → exec25-deckfix @ ba00e0e (via git fetch + reset or manual)
```

---

## §9 — Open Items Carried Forward

All OPEN-01..06 items from EXEC_23_REPORT.md remain open. No new open items from EXEC_25 Phase 1.

**OPEN-07 (new):** SHiP_timing_first_principles.pptx (575 KB, 2026-08-19) is untracked on HOST.
Not a build artifact; manually created PowerPoint. Decision: do not add to .gitignore.
Future: user decides whether to track or keep untracked.

---

## §10 — Synchronization State

| Location | Branch | HEAD | Ahead of origin |
|----------|--------|------|----------------|
| HOST main | main | ba00e0e | 24 commits |
| MSI exec25-deckfix | exec25-deckfix | ba00e0e | — |
| GitHub origin/main | main | 5e90025 | — |

Push commands (NO push yet — user approval required):
```
# Push HOST main to GitHub:
git push origin main

# Print only — no auto-push:
# git --no-pager log --oneline origin/main..main | wc -l
```

No push executed. No merge. No rebase.
