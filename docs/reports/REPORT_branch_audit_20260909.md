# ej200 forensic branch audit — 2026-09-09
Scope: HOST `t0minidaq`, the discovered additional clone, MSI through the existing reverse tunnel, and live GitHub origin. Audit time: 2026-09-09T16:27:34+02:00. Report and evidence are outside both repositories.
## 1. Verdict
**`main` at `8349041140958226a0ac1cb3bb3e30aff2303435` is the branch containing the specified current deliverable, and matches live `origin/main`. No inspected branch satisfies all requested scientific checks.** The exact PDF passes criterion E, but main fails the strict sidecar completeness check and the requested `69.2 / TOP_SUM4_N1` marker check; design documentation still contains a suspect `701.3` value. A `polishedbackpainted` active reflector was not found. These are observed gaps, not grounds to silently substitute an older branch.
Numbered evidence against criteria A–E:
1. **A — available clones:** both local working trees are clean, HOST is on main, and MSI authentication failed. W1/W2 were not identified in the bounded filesystem search. [E008](#e008) [E006](#e006) [E005](#e005)
2. **B — topology:** HOST has **11 local branches**, while live origin has **18 branch names**. Only local main and origin/main contain `8349041`; origin/HEAD is a symbolic alias. Only `feat/bar-end-vikuiti` is another local branch fully contained in main. [E036](#e036) [E030](#e030) [E035](#e035) [E033](#e033)
3. **C — optical configuration:** main uses `CreateBarSkinReflector()` on the air–wrap border (DetectorConstruction.cc:313,391–396), with `dielectric_dielectric`, `polished`, and `REFLECTIVITY={0.98,0.98}` (Materials.cc:347–357). This is **not** `polishedbackpainted`. `dielectric_metal` also appears in SiPM and another reflector factory; a global text hit alone cannot identify the active reflector. [E211](#e211) [E212](#e212)
4. **C — EXEC_25:** main deck contains hybrid **514.9** at line 787 and **~1300** pooled photons at line 533. Its only literal `0.95` is the `\Rvikuiti` definition at line 33. No `TOP_SUM4_N1` or contextually relevant `69.2` was found in the inspected tracked trees. However, `presentations/v6/FINAL_NUMBERS.md` declares `N_TOP=20 unless noted` at line 6 and still gives **701.3** at lines 24 and 38 without an END-only exception on those rows. **SUSPICIOUS** configuration inconsistency. Explicit END-only 701.3 uses in the deck are retained as historical values, not automatically treated as errors. [E244](#e244) [E246](#e246) [E309](#e309)
5. **C — sidecars:** the main TeX references 14 unique figure stems. **0/14** have all three adjacent `.root`, `.csv`, `.meta.json` files. Two have CSV and `_meta.json` (a different suffix), but no adjacent ROOT. A separate `blue_wscan_x0.root` exists under analysis; it does not complete the figure-specific trios. [E080](#e080) [E244](#e244)
6. **D — rollback:** `pre-exec23-260902`, `pre-exec24-260902`, and `pre-exec25-260902` are contained in main. The first two are absent from the live origin tag listing. June `v-exec2*` tags are separate campaigns and many tag tips have no containing branch; keep them. [E092](#e092) [E031](#e031)
7. **E — exact artifact:** `presentations/v6/talk_v6.pdf` at main, origin/main, and 8349041 is **472931 bytes**, SHA-256 **`45784ae6a46ab2306004ade963b9e897611a21bd6f51877975b4347dc11a55de`**. Exact size and supplied hash-prefix match. No full expected SHA-256 was supplied, so only the supplied prefix can be compared with an expectation. [E238](#e238) [E239](#e239)
## 2. Clone inventory and divergence
| Clone | Path / endpoint | Branch and full HEAD | HEAD commit date | Working tree / worktrees / stash | Submodules / LFS |
|---|---|---|---|---|---|
| HOST | `/home/rrios/ej200` | `main`<br>`8349041140958226a0ac1cb3bb3e30aff2303435` (`8349041`) [E010](#e010) | 2026-09-03T00:37:35+02:00 | Clean; one worktree at this path; 1 stash(es). [E008](#e008) [E011](#e011) [E012](#e012) | No submodule status entries / no index gitlinks; no tracked HEAD LFS pointers or LFS attributes found. LFS executable unavailable; LFS object-store state UNVERIFIED. [E013](#e013) [E014](#e014) [E233](#e233) |
| Additional local clone (not identified as W1/W2) | `/home/rrios/ej200_end` | `feat/endonly-mylar`<br>`fb3749def29716dc84a33fcad53a21086bc96822` (`fb3749d`) [E021](#e021) | 2026-08-14T16:24:33+02:00 | Clean; one worktree at this path; 0 stash(es). [E019](#e019) [E022](#e022) [E023](#e023) | No submodule status entries / no index gitlinks; no tracked HEAD LFS pointers or LFS attributes found. LFS executable unavailable; LFS object-store state UNVERIFIED. [E024](#e024) [E025](#e025) [E236](#e236) |
| MSI | Requested `/mnt/d/SHiP/ej200`, `rrios@localhost:9022` | UNVERIFIED | UNVERIFIED | UNVERIFIED — SSH exit 255, permission denied | UNVERIFIED |
| W1 / W2 | Not found in searched roots | UNVERIFIED | UNVERIFIED | Not reachable / identity not established | UNVERIFIED |
| origin | `git@github.com:dowiyogo/ej200.git` | default main; `8349041140958226a0ac1cb3bb3e30aff2303435` | 2026-09-03T00:37:35+02:00 (same commit as HOST) | Working tree, worktrees, stashes: N/A to Git remote | Server LFS storage UNVERIFIED |
MSI was tried using the observed local username `rrios`; whether MSI expects a different user is UNVERIFIED. No new tunnel or alternate user login was attempted. Discovery searched `/home` to depth 4 and `/home /mnt /media /opt /srv` to depth 5. The broader find exited 1 with stderr suppressed; it does not prove absence beyond accessible searched paths. The original requested third physical clone is not substituted with an invented W1/W2 identity. [E003](#e003) [E005](#e005) [E006](#e006) [E231](#e231)
| Left tip | Right tip | Left-only / right-only commits | Evidence |
|---|---|---:|---|
| HOST HEAD | origin/main | 0 / 0 | [E228](#e228) |
| HOST HEAD | additional clone HEAD | 110 / 23 | [E229](#e229) |
| additional clone HEAD | origin/main | 23 / 110 | [E230](#e230) |
| MSI HEAD | HOST / origin/main / additional clone | UNVERIFIED / UNVERIFIED | Authentication failure |
**Origin reachability:** all inspected local branch commits on HOST are reachable from live origin branch tips; the additional clone’s HEAD is exactly origin/feat/endonly-mylar. Its second local branch, main at `84e902c522d8af54ace11e5bb38202bf70e8c94d`, also has no commits absent from origin remote-tracking branches after HOST fetch. Being unmerged into main does not mean a commit is absent from origin. [E235](#e235) [E261](#e261)
The stash reachability test was also repeated against **every live origin-advertised ref**, including tags: the same two commits remain absent. [E332](#e332)
**Preserve the HOST stash.** The stash and index-parent commits `927c006` and `74d0cd1` are not reachable from origin remote-tracking branches. The stash patch adds 57 lines to `docs/branch_diagnosis/DATA_AUDIT_2026-08-31.md`. No claim is made about private/unadvertised GitHub refs or inaccessible clones. [E227](#e227) [E262](#e262)
`feat/bar-end-vikuiti` is 13 commits ahead of its same-named origin branch, but all 13 are already reachable from origin/main; they are not origin-missing commits. [E258](#e258) [E260](#e260)
## 3. Branch table — 18 observed branch names
There are 11 HOST-local heads and seven names represented only by origin refs. Each row uses the local tip when present; otherwise the verified origin tip. Ahead / behind is **branch-only / main-only** (the reverse order of the raw `main...branch` command). “Exclusive” means reachable from branch but not main; it does not establish patch uniqueness or uniqueness against every other feature branch. File stats are the requested merge-base-to-branch diff. Full commits, file lists, markers, and evidence are in section 7.
| Branch / scope / tip | Ahead / behind main | Merged / contains 8349041 | Origin exists / same tip | Exclusive commits / files | Optical surface / reflectivity | EXEC_25: 514.9 / 1300 / 69.2 | Deck sidecars |
|---|---:|---|---|---|---|---|---|
| `diag/phase7-delta-2026-08-31`<br>HOST local; `a18e863` | 6 / 93 | No / No | Yes / Yes | 6 commits; 41 files changed, 13507 insertions(+); [details](#b1) | Air-gap + polished dielectric_dielectric reflector border; 0.95 active | MISS / MISS / MISS | 0/14 complete trios |
| `docs/branch-diagnosis-2026-08-31`<br>HOST local; `2a9645f` | 4 / 93 | No / No | Yes / Yes | 4 commits; 3 files changed, 1054 insertions(+); [details](#b2) | Air-gap + polished dielectric_dielectric reflector border; 0.95 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `exp/pair-scan-2026-06-11`<br>HOST local; `f19c093` | 30 / 115 | No / No | Yes / Yes | 30 commits; 444 files changed, 13222 insertions(+), 128 deletions(-); [details](#b3) | dielectric_metal polished bar skin; 0.98 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/bar-end-vikuiti`<br>HOST local; `50cf02e` | 0 / 74 | Yes / No | Yes / No: 1bfb827 | 0 commits; 0 files; [details](#b4) | Air-gap + polished dielectric_dielectric reflector border; 0.95 active | MISS / MISS / MISS | 0/14 complete trios |
| `feat/bar-vikuiti`<br>HOST local; `219fbe3` | 36 / 115 | No / No | Yes / Yes | 36 commits; 452 files changed, 15402 insertions(+), 135 deletions(-); [details](#b5) | dielectric_metal polished bar skin; 0.98 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej204-bar-tir-only`<br>HOST local; `09f8b18` | 32 / 115 | No / No | Yes / Yes | 32 commits; 446 files changed, 13647 insertions(+), 131 deletions(-); [details](#b6) | polished dielectric_dielectric TIR-only bar skin; 0.98 factory retained, not used by bar skin | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej204-event-display-tracks`<br>origin only; `47a9a4f` | 4 / 110 | No / No | Yes / Yes | 4 commits; 10 files changed, 425 insertions(+), 113 deletions(-); [details](#b7) | dielectric_metal polished bar skin; 0.98 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej228-cylinder`<br>origin only; `66b674c` | 3 / 103 | No / No | Yes / Yes | 3 commits; 18 files changed, 766 insertions(+), 575 deletions(-); [details](#b8) | Cylinder air-gap + dielectric_metal polished outer border; 0.98 active; 0.95 unused bar factory retained | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej228-tir-only`<br>origin only; `0006919` | 4 / 103 | No / No | Yes / Yes | 4 commits; 23 files changed, 2523 insertions(+), 575 deletions(-); [details](#b9) | Cylinder polished dielectric_dielectric TIR boundary; 0.98/0.95 factories retained, not used by cylinder reflector | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej230-bar-tir-only`<br>HOST local; `b281aea` | 34 / 115 | No / No | Yes / Yes | 34 commits; 448 files changed, 14226 insertions(+), 135 deletions(-); [details](#b10) | polished dielectric_dielectric TIR-only bar skin; 0.98 factory retained, not used by bar skin | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej230-endonly-mylar`<br>origin only; `04a8047` | 42 / 115 | No / No | Yes / Yes | 42 commits; 440 files changed, 30829 insertions(+), 838 deletions(-); [details](#b11) | dielectric_metal ground Mylar surface; 0.90 active; 0.98 unused bar factory retained | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/ej230-sslg4`<br>origin only; `5b93f4c` | 34 / 115 | No / No | Yes / Yes | 34 commits; 395 files changed, 27341 insertions(+), 192 deletions(-); [details](#b12) | dielectric_metal polished bar skin; 0.98 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/endonly-mylar`<br>origin only; `fb3749d` | 23 / 110 | No / No | Yes / Yes | 23 commits; 76 files changed, 14035 insertions(+), 153 deletions(-); [details](#b13) | dielectric_metal ground/polished configurable bar skin; 0.90 Mylar default; 0.98 fallback; runtime overrides UNVERIFIED | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feat/endtop-sslg4`<br>HOST local; `5576687` | 2 / 103 | No / No | Yes / Yes | 2 commits; 10 files changed, 451 insertions(+), 5 deletions(-); [details](#b14) | Air-gap + dielectric_metal polished outer reflector border; 0.95 active default | MISS / MISS / MISS | MISS: no talk_v6 source |
| `feature/sipm-electronics-response`<br>origin only; `bb4cf6a` | 15 / 163 | No / No | Yes / Yes | 15 commits; 60 files changed, 8040 insertions(+), 150 deletions(-); [details](#b15) | Passive Mylar volume; polished dielectric boundaries; No explicit 0.95/0.98 REFLECTIVITY assignment | MISS / MISS / MISS | MISS: no talk_v6 source |
| `main`<br>HOST local; `8349041` | 0 / 0 | Yes / Yes | Yes / Yes | 0 commits; 0 files; [details](#b16) | Air-gap + polished dielectric_dielectric reflector border; 0.98 active; historical deck data 0.95 | HIT / HIT / MISS | 0/14 complete trios |
| `wip/host-stash-endtop-junio`<br>HOST local; `9710f34` | 1 / 136 | No / No | Yes / Yes | 1 commits; 2 files changed, 16 insertions(+), 13 deletions(-); [details](#b17) | dielectric_metal polished explicit sibling-panel border; 0.98 active | MISS / MISS / MISS | MISS: no talk_v6 source |
| `wip/host-uncommitted-2026-08-31`<br>HOST local; `d2c8d4c` | 2 / 93 | No / No | Yes / Yes | 2 commits; 4 files changed, 192 insertions(+), 15 deletions(-); [details](#b18) | Air-gap + polished dielectric_dielectric reflector border; 0.95 active | MISS / MISS / MISS | MISS: no talk_v6 source |
All rows have **MISS for active `polishedbackpainted`**. All except the passive-wrap electronics branch have `dielectric_metal` factory/source hits, but their active reflector differs as described. TIR-only and cylinder branches are not automatically called invalid because another unused factory is metallic. For branches without a talk_v6 tree, EXEC_25 marker MISS means the requested deliverable context is absent, not that a coincidental raw-data number was impossible.
The differing origin/feat/bar-end-vikuiti tip is `1bfb82743270dcc3ea9b17be03a6af5704cc6850` and retains `analysis/presentation_v6/talk_v6.tex`; its local tip has the relocated `presentations/v6/talk_v6.tex`. Both retain R=0.95 code and lack the final EXEC_25 corrections. [E256](#e256) [E257](#e257) [E328](#e328)
## 4. Taxonomy
These groups overlap: scientific obsolescence is not permission to discard unmerged work.
**(a) Merged, eligible for local branch-label cleanup:** `feat/bar-end-vikuiti` only, at `50cf02e...`. It has 0 commits exclusive to main, and HOST worktree inventory shows main checked out. Main itself is retained. Existing same-named origin is older, so create the proposed exact-tip backup tag before any manual deletion. Remote deletion is deferred while MSI/W1/W2 remain unverified.
**(b) Exclusive, unmerged work — retain:** `diag/phase7-delta-2026-08-31` (6 commits), `docs/branch-diagnosis-2026-08-31` (4 commits), `exp/pair-scan-2026-06-11` (30 commits), `feat/bar-vikuiti` (36 commits), `feat/ej204-bar-tir-only` (32 commits), `feat/ej204-event-display-tracks` (4 commits), `feat/ej228-cylinder` (3 commits), `feat/ej228-tir-only` (4 commits), `feat/ej230-bar-tir-only` (34 commits), `feat/ej230-endonly-mylar` (42 commits), `feat/ej230-sslg4` (34 commits), `feat/endonly-mylar` (23 commits), `feat/endtop-sslg4` (2 commits), `feature/sipm-electronics-response` (15 commits), `wip/host-stash-endtop-junio` (1 commits), `wip/host-uncommitted-2026-08-31` (2 commits). The complete commit lists are printed for each branch in section 7; none is shortened to the first 20. Shared commits may appear in multiple lists. No patch-equivalence test was used to justify deletion.
**(c) Obsolete for the target Phase-7 bar-reflector baseline:** `exp/pair-scan-2026-06-11`, `feat/bar-vikuiti`, `feat/ej204-event-display-tracks`, `feat/ej230-endonly-mylar`, `feat/ej230-sslg4`, `feat/endonly-mylar`, `feat/endtop-sslg4`, `wip/host-stash-endtop-junio`. Their inspected geometry uses a metallic bar skin, older sibling-panel boundary, or the different metallic air-gap reflector implementation. See source excerpts and call sites in section 7. “Obsolete” here means different from the target main configuration; this audit did not rerun physics or invalidate every historical study.
**(d) Indeterminate or incomplete scientifically:** TIR-only EJ-204/EJ-230, both EJ-228 cylinder studies, and the electronics-response branch have different purposes; equivalence to the requested scientific deliverable is UNVERIFIED. Diagnosis/WIP branches with historical dielectric air-gap code lack the final deliverable. Main remains the current specified artifact but is incomplete against the full audit checklist. MSI/W1/W2 are UNVERIFIED.
## 5. Tags and rollback mapping
Dates below are tag creator dates (commit dates for lightweight tags), not dates inferred from tag-name strings. Campaign labels are explicitly taken from names and stored tag/commit messages. June optical EXEC_25 and September presentation EXEC_25 are different contexts. A containing branch is a reachability result, not evidence of the branch on which the tag was originally created.
| Date | Tag | Tag object / peeled commit | Campaign mapping | Containing branches | On origin |
|---|---|---|---|---|---|
| 2026-09-02 | `pre-exec25-260902` | `1c5a4c5` / `1c5a4c5142efa555f49ba9db1c8e9a2eeb411035` | Before EXEC_25 | G1 [E107](#e107) | Yes |
| 2026-09-02 | `pre-exec24-260902` | `09dc9e2` / `09dc9e21742340cd277e83f8c27f4c79cb04c831` | Before EXEC_24 | G1 [E106](#e106) | No |
| 2026-09-01 | `pre-exec23-260902` | `5e90025` / `5e90025bc4038835f1840f6f11dca3ca27778d4f` | Before EXEC_23 | G1 [E105](#e105) | No |
| 2026-06-21 | `v-exec27-edge-estimator-coupling` | `bd24951` / `bd2495183420a8097c9742069a0ef556cfc831a7` | EXEC_27 | G2 [E114](#e114) | Yes |
| 2026-06-21 | `v-exec26-scint-air-surface-realism` | `3148da1` / `3148da19a3d42b300f71d87d32f55d2064ba1b93` | EXEC_26 | G2 [E113](#e113) | Yes |
| 2026-06-21 | `v-exec25-optical-realism-bracket` | `7c55326` / `7c553263bbe40302b4079cf40bcbac74f72e5544` | EXEC_25 | G2 [E112](#e112) | Yes |
| 2026-06-21 | `v-exec24-pe-budget-audit` | `e20b47d` / `e20b47d492e11f099aefb1074b1d8bd26cf8a88a` | EXEC_24 | G2 [E111](#e111) | Yes |
| 2026-06-21 | `v-exec23-explicit-airgap` | `bd78211` / `bd7821170c1c6997a733d0f279f3eac85821c10c` | EXEC_23 | G2 [E110](#e110) | Yes |
| 2026-06-21 | `v-exec22-endtop-optfix` | `76b582c` / `76b582c5cf982706d2b72d3b82a3f879696957f8` | EXEC_22 | G2 [E109](#e109) | Yes |
| 2026-06-21 | `v-exec21-optfix` | `9b8361f` / `9b8361f189712fd8dba332a47e2ed685f194f1c1` | EXEC_21 label; stored commit message says exec22 follow-up | G2 [E108](#e108) | Yes |
| 2026-06-18 | `diag-photon-budget-v1` | `0ca6855` / `0ca685506faa4d309a7dc9ab13b18c85e127b807` | Photon-budget diagnostic; EXEC_N UNVERIFIED | G2 [E102](#e102) | Yes |
| 2026-06-18 | `physics-baseline-v1` | `26f4767` / `26f4767a3643b5ab4faffccd9cfe23ac1b3042a9` | Physics baseline; EXEC_N UNVERIFIED | G2 [E104](#e104) | Yes |
| 2026-06-18 | `exec07-09-analysis` | `c9e61af` / `c9e61afd8a85182ac8ec23abf8f0e020b1972930` | EXEC_07–09 (tag name) | G2 [E103](#e103) | Yes |
| 2026-06-11 | `checkpoint/pre-pairscan-2026-06-11` | `5783a0d` / `5783a0d2914b01fac6f9f3f48cdb02879384aa86` | Before pair scan; preserves EXEC_12b | G3 [E100](#e100) | Yes |
| 2026-06-11 | `checkpoint/pre-exec12b-2026-06-11` | `f767be8` / `f767be8709a424b042361014fd016673f00d6b69` | Before EXEC_12b | G3 [E099](#e099) | Yes |
| 2026-06-11 | `checkpoint/pre-exec12-beamer-2026-06-11` | `6e8d3d6` / `6e8d3d6ed22073185ce57ca14ed34be1c56af463` | Before EXEC_12 | G3 [E098](#e098) | Yes |
| 2026-06-11 | `checkpoint/pre-exec11-2026-06-11` | `7596697` / `759669750e8d7396c6d1e4afc2e1922fc119e682` | Before EXEC_11 | G3 [E096](#e096) | Yes |
| 2026-06-11 | `checkpoint/pre-exec11b-2026-06-11` | `7596697` / `759669750e8d7396c6d1e4afc2e1922fc119e682` | Before EXEC_11b | G3 [E097](#e097) | Yes |
| 2026-06-10 | `checkpoint/pre-endtop-sslg4-2026-06-10` | `2a9d57b` / `2a9d57b454cfb2e659c8551b6c8d12c1b2a34b2b` | Before EndTop fork; target commit says CODEX_EXEC_09 | G4 [E095](#e095) | Yes |
| 2026-06-08 | `checkpoint/pre-physics-baseline-2026-05-08` | `e9245bb` / `4e469599670608f7a5b74e56ea0289187b176356` | Before physics baseline; EXEC_N UNVERIFIED | G4 [E101](#e101) | Yes |
- **G1**: `main`, `origin/main`
- **G2**: **No local or origin-tracking branch contains this tag tip. Preserve the tag itself.**
- **G3**: `diag/phase7-delta-2026-08-31`, `docs/branch-diagnosis-2026-08-31`, `exp/pair-scan-2026-06-11`, `feat/bar-end-vikuiti`, `feat/bar-vikuiti`, `feat/ej204-bar-tir-only`, `feat/ej230-bar-tir-only`, `feat/endtop-sslg4`, `main`, `wip/host-uncommitted-2026-08-31`, `origin/diag/phase7-delta-2026-08-31`, `origin/docs/branch-diagnosis-2026-08-31`, `origin/exp/pair-scan-2026-06-11`, `origin/feat/bar-end-vikuiti`, `origin/feat/bar-vikuiti`, `origin/feat/ej204-bar-tir-only`, `origin/feat/ej204-event-display-tracks`, `origin/feat/ej228-cylinder`, `origin/feat/ej228-tir-only`, `origin/feat/ej230-bar-tir-only`, `origin/feat/ej230-endonly-mylar`, `origin/feat/ej230-sslg4`, `origin/feat/endonly-mylar`, `origin/feat/endtop-sslg4`, `origin/main`, `origin/wip/host-uncommitted-2026-08-31`
- **G4**: `diag/phase7-delta-2026-08-31`, `docs/branch-diagnosis-2026-08-31`, `exp/pair-scan-2026-06-11`, `feat/bar-end-vikuiti`, `feat/bar-vikuiti`, `feat/ej204-bar-tir-only`, `feat/ej230-bar-tir-only`, `feat/endtop-sslg4`, `main`, `wip/host-stash-endtop-junio`, `wip/host-uncommitted-2026-08-31`, `origin/diag/phase7-delta-2026-08-31`, `origin/docs/branch-diagnosis-2026-08-31`, `origin/exp/pair-scan-2026-06-11`, `origin/feat/bar-end-vikuiti`, `origin/feat/bar-vikuiti`, `origin/feat/ej204-bar-tir-only`, `origin/feat/ej204-event-display-tracks`, `origin/feat/ej228-cylinder`, `origin/feat/ej228-tir-only`, `origin/feat/ej230-bar-tir-only`, `origin/feat/ej230-endonly-mylar`, `origin/feat/ej230-sslg4`, `origin/feat/endonly-mylar`, `origin/feat/endtop-sslg4`, `origin/main`, `origin/wip/host-stash-endtop-junio`, `origin/wip/host-uncommitted-2026-08-31`
Tag dates and subjects: [E092](#e092). Tag messages: [E093](#e093). Live remote existence: [E031](#e031). Do not treat `pre-exec25-260902` as a backup of the final 8349041 deliverable: it precedes those corrections.
## 6. Proposed cleanup — NOT EXECUTED
Only one local branch-label deletion is proposed. No source edits, remote deletions, tag deletions, stash drops, or clone removals are proposed. All commands in the following block are for René to review and execute manually. Recheck live main, origin/main, worktrees and the exact branch tip before use; this is a dated snapshot.
```bash
# Preserve the final deliverable; proposed new tag, NOT created by this audit.
git -C /home/rrios/ej200 tag audit/20260909/current-main 8349041140958226a0ac1cb3bb3e30aff2303435

# Preserve the exact merged local branch tip; origin of the same name is 13 commits older.
git -C /home/rrios/ej200 tag audit/20260909/feat-bar-end-vikuiti 50cf02e34cd7dddfb847739f6249afd1a600478f

# Preserve the stash commit and its parents; keep the stash itself.
git -C /home/rrios/ej200 tag audit/20260909/host-data-audit-stash 927c0069a238a31a68db796355326242fe51abd4

# DESTRUCTIVE, MANUAL ONLY: remove this local branch label after verifying the backup tag.
# Justification: exact tip is contained in main and origin/main; zero main-exclusive commits.
# Backup: audit/20260909/feat-bar-end-vikuiti. No worktree here has this branch checked out.
git -C /home/rrios/ej200 branch -D feat/bar-end-vikuiti
```
The force spelling is explicit because the branch tracks an older origin counterpart even though its tip is contained in main. This is not permission to use the same command for any unmerged branch. Backup tags above are local proposals, not published remote backups. Tag creation and branch deletion were **not executed**.
## 7. Per-branch evidence: science, complete exclusive commits, and file stats
<a id="b1"></a>

### 1. diag/phase7-delta-2026-08-31
Inspected ref `diag/phase7-delta-2026-08-31`, full SHA `a18e863bab76cb54ca96a86c9e871ff3501d434c`. Classification: **Historical Phase-7 geometry; pre-EXEC_25 deck**. [E032](#e032)
**B:** ahead 6, behind 93; merged NO; contains 8349041 NO. [E037](#e037) [E223](#e223) [E035](#e035)
Upstream: `origin/diag/phase7-delta-2026-08-31` (no ahead/behind annotation). Origin counterpart: `a18e863bab76cb54ca96a86c9e871ff3501d434c`. [E030](#e030)
Tip date / author / subject: 2026-09-01 00:01:29 +0200 / rrios / feat(analysis): track talk_v6 figure generation scripts and configs. [E032](#e032)
**C1/C2:** active configuration: Air-gap + polished dielectric_dielectric reflector border. Reflectivity: 0.95 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E116](#e116) [E117](#e117) [E290](#e290)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.95, 0.95};
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 HIT**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **90**. Full paths, line numbers and contents are preserved in [E118](#e118) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E117](#e117)
```text
311:     {
312:         auto* scintAirSurface     = Materials::CreateBarSurface();
313:         auto* airReflectorSurface = Materials::CreateBarSkinReflector();
314:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E116](#e116) [E290](#e290)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
332:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
354:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E291](#e291). No matches outside the bundled external libraries.
Deck `analysis/presentation_v6/talk_v6.tex`: **514.9 MISS**, **1300 MISS**. [E289](#e289)

```text
33: \newcommand{\Rvikuiti}{0.95}
35: \newcommand{\NpeG}{701.3}
82:     \textbf{Reflector:} Vikuiti ESR on all non-SiPM surfaces ($R = \Rvikuiti$).\\[5pt]
415:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
756:       $N_\text{pe}$/end at $x=0$ & $701.3$ \\
916:     $N_\text{pe}$ & $\sim\!0.37$ (BUG) & — & 701.3 \\
963:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
995:       \textbf{G4/end} & \textbf{701.3} & +refl.-rec. \\
```

Exact loose `0.95` occurrences in the deck: **16**. **OPEN-07-style macro gap:** historical literals remain; the requested final macro cleanup is absent. This label describes the requested audit check, not a verified issue identifier.
```text
33: \newcommand{\Rvikuiti}{0.95}
200:       \item Reflected ($R=0.95$ per bounce) back into bar
251:     At $L/2=700\mm$: survival $= 0.95^{85.6} \approx 1.2\%$ (code $R=0.95$). \\[4pt]
363:       \item Air$\to$Mylar: \texttt{CreateBarSkinReflector()} as border surface — \texttt{dielectric\_dielectric}, $R=0.95$ constant
741:       \item Reflector: \textbf{Vikuiti ESR} ($R=0.95$ constant in code)
790:       \item Vikuiti (R=0.95) recovers non-TIR photons (+23\%)
845:       \item $R = 0.95$ constant (wavelength-independent)
857:     \textbf{A2 — R=0.95 vs R=0.98}\\[4pt]
858:     Code value (Materials.cc): $R = 0.95$\\
860:     The code comment says ``R=0.98 Vikuiti ESR'' but the actual constant set is 0.95.
861:     All simulation results in this analysis use \alert{$R=0.95$}.\\[6pt]
864:     With R=0.95: $\Lambda_\text{refl}^H = 240\mm$ (shorter — more loss).\\[4pt]
887:       \item Vikuiti reflectivity at air–Mylar interface ($R=0.95$)
914:     Reflector & \texttt{dielectric\_metal skin} & skin surface, R=0.95 & border surface, R=0.95 \\
939:     \textbf{Fix (exec21-optfix):} Changed to \texttt{dielectric\_dielectric} border surface with $R=0.95$, preserving natural TIR at bar–air interface.
1192:       Reflector & Vikuiti ESR ($R=0.95$) \\
```

701.3 contextual caution: **SUSPICIOUS**: the old deck retains 701.3 before the hybrid correction; inspect the design frame below.

```text
32: \newcommand{\LambdaH}{405\mm}
33: \newcommand{\Rvikuiti}{0.95}
34: \newcommand{\NpeNapkin}{570}
35: \newcommand{\NpeG}{701.3}
36: \newcommand{\sigmaENDval}{53.68\ps}
37: \newcommand{\sigmaTOPval}{15.20\ps}
38: \newcommand{\sigmaBLUEval}{15.21\ps}
412:       \midrule
413:       TIR-guided / end & 570 & — \\
414:       Reflector-recovered & 0 & — \\
415:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
416:       \bottomrule
417:     \end{tabular}\\[6pt]
418:     $\Delta N_\text{pe}/N_\text{nap} = +23\%$.\\[4pt]
732: \section{Design Summary}
733: % ====================================================================
734: 
735: \begin{frame}{Design Decision: EJ-230 + 20 TOP + 16 END}
736:   \begin{columns}[T]
737:     \column{0.52\linewidth}
738:     \textbf{Selected configuration:}
753:       $\sigma_t$ at $x=0$ & $\mathbf{15.2\ps}$ \\
754:       Mean $\sigma_t$ across bar & $19.5\ps$ \\
755:       $\sigma_x$ (END $\Delta t$) & $\mathbf{7.9\mm}$ \\
756:       $N_\text{pe}$/end at $x=0$ & $701.3$ \\
757:       \bottomrule
758:     \end{tabular}
759:     \end{center}
913:     Air gap & No & No & Yes (0.10 mm) \\
914:     Reflector & \texttt{dielectric\_metal skin} & skin surface, R=0.95 & border surface, R=0.95 \\
915:     TIR & Eliminated by skin & Natural (no bar surface) & Natural Fresnel \\
916:     $N_\text{pe}$ & $\sim\!0.37$ (BUG) & — & 701.3 \\
917:     $\sigma_\text{END}$ & N/A (insufficient pe) & 50.5 ps & 53.7 ps \\
918:     TOP & No & No & Yes ($N=4,8,14,20$) \\
919:     Status & \textcolor{red}{\textbf{INVALID}} & superseded & \textcolor{green!60!black}{\textbf{current}} \\
960:     Yield & 10000 ph/MeV & 10400 ph/MeV & 9700 ph/MeV \\
961:     $n$ & 1.58 & 1.58 & 1.58 \\
962:     \midrule
963:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
964:     $N_\text{pe}$/end (napkin) & — & — & 570 \\
965:     $\sigma_t$ END-only (G4) & 47.6 ps & 52.1 ps & 49.6 ps \\
966:     Napkin bulk survival & 0.832 & 0.646 & \textbf{0.558} \\
992:       PDE & 0.400 & assumed \\
993:       \textbf{Napkin/end} & \textbf{570} & TIR-only \\
994:       \midrule
995:       \textbf{G4/end} & \textbf{701.3} & +refl.-rec. \\
996:       Surplus & $+23\%$ & Pop.\ II \\
997:       \bottomrule
998:     \end{tabular}}
```

**C4:** **0/14 complete exact-suffix figure trios.** [E040](#e040)
| Figure stem (relative to deck figs/) | Figure PDF | .root | .csv | .meta.json | _meta.json alternative | Matching sidecars elsewhere in tree |
|---|---|---|---|---|---|---|
| `fig1_kscan` | MISS | MISS | MISS | MISS | MISS | none |
| `figM1_mat_sigma_end` | MISS | MISS | MISS | MISS | MISS | none |
| `figM2_mat_npe` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_bulk_survival` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_end_mscan` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_npe_x` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_ntop_scan_full` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_sigma_t_x` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_survival_RN` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_pareto` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_sigma_vs_x` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_top_position_loo` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_veff_fit` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_veff_residual` | MISS | MISS | MISS | MISS | MISS | none |
Figure list is derived from literal `\anafig{...}` and `\includegraphics{figs/...}` calls in the tracked TeX. No dynamic figure paths were assumed; TeX build execution and contents inside figure PDFs were not tested.
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E038](#e038)
```text
a18e863 feat(analysis): track talk_v6 figure generation scripts and configs
6a29b09 docs(diagnosis): Phase 4 vs Phase 7 delta diagnostic — regeneration scope
2a9645f docs(diagnosis): data audit — dataset-to-optical-phase mapping for t0minidaq runs
3af336c docs(diagnosis): surface model impact analysis for feat/endtop-sslg4
e67d57f docs(diagnosis): V1/V2/V3 corrections + §11 Recomendación sobre main
822f06f docs(diagnosis): branch content diagnosis 2026-08-31
```

**Files touched on the branch since its merge base with main:** [E039](#e039)
```text
 analysis/optim/make_summary_plots.py               |  210 ++++
 analysis/optim/phase_ab.py                         |  202 ++++
 analysis/optim/phase_cd.py                         |  197 ++++
 analysis/optim/phase_e_pareto.py                   |  116 ++
 analysis/optim/phase_sparse_top.py                 |  203 ++++
 analysis/optim/root_best_est/README.md             |  105 ++
 analysis/optim/root_best_est/REVISION_NOTES.md     |  177 +++
 analysis/optim/root_best_est/best_est_analysis.py  |  721 ++++++++++++
 analysis/optim/root_best_est/new_analysis_plots.C  |  355 ++++++
 analysis/optim/root_best_est/verification.json     |  250 ++++
 analysis/presentation_v4/scripts/analysis_v4.py    | 1118 ++++++++++++++++++
 analysis/presentation_v4/scripts/fig_materials.C   |   65 ++
 analysis/presentation_v4/scripts/fig_materials.py  |  171 +++
 .../presentation_v4/tables/summary_numbers.json    |   26 +
 analysis/presentation_v4/talk_v4.tex               | 1067 +++++++++++++++++
 analysis/presentation_v5/BUILD.md                  |   66 ++
 analysis/presentation_v5/FINAL_NUMBERS.md          |  155 +++
 analysis/presentation_v5/REVISION_NOTES.md         |  101 ++
 analysis/presentation_v5/scripts/analysis_v5.py    |  780 +++++++++++++
 analysis/presentation_v5/talk_v5.tex               |  758 ++++++++++++
 analysis/presentation_v6/BUILD.md                  |   89 ++
 analysis/presentation_v6/CONFIGURATION_AUDIT.md    |  112 ++
 analysis/presentation_v6/CONTENT_AUDIT.md          |  123 ++
 analysis/presentation_v6/FINAL_NUMBERS.md          |  143 +++
 analysis/presentation_v6/REVISION_NOTES.md         |  102 ++
 analysis/presentation_v6/talk_v6.tex               | 1224 ++++++++++++++++++++
 docs/branch_diagnosis/DATA_AUDIT_2026-08-31.md     |  298 +++++
 docs/branch_diagnosis/DIAGNOSIS_2026-08-31.md      |  543 +++++++++
 .../IMPACT_surface_model_2026-08-31.md             |  213 ++++
 docs/branch_diagnosis/PHASE7_DELTA_2026-08-31.md   |  236 ++++
 presentations/best_est_2026-08-17/talk.tex         |  659 +++++++++++
 presentations/optim_2026-08-16/talk.tex            |  293 +++++
 presentations/optim_2026-08-17/talk.tex            |  301 +++++
 talks/napkin_first_principles/fig_gen.C            |  247 ++++
 talks/napkin_first_principles/fig_gen.py           |  194 ++++
 talks/napkin_first_principles/napkin.py            |  342 ++++++
 talks/napkin_first_principles/tex/talk.tex         |  556 +++++++++
 talks/napkin_first_principles/tex/talk_v3.tex      |  840 ++++++++++++++
 .../napkin_first_principles/values/convention.json |   30 +
 .../napkin_first_principles/values/materials.yaml  |   37 +
 .../values/napkin_macros.tex                       |   82 ++
 41 files changed, 13507 insertions(+)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E119](#e119)
<a id="b2"></a>

### 2. docs/branch-diagnosis-2026-08-31
Inspected ref `docs/branch-diagnosis-2026-08-31`, full SHA `2a9645fde563abfd20cfdde810124038e312f85e`. Classification: **Historical geometry; diagnosis work**. [E032](#e032)
**B:** ahead 4, behind 93; merged NO; contains 8349041 NO. [E042](#e042) [E223](#e223) [E035](#e035)
Upstream: none configured. Origin counterpart: `2a9645fde563abfd20cfdde810124038e312f85e`. [E030](#e030)
Tip date / author / subject: 2026-08-31 22:41:13 +0200 / rrios / docs(diagnosis): data audit — dataset-to-optical-phase mapping for t0minidaq runs. [E032](#e032)
**C1/C2:** active configuration: Air-gap + polished dielectric_dielectric reflector border. Reflectivity: 0.95 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E120](#e120) [E121](#e121) [E292](#e292)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.95, 0.95};
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 HIT**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **37**. Full paths, line numbers and contents are preserved in [E122](#e122) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E121](#e121)
```text
311:     {
312:         auto* scintAirSurface     = Materials::CreateBarSurface();
313:         auto* airReflectorSurface = Materials::CreateBarSkinReflector();
314:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E120](#e120) [E292](#e292)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
332:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
354:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E293](#e293). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E045](#e045)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E045](#e045)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E043](#e043)
```text
2a9645f docs(diagnosis): data audit — dataset-to-optical-phase mapping for t0minidaq runs
3af336c docs(diagnosis): surface model impact analysis for feat/endtop-sslg4
e67d57f docs(diagnosis): V1/V2/V3 corrections + §11 Recomendación sobre main
822f06f docs(diagnosis): branch content diagnosis 2026-08-31
```

**Files touched on the branch since its merge base with main:** [E044](#e044)
```text
 docs/branch_diagnosis/DATA_AUDIT_2026-08-31.md     | 298 +++++++++++
 docs/branch_diagnosis/DIAGNOSIS_2026-08-31.md      | 543 +++++++++++++++++++++
 .../IMPACT_surface_model_2026-08-31.md             | 213 ++++++++
 3 files changed, 1054 insertions(+)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E123](#e123)
<a id="b3"></a>

### 3. exp/pair-scan-2026-06-11
Inspected ref `exp/pair-scan-2026-06-11`, full SHA `f19c0933133c6ba25e2c4d968d03361fd804fd90`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 30, behind 115; merged NO; contains 8349041 NO. [E047](#e047) [E223](#e223) [E035](#e035)
Upstream: `origin/exp/pair-scan-2026-06-11` (no ahead/behind annotation). Origin counterpart: `f19c0933133c6ba25e2c4d968d03361fd804fd90`. [E030](#e030)
Tip date / author / subject: 2026-08-11 22:00:44 +0200 / dowiyogo / feat(pairscan): CMakeLists refactor, pair scan macros and run script. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal polished bar skin. Reflectivity: 0.98 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E124](#e124) [E125](#e125) [E294](#e294)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **6**. Full paths, line numbers and contents are preserved in [E126](#e126) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E125](#e125)
```text
254:     // Apply branch-specific reflector properties directly to the bar.
255:     auto* reflector = Materials::CreateBarSkinReflector();
256:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
257:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E124](#e124) [E294](#e294)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E295](#e295). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E050](#e050)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E050](#e050)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E048](#e048)
```text
f19c093 feat(pairscan): CMakeLists refactor, pair scan macros and run script
14ae395 fix(tests): replace escape-fraction guard with sipm-entry guard
c7a627e fix(optics): replace reflector volumes with bar skin surface
021d5f4 EXEC_12TB: add rebuilt self-contained beamer and speaker notes
1ac830a EXEC_12TB: add beamer manifest and Makefile exec12tb targets
d764892 EXEC_12TB: add generated tables and TeX number macros
9cfe042 EXEC_12TB: add figure style module and regenerated deck figures
66ac24c EXEC_12T: add complete Beamer presentation and reproducible build
b0aaac1 EXEC_12T: add integrated timing-position technical report
a2563d9 EXEC_12T: add window-dip timing and collection study
aaabe89 EXEC_12T: add threshold sweep and covariance timing decomposition
0a41034 EXEC_12T: add 4PE versus 20PE temporal and position analysis
f38c9e1 EXEC_12T: add cached fourth-to-thirtieth hit order statistics
003d106 EXEC_12: add staged XY scan proposal
6e01606 EXEC_12: add EXEC11 and EndTop technical reports
30166cd EXEC_12: add y-zero observability feasibility study
a9571fd EXEC_12: add covariance-aware EndTop estimator combination
bd21563 EXEC_12: add leave-one-position-out global X reconstruction
8332733 EXEC_12: add leave-one-position-out global X reconstruction
17d9e7d EXEC_12: add EndTop event observable builder and geometry tests
4dba735 EXEC_11: add reproducible report tables and README
a103ed6 EXEC_11: add temporal ratio and covariance-aware position reconstruction
1e95ac6 EXEC_11: add detailed two-position timing analysis
61bca2e EXEC_11: regenerate pair-scan fit QA and summary v2
bf0d2d5 EXEC_11: add per-event pair observable builder and tests
f431c01 EXEC_08/Step7: CTest guardrails for pair-scan geometry and macros
20bed82 EXEC_08/Step6: batch runner scripts/run_pair_scan.sh
549f3a6 EXEC_08/Step5: ROOT analysis macro analyze_pair_scan.C
30f3da3 EXEC_08/Step3: generate 41 pair-scan macros in macros/pairscan/
6bb73d5 EXEC_08/Step2: geometry analysis + pair selection
```

**Files touched on the branch since its merge base with main:** [E049](#e049)
```text
 .gitignore                                         |    6 +
 CMakeLists.txt                                     |   76 +-
 Makefile                                           |   40 +
 analysis/analyze_pair_scan.C                       |  496 +++++++
 analysis/exec11_pair_analysis.py                   |  934 +++++++++++++
 analysis/exec12_endtop_position.py                 |  161 +++
 analysis/exec12_make_reports.py                    |   99 ++
 analysis/exec12t_make_products.py                  |  151 +++
 analysis/exec12t_timing_threshold_analysis.py      |  136 ++
 analysis/exec12tb_figstyle.py                      |   56 +
 analysis/exec12tb_figures.py                       |  823 ++++++++++++
 analysis/exec12tb_tables.py                        |  487 +++++++
 include/DetectorConstruction.hh                    |    8 +-
 macros/pairscan/pairscan_x-422.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-423.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-424.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-425.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-426.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-427.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-428.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-429.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-430.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-431.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-432.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-433.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-434.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-435.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-436.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-437.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-438.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-439.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-440.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-441.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-442.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-443.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-444.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-445.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-446.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-447.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-448.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-449.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-450.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-451.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-452.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-453.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-454.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-455.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-456.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-457.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-458.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-459.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-460.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-461.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-462.0mm.mac             |   20 +
 pair_scan_config.json                              |   71 +
 report/exec13_xy_scan_plan.md                      |   30 +
 results/exec11_20260612_182454/README.md           |  124 ++
 .../analysis/calibration_ratio.csv                 |    3 +
 .../analysis/calibration_temporal.csv              |    2 +
 .../analysis/calibration_v1_v2.csv                 |    3 +
 .../analysis/data_inventory.csv                    |   42 +
 results/exec11_20260612_182454/analysis/fit_qa.csv |   42 +
 .../analysis/fit_qa_v1_v2_focus.csv                |    4 +
 .../exec11_20260612_182454/analysis/metadata.json  |   13 +
 .../analysis/pairscan_summary_v2.csv               |   42 +
 .../analysis/reconstruction_summary.csv            |   11 +
 .../analysis/reference_comparison.csv              |    3 +
 .../analysis/reference_positions.csv               |    3 +
 .../derived/pair_events_x-422.0mm.npz              |  Bin 0 -> 123534 bytes
 .../derived/pair_events_x-423.0mm.npz              |  Bin 0 -> 122740 bytes
 .../derived/pair_events_x-424.0mm.npz              |  Bin 0 -> 123199 bytes
 .../derived/pair_events_x-425.0mm.npz              |  Bin 0 -> 123147 bytes
 .../derived/pair_events_x-426.0mm.npz              |  Bin 0 -> 123368 bytes
 .../derived/pair_events_x-427.0mm.npz              |  Bin 0 -> 122917 bytes
 .../derived/pair_events_x-428.0mm.npz              |  Bin 0 -> 122839 bytes
 .../derived/pair_events_x-429.0mm.npz              |  Bin 0 -> 122318 bytes
 .../derived/pair_events_x-430.0mm.npz              |  Bin 0 -> 121727 bytes
 .../derived/pair_events_x-431.0mm.npz              |  Bin 0 -> 120954 bytes
 .../derived/pair_events_x-432.0mm.npz              |  Bin 0 -> 120412 bytes
 .../derived/pair_events_x-433.0mm.npz              |  Bin 0 -> 121077 bytes
 .../derived/pair_events_x-434.0mm.npz              |  Bin 0 -> 122386 bytes
 .../derived/pair_events_x-435.0mm.npz              |  Bin 0 -> 122900 bytes
 .../derived/pair_events_x-436.0mm.npz              |  Bin 0 -> 122357 bytes
 .../derived/pair_events_x-437.0mm.npz              |  Bin 0 -> 122009 bytes
 .../derived/pair_events_x-438.0mm.npz              |  Bin 0 -> 122595 bytes
 .../derived/pair_events_x-439.0mm.npz              |  Bin 0 -> 122198 bytes
 .../derived/pair_events_x-440.0mm.npz              |  Bin 0 -> 122232 bytes
 .../derived/pair_events_x-441.0mm.npz              |  Bin 0 -> 122736 bytes
 .../derived/pair_events_x-442.0mm.npz              |  Bin 0 -> 122672 bytes
 .../derived/pair_events_x-443.0mm.npz              |  Bin 0 -> 122553 bytes
 .../derived/pair_events_x-444.0mm.npz              |  Bin 0 -> 122542 bytes
 .../derived/pair_events_x-445.0mm.npz              |  Bin 0 -> 122293 bytes
 .../derived/pair_events_x-446.0mm.npz              |  Bin 0 -> 122306 bytes
 .../derived/pair_events_x-447.0mm.npz              |  Bin 0 -> 121906 bytes
 .../derived/pair_events_x-448.0mm.npz              |  Bin 0 -> 122022 bytes
 .../derived/pair_events_x-449.0mm.npz              |  Bin 0 -> 122871 bytes
 .../derived/pair_events_x-450.0mm.npz              |  Bin 0 -> 122284 bytes
 .../derived/pair_events_x-451.0mm.npz              |  Bin 0 -> 120956 bytes
 .../derived/pair_events_x-452.0mm.npz              |  Bin 0 -> 120882 bytes
 .../derived/pair_events_x-453.0mm.npz              |  Bin 0 -> 120856 bytes
 .../derived/pair_events_x-454.0mm.npz              |  Bin 0 -> 121828 bytes
 .../derived/pair_events_x-455.0mm.npz              |  Bin 0 -> 122035 bytes
 .../derived/pair_events_x-456.0mm.npz              |  Bin 0 -> 122679 bytes
 .../derived/pair_events_x-457.0mm.npz              |  Bin 0 -> 123063 bytes
 .../derived/pair_events_x-458.0mm.npz              |  Bin 0 -> 123181 bytes
 .../derived/pair_events_x-459.0mm.npz              |  Bin 0 -> 123238 bytes
 .../derived/pair_events_x-460.0mm.npz              |  Bin 0 -> 123352 bytes
 .../derived/pair_events_x-461.0mm.npz              |  Bin 0 -> 123315 bytes
 .../derived/pair_events_x-462.0mm.npz              |  Bin 0 -> 123329 bytes
 .../figures/calibrations_and_residuals.pdf         |  Bin 0 -> 25616 bytes
 .../figures/calibrations_and_residuals.png         |  Bin 0 -> 284318 bytes
 .../figures/pos_ref_1_all_hit_times.pdf            |  Bin 0 -> 23571 bytes
 .../figures/pos_ref_1_all_hit_times.png            |  Bin 0 -> 146556 bytes
 .../figures/pos_ref_1_correlations.pdf             |  Bin 0 -> 25033 bytes
 .../figures/pos_ref_1_correlations.png             |  Bin 0 -> 298577 bytes
 .../figures/pos_ref_1_delta_t.pdf                  |  Bin 0 -> 20414 bytes
 .../figures/pos_ref_1_delta_t.png                  |  Bin 0 -> 113077 bytes
 .../figures/pos_ref_1_event_times.pdf              |  Bin 0 -> 20672 bytes
 .../figures/pos_ref_1_event_times.png              |  Bin 0 -> 89908 bytes
 .../figures/pos_ref_1_npe_moyal.pdf                |  Bin 0 -> 26752 bytes
 .../figures/pos_ref_1_npe_moyal.png                |  Bin 0 -> 180911 bytes
 .../figures/pos_ref_2_all_hit_times.pdf            |  Bin 0 -> 23694 bytes
 .../figures/pos_ref_2_all_hit_times.png            |  Bin 0 -> 149031 bytes
 .../figures/pos_ref_2_correlations.pdf             |  Bin 0 -> 24595 bytes
 .../figures/pos_ref_2_correlations.png             |  Bin 0 -> 275444 bytes
 .../figures/pos_ref_2_delta_t.pdf                  |  Bin 0 -> 19869 bytes
 .../figures/pos_ref_2_delta_t.png                  |  Bin 0 -> 114187 bytes
 .../figures/pos_ref_2_event_times.pdf              |  Bin 0 -> 20491 bytes
 .../figures/pos_ref_2_event_times.png              |  Bin 0 -> 88939 bytes
 .../figures/pos_ref_2_npe_moyal.pdf                |  Bin 0 -> 26630 bytes
 .../figures/pos_ref_2_npe_moyal.png                |  Bin 0 -> 177596 bytes
 .../figures/position_reconstruction.pdf            |  Bin 0 -> 23751 bytes
 .../figures/position_reconstruction.png            |  Bin 0 -> 180735 bytes
 results/exec11_20260612_182454/logs/derive.log     |   41 +
 results/exec11_20260612_182454/logs/detail.log     |    0
 results/exec11_20260612_182454/logs/qa.log         |    0
 .../exec11_20260612_182454/logs/reconstruct.log    |    0
 results/exec11_20260612_182454/logs/report.log     |    0
 .../exec11_20260612_182454/tables/fit_qa_v1_v2.tex |    9 +
 .../tables/reconstruction_summary.tex              |   16 +
 .../tables/reference_comparison.tex                |    8 +
 results/exec12_20260612_191000/README.md           |   62 +
 .../analysis/blue_summary.csv                      |   32 +
 .../analysis/blue_weights.csv                      |   32 +
 .../analysis/configuration_provenance.json         |   29 +
 .../analysis/cv_calibrations.csv                   |  187 +++
 .../analysis/cv_model_selection.csv                |    2 +
 .../analysis/cv_predictions.csv.gz                 |  Bin 0 -> 4316374 bytes
 .../analysis/data_inventory.csv                    |   32 +
 .../analysis/x_reconstruction_summary.csv          |  187 +++
 .../analysis/y0_feasibility_summary.csv            |   94 ++
 .../derived/events/events_x+0mm.npz                |  Bin 0 -> 1673253 bytes
 .../derived/events/events_x+100mm.npz              |  Bin 0 -> 1690764 bytes
 .../derived/events/events_x+150mm.npz              |  Bin 0 -> 1686376 bytes
 .../derived/events/events_x+200mm.npz              |  Bin 0 -> 1675705 bytes
 .../derived/events/events_x+250mm.npz              |  Bin 0 -> 1660646 bytes
 .../derived/events/events_x+300mm.npz              |  Bin 0 -> 1641542 bytes
 .../derived/events/events_x+350mm.npz              |  Bin 0 -> 1619395 bytes
 .../derived/events/events_x+400mm.npz              |  Bin 0 -> 1596329 bytes
 .../derived/events/events_x+450mm.npz              |  Bin 0 -> 1564409 bytes
 .../derived/events/events_x+500mm.npz              |  Bin 0 -> 1525987 bytes
 .../derived/events/events_x+50mm.npz               |  Bin 0 -> 1693544 bytes
 .../derived/events/events_x+550mm.npz              |  Bin 0 -> 1488016 bytes
 .../derived/events/events_x+600mm.npz              |  Bin 0 -> 1444699 bytes
 .../derived/events/events_x+650mm.npz              |  Bin 0 -> 1398577 bytes
 .../derived/events/events_x+670mm.npz              |  Bin 0 -> 1376021 bytes
 .../derived/events/events_x+690mm.npz              |  Bin 0 -> 1364036 bytes
 .../derived/events/events_x-100mm.npz              |  Bin 0 -> 1688855 bytes
 .../derived/events/events_x-150mm.npz              |  Bin 0 -> 1685757 bytes
 .../derived/events/events_x-200mm.npz              |  Bin 0 -> 1674171 bytes
 .../derived/events/events_x-250mm.npz              |  Bin 0 -> 1659670 bytes
 .../derived/events/events_x-300mm.npz              |  Bin 0 -> 1643825 bytes
 .../derived/events/events_x-350mm.npz              |  Bin 0 -> 1619374 bytes
 .../derived/events/events_x-400mm.npz              |  Bin 0 -> 1594561 bytes
 .../derived/events/events_x-450mm.npz              |  Bin 0 -> 1560417 bytes
 .../derived/events/events_x-500mm.npz              |  Bin 0 -> 1526242 bytes
 .../derived/events/events_x-50mm.npz               |  Bin 0 -> 1693474 bytes
 .../derived/events/events_x-550mm.npz              |  Bin 0 -> 1487890 bytes
 .../derived/events/events_x-600mm.npz              |  Bin 0 -> 1439585 bytes
 .../derived/events/events_x-650mm.npz              |  Bin 0 -> 1395798 bytes
 .../derived/events/events_x-670mm.npz              |  Bin 0 -> 1377928 bytes
 .../derived/events/events_x-690mm.npz              |  Bin 0 -> 1363389 bytes
 .../exec12_20260612_191000/figures/bias_vs_x.pdf   |  Bin 0 -> 22417 bytes
 .../figures/calibration_end_ratio.pdf              |  Bin 0 -> 15355 bytes
 .../figures/calibration_end_timing.pdf             |  Bin 0 -> 15330 bytes
 .../figures/calibration_top_centroid.pdf           |  Bin 0 -> 15180 bytes
 .../figures/end_channel_symmetry.pdf               |  Bin 0 -> 13810 bytes
 .../residual_distributions_selected_positions.pdf  |  Bin 0 -> 18551 bytes
 .../exec12_20260612_191000/figures/rms68_vs_x.pdf  |  Bin 0 -> 24309 bytes
 .../figures/sigma_core_vs_x.pdf                    |  Bin 0 -> 24945 bytes
 .../figures/valid_fraction_vs_x.pdf                |  Bin 0 -> 19674 bytes
 .../figures/y_centroid_mean_vs_x.pdf               |  Bin 0 -> 16926 bytes
 .../figures/y_centroid_width_vs_x.pdf              |  Bin 0 -> 17309 bytes
 .../figures/y_left_vs_y_right.pdf                  |  Bin 0 -> 18062 bytes
 results/exec12_20260612_191000/logs/analysis.log   |    0
 results/exec12_20260612_191000/logs/inventory.log  |   31 +
 .../report/endtop_position_reconstruction.md       |   59 +
 .../report/endtop_position_reconstruction.tex      |   60 +
 .../report/exec11_technical_note.md                |   42 +
 .../report/exec11_technical_note.tex               |   43 +
 .../tables/x_reconstruction_summary.tex            |  190 +++
 results/exec12t_20260612_195426/README.md          |   48 +
 .../analysis/calibration_20pe.csv                  |    2 +
 .../analysis/calibration_4pe.csv                   |    2 +
 .../analysis/configuration_provenance.json         |   10 +
 .../analysis/data_inventory.csv                    |   42 +
 .../analysis/exec11_reproduction_check.csv         |   42 +
 .../analysis/loo_predictions_20pe.npz              |  Bin 0 -> 892931 bytes
 .../analysis/loo_predictions_4pe.npz               |  Bin 0 -> 892788 bytes
 .../analysis/position_reconstruction_20pe.csv      |   42 +
 .../analysis/position_reconstruction_4pe.csv       |   42 +
 .../analysis/temporal_position_summary.csv         |   83 ++
 .../analysis/threshold_4_20_summary.csv            |    5 +
 .../analysis/threshold_sweep_summary.csv           |   31 +
 .../analysis/window_dip_summary.csv                |    9 +
 .../beamer/beamer_contact_sheet.png                |  Bin 0 -> 473675 bytes
 .../beamer/exec12t_timing_position_beamer.aux      |  128 ++
 .../beamer/exec12t_timing_position_beamer.log      | 1406 ++++++++++++++++++++
 .../beamer/exec12t_timing_position_beamer.nav      |  109 ++
 .../beamer/exec12t_timing_position_beamer.out      |    1 +
 .../beamer/exec12t_timing_position_beamer.pdf      |  Bin 0 -> 234763 bytes
 .../beamer/exec12t_timing_position_beamer.snm      |    0
 .../beamer/exec12t_timing_position_beamer.tex      |  166 +++
 .../beamer/exec12t_timing_position_beamer.toc      |    0
 .../exec12t_20260612_195426/beamer/references.bib  |    1 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 29095 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 47777 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 19663 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 27031 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 15994 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 20580 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 32086 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 26819 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 26911 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 20004 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 20223 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 23610 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 24670 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 22094 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 30569 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 30458 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 18922 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 29235 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 15229 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 19449 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 20846 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 25617 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 24429 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 19255 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 19840 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 24193 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 24012 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 23457 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 31499 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 33489 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 19206 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 27919 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 15664 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 18758 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 21653 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 25116 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 46153 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 21146 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 25827 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 26188 bytes
 .../beamer/rendered/slide-41.png                   |  Bin 0 -> 18935 bytes
 .../beamer/rendered/slide-42.png                   |  Bin 0 -> 19827 bytes
 .../beamer/rendered/slide-43.png                   |  Bin 0 -> 25075 bytes
 .../beamer/rendered/slide-44.png                   |  Bin 0 -> 25549 bytes
 .../beamer/rendered/slide-45.png                   |  Bin 0 -> 22175 bytes
 .../beamer/rendered/slide-46.png                   |  Bin 0 -> 33040 bytes
 .../beamer/rendered/slide-47.png                   |  Bin 0 -> 32296 bytes
 .../beamer/rendered/slide-48.png                   |  Bin 0 -> 20505 bytes
 .../beamer/rendered/slide-49.png                   |  Bin 0 -> 28637 bytes
 .../beamer/speaker_notes.md                        |  107 ++
 .../derived/order_statistics/pair_x-422.0mm.npz    |  Bin 0 -> 1362830 bytes
 .../derived/order_statistics/pair_x-423.0mm.npz    |  Bin 0 -> 1362991 bytes
 .../derived/order_statistics/pair_x-424.0mm.npz    |  Bin 0 -> 1362994 bytes
 .../derived/order_statistics/pair_x-425.0mm.npz    |  Bin 0 -> 1362676 bytes
 .../derived/order_statistics/pair_x-426.0mm.npz    |  Bin 0 -> 1362617 bytes
 .../derived/order_statistics/pair_x-427.0mm.npz    |  Bin 0 -> 1362456 bytes
 .../derived/order_statistics/pair_x-428.0mm.npz    |  Bin 0 -> 1362417 bytes
 .../derived/order_statistics/pair_x-429.0mm.npz    |  Bin 0 -> 1362477 bytes
 .../derived/order_statistics/pair_x-430.0mm.npz    |  Bin 0 -> 1362246 bytes
 .../derived/order_statistics/pair_x-431.0mm.npz    |  Bin 0 -> 1361625 bytes
 .../derived/order_statistics/pair_x-432.0mm.npz    |  Bin 0 -> 1361227 bytes
 .../derived/order_statistics/pair_x-433.0mm.npz    |  Bin 0 -> 1361512 bytes
 .../derived/order_statistics/pair_x-434.0mm.npz    |  Bin 0 -> 1361634 bytes
 .../derived/order_statistics/pair_x-435.0mm.npz    |  Bin 0 -> 1361842 bytes
 .../derived/order_statistics/pair_x-436.0mm.npz    |  Bin 0 -> 1361904 bytes
 .../derived/order_statistics/pair_x-437.0mm.npz    |  Bin 0 -> 1362200 bytes
 .../derived/order_statistics/pair_x-438.0mm.npz    |  Bin 0 -> 1362269 bytes
 .../derived/order_statistics/pair_x-439.0mm.npz    |  Bin 0 -> 1362081 bytes
 .../derived/order_statistics/pair_x-440.0mm.npz    |  Bin 0 -> 1361731 bytes
 .../derived/order_statistics/pair_x-441.0mm.npz    |  Bin 0 -> 1361645 bytes
 .../derived/order_statistics/pair_x-442.0mm.npz    |  Bin 0 -> 1361574 bytes
 .../derived/order_statistics/pair_x-443.0mm.npz    |  Bin 0 -> 1361936 bytes
 .../derived/order_statistics/pair_x-444.0mm.npz    |  Bin 0 -> 1361808 bytes
 .../derived/order_statistics/pair_x-445.0mm.npz    |  Bin 0 -> 1362189 bytes
 .../derived/order_statistics/pair_x-446.0mm.npz    |  Bin 0 -> 1362266 bytes
 .../derived/order_statistics/pair_x-447.0mm.npz    |  Bin 0 -> 1362136 bytes
 .../derived/order_statistics/pair_x-448.0mm.npz    |  Bin 0 -> 1362063 bytes
 .../derived/order_statistics/pair_x-449.0mm.npz    |  Bin 0 -> 1362102 bytes
 .../derived/order_statistics/pair_x-450.0mm.npz    |  Bin 0 -> 1361698 bytes
 .../derived/order_statistics/pair_x-451.0mm.npz    |  Bin 0 -> 1361068 bytes
 .../derived/order_statistics/pair_x-452.0mm.npz    |  Bin 0 -> 1361575 bytes
 .../derived/order_statistics/pair_x-453.0mm.npz    |  Bin 0 -> 1361701 bytes
 .../derived/order_statistics/pair_x-454.0mm.npz    |  Bin 0 -> 1362028 bytes
 .../derived/order_statistics/pair_x-455.0mm.npz    |  Bin 0 -> 1362411 bytes
 .../derived/order_statistics/pair_x-456.0mm.npz    |  Bin 0 -> 1362338 bytes
 .../derived/order_statistics/pair_x-457.0mm.npz    |  Bin 0 -> 1362651 bytes
 .../derived/order_statistics/pair_x-458.0mm.npz    |  Bin 0 -> 1362862 bytes
 .../derived/order_statistics/pair_x-459.0mm.npz    |  Bin 0 -> 1362956 bytes
 .../derived/order_statistics/pair_x-460.0mm.npz    |  Bin 0 -> 1362783 bytes
 .../derived/order_statistics/pair_x-461.0mm.npz    |  Bin 0 -> 1362770 bytes
 .../derived/order_statistics/pair_x-462.0mm.npz    |  Bin 0 -> 1362960 bytes
 .../environment_exec12t.txt                        |   99 ++
 .../figures/correlation_ab_vs_x.pdf                |  Bin 0 -> 12608 bytes
 .../figures/efficiency_4_20_vs_x.pdf               |  Bin 0 -> 12461 bytes
 .../figures/mean_delta_4_20_vs_x.pdf               |  Bin 0 -> 13853 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 9467 bytes
 .../figures/pos_ref_1_delta_t_4_20.pdf             |  Bin 0 -> 15320 bytes
 .../figures/pos_ref_2_delta_t_4_20.pdf             |  Bin 0 -> 15348 bytes
 .../figures/threshold_bias.pdf                     |  Bin 0 -> 15585 bytes
 .../figures/threshold_chi2.pdf                     |  Bin 0 -> 15384 bytes
 .../figures/threshold_efficiency.pdf               |  Bin 0 -> 14268 bytes
 .../figures/threshold_pareto.pdf                   |  Bin 0 -> 14148 bytes
 .../figures/threshold_sigma_dt.pdf                 |  Bin 0 -> 15598 bytes
 .../figures/threshold_sigma_x.pdf                  |  Bin 0 -> 15710 bytes
 .../figures/threshold_slope.pdf                    |  Bin 0 -> 15487 bytes
 .../figures/timing_width_4_20_vs_x.pdf             |  Bin 0 -> 15411 bytes
 results/exec12t_20260612_195426/logs/analysis.log  |    0
 .../logs/beamer_pdfinfo.txt                        |   20 +
 .../logs/beamer_warnings.txt                       |    1 +
 .../logs/exec11_reproduction.log                   |    3 +
 results/exec12t_20260612_195426/logs/inventory.log |   41 +
 results/exec12t_20260612_195426/logs/window.log    |    0
 .../report/exec12t_timing_position_report.aux      |   26 +
 .../report/exec12t_timing_position_report.log      |  376 ++++++
 .../report/exec12t_timing_position_report.md       |   45 +
 .../report/exec12t_timing_position_report.out      |    7 +
 .../report/exec12t_timing_position_report.pdf      |  Bin 0 -> 35001 bytes
 .../report/exec12t_timing_position_report.tex      |   63 +
 .../tables/generated_numbers.tex                   |    6 +
 .../tables/global_context.csv                      |    7 +
 .../tables/global_context.tex                      |   10 +
 .../tables/threshold_4_20_summary.csv              |    7 +
 .../tables/threshold_4_20_summary.tex              |   10 +
 .../tables/threshold_sweep.csv                     |   31 +
 .../tables/threshold_sweep.tex                     |   34 +
 .../tables/window_dip_summary.csv                  |    9 +
 .../tables/window_dip_summary.tex                  |   12 +
 .../beamer/exec12tb_beamer.pdf                     |  Bin 0 -> 544006 bytes
 .../beamer/exec12tb_beamer.tex                     |  971 ++++++++++++++
 .../beamer/generated_numbers.tex                   |   24 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 37541 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 55903 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 39821 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 44919 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 37183 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 50484 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 53853 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 54725 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 55833 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 57120 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 55518 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 50060 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 46531 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 53957 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 45247 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 58873 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 56682 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 50931 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 57910 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 49913 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 63029 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 66176 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 49963 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 45261 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 56098 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 68087 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 72976 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 71762 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 45414 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 8355 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 54362 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 54841 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 40412 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 46539 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 61442 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 54635 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 47345 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 50838 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 47245 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 44446 bytes
 .../beamer/speaker_notes.md                        |  227 ++++
 .../figures/cfd_vs_orderstat_schematic.pdf         |  Bin 0 -> 32384 bytes
 .../figures/fine_scan_geometry.pdf                 |  Bin 0 -> 18315 bytes
 .../figures/global_context_sigma_bias.pdf          |  Bin 0 -> 21642 bytes
 .../figures/mean_dt_vs_x_with_residuals.pdf        |  Bin 0 -> 25199 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 21280 bytes
 .../figures/position_resolution_4_20_vs_x.pdf      |  Bin 0 -> 25708 bytes
 .../figures/ref1_dt_overlay.pdf                    |  Bin 0 -> 25821 bytes
 .../figures/ref1_xrec_overlay.pdf                  |  Bin 0 -> 26514 bytes
 .../figures/ref2_dt_overlay.pdf                    |  Bin 0 -> 25730 bytes
 .../figures/ref2_xrec_overlay.pdf                  |  Bin 0 -> 26561 bytes
 .../figures/rho_ab_vs_x_4_20.pdf                   |  Bin 0 -> 22331 bytes
 .../figures/rms68_dt_vs_x_4_20.pdf                 |  Bin 0 -> 21590 bytes
 .../figures/sha_manifest.json                      |   22 +
 .../figures/sweep_bias_vs_k.pdf                    |  Bin 0 -> 17260 bytes
 .../figures/sweep_chi2_vs_k.pdf                    |  Bin 0 -> 19883 bytes
 .../figures/sweep_rms68x_vs_k.pdf                  |  Bin 0 -> 20698 bytes
 .../figures/sweep_slope_vs_k.pdf                   |  Bin 0 -> 17743 bytes
 .../figures/tplus_rms68_vs_x_4_20.pdf              |  Bin 0 -> 22163 bytes
 .../figures/tradeoff_resolution_vs_bias.pdf        |  Bin 0 -> 27883 bytes
 .../figures/window_dip_counts.pdf                  |  Bin 0 -> 16254 bytes
 .../figures/window_dip_t4_t20.pdf                  |  Bin 0 -> 18771 bytes
 .../logs/beamer_compile1.log                       |  255 ++++
 .../logs/beamer_compile2.log                       |  265 ++++
 .../logs/beamer_compile3.log                       |  275 ++++
 .../logs/beamer_compile4.log                       |  618 +++++++++
 .../logs/beamer_compile5.log                       |  598 +++++++++
 results/exec12tb_20260612_204216/logs/visual_qa.md |   60 +
 .../manifest/beamer_manifest.csv                   |   41 +
 .../tables/dataset_inventory.tex                   |   15 +
 .../tables/global_context.tex                      |   18 +
 .../tables/lv_comparison.tex                       |   23 +
 .../tables/reference_positions.tex                 |   14 +
 .../tables/threshold_4_20_comparison.tex           |   22 +
 .../tables/threshold_sweep.tex                     |   22 +
 .../tables/window_dip_summary.tex                  |   16 +
 scripts/analyze_geometry.py                        |  172 +++
 scripts/build_exec12t.sh                           |   13 +
 scripts/gen_pair_scan_macros.py                    |  119 ++
 scripts/make_contact_sheet.py                      |   36 +
 scripts/run_pair_scan.sh                           |  156 +++
 src/DetectorConstruction.cc                        |   73 +-
 tests/check_endtop_balance.py                      |   14 +-
 tests/check_pairscan_geometry.py                   |  123 ++
 tests/check_pairscan_macros.py                     |   81 ++
 tests/readout_config_check.cc                      |   53 +-
 tests/test_exec11_pair_analysis.py                 |   64 +
 tests/test_exec12_endtop_position.py               |   20 +
 tests/test_exec12t_timing_threshold_analysis.py    |   15 +
 444 files changed, 13222 insertions(+), 128 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E127](#e127)
<a id="b4"></a>

### 4. feat/bar-end-vikuiti
Inspected ref `feat/bar-end-vikuiti`, full SHA `50cf02e34cd7dddfb847739f6249afd1a600478f`. Classification: **Merged ancestor; pre-EXEC_25 deck**. [E032](#e032)
**B:** ahead 0, behind 74; merged YES; contains 8349041 NO. [E052](#e052) [E223](#e223) [E035](#e035)
Upstream: `origin/feat/bar-end-vikuiti` [ahead 13]. Origin counterpart: `1bfb82743270dcc3ea9b17be03a6af5704cc6850`. [E030](#e030)
Tip date / author / subject: 2026-09-01 12:27:32 +0200 / rrios / docs(presentations): add master index presentations/README.md. [E032](#e032)
**C1/C2:** active configuration: Air-gap + polished dielectric_dielectric reflector border. Reflectivity: 0.95 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E128](#e128) [E129](#e129) [E296](#e296)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.95, 0.95};
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 HIT**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **118**. Full paths, line numbers and contents are preserved in [E130](#e130) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E129](#e129)
```text
311:     {
312:         auto* scintAirSurface     = Materials::CreateBarSurface();
313:         auto* airReflectorSurface = Materials::CreateBarSkinReflector();
314:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E128](#e128) [E296](#e296)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
332:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
354:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E298](#e298). No matches outside the bundled external libraries.
Deck `presentations/v6/talk_v6.tex`: **514.9 MISS**, **1300 MISS**. [E297](#e297)

```text
33: \newcommand{\Rvikuiti}{0.95}
35: \newcommand{\NpeG}{701.3}
82:     \textbf{Reflector:} Vikuiti ESR on all non-SiPM surfaces ($R = \Rvikuiti$).\\[5pt]
415:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
756:       $N_\text{pe}$/end at $x=0$ & $701.3$ \\
916:     $N_\text{pe}$ & $\sim\!0.37$ (BUG) & — & 701.3 \\
963:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
995:       \textbf{G4/end} & \textbf{701.3} & +refl.-rec. \\
```

Exact loose `0.95` occurrences in the deck: **16**. **OPEN-07-style macro gap:** historical literals remain; the requested final macro cleanup is absent. This label describes the requested audit check, not a verified issue identifier.
```text
33: \newcommand{\Rvikuiti}{0.95}
200:       \item Reflected ($R=0.95$ per bounce) back into bar
251:     At $L/2=700\mm$: survival $= 0.95^{85.6} \approx 1.2\%$ (code $R=0.95$). \\[4pt]
363:       \item Air$\to$Mylar: \texttt{CreateBarSkinReflector()} as border surface — \texttt{dielectric\_dielectric}, $R=0.95$ constant
741:       \item Reflector: \textbf{Vikuiti ESR} ($R=0.95$ constant in code)
790:       \item Vikuiti (R=0.95) recovers non-TIR photons (+23\%)
845:       \item $R = 0.95$ constant (wavelength-independent)
857:     \textbf{A2 — R=0.95 vs R=0.98}\\[4pt]
858:     Code value (Materials.cc): $R = 0.95$\\
860:     The code comment says ``R=0.98 Vikuiti ESR'' but the actual constant set is 0.95.
861:     All simulation results in this analysis use \alert{$R=0.95$}.\\[6pt]
864:     With R=0.95: $\Lambda_\text{refl}^H = 240\mm$ (shorter — more loss).\\[4pt]
887:       \item Vikuiti reflectivity at air–Mylar interface ($R=0.95$)
914:     Reflector & \texttt{dielectric\_metal skin} & skin surface, R=0.95 & border surface, R=0.95 \\
939:     \textbf{Fix (exec21-optfix):} Changed to \texttt{dielectric\_dielectric} border surface with $R=0.95$, preserving natural TIR at bar–air interface.
1192:       Reflector & Vikuiti ESR ($R=0.95$) \\
```

701.3 contextual caution: **SUSPICIOUS**: the old deck retains 701.3 before the hybrid correction; inspect the design frame below.

```text
32: \newcommand{\LambdaH}{405\mm}
33: \newcommand{\Rvikuiti}{0.95}
34: \newcommand{\NpeNapkin}{570}
35: \newcommand{\NpeG}{701.3}
36: \newcommand{\sigmaENDval}{53.68\ps}
37: \newcommand{\sigmaTOPval}{15.20\ps}
38: \newcommand{\sigmaBLUEval}{15.21\ps}
412:       \midrule
413:       TIR-guided / end & 570 & — \\
414:       Reflector-recovered & 0 & — \\
415:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
416:       \bottomrule
417:     \end{tabular}\\[6pt]
418:     $\Delta N_\text{pe}/N_\text{nap} = +23\%$.\\[4pt]
732: \section{Design Summary}
733: % ====================================================================
734: 
735: \begin{frame}{Design Decision: EJ-230 + 20 TOP + 16 END}
736:   \begin{columns}[T]
737:     \column{0.52\linewidth}
738:     \textbf{Selected configuration:}
753:       $\sigma_t$ at $x=0$ & $\mathbf{15.2\ps}$ \\
754:       Mean $\sigma_t$ across bar & $19.5\ps$ \\
755:       $\sigma_x$ (END $\Delta t$) & $\mathbf{7.9\mm}$ \\
756:       $N_\text{pe}$/end at $x=0$ & $701.3$ \\
757:       \bottomrule
758:     \end{tabular}
759:     \end{center}
913:     Air gap & No & No & Yes (0.10 mm) \\
914:     Reflector & \texttt{dielectric\_metal skin} & skin surface, R=0.95 & border surface, R=0.95 \\
915:     TIR & Eliminated by skin & Natural (no bar surface) & Natural Fresnel \\
916:     $N_\text{pe}$ & $\sim\!0.37$ (BUG) & — & 701.3 \\
917:     $\sigma_\text{END}$ & N/A (insufficient pe) & 50.5 ps & 53.7 ps \\
918:     TOP & No & No & Yes ($N=4,8,14,20$) \\
919:     Status & \textcolor{red}{\textbf{INVALID}} & superseded & \textcolor{green!60!black}{\textbf{current}} \\
960:     Yield & 10000 ph/MeV & 10400 ph/MeV & 9700 ph/MeV \\
961:     $n$ & 1.58 & 1.58 & 1.58 \\
962:     \midrule
963:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
964:     $N_\text{pe}$/end (napkin) & — & — & 570 \\
965:     $\sigma_t$ END-only (G4) & 47.6 ps & 52.1 ps & 49.6 ps \\
966:     Napkin bulk survival & 0.832 & 0.646 & \textbf{0.558} \\
992:       PDE & 0.400 & assumed \\
993:       \textbf{Napkin/end} & \textbf{570} & TIR-only \\
994:       \midrule
995:       \textbf{G4/end} & \textbf{701.3} & +refl.-rec. \\
996:       Surplus & $+23\%$ & Pop.\ II \\
997:       \bottomrule
998:     \end{tabular}}
```

**C4:** **0/14 complete exact-suffix figure trios.** [E055](#e055)
| Figure stem (relative to deck figs/) | Figure PDF | .root | .csv | .meta.json | _meta.json alternative | Matching sidecars elsewhere in tree |
|---|---|---|---|---|---|---|
| `fig1_kscan` | MISS | MISS | MISS | MISS | MISS | none |
| `figM1_mat_sigma_end` | MISS | MISS | MISS | MISS | MISS | none |
| `figM2_mat_npe` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_bulk_survival` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_end_mscan` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_npe_x` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_ntop_scan_full` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_sigma_t_x` | MISS | MISS | MISS | MISS | MISS | none |
| `fig_survival_RN` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_pareto` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_sigma_vs_x` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_top_position_loo` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_veff_fit` | MISS | MISS | MISS | MISS | MISS | none |
| `v5_veff_residual` | MISS | MISS | MISS | MISS | MISS | none |
Figure list is derived from literal `\anafig{...}` and `\includegraphics{figs/...}` calls in the tracked TeX. No dynamic figure paths were assumed; TeX build execution and contents inside figure PDFs were not tested.
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E053](#e053)
```text
(none)
```

**Files touched on the branch since its merge base with main:** [E054](#e054)
```text
(none)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E131](#e131)
<a id="b5"></a>

### 5. feat/bar-vikuiti
Inspected ref `feat/bar-vikuiti`, full SHA `219fbe385d43629c523c29b4a8d550783b469f35`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 36, behind 115; merged NO; contains 8349041 NO. [E057](#e057) [E223](#e223) [E035](#e035)
Upstream: none configured. Origin counterpart: `219fbe385d43629c523c29b4a8d550783b469f35`. [E030](#e030)
Tip date / author / subject: 2026-08-15 23:53:49 +0200 / rrios / feat(scan): add comparative scan analysis and Beamer for 4 bar configs. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal polished bar skin. Reflectivity: 0.98 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E132](#e132) [E133](#e133) [E299](#e299)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **6**. Full paths, line numbers and contents are preserved in [E134](#e134) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E133](#e133)
```text
292:     // Reflects 98% of photons at all angles (no TIR threshold).
293:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV,
294:                                               Materials::CreateBarSkinReflector());
295:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E132](#e132) [E299](#e299)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E300](#e300). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E060](#e060)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E060](#e060)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E058](#e058)
```text
219fbe3 feat(scan): add comparative scan analysis and Beamer for 4 bar configs
8bf2215 feat(vikuiti): add Vikuiti 3M reflector branch, scan script, and analysis macro
b281aea feat(beamer): add EJ-204 vs EJ-230 TIR-only bar timing comparison report
a95f610 feat(analysis): add bar timing resolution analysis and run macro
5f30646 feat(ej230-bar-tir-only): EJ-230 scintillator, TIR-only bar
234e0cb feat(ej204-bar-tir-only): EJ-204 bar TIR-only — remove Mylar skin reflector
f19c093 feat(pairscan): CMakeLists refactor, pair scan macros and run script
14ae395 fix(tests): replace escape-fraction guard with sipm-entry guard
c7a627e fix(optics): replace reflector volumes with bar skin surface
021d5f4 EXEC_12TB: add rebuilt self-contained beamer and speaker notes
1ac830a EXEC_12TB: add beamer manifest and Makefile exec12tb targets
d764892 EXEC_12TB: add generated tables and TeX number macros
9cfe042 EXEC_12TB: add figure style module and regenerated deck figures
66ac24c EXEC_12T: add complete Beamer presentation and reproducible build
b0aaac1 EXEC_12T: add integrated timing-position technical report
a2563d9 EXEC_12T: add window-dip timing and collection study
aaabe89 EXEC_12T: add threshold sweep and covariance timing decomposition
0a41034 EXEC_12T: add 4PE versus 20PE temporal and position analysis
f38c9e1 EXEC_12T: add cached fourth-to-thirtieth hit order statistics
003d106 EXEC_12: add staged XY scan proposal
6e01606 EXEC_12: add EXEC11 and EndTop technical reports
30166cd EXEC_12: add y-zero observability feasibility study
a9571fd EXEC_12: add covariance-aware EndTop estimator combination
bd21563 EXEC_12: add leave-one-position-out global X reconstruction
8332733 EXEC_12: add leave-one-position-out global X reconstruction
17d9e7d EXEC_12: add EndTop event observable builder and geometry tests
4dba735 EXEC_11: add reproducible report tables and README
a103ed6 EXEC_11: add temporal ratio and covariance-aware position reconstruction
1e95ac6 EXEC_11: add detailed two-position timing analysis
61bca2e EXEC_11: regenerate pair-scan fit QA and summary v2
bf0d2d5 EXEC_11: add per-event pair observable builder and tests
f431c01 EXEC_08/Step7: CTest guardrails for pair-scan geometry and macros
20bed82 EXEC_08/Step6: batch runner scripts/run_pair_scan.sh
549f3a6 EXEC_08/Step5: ROOT analysis macro analyze_pair_scan.C
30f3da3 EXEC_08/Step3: generate 41 pair-scan macros in macros/pairscan/
6bb73d5 EXEC_08/Step2: geometry analysis + pair selection
```

**Files touched on the branch since its merge base with main:** [E059](#e059)
```text
 .gitignore                                         |    7 +
 CMakeLists.txt                                     |   76 +-
 Makefile                                           |   40 +
 analysis/analyze_pair_scan.C                       |  496 +++++++
 analysis/bar_resolution_scan.cxx                   |  396 ++++++
 analysis/bar_timing_resolution.cxx                 |  387 ++++++
 analysis/exec11_pair_analysis.py                   |  934 +++++++++++++
 analysis/exec12_endtop_position.py                 |  161 +++
 analysis/exec12_make_reports.py                    |   99 ++
 analysis/exec12t_make_products.py                  |  151 +++
 analysis/exec12t_timing_threshold_analysis.py      |  136 ++
 analysis/exec12tb_figstyle.py                      |   56 +
 analysis/exec12tb_figures.py                       |  823 ++++++++++++
 analysis/exec12tb_tables.py                        |  487 +++++++
 beamer/bar_scan_comparative.pdf                    |  Bin 0 -> 282454 bytes
 beamer/bar_scan_comparative.tex                    |  637 +++++++++
 beamer/ej204_vs_ej230_bar_tir.pdf                  |  Bin 0 -> 412748 bytes
 beamer/ej204_vs_ej230_bar_tir.tex                  |  572 ++++++++
 include/DetectorConstruction.hh                    |   10 +-
 macros/pairscan/pairscan_x-422.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-423.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-424.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-425.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-426.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-427.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-428.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-429.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-430.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-431.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-432.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-433.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-434.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-435.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-436.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-437.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-438.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-439.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-440.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-441.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-442.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-443.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-444.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-445.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-446.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-447.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-448.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-449.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-450.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-451.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-452.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-453.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-454.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-455.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-456.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-457.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-458.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-459.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-460.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-461.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-462.0mm.mac             |   20 +
 macros/timing_tir_center.mac                       |   27 +
 pair_scan_config.json                              |   71 +
 report/exec13_xy_scan_plan.md                      |   30 +
 results/exec11_20260612_182454/README.md           |  124 ++
 .../analysis/calibration_ratio.csv                 |    3 +
 .../analysis/calibration_temporal.csv              |    2 +
 .../analysis/calibration_v1_v2.csv                 |    3 +
 .../analysis/data_inventory.csv                    |   42 +
 results/exec11_20260612_182454/analysis/fit_qa.csv |   42 +
 .../analysis/fit_qa_v1_v2_focus.csv                |    4 +
 .../exec11_20260612_182454/analysis/metadata.json  |   13 +
 .../analysis/pairscan_summary_v2.csv               |   42 +
 .../analysis/reconstruction_summary.csv            |   11 +
 .../analysis/reference_comparison.csv              |    3 +
 .../analysis/reference_positions.csv               |    3 +
 .../derived/pair_events_x-422.0mm.npz              |  Bin 0 -> 123534 bytes
 .../derived/pair_events_x-423.0mm.npz              |  Bin 0 -> 122740 bytes
 .../derived/pair_events_x-424.0mm.npz              |  Bin 0 -> 123199 bytes
 .../derived/pair_events_x-425.0mm.npz              |  Bin 0 -> 123147 bytes
 .../derived/pair_events_x-426.0mm.npz              |  Bin 0 -> 123368 bytes
 .../derived/pair_events_x-427.0mm.npz              |  Bin 0 -> 122917 bytes
 .../derived/pair_events_x-428.0mm.npz              |  Bin 0 -> 122839 bytes
 .../derived/pair_events_x-429.0mm.npz              |  Bin 0 -> 122318 bytes
 .../derived/pair_events_x-430.0mm.npz              |  Bin 0 -> 121727 bytes
 .../derived/pair_events_x-431.0mm.npz              |  Bin 0 -> 120954 bytes
 .../derived/pair_events_x-432.0mm.npz              |  Bin 0 -> 120412 bytes
 .../derived/pair_events_x-433.0mm.npz              |  Bin 0 -> 121077 bytes
 .../derived/pair_events_x-434.0mm.npz              |  Bin 0 -> 122386 bytes
 .../derived/pair_events_x-435.0mm.npz              |  Bin 0 -> 122900 bytes
 .../derived/pair_events_x-436.0mm.npz              |  Bin 0 -> 122357 bytes
 .../derived/pair_events_x-437.0mm.npz              |  Bin 0 -> 122009 bytes
 .../derived/pair_events_x-438.0mm.npz              |  Bin 0 -> 122595 bytes
 .../derived/pair_events_x-439.0mm.npz              |  Bin 0 -> 122198 bytes
 .../derived/pair_events_x-440.0mm.npz              |  Bin 0 -> 122232 bytes
 .../derived/pair_events_x-441.0mm.npz              |  Bin 0 -> 122736 bytes
 .../derived/pair_events_x-442.0mm.npz              |  Bin 0 -> 122672 bytes
 .../derived/pair_events_x-443.0mm.npz              |  Bin 0 -> 122553 bytes
 .../derived/pair_events_x-444.0mm.npz              |  Bin 0 -> 122542 bytes
 .../derived/pair_events_x-445.0mm.npz              |  Bin 0 -> 122293 bytes
 .../derived/pair_events_x-446.0mm.npz              |  Bin 0 -> 122306 bytes
 .../derived/pair_events_x-447.0mm.npz              |  Bin 0 -> 121906 bytes
 .../derived/pair_events_x-448.0mm.npz              |  Bin 0 -> 122022 bytes
 .../derived/pair_events_x-449.0mm.npz              |  Bin 0 -> 122871 bytes
 .../derived/pair_events_x-450.0mm.npz              |  Bin 0 -> 122284 bytes
 .../derived/pair_events_x-451.0mm.npz              |  Bin 0 -> 120956 bytes
 .../derived/pair_events_x-452.0mm.npz              |  Bin 0 -> 120882 bytes
 .../derived/pair_events_x-453.0mm.npz              |  Bin 0 -> 120856 bytes
 .../derived/pair_events_x-454.0mm.npz              |  Bin 0 -> 121828 bytes
 .../derived/pair_events_x-455.0mm.npz              |  Bin 0 -> 122035 bytes
 .../derived/pair_events_x-456.0mm.npz              |  Bin 0 -> 122679 bytes
 .../derived/pair_events_x-457.0mm.npz              |  Bin 0 -> 123063 bytes
 .../derived/pair_events_x-458.0mm.npz              |  Bin 0 -> 123181 bytes
 .../derived/pair_events_x-459.0mm.npz              |  Bin 0 -> 123238 bytes
 .../derived/pair_events_x-460.0mm.npz              |  Bin 0 -> 123352 bytes
 .../derived/pair_events_x-461.0mm.npz              |  Bin 0 -> 123315 bytes
 .../derived/pair_events_x-462.0mm.npz              |  Bin 0 -> 123329 bytes
 .../figures/calibrations_and_residuals.pdf         |  Bin 0 -> 25616 bytes
 .../figures/calibrations_and_residuals.png         |  Bin 0 -> 284318 bytes
 .../figures/pos_ref_1_all_hit_times.pdf            |  Bin 0 -> 23571 bytes
 .../figures/pos_ref_1_all_hit_times.png            |  Bin 0 -> 146556 bytes
 .../figures/pos_ref_1_correlations.pdf             |  Bin 0 -> 25033 bytes
 .../figures/pos_ref_1_correlations.png             |  Bin 0 -> 298577 bytes
 .../figures/pos_ref_1_delta_t.pdf                  |  Bin 0 -> 20414 bytes
 .../figures/pos_ref_1_delta_t.png                  |  Bin 0 -> 113077 bytes
 .../figures/pos_ref_1_event_times.pdf              |  Bin 0 -> 20672 bytes
 .../figures/pos_ref_1_event_times.png              |  Bin 0 -> 89908 bytes
 .../figures/pos_ref_1_npe_moyal.pdf                |  Bin 0 -> 26752 bytes
 .../figures/pos_ref_1_npe_moyal.png                |  Bin 0 -> 180911 bytes
 .../figures/pos_ref_2_all_hit_times.pdf            |  Bin 0 -> 23694 bytes
 .../figures/pos_ref_2_all_hit_times.png            |  Bin 0 -> 149031 bytes
 .../figures/pos_ref_2_correlations.pdf             |  Bin 0 -> 24595 bytes
 .../figures/pos_ref_2_correlations.png             |  Bin 0 -> 275444 bytes
 .../figures/pos_ref_2_delta_t.pdf                  |  Bin 0 -> 19869 bytes
 .../figures/pos_ref_2_delta_t.png                  |  Bin 0 -> 114187 bytes
 .../figures/pos_ref_2_event_times.pdf              |  Bin 0 -> 20491 bytes
 .../figures/pos_ref_2_event_times.png              |  Bin 0 -> 88939 bytes
 .../figures/pos_ref_2_npe_moyal.pdf                |  Bin 0 -> 26630 bytes
 .../figures/pos_ref_2_npe_moyal.png                |  Bin 0 -> 177596 bytes
 .../figures/position_reconstruction.pdf            |  Bin 0 -> 23751 bytes
 .../figures/position_reconstruction.png            |  Bin 0 -> 180735 bytes
 results/exec11_20260612_182454/logs/derive.log     |   41 +
 results/exec11_20260612_182454/logs/detail.log     |    0
 results/exec11_20260612_182454/logs/qa.log         |    0
 .../exec11_20260612_182454/logs/reconstruct.log    |    0
 results/exec11_20260612_182454/logs/report.log     |    0
 .../exec11_20260612_182454/tables/fit_qa_v1_v2.tex |    9 +
 .../tables/reconstruction_summary.tex              |   16 +
 .../tables/reference_comparison.tex                |    8 +
 results/exec12_20260612_191000/README.md           |   62 +
 .../analysis/blue_summary.csv                      |   32 +
 .../analysis/blue_weights.csv                      |   32 +
 .../analysis/configuration_provenance.json         |   29 +
 .../analysis/cv_calibrations.csv                   |  187 +++
 .../analysis/cv_model_selection.csv                |    2 +
 .../analysis/cv_predictions.csv.gz                 |  Bin 0 -> 4316374 bytes
 .../analysis/data_inventory.csv                    |   32 +
 .../analysis/x_reconstruction_summary.csv          |  187 +++
 .../analysis/y0_feasibility_summary.csv            |   94 ++
 .../derived/events/events_x+0mm.npz                |  Bin 0 -> 1673253 bytes
 .../derived/events/events_x+100mm.npz              |  Bin 0 -> 1690764 bytes
 .../derived/events/events_x+150mm.npz              |  Bin 0 -> 1686376 bytes
 .../derived/events/events_x+200mm.npz              |  Bin 0 -> 1675705 bytes
 .../derived/events/events_x+250mm.npz              |  Bin 0 -> 1660646 bytes
 .../derived/events/events_x+300mm.npz              |  Bin 0 -> 1641542 bytes
 .../derived/events/events_x+350mm.npz              |  Bin 0 -> 1619395 bytes
 .../derived/events/events_x+400mm.npz              |  Bin 0 -> 1596329 bytes
 .../derived/events/events_x+450mm.npz              |  Bin 0 -> 1564409 bytes
 .../derived/events/events_x+500mm.npz              |  Bin 0 -> 1525987 bytes
 .../derived/events/events_x+50mm.npz               |  Bin 0 -> 1693544 bytes
 .../derived/events/events_x+550mm.npz              |  Bin 0 -> 1488016 bytes
 .../derived/events/events_x+600mm.npz              |  Bin 0 -> 1444699 bytes
 .../derived/events/events_x+650mm.npz              |  Bin 0 -> 1398577 bytes
 .../derived/events/events_x+670mm.npz              |  Bin 0 -> 1376021 bytes
 .../derived/events/events_x+690mm.npz              |  Bin 0 -> 1364036 bytes
 .../derived/events/events_x-100mm.npz              |  Bin 0 -> 1688855 bytes
 .../derived/events/events_x-150mm.npz              |  Bin 0 -> 1685757 bytes
 .../derived/events/events_x-200mm.npz              |  Bin 0 -> 1674171 bytes
 .../derived/events/events_x-250mm.npz              |  Bin 0 -> 1659670 bytes
 .../derived/events/events_x-300mm.npz              |  Bin 0 -> 1643825 bytes
 .../derived/events/events_x-350mm.npz              |  Bin 0 -> 1619374 bytes
 .../derived/events/events_x-400mm.npz              |  Bin 0 -> 1594561 bytes
 .../derived/events/events_x-450mm.npz              |  Bin 0 -> 1560417 bytes
 .../derived/events/events_x-500mm.npz              |  Bin 0 -> 1526242 bytes
 .../derived/events/events_x-50mm.npz               |  Bin 0 -> 1693474 bytes
 .../derived/events/events_x-550mm.npz              |  Bin 0 -> 1487890 bytes
 .../derived/events/events_x-600mm.npz              |  Bin 0 -> 1439585 bytes
 .../derived/events/events_x-650mm.npz              |  Bin 0 -> 1395798 bytes
 .../derived/events/events_x-670mm.npz              |  Bin 0 -> 1377928 bytes
 .../derived/events/events_x-690mm.npz              |  Bin 0 -> 1363389 bytes
 .../exec12_20260612_191000/figures/bias_vs_x.pdf   |  Bin 0 -> 22417 bytes
 .../figures/calibration_end_ratio.pdf              |  Bin 0 -> 15355 bytes
 .../figures/calibration_end_timing.pdf             |  Bin 0 -> 15330 bytes
 .../figures/calibration_top_centroid.pdf           |  Bin 0 -> 15180 bytes
 .../figures/end_channel_symmetry.pdf               |  Bin 0 -> 13810 bytes
 .../residual_distributions_selected_positions.pdf  |  Bin 0 -> 18551 bytes
 .../exec12_20260612_191000/figures/rms68_vs_x.pdf  |  Bin 0 -> 24309 bytes
 .../figures/sigma_core_vs_x.pdf                    |  Bin 0 -> 24945 bytes
 .../figures/valid_fraction_vs_x.pdf                |  Bin 0 -> 19674 bytes
 .../figures/y_centroid_mean_vs_x.pdf               |  Bin 0 -> 16926 bytes
 .../figures/y_centroid_width_vs_x.pdf              |  Bin 0 -> 17309 bytes
 .../figures/y_left_vs_y_right.pdf                  |  Bin 0 -> 18062 bytes
 results/exec12_20260612_191000/logs/analysis.log   |    0
 results/exec12_20260612_191000/logs/inventory.log  |   31 +
 .../report/endtop_position_reconstruction.md       |   59 +
 .../report/endtop_position_reconstruction.tex      |   60 +
 .../report/exec11_technical_note.md                |   42 +
 .../report/exec11_technical_note.tex               |   43 +
 .../tables/x_reconstruction_summary.tex            |  190 +++
 results/exec12t_20260612_195426/README.md          |   48 +
 .../analysis/calibration_20pe.csv                  |    2 +
 .../analysis/calibration_4pe.csv                   |    2 +
 .../analysis/configuration_provenance.json         |   10 +
 .../analysis/data_inventory.csv                    |   42 +
 .../analysis/exec11_reproduction_check.csv         |   42 +
 .../analysis/loo_predictions_20pe.npz              |  Bin 0 -> 892931 bytes
 .../analysis/loo_predictions_4pe.npz               |  Bin 0 -> 892788 bytes
 .../analysis/position_reconstruction_20pe.csv      |   42 +
 .../analysis/position_reconstruction_4pe.csv       |   42 +
 .../analysis/temporal_position_summary.csv         |   83 ++
 .../analysis/threshold_4_20_summary.csv            |    5 +
 .../analysis/threshold_sweep_summary.csv           |   31 +
 .../analysis/window_dip_summary.csv                |    9 +
 .../beamer/beamer_contact_sheet.png                |  Bin 0 -> 473675 bytes
 .../beamer/exec12t_timing_position_beamer.aux      |  128 ++
 .../beamer/exec12t_timing_position_beamer.log      | 1406 ++++++++++++++++++++
 .../beamer/exec12t_timing_position_beamer.nav      |  109 ++
 .../beamer/exec12t_timing_position_beamer.out      |    1 +
 .../beamer/exec12t_timing_position_beamer.pdf      |  Bin 0 -> 234763 bytes
 .../beamer/exec12t_timing_position_beamer.snm      |    0
 .../beamer/exec12t_timing_position_beamer.tex      |  166 +++
 .../beamer/exec12t_timing_position_beamer.toc      |    0
 .../exec12t_20260612_195426/beamer/references.bib  |    1 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 29095 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 47777 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 19663 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 27031 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 15994 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 20580 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 32086 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 26819 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 26911 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 20004 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 20223 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 23610 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 24670 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 22094 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 30569 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 30458 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 18922 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 29235 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 15229 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 19449 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 20846 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 25617 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 24429 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 19255 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 19840 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 24193 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 24012 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 23457 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 31499 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 33489 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 19206 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 27919 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 15664 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 18758 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 21653 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 25116 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 46153 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 21146 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 25827 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 26188 bytes
 .../beamer/rendered/slide-41.png                   |  Bin 0 -> 18935 bytes
 .../beamer/rendered/slide-42.png                   |  Bin 0 -> 19827 bytes
 .../beamer/rendered/slide-43.png                   |  Bin 0 -> 25075 bytes
 .../beamer/rendered/slide-44.png                   |  Bin 0 -> 25549 bytes
 .../beamer/rendered/slide-45.png                   |  Bin 0 -> 22175 bytes
 .../beamer/rendered/slide-46.png                   |  Bin 0 -> 33040 bytes
 .../beamer/rendered/slide-47.png                   |  Bin 0 -> 32296 bytes
 .../beamer/rendered/slide-48.png                   |  Bin 0 -> 20505 bytes
 .../beamer/rendered/slide-49.png                   |  Bin 0 -> 28637 bytes
 .../beamer/speaker_notes.md                        |  107 ++
 .../derived/order_statistics/pair_x-422.0mm.npz    |  Bin 0 -> 1362830 bytes
 .../derived/order_statistics/pair_x-423.0mm.npz    |  Bin 0 -> 1362991 bytes
 .../derived/order_statistics/pair_x-424.0mm.npz    |  Bin 0 -> 1362994 bytes
 .../derived/order_statistics/pair_x-425.0mm.npz    |  Bin 0 -> 1362676 bytes
 .../derived/order_statistics/pair_x-426.0mm.npz    |  Bin 0 -> 1362617 bytes
 .../derived/order_statistics/pair_x-427.0mm.npz    |  Bin 0 -> 1362456 bytes
 .../derived/order_statistics/pair_x-428.0mm.npz    |  Bin 0 -> 1362417 bytes
 .../derived/order_statistics/pair_x-429.0mm.npz    |  Bin 0 -> 1362477 bytes
 .../derived/order_statistics/pair_x-430.0mm.npz    |  Bin 0 -> 1362246 bytes
 .../derived/order_statistics/pair_x-431.0mm.npz    |  Bin 0 -> 1361625 bytes
 .../derived/order_statistics/pair_x-432.0mm.npz    |  Bin 0 -> 1361227 bytes
 .../derived/order_statistics/pair_x-433.0mm.npz    |  Bin 0 -> 1361512 bytes
 .../derived/order_statistics/pair_x-434.0mm.npz    |  Bin 0 -> 1361634 bytes
 .../derived/order_statistics/pair_x-435.0mm.npz    |  Bin 0 -> 1361842 bytes
 .../derived/order_statistics/pair_x-436.0mm.npz    |  Bin 0 -> 1361904 bytes
 .../derived/order_statistics/pair_x-437.0mm.npz    |  Bin 0 -> 1362200 bytes
 .../derived/order_statistics/pair_x-438.0mm.npz    |  Bin 0 -> 1362269 bytes
 .../derived/order_statistics/pair_x-439.0mm.npz    |  Bin 0 -> 1362081 bytes
 .../derived/order_statistics/pair_x-440.0mm.npz    |  Bin 0 -> 1361731 bytes
 .../derived/order_statistics/pair_x-441.0mm.npz    |  Bin 0 -> 1361645 bytes
 .../derived/order_statistics/pair_x-442.0mm.npz    |  Bin 0 -> 1361574 bytes
 .../derived/order_statistics/pair_x-443.0mm.npz    |  Bin 0 -> 1361936 bytes
 .../derived/order_statistics/pair_x-444.0mm.npz    |  Bin 0 -> 1361808 bytes
 .../derived/order_statistics/pair_x-445.0mm.npz    |  Bin 0 -> 1362189 bytes
 .../derived/order_statistics/pair_x-446.0mm.npz    |  Bin 0 -> 1362266 bytes
 .../derived/order_statistics/pair_x-447.0mm.npz    |  Bin 0 -> 1362136 bytes
 .../derived/order_statistics/pair_x-448.0mm.npz    |  Bin 0 -> 1362063 bytes
 .../derived/order_statistics/pair_x-449.0mm.npz    |  Bin 0 -> 1362102 bytes
 .../derived/order_statistics/pair_x-450.0mm.npz    |  Bin 0 -> 1361698 bytes
 .../derived/order_statistics/pair_x-451.0mm.npz    |  Bin 0 -> 1361068 bytes
 .../derived/order_statistics/pair_x-452.0mm.npz    |  Bin 0 -> 1361575 bytes
 .../derived/order_statistics/pair_x-453.0mm.npz    |  Bin 0 -> 1361701 bytes
 .../derived/order_statistics/pair_x-454.0mm.npz    |  Bin 0 -> 1362028 bytes
 .../derived/order_statistics/pair_x-455.0mm.npz    |  Bin 0 -> 1362411 bytes
 .../derived/order_statistics/pair_x-456.0mm.npz    |  Bin 0 -> 1362338 bytes
 .../derived/order_statistics/pair_x-457.0mm.npz    |  Bin 0 -> 1362651 bytes
 .../derived/order_statistics/pair_x-458.0mm.npz    |  Bin 0 -> 1362862 bytes
 .../derived/order_statistics/pair_x-459.0mm.npz    |  Bin 0 -> 1362956 bytes
 .../derived/order_statistics/pair_x-460.0mm.npz    |  Bin 0 -> 1362783 bytes
 .../derived/order_statistics/pair_x-461.0mm.npz    |  Bin 0 -> 1362770 bytes
 .../derived/order_statistics/pair_x-462.0mm.npz    |  Bin 0 -> 1362960 bytes
 .../environment_exec12t.txt                        |   99 ++
 .../figures/correlation_ab_vs_x.pdf                |  Bin 0 -> 12608 bytes
 .../figures/efficiency_4_20_vs_x.pdf               |  Bin 0 -> 12461 bytes
 .../figures/mean_delta_4_20_vs_x.pdf               |  Bin 0 -> 13853 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 9467 bytes
 .../figures/pos_ref_1_delta_t_4_20.pdf             |  Bin 0 -> 15320 bytes
 .../figures/pos_ref_2_delta_t_4_20.pdf             |  Bin 0 -> 15348 bytes
 .../figures/threshold_bias.pdf                     |  Bin 0 -> 15585 bytes
 .../figures/threshold_chi2.pdf                     |  Bin 0 -> 15384 bytes
 .../figures/threshold_efficiency.pdf               |  Bin 0 -> 14268 bytes
 .../figures/threshold_pareto.pdf                   |  Bin 0 -> 14148 bytes
 .../figures/threshold_sigma_dt.pdf                 |  Bin 0 -> 15598 bytes
 .../figures/threshold_sigma_x.pdf                  |  Bin 0 -> 15710 bytes
 .../figures/threshold_slope.pdf                    |  Bin 0 -> 15487 bytes
 .../figures/timing_width_4_20_vs_x.pdf             |  Bin 0 -> 15411 bytes
 results/exec12t_20260612_195426/logs/analysis.log  |    0
 .../logs/beamer_pdfinfo.txt                        |   20 +
 .../logs/beamer_warnings.txt                       |    1 +
 .../logs/exec11_reproduction.log                   |    3 +
 results/exec12t_20260612_195426/logs/inventory.log |   41 +
 results/exec12t_20260612_195426/logs/window.log    |    0
 .../report/exec12t_timing_position_report.aux      |   26 +
 .../report/exec12t_timing_position_report.log      |  376 ++++++
 .../report/exec12t_timing_position_report.md       |   45 +
 .../report/exec12t_timing_position_report.out      |    7 +
 .../report/exec12t_timing_position_report.pdf      |  Bin 0 -> 35001 bytes
 .../report/exec12t_timing_position_report.tex      |   63 +
 .../tables/generated_numbers.tex                   |    6 +
 .../tables/global_context.csv                      |    7 +
 .../tables/global_context.tex                      |   10 +
 .../tables/threshold_4_20_summary.csv              |    7 +
 .../tables/threshold_4_20_summary.tex              |   10 +
 .../tables/threshold_sweep.csv                     |   31 +
 .../tables/threshold_sweep.tex                     |   34 +
 .../tables/window_dip_summary.csv                  |    9 +
 .../tables/window_dip_summary.tex                  |   12 +
 .../beamer/exec12tb_beamer.pdf                     |  Bin 0 -> 544006 bytes
 .../beamer/exec12tb_beamer.tex                     |  971 ++++++++++++++
 .../beamer/generated_numbers.tex                   |   24 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 37541 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 55903 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 39821 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 44919 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 37183 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 50484 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 53853 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 54725 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 55833 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 57120 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 55518 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 50060 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 46531 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 53957 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 45247 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 58873 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 56682 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 50931 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 57910 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 49913 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 63029 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 66176 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 49963 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 45261 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 56098 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 68087 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 72976 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 71762 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 45414 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 8355 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 54362 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 54841 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 40412 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 46539 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 61442 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 54635 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 47345 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 50838 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 47245 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 44446 bytes
 .../beamer/speaker_notes.md                        |  227 ++++
 .../figures/cfd_vs_orderstat_schematic.pdf         |  Bin 0 -> 32384 bytes
 .../figures/fine_scan_geometry.pdf                 |  Bin 0 -> 18315 bytes
 .../figures/global_context_sigma_bias.pdf          |  Bin 0 -> 21642 bytes
 .../figures/mean_dt_vs_x_with_residuals.pdf        |  Bin 0 -> 25199 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 21280 bytes
 .../figures/position_resolution_4_20_vs_x.pdf      |  Bin 0 -> 25708 bytes
 .../figures/ref1_dt_overlay.pdf                    |  Bin 0 -> 25821 bytes
 .../figures/ref1_xrec_overlay.pdf                  |  Bin 0 -> 26514 bytes
 .../figures/ref2_dt_overlay.pdf                    |  Bin 0 -> 25730 bytes
 .../figures/ref2_xrec_overlay.pdf                  |  Bin 0 -> 26561 bytes
 .../figures/rho_ab_vs_x_4_20.pdf                   |  Bin 0 -> 22331 bytes
 .../figures/rms68_dt_vs_x_4_20.pdf                 |  Bin 0 -> 21590 bytes
 .../figures/sha_manifest.json                      |   22 +
 .../figures/sweep_bias_vs_k.pdf                    |  Bin 0 -> 17260 bytes
 .../figures/sweep_chi2_vs_k.pdf                    |  Bin 0 -> 19883 bytes
 .../figures/sweep_rms68x_vs_k.pdf                  |  Bin 0 -> 20698 bytes
 .../figures/sweep_slope_vs_k.pdf                   |  Bin 0 -> 17743 bytes
 .../figures/tplus_rms68_vs_x_4_20.pdf              |  Bin 0 -> 22163 bytes
 .../figures/tradeoff_resolution_vs_bias.pdf        |  Bin 0 -> 27883 bytes
 .../figures/window_dip_counts.pdf                  |  Bin 0 -> 16254 bytes
 .../figures/window_dip_t4_t20.pdf                  |  Bin 0 -> 18771 bytes
 .../logs/beamer_compile1.log                       |  255 ++++
 .../logs/beamer_compile2.log                       |  265 ++++
 .../logs/beamer_compile3.log                       |  275 ++++
 .../logs/beamer_compile4.log                       |  618 +++++++++
 .../logs/beamer_compile5.log                       |  598 +++++++++
 results/exec12tb_20260612_204216/logs/visual_qa.md |   60 +
 .../manifest/beamer_manifest.csv                   |   41 +
 .../tables/dataset_inventory.tex                   |   15 +
 .../tables/global_context.tex                      |   18 +
 .../tables/lv_comparison.tex                       |   23 +
 .../tables/reference_positions.tex                 |   14 +
 .../tables/threshold_4_20_comparison.tex           |   22 +
 .../tables/threshold_sweep.tex                     |   22 +
 .../tables/window_dip_summary.tex                  |   16 +
 scripts/analyze_geometry.py                        |  172 +++
 scripts/build_exec12t.sh                           |   13 +
 scripts/gen_pair_scan_macros.py                    |  119 ++
 scripts/make_contact_sheet.py                      |   36 +
 scripts/run_pair_scan.sh                           |  156 +++
 scripts/run_resolution_scan.sh                     |  115 ++
 src/DetectorConstruction.cc                        |  123 +-
 tests/check_endtop_balance.py                      |   14 +-
 tests/check_pairscan_geometry.py                   |  123 ++
 tests/check_pairscan_macros.py                     |   81 ++
 tests/readout_config_check.cc                      |   53 +-
 tests/test_exec11_pair_analysis.py                 |   64 +
 tests/test_exec12_endtop_position.py               |   20 +
 tests/test_exec12t_timing_threshold_analysis.py    |   15 +
 452 files changed, 15402 insertions(+), 135 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E135](#e135)
<a id="b6"></a>

### 6. feat/ej204-bar-tir-only
Inspected ref `feat/ej204-bar-tir-only`, full SHA `09f8b1828889f50881639b6a44a137b4a8b27ef5`. Classification: **Alternative TIR configuration; scientific suitability UNVERIFIED**. [E032](#e032)
**B:** ahead 32, behind 115; merged NO; contains 8349041 NO. [E062](#e062) [E223](#e223) [E035](#e035)
Upstream: none configured. Origin counterpart: `09f8b1828889f50881639b6a44a137b4a8b27ef5`. [E030](#e030)
Tip date / author / subject: 2026-08-15 18:50:16 +0200 / rrios / feat(analysis): add bar timing resolution analysis and run macro. [E032](#e032)
**C1/C2:** active configuration: polished dielectric_dielectric TIR-only bar skin. Reflectivity: 0.98 factory retained, not used by bar skin. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E136](#e136) [E137](#e137) [E301](#e301)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **6**. Full paths, line numbers and contents are preserved in [E138](#e138) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E137](#e137)
```text
258: 
259:     // TIR-only: polished dielectric-dielectric surface on all bar faces.
260:     // Geant4 Fresnel equations give 100% TIR for theta > arcsin(1/1.58) = 39.3 deg.
262:     // No reflective coating — lateral faces are bare air.
263:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV,
264:                                               Materials::CreateBarSurface());
265:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E136](#e136) [E301](#e301)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E302](#e302). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E065](#e065)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E065](#e065)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E063](#e063)
```text
09f8b18 feat(analysis): add bar timing resolution analysis and run macro
234e0cb feat(ej204-bar-tir-only): EJ-204 bar TIR-only — remove Mylar skin reflector
f19c093 feat(pairscan): CMakeLists refactor, pair scan macros and run script
14ae395 fix(tests): replace escape-fraction guard with sipm-entry guard
c7a627e fix(optics): replace reflector volumes with bar skin surface
021d5f4 EXEC_12TB: add rebuilt self-contained beamer and speaker notes
1ac830a EXEC_12TB: add beamer manifest and Makefile exec12tb targets
d764892 EXEC_12TB: add generated tables and TeX number macros
9cfe042 EXEC_12TB: add figure style module and regenerated deck figures
66ac24c EXEC_12T: add complete Beamer presentation and reproducible build
b0aaac1 EXEC_12T: add integrated timing-position technical report
a2563d9 EXEC_12T: add window-dip timing and collection study
aaabe89 EXEC_12T: add threshold sweep and covariance timing decomposition
0a41034 EXEC_12T: add 4PE versus 20PE temporal and position analysis
f38c9e1 EXEC_12T: add cached fourth-to-thirtieth hit order statistics
003d106 EXEC_12: add staged XY scan proposal
6e01606 EXEC_12: add EXEC11 and EndTop technical reports
30166cd EXEC_12: add y-zero observability feasibility study
a9571fd EXEC_12: add covariance-aware EndTop estimator combination
bd21563 EXEC_12: add leave-one-position-out global X reconstruction
8332733 EXEC_12: add leave-one-position-out global X reconstruction
17d9e7d EXEC_12: add EndTop event observable builder and geometry tests
4dba735 EXEC_11: add reproducible report tables and README
a103ed6 EXEC_11: add temporal ratio and covariance-aware position reconstruction
1e95ac6 EXEC_11: add detailed two-position timing analysis
61bca2e EXEC_11: regenerate pair-scan fit QA and summary v2
bf0d2d5 EXEC_11: add per-event pair observable builder and tests
f431c01 EXEC_08/Step7: CTest guardrails for pair-scan geometry and macros
20bed82 EXEC_08/Step6: batch runner scripts/run_pair_scan.sh
549f3a6 EXEC_08/Step5: ROOT analysis macro analyze_pair_scan.C
30f3da3 EXEC_08/Step3: generate 41 pair-scan macros in macros/pairscan/
6bb73d5 EXEC_08/Step2: geometry analysis + pair selection
```

**Files touched on the branch since its merge base with main:** [E064](#e064)
```text
 .gitignore                                         |    6 +
 CMakeLists.txt                                     |   76 +-
 Makefile                                           |   40 +
 analysis/analyze_pair_scan.C                       |  496 +++++++
 analysis/bar_timing_resolution.cxx                 |  387 ++++++
 analysis/exec11_pair_analysis.py                   |  934 +++++++++++++
 analysis/exec12_endtop_position.py                 |  161 +++
 analysis/exec12_make_reports.py                    |   99 ++
 analysis/exec12t_make_products.py                  |  151 +++
 analysis/exec12t_timing_threshold_analysis.py      |  136 ++
 analysis/exec12tb_figstyle.py                      |   56 +
 analysis/exec12tb_figures.py                       |  823 ++++++++++++
 analysis/exec12tb_tables.py                        |  487 +++++++
 include/DetectorConstruction.hh                    |    8 +-
 macros/pairscan/pairscan_x-422.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-423.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-424.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-425.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-426.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-427.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-428.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-429.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-430.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-431.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-432.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-433.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-434.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-435.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-436.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-437.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-438.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-439.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-440.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-441.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-442.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-443.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-444.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-445.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-446.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-447.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-448.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-449.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-450.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-451.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-452.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-453.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-454.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-455.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-456.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-457.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-458.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-459.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-460.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-461.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-462.0mm.mac             |   20 +
 macros/timing_tir_center.mac                       |   27 +
 pair_scan_config.json                              |   71 +
 report/exec13_xy_scan_plan.md                      |   30 +
 results/exec11_20260612_182454/README.md           |  124 ++
 .../analysis/calibration_ratio.csv                 |    3 +
 .../analysis/calibration_temporal.csv              |    2 +
 .../analysis/calibration_v1_v2.csv                 |    3 +
 .../analysis/data_inventory.csv                    |   42 +
 results/exec11_20260612_182454/analysis/fit_qa.csv |   42 +
 .../analysis/fit_qa_v1_v2_focus.csv                |    4 +
 .../exec11_20260612_182454/analysis/metadata.json  |   13 +
 .../analysis/pairscan_summary_v2.csv               |   42 +
 .../analysis/reconstruction_summary.csv            |   11 +
 .../analysis/reference_comparison.csv              |    3 +
 .../analysis/reference_positions.csv               |    3 +
 .../derived/pair_events_x-422.0mm.npz              |  Bin 0 -> 123534 bytes
 .../derived/pair_events_x-423.0mm.npz              |  Bin 0 -> 122740 bytes
 .../derived/pair_events_x-424.0mm.npz              |  Bin 0 -> 123199 bytes
 .../derived/pair_events_x-425.0mm.npz              |  Bin 0 -> 123147 bytes
 .../derived/pair_events_x-426.0mm.npz              |  Bin 0 -> 123368 bytes
 .../derived/pair_events_x-427.0mm.npz              |  Bin 0 -> 122917 bytes
 .../derived/pair_events_x-428.0mm.npz              |  Bin 0 -> 122839 bytes
 .../derived/pair_events_x-429.0mm.npz              |  Bin 0 -> 122318 bytes
 .../derived/pair_events_x-430.0mm.npz              |  Bin 0 -> 121727 bytes
 .../derived/pair_events_x-431.0mm.npz              |  Bin 0 -> 120954 bytes
 .../derived/pair_events_x-432.0mm.npz              |  Bin 0 -> 120412 bytes
 .../derived/pair_events_x-433.0mm.npz              |  Bin 0 -> 121077 bytes
 .../derived/pair_events_x-434.0mm.npz              |  Bin 0 -> 122386 bytes
 .../derived/pair_events_x-435.0mm.npz              |  Bin 0 -> 122900 bytes
 .../derived/pair_events_x-436.0mm.npz              |  Bin 0 -> 122357 bytes
 .../derived/pair_events_x-437.0mm.npz              |  Bin 0 -> 122009 bytes
 .../derived/pair_events_x-438.0mm.npz              |  Bin 0 -> 122595 bytes
 .../derived/pair_events_x-439.0mm.npz              |  Bin 0 -> 122198 bytes
 .../derived/pair_events_x-440.0mm.npz              |  Bin 0 -> 122232 bytes
 .../derived/pair_events_x-441.0mm.npz              |  Bin 0 -> 122736 bytes
 .../derived/pair_events_x-442.0mm.npz              |  Bin 0 -> 122672 bytes
 .../derived/pair_events_x-443.0mm.npz              |  Bin 0 -> 122553 bytes
 .../derived/pair_events_x-444.0mm.npz              |  Bin 0 -> 122542 bytes
 .../derived/pair_events_x-445.0mm.npz              |  Bin 0 -> 122293 bytes
 .../derived/pair_events_x-446.0mm.npz              |  Bin 0 -> 122306 bytes
 .../derived/pair_events_x-447.0mm.npz              |  Bin 0 -> 121906 bytes
 .../derived/pair_events_x-448.0mm.npz              |  Bin 0 -> 122022 bytes
 .../derived/pair_events_x-449.0mm.npz              |  Bin 0 -> 122871 bytes
 .../derived/pair_events_x-450.0mm.npz              |  Bin 0 -> 122284 bytes
 .../derived/pair_events_x-451.0mm.npz              |  Bin 0 -> 120956 bytes
 .../derived/pair_events_x-452.0mm.npz              |  Bin 0 -> 120882 bytes
 .../derived/pair_events_x-453.0mm.npz              |  Bin 0 -> 120856 bytes
 .../derived/pair_events_x-454.0mm.npz              |  Bin 0 -> 121828 bytes
 .../derived/pair_events_x-455.0mm.npz              |  Bin 0 -> 122035 bytes
 .../derived/pair_events_x-456.0mm.npz              |  Bin 0 -> 122679 bytes
 .../derived/pair_events_x-457.0mm.npz              |  Bin 0 -> 123063 bytes
 .../derived/pair_events_x-458.0mm.npz              |  Bin 0 -> 123181 bytes
 .../derived/pair_events_x-459.0mm.npz              |  Bin 0 -> 123238 bytes
 .../derived/pair_events_x-460.0mm.npz              |  Bin 0 -> 123352 bytes
 .../derived/pair_events_x-461.0mm.npz              |  Bin 0 -> 123315 bytes
 .../derived/pair_events_x-462.0mm.npz              |  Bin 0 -> 123329 bytes
 .../figures/calibrations_and_residuals.pdf         |  Bin 0 -> 25616 bytes
 .../figures/calibrations_and_residuals.png         |  Bin 0 -> 284318 bytes
 .../figures/pos_ref_1_all_hit_times.pdf            |  Bin 0 -> 23571 bytes
 .../figures/pos_ref_1_all_hit_times.png            |  Bin 0 -> 146556 bytes
 .../figures/pos_ref_1_correlations.pdf             |  Bin 0 -> 25033 bytes
 .../figures/pos_ref_1_correlations.png             |  Bin 0 -> 298577 bytes
 .../figures/pos_ref_1_delta_t.pdf                  |  Bin 0 -> 20414 bytes
 .../figures/pos_ref_1_delta_t.png                  |  Bin 0 -> 113077 bytes
 .../figures/pos_ref_1_event_times.pdf              |  Bin 0 -> 20672 bytes
 .../figures/pos_ref_1_event_times.png              |  Bin 0 -> 89908 bytes
 .../figures/pos_ref_1_npe_moyal.pdf                |  Bin 0 -> 26752 bytes
 .../figures/pos_ref_1_npe_moyal.png                |  Bin 0 -> 180911 bytes
 .../figures/pos_ref_2_all_hit_times.pdf            |  Bin 0 -> 23694 bytes
 .../figures/pos_ref_2_all_hit_times.png            |  Bin 0 -> 149031 bytes
 .../figures/pos_ref_2_correlations.pdf             |  Bin 0 -> 24595 bytes
 .../figures/pos_ref_2_correlations.png             |  Bin 0 -> 275444 bytes
 .../figures/pos_ref_2_delta_t.pdf                  |  Bin 0 -> 19869 bytes
 .../figures/pos_ref_2_delta_t.png                  |  Bin 0 -> 114187 bytes
 .../figures/pos_ref_2_event_times.pdf              |  Bin 0 -> 20491 bytes
 .../figures/pos_ref_2_event_times.png              |  Bin 0 -> 88939 bytes
 .../figures/pos_ref_2_npe_moyal.pdf                |  Bin 0 -> 26630 bytes
 .../figures/pos_ref_2_npe_moyal.png                |  Bin 0 -> 177596 bytes
 .../figures/position_reconstruction.pdf            |  Bin 0 -> 23751 bytes
 .../figures/position_reconstruction.png            |  Bin 0 -> 180735 bytes
 results/exec11_20260612_182454/logs/derive.log     |   41 +
 results/exec11_20260612_182454/logs/detail.log     |    0
 results/exec11_20260612_182454/logs/qa.log         |    0
 .../exec11_20260612_182454/logs/reconstruct.log    |    0
 results/exec11_20260612_182454/logs/report.log     |    0
 .../exec11_20260612_182454/tables/fit_qa_v1_v2.tex |    9 +
 .../tables/reconstruction_summary.tex              |   16 +
 .../tables/reference_comparison.tex                |    8 +
 results/exec12_20260612_191000/README.md           |   62 +
 .../analysis/blue_summary.csv                      |   32 +
 .../analysis/blue_weights.csv                      |   32 +
 .../analysis/configuration_provenance.json         |   29 +
 .../analysis/cv_calibrations.csv                   |  187 +++
 .../analysis/cv_model_selection.csv                |    2 +
 .../analysis/cv_predictions.csv.gz                 |  Bin 0 -> 4316374 bytes
 .../analysis/data_inventory.csv                    |   32 +
 .../analysis/x_reconstruction_summary.csv          |  187 +++
 .../analysis/y0_feasibility_summary.csv            |   94 ++
 .../derived/events/events_x+0mm.npz                |  Bin 0 -> 1673253 bytes
 .../derived/events/events_x+100mm.npz              |  Bin 0 -> 1690764 bytes
 .../derived/events/events_x+150mm.npz              |  Bin 0 -> 1686376 bytes
 .../derived/events/events_x+200mm.npz              |  Bin 0 -> 1675705 bytes
 .../derived/events/events_x+250mm.npz              |  Bin 0 -> 1660646 bytes
 .../derived/events/events_x+300mm.npz              |  Bin 0 -> 1641542 bytes
 .../derived/events/events_x+350mm.npz              |  Bin 0 -> 1619395 bytes
 .../derived/events/events_x+400mm.npz              |  Bin 0 -> 1596329 bytes
 .../derived/events/events_x+450mm.npz              |  Bin 0 -> 1564409 bytes
 .../derived/events/events_x+500mm.npz              |  Bin 0 -> 1525987 bytes
 .../derived/events/events_x+50mm.npz               |  Bin 0 -> 1693544 bytes
 .../derived/events/events_x+550mm.npz              |  Bin 0 -> 1488016 bytes
 .../derived/events/events_x+600mm.npz              |  Bin 0 -> 1444699 bytes
 .../derived/events/events_x+650mm.npz              |  Bin 0 -> 1398577 bytes
 .../derived/events/events_x+670mm.npz              |  Bin 0 -> 1376021 bytes
 .../derived/events/events_x+690mm.npz              |  Bin 0 -> 1364036 bytes
 .../derived/events/events_x-100mm.npz              |  Bin 0 -> 1688855 bytes
 .../derived/events/events_x-150mm.npz              |  Bin 0 -> 1685757 bytes
 .../derived/events/events_x-200mm.npz              |  Bin 0 -> 1674171 bytes
 .../derived/events/events_x-250mm.npz              |  Bin 0 -> 1659670 bytes
 .../derived/events/events_x-300mm.npz              |  Bin 0 -> 1643825 bytes
 .../derived/events/events_x-350mm.npz              |  Bin 0 -> 1619374 bytes
 .../derived/events/events_x-400mm.npz              |  Bin 0 -> 1594561 bytes
 .../derived/events/events_x-450mm.npz              |  Bin 0 -> 1560417 bytes
 .../derived/events/events_x-500mm.npz              |  Bin 0 -> 1526242 bytes
 .../derived/events/events_x-50mm.npz               |  Bin 0 -> 1693474 bytes
 .../derived/events/events_x-550mm.npz              |  Bin 0 -> 1487890 bytes
 .../derived/events/events_x-600mm.npz              |  Bin 0 -> 1439585 bytes
 .../derived/events/events_x-650mm.npz              |  Bin 0 -> 1395798 bytes
 .../derived/events/events_x-670mm.npz              |  Bin 0 -> 1377928 bytes
 .../derived/events/events_x-690mm.npz              |  Bin 0 -> 1363389 bytes
 .../exec12_20260612_191000/figures/bias_vs_x.pdf   |  Bin 0 -> 22417 bytes
 .../figures/calibration_end_ratio.pdf              |  Bin 0 -> 15355 bytes
 .../figures/calibration_end_timing.pdf             |  Bin 0 -> 15330 bytes
 .../figures/calibration_top_centroid.pdf           |  Bin 0 -> 15180 bytes
 .../figures/end_channel_symmetry.pdf               |  Bin 0 -> 13810 bytes
 .../residual_distributions_selected_positions.pdf  |  Bin 0 -> 18551 bytes
 .../exec12_20260612_191000/figures/rms68_vs_x.pdf  |  Bin 0 -> 24309 bytes
 .../figures/sigma_core_vs_x.pdf                    |  Bin 0 -> 24945 bytes
 .../figures/valid_fraction_vs_x.pdf                |  Bin 0 -> 19674 bytes
 .../figures/y_centroid_mean_vs_x.pdf               |  Bin 0 -> 16926 bytes
 .../figures/y_centroid_width_vs_x.pdf              |  Bin 0 -> 17309 bytes
 .../figures/y_left_vs_y_right.pdf                  |  Bin 0 -> 18062 bytes
 results/exec12_20260612_191000/logs/analysis.log   |    0
 results/exec12_20260612_191000/logs/inventory.log  |   31 +
 .../report/endtop_position_reconstruction.md       |   59 +
 .../report/endtop_position_reconstruction.tex      |   60 +
 .../report/exec11_technical_note.md                |   42 +
 .../report/exec11_technical_note.tex               |   43 +
 .../tables/x_reconstruction_summary.tex            |  190 +++
 results/exec12t_20260612_195426/README.md          |   48 +
 .../analysis/calibration_20pe.csv                  |    2 +
 .../analysis/calibration_4pe.csv                   |    2 +
 .../analysis/configuration_provenance.json         |   10 +
 .../analysis/data_inventory.csv                    |   42 +
 .../analysis/exec11_reproduction_check.csv         |   42 +
 .../analysis/loo_predictions_20pe.npz              |  Bin 0 -> 892931 bytes
 .../analysis/loo_predictions_4pe.npz               |  Bin 0 -> 892788 bytes
 .../analysis/position_reconstruction_20pe.csv      |   42 +
 .../analysis/position_reconstruction_4pe.csv       |   42 +
 .../analysis/temporal_position_summary.csv         |   83 ++
 .../analysis/threshold_4_20_summary.csv            |    5 +
 .../analysis/threshold_sweep_summary.csv           |   31 +
 .../analysis/window_dip_summary.csv                |    9 +
 .../beamer/beamer_contact_sheet.png                |  Bin 0 -> 473675 bytes
 .../beamer/exec12t_timing_position_beamer.aux      |  128 ++
 .../beamer/exec12t_timing_position_beamer.log      | 1406 ++++++++++++++++++++
 .../beamer/exec12t_timing_position_beamer.nav      |  109 ++
 .../beamer/exec12t_timing_position_beamer.out      |    1 +
 .../beamer/exec12t_timing_position_beamer.pdf      |  Bin 0 -> 234763 bytes
 .../beamer/exec12t_timing_position_beamer.snm      |    0
 .../beamer/exec12t_timing_position_beamer.tex      |  166 +++
 .../beamer/exec12t_timing_position_beamer.toc      |    0
 .../exec12t_20260612_195426/beamer/references.bib  |    1 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 29095 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 47777 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 19663 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 27031 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 15994 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 20580 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 32086 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 26819 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 26911 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 20004 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 20223 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 23610 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 24670 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 22094 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 30569 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 30458 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 18922 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 29235 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 15229 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 19449 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 20846 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 25617 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 24429 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 19255 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 19840 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 24193 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 24012 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 23457 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 31499 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 33489 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 19206 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 27919 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 15664 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 18758 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 21653 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 25116 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 46153 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 21146 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 25827 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 26188 bytes
 .../beamer/rendered/slide-41.png                   |  Bin 0 -> 18935 bytes
 .../beamer/rendered/slide-42.png                   |  Bin 0 -> 19827 bytes
 .../beamer/rendered/slide-43.png                   |  Bin 0 -> 25075 bytes
 .../beamer/rendered/slide-44.png                   |  Bin 0 -> 25549 bytes
 .../beamer/rendered/slide-45.png                   |  Bin 0 -> 22175 bytes
 .../beamer/rendered/slide-46.png                   |  Bin 0 -> 33040 bytes
 .../beamer/rendered/slide-47.png                   |  Bin 0 -> 32296 bytes
 .../beamer/rendered/slide-48.png                   |  Bin 0 -> 20505 bytes
 .../beamer/rendered/slide-49.png                   |  Bin 0 -> 28637 bytes
 .../beamer/speaker_notes.md                        |  107 ++
 .../derived/order_statistics/pair_x-422.0mm.npz    |  Bin 0 -> 1362830 bytes
 .../derived/order_statistics/pair_x-423.0mm.npz    |  Bin 0 -> 1362991 bytes
 .../derived/order_statistics/pair_x-424.0mm.npz    |  Bin 0 -> 1362994 bytes
 .../derived/order_statistics/pair_x-425.0mm.npz    |  Bin 0 -> 1362676 bytes
 .../derived/order_statistics/pair_x-426.0mm.npz    |  Bin 0 -> 1362617 bytes
 .../derived/order_statistics/pair_x-427.0mm.npz    |  Bin 0 -> 1362456 bytes
 .../derived/order_statistics/pair_x-428.0mm.npz    |  Bin 0 -> 1362417 bytes
 .../derived/order_statistics/pair_x-429.0mm.npz    |  Bin 0 -> 1362477 bytes
 .../derived/order_statistics/pair_x-430.0mm.npz    |  Bin 0 -> 1362246 bytes
 .../derived/order_statistics/pair_x-431.0mm.npz    |  Bin 0 -> 1361625 bytes
 .../derived/order_statistics/pair_x-432.0mm.npz    |  Bin 0 -> 1361227 bytes
 .../derived/order_statistics/pair_x-433.0mm.npz    |  Bin 0 -> 1361512 bytes
 .../derived/order_statistics/pair_x-434.0mm.npz    |  Bin 0 -> 1361634 bytes
 .../derived/order_statistics/pair_x-435.0mm.npz    |  Bin 0 -> 1361842 bytes
 .../derived/order_statistics/pair_x-436.0mm.npz    |  Bin 0 -> 1361904 bytes
 .../derived/order_statistics/pair_x-437.0mm.npz    |  Bin 0 -> 1362200 bytes
 .../derived/order_statistics/pair_x-438.0mm.npz    |  Bin 0 -> 1362269 bytes
 .../derived/order_statistics/pair_x-439.0mm.npz    |  Bin 0 -> 1362081 bytes
 .../derived/order_statistics/pair_x-440.0mm.npz    |  Bin 0 -> 1361731 bytes
 .../derived/order_statistics/pair_x-441.0mm.npz    |  Bin 0 -> 1361645 bytes
 .../derived/order_statistics/pair_x-442.0mm.npz    |  Bin 0 -> 1361574 bytes
 .../derived/order_statistics/pair_x-443.0mm.npz    |  Bin 0 -> 1361936 bytes
 .../derived/order_statistics/pair_x-444.0mm.npz    |  Bin 0 -> 1361808 bytes
 .../derived/order_statistics/pair_x-445.0mm.npz    |  Bin 0 -> 1362189 bytes
 .../derived/order_statistics/pair_x-446.0mm.npz    |  Bin 0 -> 1362266 bytes
 .../derived/order_statistics/pair_x-447.0mm.npz    |  Bin 0 -> 1362136 bytes
 .../derived/order_statistics/pair_x-448.0mm.npz    |  Bin 0 -> 1362063 bytes
 .../derived/order_statistics/pair_x-449.0mm.npz    |  Bin 0 -> 1362102 bytes
 .../derived/order_statistics/pair_x-450.0mm.npz    |  Bin 0 -> 1361698 bytes
 .../derived/order_statistics/pair_x-451.0mm.npz    |  Bin 0 -> 1361068 bytes
 .../derived/order_statistics/pair_x-452.0mm.npz    |  Bin 0 -> 1361575 bytes
 .../derived/order_statistics/pair_x-453.0mm.npz    |  Bin 0 -> 1361701 bytes
 .../derived/order_statistics/pair_x-454.0mm.npz    |  Bin 0 -> 1362028 bytes
 .../derived/order_statistics/pair_x-455.0mm.npz    |  Bin 0 -> 1362411 bytes
 .../derived/order_statistics/pair_x-456.0mm.npz    |  Bin 0 -> 1362338 bytes
 .../derived/order_statistics/pair_x-457.0mm.npz    |  Bin 0 -> 1362651 bytes
 .../derived/order_statistics/pair_x-458.0mm.npz    |  Bin 0 -> 1362862 bytes
 .../derived/order_statistics/pair_x-459.0mm.npz    |  Bin 0 -> 1362956 bytes
 .../derived/order_statistics/pair_x-460.0mm.npz    |  Bin 0 -> 1362783 bytes
 .../derived/order_statistics/pair_x-461.0mm.npz    |  Bin 0 -> 1362770 bytes
 .../derived/order_statistics/pair_x-462.0mm.npz    |  Bin 0 -> 1362960 bytes
 .../environment_exec12t.txt                        |   99 ++
 .../figures/correlation_ab_vs_x.pdf                |  Bin 0 -> 12608 bytes
 .../figures/efficiency_4_20_vs_x.pdf               |  Bin 0 -> 12461 bytes
 .../figures/mean_delta_4_20_vs_x.pdf               |  Bin 0 -> 13853 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 9467 bytes
 .../figures/pos_ref_1_delta_t_4_20.pdf             |  Bin 0 -> 15320 bytes
 .../figures/pos_ref_2_delta_t_4_20.pdf             |  Bin 0 -> 15348 bytes
 .../figures/threshold_bias.pdf                     |  Bin 0 -> 15585 bytes
 .../figures/threshold_chi2.pdf                     |  Bin 0 -> 15384 bytes
 .../figures/threshold_efficiency.pdf               |  Bin 0 -> 14268 bytes
 .../figures/threshold_pareto.pdf                   |  Bin 0 -> 14148 bytes
 .../figures/threshold_sigma_dt.pdf                 |  Bin 0 -> 15598 bytes
 .../figures/threshold_sigma_x.pdf                  |  Bin 0 -> 15710 bytes
 .../figures/threshold_slope.pdf                    |  Bin 0 -> 15487 bytes
 .../figures/timing_width_4_20_vs_x.pdf             |  Bin 0 -> 15411 bytes
 results/exec12t_20260612_195426/logs/analysis.log  |    0
 .../logs/beamer_pdfinfo.txt                        |   20 +
 .../logs/beamer_warnings.txt                       |    1 +
 .../logs/exec11_reproduction.log                   |    3 +
 results/exec12t_20260612_195426/logs/inventory.log |   41 +
 results/exec12t_20260612_195426/logs/window.log    |    0
 .../report/exec12t_timing_position_report.aux      |   26 +
 .../report/exec12t_timing_position_report.log      |  376 ++++++
 .../report/exec12t_timing_position_report.md       |   45 +
 .../report/exec12t_timing_position_report.out      |    7 +
 .../report/exec12t_timing_position_report.pdf      |  Bin 0 -> 35001 bytes
 .../report/exec12t_timing_position_report.tex      |   63 +
 .../tables/generated_numbers.tex                   |    6 +
 .../tables/global_context.csv                      |    7 +
 .../tables/global_context.tex                      |   10 +
 .../tables/threshold_4_20_summary.csv              |    7 +
 .../tables/threshold_4_20_summary.tex              |   10 +
 .../tables/threshold_sweep.csv                     |   31 +
 .../tables/threshold_sweep.tex                     |   34 +
 .../tables/window_dip_summary.csv                  |    9 +
 .../tables/window_dip_summary.tex                  |   12 +
 .../beamer/exec12tb_beamer.pdf                     |  Bin 0 -> 544006 bytes
 .../beamer/exec12tb_beamer.tex                     |  971 ++++++++++++++
 .../beamer/generated_numbers.tex                   |   24 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 37541 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 55903 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 39821 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 44919 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 37183 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 50484 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 53853 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 54725 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 55833 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 57120 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 55518 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 50060 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 46531 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 53957 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 45247 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 58873 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 56682 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 50931 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 57910 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 49913 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 63029 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 66176 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 49963 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 45261 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 56098 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 68087 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 72976 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 71762 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 45414 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 8355 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 54362 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 54841 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 40412 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 46539 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 61442 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 54635 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 47345 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 50838 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 47245 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 44446 bytes
 .../beamer/speaker_notes.md                        |  227 ++++
 .../figures/cfd_vs_orderstat_schematic.pdf         |  Bin 0 -> 32384 bytes
 .../figures/fine_scan_geometry.pdf                 |  Bin 0 -> 18315 bytes
 .../figures/global_context_sigma_bias.pdf          |  Bin 0 -> 21642 bytes
 .../figures/mean_dt_vs_x_with_residuals.pdf        |  Bin 0 -> 25199 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 21280 bytes
 .../figures/position_resolution_4_20_vs_x.pdf      |  Bin 0 -> 25708 bytes
 .../figures/ref1_dt_overlay.pdf                    |  Bin 0 -> 25821 bytes
 .../figures/ref1_xrec_overlay.pdf                  |  Bin 0 -> 26514 bytes
 .../figures/ref2_dt_overlay.pdf                    |  Bin 0 -> 25730 bytes
 .../figures/ref2_xrec_overlay.pdf                  |  Bin 0 -> 26561 bytes
 .../figures/rho_ab_vs_x_4_20.pdf                   |  Bin 0 -> 22331 bytes
 .../figures/rms68_dt_vs_x_4_20.pdf                 |  Bin 0 -> 21590 bytes
 .../figures/sha_manifest.json                      |   22 +
 .../figures/sweep_bias_vs_k.pdf                    |  Bin 0 -> 17260 bytes
 .../figures/sweep_chi2_vs_k.pdf                    |  Bin 0 -> 19883 bytes
 .../figures/sweep_rms68x_vs_k.pdf                  |  Bin 0 -> 20698 bytes
 .../figures/sweep_slope_vs_k.pdf                   |  Bin 0 -> 17743 bytes
 .../figures/tplus_rms68_vs_x_4_20.pdf              |  Bin 0 -> 22163 bytes
 .../figures/tradeoff_resolution_vs_bias.pdf        |  Bin 0 -> 27883 bytes
 .../figures/window_dip_counts.pdf                  |  Bin 0 -> 16254 bytes
 .../figures/window_dip_t4_t20.pdf                  |  Bin 0 -> 18771 bytes
 .../logs/beamer_compile1.log                       |  255 ++++
 .../logs/beamer_compile2.log                       |  265 ++++
 .../logs/beamer_compile3.log                       |  275 ++++
 .../logs/beamer_compile4.log                       |  618 +++++++++
 .../logs/beamer_compile5.log                       |  598 +++++++++
 results/exec12tb_20260612_204216/logs/visual_qa.md |   60 +
 .../manifest/beamer_manifest.csv                   |   41 +
 .../tables/dataset_inventory.tex                   |   15 +
 .../tables/global_context.tex                      |   18 +
 .../tables/lv_comparison.tex                       |   23 +
 .../tables/reference_positions.tex                 |   14 +
 .../tables/threshold_4_20_comparison.tex           |   22 +
 .../tables/threshold_sweep.tex                     |   22 +
 .../tables/window_dip_summary.tex                  |   16 +
 scripts/analyze_geometry.py                        |  172 +++
 scripts/build_exec12t.sh                           |   13 +
 scripts/gen_pair_scan_macros.py                    |  119 ++
 scripts/make_contact_sheet.py                      |   36 +
 scripts/run_pair_scan.sh                           |  156 +++
 src/DetectorConstruction.cc                        |   87 +-
 tests/check_endtop_balance.py                      |   14 +-
 tests/check_pairscan_geometry.py                   |  123 ++
 tests/check_pairscan_macros.py                     |   81 ++
 tests/readout_config_check.cc                      |   53 +-
 tests/test_exec11_pair_analysis.py                 |   64 +
 tests/test_exec12_endtop_position.py               |   20 +
 tests/test_exec12t_timing_threshold_analysis.py    |   15 +
 446 files changed, 13647 insertions(+), 131 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E139](#e139)
<a id="b7"></a>

### 7. feat/ej204-event-display-tracks
Inspected ref `origin/feat/ej204-event-display-tracks`, full SHA `47a9a4fb31118ecb6f68cb658d40ec282c6da99f`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 4, behind 110; merged NO; contains 8349041 NO. [E140](#e140) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-14 16:25:57 +0200 / dowiyogo / fix(tests): update check_endtop_balance smoke test. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal polished bar skin. Reflectivity: 0.98 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E145](#e145) [E146](#e146) [E314](#e314)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **11**. Full paths, line numbers and contents are preserved in [E147](#e147) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E146](#e146)
```text
254:     // Apply branch-specific reflector properties directly to the bar.
255:     auto* reflector = Materials::CreateBarSkinReflector();
256:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
257:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E145](#e145) [E314](#e314)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E315](#e315). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E143](#e143)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E143](#e143)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E141](#e141)
```text
47a9a4f fix(tests): update check_endtop_balance smoke test
0a42656 fix(tests): replace escape-fraction guard with sipm-entry guard
b2718c4 fix(optics): replace reflector volumes with bar skin surface
c93502b Add selective optical trajectory storage for EJ-204 event display
```

**Files touched on the branch since its merge base with main:** [E142](#e142)
```text
 include/DetectorConstruction.hh     |   8 +-
 include/DisplayTrackingAction.hh    |  35 +++++
 macros/ej204_x0_first100_tracks.mac |  34 +++++
 src/ActionInitialization.cc         |   2 +
 src/DetectorConstruction.cc         |  73 +---------
 src/DisplayTrackingAction.cc        | 275 ++++++++++++++++++++++++++++++++++++
 src/EventAction.cc                  |  17 ++-
 src/SteppingAction.cc               |  34 +++++
 tests/check_endtop_balance.py       |   7 +-
 tests/readout_config_check.cc       |  53 +++----
 10 files changed, 425 insertions(+), 113 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E148](#e148)
<a id="b8"></a>

### 8. feat/ej228-cylinder
Inspected ref `origin/feat/ej228-cylinder`, full SHA `66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e`. Classification: **Alternative cylindrical geometry; suitability UNVERIFIED**. [E032](#e032)
**B:** ahead 3, behind 103; merged NO; contains 8349041 NO. [E149](#e149) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-15 00:03:35 +0200 / rrios / feat(ej228): EJ-228 cylinder simulation — 25mm diam × 25mm height. [E032](#e032)
**C1/C2:** active configuration: Cylinder air-gap + dielectric_metal polished outer border. Reflectivity: 0.98 active; 0.95 unused bar factory retained. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E154](#e154) [E155](#e155) [E316](#e316)

```text
211: G4OpticalSurface* CreateVikuitiSurface() {
216:     surf->SetType(dielectric_metal);
218:     surf->SetFinish(polished);
222:     const std::vector<G4double> refl   = {0.98, 0.98};
226:     mpt->AddProperty("REFLECTIVITY", energy, refl);
251: G4OpticalSurface* CreateBarSurface() {
261:     surf->SetType(dielectric_dielectric);
263:     surf->SetFinish(polished);
269: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
273:     surf->SetType(dielectric_metal);
275:     surf->SetFinish(polished);
294:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
328: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
335:     surf->SetType(dielectric_metal);
337:     surf->SetFinish(polished);
341:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
346:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
371: G4OpticalSurface* CreateBarSkinReflector() {
394:     surf->SetType(dielectric_dielectric);
396:     surf->SetFinish(polished);
401:     const std::vector<G4double> refl   = {0.95, 0.95};
404:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **16**. Full paths, line numbers and contents are preserved in [E156](#e156) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E155](#e155)
```text
95:     new G4LogicalBorderSurface("CylAirSurf", cylPV, airPV,
96:                                Materials::CreateBarSurface());
97: 
99:     new G4LogicalBorderSurface("AirVikuitiSurf", airPV, wrapPV,
100:                                Materials::CreateVikuitiSurface());
101:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E154](#e154) [E316](#e316)
```text
src/Materials.cc
316:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
379:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
401:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
50: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
51: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
66: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E317](#e317). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E152](#e152)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E152](#e152)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E150](#e150)
```text
66b674c feat(ej228): EJ-228 cylinder simulation — 25mm diam × 25mm height
5576687 fix(optics): eliminate group-velocity aliasing bug in EXEC_23 air-gap geometry
610b189 feat(validation): physics validation scan — 7 pos × 500 events
```

**Files touched on the branch since its merge base with main:** [E151](#e151)
```text
 include/DetectorConstruction.hh                  | 114 +++---
 include/Materials.hh                             |  11 +
 include/PrimaryGeneratorAction.hh                |  51 +--
 macros/ej228_run.mac                             |  25 ++
 macros/ej228_vis.mac                             |  26 ++
 macros/validation_scan/validation_01_x-600mm.mac |  18 +
 macros/validation_scan/validation_02_x-400mm.mac |  18 +
 macros/validation_scan/validation_03_x-200mm.mac |  18 +
 macros/validation_scan/validation_04_x0mm.mac    |  18 +
 macros/validation_scan/validation_05_x+200mm.mac |  18 +
 macros/validation_scan/validation_06_x+400mm.mac |  18 +
 macros/validation_scan/validation_07_x+600mm.mac |  18 +
 scripts/analyze_validation.py                    | 211 ++++++++++
 scripts/run_validation_scan.sh                   |  97 +++++
 src/DetectorConstruction.cc                      | 499 ++++++-----------------
 src/Materials.cc                                 |  47 +++
 src/PrimaryGeneratorAction.cc                    | 131 ++----
 src/SiPMSD.cc                                    |   3 +-
 18 files changed, 766 insertions(+), 575 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E157](#e157)
<a id="b9"></a>

### 9. feat/ej228-tir-only
Inspected ref `origin/feat/ej228-tir-only`, full SHA `0006919fe047258228e4eb81e8f96e301458a848`. Classification: **Alternative TIR cylinder; suitability UNVERIFIED**. [E032](#e032)
**B:** ahead 4, behind 103; merged NO; contains 8349041 NO. [E158](#e158) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-15 18:32:40 +0200 / rrios / feat(ej228-tir-only): TIR-only cylinder + Vikuiti comparison analysis and Beamer report. [E032](#e032)
**C1/C2:** active configuration: Cylinder polished dielectric_dielectric TIR boundary. Reflectivity: 0.98/0.95 factories retained, not used by cylinder reflector. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E163](#e163) [E164](#e164) [E318](#e318)

```text
211: G4OpticalSurface* CreateVikuitiSurface() {
216:     surf->SetType(dielectric_metal);
218:     surf->SetFinish(polished);
222:     const std::vector<G4double> refl   = {0.98, 0.98};
226:     mpt->AddProperty("REFLECTIVITY", energy, refl);
251: G4OpticalSurface* CreateBarSurface() {
261:     surf->SetType(dielectric_dielectric);
263:     surf->SetFinish(polished);
269: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
273:     surf->SetType(dielectric_metal);
275:     surf->SetFinish(polished);
294:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
328: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
335:     surf->SetType(dielectric_metal);
337:     surf->SetFinish(polished);
341:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
346:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
371: G4OpticalSurface* CreateBarSkinReflector() {
394:     surf->SetType(dielectric_dielectric);
396:     surf->SetFinish(polished);
401:     const std::vector<G4double> refl   = {0.95, 0.95};
404:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **16**. Full paths, line numbers and contents are preserved in [E165](#e165) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E164](#e164)
```text
78: 
79:     // ── Optical surface: scintillator ↔ air gap (polished → TIR only) ─────────
80:     // No Vikuiti wrap — photons escaping TIR are lost to the world.
83:     new G4LogicalBorderSurface("CylAirSurf", cylPV, airPV,
84:                                Materials::CreateBarSurface());
85: 
129: 
130:     G4cout << "  TIR-only mantle (no Vikuiti), air gap " << kAirGapThick/mm << " mm\n"
131:            << "  TIR angle θ_c = arcsin(1/1.58) = 39.3° on mantle\n\n";
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E163](#e163) [E318](#e318)
```text
src/Materials.cc
316:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
379:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
401:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
50: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
51: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
66: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E319](#e319). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E161](#e161)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E161](#e161)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E159](#e159)
```text
0006919 feat(ej228-tir-only): TIR-only cylinder + Vikuiti comparison analysis and Beamer report
66b674c feat(ej228): EJ-228 cylinder simulation — 25mm diam × 25mm height
5576687 fix(optics): eliminate group-velocity aliasing bug in EXEC_23 air-gap geometry
610b189 feat(validation): physics validation scan — 7 pos × 500 events
```

**Files touched on the branch since its merge base with main:** [E160](#e160)
```text
 analysis/Makefile.analyze                        |  12 +
 analysis/analyze_scan.cxx                        | 419 +++++++++++++++++
 analysis/compare_vikuiti_tir.cxx                 | 330 +++++++++++++
 analysis/timing_resolution.cxx                   | 445 ++++++++++++++++++
 beamer/ej228_vikuiti_vs_tir.tex                  | 568 +++++++++++++++++++++++
 include/DetectorConstruction.hh                  | 114 ++---
 include/Materials.hh                             |  11 +
 include/PrimaryGeneratorAction.hh                |  51 +-
 macros/ej228_run.mac                             |  25 +
 macros/ej228_vis.mac                             |  26 ++
 macros/validation_scan/validation_01_x-600mm.mac |  18 +
 macros/validation_scan/validation_02_x-400mm.mac |  18 +
 macros/validation_scan/validation_03_x-200mm.mac |  18 +
 macros/validation_scan/validation_04_x0mm.mac    |  18 +
 macros/validation_scan/validation_05_x+200mm.mac |  18 +
 macros/validation_scan/validation_06_x+400mm.mac |  18 +
 macros/validation_scan/validation_07_x+600mm.mac |  18 +
 scripts/analyze_validation.py                    | 211 +++++++++
 scripts/run_validation_scan.sh                   |  97 ++++
 src/DetectorConstruction.cc                      | 482 +++++--------------
 src/Materials.cc                                 |  47 ++
 src/PrimaryGeneratorAction.cc                    | 131 ++----
 src/SiPMSD.cc                                    |   3 +-
 23 files changed, 2523 insertions(+), 575 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E166](#e166)
<a id="b10"></a>

### 10. feat/ej230-bar-tir-only
Inspected ref `feat/ej230-bar-tir-only`, full SHA `b281aea71845a115290458928c542dad5c24b5ac`. Classification: **Alternative TIR configuration; suitability UNVERIFIED**. [E032](#e032)
**B:** ahead 34, behind 115; merged NO; contains 8349041 NO. [E067](#e067) [E223](#e223) [E035](#e035)
Upstream: none configured. Origin counterpart: `b281aea71845a115290458928c542dad5c24b5ac`. [E030](#e030)
Tip date / author / subject: 2026-08-15 22:12:30 +0200 / rrios / feat(beamer): add EJ-204 vs EJ-230 TIR-only bar timing comparison report. [E032](#e032)
**C1/C2:** active configuration: polished dielectric_dielectric TIR-only bar skin. Reflectivity: 0.98 factory retained, not used by bar skin. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E167](#e167) [E168](#e168) [E303](#e303)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **6**. Full paths, line numbers and contents are preserved in [E169](#e169) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E168](#e168)
```text
290: 
291:     // TIR-only: polished dielectric-dielectric surface on all bar faces.
292:     // Geant4 Fresnel equations give 100% TIR for theta > arcsin(1/1.58) = 39.3 deg.
294:     // No reflective coating — lateral faces are bare air.
295:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV,
296:                                               Materials::CreateBarSurface());
297:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E167](#e167) [E303](#e303)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E304](#e304). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E070](#e070)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E070](#e070)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E068](#e068)
```text
b281aea feat(beamer): add EJ-204 vs EJ-230 TIR-only bar timing comparison report
a95f610 feat(analysis): add bar timing resolution analysis and run macro
5f30646 feat(ej230-bar-tir-only): EJ-230 scintillator, TIR-only bar
234e0cb feat(ej204-bar-tir-only): EJ-204 bar TIR-only — remove Mylar skin reflector
f19c093 feat(pairscan): CMakeLists refactor, pair scan macros and run script
14ae395 fix(tests): replace escape-fraction guard with sipm-entry guard
c7a627e fix(optics): replace reflector volumes with bar skin surface
021d5f4 EXEC_12TB: add rebuilt self-contained beamer and speaker notes
1ac830a EXEC_12TB: add beamer manifest and Makefile exec12tb targets
d764892 EXEC_12TB: add generated tables and TeX number macros
9cfe042 EXEC_12TB: add figure style module and regenerated deck figures
66ac24c EXEC_12T: add complete Beamer presentation and reproducible build
b0aaac1 EXEC_12T: add integrated timing-position technical report
a2563d9 EXEC_12T: add window-dip timing and collection study
aaabe89 EXEC_12T: add threshold sweep and covariance timing decomposition
0a41034 EXEC_12T: add 4PE versus 20PE temporal and position analysis
f38c9e1 EXEC_12T: add cached fourth-to-thirtieth hit order statistics
003d106 EXEC_12: add staged XY scan proposal
6e01606 EXEC_12: add EXEC11 and EndTop technical reports
30166cd EXEC_12: add y-zero observability feasibility study
a9571fd EXEC_12: add covariance-aware EndTop estimator combination
bd21563 EXEC_12: add leave-one-position-out global X reconstruction
8332733 EXEC_12: add leave-one-position-out global X reconstruction
17d9e7d EXEC_12: add EndTop event observable builder and geometry tests
4dba735 EXEC_11: add reproducible report tables and README
a103ed6 EXEC_11: add temporal ratio and covariance-aware position reconstruction
1e95ac6 EXEC_11: add detailed two-position timing analysis
61bca2e EXEC_11: regenerate pair-scan fit QA and summary v2
bf0d2d5 EXEC_11: add per-event pair observable builder and tests
f431c01 EXEC_08/Step7: CTest guardrails for pair-scan geometry and macros
20bed82 EXEC_08/Step6: batch runner scripts/run_pair_scan.sh
549f3a6 EXEC_08/Step5: ROOT analysis macro analyze_pair_scan.C
30f3da3 EXEC_08/Step3: generate 41 pair-scan macros in macros/pairscan/
6bb73d5 EXEC_08/Step2: geometry analysis + pair selection
```

**Files touched on the branch since its merge base with main:** [E069](#e069)
```text
 .gitignore                                         |    7 +
 CMakeLists.txt                                     |   76 +-
 Makefile                                           |   40 +
 analysis/analyze_pair_scan.C                       |  496 +++++++
 analysis/bar_timing_resolution.cxx                 |  387 ++++++
 analysis/exec11_pair_analysis.py                   |  934 +++++++++++++
 analysis/exec12_endtop_position.py                 |  161 +++
 analysis/exec12_make_reports.py                    |   99 ++
 analysis/exec12t_make_products.py                  |  151 +++
 analysis/exec12t_timing_threshold_analysis.py      |  136 ++
 analysis/exec12tb_figstyle.py                      |   56 +
 analysis/exec12tb_figures.py                       |  823 ++++++++++++
 analysis/exec12tb_tables.py                        |  487 +++++++
 beamer/ej204_vs_ej230_bar_tir.pdf                  |  Bin 0 -> 413342 bytes
 beamer/ej204_vs_ej230_bar_tir.tex                  |  542 ++++++++
 include/DetectorConstruction.hh                    |   10 +-
 macros/pairscan/pairscan_x-422.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-423.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-424.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-425.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-426.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-427.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-428.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-429.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-430.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-431.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-432.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-433.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-434.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-435.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-436.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-437.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-438.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-439.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-440.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-441.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-442.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-443.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-444.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-445.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-446.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-447.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-448.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-449.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-450.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-451.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-452.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-453.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-454.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-455.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-456.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-457.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-458.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-459.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-460.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-461.0mm.mac             |   20 +
 macros/pairscan/pairscan_x-462.0mm.mac             |   20 +
 macros/timing_tir_center.mac                       |   27 +
 pair_scan_config.json                              |   71 +
 report/exec13_xy_scan_plan.md                      |   30 +
 results/exec11_20260612_182454/README.md           |  124 ++
 .../analysis/calibration_ratio.csv                 |    3 +
 .../analysis/calibration_temporal.csv              |    2 +
 .../analysis/calibration_v1_v2.csv                 |    3 +
 .../analysis/data_inventory.csv                    |   42 +
 results/exec11_20260612_182454/analysis/fit_qa.csv |   42 +
 .../analysis/fit_qa_v1_v2_focus.csv                |    4 +
 .../exec11_20260612_182454/analysis/metadata.json  |   13 +
 .../analysis/pairscan_summary_v2.csv               |   42 +
 .../analysis/reconstruction_summary.csv            |   11 +
 .../analysis/reference_comparison.csv              |    3 +
 .../analysis/reference_positions.csv               |    3 +
 .../derived/pair_events_x-422.0mm.npz              |  Bin 0 -> 123534 bytes
 .../derived/pair_events_x-423.0mm.npz              |  Bin 0 -> 122740 bytes
 .../derived/pair_events_x-424.0mm.npz              |  Bin 0 -> 123199 bytes
 .../derived/pair_events_x-425.0mm.npz              |  Bin 0 -> 123147 bytes
 .../derived/pair_events_x-426.0mm.npz              |  Bin 0 -> 123368 bytes
 .../derived/pair_events_x-427.0mm.npz              |  Bin 0 -> 122917 bytes
 .../derived/pair_events_x-428.0mm.npz              |  Bin 0 -> 122839 bytes
 .../derived/pair_events_x-429.0mm.npz              |  Bin 0 -> 122318 bytes
 .../derived/pair_events_x-430.0mm.npz              |  Bin 0 -> 121727 bytes
 .../derived/pair_events_x-431.0mm.npz              |  Bin 0 -> 120954 bytes
 .../derived/pair_events_x-432.0mm.npz              |  Bin 0 -> 120412 bytes
 .../derived/pair_events_x-433.0mm.npz              |  Bin 0 -> 121077 bytes
 .../derived/pair_events_x-434.0mm.npz              |  Bin 0 -> 122386 bytes
 .../derived/pair_events_x-435.0mm.npz              |  Bin 0 -> 122900 bytes
 .../derived/pair_events_x-436.0mm.npz              |  Bin 0 -> 122357 bytes
 .../derived/pair_events_x-437.0mm.npz              |  Bin 0 -> 122009 bytes
 .../derived/pair_events_x-438.0mm.npz              |  Bin 0 -> 122595 bytes
 .../derived/pair_events_x-439.0mm.npz              |  Bin 0 -> 122198 bytes
 .../derived/pair_events_x-440.0mm.npz              |  Bin 0 -> 122232 bytes
 .../derived/pair_events_x-441.0mm.npz              |  Bin 0 -> 122736 bytes
 .../derived/pair_events_x-442.0mm.npz              |  Bin 0 -> 122672 bytes
 .../derived/pair_events_x-443.0mm.npz              |  Bin 0 -> 122553 bytes
 .../derived/pair_events_x-444.0mm.npz              |  Bin 0 -> 122542 bytes
 .../derived/pair_events_x-445.0mm.npz              |  Bin 0 -> 122293 bytes
 .../derived/pair_events_x-446.0mm.npz              |  Bin 0 -> 122306 bytes
 .../derived/pair_events_x-447.0mm.npz              |  Bin 0 -> 121906 bytes
 .../derived/pair_events_x-448.0mm.npz              |  Bin 0 -> 122022 bytes
 .../derived/pair_events_x-449.0mm.npz              |  Bin 0 -> 122871 bytes
 .../derived/pair_events_x-450.0mm.npz              |  Bin 0 -> 122284 bytes
 .../derived/pair_events_x-451.0mm.npz              |  Bin 0 -> 120956 bytes
 .../derived/pair_events_x-452.0mm.npz              |  Bin 0 -> 120882 bytes
 .../derived/pair_events_x-453.0mm.npz              |  Bin 0 -> 120856 bytes
 .../derived/pair_events_x-454.0mm.npz              |  Bin 0 -> 121828 bytes
 .../derived/pair_events_x-455.0mm.npz              |  Bin 0 -> 122035 bytes
 .../derived/pair_events_x-456.0mm.npz              |  Bin 0 -> 122679 bytes
 .../derived/pair_events_x-457.0mm.npz              |  Bin 0 -> 123063 bytes
 .../derived/pair_events_x-458.0mm.npz              |  Bin 0 -> 123181 bytes
 .../derived/pair_events_x-459.0mm.npz              |  Bin 0 -> 123238 bytes
 .../derived/pair_events_x-460.0mm.npz              |  Bin 0 -> 123352 bytes
 .../derived/pair_events_x-461.0mm.npz              |  Bin 0 -> 123315 bytes
 .../derived/pair_events_x-462.0mm.npz              |  Bin 0 -> 123329 bytes
 .../figures/calibrations_and_residuals.pdf         |  Bin 0 -> 25616 bytes
 .../figures/calibrations_and_residuals.png         |  Bin 0 -> 284318 bytes
 .../figures/pos_ref_1_all_hit_times.pdf            |  Bin 0 -> 23571 bytes
 .../figures/pos_ref_1_all_hit_times.png            |  Bin 0 -> 146556 bytes
 .../figures/pos_ref_1_correlations.pdf             |  Bin 0 -> 25033 bytes
 .../figures/pos_ref_1_correlations.png             |  Bin 0 -> 298577 bytes
 .../figures/pos_ref_1_delta_t.pdf                  |  Bin 0 -> 20414 bytes
 .../figures/pos_ref_1_delta_t.png                  |  Bin 0 -> 113077 bytes
 .../figures/pos_ref_1_event_times.pdf              |  Bin 0 -> 20672 bytes
 .../figures/pos_ref_1_event_times.png              |  Bin 0 -> 89908 bytes
 .../figures/pos_ref_1_npe_moyal.pdf                |  Bin 0 -> 26752 bytes
 .../figures/pos_ref_1_npe_moyal.png                |  Bin 0 -> 180911 bytes
 .../figures/pos_ref_2_all_hit_times.pdf            |  Bin 0 -> 23694 bytes
 .../figures/pos_ref_2_all_hit_times.png            |  Bin 0 -> 149031 bytes
 .../figures/pos_ref_2_correlations.pdf             |  Bin 0 -> 24595 bytes
 .../figures/pos_ref_2_correlations.png             |  Bin 0 -> 275444 bytes
 .../figures/pos_ref_2_delta_t.pdf                  |  Bin 0 -> 19869 bytes
 .../figures/pos_ref_2_delta_t.png                  |  Bin 0 -> 114187 bytes
 .../figures/pos_ref_2_event_times.pdf              |  Bin 0 -> 20491 bytes
 .../figures/pos_ref_2_event_times.png              |  Bin 0 -> 88939 bytes
 .../figures/pos_ref_2_npe_moyal.pdf                |  Bin 0 -> 26630 bytes
 .../figures/pos_ref_2_npe_moyal.png                |  Bin 0 -> 177596 bytes
 .../figures/position_reconstruction.pdf            |  Bin 0 -> 23751 bytes
 .../figures/position_reconstruction.png            |  Bin 0 -> 180735 bytes
 results/exec11_20260612_182454/logs/derive.log     |   41 +
 results/exec11_20260612_182454/logs/detail.log     |    0
 results/exec11_20260612_182454/logs/qa.log         |    0
 .../exec11_20260612_182454/logs/reconstruct.log    |    0
 results/exec11_20260612_182454/logs/report.log     |    0
 .../exec11_20260612_182454/tables/fit_qa_v1_v2.tex |    9 +
 .../tables/reconstruction_summary.tex              |   16 +
 .../tables/reference_comparison.tex                |    8 +
 results/exec12_20260612_191000/README.md           |   62 +
 .../analysis/blue_summary.csv                      |   32 +
 .../analysis/blue_weights.csv                      |   32 +
 .../analysis/configuration_provenance.json         |   29 +
 .../analysis/cv_calibrations.csv                   |  187 +++
 .../analysis/cv_model_selection.csv                |    2 +
 .../analysis/cv_predictions.csv.gz                 |  Bin 0 -> 4316374 bytes
 .../analysis/data_inventory.csv                    |   32 +
 .../analysis/x_reconstruction_summary.csv          |  187 +++
 .../analysis/y0_feasibility_summary.csv            |   94 ++
 .../derived/events/events_x+0mm.npz                |  Bin 0 -> 1673253 bytes
 .../derived/events/events_x+100mm.npz              |  Bin 0 -> 1690764 bytes
 .../derived/events/events_x+150mm.npz              |  Bin 0 -> 1686376 bytes
 .../derived/events/events_x+200mm.npz              |  Bin 0 -> 1675705 bytes
 .../derived/events/events_x+250mm.npz              |  Bin 0 -> 1660646 bytes
 .../derived/events/events_x+300mm.npz              |  Bin 0 -> 1641542 bytes
 .../derived/events/events_x+350mm.npz              |  Bin 0 -> 1619395 bytes
 .../derived/events/events_x+400mm.npz              |  Bin 0 -> 1596329 bytes
 .../derived/events/events_x+450mm.npz              |  Bin 0 -> 1564409 bytes
 .../derived/events/events_x+500mm.npz              |  Bin 0 -> 1525987 bytes
 .../derived/events/events_x+50mm.npz               |  Bin 0 -> 1693544 bytes
 .../derived/events/events_x+550mm.npz              |  Bin 0 -> 1488016 bytes
 .../derived/events/events_x+600mm.npz              |  Bin 0 -> 1444699 bytes
 .../derived/events/events_x+650mm.npz              |  Bin 0 -> 1398577 bytes
 .../derived/events/events_x+670mm.npz              |  Bin 0 -> 1376021 bytes
 .../derived/events/events_x+690mm.npz              |  Bin 0 -> 1364036 bytes
 .../derived/events/events_x-100mm.npz              |  Bin 0 -> 1688855 bytes
 .../derived/events/events_x-150mm.npz              |  Bin 0 -> 1685757 bytes
 .../derived/events/events_x-200mm.npz              |  Bin 0 -> 1674171 bytes
 .../derived/events/events_x-250mm.npz              |  Bin 0 -> 1659670 bytes
 .../derived/events/events_x-300mm.npz              |  Bin 0 -> 1643825 bytes
 .../derived/events/events_x-350mm.npz              |  Bin 0 -> 1619374 bytes
 .../derived/events/events_x-400mm.npz              |  Bin 0 -> 1594561 bytes
 .../derived/events/events_x-450mm.npz              |  Bin 0 -> 1560417 bytes
 .../derived/events/events_x-500mm.npz              |  Bin 0 -> 1526242 bytes
 .../derived/events/events_x-50mm.npz               |  Bin 0 -> 1693474 bytes
 .../derived/events/events_x-550mm.npz              |  Bin 0 -> 1487890 bytes
 .../derived/events/events_x-600mm.npz              |  Bin 0 -> 1439585 bytes
 .../derived/events/events_x-650mm.npz              |  Bin 0 -> 1395798 bytes
 .../derived/events/events_x-670mm.npz              |  Bin 0 -> 1377928 bytes
 .../derived/events/events_x-690mm.npz              |  Bin 0 -> 1363389 bytes
 .../exec12_20260612_191000/figures/bias_vs_x.pdf   |  Bin 0 -> 22417 bytes
 .../figures/calibration_end_ratio.pdf              |  Bin 0 -> 15355 bytes
 .../figures/calibration_end_timing.pdf             |  Bin 0 -> 15330 bytes
 .../figures/calibration_top_centroid.pdf           |  Bin 0 -> 15180 bytes
 .../figures/end_channel_symmetry.pdf               |  Bin 0 -> 13810 bytes
 .../residual_distributions_selected_positions.pdf  |  Bin 0 -> 18551 bytes
 .../exec12_20260612_191000/figures/rms68_vs_x.pdf  |  Bin 0 -> 24309 bytes
 .../figures/sigma_core_vs_x.pdf                    |  Bin 0 -> 24945 bytes
 .../figures/valid_fraction_vs_x.pdf                |  Bin 0 -> 19674 bytes
 .../figures/y_centroid_mean_vs_x.pdf               |  Bin 0 -> 16926 bytes
 .../figures/y_centroid_width_vs_x.pdf              |  Bin 0 -> 17309 bytes
 .../figures/y_left_vs_y_right.pdf                  |  Bin 0 -> 18062 bytes
 results/exec12_20260612_191000/logs/analysis.log   |    0
 results/exec12_20260612_191000/logs/inventory.log  |   31 +
 .../report/endtop_position_reconstruction.md       |   59 +
 .../report/endtop_position_reconstruction.tex      |   60 +
 .../report/exec11_technical_note.md                |   42 +
 .../report/exec11_technical_note.tex               |   43 +
 .../tables/x_reconstruction_summary.tex            |  190 +++
 results/exec12t_20260612_195426/README.md          |   48 +
 .../analysis/calibration_20pe.csv                  |    2 +
 .../analysis/calibration_4pe.csv                   |    2 +
 .../analysis/configuration_provenance.json         |   10 +
 .../analysis/data_inventory.csv                    |   42 +
 .../analysis/exec11_reproduction_check.csv         |   42 +
 .../analysis/loo_predictions_20pe.npz              |  Bin 0 -> 892931 bytes
 .../analysis/loo_predictions_4pe.npz               |  Bin 0 -> 892788 bytes
 .../analysis/position_reconstruction_20pe.csv      |   42 +
 .../analysis/position_reconstruction_4pe.csv       |   42 +
 .../analysis/temporal_position_summary.csv         |   83 ++
 .../analysis/threshold_4_20_summary.csv            |    5 +
 .../analysis/threshold_sweep_summary.csv           |   31 +
 .../analysis/window_dip_summary.csv                |    9 +
 .../beamer/beamer_contact_sheet.png                |  Bin 0 -> 473675 bytes
 .../beamer/exec12t_timing_position_beamer.aux      |  128 ++
 .../beamer/exec12t_timing_position_beamer.log      | 1406 ++++++++++++++++++++
 .../beamer/exec12t_timing_position_beamer.nav      |  109 ++
 .../beamer/exec12t_timing_position_beamer.out      |    1 +
 .../beamer/exec12t_timing_position_beamer.pdf      |  Bin 0 -> 234763 bytes
 .../beamer/exec12t_timing_position_beamer.snm      |    0
 .../beamer/exec12t_timing_position_beamer.tex      |  166 +++
 .../beamer/exec12t_timing_position_beamer.toc      |    0
 .../exec12t_20260612_195426/beamer/references.bib  |    1 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 29095 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 47777 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 19663 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 27031 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 15994 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 20580 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 32086 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 26819 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 26911 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 20004 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 20223 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 23610 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 24670 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 22094 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 30569 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 30458 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 18922 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 29235 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 15229 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 19449 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 20846 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 25617 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 24429 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 19255 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 19840 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 24193 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 24012 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 23457 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 31499 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 33489 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 19206 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 27919 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 15664 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 18758 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 21653 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 25116 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 46153 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 21146 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 25827 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 26188 bytes
 .../beamer/rendered/slide-41.png                   |  Bin 0 -> 18935 bytes
 .../beamer/rendered/slide-42.png                   |  Bin 0 -> 19827 bytes
 .../beamer/rendered/slide-43.png                   |  Bin 0 -> 25075 bytes
 .../beamer/rendered/slide-44.png                   |  Bin 0 -> 25549 bytes
 .../beamer/rendered/slide-45.png                   |  Bin 0 -> 22175 bytes
 .../beamer/rendered/slide-46.png                   |  Bin 0 -> 33040 bytes
 .../beamer/rendered/slide-47.png                   |  Bin 0 -> 32296 bytes
 .../beamer/rendered/slide-48.png                   |  Bin 0 -> 20505 bytes
 .../beamer/rendered/slide-49.png                   |  Bin 0 -> 28637 bytes
 .../beamer/speaker_notes.md                        |  107 ++
 .../derived/order_statistics/pair_x-422.0mm.npz    |  Bin 0 -> 1362830 bytes
 .../derived/order_statistics/pair_x-423.0mm.npz    |  Bin 0 -> 1362991 bytes
 .../derived/order_statistics/pair_x-424.0mm.npz    |  Bin 0 -> 1362994 bytes
 .../derived/order_statistics/pair_x-425.0mm.npz    |  Bin 0 -> 1362676 bytes
 .../derived/order_statistics/pair_x-426.0mm.npz    |  Bin 0 -> 1362617 bytes
 .../derived/order_statistics/pair_x-427.0mm.npz    |  Bin 0 -> 1362456 bytes
 .../derived/order_statistics/pair_x-428.0mm.npz    |  Bin 0 -> 1362417 bytes
 .../derived/order_statistics/pair_x-429.0mm.npz    |  Bin 0 -> 1362477 bytes
 .../derived/order_statistics/pair_x-430.0mm.npz    |  Bin 0 -> 1362246 bytes
 .../derived/order_statistics/pair_x-431.0mm.npz    |  Bin 0 -> 1361625 bytes
 .../derived/order_statistics/pair_x-432.0mm.npz    |  Bin 0 -> 1361227 bytes
 .../derived/order_statistics/pair_x-433.0mm.npz    |  Bin 0 -> 1361512 bytes
 .../derived/order_statistics/pair_x-434.0mm.npz    |  Bin 0 -> 1361634 bytes
 .../derived/order_statistics/pair_x-435.0mm.npz    |  Bin 0 -> 1361842 bytes
 .../derived/order_statistics/pair_x-436.0mm.npz    |  Bin 0 -> 1361904 bytes
 .../derived/order_statistics/pair_x-437.0mm.npz    |  Bin 0 -> 1362200 bytes
 .../derived/order_statistics/pair_x-438.0mm.npz    |  Bin 0 -> 1362269 bytes
 .../derived/order_statistics/pair_x-439.0mm.npz    |  Bin 0 -> 1362081 bytes
 .../derived/order_statistics/pair_x-440.0mm.npz    |  Bin 0 -> 1361731 bytes
 .../derived/order_statistics/pair_x-441.0mm.npz    |  Bin 0 -> 1361645 bytes
 .../derived/order_statistics/pair_x-442.0mm.npz    |  Bin 0 -> 1361574 bytes
 .../derived/order_statistics/pair_x-443.0mm.npz    |  Bin 0 -> 1361936 bytes
 .../derived/order_statistics/pair_x-444.0mm.npz    |  Bin 0 -> 1361808 bytes
 .../derived/order_statistics/pair_x-445.0mm.npz    |  Bin 0 -> 1362189 bytes
 .../derived/order_statistics/pair_x-446.0mm.npz    |  Bin 0 -> 1362266 bytes
 .../derived/order_statistics/pair_x-447.0mm.npz    |  Bin 0 -> 1362136 bytes
 .../derived/order_statistics/pair_x-448.0mm.npz    |  Bin 0 -> 1362063 bytes
 .../derived/order_statistics/pair_x-449.0mm.npz    |  Bin 0 -> 1362102 bytes
 .../derived/order_statistics/pair_x-450.0mm.npz    |  Bin 0 -> 1361698 bytes
 .../derived/order_statistics/pair_x-451.0mm.npz    |  Bin 0 -> 1361068 bytes
 .../derived/order_statistics/pair_x-452.0mm.npz    |  Bin 0 -> 1361575 bytes
 .../derived/order_statistics/pair_x-453.0mm.npz    |  Bin 0 -> 1361701 bytes
 .../derived/order_statistics/pair_x-454.0mm.npz    |  Bin 0 -> 1362028 bytes
 .../derived/order_statistics/pair_x-455.0mm.npz    |  Bin 0 -> 1362411 bytes
 .../derived/order_statistics/pair_x-456.0mm.npz    |  Bin 0 -> 1362338 bytes
 .../derived/order_statistics/pair_x-457.0mm.npz    |  Bin 0 -> 1362651 bytes
 .../derived/order_statistics/pair_x-458.0mm.npz    |  Bin 0 -> 1362862 bytes
 .../derived/order_statistics/pair_x-459.0mm.npz    |  Bin 0 -> 1362956 bytes
 .../derived/order_statistics/pair_x-460.0mm.npz    |  Bin 0 -> 1362783 bytes
 .../derived/order_statistics/pair_x-461.0mm.npz    |  Bin 0 -> 1362770 bytes
 .../derived/order_statistics/pair_x-462.0mm.npz    |  Bin 0 -> 1362960 bytes
 .../environment_exec12t.txt                        |   99 ++
 .../figures/correlation_ab_vs_x.pdf                |  Bin 0 -> 12608 bytes
 .../figures/efficiency_4_20_vs_x.pdf               |  Bin 0 -> 12461 bytes
 .../figures/mean_delta_4_20_vs_x.pdf               |  Bin 0 -> 13853 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 9467 bytes
 .../figures/pos_ref_1_delta_t_4_20.pdf             |  Bin 0 -> 15320 bytes
 .../figures/pos_ref_2_delta_t_4_20.pdf             |  Bin 0 -> 15348 bytes
 .../figures/threshold_bias.pdf                     |  Bin 0 -> 15585 bytes
 .../figures/threshold_chi2.pdf                     |  Bin 0 -> 15384 bytes
 .../figures/threshold_efficiency.pdf               |  Bin 0 -> 14268 bytes
 .../figures/threshold_pareto.pdf                   |  Bin 0 -> 14148 bytes
 .../figures/threshold_sigma_dt.pdf                 |  Bin 0 -> 15598 bytes
 .../figures/threshold_sigma_x.pdf                  |  Bin 0 -> 15710 bytes
 .../figures/threshold_slope.pdf                    |  Bin 0 -> 15487 bytes
 .../figures/timing_width_4_20_vs_x.pdf             |  Bin 0 -> 15411 bytes
 results/exec12t_20260612_195426/logs/analysis.log  |    0
 .../logs/beamer_pdfinfo.txt                        |   20 +
 .../logs/beamer_warnings.txt                       |    1 +
 .../logs/exec11_reproduction.log                   |    3 +
 results/exec12t_20260612_195426/logs/inventory.log |   41 +
 results/exec12t_20260612_195426/logs/window.log    |    0
 .../report/exec12t_timing_position_report.aux      |   26 +
 .../report/exec12t_timing_position_report.log      |  376 ++++++
 .../report/exec12t_timing_position_report.md       |   45 +
 .../report/exec12t_timing_position_report.out      |    7 +
 .../report/exec12t_timing_position_report.pdf      |  Bin 0 -> 35001 bytes
 .../report/exec12t_timing_position_report.tex      |   63 +
 .../tables/generated_numbers.tex                   |    6 +
 .../tables/global_context.csv                      |    7 +
 .../tables/global_context.tex                      |   10 +
 .../tables/threshold_4_20_summary.csv              |    7 +
 .../tables/threshold_4_20_summary.tex              |   10 +
 .../tables/threshold_sweep.csv                     |   31 +
 .../tables/threshold_sweep.tex                     |   34 +
 .../tables/window_dip_summary.csv                  |    9 +
 .../tables/window_dip_summary.tex                  |   12 +
 .../beamer/exec12tb_beamer.pdf                     |  Bin 0 -> 544006 bytes
 .../beamer/exec12tb_beamer.tex                     |  971 ++++++++++++++
 .../beamer/generated_numbers.tex                   |   24 +
 .../beamer/rendered/slide-01.png                   |  Bin 0 -> 37541 bytes
 .../beamer/rendered/slide-02.png                   |  Bin 0 -> 55903 bytes
 .../beamer/rendered/slide-03.png                   |  Bin 0 -> 39821 bytes
 .../beamer/rendered/slide-04.png                   |  Bin 0 -> 44919 bytes
 .../beamer/rendered/slide-05.png                   |  Bin 0 -> 37183 bytes
 .../beamer/rendered/slide-06.png                   |  Bin 0 -> 50484 bytes
 .../beamer/rendered/slide-07.png                   |  Bin 0 -> 53853 bytes
 .../beamer/rendered/slide-08.png                   |  Bin 0 -> 54725 bytes
 .../beamer/rendered/slide-09.png                   |  Bin 0 -> 55833 bytes
 .../beamer/rendered/slide-10.png                   |  Bin 0 -> 57120 bytes
 .../beamer/rendered/slide-11.png                   |  Bin 0 -> 55518 bytes
 .../beamer/rendered/slide-12.png                   |  Bin 0 -> 50060 bytes
 .../beamer/rendered/slide-13.png                   |  Bin 0 -> 46531 bytes
 .../beamer/rendered/slide-14.png                   |  Bin 0 -> 53957 bytes
 .../beamer/rendered/slide-15.png                   |  Bin 0 -> 45247 bytes
 .../beamer/rendered/slide-16.png                   |  Bin 0 -> 58873 bytes
 .../beamer/rendered/slide-17.png                   |  Bin 0 -> 56682 bytes
 .../beamer/rendered/slide-18.png                   |  Bin 0 -> 50931 bytes
 .../beamer/rendered/slide-19.png                   |  Bin 0 -> 57910 bytes
 .../beamer/rendered/slide-20.png                   |  Bin 0 -> 49913 bytes
 .../beamer/rendered/slide-21.png                   |  Bin 0 -> 63029 bytes
 .../beamer/rendered/slide-22.png                   |  Bin 0 -> 66176 bytes
 .../beamer/rendered/slide-23.png                   |  Bin 0 -> 49963 bytes
 .../beamer/rendered/slide-24.png                   |  Bin 0 -> 45261 bytes
 .../beamer/rendered/slide-25.png                   |  Bin 0 -> 56098 bytes
 .../beamer/rendered/slide-26.png                   |  Bin 0 -> 68087 bytes
 .../beamer/rendered/slide-27.png                   |  Bin 0 -> 72976 bytes
 .../beamer/rendered/slide-28.png                   |  Bin 0 -> 71762 bytes
 .../beamer/rendered/slide-29.png                   |  Bin 0 -> 45414 bytes
 .../beamer/rendered/slide-30.png                   |  Bin 0 -> 8355 bytes
 .../beamer/rendered/slide-31.png                   |  Bin 0 -> 54362 bytes
 .../beamer/rendered/slide-32.png                   |  Bin 0 -> 54841 bytes
 .../beamer/rendered/slide-33.png                   |  Bin 0 -> 40412 bytes
 .../beamer/rendered/slide-34.png                   |  Bin 0 -> 46539 bytes
 .../beamer/rendered/slide-35.png                   |  Bin 0 -> 61442 bytes
 .../beamer/rendered/slide-36.png                   |  Bin 0 -> 54635 bytes
 .../beamer/rendered/slide-37.png                   |  Bin 0 -> 47345 bytes
 .../beamer/rendered/slide-38.png                   |  Bin 0 -> 50838 bytes
 .../beamer/rendered/slide-39.png                   |  Bin 0 -> 47245 bytes
 .../beamer/rendered/slide-40.png                   |  Bin 0 -> 44446 bytes
 .../beamer/speaker_notes.md                        |  227 ++++
 .../figures/cfd_vs_orderstat_schematic.pdf         |  Bin 0 -> 32384 bytes
 .../figures/fine_scan_geometry.pdf                 |  Bin 0 -> 18315 bytes
 .../figures/global_context_sigma_bias.pdf          |  Bin 0 -> 21642 bytes
 .../figures/mean_dt_vs_x_with_residuals.pdf        |  Bin 0 -> 25199 bytes
 .../figures/order_statistic_schematic.pdf          |  Bin 0 -> 21280 bytes
 .../figures/position_resolution_4_20_vs_x.pdf      |  Bin 0 -> 25708 bytes
 .../figures/ref1_dt_overlay.pdf                    |  Bin 0 -> 25821 bytes
 .../figures/ref1_xrec_overlay.pdf                  |  Bin 0 -> 26514 bytes
 .../figures/ref2_dt_overlay.pdf                    |  Bin 0 -> 25730 bytes
 .../figures/ref2_xrec_overlay.pdf                  |  Bin 0 -> 26561 bytes
 .../figures/rho_ab_vs_x_4_20.pdf                   |  Bin 0 -> 22331 bytes
 .../figures/rms68_dt_vs_x_4_20.pdf                 |  Bin 0 -> 21590 bytes
 .../figures/sha_manifest.json                      |   22 +
 .../figures/sweep_bias_vs_k.pdf                    |  Bin 0 -> 17260 bytes
 .../figures/sweep_chi2_vs_k.pdf                    |  Bin 0 -> 19883 bytes
 .../figures/sweep_rms68x_vs_k.pdf                  |  Bin 0 -> 20698 bytes
 .../figures/sweep_slope_vs_k.pdf                   |  Bin 0 -> 17743 bytes
 .../figures/tplus_rms68_vs_x_4_20.pdf              |  Bin 0 -> 22163 bytes
 .../figures/tradeoff_resolution_vs_bias.pdf        |  Bin 0 -> 27883 bytes
 .../figures/window_dip_counts.pdf                  |  Bin 0 -> 16254 bytes
 .../figures/window_dip_t4_t20.pdf                  |  Bin 0 -> 18771 bytes
 .../logs/beamer_compile1.log                       |  255 ++++
 .../logs/beamer_compile2.log                       |  265 ++++
 .../logs/beamer_compile3.log                       |  275 ++++
 .../logs/beamer_compile4.log                       |  618 +++++++++
 .../logs/beamer_compile5.log                       |  598 +++++++++
 results/exec12tb_20260612_204216/logs/visual_qa.md |   60 +
 .../manifest/beamer_manifest.csv                   |   41 +
 .../tables/dataset_inventory.tex                   |   15 +
 .../tables/global_context.tex                      |   18 +
 .../tables/lv_comparison.tex                       |   23 +
 .../tables/reference_positions.tex                 |   14 +
 .../tables/threshold_4_20_comparison.tex           |   22 +
 .../tables/threshold_sweep.tex                     |   22 +
 .../tables/window_dip_summary.tex                  |   16 +
 scripts/analyze_geometry.py                        |  172 +++
 scripts/build_exec12t.sh                           |   13 +
 scripts/gen_pair_scan_macros.py                    |  119 ++
 scripts/make_contact_sheet.py                      |   36 +
 scripts/run_pair_scan.sh                           |  156 +++
 src/DetectorConstruction.cc                        |  125 +-
 tests/check_endtop_balance.py                      |   14 +-
 tests/check_pairscan_geometry.py                   |  123 ++
 tests/check_pairscan_macros.py                     |   81 ++
 tests/readout_config_check.cc                      |   53 +-
 tests/test_exec11_pair_analysis.py                 |   64 +
 tests/test_exec12_endtop_position.py               |   20 +
 tests/test_exec12t_timing_threshold_analysis.py    |   15 +
 448 files changed, 14226 insertions(+), 135 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E170](#e170)
<a id="b11"></a>

### 11. feat/ej230-endonly-mylar
Inspected ref `origin/feat/ej230-endonly-mylar`, full SHA `04a80473a2ed9aec206af9e0eb1c70a7e25f84d1`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 42, behind 115; merged NO; contains 8349041 NO. [E171](#e171) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-31 20:32:44 +0200 / dowiyogo / fix(beamer): manual layout corrections to ej230 endonly-mylar deck. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal ground Mylar surface. Reflectivity: 0.90 active; 0.98 unused bar factory retained. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E176](#e176) [E177](#e177) [E320](#e320)

```text
205: G4OpticalSurface* CreateBarSurface() {
215:     surf->SetType(dielectric_dielectric);
217:     surf->SetFinish(polished);
223: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
229:     surf->SetType(dielectric_metal);
231:     surf->SetFinish(polished);
250:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
296: G4OpticalSurface* CreateBarSkinReflector() {
315:     surf->SetType(dielectric_metal);
317:     surf->SetFinish(polished);
324:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
333:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
339: G4OpticalSurface* CreateMylarReflector() {
341:     surf->SetType(dielectric_metal);
343:     surf->SetFinish(ground);
347:     const std::vector<G4double> uniformReflectivity = {0.90, 0.90};
352:     mpt->AddProperty("REFLECTIVITY", energy, uniformReflectivity);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **195**. Full paths, line numbers and contents are preserved in [E178](#e178) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E177](#e177)
```text
285:     // volumes avoid global skin surfaces shadowing active SiPM faces.
286:     auto* reflector = Materials::CreateMylarReflector();
287:     const G4double foilHalfT = 0.5 * um;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E176](#e176) [E320](#e320)
```text
src/Materials.cc
308:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
53: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E321](#e321). Raw matches below are historical finish-test documentation or percentages, not the requested configuration/result marker.
```text
origin/feat/ej230-endonly-mylar:docs/exec14d_adaptive_tN.md:21:| -450 | 20 | 14 | 23.69 | 69.2% | 95.7% | reduced |
origin/feat/ej230-endonly-mylar:results_ej230_analysis/logs/exec14b_pdftotext.txt:1784:69.2%
origin/feat/ej230-endonly-mylar:results_ej230_analysis/tables/adaptive_tN.tex:13:-450 & 4 & 4 & 14 & 69.2\% \\
```

**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E174](#e174)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E174](#e174)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E172](#e172)
```text
04a8047 fix(beamer): manual layout corrections to ej230 endonly-mylar deck
98af271 feat(exec14b/14d): add EXEC14b/14d analysis scripts and t0minidaq adaptations
e4110ce DIAGNOSTIC: Add PrintGroupVelDiagnostic() to investigate GROUPVEL propagation velocity anomaly
f22de30 EXEC_14: análisis completo EJ-230 End-only+Mylar (ADDENDUM FINAL aplicado)
7e18117 fix(readout): elimina selector residual de midpoint TOP
5c85a40 docs: README + nota AFBR P024M/P014M
ef8c4de chore(scripts): run_msi.sh (16t) y run_t0minidaq.sh (24t)
0c02e00 build: GDML opcional (compila sin GDML)
b1c738d refactor(readout): canales solo END; limpia referencias TOP
5668bab feat(geom): Mylar en todas las caras salvo extremos (misma def. que feat/endonly-mylar)
8e2b351 feat(geom): elimina TOP SiPMs y sus ventanas (réplica feat/endonly-mylar, EJ230)
97bdc8e chore: worktree/branch feat/ej230-endonly-mylar desde feat/ej230-sslg4
79e701d EXEC_14E/T5: compile and validate final report
5ddc3e0 EXEC_14E/T4: append legible adaptive t_N backup frames
144ed92 EXEC_14E/T3: validate appended adaptive backup frames
153be3c EXEC_14E/T2: move adaptive t_N material to appended backup
ebbb630 EXEC_14E/T1: restore fixed-bin six-panel t_N displays
8a4d338 EXEC_14D/T6: recompilar y validar reporte final
80fd1fd EXEC_14D/T5: pulir coeficiente analítico y errores Npe
2ae390c EXEC_14D/T4: reescribir diagnóstico y restaurar cajas pedagógicas
3333a8f EXEC_14D/T3: aplicar t_N adaptativo y anotaciones
1379d4b EXEC_14D/T2: unificar fitted sigma(t4)
2a6017b EXEC_14D/T1: diagnosticar EndTop contra End-only
cb4b912 EXEC_14C/T5: preflight + render + documentación
0e70266 EXEC_14C/T4: tablas, reconstrucción y compilación estricta del Beamer
0b8be52 EXEC_14C/T3: regenerar EXEC_09 tail (figura/CSV/tablas EJ-230 auténticas)
34e67a5 EXEC_14C/T2: verificar propagación de constantes 0.5/1.5 ns
41b25d9 EXEC_14C/T1: completar End-only x=0,+400 y validar ROOT
533b360 EXEC_14B WIP: estado interrumpido (CODEX), pre-reanudación
4842079 EXEC_14B/T1: audit Beamer assets and identify broken paths
be6743e EXEC_14/T5: hallazgos documentation
f6aa121 EXEC_14/T4: Beamer compiled + frame parity verified
0b26c5f EXEC_14/T3: EJ-230 analysis products — 161 figs, CSVs, tables
c16f679 EXEC_14/T2: analysis route config — exec12b --tau-d arg + EJ-230 audit files
5dad9e3 EXEC_14/T1: data inventory — 31/31 EJ-230 ROOT files validated
192c4fe EXEC_13/T7: RUNBOOK_EJ230.md — per-machine copy-paste instructions
09c7b8c EXEC_13/T6: exec13_ej230_report_full.tex — 1:1 Beamer mirror (119 frames)
e162fde EXEC_13/T5: headless analysis pipeline for t0minidaq
b53d688 EXEC_13/T4: machine-specific scan and center scripts (24t/16t)
7e7e00e EXEC_13/T3: GDML test conditional via ENABLE_GDML_TEST CMake option
4f4f8fe EXEC_13/T2: tee sim output to terminal in run script
4be50d8 EXEC_13/T1: switch scintillator OPSC-101→OPSC-106 (EJ-204→EJ-230)
```

**Files touched on the branch since its merge base with main:** [E173](#e173)
```text
 .gitignore                                         |    7 +
 CMakeLists.txt                                     |   28 +-
 README.md                                          |  160 +-
 RUNBOOK_EJ230.md                                   |  123 +
 analysis/endonly_sum4.py                           |  432 +++
 analysis/exec07/common.py                          |    4 +-
 analysis/exec07/exec09_timing_mechanism.py         |    6 +-
 analysis/exec07/exec10_landau_analysis.py          |    2 +
 analysis/exec07/exec10_time_figures.py             |    2 +
 analysis/exec07/exec10_veff_diagnosis.py           |    2 +
 analysis/exec07/exec11_time_arrival.py             |    2 +-
 analysis/exec07/exec12b_tn_dispersion.py           |    8 +
 analysis/exec07/exec13_ej230_report_full.tex       | 1696 +++++++++++
 analysis/exec07_photon_budget.py                   |    8 +-
 analysis/exec13/__init__.py                        |    0
 analysis/exec13/common13.py                        |   71 +
 analysis/exec13/exec13_tN_analysis.py              |  906 ++++++
 analysis_ej230/DATASET_INVENTORY.md                |  110 +
 analysis_ej230/EJ230_MATERIAL_AUDIT.md             |  116 +
 analysis_ej230/FINAL_PHYSICS_AUDIT.md              |  235 ++
 analysis_ej230/csv/attenuation_bootstrap_left.csv  |    8 +
 .../csv/attenuation_bootstrap_replicates.csv       |  201 ++
 .../csv/attenuation_bootstrap_replicates_left.csv  |  201 ++
 .../csv/attenuation_bootstrap_replicates_right.csv |  201 ++
 analysis_ej230/csv/attenuation_bootstrap_right.csv |    8 +
 .../csv/attenuation_bootstrap_summary.csv          |    9 +
 .../csv/attenuation_bootstrap_summary_left.csv     |    9 +
 .../csv/attenuation_bootstrap_summary_right.csv    |    9 +
 analysis_ej230/csv/attenuation_cm.csv              |   32 +
 .../csv/attenuation_model_comparison.csv           |   11 +
 .../csv/attenuation_points_endonly_mylar_ej230.csv |   32 +
 .../csv/attenuation_points_endtop_ej230.csv        |   32 +
 .../csv/attenuation_range_sensitivity.csv          |    7 +
 analysis_ej230/csv/far_end_light_yield.csv         |    3 +
 analysis_ej230/csv/group_velocity_validation.csv   |    6 +
 .../csv/photon_budget_sensor_conversion.csv        |    6 +
 .../csv/photon_budget_terminal_fates.csv           |    6 +
 analysis_ej230/csv/timing_estimator_slopes.csv     |    5 +
 .../csv/timing_resolution_by_position.csv          |   32 +
 analysis_ej230/logs/analysis_run.log               |  104 +
 analysis_ej230/logs/bootstrap_cm_runtime.log       |   13 +
 analysis_ej230/logs/bootstrap_runtime.log          |   13 +
 .../scripts/analysis_ej230_endonly_mylar.py        | 1610 +++++++++++
 analysis_ej230/scripts/bootstrap_attenuation_ej230 |  Bin 0 -> 65648 bytes
 .../scripts/bootstrap_attenuation_openmp.cpp       |  456 +++
 .../scripts/validate_group_velocity_ej230.py       |  100 +
 audit/exec10_landau_analysis.md                    |    6 +-
 audit/exec10_veff_diagnosis.md                     |    8 +-
 audit/exec11_arrival.md                            |   44 +-
 docs/EJ230_material_audit.md                       |  235 ++
 docs/exec14_data_inventory.md                      |   56 +
 docs/exec14_fit_failures.md                        |   63 +
 docs/exec14b_beamer_audit.md                       |   52 +
 docs/exec14b_missing_analysis.md                   |   15 +
 docs/exec14c_state.md                              |   84 +
 docs/exec14d_adaptive_tN.md                        |   49 +
 docs/exec14d_endtop_diagnosis.md                   |   40 +
 docs/exec14d_title_exceptions.md                   |   16 +
 docs/exec14e_tN_revert.md                          |   82 +
 include/DetectorConstruction.hh                    |   32 +-
 include/EventAction.hh                             |    2 -
 include/Materials.hh                               |    3 +
 include/PrimaryGeneratorAction.hh                  |   12 +-
 include/RunAction.hh                               |    2 -
 include/SiPMSD.hh                                  |    3 +-
 macros/edge_scan_smoke.mac                         |    3 +-
 macros/endtop_smoke_center.mac                     |   17 -
 macros/endtop_smoke_edge.mac                       |   17 -
 macros/event_display_top_3d.mac                    |   38 -
 macros/event_display_top_batch.mac                 |   20 -
 macros/event_display_top_lateral.mac               |   39 -
 macros/event_display_top_midpoint.mac              |   49 -
 macros/exec08b_run_a_window_center.mac             |   16 -
 macros/exec08b_run_b_window_midpoint.mac           |   16 -
 macros/exec08b_run_c1_window_mirror.mac            |   16 -
 macros/exec08b_run_c2_window_exact_mirror.mac      |   16 -
 macros/run.mac                                     |    3 -
 macros/scan.mac                                    |    1 -
 macros/scan_angle.mac                              |   24 -
 main.cc                                            |   68 +-
 .../endonly_mylar_ej230/build_deck_ej230.py        | 1063 +++++++
 presentation/endonly_mylar_ej230/deck_values.json  |  607 ++++
 .../figures/fig_attenuation.png                    |  Bin 0 -> 123804 bytes
 .../figures/fig_attenuation.txt                    |    1 +
 .../figures/fig_attenuation_twocomp.png            |  Bin 0 -> 106944 bytes
 .../figures/fig_attenuation_twocomp.txt            |    1 +
 .../figures/fig_endtop_comparison.png              |  Bin 0 -> 95161 bytes
 .../figures/fig_endtop_comparison.txt              |    1 +
 .../figures/fig_far_end_yield.png                  |  Bin 0 -> 39947 bytes
 .../figures/fig_far_end_yield.txt                  |    1 +
 .../endonly_mylar_ej230/figures/fig_fpt_vs_t50.png |  Bin 0 -> 188299 bytes
 .../endonly_mylar_ej230/figures/fig_fpt_vs_t50.txt |    1 +
 .../figures/fig_mean_arrival_time.png              |  Bin 0 -> 113282 bytes
 .../figures/fig_mean_arrival_time.txt              |    1 +
 .../endonly_mylar_ej230/figures/fig_npe_vs_x.png   |  Bin 0 -> 73036 bytes
 .../endonly_mylar_ej230/figures/fig_npe_vs_x.txt   |    1 +
 .../figures/fig_photon_budget.png                  |  Bin 0 -> 89586 bytes
 .../figures/fig_photon_budget.txt                  |    1 +
 .../endonly_mylar_ej230/figures/fig_sigma_t.png    |  Bin 0 -> 67685 bytes
 .../endonly_mylar_ej230/figures/fig_sigma_t.txt    |    1 +
 .../figures/fig_timing_histograms.png              |  Bin 0 -> 151217 bytes
 .../figures/fig_timing_histograms.txt              |    1 +
 presentation/endonly_mylar_ej230/main.tex          |  770 +++++
 .../tables/tab_attenuation_fits.tex                |   10 +
 .../endonly_mylar_ej230/tables/tab_sigma_t.tex     |    9 +
 .../endonly_mylar_ej230/tables/tab_sim_params.tex  |   24 +
 .../endonly_mylar_ej230/verify_deck_ej230.py       |  229 ++
 results_ej230_analysis/conclusions_exec07.md       |   19 +
 results_ej230_analysis/csv/exec08b_timing_gate.csv |    7 +
 .../csv/exec08b_timing_gate_raw.csv                |   13 +
 .../csv/exec08b_window_dip_profiles.csv            |   56 +
 .../csv/exec09_tail_comparison.csv                 |   16 +
 results_ej230_analysis/csv/exec09_tail_metrics.csv |   31 +
 .../csv/exec09_timing_verdict.txt                  |    6 +
 results_ej230_analysis/csv/exec13_tN_summary.csv   |   94 +
 results_ej230_analysis/csv/exec14d_adaptive_tN.csv |   94 +
 results_ej230_analysis/csv/exec14d_endtop_diag.csv |    7 +
 .../csv/exec14d_nominal_parameters.csv             |    2 +
 .../csv/exec14e_fixed_tN_summary.csv               |  187 ++
 results_ej230_analysis/exec10_fano_by_channel.csv  | 2667 +++++++++++++++++
 results_ej230_analysis/exec10_fano_fit.csv         |    2 +
 results_ej230_analysis/exec10_landau_mpv.csv       |  603 ++++
 .../exec10_late_fraction_representative.csv        |    9 +
 results_ej230_analysis/exec10_velocity_fits.csv    |   10 +
 results_ej230_analysis/exec10_velocity_metrics.csv |   63 +
 results_ej230_analysis/exec11_arrival_metrics.csv  |   94 +
 results_ej230_analysis/figs/P1_npe_vs_x.png        |  Bin 0 -> 233001 bytes
 results_ej230_analysis/figs/P2_npe_heatmap_top.png |  Bin 0 -> 118849 bytes
 results_ej230_analysis/figs/P3_fano_vs_x.png       |  Bin 0 -> 318907 bytes
 .../figs/P4_poisson_check_x-400.png                |  Bin 0 -> 199044 bytes
 .../figs/P4_poisson_check_x-690.png                |  Bin 0 -> 203861 bytes
 .../figs/P4_poisson_check_x0.png                   |  Bin 0 -> 229756 bytes
 .../figs/P4_poisson_check_x400.png                 |  Bin 0 -> 202640 bytes
 .../figs/P4_poisson_check_x690.png                 |  Bin 0 -> 198169 bytes
 results_ej230_analysis/figs/P5_tdist_examples.png  |  Bin 0 -> 255944 bytes
 results_ej230_analysis/figs/P6_tmean_vs_x.png      |  Bin 0 -> 258753 bytes
 results_ej230_analysis/figs/P7_deltaT_end.png      |  Bin 0 -> 270072 bytes
 .../figs/exec08b_id18_impact_maps.png              |  Bin 0 -> 374058 bytes
 .../figs/exec08b_window_dip_profiles.png           |  Bin 0 -> 489245 bytes
 .../figs/exec09_tail_comparison.png                |  Bin 0 -> 340011 bytes
 .../figs/exec10_fano_vs_mean.png                   |  Bin 0 -> 233930 bytes
 .../figs/exec10_landau_fit_example.png             |  Bin 0 -> 100479 bytes
 .../figs/exec10_velocity_estimators.png            |  Bin 0 -> 280181 bytes
 .../figs/exec11_arrival_-100mm.png                 |  Bin 0 -> 204685 bytes
 .../figs/exec11_arrival_-150mm.png                 |  Bin 0 -> 214515 bytes
 .../figs/exec11_arrival_-200mm.png                 |  Bin 0 -> 207191 bytes
 .../figs/exec11_arrival_-250mm.png                 |  Bin 0 -> 211922 bytes
 .../figs/exec11_arrival_-300mm.png                 |  Bin 0 -> 204030 bytes
 .../figs/exec11_arrival_-350mm.png                 |  Bin 0 -> 210514 bytes
 .../figs/exec11_arrival_-400mm.png                 |  Bin 0 -> 207310 bytes
 .../figs/exec11_arrival_-450mm.png                 |  Bin 0 -> 210300 bytes
 .../figs/exec11_arrival_-500mm.png                 |  Bin 0 -> 206720 bytes
 .../figs/exec11_arrival_-50mm.png                  |  Bin 0 -> 208467 bytes
 .../figs/exec11_arrival_-550mm.png                 |  Bin 0 -> 219071 bytes
 .../figs/exec11_arrival_-600mm.png                 |  Bin 0 -> 203660 bytes
 .../figs/exec11_arrival_-650mm.png                 |  Bin 0 -> 216676 bytes
 .../figs/exec11_arrival_-670mm.png                 |  Bin 0 -> 225139 bytes
 .../figs/exec11_arrival_-690mm.png                 |  Bin 0 -> 229147 bytes
 results_ej230_analysis/figs/exec11_arrival_0mm.png |  Bin 0 -> 191246 bytes
 .../figs/exec11_arrival_100mm.png                  |  Bin 0 -> 203459 bytes
 .../figs/exec11_arrival_150mm.png                  |  Bin 0 -> 212372 bytes
 .../figs/exec11_arrival_200mm.png                  |  Bin 0 -> 207820 bytes
 .../figs/exec11_arrival_250mm.png                  |  Bin 0 -> 211887 bytes
 .../figs/exec11_arrival_300mm.png                  |  Bin 0 -> 204235 bytes
 .../figs/exec11_arrival_350mm.png                  |  Bin 0 -> 210990 bytes
 .../figs/exec11_arrival_400mm.png                  |  Bin 0 -> 208226 bytes
 .../figs/exec11_arrival_450mm.png                  |  Bin 0 -> 209545 bytes
 .../figs/exec11_arrival_500mm.png                  |  Bin 0 -> 203660 bytes
 .../figs/exec11_arrival_50mm.png                   |  Bin 0 -> 207820 bytes
 .../figs/exec11_arrival_550mm.png                  |  Bin 0 -> 218482 bytes
 .../figs/exec11_arrival_600mm.png                  |  Bin 0 -> 203211 bytes
 .../figs/exec11_arrival_650mm.png                  |  Bin 0 -> 214913 bytes
 .../figs/exec11_arrival_670mm.png                  |  Bin 0 -> 221144 bytes
 .../figs/exec11_arrival_690mm.png                  |  Bin 0 -> 226340 bytes
 results_ej230_analysis/figs/exec13_tN_-400mm.png   |  Bin 0 -> 227690 bytes
 results_ej230_analysis/figs/exec13_tN_-650mm.png   |  Bin 0 -> 241442 bytes
 results_ej230_analysis/figs/exec13_tN_-690mm.png   |  Bin 0 -> 235239 bytes
 results_ej230_analysis/figs/exec13_tN_0mm.png      |  Bin 0 -> 246137 bytes
 results_ej230_analysis/figs/exec13_tN_400mm.png    |  Bin 0 -> 235145 bytes
 results_ej230_analysis/figs/exec13_tN_650mm.png    |  Bin 0 -> 241645 bytes
 results_ej230_analysis/figs/exec13_tN_690mm.png    |  Bin 0 -> 238605 bytes
 results_ej230_analysis/figs/exec13_tN_summary.png  |  Bin 0 -> 105222 bytes
 .../figs/exec13_tn_-400mm_endL.png                 |  Bin 0 -> 164703 bytes
 .../figs/exec13_tn_-400mm_endR.png                 |  Bin 0 -> 161097 bytes
 .../figs/exec13_tn_-400mm_top.png                  |  Bin 0 -> 212101 bytes
 .../figs/exec13_tn_-650mm_endL.png                 |  Bin 0 -> 168547 bytes
 .../figs/exec13_tn_-650mm_endR.png                 |  Bin 0 -> 131747 bytes
 .../figs/exec13_tn_-650mm_top.png                  |  Bin 0 -> 213667 bytes
 .../figs/exec13_tn_-690mm_endL.png                 |  Bin 0 -> 184723 bytes
 .../figs/exec13_tn_-690mm_endR.png                 |  Bin 0 -> 128295 bytes
 .../figs/exec13_tn_-690mm_top.png                  |  Bin 0 -> 212919 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_endL.png |  Bin 0 -> 153867 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_endR.png |  Bin 0 -> 156515 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_top.png  |  Bin 0 -> 210186 bytes
 .../figs/exec13_tn_400mm_endL.png                  |  Bin 0 -> 156488 bytes
 .../figs/exec13_tn_400mm_endR.png                  |  Bin 0 -> 165674 bytes
 .../figs/exec13_tn_400mm_top.png                   |  Bin 0 -> 211626 bytes
 .../figs/exec13_tn_650mm_endL.png                  |  Bin 0 -> 136647 bytes
 .../figs/exec13_tn_650mm_endR.png                  |  Bin 0 -> 170238 bytes
 .../figs/exec13_tn_650mm_top.png                   |  Bin 0 -> 215272 bytes
 .../figs/exec13_tn_690mm_endL.png                  |  Bin 0 -> 128415 bytes
 .../figs/exec13_tn_690mm_endR.png                  |  Bin 0 -> 185895 bytes
 .../figs/exec13_tn_690mm_top.png                   |  Bin 0 -> 213888 bytes
 .../figs/exec14d_adaptive_tN_-400mm.png            |  Bin 0 -> 161061 bytes
 .../figs/exec14d_adaptive_tN_-650mm.png            |  Bin 0 -> 167773 bytes
 .../figs/exec14d_adaptive_tN_-690mm.png            |  Bin 0 -> 162200 bytes
 .../figs/exec14d_adaptive_tN_0mm.png               |  Bin 0 -> 157373 bytes
 .../figs/exec14d_adaptive_tN_400mm.png             |  Bin 0 -> 168065 bytes
 .../figs/exec14d_adaptive_tN_650mm.png             |  Bin 0 -> 170074 bytes
 .../figs/exec14d_adaptive_tN_690mm.png             |  Bin 0 -> 166546 bytes
 .../figs/exec14d_adaptive_tN_summary.png           |  Bin 0 -> 124988 bytes
 .../figs/fits_attenuation_velocity.png             |  Bin 0 -> 239898 bytes
 .../figs/muon_-100mm_geometry.png                  |  Bin 0 -> 86814 bytes
 .../figs/muon_-100mm_top_profile.png               |  Bin 0 -> 130050 bytes
 .../figs/muon_-150mm_geometry.png                  |  Bin 0 -> 86584 bytes
 .../figs/muon_-150mm_top_profile.png               |  Bin 0 -> 130972 bytes
 .../figs/muon_-200mm_geometry.png                  |  Bin 0 -> 86031 bytes
 .../figs/muon_-200mm_top_profile.png               |  Bin 0 -> 128479 bytes
 .../figs/muon_-250mm_geometry.png                  |  Bin 0 -> 88413 bytes
 .../figs/muon_-250mm_top_profile.png               |  Bin 0 -> 131490 bytes
 .../figs/muon_-300mm_geometry.png                  |  Bin 0 -> 87067 bytes
 .../figs/muon_-300mm_top_profile.png               |  Bin 0 -> 130048 bytes
 .../figs/muon_-350mm_geometry.png                  |  Bin 0 -> 86576 bytes
 .../figs/muon_-350mm_top_profile.png               |  Bin 0 -> 130743 bytes
 .../figs/muon_-400mm_geometry.png                  |  Bin 0 -> 86364 bytes
 .../figs/muon_-400mm_top_profile.png               |  Bin 0 -> 126950 bytes
 .../figs/muon_-450mm_geometry.png                  |  Bin 0 -> 85343 bytes
 .../figs/muon_-450mm_top_profile.png               |  Bin 0 -> 127479 bytes
 .../figs/muon_-500mm_geometry.png                  |  Bin 0 -> 86759 bytes
 .../figs/muon_-500mm_top_profile.png               |  Bin 0 -> 126793 bytes
 .../figs/muon_-50mm_geometry.png                   |  Bin 0 -> 87030 bytes
 .../figs/muon_-50mm_top_profile.png                |  Bin 0 -> 130249 bytes
 .../figs/muon_-550mm_geometry.png                  |  Bin 0 -> 85110 bytes
 .../figs/muon_-550mm_top_profile.png               |  Bin 0 -> 124255 bytes
 .../figs/muon_-600mm_geometry.png                  |  Bin 0 -> 86326 bytes
 .../figs/muon_-600mm_top_profile.png               |  Bin 0 -> 122760 bytes
 .../figs/muon_-650mm_geometry.png                  |  Bin 0 -> 85418 bytes
 .../figs/muon_-650mm_top_profile.png               |  Bin 0 -> 118836 bytes
 .../figs/muon_-670mm_geometry.png                  |  Bin 0 -> 85573 bytes
 .../figs/muon_-670mm_top_profile.png               |  Bin 0 -> 115700 bytes
 .../figs/muon_-690mm_geometry.png                  |  Bin 0 -> 86032 bytes
 .../figs/muon_-690mm_top_profile.png               |  Bin 0 -> 115351 bytes
 results_ej230_analysis/figs/muon_0mm_geometry.png  |  Bin 0 -> 84816 bytes
 .../figs/muon_0mm_top_profile.png                  |  Bin 0 -> 128228 bytes
 .../figs/muon_100mm_geometry.png                   |  Bin 0 -> 85800 bytes
 .../figs/muon_100mm_top_profile.png                |  Bin 0 -> 129881 bytes
 .../figs/muon_150mm_geometry.png                   |  Bin 0 -> 87145 bytes
 .../figs/muon_150mm_top_profile.png                |  Bin 0 -> 132085 bytes
 .../figs/muon_200mm_geometry.png                   |  Bin 0 -> 88515 bytes
 .../figs/muon_200mm_top_profile.png                |  Bin 0 -> 130820 bytes
 .../figs/muon_250mm_geometry.png                   |  Bin 0 -> 88327 bytes
 .../figs/muon_250mm_top_profile.png                |  Bin 0 -> 131353 bytes
 .../figs/muon_300mm_geometry.png                   |  Bin 0 -> 87084 bytes
 .../figs/muon_300mm_top_profile.png                |  Bin 0 -> 129045 bytes
 .../figs/muon_350mm_geometry.png                   |  Bin 0 -> 88801 bytes
 .../figs/muon_350mm_top_profile.png                |  Bin 0 -> 131357 bytes
 .../figs/muon_400mm_geometry.png                   |  Bin 0 -> 86811 bytes
 .../figs/muon_400mm_top_profile.png                |  Bin 0 -> 127528 bytes
 .../figs/muon_450mm_geometry.png                   |  Bin 0 -> 85908 bytes
 .../figs/muon_450mm_top_profile.png                |  Bin 0 -> 127149 bytes
 .../figs/muon_500mm_geometry.png                   |  Bin 0 -> 86543 bytes
 .../figs/muon_500mm_top_profile.png                |  Bin 0 -> 124196 bytes
 results_ej230_analysis/figs/muon_50mm_geometry.png |  Bin 0 -> 86568 bytes
 .../figs/muon_50mm_top_profile.png                 |  Bin 0 -> 129626 bytes
 .../figs/muon_550mm_geometry.png                   |  Bin 0 -> 87039 bytes
 .../figs/muon_550mm_top_profile.png                |  Bin 0 -> 124104 bytes
 .../figs/muon_600mm_geometry.png                   |  Bin 0 -> 88127 bytes
 .../figs/muon_600mm_top_profile.png                |  Bin 0 -> 123025 bytes
 .../figs/muon_650mm_geometry.png                   |  Bin 0 -> 86700 bytes
 .../figs/muon_650mm_top_profile.png                |  Bin 0 -> 118119 bytes
 .../figs/muon_670mm_geometry.png                   |  Bin 0 -> 87021 bytes
 .../figs/muon_670mm_top_profile.png                |  Bin 0 -> 115175 bytes
 .../figs/muon_690mm_geometry.png                   |  Bin 0 -> 87513 bytes
 .../figs/muon_690mm_top_profile.png                |  Bin 0 -> 114202 bytes
 results_ej230_analysis/figs/tn_-400mm_endL.png     |  Bin 0 -> 164703 bytes
 results_ej230_analysis/figs/tn_-400mm_endR.png     |  Bin 0 -> 161097 bytes
 results_ej230_analysis/figs/tn_-400mm_top.png      |  Bin 0 -> 212101 bytes
 results_ej230_analysis/figs/tn_-650mm_endL.png     |  Bin 0 -> 168547 bytes
 results_ej230_analysis/figs/tn_-650mm_endR.png     |  Bin 0 -> 131747 bytes
 results_ej230_analysis/figs/tn_-650mm_top.png      |  Bin 0 -> 213667 bytes
 results_ej230_analysis/figs/tn_-690mm_endL.png     |  Bin 0 -> 184723 bytes
 results_ej230_analysis/figs/tn_-690mm_endR.png     |  Bin 0 -> 128295 bytes
 results_ej230_analysis/figs/tn_-690mm_top.png      |  Bin 0 -> 212919 bytes
 results_ej230_analysis/figs/tn_0mm_endL.png        |  Bin 0 -> 153867 bytes
 results_ej230_analysis/figs/tn_0mm_endR.png        |  Bin 0 -> 156515 bytes
 results_ej230_analysis/figs/tn_0mm_top.png         |  Bin 0 -> 210186 bytes
 results_ej230_analysis/figs/tn_400mm_endL.png      |  Bin 0 -> 156488 bytes
 results_ej230_analysis/figs/tn_400mm_endR.png      |  Bin 0 -> 165674 bytes
 results_ej230_analysis/figs/tn_400mm_top.png       |  Bin 0 -> 211626 bytes
 results_ej230_analysis/figs/tn_650mm_endL.png      |  Bin 0 -> 136647 bytes
 results_ej230_analysis/figs/tn_650mm_endR.png      |  Bin 0 -> 170238 bytes
 results_ej230_analysis/figs/tn_650mm_top.png       |  Bin 0 -> 215272 bytes
 results_ej230_analysis/figs/tn_690mm_endL.png      |  Bin 0 -> 128415 bytes
 results_ej230_analysis/figs/tn_690mm_endR.png      |  Bin 0 -> 185895 bytes
 results_ej230_analysis/figs/tn_690mm_top.png       |  Bin 0 -> 213888 bytes
 results_ej230_analysis/fit_results_exec07.csv      |    3 +
 .../logs/exec07_photon_budget.log                  |   34 +
 results_ej230_analysis/logs/exec10_landau.log      |   42 +
 .../logs/exec10_time_figures.log                   |    1 +
 results_ej230_analysis/logs/exec10_veff.log        |   47 +
 results_ej230_analysis/logs/exec11_arrival.log     |   94 +
 results_ej230_analysis/logs/exec12b_dispersion.log |   29 +
 results_ej230_analysis/logs/exec13_tN.log          |   24 +
 .../logs/exec14b_asset_audit_after.log             |    1 +
 .../logs/exec14b_asset_audit_before.log            |   28 +
 .../logs/exec14b_asset_parity.log                  |    1 +
 .../logs/exec14b_figure_frame_consistency.log      |    1 +
 .../logs/exec14b_frame_parity.log                  |    1 +
 results_ej230_analysis/logs/exec14b_latexmk.log    |  792 +++++
 .../logs/exec14b_number_audit.log                  |    1 +
 .../logs/exec14b_pdf_preflight.log                 |   10 +
 results_ej230_analysis/logs/exec14b_pdfimages.txt  |  330 +++
 results_ej230_analysis/logs/exec14b_pdftotext.txt  | 1894 ++++++++++++
 .../logs/exec14b_raster_text_audit.log             |    1 +
 .../logs/exec14b_timing_gate.log                   |   13 +
 .../logs/exec14b_timing_mechanism.log              |   12 +
 results_ej230_analysis/logs/exec14b_window_dip.log |    5 +
 .../logs/exec14c_asset_audit.log                   |    1 +
 .../logs/exec14c_asset_parity.log                  |    1 +
 .../logs/exec14c_duplicate_macros.log              |    0
 .../logs/exec14c_figure_frame_consistency.log      |    1 +
 .../logs/exec14c_frame_parity.log                  |    1 +
 .../logs/exec14c_generate_tables.log               |   50 +
 results_ej230_analysis/logs/exec14c_latexmk.log    |  393 +++
 .../logs/exec14c_number_audit.log                  |    1 +
 .../logs/exec14c_raster_text_audit.log             |    1 +
 .../logs/exec14c_rebuild_report.log                |    1 +
 .../logs/exec14c_root_validation.log               |    1 +
 .../logs/exec14c_timing_mechanism.log              |   12 +
 results_ej230_analysis/per_position_exec07.csv     |   32 +
 results_ej230_analysis/poisson_diagnostics.csv     |  125 +
 .../report/exec13_ej230_report_full.aux            |  302 ++
 .../report/exec13_ej230_report_full.log            | 2512 ++++++++++++++++
 .../report/exec13_ej230_report_full.nav            |  279 ++
 .../report/exec13_ej230_report_full.out            |    5 +
 .../report/exec13_ej230_report_full.pdf            |  Bin 0 -> 20772789 bytes
 .../report/exec13_ej230_report_full.snm            |    0
 .../report/exec13_ej230_report_full.tex            | 1251 ++++++++
 .../report/exec13_ej230_report_full.toc            |    4 +
 results_ej230_analysis/report/figure_manifest.csv  |  165 ++
 results_ej230_analysis/report/figure_manifest.md   |  168 ++
 results_ej230_analysis/report/pdf_preflight.json   |  139 +
 results_ej230_analysis/root_validation_exec14b.csv |   39 +
 results_ej230_analysis/summary_exec07.csv          | 3039 ++++++++++++++++++++
 results_ej230_analysis/tables/adaptive_tN.tex      |   49 +
 results_ej230_analysis/tables/attenuation_fit.tex  |    3 +
 .../tables/endtop_endonly_ratio.tex                |   11 +
 .../tables/endtop_endonly_tails.tex                |    3 +
 .../tables/endtop_physical_diagnosis.tex           |   14 +
 results_ej230_analysis/tables/exec14_macros.tex    |  292 ++
 results_ej230_analysis/tables/fano_fit.tex         |    3 +
 .../tables/key_position_-400.tex                   |   14 +
 .../tables/key_position_-650.tex                   |   14 +
 .../tables/key_position_-690.tex                   |   14 +
 results_ej230_analysis/tables/key_position_0.tex   |   14 +
 results_ej230_analysis/tables/key_position_400.tex |   14 +
 results_ej230_analysis/tables/key_position_650.tex |   14 +
 results_ej230_analysis/tables/key_position_690.tex |   14 +
 .../tables/landau_key_positions.tex                |   18 +
 .../tables/numerical_conclusions.tex               |   26 +
 results_ej230_analysis/tables/position_-100.tex    |   14 +
 results_ej230_analysis/tables/position_-150.tex    |   14 +
 results_ej230_analysis/tables/position_-200.tex    |   14 +
 results_ej230_analysis/tables/position_-250.tex    |   14 +
 results_ej230_analysis/tables/position_-300.tex    |   14 +
 results_ej230_analysis/tables/position_-350.tex    |   14 +
 results_ej230_analysis/tables/position_-450.tex    |   14 +
 results_ej230_analysis/tables/position_-50.tex     |   14 +
 results_ej230_analysis/tables/position_-500.tex    |   14 +
 results_ej230_analysis/tables/position_-550.tex    |   14 +
 results_ej230_analysis/tables/position_-600.tex    |   14 +
 results_ej230_analysis/tables/position_-670.tex    |   14 +
 results_ej230_analysis/tables/position_100.tex     |   14 +
 results_ej230_analysis/tables/position_150.tex     |   14 +
 results_ej230_analysis/tables/position_200.tex     |   14 +
 results_ej230_analysis/tables/position_250.tex     |   14 +
 results_ej230_analysis/tables/position_300.tex     |   14 +
 results_ej230_analysis/tables/position_350.tex     |   14 +
 results_ej230_analysis/tables/position_450.tex     |   14 +
 results_ej230_analysis/tables/position_50.tex      |   14 +
 results_ej230_analysis/tables/position_500.tex     |   14 +
 results_ej230_analysis/tables/position_550.tex     |   14 +
 results_ej230_analysis/tables/position_600.tex     |   14 +
 results_ej230_analysis/tables/position_670.tex     |   14 +
 results_ej230_analysis/tables/provenance.tex       |   31 +
 results_ej230_analysis/tables/sum4_timing.tex      |    3 +
 results_ej230_analysis/tables/tN_summary.tex       |   15 +
 .../tables/threshold_rationale.tex                 |    3 +
 .../tables/top_localization_summary.tex            |    5 +
 .../tables/top_timing_comparison.tex               |   13 +
 .../tables/top_timing_definition.tex               |    9 +
 .../tables/top_timing_estimates.tex                |   13 +
 .../tables/validation_summary.tex                  |   11 +
 .../tables/velocity_estimators.tex                 |    8 +
 .../tables/window_dip_mechanism.tex                |    3 +
 .../tables/window_dip_summary.tex                  |    3 +
 results_ej230_analysis/tables/window_dip_test.tex  |   10 +
 results_ej230_analysis/top_localization_gate.csv   |   32 +
 results_ej230_analysis/validation_exec07.txt       |    1 +
 resume_scan_2.sh                                   |    0
 scripts/audit_beamer_assets.py                     |  260 ++
 scripts/audit_beamer_numbers.py                    |  168 ++
 scripts/audit_exec14b_raster_text.py               |   73 +
 scripts/build_report_msi.sh                        |   28 +
 scripts/check_exec14b_asset_parity.py              |   62 +
 scripts/check_exec14b_figure_frame_consistency.py  |  149 +
 scripts/check_exec14b_frame_parity.py              |  102 +
 scripts/diag_exec14d_endtop_ratio.py               |  353 +++
 scripts/endonly_positions.txt                      |   32 +
 scripts/generate_exec14b_tables.py                 |  615 ++++
 scripts/generate_exec14d_nominal_parameters.py     |   44 +
 scripts/make_beamer_ej230.py                       |  161 ++
 scripts/preflight_exec14b_pdf.py                   |   80 +
 scripts/rebuild_exec14b_report.py                  |  372 +++
 scripts/run_analysis_t0minidaq.sh                  |  170 ++
 scripts/run_center.sh                              |   58 +
 scripts/run_center_msi_16t.sh                      |    3 +
 scripts/run_center_t0minidaq_24t.sh                |    3 +
 scripts/run_exec07_scan.sh                         |    4 +-
 scripts/run_exec14b_main_analysis_repair.sh        |   30 +
 scripts/run_exec14b_special_analysis.sh            |   68 +
 scripts/run_exec14b_special_sims.sh                |  123 +
 scripts/run_exec14b_strict_preflight.sh            |   71 +
 scripts/run_msi.sh                                 |    6 +
 scripts/run_scan.sh                                |  152 +
 scripts/run_scan_msi_16t.sh                        |    3 +
 scripts/run_scan_t0minidaq_24t.sh                  |    3 +
 scripts/run_t0minidaq.sh                           |    6 +
 scripts/validate_exec14b_roots.py                  |  163 ++
 src/DetectorConstruction.cc                        |  202 +-
 src/EventAction.cc                                 |    4 +-
 src/Materials.cc                                   |   25 +
 src/PrimaryGeneratorAction.cc                      |   50 +-
 src/RunAction.cc                                   |   15 +-
 src/SiPMSD.cc                                      |   13 +-
 src/SteppingAction.cc                              |    7 +-
 tests/check_endtop_gdml.py                         |   72 -
 ...xport_endtop_gdml.cc => export_endonly_gdml.cc} |   13 +-
 tests/physics_baseline_check.cc                    |   34 +-
 tests/readout_config_check.cc                      |   89 +-
 440 files changed, 30829 insertions(+), 838 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E179](#e179)
<a id="b12"></a>

### 12. feat/ej230-sslg4
Inspected ref `origin/feat/ej230-sslg4`, full SHA `5b93f4cb6b3b0d002a1fd6739cdf770091877005`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 34, behind 115; merged NO; contains 8349041 NO. [E180](#e180) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-12 09:35:03 +0200 / dowiyogo / EXEC_13-230: tN escala-fija (10 ps), dispersión clúster nearest±2, perfiles Top Npe/timing, núcleo-vs-cola (MAD) — datos EJ-230. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal polished bar skin. Reflectivity: 0.98 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E185](#e185) [E186](#e186) [E322](#e322)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **195**. Full paths, line numbers and contents are preserved in [E187](#e187) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E186](#e186)
```text
288:     // Apply branch-specific reflector properties directly to the bar.
289:     auto* reflector = Materials::CreateBarSkinReflector();
290:     auto* barSkin = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
291:     (void)barSkin;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E185](#e185) [E322](#e322)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E323](#e323). Raw matches below are historical finish-test documentation or percentages, not the requested configuration/result marker.
```text
origin/feat/ej230-sslg4:docs/exec14d_adaptive_tN.md:21:| -450 | 20 | 14 | 23.69 | 69.2% | 95.7% | reduced |
origin/feat/ej230-sslg4:results_ej230_analysis/logs/exec14b_pdftotext.txt:1784:69.2%
origin/feat/ej230-sslg4:results_ej230_analysis/tables/adaptive_tN.tex:13:-450 & 4 & 4 & 14 & 69.2\% \\
```

**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E183](#e183)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E183](#e183)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E181](#e181)
```text
5b93f4c EXEC_13-230: tN escala-fija (10 ps), dispersión clúster nearest±2, perfiles Top Npe/timing, núcleo-vs-cola (MAD) — datos EJ-230
3eb0928 feat(t0minidaq): add t0minidaq-24t and MSI-16t script variants
ca2f1c3 fix(tests): replace escape-fraction guard with sipm-entry guard
264d263 fix(optics): replace reflector volumes with bar skin surface
79e701d EXEC_14E/T5: compile and validate final report
5ddc3e0 EXEC_14E/T4: append legible adaptive t_N backup frames
144ed92 EXEC_14E/T3: validate appended adaptive backup frames
153be3c EXEC_14E/T2: move adaptive t_N material to appended backup
ebbb630 EXEC_14E/T1: restore fixed-bin six-panel t_N displays
8a4d338 EXEC_14D/T6: recompilar y validar reporte final
80fd1fd EXEC_14D/T5: pulir coeficiente analítico y errores Npe
2ae390c EXEC_14D/T4: reescribir diagnóstico y restaurar cajas pedagógicas
3333a8f EXEC_14D/T3: aplicar t_N adaptativo y anotaciones
1379d4b EXEC_14D/T2: unificar fitted sigma(t4)
2a6017b EXEC_14D/T1: diagnosticar EndTop contra End-only
cb4b912 EXEC_14C/T5: preflight + render + documentación
0e70266 EXEC_14C/T4: tablas, reconstrucción y compilación estricta del Beamer
0b8be52 EXEC_14C/T3: regenerar EXEC_09 tail (figura/CSV/tablas EJ-230 auténticas)
34e67a5 EXEC_14C/T2: verificar propagación de constantes 0.5/1.5 ns
41b25d9 EXEC_14C/T1: completar End-only x=0,+400 y validar ROOT
533b360 EXEC_14B WIP: estado interrumpido (CODEX), pre-reanudación
4842079 EXEC_14B/T1: audit Beamer assets and identify broken paths
be6743e EXEC_14/T5: hallazgos documentation
f6aa121 EXEC_14/T4: Beamer compiled + frame parity verified
0b26c5f EXEC_14/T3: EJ-230 analysis products — 161 figs, CSVs, tables
c16f679 EXEC_14/T2: analysis route config — exec12b --tau-d arg + EJ-230 audit files
5dad9e3 EXEC_14/T1: data inventory — 31/31 EJ-230 ROOT files validated
192c4fe EXEC_13/T7: RUNBOOK_EJ230.md — per-machine copy-paste instructions
09c7b8c EXEC_13/T6: exec13_ej230_report_full.tex — 1:1 Beamer mirror (119 frames)
e162fde EXEC_13/T5: headless analysis pipeline for t0minidaq
b53d688 EXEC_13/T4: machine-specific scan and center scripts (24t/16t)
7e7e00e EXEC_13/T3: GDML test conditional via ENABLE_GDML_TEST CMake option
4f4f8fe EXEC_13/T2: tee sim output to terminal in run script
4be50d8 EXEC_13/T1: switch scintillator OPSC-101→OPSC-106 (EJ-204→EJ-230)
```

**Files touched on the branch since its merge base with main:** [E182](#e182)
```text
 .gitignore                                         |    7 +
 CMakeLists.txt                                     |   30 +-
 RUNBOOK_EJ230.md                                   |  123 +
 analysis/exec07/common.py                          |    4 +-
 analysis/exec07/exec09_timing_mechanism.py         |    6 +-
 analysis/exec07/exec10_landau_analysis.py          |    2 +
 analysis/exec07/exec10_time_figures.py             |    2 +
 analysis/exec07/exec10_veff_diagnosis.py           |    2 +
 analysis/exec07/exec11_time_arrival.py             |    2 +-
 analysis/exec07/exec12b_tn_dispersion.py           |    8 +
 analysis/exec07/exec13_ej230_report_full.tex       | 1696 +++++++++++
 analysis/exec07_photon_budget.py                   |    8 +-
 analysis/exec13/__init__.py                        |    0
 analysis/exec13/common13.py                        |   71 +
 analysis/exec13/exec13_230_core_resolution.csv     |    7 +
 analysis/exec13/exec13_230_f1_tn_histograms.csv    |  601 ++++
 .../exec13/exec13_230_f2_cluster_dispersion.csv    |  261 ++
 analysis/exec13/exec13_230_f3_npe_profile.csv      |  211 ++
 analysis/exec13/exec13_230_f4_timing_profile.csv   |  421 +++
 analysis/exec13/exec13_230_fixed_scale.py          | 1675 +++++++++++
 analysis/exec13/exec13_230_gates.csv               |   24 +
 analysis/exec13/exec13_230_report.pdf              |  Bin 0 -> 573636 bytes
 analysis/exec13/exec13_230_report.tex              |  233 ++
 analysis/exec13/exec13_230_resolution_macros.tex   |   58 +
 analysis/exec13/exec13_230_tN_table.tex            |   22 +
 analysis/exec13/exec13_tN_analysis.py              |  906 ++++++
 .../figs/exec13_230_f1_t20pe_overlay_abs.pdf       |  Bin 0 -> 22319 bytes
 .../figs/exec13_230_f1_t20pe_overlay_abs.png       |  Bin 0 -> 62320 bytes
 .../figs/exec13_230_f1_t20pe_overlay_norm.pdf      |  Bin 0 -> 23609 bytes
 .../figs/exec13_230_f1_t20pe_overlay_norm.png      |  Bin 0 -> 65066 bytes
 .../exec13/figs/exec13_230_f1_t20pe_panels.pdf     |  Bin 0 -> 28124 bytes
 .../exec13/figs/exec13_230_f1_t20pe_panels.png     |  Bin 0 -> 90318 bytes
 .../exec13/figs/exec13_230_f1_t4pe_overlay_abs.pdf |  Bin 0 -> 21883 bytes
 .../exec13/figs/exec13_230_f1_t4pe_overlay_abs.png |  Bin 0 -> 59425 bytes
 .../figs/exec13_230_f1_t4pe_overlay_norm.pdf       |  Bin 0 -> 23120 bytes
 .../figs/exec13_230_f1_t4pe_overlay_norm.png       |  Bin 0 -> 59538 bytes
 analysis/exec13/figs/exec13_230_f1_t4pe_panels.pdf |  Bin 0 -> 26852 bytes
 analysis/exec13/figs/exec13_230_f1_t4pe_panels.png |  Bin 0 -> 82026 bytes
 .../exec13/figs/exec13_230_f2_cluster_+0mm.pdf     |  Bin 0 -> 33062 bytes
 .../exec13/figs/exec13_230_f2_cluster_+0mm.png     |  Bin 0 -> 108982 bytes
 .../exec13/figs/exec13_230_f2_cluster_-400mm.pdf   |  Bin 0 -> 32975 bytes
 .../exec13/figs/exec13_230_f2_cluster_-400mm.png   |  Bin 0 -> 125892 bytes
 .../exec13/figs/exec13_230_f2_cluster_-690mm.pdf   |  Bin 0 -> 29995 bytes
 .../exec13/figs/exec13_230_f2_cluster_-690mm.png   |  Bin 0 -> 94484 bytes
 analysis/exec13/figs/exec13_230_f3_npe_profile.pdf |  Bin 0 -> 36343 bytes
 analysis/exec13/figs/exec13_230_f3_npe_profile.png |  Bin 0 -> 72583 bytes
 .../exec13/figs/exec13_230_f4_timing_t20pe.pdf     |  Bin 0 -> 29399 bytes
 .../exec13/figs/exec13_230_f4_timing_t20pe.png     |  Bin 0 -> 82023 bytes
 analysis/exec13/figs/exec13_230_f4_timing_t4pe.pdf |  Bin 0 -> 31849 bytes
 analysis/exec13/figs/exec13_230_f4_timing_t4pe.png |  Bin 0 -> 86999 bytes
 analysis/exec13/figs/exec13_230_f5_t4_core.pdf     |  Bin 0 -> 35931 bytes
 analysis/exec13/figs/exec13_230_f5_t4_core.png     |  Bin 0 -> 102846 bytes
 audit/exec10_landau_analysis.md                    |    6 +-
 audit/exec10_veff_diagnosis.md                     |    8 +-
 audit/exec11_arrival.md                            |   44 +-
 docs/EJ230_material_audit.md                       |  235 ++
 docs/exec14_data_inventory.md                      |   56 +
 docs/exec14_fit_failures.md                        |   63 +
 docs/exec14b_beamer_audit.md                       |   52 +
 docs/exec14b_missing_analysis.md                   |   15 +
 docs/exec14c_state.md                              |   84 +
 docs/exec14d_adaptive_tN.md                        |   49 +
 docs/exec14d_endtop_diagnosis.md                   |   40 +
 docs/exec14d_title_exceptions.md                   |   16 +
 docs/exec14e_tN_revert.md                          |   82 +
 include/DetectorConstruction.hh                    |   14 +-
 macros/edge_scan_smoke.mac                         |    2 +-
 macros/endtop_smoke_center.mac                     |    2 +-
 macros/endtop_smoke_edge.mac                       |    2 +-
 macros/exec08b_run_a_window_center.mac             |    2 +-
 macros/exec08b_run_b_window_midpoint.mac           |    2 +-
 macros/exec08b_run_c1_window_mirror.mac            |    2 +-
 macros/exec08b_run_c2_window_exact_mirror.mac      |    2 +-
 results_ej230_analysis/conclusions_exec07.md       |   19 +
 results_ej230_analysis/csv/exec08b_timing_gate.csv |    7 +
 .../csv/exec08b_timing_gate_raw.csv                |   13 +
 .../csv/exec08b_window_dip_profiles.csv            |   56 +
 .../csv/exec09_tail_comparison.csv                 |   16 +
 results_ej230_analysis/csv/exec09_tail_metrics.csv |   31 +
 .../csv/exec09_timing_verdict.txt                  |    6 +
 results_ej230_analysis/csv/exec13_tN_summary.csv   |   94 +
 results_ej230_analysis/csv/exec14d_adaptive_tN.csv |   94 +
 results_ej230_analysis/csv/exec14d_endtop_diag.csv |    7 +
 .../csv/exec14d_nominal_parameters.csv             |    2 +
 .../csv/exec14e_fixed_tN_summary.csv               |  187 ++
 results_ej230_analysis/exec10_fano_by_channel.csv  | 2667 +++++++++++++++++
 results_ej230_analysis/exec10_fano_fit.csv         |    2 +
 results_ej230_analysis/exec10_landau_mpv.csv       |  603 ++++
 .../exec10_late_fraction_representative.csv        |    9 +
 results_ej230_analysis/exec10_velocity_fits.csv    |   10 +
 results_ej230_analysis/exec10_velocity_metrics.csv |   63 +
 results_ej230_analysis/exec11_arrival_metrics.csv  |   94 +
 results_ej230_analysis/figs/P1_npe_vs_x.png        |  Bin 0 -> 233001 bytes
 results_ej230_analysis/figs/P2_npe_heatmap_top.png |  Bin 0 -> 118849 bytes
 results_ej230_analysis/figs/P3_fano_vs_x.png       |  Bin 0 -> 318907 bytes
 .../figs/P4_poisson_check_x-400.png                |  Bin 0 -> 199044 bytes
 .../figs/P4_poisson_check_x-690.png                |  Bin 0 -> 203861 bytes
 .../figs/P4_poisson_check_x0.png                   |  Bin 0 -> 229756 bytes
 .../figs/P4_poisson_check_x400.png                 |  Bin 0 -> 202640 bytes
 .../figs/P4_poisson_check_x690.png                 |  Bin 0 -> 198169 bytes
 results_ej230_analysis/figs/P5_tdist_examples.png  |  Bin 0 -> 255944 bytes
 results_ej230_analysis/figs/P6_tmean_vs_x.png      |  Bin 0 -> 258753 bytes
 results_ej230_analysis/figs/P7_deltaT_end.png      |  Bin 0 -> 270072 bytes
 .../figs/exec08b_id18_impact_maps.png              |  Bin 0 -> 374058 bytes
 .../figs/exec08b_window_dip_profiles.png           |  Bin 0 -> 489245 bytes
 .../figs/exec09_tail_comparison.png                |  Bin 0 -> 340011 bytes
 .../figs/exec10_fano_vs_mean.png                   |  Bin 0 -> 233930 bytes
 .../figs/exec10_landau_fit_example.png             |  Bin 0 -> 100479 bytes
 .../figs/exec10_velocity_estimators.png            |  Bin 0 -> 280181 bytes
 .../figs/exec11_arrival_-100mm.png                 |  Bin 0 -> 204685 bytes
 .../figs/exec11_arrival_-150mm.png                 |  Bin 0 -> 214515 bytes
 .../figs/exec11_arrival_-200mm.png                 |  Bin 0 -> 207191 bytes
 .../figs/exec11_arrival_-250mm.png                 |  Bin 0 -> 211922 bytes
 .../figs/exec11_arrival_-300mm.png                 |  Bin 0 -> 204030 bytes
 .../figs/exec11_arrival_-350mm.png                 |  Bin 0 -> 210514 bytes
 .../figs/exec11_arrival_-400mm.png                 |  Bin 0 -> 207310 bytes
 .../figs/exec11_arrival_-450mm.png                 |  Bin 0 -> 210300 bytes
 .../figs/exec11_arrival_-500mm.png                 |  Bin 0 -> 206720 bytes
 .../figs/exec11_arrival_-50mm.png                  |  Bin 0 -> 208467 bytes
 .../figs/exec11_arrival_-550mm.png                 |  Bin 0 -> 219071 bytes
 .../figs/exec11_arrival_-600mm.png                 |  Bin 0 -> 203660 bytes
 .../figs/exec11_arrival_-650mm.png                 |  Bin 0 -> 216676 bytes
 .../figs/exec11_arrival_-670mm.png                 |  Bin 0 -> 225139 bytes
 .../figs/exec11_arrival_-690mm.png                 |  Bin 0 -> 229147 bytes
 results_ej230_analysis/figs/exec11_arrival_0mm.png |  Bin 0 -> 191246 bytes
 .../figs/exec11_arrival_100mm.png                  |  Bin 0 -> 203459 bytes
 .../figs/exec11_arrival_150mm.png                  |  Bin 0 -> 212372 bytes
 .../figs/exec11_arrival_200mm.png                  |  Bin 0 -> 207820 bytes
 .../figs/exec11_arrival_250mm.png                  |  Bin 0 -> 211887 bytes
 .../figs/exec11_arrival_300mm.png                  |  Bin 0 -> 204235 bytes
 .../figs/exec11_arrival_350mm.png                  |  Bin 0 -> 210990 bytes
 .../figs/exec11_arrival_400mm.png                  |  Bin 0 -> 208226 bytes
 .../figs/exec11_arrival_450mm.png                  |  Bin 0 -> 209545 bytes
 .../figs/exec11_arrival_500mm.png                  |  Bin 0 -> 203660 bytes
 .../figs/exec11_arrival_50mm.png                   |  Bin 0 -> 207820 bytes
 .../figs/exec11_arrival_550mm.png                  |  Bin 0 -> 218482 bytes
 .../figs/exec11_arrival_600mm.png                  |  Bin 0 -> 203211 bytes
 .../figs/exec11_arrival_650mm.png                  |  Bin 0 -> 214913 bytes
 .../figs/exec11_arrival_670mm.png                  |  Bin 0 -> 221144 bytes
 .../figs/exec11_arrival_690mm.png                  |  Bin 0 -> 226340 bytes
 results_ej230_analysis/figs/exec13_tN_-400mm.png   |  Bin 0 -> 227690 bytes
 results_ej230_analysis/figs/exec13_tN_-650mm.png   |  Bin 0 -> 241442 bytes
 results_ej230_analysis/figs/exec13_tN_-690mm.png   |  Bin 0 -> 235239 bytes
 results_ej230_analysis/figs/exec13_tN_0mm.png      |  Bin 0 -> 246137 bytes
 results_ej230_analysis/figs/exec13_tN_400mm.png    |  Bin 0 -> 235145 bytes
 results_ej230_analysis/figs/exec13_tN_650mm.png    |  Bin 0 -> 241645 bytes
 results_ej230_analysis/figs/exec13_tN_690mm.png    |  Bin 0 -> 238605 bytes
 results_ej230_analysis/figs/exec13_tN_summary.png  |  Bin 0 -> 105222 bytes
 .../figs/exec13_tn_-400mm_endL.png                 |  Bin 0 -> 164703 bytes
 .../figs/exec13_tn_-400mm_endR.png                 |  Bin 0 -> 161097 bytes
 .../figs/exec13_tn_-400mm_top.png                  |  Bin 0 -> 212101 bytes
 .../figs/exec13_tn_-650mm_endL.png                 |  Bin 0 -> 168547 bytes
 .../figs/exec13_tn_-650mm_endR.png                 |  Bin 0 -> 131747 bytes
 .../figs/exec13_tn_-650mm_top.png                  |  Bin 0 -> 213667 bytes
 .../figs/exec13_tn_-690mm_endL.png                 |  Bin 0 -> 184723 bytes
 .../figs/exec13_tn_-690mm_endR.png                 |  Bin 0 -> 128295 bytes
 .../figs/exec13_tn_-690mm_top.png                  |  Bin 0 -> 212919 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_endL.png |  Bin 0 -> 153867 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_endR.png |  Bin 0 -> 156515 bytes
 results_ej230_analysis/figs/exec13_tn_0mm_top.png  |  Bin 0 -> 210186 bytes
 .../figs/exec13_tn_400mm_endL.png                  |  Bin 0 -> 156488 bytes
 .../figs/exec13_tn_400mm_endR.png                  |  Bin 0 -> 165674 bytes
 .../figs/exec13_tn_400mm_top.png                   |  Bin 0 -> 211626 bytes
 .../figs/exec13_tn_650mm_endL.png                  |  Bin 0 -> 136647 bytes
 .../figs/exec13_tn_650mm_endR.png                  |  Bin 0 -> 170238 bytes
 .../figs/exec13_tn_650mm_top.png                   |  Bin 0 -> 215272 bytes
 .../figs/exec13_tn_690mm_endL.png                  |  Bin 0 -> 128415 bytes
 .../figs/exec13_tn_690mm_endR.png                  |  Bin 0 -> 185895 bytes
 .../figs/exec13_tn_690mm_top.png                   |  Bin 0 -> 213888 bytes
 .../figs/exec14d_adaptive_tN_-400mm.png            |  Bin 0 -> 161061 bytes
 .../figs/exec14d_adaptive_tN_-650mm.png            |  Bin 0 -> 167773 bytes
 .../figs/exec14d_adaptive_tN_-690mm.png            |  Bin 0 -> 162200 bytes
 .../figs/exec14d_adaptive_tN_0mm.png               |  Bin 0 -> 157373 bytes
 .../figs/exec14d_adaptive_tN_400mm.png             |  Bin 0 -> 168065 bytes
 .../figs/exec14d_adaptive_tN_650mm.png             |  Bin 0 -> 170074 bytes
 .../figs/exec14d_adaptive_tN_690mm.png             |  Bin 0 -> 166546 bytes
 .../figs/exec14d_adaptive_tN_summary.png           |  Bin 0 -> 124988 bytes
 .../figs/fits_attenuation_velocity.png             |  Bin 0 -> 239898 bytes
 .../figs/muon_-100mm_geometry.png                  |  Bin 0 -> 86814 bytes
 .../figs/muon_-100mm_top_profile.png               |  Bin 0 -> 130050 bytes
 .../figs/muon_-150mm_geometry.png                  |  Bin 0 -> 86584 bytes
 .../figs/muon_-150mm_top_profile.png               |  Bin 0 -> 130972 bytes
 .../figs/muon_-200mm_geometry.png                  |  Bin 0 -> 86031 bytes
 .../figs/muon_-200mm_top_profile.png               |  Bin 0 -> 128479 bytes
 .../figs/muon_-250mm_geometry.png                  |  Bin 0 -> 88413 bytes
 .../figs/muon_-250mm_top_profile.png               |  Bin 0 -> 131490 bytes
 .../figs/muon_-300mm_geometry.png                  |  Bin 0 -> 87067 bytes
 .../figs/muon_-300mm_top_profile.png               |  Bin 0 -> 130048 bytes
 .../figs/muon_-350mm_geometry.png                  |  Bin 0 -> 86576 bytes
 .../figs/muon_-350mm_top_profile.png               |  Bin 0 -> 130743 bytes
 .../figs/muon_-400mm_geometry.png                  |  Bin 0 -> 86364 bytes
 .../figs/muon_-400mm_top_profile.png               |  Bin 0 -> 126950 bytes
 .../figs/muon_-450mm_geometry.png                  |  Bin 0 -> 85343 bytes
 .../figs/muon_-450mm_top_profile.png               |  Bin 0 -> 127479 bytes
 .../figs/muon_-500mm_geometry.png                  |  Bin 0 -> 86759 bytes
 .../figs/muon_-500mm_top_profile.png               |  Bin 0 -> 126793 bytes
 .../figs/muon_-50mm_geometry.png                   |  Bin 0 -> 87030 bytes
 .../figs/muon_-50mm_top_profile.png                |  Bin 0 -> 130249 bytes
 .../figs/muon_-550mm_geometry.png                  |  Bin 0 -> 85110 bytes
 .../figs/muon_-550mm_top_profile.png               |  Bin 0 -> 124255 bytes
 .../figs/muon_-600mm_geometry.png                  |  Bin 0 -> 86326 bytes
 .../figs/muon_-600mm_top_profile.png               |  Bin 0 -> 122760 bytes
 .../figs/muon_-650mm_geometry.png                  |  Bin 0 -> 85418 bytes
 .../figs/muon_-650mm_top_profile.png               |  Bin 0 -> 118836 bytes
 .../figs/muon_-670mm_geometry.png                  |  Bin 0 -> 85573 bytes
 .../figs/muon_-670mm_top_profile.png               |  Bin 0 -> 115700 bytes
 .../figs/muon_-690mm_geometry.png                  |  Bin 0 -> 86032 bytes
 .../figs/muon_-690mm_top_profile.png               |  Bin 0 -> 115351 bytes
 results_ej230_analysis/figs/muon_0mm_geometry.png  |  Bin 0 -> 84816 bytes
 .../figs/muon_0mm_top_profile.png                  |  Bin 0 -> 128228 bytes
 .../figs/muon_100mm_geometry.png                   |  Bin 0 -> 85800 bytes
 .../figs/muon_100mm_top_profile.png                |  Bin 0 -> 129881 bytes
 .../figs/muon_150mm_geometry.png                   |  Bin 0 -> 87145 bytes
 .../figs/muon_150mm_top_profile.png                |  Bin 0 -> 132085 bytes
 .../figs/muon_200mm_geometry.png                   |  Bin 0 -> 88515 bytes
 .../figs/muon_200mm_top_profile.png                |  Bin 0 -> 130820 bytes
 .../figs/muon_250mm_geometry.png                   |  Bin 0 -> 88327 bytes
 .../figs/muon_250mm_top_profile.png                |  Bin 0 -> 131353 bytes
 .../figs/muon_300mm_geometry.png                   |  Bin 0 -> 87084 bytes
 .../figs/muon_300mm_top_profile.png                |  Bin 0 -> 129045 bytes
 .../figs/muon_350mm_geometry.png                   |  Bin 0 -> 88801 bytes
 .../figs/muon_350mm_top_profile.png                |  Bin 0 -> 131357 bytes
 .../figs/muon_400mm_geometry.png                   |  Bin 0 -> 86811 bytes
 .../figs/muon_400mm_top_profile.png                |  Bin 0 -> 127528 bytes
 .../figs/muon_450mm_geometry.png                   |  Bin 0 -> 85908 bytes
 .../figs/muon_450mm_top_profile.png                |  Bin 0 -> 127149 bytes
 .../figs/muon_500mm_geometry.png                   |  Bin 0 -> 86543 bytes
 .../figs/muon_500mm_top_profile.png                |  Bin 0 -> 124196 bytes
 results_ej230_analysis/figs/muon_50mm_geometry.png |  Bin 0 -> 86568 bytes
 .../figs/muon_50mm_top_profile.png                 |  Bin 0 -> 129626 bytes
 .../figs/muon_550mm_geometry.png                   |  Bin 0 -> 87039 bytes
 .../figs/muon_550mm_top_profile.png                |  Bin 0 -> 124104 bytes
 .../figs/muon_600mm_geometry.png                   |  Bin 0 -> 88127 bytes
 .../figs/muon_600mm_top_profile.png                |  Bin 0 -> 123025 bytes
 .../figs/muon_650mm_geometry.png                   |  Bin 0 -> 86700 bytes
 .../figs/muon_650mm_top_profile.png                |  Bin 0 -> 118119 bytes
 .../figs/muon_670mm_geometry.png                   |  Bin 0 -> 87021 bytes
 .../figs/muon_670mm_top_profile.png                |  Bin 0 -> 115175 bytes
 .../figs/muon_690mm_geometry.png                   |  Bin 0 -> 87513 bytes
 .../figs/muon_690mm_top_profile.png                |  Bin 0 -> 114202 bytes
 results_ej230_analysis/figs/tn_-400mm_endL.png     |  Bin 0 -> 164703 bytes
 results_ej230_analysis/figs/tn_-400mm_endR.png     |  Bin 0 -> 161097 bytes
 results_ej230_analysis/figs/tn_-400mm_top.png      |  Bin 0 -> 212101 bytes
 results_ej230_analysis/figs/tn_-650mm_endL.png     |  Bin 0 -> 168547 bytes
 results_ej230_analysis/figs/tn_-650mm_endR.png     |  Bin 0 -> 131747 bytes
 results_ej230_analysis/figs/tn_-650mm_top.png      |  Bin 0 -> 213667 bytes
 results_ej230_analysis/figs/tn_-690mm_endL.png     |  Bin 0 -> 184723 bytes
 results_ej230_analysis/figs/tn_-690mm_endR.png     |  Bin 0 -> 128295 bytes
 results_ej230_analysis/figs/tn_-690mm_top.png      |  Bin 0 -> 212919 bytes
 results_ej230_analysis/figs/tn_0mm_endL.png        |  Bin 0 -> 153867 bytes
 results_ej230_analysis/figs/tn_0mm_endR.png        |  Bin 0 -> 156515 bytes
 results_ej230_analysis/figs/tn_0mm_top.png         |  Bin 0 -> 210186 bytes
 results_ej230_analysis/figs/tn_400mm_endL.png      |  Bin 0 -> 156488 bytes
 results_ej230_analysis/figs/tn_400mm_endR.png      |  Bin 0 -> 165674 bytes
 results_ej230_analysis/figs/tn_400mm_top.png       |  Bin 0 -> 211626 bytes
 results_ej230_analysis/figs/tn_650mm_endL.png      |  Bin 0 -> 136647 bytes
 results_ej230_analysis/figs/tn_650mm_endR.png      |  Bin 0 -> 170238 bytes
 results_ej230_analysis/figs/tn_650mm_top.png       |  Bin 0 -> 215272 bytes
 results_ej230_analysis/figs/tn_690mm_endL.png      |  Bin 0 -> 128415 bytes
 results_ej230_analysis/figs/tn_690mm_endR.png      |  Bin 0 -> 185895 bytes
 results_ej230_analysis/figs/tn_690mm_top.png       |  Bin 0 -> 213888 bytes
 results_ej230_analysis/fit_results_exec07.csv      |    3 +
 .../logs/exec07_photon_budget.log                  |   34 +
 results_ej230_analysis/logs/exec10_landau.log      |   42 +
 .../logs/exec10_time_figures.log                   |    1 +
 results_ej230_analysis/logs/exec10_veff.log        |   47 +
 results_ej230_analysis/logs/exec11_arrival.log     |   94 +
 results_ej230_analysis/logs/exec12b_dispersion.log |   29 +
 results_ej230_analysis/logs/exec13_tN.log          |   24 +
 .../logs/exec14b_asset_audit_after.log             |    1 +
 .../logs/exec14b_asset_audit_before.log            |   28 +
 .../logs/exec14b_asset_parity.log                  |    1 +
 .../logs/exec14b_figure_frame_consistency.log      |    1 +
 .../logs/exec14b_frame_parity.log                  |    1 +
 results_ej230_analysis/logs/exec14b_latexmk.log    |  792 +++++
 .../logs/exec14b_number_audit.log                  |    1 +
 .../logs/exec14b_pdf_preflight.log                 |   10 +
 results_ej230_analysis/logs/exec14b_pdfimages.txt  |  330 +++
 results_ej230_analysis/logs/exec14b_pdftotext.txt  | 1894 ++++++++++++
 .../logs/exec14b_raster_text_audit.log             |    1 +
 .../logs/exec14b_timing_gate.log                   |   13 +
 .../logs/exec14b_timing_mechanism.log              |   12 +
 results_ej230_analysis/logs/exec14b_window_dip.log |    5 +
 .../logs/exec14c_asset_audit.log                   |    1 +
 .../logs/exec14c_asset_parity.log                  |    1 +
 .../logs/exec14c_duplicate_macros.log              |    0
 .../logs/exec14c_figure_frame_consistency.log      |    1 +
 .../logs/exec14c_frame_parity.log                  |    1 +
 .../logs/exec14c_generate_tables.log               |   50 +
 results_ej230_analysis/logs/exec14c_latexmk.log    |  393 +++
 .../logs/exec14c_number_audit.log                  |    1 +
 .../logs/exec14c_raster_text_audit.log             |    1 +
 .../logs/exec14c_rebuild_report.log                |    1 +
 .../logs/exec14c_root_validation.log               |    1 +
 .../logs/exec14c_timing_mechanism.log              |   12 +
 results_ej230_analysis/per_position_exec07.csv     |   32 +
 results_ej230_analysis/poisson_diagnostics.csv     |  125 +
 .../report/exec13_ej230_report_full.aux            |  302 ++
 .../report/exec13_ej230_report_full.log            | 2512 ++++++++++++++++
 .../report/exec13_ej230_report_full.nav            |  279 ++
 .../report/exec13_ej230_report_full.out            |    5 +
 .../report/exec13_ej230_report_full.pdf            |  Bin 0 -> 20772789 bytes
 .../report/exec13_ej230_report_full.snm            |    0
 .../report/exec13_ej230_report_full.tex            | 1251 ++++++++
 .../report/exec13_ej230_report_full.toc            |    4 +
 results_ej230_analysis/report/figure_manifest.csv  |  165 ++
 results_ej230_analysis/report/figure_manifest.md   |  168 ++
 results_ej230_analysis/report/pdf_preflight.json   |  139 +
 results_ej230_analysis/root_validation_exec14b.csv |   39 +
 results_ej230_analysis/summary_exec07.csv          | 3039 ++++++++++++++++++++
 results_ej230_analysis/tables/adaptive_tN.tex      |   49 +
 results_ej230_analysis/tables/attenuation_fit.tex  |    3 +
 .../tables/endtop_endonly_ratio.tex                |   11 +
 .../tables/endtop_endonly_tails.tex                |    3 +
 .../tables/endtop_physical_diagnosis.tex           |   14 +
 results_ej230_analysis/tables/exec14_macros.tex    |  292 ++
 results_ej230_analysis/tables/fano_fit.tex         |    3 +
 .../tables/key_position_-400.tex                   |   14 +
 .../tables/key_position_-650.tex                   |   14 +
 .../tables/key_position_-690.tex                   |   14 +
 results_ej230_analysis/tables/key_position_0.tex   |   14 +
 results_ej230_analysis/tables/key_position_400.tex |   14 +
 results_ej230_analysis/tables/key_position_650.tex |   14 +
 results_ej230_analysis/tables/key_position_690.tex |   14 +
 .../tables/landau_key_positions.tex                |   18 +
 .../tables/numerical_conclusions.tex               |   26 +
 results_ej230_analysis/tables/position_-100.tex    |   14 +
 results_ej230_analysis/tables/position_-150.tex    |   14 +
 results_ej230_analysis/tables/position_-200.tex    |   14 +
 results_ej230_analysis/tables/position_-250.tex    |   14 +
 results_ej230_analysis/tables/position_-300.tex    |   14 +
 results_ej230_analysis/tables/position_-350.tex    |   14 +
 results_ej230_analysis/tables/position_-450.tex    |   14 +
 results_ej230_analysis/tables/position_-50.tex     |   14 +
 results_ej230_analysis/tables/position_-500.tex    |   14 +
 results_ej230_analysis/tables/position_-550.tex    |   14 +
 results_ej230_analysis/tables/position_-600.tex    |   14 +
 results_ej230_analysis/tables/position_-670.tex    |   14 +
 results_ej230_analysis/tables/position_100.tex     |   14 +
 results_ej230_analysis/tables/position_150.tex     |   14 +
 results_ej230_analysis/tables/position_200.tex     |   14 +
 results_ej230_analysis/tables/position_250.tex     |   14 +
 results_ej230_analysis/tables/position_300.tex     |   14 +
 results_ej230_analysis/tables/position_350.tex     |   14 +
 results_ej230_analysis/tables/position_450.tex     |   14 +
 results_ej230_analysis/tables/position_50.tex      |   14 +
 results_ej230_analysis/tables/position_500.tex     |   14 +
 results_ej230_analysis/tables/position_550.tex     |   14 +
 results_ej230_analysis/tables/position_600.tex     |   14 +
 results_ej230_analysis/tables/position_670.tex     |   14 +
 results_ej230_analysis/tables/provenance.tex       |   31 +
 results_ej230_analysis/tables/sum4_timing.tex      |    3 +
 results_ej230_analysis/tables/tN_summary.tex       |   15 +
 .../tables/threshold_rationale.tex                 |    3 +
 .../tables/top_localization_summary.tex            |    5 +
 .../tables/top_timing_comparison.tex               |   13 +
 .../tables/top_timing_definition.tex               |    9 +
 .../tables/top_timing_estimates.tex                |   13 +
 .../tables/validation_summary.tex                  |   11 +
 .../tables/velocity_estimators.tex                 |    8 +
 .../tables/window_dip_mechanism.tex                |    3 +
 .../tables/window_dip_summary.tex                  |    3 +
 results_ej230_analysis/tables/window_dip_test.tex  |   10 +
 results_ej230_analysis/top_localization_gate.csv   |   32 +
 results_ej230_analysis/validation_exec07.txt       |    1 +
 scripts/audit_beamer_assets.py                     |  260 ++
 scripts/audit_beamer_numbers.py                    |  168 ++
 scripts/audit_exec14b_raster_text.py               |   73 +
 scripts/build_report_msi.sh                        |   28 +
 scripts/check_exec14b_asset_parity.py              |   62 +
 scripts/check_exec14b_figure_frame_consistency.py  |  149 +
 scripts/check_exec14b_frame_parity.py              |  102 +
 scripts/diag_exec14d_endtop_ratio.py               |  353 +++
 scripts/generate_exec14b_tables.py                 |  615 ++++
 scripts/generate_exec14d_nominal_parameters.py     |   44 +
 scripts/make_beamer_ej230.py                       |  161 ++
 scripts/preflight_exec14b_pdf.py                   |   80 +
 scripts/rebuild_exec14b_report.py                  |  372 +++
 scripts/run_analysis_t0minidaq.sh                  |  170 ++
 scripts/run_center.sh                              |   58 +
 scripts/run_center_msi_16t.sh                      |    3 +
 scripts/run_center_t0minidaq_24t.sh                |    3 +
 scripts/run_exec07_scan.sh                         |    4 +-
 scripts/run_exec14b_main_analysis_repair.sh        |   30 +
 scripts/run_exec14b_special_analysis.sh            |   68 +
 scripts/run_exec14b_special_sims.sh                |  123 +
 scripts/run_exec14b_strict_preflight.sh            |   71 +
 scripts/run_scan.sh                                |  115 +
 scripts/run_scan_msi_16t.sh                        |    3 +
 scripts/run_scan_t0minidaq_24t.sh                  |    3 +
 scripts/validate_exec14b_roots.py                  |  163 ++
 src/DetectorConstruction.cc                        |  113 +-
 tests/check_endtop_balance.py                      |   14 +-
 tests/physics_baseline_check.cc                    |   34 +-
 tests/readout_config_check.cc                      |   53 +-
 395 files changed, 27341 insertions(+), 192 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E188](#e188)
<a id="b13"></a>

### 13. feat/endonly-mylar
Inspected ref `origin/feat/endonly-mylar`, full SHA `fb3749def29716dc84a33fcad53a21086bc96822`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 23, behind 110; merged NO; contains 8349041 NO. [E189](#e189) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-08-14 16:24:33 +0200 / dowiyogo / docs(beamer): update endonly-mylar presentation with corrections. [E032](#e032)
**C1/C2:** active configuration: dielectric_metal ground/polished configurable bar skin. Reflectivity: 0.90 Mylar default; 0.98 fallback; runtime overrides UNVERIFIED. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E194](#e194) [E195](#e195) [E324](#e324)

```text
205: G4OpticalSurface* CreateBarSurface() {
215:     surf->SetType(dielectric_dielectric);
217:     surf->SetFinish(polished);
223: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
227:     surf->SetType(dielectric_metal);
229:     surf->SetFinish(polished);
248:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
294: G4OpticalSurface* CreateBarSkinReflector() {
313:     surf->SetType(dielectric_metal);
315:     surf->SetFinish(polished);
322:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
331:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
337: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
341:     surf->SetType(dielectric_metal);
343:     surf->SetFinish(ground);
347:     const std::vector<G4double> uniformReflectivity = {reflectivity, reflectivity};
352:     mpt->AddProperty("REFLECTIVITY", energy, uniformReflectivity);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
Default Mylar mode and R=0.90, with configurable parameters: [E333](#e333)
```text
111:     G4String           fTopSurface = "mylar";
112:     G4double           fMylarReflectivity = 0.90;
113:     G4double           fMylarSpecularLobe = 1.0;
114:     G4double           fMylarSigmaAlpha = 0.1 * CLHEP::deg;
```

All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **20**. Full paths, line numbers and contents are preserved in [E196](#e196) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E195](#e195)
```text
323:     auto* reflector = fTopSurface == "mylar"
324:         ? Materials::CreateMylarReflector(
325:               fMylarReflectivity, fMylarSpecularLobe, fMylarSigmaAlpha)
326:         : Materials::CreateBarSkinReflector();
327:     fBarSkinSurface = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
328:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E194](#e194) [E324](#e324)
```text
src/Materials.cc
306:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
56: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E325](#e325). Raw matches below are historical finish-test documentation or percentages, not the requested configuration/result marker.
```text
origin/feat/endonly-mylar:EXEC_22b_REPORT.md:1:# EXEC_22b - polishedbackpainted quick test
origin/feat/endonly-mylar:EXEC_22b_REPORT.md:9:- finish nuevo: `polishedbackpainted`
origin/feat/endonly-mylar:analysis/exec22b_quick.py:376:        "finish_tested": "polishedbackpainted",
```

**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E192](#e192)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E192](#e192)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E190](#e190)
```text
fb3749d docs(beamer): update endonly-mylar presentation with corrections
b6874fe analysis(exec28): scan11 weighted/GLS estimator analysis script
fad89b8 analysis(exec22b): backpainted surface test — script and verdict report
fee455f fix(gitignore): add out/ and LaTeX build artifacts
5bb2043 fix(t0minidaq): update scan scripts and CMakeLists for t0minidaq server
6942ee6 fix(optics): replace reflector volumes with bar skin surface
7347db6 audit: add FINAL_PHYSICS_AUDIT.md, C++/OpenMP bootstrap skeleton, group-velocity validator
8b65786 verify: extend checks — range sensitivity, budget closure, 1-R, σ_eq, no 'same attenuation', M2 bootstrap, v_app
20653d8 regen: regenerate End-only Mylar deck from corrected pipeline; all verify_deck checks pass
1ef8ddf deck: fix physics framework, NPE display, photon budget, velocity/timing interpretation
6e128af analysis: constrained M2 fit, AIC/BIC, bootstrap, photon-budget and group-velocity corrections
bf41a16 pre-audit: snapshot of End-only Mylar files before EXEC FINAL corrections
9a97f16 verify: extend checks for three-model fits, datasheet properties, no ±0.0
4ce2a07 deck: rebuild with 4 new frames, λ_bulk vs λ_eff framework, EndTop comparison
c120c51 analysis: fix yield 10000→10400; weighted three-model attenuation fits; EndTop comparison
489d238 feat(presentation): add end-only + Mylar EJ-204 Beamer deck pipeline
3ae135f test: strengthen base parity guard for SiPM coupling and scintillator
7f7fa4f docs: document end-only Mylar workflow and verification
be33a88 test: add end-only guardrails (no top channels, +Y closed, physical ordering, anti-artifact parity)
898c9eb feat(analysis): make end readout pipeline robust to top-channel absence; emit attenuation + sigma_t CSVs
ccc256d feat(scripts): add multi-host run scripts (msi 16t / t0minidaq 24t)
1c8ef21 feat(geom): add topSurface hook (mylar|sipm) with tunable Mylar reflector
5266c33 docs: audit of EJ204 geometry and readout before end-only refactor
```

**Files touched on the branch since its merge base with main:** [E191](#e191)
```text
 .gitignore                                         |   19 +
 CMakeLists.txt                                     |   51 +-
 EXEC_22b_REPORT.md                                 |   72 +
 FINAL_PHYSICS_AUDIT.md                             |  362 +++++
 README_endonly.md                                  |  129 ++
 analysis/bootstrap_attenuation_openmp.cpp          |  456 ++++++
 analysis/endonly_sum4.py                           |  406 ++++++
 analysis/exec22b_quick.py                          |  402 +++++
 analysis/exec28_scan11_weighted.py                 |  930 ++++++++++++
 analysis/validate_group_velocity.py                |  170 +++
 audit_backup/endonly_mylar_before_final/README.md  |   69 +
 .../analysis_endonly_mylar.py                      | 1252 ++++++++++++++++
 .../endonly_mylar_before_final/build_deck.py       |  840 +++++++++++
 .../endonly_mylar_before_final/deck_values.json    |  202 +++
 .../figures/fig_attenuation.txt                    |    1 +
 .../figures/fig_attenuation_twocomp.txt            |    1 +
 .../figures/fig_endtop_comparison.txt              |    1 +
 .../figures/fig_far_end_yield.txt                  |    1 +
 .../figures/fig_fpt_vs_t50.txt                     |    1 +
 .../figures/fig_mean_arrival_time.txt              |    1 +
 .../figures/fig_npe_vs_x.txt                       |    1 +
 .../figures/fig_photon_budget.txt                  |    1 +
 .../figures/fig_sigma_t.txt                        |    1 +
 .../figures/fig_timing_histograms.txt              |    1 +
 .../endonly_mylar_before_final/main.fdb_latexmk    |  290 ++++
 audit_backup/endonly_mylar_before_final/main.fls   | 1463 +++++++++++++++++++
 .../endonly_mylar_before_final/main.synctex.gz     |  Bin 0 -> 52764 bytes
 audit_backup/endonly_mylar_before_final/main.tex   |  620 ++++++++
 .../tables/tab_attenuation_fits.tex                |   10 +
 .../tables/tab_sigma_t.tex                         |    9 +
 .../tables/tab_sim_params.tex                      |   24 +
 .../endonly_mylar_before_final/verify_deck.py      |  259 ++++
 audit_backup/pre_final_changes.patch               |   30 +
 docs/AUDIT_endonly.md                              |  243 ++++
 include/DetectorConstruction.hh                    |   47 +-
 include/Materials.hh                               |    6 +
 macros/endonly_guard.mac                           |   16 +
 macros/endtop_smoke_center.mac                     |    1 +
 macros/endtop_smoke_edge.mac                       |    1 +
 main.cc                                            |   68 +-
 presentation/endonly_mylar/README.md               |   69 +
 .../endonly_mylar/analysis_endonly_mylar.py        | 1532 ++++++++++++++++++++
 presentation/endonly_mylar/build_deck.py           | 1038 +++++++++++++
 presentation/endonly_mylar/deck_values.json        |  442 ++++++
 .../endonly_mylar/figures/fig_attenuation.txt      |    1 +
 .../figures/fig_attenuation_twocomp.txt            |    1 +
 .../figures/fig_endtop_comparison.txt              |    1 +
 .../endonly_mylar/figures/fig_far_end_yield.txt    |    1 +
 .../endonly_mylar/figures/fig_fpt_vs_t50.txt       |    1 +
 .../figures/fig_mean_arrival_time.txt              |    1 +
 .../endonly_mylar/figures/fig_npe_vs_x.txt         |    1 +
 .../endonly_mylar/figures/fig_photon_budget.txt    |    1 +
 presentation/endonly_mylar/figures/fig_sigma_t.txt |    1 +
 .../figures/fig_timing_histograms.txt              |    1 +
 presentation/endonly_mylar/main.tex                |  746 ++++++++++
 presentation/endonly_mylar/main_corrected.tex      |  760 ++++++++++
 .../endonly_mylar/tables/tab_attenuation_fits.tex  |   10 +
 presentation/endonly_mylar/tables/tab_sigma_t.tex  |    9 +
 .../endonly_mylar/tables/tab_sim_params.tex        |   24 +
 presentation/endonly_mylar/verify_deck.py          |  346 +++++
 resume_scan_2.sh                                   |    0
 scripts/endonly_positions.txt                      |   32 +
 scripts/run_exec07_scan.sh                         |    1 +
 scripts/run_msi.sh                                 |    5 +
 scripts/run_scan.sh                                |  164 +++
 scripts/run_t0minidaq.sh                           |    5 +
 src/DetectorConstruction.cc                        |  149 +-
 src/Materials.cc                                   |   25 +
 src/RunAction.cc                                   |    3 +
 tests/check_anti_artifact_parity.py                |  127 ++
 tests/check_endonly_photon_budget.py               |   56 +
 tests/check_endtop_balance.py                      |   14 +-
 tests/check_physical_ordering.py                   |   52 +
 tests/endonly_geometry_check.cc                    |   48 +
 tests/export_endtop_gdml.cc                        |    1 +
 tests/readout_config_check.cc                      |   63 +-
 76 files changed, 14035 insertions(+), 153 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E197](#e197)
<a id="b14"></a>

### 14. feat/endtop-sslg4
Inspected ref `feat/endtop-sslg4`, full SHA `55766876ffa66e44e8a42461029d6276ff14d164`. Classification: **Different pre-Phase-7 reflector implementation**. [E032](#e032)
**B:** ahead 2, behind 103; merged NO; contains 8349041 NO. [E072](#e072) [E223](#e223) [E035](#e035)
Upstream: `origin/feat/endtop-sslg4` (no ahead/behind annotation). Origin counterpart: `55766876ffa66e44e8a42461029d6276ff14d164`. [E030](#e030)
Tip date / author / subject: 2026-08-14 17:31:26 +0200 / dowiyogo / fix(optics): eliminate group-velocity aliasing bug in EXEC_23 air-gap geometry. [E032](#e032)
**C1/C2:** active configuration: Air-gap + dielectric_metal polished outer reflector border. Reflectivity: 0.95 active default. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E198](#e198) [E199](#e199) [E305](#e305)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.95, 0.95};
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 HIT**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **17**. Full paths, line numbers and contents are preserved in [E200](#e200) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E199](#e199)
```text
284:     {
285:         auto* scintAirSurface     = Materials::CreateBarSurface();
286:         auto* airReflectorSurface = Materials::CreateMylarReflector();
287:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E198](#e198) [E305](#e305)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
332:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
354:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E306](#e306). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E075](#e075)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E075](#e075)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E073](#e073)
```text
5576687 fix(optics): eliminate group-velocity aliasing bug in EXEC_23 air-gap geometry
610b189 feat(validation): physics validation scan — 7 pos × 500 events
```

**Files touched on the branch since its merge base with main:** [E074](#e074)
```text
 macros/validation_scan/validation_01_x-600mm.mac |  18 ++
 macros/validation_scan/validation_02_x-400mm.mac |  18 ++
 macros/validation_scan/validation_03_x-200mm.mac |  18 ++
 macros/validation_scan/validation_04_x0mm.mac    |  18 ++
 macros/validation_scan/validation_05_x+200mm.mac |  18 ++
 macros/validation_scan/validation_06_x+400mm.mac |  18 ++
 macros/validation_scan/validation_07_x+600mm.mac |  18 ++
 scripts/analyze_validation.py                    | 211 +++++++++++++++++++++++
 scripts/run_validation_scan.sh                   |  97 +++++++++++
 src/DetectorConstruction.cc                      |  22 ++-
 10 files changed, 451 insertions(+), 5 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E201](#e201)
<a id="b15"></a>

### 15. feature/sipm-electronics-response
Inspected ref `origin/feature/sipm-electronics-response`, full SHA `bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf`. Classification: **Older alternative electronics/geometry; suitability UNVERIFIED**. [E032](#e032)
**B:** ahead 15, behind 163; merged NO; contains 8349041 NO. [E202](#e202) [E223](#e223) [E035](#e035)
Tip date / author / subject: 2026-06-09 14:35:34 +0200 / dowiyogo / Document and ignore generated simulation runs. [E032](#e032)
**C1/C2:** active configuration: Passive Mylar volume; polished dielectric boundaries. Reflectivity: No explicit 0.95/0.98 REFLECTIVITY assignment. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **MISS**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E207](#e207) [E208](#e208) [E326](#e326)

```text
211: G4OpticalSurface* CreateBarSurface() {
217:     surf->SetType(dielectric_dielectric);
219:     surf->SetFinish(polished);
225: G4OpticalSurface* CreateSiPMSurface() {
232:     surf->SetType(dielectric_dielectric);
234:     surf->SetFinish(polished);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **2**. Full paths, line numbers and contents are preserved in [E209](#e209) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E208](#e208)
```text
214:     // ── Segmented wrap — centre Mylar plus configurable 50 mm edge caps ─────
215:     // The 25 µm Mylar layer (n=1.65) acts as a passive reflector:
216:     //   • Photons from bar (n=1.58) → Mylar: small Fresnel reflection, mostly transmitted
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E207](#e207) [E326](#e326)
```text
src/Materials.cc

include/Materials.hh
38: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E327](#e327). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E205](#e205)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E205](#e205)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E203](#e203)
```text
bb4cf6a Document and ignore generated simulation runs
05e79e2 Add trailing newline to src/SiPMSD.cc
3546c52 Fix missing closing brace in SiPMSD
74ff2e1 Fix SiPMSD: remove stray closing brace and tidy EOF newline
d870974 fix(analysis): handle sparse black edge scans
a373b37 docs(slides): add EJ230 timing detector presentation source
26eaf07 fix(analysis): keep fully dead edge-scan positions in summaries
0c95a68 fix(analysis): support current pandas plotting semantics
3164be5 feat(analysis): emulate FastIC+ Sum-of-N channel grouping for timing resolution
974df48 feat(vis): event displays for top-readout midpoint and cross-talk analysis
af32e07 feat(geom): segmented Mylar wrap with configurable edge caps (mylar|air|black)
c53ea53 feat(scan): add 10mm edge scan and dedicated edge analysis
9683346 feat(sipm): checkpoint electronics response and EJ230 analysis updates
e720244 Agregando generador beamer
9f54f7d Add SiPM electronics response analysis
```

**Files touched on the branch since its merge base with main:** [E204](#e204)
```text
 .claude/settings.json                              |  12 +
 .codex                                             |   0
 .gitignore                                         |   1 +
 CMakeLists.txt                                     |  11 +-
 README.md                                          |  79 ++-
 analysis/ResolutionScan_v2.C                       |  89 ++-
 analysis/SiPMRankingScan_RMS.C                     | 365 ++++++++++++
 analysis/SiPMRankingScan_coreSigma.C               | 474 +++++++++++++++
 analysis/SiPMRankingScan_v2.C                      | 613 +++++++++++++++++++
 analysis/TimeMarkScan.C                            | 250 ++++++++
 analysis/{analyze.py => analyze_basic.py}          |   6 +-
 analysis/analyze_dCFD.py                           | 323 ++++++++++
 .../{analyze_dCFD.C => analyze_dCFD_5thPhoton.C}   |  12 +-
 analysis/analyze_dCFD_fraction.C                   | 110 ++++
 analysis/compare_edge_wraps.py                     |  78 +++
 analysis/convert_vis_exports.py                    |  32 +
 analysis/edge_resolution.C                         | 122 ++++
 analysis/edge_resolution.py                        | 308 ++++++++++
 analysis/fpt_vs_n_profile.C                        | 170 ++++++
 analysis/fpt_vs_n_profile_batch.C                  | 137 +++++
 analysis/grouped_resolution.C                      | 205 +++++++
 analysis/grouped_resolution.py                     | 293 +++++++++
 analysis/merge_runs.py                             |   4 +-
 .../{resolution_vs_x.py => resolution_vs_x_FPT.py} |   8 +-
 analysis/resolution_vs_x_dCFD.py                   | 470 +++++++++++++++
 analysis/sipm_waveform_dcfd                        | Bin 0 -> 83320 bytes
 analysis/sipm_waveform_dcfd.cpp                    | 533 +++++++++++++++++
 analysis/sipm_waveform_dcfd.py                     | 592 +++++++++++++++++++
 analysis/topreadout_crosstalk.py                   | 131 ++++
 audit/runs_file_inventory.txt                      | 137 +++++
 fpt_manifest_to_beamer.py                          | 202 +++++++
 fpt_vs_n_profile.C                                 | 170 ++++++
 fpt_vs_n_profile_batch.C                           | 137 +++++
 fpt_vs_n_profile_batch_slides.C                    | 235 ++++++++
 include/DetectorConstruction.hh                    |  35 +-
 include/Materials.hh                               |  15 +-
 include/SiPMSD.hh                                  |  40 +-
 macros/event_display_top_3d.mac                    |  38 ++
 macros/event_display_top_batch.mac                 |  20 +
 macros/event_display_top_lateral.mac               |  39 ++
 macros/event_display_top_midpoint.mac              |  49 ++
 macros/run.mac                                     |  13 +-
 macros/scan.mac                                    |  12 +-
 macros/scan_edge.mac                               |  22 +
 macros/scan_edge_air.mac                           |  18 +
 macros/scan_edge_black.mac                         |  18 +
 macros/scan_edge_neg.mac                           |  21 +
 macros/scan_edge_step.mac                          |   3 +
 macros/scan_step.mac                               |   2 +-
 macros/vis_scan_png.mac                            |  58 ++
 macros/vis_scan_png_lite.mac                       |  47 ++
 macros/vis_scan_png_step.mac                       |   9 +
 main.cc                                            |   2 +
 png_to_beamer.py                                   | 177 ++++++
 presentacion_timing_detector_scan_ej-230.tex       | 657 +++++++++++++++++++++
 src/DetectorConstruction.cc                        | 190 +++++-
 src/Materials.cc                                   | 126 +++-
 src/RunAction.cc                                   |   8 +-
 src/SiPMSD.cc                                      | 253 +++++++-
 src/SteppingAction.cc                              |   9 +-
 60 files changed, 8040 insertions(+), 150 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E210](#e210)
<a id="b16"></a>

### 16. main
Inspected ref `main`, full SHA `8349041140958226a0ac1cb3bb3e30aff2303435`. Classification: **Current specified artifact, with audit gaps**. [E032](#e032)
**B:** ahead 0, behind 0; merged YES; contains 8349041 YES. [E077](#e077) [E223](#e223) [E035](#e035)
Upstream: `origin/main` (no ahead/behind annotation). Origin counterpart: `8349041140958226a0ac1cb3bb3e30aff2303435`. [E030](#e030)
Tip date / author / subject: 2026-09-03 00:37:35 +0200 / rrios / EXEC_25: EXEC_25_REPORT.md — CHECKPOINT 5 (FIX-27/27b final). [E032](#e032)
**C1/C2:** active configuration: Air-gap + polished dielectric_dielectric reflector border. Reflectivity: 0.98 active; historical deck data 0.95. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E211](#e211) [E212](#e212) [E307](#e307)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.98, 0.98};  // Vikuiti 3M ESR spec (R≈0.98, 3M product sheet)
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **165**. Full paths, line numbers and contents are preserved in [E213](#e213) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E212](#e212)
```text
311:     {
312:         auto* scintAirSurface     = Materials::CreateBarSurface();
313:         auto* airReflectorSurface = Materials::CreateBarSkinReflector();
314:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E211](#e211) [E307](#e307)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E309](#e309). No matches outside the bundled external libraries.
Deck `presentations/v6/talk_v6.tex`: **514.9 HIT**, **1300 HIT**. [E244](#e244)

```text
33: \newcommand{\Rvikuiti}{0.95}
35: \newcommand{\NpeG}{701.3}
107:     \textbf{Reflector:} Vikuiti ESR on all non-SiPM surfaces ($R = \Rvikuiti$).\\[5pt]
225:       \item Reflected ($R=\Rvikuiti$ per bounce) back into bar
276:     Sim.\ data ($R=\Rvikuiti$): $\Lambda_\text{refl}^H \approx \NapLambdareflHmmsim$.}\\[2pt]
277:     At $L/2=700\mm$: survival $= \Rvikuiti^{85.6} \approx 1.2\%$.\\[4pt]
389:       \item Air$\to$Mylar: \texttt{CreateBarSkinReflector()} as border surface — \texttt{dielectric\_dielectric}, $R=\Rvikuiti$ constant
441:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
533:     The 3rd-order statistic (out of $\sim\!1300$ pooled photons) suppresses random single-photon fluctuations.\\[6pt]
772:       \item Reflector: \textbf{Vikuiti ESR} ($R=\Rvikuiti$ in data; code $0.98$ since \texttt{c7acb7a})
787:       $N_\text{pe}$/end at $x=0$ & $514.9$ (hybrid, $N_\text{TOP}=20$) \\
821:       \item Vikuiti (R=\Rvikuiti) recovers non-TIR photons (+23\%)
876:       \item $R = \Rvikuiti$ constant (wavelength-independent)
888:     \textbf{A2 — R=\Rvikuiti{} vs R=0.98}\\[4pt]
889:     All simulation data shown: $R=\Rvikuiti$ (Materials.cc at data generation).\\[4pt]
893:     With $R=\Rvikuiti$ (simulation data): $\Lambda_\text{refl}^H = \NapLambdareflHmmsim$ (shorter --- more loss).
915:       \item Vikuiti reflectivity at air–Mylar interface ($R=\Rvikuiti$)
942:     Reflector & \texttt{dielectric\_metal skin} & border surface, R=\Rvikuiti & border surface, R=\Rvikuiti \\
944:     $N_\text{pe}$/end & $\sim\!0.37$ (BUG) & 701.3 & 514.9 \\
967:     \textbf{Fix (exec21-optfix):} Changed to \texttt{dielectric\_dielectric} border surface with $R=\Rvikuiti$, preserving natural TIR at bar–air interface.
991:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
1021:       \textbf{G4/end (END-only)} & \textbf{701.3} & +refl.-rec. \\
1022:       G4/end (hybrid, $N_\text{TOP}{=}20$) & 514.9 & TOP intercept \\
1258:       Reflector & Vikuiti ESR ($R=\Rvikuiti$) \\
```

Exact loose `0.95` occurrences in the deck: **1**. Only the macro definition; no unapplied literal remains in this deck.
```text
33: \newcommand{\Rvikuiti}{0.95}
```

701.3 contextual caution: **SUSPICIOUS** in companion FINAL_NUMBERS.md under its N_TOP=20 default; deck S20 now correctly uses 514.9. END-only-labeled 701.3 occurrences alone are not errors.

```text
32: \newcommand{\LambdaH}{405\mm}
33: \newcommand{\Rvikuiti}{0.95}
34: \newcommand{\NpeNapkin}{570}
35: \newcommand{\NpeG}{701.3}
36: \newcommand{\sigmaENDval}{53.68\ps}
37: \newcommand{\sigmaTOPval}{15.20\ps}
38: \newcommand{\sigmaBLUEval}{15.21\ps}
438:       \midrule
439:       TIR-guided / end & 570 & — \\
440:       Reflector-recovered & 0 & — \\
441:       \textbf{Total / end} & \textbf{570} & \textbf{701.3} \\
442:       \bottomrule
443:     \end{tabular}\\[1pt]
444:     $\Delta N_\text{pe}/N_\text{nap} = +23\%$ (END-only).\\[4pt]
763: \section{Design Summary}
764: % ====================================================================
765: 
766: \begin{frame}{Design Decision: EJ-230 + 20 TOP + 16 END}
767:   \begin{columns}[T]
768:     \column{0.52\linewidth}
769:     \textbf{Selected configuration:}
941:     Air gap & No & Yes (0.10 mm) & Yes (0.10 mm) \\
942:     Reflector & \texttt{dielectric\_metal skin} & border surface, R=\Rvikuiti & border surface, R=\Rvikuiti \\
943:     TIR & Eliminated by skin & Natural Fresnel & Natural Fresnel \\
944:     $N_\text{pe}$/end & $\sim\!0.37$ (BUG) & 701.3 & 514.9 \\
945:     $\sigma_\text{END}$ & N/A (insufficient pe) & \sigmaENDonlyval\ ($m^*=7$) & \sigmaENDval\ ($m=8$) \\
946:     TOP & No & No ($N=0$) & Yes ($N=4,8,14,20$) \\
947:     Status & \textcolor{red}{\textbf{INVALID}} & superseded & \textcolor{green!60!black}{\textbf{current}} \\
988:     Yield & 10000 ph/MeV & 10400 ph/MeV & 9700 ph/MeV \\
989:     $n$ & 1.58 & 1.58 & 1.58 \\
990:     \midrule
991:     $N_\text{pe}$/end (G4, $x=0$) & 1311 & 941 & \textbf{701.3} \\
992:     $N_\text{pe}$/end (napkin) & — & — & 570 \\
993:     $\sigma_t$ (G4, $x{=}0$, $m^*$-opt) & \RessigmaxzEROejA & \RessigmaxzEROejB & \RessigmaxzEROejC \\
994:     Napkin bulk survival & 0.832 & 0.646 & \textbf{0.558} \\
1018:       PDE & 0.400 & assumed \\
1019:       \textbf{Napkin/end} & \textbf{570} & TIR-only \\
1020:       \midrule
1021:       \textbf{G4/end (END-only)} & \textbf{701.3} & +refl.-rec. \\
1022:       G4/end (hybrid, $N_\text{TOP}{=}20$) & 514.9 & TOP intercept \\
1023:       TOP interception & $-26.6\%$ & vs END-only \\
1024:       Surplus (END-only) & $+23\%$ & Pop.\ II \\
```

**C4:** **0/14 complete exact-suffix figure trios.** [E080](#e080)
| Figure stem (relative to deck figs/) | Figure PDF | .root | .csv | .meta.json | _meta.json alternative | Matching sidecars elsewhere in tree |
|---|---|---|---|---|---|---|
| `fig1_kscan` | HIT | MISS | MISS | MISS | MISS | none |
| `figM1_mat_sigma_end` | HIT | MISS | MISS | MISS | MISS | none |
| `figM2_mat_npe` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_bulk_survival` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_end_mscan` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_npe_x` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_ntop_scan_full` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_sigma_t_x` | HIT | MISS | MISS | MISS | MISS | none |
| `fig_survival_RN` | HIT | MISS | MISS | MISS | MISS | none |
| `v5_pareto` | HIT | MISS | MISS | MISS | MISS | none |
| `v5_sigma_vs_x` | HIT | MISS | MISS | MISS | MISS | none |
| `v5_top_position_loo` | HIT | MISS | HIT | MISS | HIT | none |
| `v5_veff_fit` | HIT | MISS | HIT | MISS | HIT | none |
| `v5_veff_residual` | HIT | MISS | MISS | MISS | MISS | none |
Figure list is derived from literal `\anafig{...}` and `\includegraphics{figs/...}` calls in the tracked TeX. No dynamic figure paths were assumed; TeX build execution and contents inside figure PDFs were not tested.
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E078](#e078)
```text
(none)
```

**Files touched on the branch since its merge base with main:** [E079](#e079)
```text
(none)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E214](#e214)
<a id="b17"></a>

### 17. wip/host-stash-endtop-junio
Inspected ref `wip/host-stash-endtop-junio`, full SHA `9710f340f239e9cc40be9ad7d7d360fd701dee76`. Classification: **Pre-Phase-7 target optics**. [E032](#e032)
**B:** ahead 1, behind 136; merged NO; contains 8349041 NO. [E082](#e082) [E223](#e223) [E035](#e035)
Upstream: `origin/wip/host-stash-endtop-junio` (no ahead/behind annotation). Origin counterpart: `9710f340f239e9cc40be9ad7d7d360fd701dee76`. [E030](#e030)
Tip date / author / subject: 2026-08-31 20:32:22 +0200 / rrios / wip: recover June stash (CMakeLists + exec07 scan script). [E032](#e032)
**C1/C2:** active configuration: dielectric_metal polished explicit sibling-panel border. Reflectivity: 0.98 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E215](#e215) [E216](#e216) [E310](#e310)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
293: G4OpticalSurface* CreateBarSkinReflector() {
312:     surf->SetType(dielectric_metal);
314:     surf->SetFinish(polished);
321:     G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};
330:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
```

Active/default reflector assignment markers: **R=0.98 HIT; R=0.95 MISS**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **5**. Full paths, line numbers and contents are preserved in [E217](#e217) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E216](#e216)
```text
255:     // volumes avoid global skin surfaces shadowing active SiPM faces.
256:     auto* reflector = Materials::CreateBarSkinReflector();
257:     const G4double foilHalfT = 0.5 * um;
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E215](#e215) [E310](#e310)
```text
src/Materials.cc
305:     // R=0.95 (aluminized Mylar / high-quality reflector).
include/Materials.hh
50: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E311](#e311). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E085](#e085)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E085](#e085)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E083](#e083)
```text
9710f34 wip: recover June stash (CMakeLists + exec07 scan script)
```

**Files touched on the branch since its merge base with main:** [E084](#e084)
```text
 CMakeLists.txt             | 23 +++++++++++++----------
 scripts/run_exec07_scan.sh |  6 +++---
 2 files changed, 16 insertions(+), 13 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E218](#e218)
<a id="b18"></a>

### 18. wip/host-uncommitted-2026-08-31
Inspected ref `wip/host-uncommitted-2026-08-31`, full SHA `d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90`. Classification: **Historical geometry; unmerged work**. [E032](#e032)
**B:** ahead 2, behind 93; merged NO; contains 8349041 NO. [E087](#e087) [E223](#e223) [E035](#e035)
Upstream: `origin/wip/host-uncommitted-2026-08-31` (no ahead/behind annotation). Origin counterpart: `d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90`. [E030](#e030)
Tip date / author / subject: 2026-08-31 20:33:09 +0200 / rrios / wip: add EndSparseTop run script. [E032](#e032)
**C1/C2:** active configuration: Air-gap + polished dielectric_dielectric reflector border. Reflectivity: 0.95 active. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **HIT**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. [E219](#e219) [E220](#e220) [E312](#e312)

```text
204: G4OpticalSurface* CreateBarSurface() {
214:     surf->SetType(dielectric_dielectric);
216:     surf->SetFinish(polished);
222: G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
226:     surf->SetType(dielectric_metal);
228:     surf->SetFinish(polished);
247:     mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
281: G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
288:     surf->SetType(dielectric_metal);
290:     surf->SetFinish(polished);
294:     const std::vector<G4double> refl   = {reflectivity, reflectivity};
299:     mpt->AddProperty("REFLECTIVITY",        energy, refl);
324: G4OpticalSurface* CreateBarSkinReflector() {
347:     surf->SetType(dielectric_dielectric);
349:     surf->SetFinish(polished);
354:     const std::vector<G4double> refl   = {0.95, 0.95};
357:     mpt->AddProperty("REFLECTIVITY", energy, refl);
```

Active/default reflector assignment markers: **R=0.98 MISS; R=0.95 HIT**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.
All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **19**. Full paths, line numbers and contents are preserved in [E221](#e221) and the complete outputs appendix. No fixed line-number assumption was used.
Geometry call sites: [E220](#e220)
```text
346:     {
347:         auto* scintAirSurface     = Materials::CreateBarSurface();
348:         auto* airReflectorSurface = Materials::CreateBarSkinReflector();
349:
```

Loose `0.95` in project surface/header sources (including historical comments; not automatically active): [E219](#e219) [E312](#e312)
```text
src/Materials.cc
269:     // transmitted; the air→Mylar surface with dielectric_metal + REFLECTIVITY=0.95
332:     //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
354:     const std::vector<G4double> refl   = {0.95, 0.95};
include/Materials.hh
39: // dielectric_metal with R=0.95 to model Mylar substrate reflectance.
40: G4OpticalSurface* CreateMylarReflector(G4double reflectivity = 0.95,
55: // dielectric_metal | groundfrontpainted | R = 0.95
```

**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search [E313](#e313). No matches outside the bundled external libraries.
**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. [E090](#e090)
**C4:** No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass. [E090](#e090)
**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): [E088](#e088)
```text
d2c8d4c wip: add EndSparseTop run script
1f75563 wip: uncommitted DetectorConstruction and SiPMSD changes from host
```

**Files touched on the branch since its merge base with main:** [E089](#e089)
```text
 include/DetectorConstruction.hh   |  30 ++++++++---
 scripts/run_end_vik_sparse_top.sh | 102 ++++++++++++++++++++++++++++++++++++++
 src/DetectorConstruction.cc       |  72 ++++++++++++++++++++++++---
 src/SiPMSD.cc                     |   3 +-
 4 files changed, 192 insertions(+), 15 deletions(-)
```

Commits absent from all 18 advertised origin branch histories: **0**. [E222](#e222)
## 8. Additional main evidence and uncertainties
The tracked `analysis/top_npe_diag/top_npe_diag.csv` and its metadata provide the observable data behind the hybrid correction. This audit read their bytes; it did not rerun the raw ROOT event analysis. [E251](#e251) [E252](#e252)
Main still has historical R=0.95 prose in Materials.cc comments, DetectorConstruction.cc comments, CONFIGURATION_AUDIT.md, FINAL_NUMBERS.md, and README/REVISION_NOTES. The single remaining literal in the current deck is the deliberately defined macro. The report does not equate every historical 0.95 statement to an active 0.95 assignment. All marker hits, including loose literals outside the deck, are preserved in the linked evidence outputs.
`OPEN-07` in the current EXEC_25 report actually labels a superseded PowerPoint artifact and is marked resolved there; it is not the same issue description as the user’s OPEN-07 shorthand. [E250](#e250)
- MSI HEAD, branches, worktrees, dirty files, stash, submodules, LFS and divergence: **UNVERIFIED** because SSH authentication failed. W1/W2 identity and reachability: **UNVERIFIED** within the bounded search. Historical reports describing MSI are not substituted for a live inspection.
- Full repository history outside refs, dangling objects, ignored/generated worktree data and unadvertised origin refs were not comprehensively audited. “Clean” means the requested porcelain status has no tracked/untracked changes; it does not mean there are no ignored datasets.
- LFS executable is unavailable. No HEAD pointer or filter=lfs marker was found in either local clone; actual LFS server/object-store content remains **UNVERIFIED**. No submodule status entries or index gitlinks were found.
- No simulation, ROOT regeneration, TeX build, timing recalculation, or PDF rendering was run. Scientific correctness beyond the requested source/marker/artifact checks is **UNVERIFIED**. The old R=0.95 dataset cannot be certified as an R=0.98 result from code text.
- No valid `69.2 / TOP_SUM4_N1` pairing was found. Incidental values in percentages, CSV numbers, or external material spectra were excluded as scientific evidence.
- Sidecar absence is verified in the committed trees. Existence of ignored copies, MSI-only files, renamed provenance mappings not explicitly encoded in these file names, or reproducible regeneration is **UNVERIFIED**. `_meta.json` is reported separately from the specified `.meta.json`.
- No branch count was forced to 18 local heads. The 18-row table is the observed branch-name union; the seven remote-only rows are explicit.
- Patch equivalence between unmerged histories was not tested. All main-exclusive commits remain preservation candidates.
- Original branch of tag creation and EXEC_N mappings not stated in tag names/messages are **UNVERIFIED**.
- Proposed backup tags do not exist as a result of this audit; no cleanup command was executed.
## Appendix A. Exact command ledger and complete outputs
Every reference E### below identifies the exact recorded command, UTC timestamp and exit status. Complete untruncated stdout/stderr are available in [evidence_outputs.md](branch_audit_20260909/evidence_outputs.md), keyed by the same identifiers, and machine-readable [evidence.jsonl](branch_audit_20260909/evidence.jsonl). `git grep` exit 1 means no matches; missing `git lfs` also exited 1 and is interpreted from stderr. The deliberately attempted `main:src/Materials.hh` path failed (128); the real header was subsequently located and read as `include/Materials.hh`. No missing-path result is treated as source evidence.
Git read commands used environment `GIT_OPTIONAL_LOCKS=0`, `GIT_TERMINAL_PROMPT=0`, `GIT_SSH_COMMAND="ssh -o BatchMode=yes -o ConnectTimeout=10"`. The only authorized repository update was HOST fetch, with `fetch.prune=false`, `fetch.pruneTags=false`, `gc.auto=0`, `maintenance.auto=false`, and `--no-prune`. No checkout, merge, rebase, push, destructive reset, cleanup, branch/tag deletion, or stash mutation was executed. Initial and final porcelain outputs for both local clones are identical. Fetch success is verified by its exit code, not inferred from stale tracking refs.
<a id="e001"></a>

**E001** — 2026-09-09T14:27:34.189521+00:00 — exit `0`

```bash
date --iso-8601=seconds
```
<a id="e002"></a>

**E002** — 2026-09-09T14:27:34.190887+00:00 — exit `0`

```bash
hostname
```
<a id="e003"></a>

**E003** — 2026-09-09T14:27:34.193658+00:00 — exit `0`

```bash
id -un
```
<a id="e004"></a>

**E004** — 2026-09-09T14:27:34.198716+00:00 — exit `0`

```bash
bash -c 'find /home -maxdepth 4 -type d -name .git 2>/dev/null | grep -i ej200'
```
<a id="e005"></a>

**E005** — 2026-09-09T14:27:34.214279+00:00 — exit `1`

```bash
bash -c 'find /home /mnt /media /opt /srv -maxdepth 5 \( -type d -name .git -o -type f -name .git \) -print 2>/dev/null'
```
<a id="e006"></a>

**E006** — 2026-09-09T14:27:34.315481+00:00 — exit `255`

```bash
ssh -p 9022 -o BatchMode=yes -o ConnectTimeout=5 rrios@localhost 'echo OK'
```
<a id="e007"></a>

**E007** — 2026-09-09T14:27:34.318209+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 remote -v
```
<a id="e008"></a>

**E008** — 2026-09-09T14:27:34.323710+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 status --porcelain=v2 -b
```
<a id="e009"></a>

**E009** — 2026-09-09T14:27:34.326013+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 symbolic-ref --short HEAD
```
<a id="e010"></a>

**E010** — 2026-09-09T14:27:34.328845+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show -s '--format=%H|%h|%cI|%aI|%s' HEAD
```
<a id="e011"></a>

**E011** — 2026-09-09T14:27:34.330584+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 worktree list --porcelain
```
<a id="e012"></a>

**E012** — 2026-09-09T14:27:34.334193+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 stash list
```
<a id="e013"></a>

**E013** — 2026-09-09T14:27:34.356409+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 submodule status --recursive
```
<a id="e014"></a>

**E014** — 2026-09-09T14:27:34.358657+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-files --stage
```
<a id="e015"></a>

**E015** — 2026-09-09T14:27:34.362611+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 lfs version
```
<a id="e016"></a>

**E016** — 2026-09-09T14:27:34.365523+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 lfs ls-files
```
<a id="e017"></a>

**E017** — 2026-09-09T14:27:34.367635+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 config --get-regexp '^(fetch\.|remote\..*\.(prune|pruneTags)|maintenance\.|gc\.)'
```
<a id="e018"></a>

**E018** — 2026-09-09T14:27:34.369543+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end remote -v
```
<a id="e019"></a>

**E019** — 2026-09-09T14:27:34.377914+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end status --porcelain=v2 -b
```
<a id="e020"></a>

**E020** — 2026-09-09T14:27:34.379996+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end symbolic-ref --short HEAD
```
<a id="e021"></a>

**E021** — 2026-09-09T14:27:34.383381+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end show -s '--format=%H|%h|%cI|%aI|%s' HEAD
```
<a id="e022"></a>

**E022** — 2026-09-09T14:27:34.385503+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end worktree list --porcelain
```
<a id="e023"></a>

**E023** — 2026-09-09T14:27:34.387531+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end stash list
```
<a id="e024"></a>

**E024** — 2026-09-09T14:27:34.408665+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end submodule status --recursive
```
<a id="e025"></a>

**E025** — 2026-09-09T14:27:34.410735+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end ls-files --stage
```
<a id="e026"></a>

**E026** — 2026-09-09T14:27:34.413816+00:00 — exit `1`

```bash
git -C /home/rrios/ej200_end lfs version
```
<a id="e027"></a>

**E027** — 2026-09-09T14:27:34.416612+00:00 — exit `1`

```bash
git -C /home/rrios/ej200_end lfs ls-files
```
<a id="e028"></a>

**E028** — 2026-09-09T14:27:34.418694+00:00 — exit `1`

```bash
git -C /home/rrios/ej200_end config --get-regexp '^(fetch\.|remote\..*\.(prune|pruneTags)|maintenance\.|gc\.)'
```
<a id="e029"></a>

**E029** — 2026-09-09T14:27:35.671200+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 -c fetch.prune=false -c fetch.pruneTags=false -c gc.auto=0 -c maintenance.auto=false fetch --all --no-prune
```
<a id="e030"></a>

**E030** — 2026-09-09T14:27:36.988171+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-remote --heads origin
```
<a id="e031"></a>

**E031** — 2026-09-09T14:27:38.312408+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-remote --tags origin
```
<a id="e032"></a>

**E032** — 2026-09-09T14:27:38.318356+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref --sort=-committerdate refs/heads refs/remotes '--format=%(refname:short)|%(committerdate:iso8601)|%(objectname)|%(upstream:short)|%(upstream:track)|%(authorname)|%(contents:subject)'
```
<a id="e033"></a>

**E033** — 2026-09-09T14:27:38.325104+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch --merged main
```
<a id="e034"></a>

**E034** — 2026-09-09T14:27:38.329521+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch --no-merged main
```
<a id="e035"></a>

**E035** — 2026-09-09T14:27:38.337293+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains 8349041
```
<a id="e036"></a>

**E036** — 2026-09-09T14:27:38.339601+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref '--format=%(refname:short)' refs/heads
```
<a id="e037"></a>

**E037** — 2026-09-09T14:27:38.343453+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...diag/phase7-delta-2026-08-31
```
<a id="e038"></a>

**E038** — 2026-09-09T14:27:38.347035+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..diag/phase7-delta-2026-08-31
```
<a id="e039"></a>

**E039** — 2026-09-09T14:27:38.358176+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...diag/phase7-delta-2026-08-31
```
<a id="e040"></a>

**E040** — 2026-09-09T14:27:38.361354+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r diag/phase7-delta-2026-08-31
```
<a id="e041"></a>

**E041** — 2026-09-09T14:27:38.432708+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' diag/phase7-delta-2026-08-31 --
```
<a id="e042"></a>

**E042** — 2026-09-09T14:27:38.437416+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...docs/branch-diagnosis-2026-08-31
```
<a id="e043"></a>

**E043** — 2026-09-09T14:27:38.441434+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..docs/branch-diagnosis-2026-08-31
```
<a id="e044"></a>

**E044** — 2026-09-09T14:27:38.445917+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...docs/branch-diagnosis-2026-08-31
```
<a id="e045"></a>

**E045** — 2026-09-09T14:27:38.448759+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r docs/branch-diagnosis-2026-08-31
```
<a id="e046"></a>

**E046** — 2026-09-09T14:27:38.513423+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' docs/branch-diagnosis-2026-08-31 --
```
<a id="e047"></a>

**E047** — 2026-09-09T14:27:38.518518+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...exp/pair-scan-2026-06-11
```
<a id="e048"></a>

**E048** — 2026-09-09T14:27:38.522982+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..exp/pair-scan-2026-06-11
```
<a id="e049"></a>

**E049** — 2026-09-09T14:27:38.650684+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...exp/pair-scan-2026-06-11
```
<a id="e050"></a>

**E050** — 2026-09-09T14:27:38.653918+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r exp/pair-scan-2026-06-11
```
<a id="e051"></a>

**E051** — 2026-09-09T14:27:38.730099+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' exp/pair-scan-2026-06-11 --
```
<a id="e052"></a>

**E052** — 2026-09-09T14:27:38.734694+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...feat/bar-end-vikuiti
```
<a id="e053"></a>

**E053** — 2026-09-09T14:27:38.738431+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..feat/bar-end-vikuiti
```
<a id="e054"></a>

**E054** — 2026-09-09T14:27:38.742621+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...feat/bar-end-vikuiti
```
<a id="e055"></a>

**E055** — 2026-09-09T14:27:38.746371+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r feat/bar-end-vikuiti
```
<a id="e056"></a>

**E056** — 2026-09-09T14:27:38.818240+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' feat/bar-end-vikuiti --
```
<a id="e057"></a>

**E057** — 2026-09-09T14:27:38.823316+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...feat/bar-vikuiti
```
<a id="e058"></a>

**E058** — 2026-09-09T14:27:38.828177+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..feat/bar-vikuiti
```
<a id="e059"></a>

**E059** — 2026-09-09T14:27:38.927286+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...feat/bar-vikuiti
```
<a id="e060"></a>

**E060** — 2026-09-09T14:27:38.930511+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r feat/bar-vikuiti
```
<a id="e061"></a>

**E061** — 2026-09-09T14:27:39.005012+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' feat/bar-vikuiti --
```
<a id="e062"></a>

**E062** — 2026-09-09T14:27:39.009943+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...feat/ej204-bar-tir-only
```
<a id="e063"></a>

**E063** — 2026-09-09T14:27:39.014337+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..feat/ej204-bar-tir-only
```
<a id="e064"></a>

**E064** — 2026-09-09T14:27:39.111626+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...feat/ej204-bar-tir-only
```
<a id="e065"></a>

**E065** — 2026-09-09T14:27:39.114831+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r feat/ej204-bar-tir-only
```
<a id="e066"></a>

**E066** — 2026-09-09T14:27:39.191611+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' feat/ej204-bar-tir-only --
```
<a id="e067"></a>

**E067** — 2026-09-09T14:27:39.196363+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...feat/ej230-bar-tir-only
```
<a id="e068"></a>

**E068** — 2026-09-09T14:27:39.200208+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..feat/ej230-bar-tir-only
```
<a id="e069"></a>

**E069** — 2026-09-09T14:27:39.297293+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...feat/ej230-bar-tir-only
```
<a id="e070"></a>

**E070** — 2026-09-09T14:27:39.300328+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r feat/ej230-bar-tir-only
```
<a id="e071"></a>

**E071** — 2026-09-09T14:27:39.378454+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' feat/ej230-bar-tir-only --
```
<a id="e072"></a>

**E072** — 2026-09-09T14:27:39.382855+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...feat/endtop-sslg4
```
<a id="e073"></a>

**E073** — 2026-09-09T14:27:39.386830+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..feat/endtop-sslg4
```
<a id="e074"></a>

**E074** — 2026-09-09T14:27:39.392511+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...feat/endtop-sslg4
```
<a id="e075"></a>

**E075** — 2026-09-09T14:27:39.395253+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r feat/endtop-sslg4
```
<a id="e076"></a>

**E076** — 2026-09-09T14:27:39.458884+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' feat/endtop-sslg4 --
```
<a id="e077"></a>

**E077** — 2026-09-09T14:27:39.461915+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...main
```
<a id="e078"></a>

**E078** — 2026-09-09T14:27:39.464455+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..main
```
<a id="e079"></a>

**E079** — 2026-09-09T14:27:39.466989+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...main
```
<a id="e080"></a>

**E080** — 2026-09-09T14:27:39.469896+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r main
```
<a id="e081"></a>

**E081** — 2026-09-09T14:27:39.546568+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' main --
```
<a id="e082"></a>

**E082** — 2026-09-09T14:27:39.551293+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...wip/host-stash-endtop-junio
```
<a id="e083"></a>

**E083** — 2026-09-09T14:27:39.555423+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..wip/host-stash-endtop-junio
```
<a id="e084"></a>

**E084** — 2026-09-09T14:27:39.559877+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...wip/host-stash-endtop-junio
```
<a id="e085"></a>

**E085** — 2026-09-09T14:27:39.562537+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r wip/host-stash-endtop-junio
```
<a id="e086"></a>

**E086** — 2026-09-09T14:27:39.569308+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' wip/host-stash-endtop-junio --
```
<a id="e087"></a>

**E087** — 2026-09-09T14:27:39.573547+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...wip/host-uncommitted-2026-08-31
```
<a id="e088"></a>

**E088** — 2026-09-09T14:27:39.577195+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..wip/host-uncommitted-2026-08-31
```
<a id="e089"></a>

**E089** — 2026-09-09T14:27:39.582184+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...wip/host-uncommitted-2026-08-31
```
<a id="e090"></a>

**E090** — 2026-09-09T14:27:39.584956+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r wip/host-uncommitted-2026-08-31
```
<a id="e091"></a>

**E091** — 2026-09-09T14:27:39.648764+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' wip/host-uncommitted-2026-08-31 --
```
<a id="e092"></a>

**E092** — 2026-09-09T14:27:39.652701+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 tag --list --sort=-creatordate '--format=%(creatordate:short) %(refname:short) %(objectname:short) %(subject)'
```
<a id="e093"></a>

**E093** — 2026-09-09T14:27:39.655949+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref refs/tags '--format=%(refname:short)|%(objectname)|%(*objectname)|%(contents)'
```
<a id="e094"></a>

**E094** — 2026-09-09T14:27:39.658837+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 tag --list
```
<a id="e095"></a>

**E095** — 2026-09-09T14:27:39.664750+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-endtop-sslg4-2026-06-10
```
<a id="e096"></a>

**E096** — 2026-09-09T14:27:39.670345+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-exec11-2026-06-11
```
<a id="e097"></a>

**E097** — 2026-09-09T14:27:39.675670+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-exec11b-2026-06-11
```
<a id="e098"></a>

**E098** — 2026-09-09T14:27:39.681306+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-exec12-beamer-2026-06-11
```
<a id="e099"></a>

**E099** — 2026-09-09T14:27:39.686488+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-exec12b-2026-06-11
```
<a id="e100"></a>

**E100** — 2026-09-09T14:27:39.691945+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-pairscan-2026-06-11
```
<a id="e101"></a>

**E101** — 2026-09-09T14:27:39.697410+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains checkpoint/pre-physics-baseline-2026-05-08
```
<a id="e102"></a>

**E102** — 2026-09-09T14:27:39.702764+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains diag-photon-budget-v1
```
<a id="e103"></a>

**E103** — 2026-09-09T14:27:39.707743+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains exec07-09-analysis
```
<a id="e104"></a>

**E104** — 2026-09-09T14:27:39.712692+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains physics-baseline-v1
```
<a id="e105"></a>

**E105** — 2026-09-09T14:27:39.718020+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains pre-exec23-260902
```
<a id="e106"></a>

**E106** — 2026-09-09T14:27:39.723309+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains pre-exec24-260902
```
<a id="e107"></a>

**E107** — 2026-09-09T14:27:39.728726+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains pre-exec25-260902
```
<a id="e108"></a>

**E108** — 2026-09-09T14:27:39.733882+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec21-optfix
```
<a id="e109"></a>

**E109** — 2026-09-09T14:27:39.738735+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec22-endtop-optfix
```
<a id="e110"></a>

**E110** — 2026-09-09T14:27:39.743564+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec23-explicit-airgap
```
<a id="e111"></a>

**E111** — 2026-09-09T14:27:39.748369+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec24-pe-budget-audit
```
<a id="e112"></a>

**E112** — 2026-09-09T14:27:39.753217+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec25-optical-realism-bracket
```
<a id="e113"></a>

**E113** — 2026-09-09T14:27:39.758221+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec26-scint-air-surface-realism
```
<a id="e114"></a>

**E114** — 2026-09-09T14:27:39.763217+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains v-exec27-edge-estimator-coupling
```
<a id="e115"></a>

**E115** — 2026-09-09T14:27:39.766100+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show --stat 8349041
```
<a id="e116"></a>

**E116** — 2026-09-09T14:28:41.270700+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show diag/phase7-delta-2026-08-31:src/Materials.cc
```
<a id="e117"></a>

**E117** — 2026-09-09T14:28:41.273580+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show diag/phase7-delta-2026-08-31:src/DetectorConstruction.cc
```
<a id="e118"></a>

**E118** — 2026-09-09T14:28:41.337423+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' diag/phase7-delta-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e119"></a>

**E119** — 2026-09-09T14:28:41.341428+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline diag/phase7-delta-2026-08-31 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e120"></a>

**E120** — 2026-09-09T14:28:41.344355+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show docs/branch-diagnosis-2026-08-31:src/Materials.cc
```
<a id="e121"></a>

**E121** — 2026-09-09T14:28:41.346714+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show docs/branch-diagnosis-2026-08-31:src/DetectorConstruction.cc
```
<a id="e122"></a>

**E122** — 2026-09-09T14:28:41.409833+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' docs/branch-diagnosis-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e123"></a>

**E123** — 2026-09-09T14:28:41.413888+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline docs/branch-diagnosis-2026-08-31 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e124"></a>

**E124** — 2026-09-09T14:28:41.416543+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show exp/pair-scan-2026-06-11:src/Materials.cc
```
<a id="e125"></a>

**E125** — 2026-09-09T14:28:41.418999+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show exp/pair-scan-2026-06-11:src/DetectorConstruction.cc
```
<a id="e126"></a>

**E126** — 2026-09-09T14:28:41.494255+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' exp/pair-scan-2026-06-11 -- ':(exclude)src/external/**'
```
<a id="e127"></a>

**E127** — 2026-09-09T14:28:41.498033+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline exp/pair-scan-2026-06-11 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e128"></a>

**E128** — 2026-09-09T14:28:41.500587+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-end-vikuiti:src/Materials.cc
```
<a id="e129"></a>

**E129** — 2026-09-09T14:28:41.503523+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-end-vikuiti:src/DetectorConstruction.cc
```
<a id="e130"></a>

**E130** — 2026-09-09T14:28:41.572662+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' feat/bar-end-vikuiti -- ':(exclude)src/external/**'
```
<a id="e131"></a>

**E131** — 2026-09-09T14:28:41.577813+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline feat/bar-end-vikuiti --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e132"></a>

**E132** — 2026-09-09T14:28:41.580256+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-vikuiti:src/Materials.cc
```
<a id="e133"></a>

**E133** — 2026-09-09T14:28:41.583175+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-vikuiti:src/DetectorConstruction.cc
```
<a id="e134"></a>

**E134** — 2026-09-09T14:28:41.658413+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' feat/bar-vikuiti -- ':(exclude)src/external/**'
```
<a id="e135"></a>

**E135** — 2026-09-09T14:28:41.662240+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline feat/bar-vikuiti --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e136"></a>

**E136** — 2026-09-09T14:28:41.664671+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej204-bar-tir-only:src/Materials.cc
```
<a id="e137"></a>

**E137** — 2026-09-09T14:28:41.667570+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej204-bar-tir-only:src/DetectorConstruction.cc
```
<a id="e138"></a>

**E138** — 2026-09-09T14:28:41.742578+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' feat/ej204-bar-tir-only -- ':(exclude)src/external/**'
```
<a id="e139"></a>

**E139** — 2026-09-09T14:28:41.746391+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline feat/ej204-bar-tir-only --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e140"></a>

**E140** — 2026-09-09T14:28:41.750532+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/ej204-event-display-tracks
```
<a id="e141"></a>

**E141** — 2026-09-09T14:28:41.754522+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/ej204-event-display-tracks
```
<a id="e142"></a>

**E142** — 2026-09-09T14:28:41.758999+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/ej204-event-display-tracks
```
<a id="e143"></a>

**E143** — 2026-09-09T14:28:41.761703+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/ej204-event-display-tracks
```
<a id="e144"></a>

**E144** — 2026-09-09T14:28:41.825538+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/ej204-event-display-tracks --
```
<a id="e145"></a>

**E145** — 2026-09-09T14:28:41.828771+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej204-event-display-tracks:src/Materials.cc
```
<a id="e146"></a>

**E146** — 2026-09-09T14:28:41.831605+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej204-event-display-tracks:src/DetectorConstruction.cc
```
<a id="e147"></a>

**E147** — 2026-09-09T14:28:41.894324+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/ej204-event-display-tracks -- ':(exclude)src/external/**'
```
<a id="e148"></a>

**E148** — 2026-09-09T14:28:41.898008+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/ej204-event-display-tracks --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e149"></a>

**E149** — 2026-09-09T14:28:41.902177+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/ej228-cylinder
```
<a id="e150"></a>

**E150** — 2026-09-09T14:28:41.906031+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/ej228-cylinder
```
<a id="e151"></a>

**E151** — 2026-09-09T14:28:41.912271+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/ej228-cylinder
```
<a id="e152"></a>

**E152** — 2026-09-09T14:28:41.915035+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/ej228-cylinder
```
<a id="e153"></a>

**E153** — 2026-09-09T14:28:41.978655+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/ej228-cylinder --
```
<a id="e154"></a>

**E154** — 2026-09-09T14:28:41.982128+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-cylinder:src/Materials.cc
```
<a id="e155"></a>

**E155** — 2026-09-09T14:28:41.984947+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-cylinder:src/DetectorConstruction.cc
```
<a id="e156"></a>

**E156** — 2026-09-09T14:28:42.046938+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/ej228-cylinder -- ':(exclude)src/external/**'
```
<a id="e157"></a>

**E157** — 2026-09-09T14:28:42.050607+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/ej228-cylinder --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e158"></a>

**E158** — 2026-09-09T14:28:42.055083+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/ej228-tir-only
```
<a id="e159"></a>

**E159** — 2026-09-09T14:28:42.059109+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/ej228-tir-only
```
<a id="e160"></a>

**E160** — 2026-09-09T14:28:42.066006+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/ej228-tir-only
```
<a id="e161"></a>

**E161** — 2026-09-09T14:28:42.068848+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/ej228-tir-only
```
<a id="e162"></a>

**E162** — 2026-09-09T14:28:42.131861+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/ej228-tir-only --
```
<a id="e163"></a>

**E163** — 2026-09-09T14:28:42.135080+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-tir-only:src/Materials.cc
```
<a id="e164"></a>

**E164** — 2026-09-09T14:28:42.137916+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-tir-only:src/DetectorConstruction.cc
```
<a id="e165"></a>

**E165** — 2026-09-09T14:28:42.200679+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/ej228-tir-only -- ':(exclude)src/external/**'
```
<a id="e166"></a>

**E166** — 2026-09-09T14:28:42.204438+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/ej228-tir-only --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e167"></a>

**E167** — 2026-09-09T14:28:42.206829+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej230-bar-tir-only:src/Materials.cc
```
<a id="e168"></a>

**E168** — 2026-09-09T14:28:42.209761+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej230-bar-tir-only:src/DetectorConstruction.cc
```
<a id="e169"></a>

**E169** — 2026-09-09T14:28:42.284169+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' feat/ej230-bar-tir-only -- ':(exclude)src/external/**'
```
<a id="e170"></a>

**E170** — 2026-09-09T14:28:42.287973+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline feat/ej230-bar-tir-only --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e171"></a>

**E171** — 2026-09-09T14:28:42.292139+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/ej230-endonly-mylar
```
<a id="e172"></a>

**E172** — 2026-09-09T14:28:42.296378+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/ej230-endonly-mylar
```
<a id="e173"></a>

**E173** — 2026-09-09T14:28:42.471927+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/ej230-endonly-mylar
```
<a id="e174"></a>

**E174** — 2026-09-09T14:28:42.475114+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/ej230-endonly-mylar
```
<a id="e175"></a>

**E175** — 2026-09-09T14:28:42.605843+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/ej230-endonly-mylar --
```
<a id="e176"></a>

**E176** — 2026-09-09T14:28:42.609241+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-endonly-mylar:src/Materials.cc
```
<a id="e177"></a>

**E177** — 2026-09-09T14:28:42.612083+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-endonly-mylar:src/DetectorConstruction.cc
```
<a id="e178"></a>

**E178** — 2026-09-09T14:28:42.737699+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/ej230-endonly-mylar -- ':(exclude)src/external/**'
```
<a id="e179"></a>

**E179** — 2026-09-09T14:28:42.741719+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/ej230-endonly-mylar --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e180"></a>

**E180** — 2026-09-09T14:28:42.745899+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/ej230-sslg4
```
<a id="e181"></a>

**E181** — 2026-09-09T14:28:42.750839+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/ej230-sslg4
```
<a id="e182"></a>

**E182** — 2026-09-09T14:28:42.901871+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/ej230-sslg4
```
<a id="e183"></a>

**E183** — 2026-09-09T14:28:42.905016+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/ej230-sslg4
```
<a id="e184"></a>

**E184** — 2026-09-09T14:28:43.032859+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/ej230-sslg4 --
```
<a id="e185"></a>

**E185** — 2026-09-09T14:28:43.036041+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-sslg4:src/Materials.cc
```
<a id="e186"></a>

**E186** — 2026-09-09T14:28:43.038919+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-sslg4:src/DetectorConstruction.cc
```
<a id="e187"></a>

**E187** — 2026-09-09T14:28:43.165186+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/ej230-sslg4 -- ':(exclude)src/external/**'
```
<a id="e188"></a>

**E188** — 2026-09-09T14:28:43.169242+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/ej230-sslg4 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e189"></a>

**E189** — 2026-09-09T14:28:43.173410+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feat/endonly-mylar
```
<a id="e190"></a>

**E190** — 2026-09-09T14:28:43.177999+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feat/endonly-mylar
```
<a id="e191"></a>

**E191** — 2026-09-09T14:28:43.186360+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feat/endonly-mylar
```
<a id="e192"></a>

**E192** — 2026-09-09T14:28:43.189175+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/endonly-mylar
```
<a id="e193"></a>

**E193** — 2026-09-09T14:28:43.252610+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feat/endonly-mylar --
```
<a id="e194"></a>

**E194** — 2026-09-09T14:28:43.255672+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/endonly-mylar:src/Materials.cc
```
<a id="e195"></a>

**E195** — 2026-09-09T14:28:43.258484+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/endonly-mylar:src/DetectorConstruction.cc
```
<a id="e196"></a>

**E196** — 2026-09-09T14:28:43.319534+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feat/endonly-mylar -- ':(exclude)src/external/**'
```
<a id="e197"></a>

**E197** — 2026-09-09T14:28:43.323461+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/endonly-mylar --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e198"></a>

**E198** — 2026-09-09T14:28:43.326107+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/endtop-sslg4:src/Materials.cc
```
<a id="e199"></a>

**E199** — 2026-09-09T14:28:43.329007+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/endtop-sslg4:src/DetectorConstruction.cc
```
<a id="e200"></a>

**E200** — 2026-09-09T14:28:43.391037+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' feat/endtop-sslg4 -- ':(exclude)src/external/**'
```
<a id="e201"></a>

**E201** — 2026-09-09T14:28:43.394679+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline feat/endtop-sslg4 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e202"></a>

**E202** — 2026-09-09T14:28:43.398678+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count main...origin/feature/sipm-electronics-response
```
<a id="e203"></a>

**E203** — 2026-09-09T14:28:43.402458+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main..origin/feature/sipm-electronics-response
```
<a id="e204"></a>

**E204** — 2026-09-09T14:28:43.410687+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 diff --stat main...origin/feature/sipm-electronics-response
```
<a id="e205"></a>

**E205** — 2026-09-09T14:28:43.413091+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feature/sipm-electronics-response
```
<a id="e206"></a>

**E206** — 2026-09-09T14:28:43.417436+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti' origin/feature/sipm-electronics-response --
```
<a id="e207"></a>

**E207** — 2026-09-09T14:28:43.419546+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feature/sipm-electronics-response:src/Materials.cc
```
<a id="e208"></a>

**E208** — 2026-09-09T14:28:43.421953+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feature/sipm-electronics-response:src/DetectorConstruction.cc
```
<a id="e209"></a>

**E209** — 2026-09-09T14:28:43.425904+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' origin/feature/sipm-electronics-response -- ':(exclude)src/external/**'
```
<a id="e210"></a>

**E210** — 2026-09-09T14:28:43.428939+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feature/sipm-electronics-response --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e211"></a>

**E211** — 2026-09-09T14:28:43.431779+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:src/Materials.cc
```
<a id="e212"></a>

**E212** — 2026-09-09T14:28:43.434350+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:src/DetectorConstruction.cc
```
<a id="e213"></a>

**E213** — 2026-09-09T14:28:43.508626+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' main -- ':(exclude)src/external/**'
```
<a id="e214"></a>

**E214** — 2026-09-09T14:28:43.512798+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline main --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e215"></a>

**E215** — 2026-09-09T14:28:43.515246+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-stash-endtop-junio:src/Materials.cc
```
<a id="e216"></a>

**E216** — 2026-09-09T14:28:43.518112+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-stash-endtop-junio:src/DetectorConstruction.cc
```
<a id="e217"></a>

**E217** — 2026-09-09T14:28:43.523565+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' wip/host-stash-endtop-junio -- ':(exclude)src/external/**'
```
<a id="e218"></a>

**E218** — 2026-09-09T14:28:43.526789+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline wip/host-stash-endtop-junio --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e219"></a>

**E219** — 2026-09-09T14:28:43.529191+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-uncommitted-2026-08-31:src/Materials.cc
```
<a id="e220"></a>

**E220** — 2026-09-09T14:28:43.531876+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-uncommitted-2026-08-31:src/DetectorConstruction.cc
```
<a id="e221"></a>

**E221** — 2026-09-09T14:28:43.594074+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1' wip/host-uncommitted-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e222"></a>

**E222** — 2026-09-09T14:28:43.597886+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline wip/host-uncommitted-2026-08-31 --not a18e863bab76cb54ca96a86c9e871ff3501d434c 2a9645fde563abfd20cfdde810124038e312f85e f19c0933133c6ba25e2c4d968d03361fd804fd90 1bfb82743270dcc3ea9b17be03a6af5704cc6850 219fbe385d43629c523c29b4a8d550783b469f35 09f8b1828889f50881639b6a44a137b4a8b27ef5 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 0006919fe047258228e4eb81e8f96e301458a848 b281aea71845a115290458928c542dad5c24b5ac 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 5b93f4cb6b3b0d002a1fd6739cdf770091877005 fb3749def29716dc84a33fcad53a21086bc96822 55766876ffa66e44e8a42461029d6276ff14d164 bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90
```
<a id="e223"></a>

**E223** — 2026-09-09T14:28:43.603052+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --merged main
```
<a id="e224"></a>

**E224** — 2026-09-09T14:28:43.605910+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref '--format=%(refname)|%(objectname)' refs/heads refs/remotes
```
<a id="e225"></a>

**E225** — 2026-09-09T14:28:43.609578+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 stash show --stat 'stash@{0}'
```
<a id="e226"></a>

**E226** — 2026-09-09T14:28:43.612358+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show -s '--format=%H|%P|%cI|%s' 'stash@{0}'
```
<a id="e227"></a>

**E227** — 2026-09-09T14:28:43.616695+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline 'stash@{0}' --not --remotes=origin
```
<a id="e228"></a>

**E228** — 2026-09-09T14:28:43.619262+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count HEAD...origin/main
```
<a id="e229"></a>

**E229** — 2026-09-09T14:28:43.623273+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count HEAD...fb3749def29716dc84a33fcad53a21086bc96822
```
<a id="e230"></a>

**E230** — 2026-09-09T14:28:43.627455+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count fb3749def29716dc84a33fcad53a21086bc96822...origin/main
```
<a id="e231"></a>

**E231** — 2026-09-09T14:28:44.879570+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-remote --symref origin HEAD
```
<a id="e232"></a>

**E232** — 2026-09-09T14:28:44.883382+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref refs/heads '--format=%(refname:short)|%(objectname)|%(upstream:short)|%(upstream:track)'
```
<a id="e233"></a>

**E233** — 2026-09-09T14:28:44.951482+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'filter=lfs|version https://git-lfs.github.com/spec/v1' HEAD --
```
<a id="e234"></a>

**E234** — 2026-09-09T14:28:44.954040+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r HEAD .gitmodules .gitattributes
```
<a id="e235"></a>

**E235** — 2026-09-09T14:28:44.957058+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end for-each-ref refs/heads '--format=%(refname:short)|%(objectname)|%(upstream:short)|%(upstream:track)'
```
<a id="e236"></a>

**E236** — 2026-09-09T14:28:45.027442+00:00 — exit `1`

```bash
git -C /home/rrios/ej200_end grep -n -I -E 'filter=lfs|version https://git-lfs.github.com/spec/v1' HEAD --
```
<a id="e237"></a>

**E237** — 2026-09-09T14:28:45.030155+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end ls-tree -r HEAD .gitmodules .gitattributes
```
<a id="e238"></a>

**E238** — 2026-09-09T14:28:45.035061+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p main:presentations/v6/talk_v6.pdf | sha256sum'
```
<a id="e239"></a>

**E239** — 2026-09-09T14:28:45.039525+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p main:presentations/v6/talk_v6.pdf | wc -c'
```
<a id="e240"></a>

**E240** — 2026-09-09T14:28:45.044081+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p origin/main:presentations/v6/talk_v6.pdf | sha256sum'
```
<a id="e241"></a>

**E241** — 2026-09-09T14:28:45.048472+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p origin/main:presentations/v6/talk_v6.pdf | wc -c'
```
<a id="e242"></a>

**E242** — 2026-09-09T14:28:45.053078+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p 8349041:presentations/v6/talk_v6.pdf | sha256sum'
```
<a id="e243"></a>

**E243** — 2026-09-09T14:28:45.057719+00:00 — exit `0`

```bash
bash -o pipefail -c 'git -C /home/rrios/ej200 cat-file -p 8349041:presentations/v6/talk_v6.pdf | wc -c'
```
<a id="e244"></a>

**E244** — 2026-09-09T14:28:45.060674+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/talk_v6.tex
```
<a id="e245"></a>

**E245** — 2026-09-09T14:28:45.063391+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/results_macros.tex
```
<a id="e246"></a>

**E246** — 2026-09-09T14:28:45.066063+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/FINAL_NUMBERS.md
```
<a id="e247"></a>

**E247** — 2026-09-09T14:28:45.068543+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/README.md
```
<a id="e248"></a>

**E248** — 2026-09-09T14:28:45.071111+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/EXEC_25_REPORT.md
```
<a id="e249"></a>

**E249** — 2026-09-09T14:28:45.073784+00:00 — exit `128`

```bash
git -C /home/rrios/ej200 show main:src/Materials.hh
```
<a id="e250"></a>

**E250** — 2026-09-09T14:28:45.142763+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E 'TOP_SUM4_N1|69[.,]2|69\.1|sigmaTOP|TOP.*N1|OPEN-07' main -- ':(exclude)src/external/**'
```
<a id="e251"></a>

**E251** — 2026-09-09T14:30:04.838827+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:analysis/top_npe_diag/top_npe_diag.csv
```
<a id="e252"></a>

**E252** — 2026-09-09T14:30:04.841372+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:analysis/top_npe_diag/top_npe_diag_meta.json
```
<a id="e253"></a>

**E253** — 2026-09-09T14:30:04.843648+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/figs/v5_veff_fit_meta.json
```
<a id="e254"></a>

**E254** — 2026-09-09T14:30:04.846698+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/figs/v5_top_position_loo_meta.json
```
<a id="e255"></a>

**E255** — 2026-09-09T14:30:04.849698+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:docs/branch_diagnosis/REFLECTIVITY_CHANGE.md
```
<a id="e256"></a>

**E256** — 2026-09-09T14:30:04.853383+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-tree -r origin/feat/bar-end-vikuiti
```
<a id="e257"></a>

**E257** — 2026-09-09T14:30:04.856014+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/bar-end-vikuiti:src/Materials.cc
```
<a id="e258"></a>

**E258** — 2026-09-09T14:30:04.858654+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-list --left-right --count origin/feat/bar-end-vikuiti...feat/bar-end-vikuiti
```
<a id="e259"></a>

**E259** — 2026-09-09T14:30:04.861647+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline origin/feat/bar-end-vikuiti..feat/bar-end-vikuiti
```
<a id="e260"></a>

**E260** — 2026-09-09T14:30:04.866955+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 branch -a --contains 50cf02e34cd7dddfb847739f6249afd1a600478f
```
<a id="e261"></a>

**E261** — 2026-09-09T14:30:04.871705+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline 84e902c522d8af54ace11e5bb38202bf70e8c94d --not --remotes=origin
```
<a id="e262"></a>

**E262** — 2026-09-09T14:30:04.875098+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 stash show -p 'stash@{0}'
```
<a id="e263"></a>

**E263** — 2026-09-09T14:30:06.636436+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-remote origin
```
<a id="e264"></a>

**E264** — 2026-09-09T14:30:06.642323+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 status --porcelain=v2 -b
```
<a id="e265"></a>

**E265** — 2026-09-09T14:30:06.645087+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref refs/heads refs/tags refs/stash '--format=%(refname)|%(objectname)'
```
<a id="e266"></a>

**E266** — 2026-09-09T14:30:06.649580+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end status --porcelain=v2 -b
```
<a id="e267"></a>

**E267** — 2026-09-09T14:30:06.652810+00:00 — exit `0`

```bash
git -C /home/rrios/ej200_end for-each-ref refs/heads refs/tags refs/stash '--format=%(refname)|%(objectname)'
```
<a id="e268"></a>

**E268** — 2026-09-09T14:30:06.655455+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 tag --list
```
<a id="e269"></a>

**E269** — 2026-09-09T14:30:06.657772+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-endtop-sslg4-2026-06-10^{commit}'
```
<a id="e270"></a>

**E270** — 2026-09-09T14:30:06.660066+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-exec11-2026-06-11^{commit}'
```
<a id="e271"></a>

**E271** — 2026-09-09T14:30:06.662418+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-exec11b-2026-06-11^{commit}'
```
<a id="e272"></a>

**E272** — 2026-09-09T14:30:06.664688+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-exec12-beamer-2026-06-11^{commit}'
```
<a id="e273"></a>

**E273** — 2026-09-09T14:30:06.666821+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-exec12b-2026-06-11^{commit}'
```
<a id="e274"></a>

**E274** — 2026-09-09T14:30:06.669613+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-pairscan-2026-06-11^{commit}'
```
<a id="e275"></a>

**E275** — 2026-09-09T14:30:06.671958+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'checkpoint/pre-physics-baseline-2026-05-08^{commit}'
```
<a id="e276"></a>

**E276** — 2026-09-09T14:30:06.674242+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'diag-photon-budget-v1^{commit}'
```
<a id="e277"></a>

**E277** — 2026-09-09T14:30:06.676314+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'exec07-09-analysis^{commit}'
```
<a id="e278"></a>

**E278** — 2026-09-09T14:30:06.678856+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'physics-baseline-v1^{commit}'
```
<a id="e279"></a>

**E279** — 2026-09-09T14:30:06.681143+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'pre-exec23-260902^{commit}'
```
<a id="e280"></a>

**E280** — 2026-09-09T14:30:06.683481+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'pre-exec24-260902^{commit}'
```
<a id="e281"></a>

**E281** — 2026-09-09T14:30:06.686264+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'pre-exec25-260902^{commit}'
```
<a id="e282"></a>

**E282** — 2026-09-09T14:30:06.688549+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec21-optfix^{commit}'
```
<a id="e283"></a>

**E283** — 2026-09-09T14:30:06.690735+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec22-endtop-optfix^{commit}'
```
<a id="e284"></a>

**E284** — 2026-09-09T14:30:06.692872+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec23-explicit-airgap^{commit}'
```
<a id="e285"></a>

**E285** — 2026-09-09T14:30:06.695575+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec24-pe-budget-audit^{commit}'
```
<a id="e286"></a>

**E286** — 2026-09-09T14:30:06.698245+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec25-optical-realism-bracket^{commit}'
```
<a id="e287"></a>

**E287** — 2026-09-09T14:30:06.700593+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec26-scint-air-surface-realism^{commit}'
```
<a id="e288"></a>

**E288** — 2026-09-09T14:30:06.703204+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 rev-parse 'v-exec27-edge-estimator-coupling^{commit}'
```
<a id="e289"></a>

**E289** — 2026-09-09T14:30:20.735398+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show diag/phase7-delta-2026-08-31:analysis/presentation_v6/talk_v6.tex
```
<a id="e290"></a>

**E290** — 2026-09-09T14:30:20.738613+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show diag/phase7-delta-2026-08-31:include/Materials.hh
```
<a id="e291"></a>

**E291** — 2026-09-09T14:30:20.805443+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' diag/phase7-delta-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e292"></a>

**E292** — 2026-09-09T14:30:20.808848+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show docs/branch-diagnosis-2026-08-31:include/Materials.hh
```
<a id="e293"></a>

**E293** — 2026-09-09T14:30:20.874144+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' docs/branch-diagnosis-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e294"></a>

**E294** — 2026-09-09T14:30:20.877574+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show exp/pair-scan-2026-06-11:include/Materials.hh
```
<a id="e295"></a>

**E295** — 2026-09-09T14:30:20.955104+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' exp/pair-scan-2026-06-11 -- ':(exclude)src/external/**'
```
<a id="e296"></a>

**E296** — 2026-09-09T14:30:20.958681+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-end-vikuiti:include/Materials.hh
```
<a id="e297"></a>

**E297** — 2026-09-09T14:30:20.961910+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-end-vikuiti:presentations/v6/talk_v6.tex
```
<a id="e298"></a>

**E298** — 2026-09-09T14:30:21.040559+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' feat/bar-end-vikuiti -- ':(exclude)src/external/**'
```
<a id="e299"></a>

**E299** — 2026-09-09T14:30:21.043832+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/bar-vikuiti:include/Materials.hh
```
<a id="e300"></a>

**E300** — 2026-09-09T14:30:21.119746+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' feat/bar-vikuiti -- ':(exclude)src/external/**'
```
<a id="e301"></a>

**E301** — 2026-09-09T14:30:21.123173+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej204-bar-tir-only:include/Materials.hh
```
<a id="e302"></a>

**E302** — 2026-09-09T14:30:21.200485+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' feat/ej204-bar-tir-only -- ':(exclude)src/external/**'
```
<a id="e303"></a>

**E303** — 2026-09-09T14:30:21.203852+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/ej230-bar-tir-only:include/Materials.hh
```
<a id="e304"></a>

**E304** — 2026-09-09T14:30:21.281904+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' feat/ej230-bar-tir-only -- ':(exclude)src/external/**'
```
<a id="e305"></a>

**E305** — 2026-09-09T14:30:21.285337+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show feat/endtop-sslg4:include/Materials.hh
```
<a id="e306"></a>

**E306** — 2026-09-09T14:30:21.349236+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' feat/endtop-sslg4 -- ':(exclude)src/external/**'
```
<a id="e307"></a>

**E307** — 2026-09-09T14:30:21.352599+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:include/Materials.hh
```
<a id="e308"></a>

**E308** — 2026-09-09T14:30:21.355851+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show main:presentations/v6/talk_v6.tex
```
<a id="e309"></a>

**E309** — 2026-09-09T14:30:21.437073+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' main -- ':(exclude)src/external/**'
```
<a id="e310"></a>

**E310** — 2026-09-09T14:30:21.440428+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-stash-endtop-junio:include/Materials.hh
```
<a id="e311"></a>

**E311** — 2026-09-09T14:30:21.447406+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' wip/host-stash-endtop-junio -- ':(exclude)src/external/**'
```
<a id="e312"></a>

**E312** — 2026-09-09T14:30:21.450619+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show wip/host-uncommitted-2026-08-31:include/Materials.hh
```
<a id="e313"></a>

**E313** — 2026-09-09T14:30:21.515319+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' wip/host-uncommitted-2026-08-31 -- ':(exclude)src/external/**'
```
<a id="e314"></a>

**E314** — 2026-09-09T14:30:21.518660+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej204-event-display-tracks:include/Materials.hh
```
<a id="e315"></a>

**E315** — 2026-09-09T14:30:21.581238+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/ej204-event-display-tracks -- ':(exclude)src/external/**'
```
<a id="e316"></a>

**E316** — 2026-09-09T14:30:21.584628+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-cylinder:include/Materials.hh
```
<a id="e317"></a>

**E317** — 2026-09-09T14:30:21.649330+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/ej228-cylinder -- ':(exclude)src/external/**'
```
<a id="e318"></a>

**E318** — 2026-09-09T14:30:21.652772+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej228-tir-only:include/Materials.hh
```
<a id="e319"></a>

**E319** — 2026-09-09T14:30:21.718158+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/ej228-tir-only -- ':(exclude)src/external/**'
```
<a id="e320"></a>

**E320** — 2026-09-09T14:30:21.721436+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-endonly-mylar:include/Materials.hh
```
<a id="e321"></a>

**E321** — 2026-09-09T14:30:21.851865+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/ej230-endonly-mylar -- ':(exclude)src/external/**'
```
<a id="e322"></a>

**E322** — 2026-09-09T14:30:21.855618+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/ej230-sslg4:include/Materials.hh
```
<a id="e323"></a>

**E323** — 2026-09-09T14:30:21.986287+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/ej230-sslg4 -- ':(exclude)src/external/**'
```
<a id="e324"></a>

**E324** — 2026-09-09T14:30:21.989779+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/endonly-mylar:include/Materials.hh
```
<a id="e325"></a>

**E325** — 2026-09-09T14:30:22.057517+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/endonly-mylar -- ':(exclude)src/external/**'
```
<a id="e326"></a>

**E326** — 2026-09-09T14:30:22.060731+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feature/sipm-electronics-response:include/Materials.hh
```
<a id="e327"></a>

**E327** — 2026-09-09T14:30:22.067951+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feature/sipm-electronics-response -- ':(exclude)src/external/**'
```
<a id="e328"></a>

**E328** — 2026-09-09T14:30:22.071250+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/bar-end-vikuiti:analysis/presentation_v6/talk_v6.tex
```
<a id="e329"></a>

**E329** — 2026-09-09T14:30:22.074195+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/bar-end-vikuiti:include/Materials.hh
```
<a id="e330"></a>

**E330** — 2026-09-09T14:30:22.141596+00:00 — exit `1`

```bash
git -C /home/rrios/ej200 grep -n -I -E '(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted' origin/feat/bar-end-vikuiti -- ':(exclude)src/external/**'
```
<a id="e331"></a>

**E331** — 2026-09-09T14:37:26.382316+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 ls-remote origin
```
<a id="e332"></a>

**E332** — 2026-09-09T14:37:26.387613+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 log --oneline 'stash@{0}' --not 0006919fe047258228e4eb81e8f96e301458a848 04a80473a2ed9aec206af9e0eb1c70a7e25f84d1 09f8b1828889f50881639b6a44a137b4a8b27ef5 0ca685506faa4d309a7dc9ab13b18c85e127b807 1bfb82743270dcc3ea9b17be03a6af5704cc6850 1c5a4c5142efa555f49ba9db1c8e9a2eeb411035 219fbe385d43629c523c29b4a8d550783b469f35 26f4767a3643b5ab4faffccd9cfe23ac1b3042a9 2a9645fde563abfd20cfdde810124038e312f85e 2a9d57b454cfb2e659c8551b6c8d12c1b2a34b2b 3148da19a3d42b300f71d87d32f55d2064ba1b93 47a9a4fb31118ecb6f68cb658d40ec282c6da99f 4e469599670608f7a5b74e56ea0289187b176356 55766876ffa66e44e8a42461029d6276ff14d164 5783a0d2914b01fac6f9f3f48cdb02879384aa86 5b93f4cb6b3b0d002a1fd6739cdf770091877005 66b674ce94c443a9e14c1c8b64a4df9f07ad5e6e 6e8d3d6ed22073185ce57ca14ed34be1c56af463 759669750e8d7396c6d1e4afc2e1922fc119e682 76b582c5cf982706d2b72d3b82a3f879696957f8 7c553263bbe40302b4079cf40bcbac74f72e5544 8349041140958226a0ac1cb3bb3e30aff2303435 9710f340f239e9cc40be9ad7d7d360fd701dee76 9b8361f189712fd8dba332a47e2ed685f194f1c1 a18e863bab76cb54ca96a86c9e871ff3501d434c b281aea71845a115290458928c542dad5c24b5ac bb4cf6a1f80d9c14a3306d98f0274213f9f8a7cf bd2495183420a8097c9742069a0ef556cfc831a7 bd7821170c1c6997a733d0f279f3eac85821c10c c9e61afd8a85182ac8ec23abf8f0e020b1972930 d2c8d4cb69e7067e2ae18f9dd498c9d603bb6c90 e20b47d492e11f099aefb1074b1d8bd26cf8a88a e9245bb9f2b728a2389d915214be988a17f82311 f19c0933133c6ba25e2c4d968d03361fd804fd90 f767be8709a424b042361014fd016673f00d6b69 fb3749def29716dc84a33fcad53a21086bc96822
```
<a id="e333"></a>

**E333** — 2026-09-09T14:37:26.390277+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 show origin/feat/endonly-mylar:include/DetectorConstruction.hh
```
<a id="e334"></a>

**E334** — 2026-09-09T14:37:26.394122+00:00 — exit `0`

```bash
git -C /home/rrios/ej200 for-each-ref --sort=-committerdate refs/heads refs/remotes '--format=%(refname:short)|%(committerdate:iso8601)|%(objectname:short)|%(upstream:short)|%(upstream:track)|%(authorname)|%(contents:subject)'
```
## Appendix B. Method and preservation
All Git commands specified `git -C <clone-path>`. Branch content was read with git show, git grep, git cat-file and git ls-tree; no branch checkout was used. The report, collection scripts, derived table and evidence were written exclusively under `/home/rrios` outside ej200 and ej200_end. Discovery and prerequisite instruction reads were filesystem reads. The first discovery/probe was repeated in the timestamped ledger for reproducibility.

Figure counts, source line numbers, branch counts and group membership are deterministic aggregations of the recorded command outputs. Line numbers in source excerpts are enumerated from the exact git-show output; discovery used patterns rather than assumed line positions. Figures use tree membership checks; no renderer or uncommitted worktree state was used to certify a sidecar. The table-generation implementation is available at [make_report.py](branch_audit_20260909/make_report.py).
