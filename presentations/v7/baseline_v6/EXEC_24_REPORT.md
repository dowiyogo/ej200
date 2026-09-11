# EXEC_24_REPORT.md — Internal Consistency Campaign for talk_v6.tex

Generated: 2026-09-02
Authoritative clone: t0minidaq:/home/rrios/ej200, branch main
Rollback tag: pre-exec24-260902 = 09dc9e2
HEAD before campaign: 09dc9e2 (37 pages)

---

## 1. FIX-11..15 Status

| FIX | Status | Files | Before → After |
|-----|--------|-------|----------------|
| FIX-11 | DONE | talk_v6.tex | Combined G1/G2 [shrink=5] frame (p. 32, truncated) → two separate frames G1 (p. 32) and G2 (p. 33). 17 preamble macros added. G1 formula literals 15.20², 53.68², 196.3, 0.013 → \sigmaTOPbval², \sigmaENDbval², \CETval, \wEND. G2 tabular: narrow p{} columns → lrrr full-width; qualitative BLUE check: analytic w_E^analytic values → empirical w* (0.065, 0.080, 0.125) with position-specific macros. |
| FIX-12 | DONE | talk_v6.tex | S5 (p. 5) items: "Path to END ≈ L/(2cosθ_c)" and "Fast: arrive ~3.7 ns" → "\NapLhalfmm to \NappathTIRlimmm" and "Arrive \NaptaxialENDns to \NaptTIRlimns after emission". Same macros as S11/S13. |
| FIX-13 | DONE | talk_v6.tex | C1 footnote (p. 27): "(see FIX-08 on TOP interception)" → "(see App.~B1)". |
| FIX-14 | DONE | talk_v6.tex | Design Decision (p. 20): "$R=0.95$ constant in code" → "$R=\Rvikuiti$ in data; code $0.98$ since \texttt{c7acb7a}". |
| FIX-15 | DONE | talk_v6.tex | S12 (p. 12): "τ_prop^END ≈ 3.7 ns" → "τ_prop^END ∈ [\NaptaxialENDns, \NaptTIRlimns]". |

---

## 2. Overflow History

### Baseline for EXEC_24 (09dc9e2, 37 pages)

| Log line | Overflow | Size | Frame |
|----------|----------|------|-------|
| — | — | Clean (EXEC_23 result) | — |

### After FIX-11..15 (first compilation attempt)

| Log line | Overflow | Size | Cause |
|----------|----------|------|-------|
| 259 | Overfull \vbox | 21.47 pt | FIX-12: longer bullets in "Two Photon Populations" frame |
| 808 | Overfull \vbox | 2.70 pt | FIX-14: longer reflector bullet in "Design Decision" frame |
| 1163 | Overfull \hbox | 30.53 pt | FIX-11: narrow right column in G2 frame (phantom alignment + long inline math) |

### After corrections

- FIX-12 bullets shortened: "Path to END: $\NapLhalfmm$ to $\NappathTIRlimmm$" (removed "(angle-dependent)"); arrival unchanged.
- FIX-14 reflector bullet shortened to one line: "$R=\Rvikuiti$ in data; code $0.98$ since \texttt{c7acb7a}".
- G2 frame restructured from two-column to single-column (full-width tabular, qualitative check as flowing paragraph — no phantom alignment).

### Final (38 pages) — CLEAN

| Log | Type | Status |
|-----|------|--------|
| Overfull \vbox | — | NONE |
| Overfull \hbox | — | NONE |
| Underfull \hbox line 953 | badness 10000 | Pre-existing; B1 frame footnote paragraph end; not caused by EXEC_24. |

---

## 3. Visual Inspection

All inspected at 80 dpi via pdftoppm.

| Page | Frame | Verdict |
|------|-------|---------|
| 5 | Two Photon Populations Drive the Timing | text complete, no clipping, items on single lines |
| 12 | TOP Estimator: k-th Order Statistic | text complete, τ_prop^END ∈ [3.689 ns, 5.829 ns] visible |
| 20 | Design Decision: EJ-230 + 20 TOP + 16 END | text complete, R=0.95 in data note visible, no clipping |
| 27 | C1 — Full Material Properties | text complete, "see App. B1" visible in footnote |
| 32 | G1 — BLUE Weight Derivation (new) | text complete, no clipping, ρ distinction column visible |
| 33 | G2 — Covariance Convention Comparison (new) | text complete, full-width table, complete qualitative BLUE check with w* values |

Page count: 37 (baseline) → 38 (EXEC_24, G1/G2 split adds one page). Pages 34–38 shift from baseline pages 33–37.

---

## 4. New Preamble Macros (FIX-11)

All added after line 43 of talk_v6.tex in the "convenience macros" block:

| Macro | Value | Source |
|-------|-------|--------|
| `\sigmaENDbval` | 53.68 | bare σ_E (ps); hybrid m=8 x=0; matches \sigmaENDval |
| `\sigmaTOPbval` | 15.20 | bare σ_T (ps); k=3 x=0; matches \sigmaTOPval |
| `\CETval` | 196.3 | C_ET best_est and empirical (ps²) |
| `\CETrobust` | 80.5 | C_ET robust IQR polarization (ps²) |
| `\wENDrobust` | 0.051 | w_E robust IQR; not canonical |
| `\wENDemp` | 0.086 | w_E empirical np.var+np.cov; not canonical |
| `\sigBLUEbestest` | 15.21 | σ_IQR^ev best_est (ps, bare) |
| `\sigBLUErobust` | 15.12 | σ_IQR^ev robust IQR (ps, bare) |
| `\sigBLUEemp` | 15.03 | σ_IQR^ev empirical (ps, bare) |
| `\sigmaformBLUE` | 15.18 | σ^formula Gaussian BLUE algebra (ps, bare) |
| `\wstarxzero` | 0.080 | w*(x=0); source: blue_wscan_x0.meta.json |
| `\wstarxptwo` | 0.125 | w*(x=+200); source: ρ(x) check |
| `\wstarxpsix` | 0.065 | w*(x=+600); source: ρ(x) check |
| `\sigmaTOPxptwo` | 24.3 | bare σ_TOP (ps) at x=+200 |
| `\sigmaTOPxpsix` | 16.8 | bare σ_TOP (ps) at x=+600 |

---

## 5. FIX-13 exec2 Audit

Grep for FIX-, EXEC, exec2 in non-comment content of talk_v6.tex (post-EXEC_24):

| Line | String | Context | Action |
|------|--------|---------|--------|
| (none) | FIX-** | — | CLEAN |
| (none) | EXEC_** | — | CLEAN |
| ~967 | exec21-optfix | `\textbf{Fix (exec21-optfix):}` in B2 frame, describing historical GEN-1 simulation bug fix | Not changed — lowercase `exec`, in-slide code name for a historical analysis campaign. Not an EXEC_ label per constraint. Noted for future reference. |

---

## 6. FIX-14 R=0.95 Audit

Other occurrences of R=0.95 in non-comment content that describe what "the code" contains (non-exhaustive):

| Line | Text | Frame | Status |
|------|------|-------|--------|
| ~371 | `$R=0.95$ constant` | Optical Surfaces (code listing description) | Describes data-generation code. Accurate for the dataset shown. Candidate for future correction if R=0.98 re-simulation is completed. |
| ~858 | `$R = 0.95$ constant (wavelength-independent)` | Simulation Model Details | Same — accurate for data. Future correction pending re-simulation. |

Lines 207, 803, 897, 924, 949 (other R=0.95 occurrences) describe data values, not code state; no correction needed.

---

## 7. MSI Worktree Status (exec24-deckfix branch)

exec24-deckfix was created on MSI (/mnt/d/SHiP/ej200) at pre-exec24 HEAD = 09dc9e2 (detached). During EXEC_24 compilation, talk_v6.tex was updated in-place via SSH pipe (not via git checkout); MSI working tree holds the compiled EXEC_24 PDF but the branch exec24-deckfix itself is at 09dc9e2. The canonical changes live in HOST main branch.

ej200_edge_scan worktree: holds main at 5e90025 (pre-exec23 state). Not touched during EXEC_24.

---

## 8. Synchronization Needed

- HOST (t0minidaq:/home/rrios/ej200) main: EXEC_24 edits made, not yet committed.
- MSI (/mnt/d/SHiP/ej200): talk_v6.tex updated in-place via SSH pipe; PDF compiled (38 pages, 472941 bytes); no git commit on MSI.
- GitHub: 12+ commits behind (pre-EXEC_23 state); push requires user approval.

Recommended next steps:
1. **Commit on HOST:** Five commits (or one combined) per FIX-11..15.
2. **Sync MSI:** Push HOST main to MSI via `git push ssh://reriosto@localhost:9022/mnt/d/SHiP/ej200 main:exec24-deckfix` (requires user approval).
3. **Push to GitHub:** `git push origin main` (requires user approval).

---

## 9. Open Items Carried Forward from EXEC_23

All OPEN-01..06 items from EXEC_23_REPORT.md remain open. No new open items introduced by EXEC_24.

---

No push. pre-exec24-260902 tag preserved.
