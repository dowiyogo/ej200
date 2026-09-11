# EXEC_23_REPORT.md — Internal Consistency Campaign for talk_v6.tex

Generated: 2026-09-02  
Authoritative clone: t0minidaq:/home/rrios/ej200, branch main  
Rollback tag: pre-exec23-260902 = 5e90025  
HEAD after campaign: cb0fd9a

---

## 1. FIX-01..10 Status

| FIX | Status | Commit | Files | Before → After |
|-----|--------|--------|-------|----------------|
| FIX-01 | DONE | fde596b | talk_v6.tex, results_macros.py, results_macros.tex, napkin.py | S8 material table: hardcoded 47.6/52.1/49.6 ps → `\RessigmabxzERO*` / `\RessigmabmEAN*` macros (sourced from phase_ab_optimal.csv). Two σ rows: x=0 and bar mean. Footnote rewritten. C1/figM1 caption updated. EJ-230 wins at x=0 (50.49 ps) AND bar mean (51.95 ps) — old claim "EJ-200 wins at x=0" was wrong. |
| FIX-01b | OPEN | 7e7180f | presentations/v6/README.md | figM1/figM2 figure sidecars missing from phase_ab.py — logged as [OPEN] in README; deferred pending R=0.98 re-simulation. |
| FIX-02 | DONE | 80feecc | talk_v6.tex | S-G4-validates frame: σ_END=53.68 ps relabeled as hybrid geometry (N_TOP=20, m=8), not END-only. σ_END=50.49 ps (END-only, m*=7) sourced from scan_end_vikuiti. |
| FIX-03 | DONE | f1591ec | talk_v6.tex, FINAL_NUMBERS.md, CONTENT_AUDIT.md, REVISION_NOTES.md, README.md, analysis/optim/root_best_est/blue_wscan_x0.{csv,meta.json,root} | w_END not resolved at n=5000. S15: ρ=0.16 removed from formula block; C_ET=196.3 ps² added; w_END → ≲0.1; physics note on TOP dominance. G1/G2 combined frame: three covariance conventions, flatness note, not-resolved closing (frame split into G1 + G2 by EXEC_24 FIX-11). S16/S21 updated. E1 footnote: ρ=0.16 and w_E=0.013 labeled v5/empirical, BLUE-consistent ρ=0.241. |
| FIX-04 | DONE | 1fdefe2 | talk_v6.tex | σ_B reconciliation footnote added to S16 and S17. |
| FIX-05 | DONE | 139a16e | talk_v6.tex, napkin_macros.tex (via napkin.py) | Λ_refl(R=0.95): hardcoded 240 mm → 159.4 mm via `\NapLambdaReflRlow` macro. Computed: −L_bar/(2·ln R) = 700/(2·ln(1/0.95)) ≈ 159.4 mm (was 240 mm which used R=0.98 value). |
| FIX-06 | DONE | 4bf53a8 | talk_v6.tex, napkin_macros.tex (via napkin.py) | Path-to-END: hardcoded "700 mm" → range "700–1106 mm" via macros `\NapLhalfmm` and `\NappathTIRlimmm`. Upper bound = L/(2 sin θ_c) = 1106 mm. |
| FIX-07 | DONE | 351f510 | talk_v6.tex | σ_x formula: spurious /√2 removed. Correct: σ_x = σ_END · v_eff / (L/2). Was: σ_x = σ_END · v_eff / (L/2 · √2). |
| FIX-08 | DONE | 540f591 | talk_v6.tex, CONFIGURATION_AUDIT.md | B1 table: "GEN-2" column relabeled "GEN-3 END-only (scan_end_vikuiti)". Air gap: corrected to Yes (0.10 mm). Reflector: corrected to border surface dielectric_dielectric R=0.95. σ_END footnote uses `\sigmaENDonlyval`. TOP-interception claim demoted to hypothesis (no controlled experiment isolating the variable). |
| FIX-09 | DONE | 76c2520 | talk_v6.tex | A2 block: reflector R discrepancy corrected. Code has R=0.98 (commit c7acb7a). Dataset scan_end_vikuiti used R=0.95. Slide now states the discrepancy explicitly. |
| FIX-10 | DONE | cb0fd9a | talk_v6.tex | Three overflow fixes (see §4). I1 split into I1/I2. Build environment: MSI WSL2 (lualatex). |

---

## 2. Full Commit Log (pre-exec23-260902..HEAD)

```
cb0fd9a EXEC_23 FIX-10: split overflowing frames (I1/I2); fix S8 and photon-budget overflows
f1591ec EXEC_23 FIX-03: w_END not resolved (CP-B w-scan); three covariance conventions in G2
540f591 EXEC_23 FIX-08: relabel B1 column GEN-2 → GEN-3 END-only (scan_end_vikuiti)
1fdefe2 EXEC_23 FIX-04: add σ_B reconciliation footnote to S16 and S17
7e7180f EXEC_23 FIX-01b: log figure sidecar defect as open item in README
fde596b EXEC_23 FIX-01: material σ table sourced from phase_ab_optimal.csv via macros
80feecc EXEC_23 FIX-02: S-G4-validates frame relabeled — 53.68 ps is hybrid geometry, not END-only
76c2520 EXEC_23 FIX-09: A2 block corrected — code now has R=0.98 (c7acb7a), data used R=0.95
351f510 EXEC_23 FIX-07: σ_x formula corrected, spurious /√2 removed
4bf53a8 EXEC_23 FIX-06: path-to-END corrected (700–1106 mm range via macros)
139a16e EXEC_23 FIX-05: Λ_refl(R=0.95) corrected 240→159.4 mm via macro
```

---

## 3. verify_deck.py

**Absent.** No verify_deck.py script exists in presentations/v6/scripts/ or anywhere in the
repository. Internal consistency was verified manually via EXEC_23 FIX-01..10 and by
inspecting each changed slide against the authoritative source in FINAL_NUMBERS.md and
CONFIGURATION_AUDIT.md.

---

## 4. Overflow Inventory

### Baseline (pre-exec23-260902 = 5e90025)

| Log line | Overflow | Size | Page | Frame |
|----------|----------|------|------|-------|
| 1362 | `\vbox` | 129.63 pt | 34 | I1 — Jitter and Timing Budget |

### HEAD before FIX-10 (f1591ec)

| Log line | Overflow | Size | Page | Frame | Cause |
|----------|----------|------|------|-------|-------|
| 1194 | `\hbox` | 10.65 pt | 8 | Material Choice: EJ-230 is the Fastest Scintillator Here | FIX-01: two σ rows with wider `($x{=}0$)` label vs old `(G4)` |
| 1219 | `\vbox` | 5.02 pt | 10 | Geant4 Validates the Napkin | FIX-08: added `\footnotesize END-only` footnote line |
| 1385 | `\vbox` | 129.63 pt | 34 | I1 — Jitter and Timing Budget | Pre-existing |

### HEAD after FIX-10 (cb0fd9a) — CLEAN

No Overfull \vbox or \hbox warnings. 37 pages.

### Per-frame inspection (overflow-affected frames)

- **page 8 inspected: text complete.** S8 material table column rebalanced (0.48→0.52 table, 0.48→0.44 figures). Table fits within column.
- **page 10 inspected: text complete.** Post-table spacer reduced 6pt→1pt. All content visible.
- **page 34 inspected (I1): text complete.** Frame retitled "I1 — Jitter: Sources and SPTR Provenance". Contains: sources itemize list (emission, propagation, SPTR with Lee 2025 measurement, electronic not-included), fig_sigma_t_x caption.
- **page 35 inspected (I2, new): text complete.** "I2 — Timing Simulation Coverage". Contains: simulation includes/does not include, alert σ_TOP=15.2 ps lower bound, order-statistic caveat.

---

## 5. Places Where Data Disagreed with a Stated Hypothesis

These are results, not corrections, and are the most important findings of EXEC_23.

**D1 — EJ-230 wins at BOTH x=0 AND bar mean (FIX-01)**  
Old slide claimed "EJ-200 wins at x=0 (47.6 vs 49.6 ps)." Source of truth (phase_ab_optimal.csv, GEN-3 END-only, m*-optimal) gives EJ-230 at x=0: 50.49 ps < EJ-204: 53.36 ps < EJ-200: 54.86 ps. The reversal is because the orphaned literal 47.6/52.1/49.6 ps came from an earlier simulation run with a different geometry or parameters than the current GEN-3 END-only dataset. EJ-230 is now the winner at both positions.

**D2 — TOP-intercept claim not isolatable (FIX-08)**  
B1 table implied that σ_END grows from 50.5 ps (END-only, N_TOP=0) to 53.68 ps (hybrid, N_TOP=20) because TOP SiPMs intercept ~3 ps worth of early photons. This is a plausible hypothesis but the comparison confounds geometry (END-only vs hybrid) and estimator configuration simultaneously. No controlled experiment holds one variable fixed. Demoted from claim to hypothesis in FIX-08.

**D3 — w_END not resolved at n=5000 (FIX-03)**  
CP-B w-scan (81 points, 1000-resample bootstrap) on the same 5000-event sample shows σ_IQR(w) flat to within one bootstrap SE over w ∈ [0, 0.10] (total variation 0.26 ps; SE ≈ 0.17 ps). Three covariance conventions give w_END ∈ {0.013, 0.051, 0.086} with σ_IQR^event within 0.18 ps. The "correction" from v5 w=0.086 to v6 w=0.013 is not resolved by the available dataset. Both are consistent with the data.

**D4 — ρ(x) does not explain empirical w* non-monotonicity (FIX-03)**  
Analytic BLUE gives monotonic w_END ordering with σ_T: x=0 (w=0.013) < x=+600 (w=0.041) < x=+200 (w=0.106), consistent with σ_T ordering. Empirical w* does not follow this ordering. This non-monotonicity is statistical noise in an unresolved argmin at n=5000, not a physics effect driven by ρ(x).

---

## 6. Open Items

**[OPEN-01] FIX-01b — figM1/figM2 figure sidecars**  
analysis/optim/phase_ab.py saves PNG figures but emits no CSV or meta.json sidecars.
figM1/figM2/fig_sigma_t_x/fig_npe_x in presentations/v6/figs/ lack provenance metadata.
Deferred until R=0.98 re-simulation campaign.

**[OPEN-02] R=0.98 re-simulation**  
All reported results (σ_END, σ_TOP, σ_BLUE, N_pe) use R=0.95 data (scan_end_vikuiti, scan_end_vik_sparse_top_v2). Code was corrected to R=0.98 in commit c7acb7a. Re-simulation required to validate whether design conclusions hold at R=0.98. Expected effect: Λ_refl grows from 159.4 mm to 404.6 mm; more reflector-recovered photons; σ_END and σ_TOP both expected to improve.

**[OPEN-03] SPTR bench measurement pending EOS data**  
README states: SPTR ≈ 85–92 ps RMS at 10 V OV (own bench, laser intensity scan). Do NOT add to I1/I2 until the formal EOS measurement is available. See analysis/timing/TODO.md.

**[OPEN-04] w_END not resolved at n=5000**  
Three covariance conventions give w_END ∈ {0.013, 0.051, 0.086}; all within 0.18 ps ≈ bootstrap SE. Resolution requires ≥50k events at x=0 under hybrid geometry. Documented in G2 appendix, FINAL_NUMBERS.md, REVISION_NOTES.md §4.
Sidecar: analysis/optim/root_best_est/blue_wscan_x0.{csv,meta.json,root}

**[OPEN-05] k-adaptive estimator**  
Global k=3 is suboptimal at x=±200, ±500 mm. Oracle gap = 3.4 ps (mean). Adaptive k*(x) would require a real-time position estimate. Not addressed in this campaign.

**[OPEN-06] TOP-intercept effect isolation**  
The hypothesis that TOP SiPMs intercept early photons and increase σ_END from 50.5 to 53.68 ps has not been tested in a controlled experiment (one variable at a time). Requires a hybrid geometry simulation with N_TOP=0 to test the hypothesis directly.

---

## 7. sim_sigma Dict and FIX-01b Status

**sim_sigma (napkin.py):** The dict is marked `# DEPRECATED (EXEC_23 FIX-01): orphaned values; source of truth is phase_ab_optimal.csv via results_macros.py` with `# ps — DO NOT USE` on the value line. It is NOT removed — it still runs and emits `\NapsigmasimejA{47.60\,\text{ps}}`, `\NapsigmasimejB`, `\NapsigmasimejC` to napkin_macros.tex. These three macros are present in napkin_macros.tex but are NOT referenced anywhere in talk_v6.tex. They cannot be silently reused in the current deck without an explicit `\NapsigmasimejX` call in the .tex source. The deprecation warning is visible to any reader of napkin.py but does not prevent the values from being generated.

**FIX-01b (figM1/figM2 sidecars):** Logged as [OPEN] in presentations/v6/README.md since commit 7e7180f. No action taken; see OPEN-01 above.

---

## 8. Build Environment

- **Compilation host:** MSI workstation, WSL2 (Linux 6.18.33.2-microsoft-standard-WSL2)
- **lualatex:** /usr/bin/lualatex (texlive)
- **latexmk:** /home/reriosto/.local/bin/latexmk
- **beamer.cls, tikz.sty:** /usr/share/texlive/texmf-dist/
- **pgfplots.sty:** Present via MiKTeX Windows installation (not used in talk_v6.tex)
- **Access:** reverse SSH tunnel; from t0minidaq: `ssh -p 9022 reriosto@localhost`
- **t0minidaq:** no TeX installation; all compilation on MSI
- **Baseline PDF:** /mnt/d/SHiP/ej200_edge_scan/presentations/v6/talk_v6.pdf (446 KB, 36 pages)
- **Final PDF:** /mnt/d/SHiP/ej200/presentations/v6/talk_v6.pdf (462 KB, 37 pages)
- **Rasters:** /tmp/rasters/baseline/ (36 PNGs), /tmp/rasters/head2/ (37 PNGs), 80 dpi

---

No push. pre-exec23-260902 tag preserved.
