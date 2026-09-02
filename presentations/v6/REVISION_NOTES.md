# REVISION_NOTES.md — Corrections and improvements from v5 to v6

---

## 1. Comprehensive first-principles narrative (new)

v5 started with simulation results. v6 opens with "Napkin first, Geant4 second":
θ_c → f_TIR → bounces → reflector attenuation → photon budget → G4 confirms.
This makes the physical reasoning explicit and the Geant4 role clear (confirmation, not discovery).

---

## 2. Two-photon-population framework (new)

v5 did not explicitly name or explain Population I (TIR-guided) vs Population II (reflector-recovered).
v6 introduces this framework in a dedicated slide (S5) and uses it throughout to explain
why G4 gives +23% over the napkin and why the arrival time spectrum has two features.

---

## 3. TikZ diagrams (new)

v6 adds five original TikZ diagrams not present in v5:
- Bar geometry with SiPMs (S3)
- TIR diagram: TIR ray + non-TIR escape + Vikuiti recovery (S4)
- Two-population arrival spectrum schematic (S5)
- Zig-zag bounce counting with N_bounce label (S6)
- GEN-3 vs GEN-1 air-gap geometry comparison (S9)

---

## 4. BLUE weights: two conventions → documented ambiguity below statistical resolution

v5 BLUE: w_END = 0.086 (from np.var() empirical variance)
v6 BLUE: w_END = 0.013 (from IQR-based σ² — consistent with best_est_analysis.py)

Both are self-consistent within their own conventions. v6 originally presented 0.013 as canonical
(σ_IQR is the reported metric, so consistency requires IQR-based weights).

**EXEC_23 FIX-03 correction (2026-09-02):** w_END is NOT resolved at n=5000. A CP-B w-scan
(81 values of w, 1000-resample bootstrap) on the same 5000-event sample shows σ_IQR(w) is
flat to within one bootstrap SE over w ∈ [0, 0.10] (total variation 0.26 ps; bootstrap SE ≈ 0.17 ps).
Three conventions give w_END ∈ {0.013, 0.051, 0.086} with σ_IQR^event within 0.18 ps.
The "correction" from 0.086 to 0.013 is now described as a documented ambiguity, not a resolved
correction. Appendix G2 presents all three and closes with the not-resolved statement.
See FINAL_NUMBERS.md §"Point performance" for the full three-convention table.
Sidecar: analysis/optim/root_best_est/blue_wscan_x0.{csv,meta.json,root}

---

## 5. G4 optical model explained, not just cited

v5 stated the model configuration without explaining the physics implications.
v6 explains WHY each surface choice matters:
- CreateBarSurface(): natural Fresnel TIR (not artificial)
- CreateBarSkinReflector() as border: R=0.95 acts only on non-TIR photons
- GEN-1 skin surface on bar: eliminates TIR → 300× photon deficit (now explained)

---

## 6. Appendix: historical GEN-1 bug fully documented

v5 mentioned the fix briefly. v6 Appendix B2 gives the complete explanation:
- dielectric_metal skin on bar volume → no TIR possible
- 0.37 pe/event claimed — 1500× below napkin prediction → physically impossible
- Napkin falsifies the GEN-1 claim before G4 fix was applied

---

## 7. R=0.95 vs R=0.98 discrepancy explicitly flagged

v5 mentioned R=0.95 in the reflector model but did not contrast with R=0.98 napkin.
v6 Appendix A2 explicitly states:
- Code (Materials.cc): R=0.95 constant
- Napkin uses R=0.98 (manufacturer Vikuiti ESR spec)
- Code comment says "R=0.98" but implements 0.95 — discrepancy should be corrected

---

## 8. Content audit and configuration audit added

New documents CONTENT_AUDIT.md and CONFIGURATION_AUDIT.md preserve the full
provenance of every number in the presentation. No results from prior presentations
are silently dropped — all are either in the main deck or explicitly moved to appendix
with their status (valid/superseded/historical/incorrect-interpretation).

---

## 9. Global vs oracle distinction reinforced

v5 introduced this distinction correctly. v6 reinforces it by:
- Calling out the 3.4 ps mean gap explicitly
- Explaining WHY spikes occur at x=±200,±500 mm (k=3 is locally suboptimal)
- Clarifying that oracle is "diagnostic lower bound, not deployable"

---

## 10. Appendix structure: 13 backup slides in blocks A–J

v5 had 5 backup slides. v6 has 13 appendix slides organized by topic:
A: Optical model (surface definitions, R discrepancy, air gap, PDE)
B: Historical simulations (GEN generations, 1-PE anomaly)
C: Materials (full property table)
D: Photon budget (detailed chain)
E: END timing (full table, v_eff discrepancy)
F: TOP estimator (k-scan, edge behavior)
G: BLUE statistics (weight derivation, two variants)
H: Position (centroid saturation, LOO results)
I: Jitter (timing budget, SPTR caveat)
J: Provenance (data files, scripts, correction table, figures index)
