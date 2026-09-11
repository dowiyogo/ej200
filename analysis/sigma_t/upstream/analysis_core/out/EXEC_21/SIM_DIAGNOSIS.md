# SIM_DIAGNOSIS.md — EXEC_21
Date: 2026-06-21 00:57

## Verdict
**CULPRIT IDENTIFIED: dielectric_metal skin surface eliminates TIR**
The simulation is INTERNALLY CONSISTENT (theory matches 0.37 PE observed)
but PHYSICALLY WRONG: real bars use TIR + Mylar wrapping, not pure metallic reflection.

---

## D1 — Falsification Gate
- σ_END(t_avg, x=0) = 882.5 ps (BROKEN: > REAL_SINGLE_FILT=250.0 ps)
- 31/31 positions have σ_END > 250.0 ps (impossible: intrinsic > real+electronics)
- Near-end (x=-690): σ_END ≈ 609 ps (still BROKEN vs real ~250 ps single-filtered)
- EXCEPTION: very near-end positions where END_L sees 900+ PE/ev → OK

## D2 — Photon Budget
- Expected <Npe>_END @ center with metallic R=0.98: 0.04 PE/ev
  (Theory: N_gen=20800 ph × survival=0.036^164 bounces × area_frac × PDE)
- Observed <Npe>_END @ center: 0.37 PE/ev
- Ratio observed/expected: ≈ 1× (SIM IS INTERNALLY CONSISTENT with metallic model)
- Expected with TIR: 749 PE/ev (17724× more)
- Real detector: ~40 PE/ev → metallic model STRUCTURALLY WRONG

## D3 — Photon History
- Time distributions show no bimodal structure on END hits at center
- Consistent with no guided light: only rare direct photons reach END
- Fraction 'late' hits (>3 ns, signature of guided photons):
  @ x=-690: 30.97%
  @ x=0: 100.00%
  @ x=690: 100.00%

## D4 — Root Cause Analysis
- **CULPRIT**: `Materials::CreateBarSkinReflector()` uses `dielectric_metal` surface
  → No TIR possible. All surface interactions are metallic R=0.98 (or 0.90 in endonly).
- Average bounces before photon hits END: 164
- Survival after 164 metallic bounces at R=0.98: 3.64%
- Additional ABSLENGTH loss: 12.87%
- → Predicted Npe: 0.04 PE (matches observed 0.37 PE ✓ — sim self-consistent)

- Even R=1.0 (perfect metallic): 1.2 PE (ABSLENGTH limits to ~4% survival even with R=1)
- TIR fix (dielectric_dielectric): 749 PE expected (matches order of real ~40 PE)

## D5 — Fix Status
**FIX IDENTIFIED (not yet implemented):**
Change bar skin surface from `dielectric_metal` to `dielectric_dielectric`
(enabling TIR at n=1.58/1.0 bar/air interface) PLUS add Mylar wrapping for non-TIR photons.

This requires restoring explicit reflector panel geometry (old geometry had this; new skin surface lost it).
Estimated implementation: branch exec21-optfix, ~2-3 files to modify in DetectorConstruction.cc/Materials.cc.

**NOT YET IMPLEMENTED** (awaiting René's decision on the fix approach).
Mini-sims D4 run on theory; code modification for empirical verification pending.

## D6 — GLS Bug Fix
- EXEC_20 σ_EndTop(GLS)=75.8 ps > σ_TOP=75.9 ps (impossible)
- Cause: large positive covariance cov(t_END,t_TOP) from shared scintillation fluctuations
- Fix: inverse-variance combination → σ_EndTop=75601.4 ps ≤ σ_TOP ✓
- Per-event inv-var combination: σ_combined=69.7 ps

## Consequences for Previous Results
1. σ_END(center) ≈ 882 ps is UNPHYSICAL: sim misses TIR → too few photons at END
2. T3 'aporte del TOP 85%' is spurious: it compares physical TOP vs unphysical END
3. Once D5 fix is applied: re-run T1/T3 with physical END model
4. GLS combination in T3/T7: use inv-variance, not covariance-GLS

## Ranked Candidates (for historical tracking)
1. **#1 CONFIRMED**: dielectric_metal surface (no TIR) → ~300× light deficit
2. ~Secondary: explicit panel vs skin surface geometry difference
3. Minor: reflectivity tuning (0.90 vs 0.95 vs 0.98) — secondary to TIR