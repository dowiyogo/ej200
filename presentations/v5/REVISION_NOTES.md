# REVISION_NOTES.md — Corrections from v4 to v5

---

## 1. Author and institution (CRITICAL)

**v4**: `R. Ríos Pino / SHiP Collaboration / CINVESTAV`  
**v5**: `René Ríos / Universidad de La Serena (ULS)`

Source: user correction during session; confirmed by `presentations/sim_status_2026-06-08/talk.tex` line 20.

---

## 2. EJ-230 attenuation length (CRITICAL)

**v4**: stated 900 mm  
**v5**: **1200 mm**

Source: `src/Materials.cc` — `AddScintillatorProperties(mat, ..., 120.0 * cm, ...)` → 120 cm = 1200 mm.  
The value 900 mm was incorrect and not present anywhere in the source code.

---

## 3. Reflector model

**v4**: described as "Vikuiti ESR spectral reflector"  
**v5**: **dielectric_dielectric, R=0.95 constant (wavelength-independent)**

Source: `src/Materials.cc`, `CreateBarSkinReflector()`. The actual Vikuiti spectral table is NOT used. The reflectivity is a constant 0.95. The code comment says "R=0.98 Vikuiti ESR" but the actual value set is 0.95.

---

## 4. fig3_sigma_vs_x: global vs oracle (CRITICAL SCIENTIFIC ERROR)

**v4**: figure described as "global k=3, m=8" but was computed with oracle k*(x), m*(x) — locally optimal per position (same-sample).  
**v5**: two separate curves:
- **Global (k=3, m=8)**: deployable, fixed parameters — mean σ_BLUE=19.5 ps
- **Oracle k*(x), m*(x)**: diagnostic lower bound only — mean σ_BLUE=16.1 ps

The gap (3.4 ps mean) represents the deployment cost.  
Spikes in global curve at x=±200, ±500 mm where k=3 is not locally optimal.

---

## 5. BLUE definition clarification

**v4**: unclear whether covariance matrix or σ_IQR used for weights.  
**v5**: explicit distinction:
- BLUE **weights** from ordinary Var/Cov matrix: `C = np.cov(T0e, T0t)` → w_END, w_TOP
- **σ_IQR** = (Q75−Q25)/1.349 is an **independent** robust performance metric — NOT derived from C

At x=0: w_END=0.086 (TOP dominates); σ_BLUE = 15.0 ps.

---

## 6. v_eff nonlinearity description

**v4**: χ²/ndf described as showing "mild" nonlinearity  
**v5**: **χ²/ndf = 1250/11, p ≈ 0 — statistically highly significant nonlinearity**

Residual RMS = 14.5 ps. S-shape pattern driven by edge photon pool asymmetry.  
"Mild" removed; replaced with "statistically significant" and explicit residual plot.

---

## 7. PDE interpretation

**v4**: "mean PDE = 0.607" stated without clarification  
**v5**: explicit note that 0.607 is the **mean over detected photons** — selection-biased.  
This is NOT the overall detection fraction. Correct budget requires N_detected/N_emitted.

---

## 8. N_TOP scan uses physical geometry files

**v4**: N_TOP scan performed by software ablation (subsampling SiPMs from ntop20 data)  
**v5**: uses **actual physical simulation files**: `photon_hits_ej230_ntop{4,8,14,20}_x0.root`

This matters because different N_TOP geometries have physically different photon distributions (different solid angles, different inter-SiPM shadowing).

Results differ slightly: ntop14 oracle k*=2 (not k*=1 as software ablation would give).

---

## 9. TOP position leave-one-out (new analysis)

Not present in v4. New in v5:
- N_pe centroid → LOO calibration on 12/13 positions → predict held-out
- Result: σ_x=17.2 mm (linear), large bias at bar ends (±29 mm at ±600 mm)
- END Δt gives σ_x=7.9 mm — 2× better, no edge bias
- Centroid saturates: ±600 mm muon → centroid at ±313 mm

---

## 10. Slide structure improvements

- One message per slide (assertion-style titles)
- Separate slides for: global results table, oracle vs global comparison, N_TOP scan, LOO position
- Consistent colour scheme: END=kAzure+2, TOP=kOrange+1, BLUE=kGreen+3, Oracle=kGray+2
- Backup slides: PDE clarification, reflector model, oracle table, train/test validation, centroid saturation
- Navigation bars suppressed (`\setbeamertemplate{headline}{}`) to prevent Warsaw overflow
