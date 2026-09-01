# CONTENT_AUDIT.md — Content Preservation Matrix for talk_v6

Every scientifically meaningful result/slide from previous presentations.
Statuses: valid | superseded | historical | config-dependent | requires-reanalysis | incorrect-interpretation | unresolved

---

## FROM: talks/napkin_first_principles/tex/talk_v3.tex (René Ríos / ULS, Aug 2026)

| Content | Configuration | Status | v6 location |
|---------|---------------|--------|-------------|
| θ_c = arcsin(1/1.58) = 39.27° | analytic | valid | Main S6 |
| TIR slab fraction f_TIR = 2cos(θ_c)−1 = 0.548 | Knoll 8.15 | valid | Main S6 |
| f_TIR per end ≈ 0.274 | analytic | valid | Main S6 |
| Two-wall constraint note (rectangle ≠ slab) | analytic | valid | Main S6 |
| v_g = c/n = 189.7 mm/ns | n=1.58 | valid | Main S8 |
| v_eff = 184.6 mm/ns (from Δt, first-photon) | GEN-3 END | config-dependent | Appendix E7 |
| σ_x ≈ 10.3 mm (from 184.6 mm/ns) | GEN-3 END | config-dependent | Appendix E7 |
| N_bounce ≈ 86 (thin face H=10 mm) | analytic | valid | Main S9 |
| N_bounce ≈ 14 (wide face W=60 mm) | analytic | valid | Main S9 |
| Λ_refl ≈ 404 mm (H=10 mm, R=0.98) | analytic, R=0.98 | config-dependent (R=0.98) | Main S9 (note R=0.95) |
| EJ-200 > EJ-204 > EJ-230 in Npe | GEN-3 | valid | Main S10 |
| EJ-230 < EJ-204 < EJ-200 in σ_t | GEN-3 | valid | Main S10 |
| Photon budget: 570 PE/end (napkin) vs 701 PE/end (G4) | GEN-3, EJ-230, x=0 | valid (20% surplus from reflector-recovered) | Main S11 |
| Per-end vs L+R distinction | analytic | valid | Main S11 |
| Two photon populations (TIR-guided, reflector-recovered) | conceptual | valid | Main S7 |
| Napkin vs Geant4 scoreboard | GEN-3 | valid | Main S10/11 |
| m-average END estimator: σ_END=53.7 ps (m*=8, x=0) | GEN-3 hybrid | valid | Main S14 |
| Top k-th order statistic: σ_TOP=15.2 ps (k*=3, N=20, x=0) | GEN-3 | valid | Main S15 |
| BLUE: w_END=0.013, w_TOP=0.987, σ_BLUE=15.21 ps | GEN-3, x=0, N=20 | valid | Main S18 |
| BLUE uncorrelated benchmark: 14.6 ps | GEN-3 | valid | Appendix G1 |
| BLUE vs N_TOP table | GEN-3 | valid | Main S19 |
| σ(x) variation 14.6-17.8 ps, peaks at x=±200,±500 mm | GEN-3, oracle k*(x) | superseded by global/oracle distinction | Main S21 + note |
| Pareto frontier: 15-channel to 36-channel | GEN-3 | valid | Main S24 |
| END provides position; TOP provides time | GEN-3 | valid | Main S26 |
| GROUPVEL audit (verified zero superluminal photons) | GEN-3 | valid | Appendix A7 |
| 1-PE anomaly (historical, before exec21-optfix) | GEN-1 (INVALID) | historical/explained | Appendix B6 |
| EJ-204 vs EJ-230: diff 0.07±0.38 ps | GEN-3, N=20, x=0 | valid | Appendix C5 |
| Specular vs diffuse: campaign not yet run | prospective | unresolved | Appendix A? |

---

## FROM: presentations/best_est_2026-08-17/talk.tex

| Content | Configuration | Status | v6 location |
|---------|---------------|--------|-------------|
| Near-field (~35 mm) vs far-field (~700 mm) path comparison | GEN-3 | valid | Main S15 |
| τ_prop END ≈ 3.7 ns vs TOP ≈ 185 ps | GEN-3 | valid | Main S15 |
| σ_END m*=8 = 53.68 ps | GEN-3, x=0 | valid | Main S14 |
| σ_TOP k*=3 = 15.20 ps | GEN-3, x=0, N=20 | valid | Main S15 |
| σ_BLUE = 15.21 ps (w_END=0.013) | GEN-3, x=0, N=20 | valid | Main S18 |
| N_TOP table: {4,8,14,20} → {45.5,28.9,18.2,15.2} ps | GEN-3 | valid | Main S19 |
| σ_BLUE(x) range 14.6–17.8 ps | GEN-3, oracle | superseded (now use global/oracle) | Appendix |
| Historical ~70 ps per-SiPM | older config | historical | Appendix B2 |
| BLUE cross-check: σ_cov=15.18 vs σ_ev=15.21 ps | GEN-3 | valid | Appendix G1 |

---

## FROM: analysis/presentation_v4/talk_v4.tex

| Content | Configuration | Status | v6 location |
|---------|---------------|--------|-------------|
| Author "R. Ríos Pino / CINVESTAV" | — | WRONG — corrected to René Ríos (ULS) | — |
| EJ-230 att = 900 mm stated | — | WRONG — actual value is 1200 mm | — |
| fig3_sigma_vs_x labeled "global k=3" but used oracle k*(x) | GEN-3 | incorrect interpretation | Appendix with note |
| Photon budget chain | GEN-3 | valid | Main S11 |
| Jitter slide (figG_jitter.pdf) | GEN-3 | valid | Main S25 |
| Estimator comparison figH | GEN-3 | valid | Appendix F2 |

---

## FROM: analysis/presentation_v5/talk_v5.tex

| Content | Configuration | Status | v6 location |
|---------|---------------|--------|-------------|
| Author "René Ríos (ULS)" | — | correct | Main S1 |
| EJ-230 att = 1200 mm | GEN-3 | correct | Main S5 |
| Global (k=3,m=8) vs oracle (k*(x),m*(x)) correctly separated | GEN-3 | valid | Main S20-21 |
| BLUE weights: w_END=0.086 (empirical variance) | GEN-3 | superseded — use best_est (0.013) | Appendix G (as comparison) |
| v_eff = 177.6 mm/ns, χ²/ndf=1250/11 | GEN-3, m=8 avg | valid | Main S16 |
| Nonlinearity residual RMS = 14.5 ps | GEN-3 | valid | Main S16 |
| Global mean σ_BLUE = 19.5 ps across x | GEN-3, k=3,m=8 | valid | Main S21 |
| Oracle mean σ_BLUE = 16.1 ps | GEN-3, oracle | valid | Main S21 |
| TOP LOO: σ_x=17.2 mm linear, END σ_x=7.9 mm | GEN-3 | valid | Main S23 |
| N_TOP physical scan: 4→45, 8→32, 14→18, 20→15 ps | GEN-3 | valid | Main S19 |
| PDE=0.607: selection-biased (mean over detected photons) | GEN-3 | valid | Appendix A8 |
| Reflector R=0.95 constant (dielectric_dielectric) | GEN-3 | valid | Main S12 + Appendix A6 |

---

## FROM: presentations/optim_2026-08-16/ and optim_2026-08-17/

| Content | Configuration | Status | v6 location |
|---------|---------------|--------|-------------|
| N_TOP ablation studies | GEN-3 with software masking | config-dependent (software ablation ≠ physical removal) | Appendix B3 |
| Pareto with ablation results | mixed | requires-reanalysis (use physical files only) | Main S24 (physical only) |

---

## ITEMS REQUIRING SPECIAL TREATMENT

### B6. Historical 1-PE anomaly
- Source: GEN-1 (dielectric_metal surface)
- The ~0.37 PE/event result is a Geant4 surface bug artifact, not physics
- Ratio Vikuiti/TIR ≈ 10³ is falsifiable and was shown to be impossible by napkin
- v6 location: Appendix B10 — labeled prominently as surface-bug artifact

### v_eff discrepancy (184.6 vs 177.6 mm/ns)
- Both are valid for their respective estimators
- 184.6 mm/ns: best-fit to first-photon Δt (napkin analysis)
- 177.6 mm/ns: best-fit to m=8 average Δt (v5 global analysis)
- v6: explain the distinction; use whichever is relevant to the estimator being discussed

### BLUE weights discrepancy (0.013 vs 0.086)
- Both are technically correct but use different variance definitions
- 0.013 (best_est): IQR-based σ² in BLUE formula → physically clean
- 0.086 (v5): empirical np.var() → includes heavy tails
- v6: use 0.013 as canonical; explain in appendix G2 that both exist

### Figure fig3_sigma_vs_x (oracle mislabeled as global)
- This was from presentation_v4 and was the key error identified in v5
- v6 uses v5_sigma_vs_x.pdf which correctly separates global and oracle
- Old figure kept in appendix with label: "Historical — oracle k*(x) mislabeled as global"
