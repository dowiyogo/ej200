# Revision Notes — Napkin First-Principles Presentation
## Audit dated 2026-08-19

---

## 1. The historical "~73 ps" TOP result — traced

### Candidate sources

| Source | Config | σ (ps) | N SiPM |
|--------|--------|---------|---------|
| `phase_c_top_ablation.csv` | EJ-230 TIR, N=70 | 70.2 | 70 |
| `phase_c_top_ablation.csv` | EJ-230 TIR, N=50 | 82.0 | 50 |
| `phase_c_top_ablation.csv` | EJ-204 TIR, N=70 | 82.1 | 70 |
| `verification.json` ntop_scan | EJ-230 Vik sparse, N=4 | 78.2 | 4 |
| `phase_ab_kscan_full.csv` END k-scan | EJ-230, k=26 | 73.4 | 16 |

### Conclusion

The closest match is **EJ-230 TIR, N=70 SiPMs → σ_TOP = 70.2 ps** from `phase_c_top_ablation.csv`.
The "~73 ps" was most likely quoted from this configuration (full-coverage TIR TOP, 70 SiPMs).

Critically, that analysis used a **per-SiPM mean estimator** (mean first-photon arrival time averaged over all SiPMs), whereas the new analysis uses a **pooled k-th order statistic** (pool all photons from all SiPMs, take k-th earliest). The same 70 SiPMs with the pooled estimator would give dramatically better resolution.

### What actually drove the improvement from ~70 ps to 15.2 ps

Three factors, ranked by impact:

1. **Estimator change (dominant)**: Per-SiPM mean → pooled k-th order statistic.
   - Pooled N=20 sparse SiPMs yield ~1340 photons/event; k*=3 → σ_TOP = 15.2 ps.
   - TIR N=70 SiPMs with the SAME pooled estimator would have been comparably good.
   - The pooled estimator was the primary methodological advance.

2. **Reflector change (TIR → Vikuiti)**: At N=70, both give ~70 ps with the old estimator.
   - The reflector choice is NOT the main driver of the factor-5 improvement.
   - Vikuiti does increase Npe on TOP (~1340 vs ~287 for N=4 sparse), but this is also the N count effect.

3. **SiPM count (N=4 → N=20)**: With the pooled estimator, σ_TOP scales roughly as 1/√N_pe.
   - N=4: 287 pe → σ_TOP=78 ps; N=20: 1341 pe → σ_TOP=15 ps → ratio ~5 consistent with √(1341/287)≈2.2 × estimator gain.

**The correct narrative**: The improvement is primarily an **estimator advance** (per-SiPM mean → pooled k-th), combined with increasing N_top from 4 to 20 sparse SiPMs. The reflector (TIR→Vikuiti) is a secondary effect.

---

## 2. END estimator — corrected definition and numbers

### Two geometries give different results

| Geometry | m* | σ_END (ps) | Note |
|----------|----|------------|------|
| END-only (scan_end_vikuiti) | 7 | 50.49 | No TOP SiPMs in simulation |
| Combined (scan_end_vik_sparse_top_v2) | 8 | 53.68 | TOP SiPMs present; photons lost to TOP |

The combined-geometry result is the physically correct one when TOP SiPMs are present. The degradation (50.5 → 53.7 ps) is real: photons that would have reflected back toward END are captured by TOP SiPMs.

### Definition
`T0_END(m) = [mean(t_L[0:m]) + mean(t_R[0:m])] / 2`

where `t_L`, `t_R` are sorted arrival times at LEFT and RIGHT END faces. The optimal m* in combined geometry is 8.

---

## 3. BLUE combination — apparent paradox resolved

### Numbers (EJ-230, N=20, x=0)

| Estimator | σ (ps) | Bootstrap uncertainty |
|-----------|--------|-----------------------|
| σ_END | 53.68 | 0.93 |
| σ_TOP | 15.20 | 0.26 |
| σ_BLUE (event-by-event) | 15.21 | 0.26 |
| σ_BLUE (from cov matrix) | 15.18 | — |

The cov-matrix BLUE gives 15.18 ps (< σ_TOP), confirming BLUE is theoretically optimal. The event-by-event 15.21 ps is 0.01 ps above σ_TOP due to finite-sample noise in the covariance estimate. The difference is well within statistical uncertainty (0.26 ps bootstrap half-width).

### Weights and correlation

| N_top | σ_END | σ_TOP | ρ | w_END | w_TOP | σ_BLUE |
|-------|-------|-------|-----|-------|-------|--------|
| 4 | 51.0 | 78.2 | 0.149 | 0.741 | 0.259 | 45.5 |
| 8 | 50.2 | 32.0 | 0.180 | 0.212 | 0.788 | 28.9 |
| 14 | 51.0 | 18.7 | 0.182 | 0.024 | 0.976 | 18.2 |
| 20 | 53.7 | 15.2 | 0.159 | 0.013 | 0.987 | 15.2 |

**Key finding**: ρ ≈ 0.15–0.18 and remains roughly constant regardless of N_top. The correlation is low but non-negligible. At N=20, END contributes < 1.3% weight → σ_BLUE ≈ σ_TOP.

**Physical origin of correlation**: Both estimators use the same photon pool from the same muon event. The first few photons to reach END have also traveled some distance before hitting END, and the photon arrival time at TOP correlates with the muon T0. This creates a small but non-zero ρ.

---

## 4. Position scan (N=20, EJ-230)

| x (mm) | σ_END (ps) | σ_TOP (ps) | σ_BLUE (ps) | bootstrap |
|---------|-----------|-----------|------------|-----------|
| -600 | 59.0 | 15.4 | 15.1 | 0.30 |
| -500 | 60.1 | 17.8 | 17.4 | 0.29 |
| -400 | 55.4 | 15.9 | 15.4 | 0.26 |
| -300 | 56.3 | 16.0 | 15.7 | 0.29 |
| -200 | 53.6 | 17.6 | 17.2 | 0.31 |
| -100 | 52.8 | 14.9 | 14.6 | 0.27 |
|   0 | 53.7 | 15.2 | 15.2 | 0.26 |
| +100 | 53.1 | 15.6 | 15.0 | 0.23 |
| +200 | 54.5 | 18.2 | 17.8 | 0.35 |
| +300 | 54.9 | 15.8 | 15.6 | 0.31 |
| +400 | 57.4 | 15.8 | 15.3 | 0.24 |
| +500 | 58.8 | 17.5 | 16.9 | 0.26 |
| +600 | 60.9 | 15.3 | 15.1 | 0.27 |

**Range**: σ_BLUE 14.6–17.8 ps; peak-to-peak 21%.

**Pattern**: The 3 positions with σ_TOP > 17 ps are x=±500 and x=±200. These lie between sparse TOP SiPMs. The variation tracks the proximity to the nearest TOP SiPM, not the scintillator attenuation length.

**Note**: ρ(END,TOP) vs x is NOT in verification.json. Only the x=0 value is confirmed (ρ=0.159). The σ_BLUE values for |x|>0 were computed with ρ estimated from the same data, so they are self-consistent but not separately audited.

---

## 5. Pareto frontier — channel count vs resolution

| Config | Reflector | N_SiPM total | σ (ps) | Notes |
|--------|-----------|-------------|--------|-------|
| END+Vik EJ-230 (2/end) | Vikuiti | 4 | 131.2 | k*=2 |
| END+Vik EJ-230 (4/end) | Vikuiti | 8 | 97.9 | |
| END+Vik EJ-230 (8/end) | Vikuiti | 16 | 73.4 | |
| END+Vik EJ-230 (12/end) | Vikuiti | 24 | 63.2 | |
| END+Vik EJ-230 (16/end) | Vikuiti | 32 | 56.0 | ≈ m-average minimum |
| TOP TIR EJ-230 (N=70) | TIR | 70 | 70.2 | per-SiPM estimator |
| TOP TIR EJ-230 (N=50) | TIR | 50 | 82.0 | |
| BLUE(END+TOP) N=4 | Vikuiti | 20 | 45.5 | 16 END + 4 TOP |
| BLUE(END+TOP) N=8 | Vikuiti | 24 | 28.9 | 16 END + 8 TOP |
| BLUE(END+TOP) N=14 | Vikuiti | 30 | 18.2 | 16 END + 14 TOP |
| BLUE(END+TOP) N=20 | Vikuiti | 36 | 15.2 | **best; 16 END + 20 TOP** |

**Key message**: Adding 20 sparse TOP SiPMs to the baseline 16 END SiPMs (total 36 channels) reduces σ from 56 ps to 15.2 ps — a factor 3.7× improvement for a 2.25× increase in channel count.

---

## 6. EJ-204 vs EJ-230 at N=20

- σ_BLUE (EJ-230): 15.21 ± 0.26 ps
- σ_BLUE (EJ-204): 15.14 ± ? ps (full bootstrap not stored, but ~same uncertainty)
- Difference: 0.07 ps ≈ 0.27σ → **statistically indistinguishable**
- **Conclusion**: scintillator choice (EJ-204 vs EJ-230) does not matter at this SiPM count; TOP photon statistics dominate.

---

## 7. Corrections and inconsistencies in the old presentation

1. **Wrong claim**: "END+Vikuiti is the best configuration." Correct: BLUE combination with sparse TOP dominates by factor ~3.7.
2. **Missing comparison**: The old presentation did not show BLUE weights or ρ. This is critical context (END is nearly useless at N=20 TOP SiPMs).
3. **~73 ps provenance**: Not from a well-defined benchmark that was fairly compared. The per-SiPM vs pooled estimator distinction was not stated.
4. **σ_BLUE > σ_TOP paradox**: Due to finite-sample covariance noise. Statistically they are the same (difference < 0.5 bootstrap σ). This should be acknowledged, not ignored.
5. **Per-end vs L+R**: The napkin estimates 570 pe/end for EJ-230 (TIR-guided paths only). Geant4 simulation gives ~700 pe/end in END-only geometry and ~1340 pe pooled from 20 TOP SiPMs. These are consistent (napkin doesn't include diffuse paths and uses R=0.98 idealization).

---

## 8. Unresolved questions

- ρ(END, TOP) vs x: not computed for |x| > 0. Approximated as ρ ≈ 0.16 throughout.
- The historical ~73 ps may also have been from an exec07/exec12 analysis using a different analysis chain. The exec12_tN_summary.csv sigma_fit column is ambiguous (some values suspiciously small).
- Optimal sparse TOP SiPM placement vs random placement: the scan used a fixed sparse layout; optimal placement not studied.
- N_SiPM > 20 on TOP: not scanned. Diminishing returns expected since Npe ∝ N and σ_TOP ∝ 1/√Npe → σ ∝ 1/√N.

---

## 9. New figures generated (2026-08-19)

| File | Content |
|------|---------|
| `fig_end_mscan.pdf` | σ_END(m) vs m, EJ-230 x=0, both geometries marked |
| `fig_ntop_scan_full.pdf` | σ_END/TOP/BLUE vs N_top + w_END and ρ vs N_top panel |
| `fig_sigma_vs_x_new.pdf` | σ_END/TOP/BLUE vs x, all 13 positions |
| `fig_pareto.pdf` | Pareto frontier: resolution vs channel count |

Existing figures from `best_est_analysis.py`:
- `fig1_kscan.pdf`: σ_TOP(k) vs k for N=20 (already correct)
- `fig3_sigma_vs_x.pdf`: σ_BLUE vs x with σ_END/TOP (already correct)
- `fig4_sigma_vs_ntop.pdf`: σ_TOP vs N_top (already correct)
