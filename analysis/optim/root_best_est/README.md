# root_best_est — EndSparseTop reanalysis (2026-08-17)

Independent reanalysis of the EndSparseTop N=20 optimal configuration.
All numbers verified directly from ROOT files in `scan_end_vik_sparse_top_v2/`.

## Commands to reproduce

```bash
cd /home/rrios/ej200/analysis/optim/root_best_est
python3 best_est_analysis.py          # ~8 min (reads 17 ROOT files + bootstrap)
python3 -c "exec(open('fig4_twopanel.py').read())"  # update fig4 (2-panel)
```

## Input files used

- `results/scan_end_vik_sparse_top_v2/photon_hits_ej230_ntop{4,8,14,20}_x{pos}.root`
  (104 files total; analysis uses 17 key files)
- `results/scan_end_vikuiti/photon_hits_ej230_x0.root`  (END+Vik reference)

## Output files

| File | Description |
|------|-------------|
| `fig1_kscan.{pdf,png}` | σ_TOP(k) for N=20 EJ-230 x=0 |
| `fig2_residuals.{pdf,png}` | T residual histograms (END/TOP/BLUE) |
| `fig3_sigma_vs_x.{pdf,png}` | σ vs position for N=20, all 13 positions |
| `fig4_sigma_vs_ntop.{pdf,png}` | σ_TOP vs N_top (two-panel) |
| `figures.root` | All ROOT objects (TGraphErrors, TH1D, TF1) |
| `verification.json` | All numbers used in presentation |

## Key verified numbers (EJ-230, N=20, x=0)

| Quantity | Value | Source |
|----------|-------|--------|
| σ_END (m*=8) | 53.68 ± 0.93 ps | ROOT files, IQR |
| σ_TOP (k*=3) | 15.20 ± 0.27 ps | ROOT files, IQR |
| σ_BLUE (event) | 15.21 ± 0.26 ps | ROOT files, IQR |
| σ_BLUE (from cov) | 15.18 ps | covariance cross-check ✓ |
| w_END | 0.013 | BLUE formula |
| w_TOP | 0.987 | BLUE formula |
| ρ(T_END, T_TOP) | 0.16 | empirical covariance |
| N_pe_top | 1341.1 | face_type=2 hits/event |

## EJ-204 vs EJ-230 at N=20

| | EJ-230 | EJ-204 |
|--|--------|--------|
| σ_BLUE | 15.21 ± 0.26 ps | 15.14 ± 0.28 ps |
| diff | -0.07 ± 0.38 ps (0.2σ) | — |
| **Conclusion** | **Statistically equivalent** | — |

## BLUE weights vs N_top

| N | σ_END | σ_TOP | ρ | w_END | w_TOP | σ_BLUE |
|---|-------|-------|---|-------|-------|--------|
| 4 | 51.0 | 78.2 | 0.15 | 0.74 | 0.26 | 45.5±0.8 |
| 8 | 50.2 | 32.0 | 0.18 | 0.21 | 0.79 | 28.9±0.5 |
| 14 | 51.0 | 18.7 | 0.18 | 0.02 | 0.98 | 18.2±0.3 |
| 20 | 53.7 | 15.2 | 0.16 | 0.013 | 0.987 | 15.2±0.3 |

## Position scan summary (EJ-230, N=20)

Range: 14.64–17.79 ps (3.15 ps peak-to-peak, 21% variation)
Constant fit: C = 15.75 ± 0.08 ps, χ²/ndf = 163/12, p < 0.001

**The variation is real** (not statistical noise). It correlates with
proximity to the nearest TOP SiPM:
- Positions at 5 mm from a SiPM: σ ≈ 14.6–15.0 ps
- Positions at 25 mm from nearest SiPM (gap midpoints): σ ≈ 17–18 ps
- x=0 at 35 mm: σ = 15.2 ps (two SiPMs at ±35 mm help)

## Gaussian fit results (Figure 2)

| | σ_IQR | σ_Gauss | χ²/ndf | p |
|--|-------|---------|--------|---|
| END (m*=8) | 53.68 ps | 53.31 ps | 74.9/60 | 0.09 |
| TOP (k*=3) | 15.20 ps | 14.93 ps | 39.5/15 | 0.001 |
| BLUE | 15.21 ps | 14.80 ps | 43.8/15 | <0.001 |

TOP/BLUE have non-Gaussian tails — IQR is the appropriate metric.
END distribution is Gaussian (sum of many photons → CLT).

## Resolution metric

σ_IQR = (Q75 − Q25) / 1.349

Used throughout. Bootstrap uncertainty: N_boot=400 with RNG seed 2026.

## What the previous analysis got right

- σ_BLUE ≈ 15.2 ps ✓ (confirmed)
- k* = 3 for N=20 ✓
- m* ≈ 8 for N=20 ✓
- w_TOP ≈ 1 for N≥14 ✓
- EJ-230 ≈ EJ-204 at N=20 (now explicitly verified: 0.2σ difference)

## What needed correction / maturation

1. σ_END in EST mode = 53.7 ps (not 50.5 ps): in EST mode some photons
   hit TOP SiPMs and are not available to END SiPMs — slight degradation.
2. Position scan: range is 14.6–17.8 ps (NOT perfectly flat); variation
   is real and must be described as "~21% peak-to-peak" not "flat".
3. EJ-204 vs EJ-230: must be stated as equivalent at N=20 (0.2σ).
4. A/sqrt(Npe) model: fails at small N (N=4: 37.6 ps predicted, 78.2 measured);
   near-field approximation only valid for N≥14.
