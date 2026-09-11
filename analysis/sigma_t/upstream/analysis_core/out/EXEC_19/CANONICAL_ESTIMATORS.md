# CANONICAL_ESTIMATORS.md — EXEC_20 update
Updated: 2026-06-20 22:34

## Convention (EXEC_20 final)
- **Binning**: sqrt_n (ALL sessions from EXEC_20 onward)
- Walk: parametric α+β/√s, anchored at median(NPE)
- Fit: MAD·1.4826 seeded, window ±2σ_MAD, N_boot=200
- **σ_int**: intrinsic ONLY (no SPTR/FastIC). SPTR/FastIC in deferred appendix.
- σ_int includes 20 ps SiPMSD jitter (set to 0 in actual macs — both datasets)

## Headline estimators @ x=0 (EndTop data, sqrt_n, walk-corr)
- END-only t_avg (all 8 per face): σ = ~882 ps (T1; center is photon-starved)
- TOP nearest N=1: σ ≈ 106 ps (T4 EXEC_19, sqrt_n, walk-corr)
- TOP_SUM4_N1: σ ≈ 69-73 ps (EXEC_16 raw/EXEC_19 walk-corr; see T5 for dynamic vs fixed)
- dynamic nearest-4: σ = 69.4 ps (sqrt_n, walk-corr)
- id//4 fixed cluster: σ = 76.4 ps (sqrt_n, walk-corr)

## Headline estimators @ x=-690 mm (EndTop data)
- END_L_SUM4 co-localized (4 END_L gids 3,4,2,5): σ ≈ 30.5 ps (EXEC_19 T4)
- TOP nearest @ x=-690: σ ≈ 118 ps (EXEC_19 T4)

## END-only baseline (T1, EXEC_20)
- t_avg = (t_L + t_R)/2 with both ends (all 8 per face), walk-corr, sqrt_n
- Position-independent by design; requires both ends to fire
- At center: ~882 ps (few photons); near end x=-690: ~609 ps (tR rarely fires)
- Single-end (inv-var / best face): ~30-44 ps near-end; NOT used for T1 baseline

## Deferred to SPTR/FastIC appendix
- σ_tot_A, σ_tot_B, N_eff, Chain A/B: see EXEC_19 for these numbers
- Not in headline (René's decision, EXEC_20 reorientation)