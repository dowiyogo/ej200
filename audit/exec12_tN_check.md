# EXEC_12 T4 consistency check: sigma_fit(t_4) vs sigma_group

sigma_group = sigma(Delta_T_AB) / sqrt(2) from exec08b_timing_gate.py.
sigma_fit(t_4) = Gaussian core sigma of the t_4 distribution for the near-End SUM4.

Physical relationship:
  t_4 = time of the 4th individual photon in the cluster (pure photon-counting jitter).
  t_LE = when the analog SUM4 waveform (sum of SPE pulses, tau_r=0.7 ns) crosses 4 PE.
  t_LE > t_4 always; the extra delay has its own jitter from the SPE pulse convolution.
  At high N_pe (x=+/-400, ~360 PE/cluster): sigma(t_4) is small (~15 ps) but
    sigma_group stays larger (~50 ps) due to residual discriminator jitter.
  At moderate N_pe (x=0, ~195 PE/cluster): photon-counting noise dominates and
    sigma(t_4) ≈ sigma_group ≈ 110 ps.

Criterion (spec: 'mismo orden, misma forma vs x'):
  (A) ratio in [0.10, 10.0] (same order of magnitude) for all positions.
  (B) same shape as function of x (both V-shaped, larger at x=0).

| x [mm] | End side | sigma_fit(t_4) [ps] | sigma_group [ps] | ratio |
|--------|----------|---------------------|------------------|-------|
| -400 | left | 14.9 | 50.0 +/- 0.8 | 0.30 [OK] |
| +0 | left | 122.9 | 107.5 +/- 1.7 | 1.14 [OK] |
| +0 | right | 122.9 | 107.2 +/- 1.7 | 1.15 [OK] |
| +400 | right | 14.8 | 51.5 +/- 0.8 | 0.29 [OK] |

## Verdict

CONSISTENCY CHECK PASSED: sigma_fit(t_4) within one order of magnitude of sigma_group at all tested positions; shape vs x is the same (V-shaped). The factor-of-3 residual at x=+/-400 is expected: sigma_group includes analog discriminator jitter from the SPE pulse convolution, which dominates over photon-counting noise at high N_pe.
