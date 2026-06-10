# EXEC_08 photon-budget conclusions

Input validation: 31/31 files passed; each has event_id 0..1999, matching gun_x_mm, and IDs 0..85.

## Numerical conclusions

- Effective attenuation: left 33.14 +/- 0.03 cm; right 33.06 +/- 0.03 cm. These include light extraction through 70 Top windows and are not expected to equal the 160 cm bulk ABSLENGTH.
- Effective propagation velocity from mean all-photon time: left 27.68 +/- 0.27 cm/ns; right 27.67 +/- 0.27 cm/ns.
- Typical nearest-Top var/mean: median 24.81, range 19.14-29.27. N_pe is not forced to Poisson.
- Intrinsic End SUM4 sigma_group= sigma(DeltaT)/sqrt(2): mean 148.28 ps, range 12.10-437.32 ps. No SPTR/electronics jitter.
- sigma(t_avg) mean 142.49 ps, range 93.20-225.36 ps.
- Top analytic estimate, nearest channel: 55.19-60.36 ps; nearest-three group: 33.61-35.55 ps. This is an analytic estimate, not simulated sigma_t.
- Dead channels: none. Mirrored End L/R asymmetry above 5%: 0/31 positions; maximum 2.05%.

## Method notes

- EJ-204 estimate uses sqrt(tau_r*tau_d)=1.122 ns from tau_r=0.7 ns and tau_d=1.8 ns.
- SUM4 leading edge is ported from `analysis/congruent_sum4_timing.C`: normalized double exponential 0.5/5.0 ns, absolute threshold 4.0 PE, no time-walk correction.
- Gaussian overlays are diagnostic only. Poisson-tail plots use log-y and report both Poisson and Gaussian chi2/ndf.
