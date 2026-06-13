# EXEC_08 photon-budget conclusions

Input validation: 31/31 files passed; each has event_id 0..1999, matching gun_x_mm, and IDs 0..85.

## Numerical conclusions

- Effective attenuation: left 30.70 +/- 0.02 cm; right 30.82 +/- 0.02 cm. These include light extraction through 70 Top windows and are not expected to equal the 120 cm bulk ABSLENGTH.
- Apparent slope-derived velocity from mean all-photon time: left 27.81 +/- 0.26 cm/ns; right 27.82 +/- 0.26 cm/ns. EXEC_10 shows this is not a direct propagation velocity.
- Typical nearest-Top var/mean: median 22.00, range 18.94-26.54. N_pe is not forced to Poisson.
- Intrinsic End SUM4 sigma_group= sigma(DeltaT)/sqrt(2): mean 176.71 ps, range 11.86-615.15 ps. No SPTR/electronics jitter.
- sigma(t_avg) mean 172.36 ps, range 93.69-325.87 ps.
- Top analytic estimate, nearest channel: 44.83-48.70 ps; nearest-three group: 27.35-28.95 ps. This is an analytic estimate, not simulated sigma_t.
- Dead channels: none. Mirrored End L/R asymmetry above 5%: 0/31 positions; maximum 1.59%.

## Method notes

- EJ-230 estimate uses sqrt(tau_r*tau_d)=0.866 ns from tau_r=0.5 ns and tau_d=1.5 ns.
- SUM4 leading edge is ported from `analysis/congruent_sum4_timing.C`: normalized double exponential 0.5/5.0 ns, absolute threshold 4.0 PE, no time-walk correction.
- Gaussian overlays are diagnostic only. Poisson-tail plots use log-y and report both Poisson and Gaussian chi2/ndf.
