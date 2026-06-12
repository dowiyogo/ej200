# EXEC_10 apparent-velocity diagnosis

All errors are delete-one-block jackknife errors over event IDs. Fits combine the mirrored left/right End samples.

| estimator | velocity [cm/ns] | chi2/ndf |
|---|---:|---:|
| All-photon mean (tail-biased) | 25.694 +/- 0.005 | 3357.53 |
| All-photon median t50 | 25.714 +/- 0.005 | 4143.25 |
| Mean FPT | 26.922 +/- 0.003 | 742.63 |

Mean FPT is not consistent with c/n ~= 19 cm/ns within 3 sigma.
The requested t>10 ns fraction test refutes the proposed trend: it increases from 0.002493 near the End to 0.076465 at 139 cm, because a fixed absolute-time threshold also tracks propagation delay.
None of the three global linear fits is statistically acceptable; the quoted slope-derived velocities are descriptive apparent velocities, not propagation measurements.

Alternative hypotheses to test in a future dedicated audit (not changed here): record optical-photon creation time and path length; inspect Geant4's effective GROUPVEL generated from the SSLG4 RINDEX table; and fit local distance ranges instead of forcing one line over 1--139 cm.
The source RINDEX table is constant at 1.58 from 200--800 nm, so the apparent FPT value above c/n is not explained by the documented phase index alone.
