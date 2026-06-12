# EXEC_10 Landau/Moyal N_pe audit

The fitted high-N_pe model is scipy.stats.moyal convolved numerically with a Gaussian.
Moyal is a practical analytic approximation to Landau, not the exact ROOT Landau density.

- Global F fit: c=0.0625472 +/- 0.0001170; chi2/ndf=1.41 over 2666 channel-position points.
- Moyal-Gaussian fits requested for mean >=20 N_pe: 257; failed: 0.
- High-N_pe fits (mean >50) with chi2/ndf >5: 31.
- Fits with sigma_Gauss at the 0.1 N_pe lower bound: 35; sigma_Gauss is not identifiable there, while MPV and Landau width remain reported.

The CSV is intended as future input for comparison with test-beam Landau MPVs through the pending ToT(NPE) calibration.
