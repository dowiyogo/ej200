# EXEC_46 Step 5 — revised chain-rule identification and localized residual

Date: 2026-09-16. Revision E1–E5 supersedes the Step 5 conclusion in commit `185a916`.

**CHAIN_RULE_WITHIN_SLOPE_REFUTED. Steps 5.3–5.5 completed; Step 5.2 cancelled; STOP BEFORE STEP 6.**

No simulation, production rerun, source-tree regeneration, push, merge, or deck edit. Input: 210,000 existing derived events, 21 cells, 10,000 per cell. The input SHA-256 is checked unchanged at exit.

Under the linear seven-point between specification, the localized positive structure is present separately in both mirrors for all materials: at ±500 mm its significance is 7.49–10.53 cell SEM, or 5.74–8.18 paired-bootstrap SE after refitting. Its even amplitudes are 6.368 ± 0.900, 7.858 ± 0.920, and 7.338 ± 0.859 ps. The same-projection source-mixture terms are +4.157 ± 0.298, +3.767 ± 0.263, and +2.800 ± 0.216 ps: a partial descriptive allocation, leaving +2.212 ± 0.864, +4.091 ± 0.877, and +4.538 ± 0.827 ps. The mechanism is not fully closed under that specification. G1 below shows that a quadratic between prediction reduces the common amplitude to a value compatible with zero.

## E1–E2: correction of the circular claim

The EJ-200 center-to-|x|=650 secant is -0.02074278 ps/pe, versus the seven-mean OLS slope -0.02105100 ps/pe. Regressing the outcome cell means against their Npe means fits the position curve that was to be explained. The former 96.997%, 97.653%, and 97.294% reductions are **reabsorption into a parameter fitted to the same data**, not explained fractions, causal evidence, or proof of an artifact. The former majority-artifact status is withdrawn.

Step 5.2 is cancelled: with cell-centered Npe, alpha_i equals the observed mean T0_i; a global reference simply rewrites the original residual. It is not an independent observable.

Define the descriptive target r_i = mean(T0)_i − mean(T0)_0 − beta_between[mean(Npe)_i − mean(Npe)_0]. Its even component is the average of the two mirrored residuals. No physical correction is implied.

## E5: conditional slopes across materials

| Material | beta within [ps/pe] | HC1 SE | within / SE magnitude | beta between [ps/pe] | OLS SE |
|---|---|---|---|---|---|
| EJ-200 | -0.056936 | 0.000756 | 75.305116 | -0.021051 | 0.001992 |
| EJ-204 | -0.054226 | 0.000802 | 67.621521 | -0.008742 | 0.002080 |
| EJ-230 | -0.056609 | 0.000894 | 63.312221 | -0.000713 | 0.001981 |

The inverse-variance common within slope is -0.055921 ± 0.000469 ps/pe. Exact equality gives chi2/ndf = 6.860/2, p = 0.03238. The range relative to the common magnitude is 4.85%. The conditional response is nearly common, but exact equality is mildly disfavored (p=0.032); do not claim exact universality. Errors are event-level HC1 fixed-effect slope errors; different material simulations are treated as independent.

The ratio of the EJ-200 to EJ-230 between-slope central magnitudes is 29.54; the EJ-230 denominator is compatible with zero, so “factor 30” is a central-value description, not a precisely measured physical factor. The between fit spans positions and source populations; the within fit conditions on position. Their discrepancy rejects using the measured within slope as the between-position response; it does not independently identify the mechanism of the discrepancy.

## E3: the localized target, both mirrors and uncertainty

The following reports r/SEM exactly as requested and a separate r/bootstrap-SE. The latter uses 2000 ordinary paired-event bootstrap replicates per cell, seed 460516, recomputing the center, seven-point between slope, and residual in each replicate. Both ENDs and all per-event counts remain paired; cells are resampled independently. These z values are descriptive standardized effects, not independent hypothesis tests after model selection. The bootstrap does not cover geometry/model systematics or between-fit misspecification.

| Material | x [mm] | r [ps] | cell SEM [ps] | r/SEM | paired SE [ps] | r/paired SE | 95% lower [ps] | 95% upper [ps] |
|---|---|---|---|---|---|---|---|---|
| EJ-200 | -500 | 6.6052 | 0.8138 | 8.1165 | 1.0544 | 6.2645 | 4.4710 | 8.5649 |
| EJ-200 | 500 | 6.1317 | 0.8186 | 7.4903 | 1.0681 | 5.7406 | 4.0926 | 8.2108 |
| EJ-204 | -500 | 7.7642 | 0.8221 | 9.4441 | 1.0782 | 7.2010 | 5.6928 | 9.8323 |
| EJ-204 | 500 | 7.9528 | 0.8213 | 9.6837 | 1.0885 | 7.3063 | 5.8236 | 10.0376 |
| EJ-230 | -500 | 6.2003 | 0.7923 | 7.8260 | 1.0114 | 6.1303 | 4.2217 | 8.1787 |
| EJ-230 | 500 | 8.4767 | 0.8052 | 10.5278 | 1.0365 | 8.1784 | 6.4042 | 10.5195 |

All 21 signed residuals and their uncertainties are in `residual_significance.csv`; the plot `residual_mirrors.pdf` retains both mirrors.

| Material | \|x\| [mm] | even r [ps] | paired SE [ps] |
|---|---|---|---|
| EJ-200 | 200 | 0.3818 | 0.9289 |
| EJ-200 | 500 | 6.3685 | 0.8996 |
| EJ-200 | 650 | 0.4743 | 0.5918 |
| EJ-204 | 200 | 1.0591 | 0.8930 |
| EJ-204 | 500 | 7.8585 | 0.9203 |
| EJ-204 | 650 | 0.9462 | 0.5724 |
| EJ-230 | 200 | 1.7327 | 0.8516 |
| EJ-230 | 500 | 7.3385 | 0.8593 |
| EJ-230 | 650 | 1.3644 | 0.5467 |

Across materials the common even amplitude is 7.183 ± 0.515 ps (chi2/ndf=1.391/2, p=0.499). This is the result under the linear seven-point between specification.

### G1: specification of the between prediction

The previous F1 table varied the within-cell beta chain and therefore quantified the registered within-chain remnant, not the approximately 7 ps between-fit residual. That table is retained below as a separate result. To test the corrected residual itself, fit the same seven cell means with unweighted linear and quadratic models of mean(T0) versus mean(Npe). In both cases subtract the fitted change relative to the x=0 cell, then average the two mirrors. The bootstrap refits the selected polynomial in every replica.

| Material | linear residual [ps] | linear bootstrap SE [ps] | pol2 residual [ps] | pol2 bootstrap SE [ps] | envelope low [ps] | envelope high [ps] | envelope contains zero |
|---|---|---|---|---|---|---|---|
| EJ-200 | 6.3685 | 0.8996 | -0.5641 | 0.5044 | -0.5641 | 6.3685 | True |
| EJ-204 | 7.8585 | 0.9203 | -0.2706 | 0.4922 | -0.2706 | 7.8585 | True |
| EJ-230 | 7.3385 | 0.8593 | 0.2458 | 0.4795 | 0.2458 | 7.3385 | False |

| between specification | common residual [ps] | common bootstrap SE [ps] | chi2/ndf | p value | common envelope low [ps] | common envelope high [ps] | envelope contains zero |
|---|---|---|---|---|---|---|---|
| linear | 7.1834 | 0.5150 | 0.6957 | 0.4987 | -0.1826 | 7.1834 | True |
| pol2 | -0.1826 | 0.2839 | 0.7011 | 0.4960 | -0.1826 | 7.1834 | True |

The quadratic between prediction reduces all three 500 mm residuals to values near zero. Accordingly, the 7.183 ps common residual does not survive its own linear-versus-quadratic specification envelope. The apparent material universality remains true conditionally within each specification, but its nonzero magnitude is not specification-stable. Neither fit is a causal explanation: both regress the outcome to be explained on a position-correlated cell mean.

### H1: nested-model and out-of-sample tests

The requested primary F-test applies the same unweighted between-fit convention to all seven signed positions. Linear has two parameters and five residual degrees of freedom; pol2 has three parameters and four residual degrees of freedom.

| Material | linear RSS [ps2] | pol2 RSS [ps2] | delta RSS [ps2] | F(1,4) | F-test p |
|---|---|---|---|---|---|
| EJ-200 | 52.577261 | 1.458408 | 51.118853 | 140.204538 | 0.000291 |
| EJ-204 | 74.074828 | 2.529924 | 71.544904 | 113.117874 | 0.000443 |
| EJ-230 | 59.797872 | 4.741181 | 55.056691 | 46.449773 | 0.002423 |

This seven-position F-test rejects the linear specification for all three materials.

The even-component nested test uses exactly the three informative shifts at |x|=200, 500, and 650 mm. The center is fixed to zero by subtraction. The linear model has one coefficient, the quadratic model two, leaving only one residual degree of freedom for pol2. The GLS covariance includes the shared x=0 uncertainty and the two mirror SEMs. The F statistic is `[(chi2_linear-chi2_pol2)/1]/(chi2_pol2/1)`; its low-denominator-dof p value is the conservative nested-model result. The delta-chi2 p value treats the supplied cell-mean uncertainties as calibrated.

| Material | linear chi2 | pol2 chi2 | delta chi2 | F(1,1) | F-test p | delta-chi2 p |
|---|---|---|---|---|---|---|
| EJ-200 | 79.667885 | 1.293377 | 78.374508 | 60.596788 | 0.081336 | 8.525e-19 |
| EJ-204 | 110.372285 | 0.336130 | 110.036155 | 327.361548 | 0.035150 | 9.622e-26 |
| EJ-230 | 90.935587 | 0.280255 | 90.655332 | 323.474383 | 0.035360 | 1.710e-21 |

The even-component F-test rejects the linear model at 5% for EJ-204 and EJ-230; EJ-200 is suggestive but does not cross 5% because the denominator has one degree of freedom. The calibrated delta-chi2 test strongly favors pol2 for all three materials.

LOO uses the seven signed positions exactly as requested. Each row fits six cell means without weights and predicts the excluded seventh; prediction errors are observed minus predicted.

| Material | excluded x [mm] | model | LOO error [ps] | absolute error [ps] | error / cell SEM |
|---|---|---|---|---|---|
| EJ-200 | -650 | linear | -3.3961 | 3.3961 | -4.2023 |
| EJ-200 | -650 | pol2 | -0.8218 | 0.8218 | -1.0169 |
| EJ-200 | -500 | linear | 5.3005 | 5.3005 | 6.5134 |
| EJ-200 | -500 | pol2 | 0.6412 | 0.6412 | 0.7879 |
| EJ-200 | -200 | linear | -1.7735 | 1.7735 | -2.2935 |
| EJ-200 | -200 | pol2 | -0.0760 | 0.0760 | -0.0983 |
| EJ-200 | 0 | linear | -2.9272 | 2.9272 | -3.8057 |
| EJ-200 | 0 | pol2 | 1.1725 | 1.1725 | 1.5244 |
| EJ-200 | 200 | linear | -2.7613 | 2.7613 | -3.5024 |
| EJ-200 | 200 | pol2 | -1.0654 | 1.0654 | -1.3514 |
| EJ-200 | 500 | linear | 4.7481 | 4.7481 | 5.8001 |
| EJ-200 | 500 | pol2 | -0.2760 | 0.2760 | -0.3371 |
| EJ-200 | 650 | linear | -2.3908 | 2.3908 | -2.9363 |
| EJ-200 | 650 | pol2 | 0.7931 | 0.7931 | 0.9741 |
| EJ-204 | -650 | linear | -2.4032 | 2.4032 | -2.8346 |
| EJ-204 | -650 | pol2 | 1.3149 | 1.3149 | 1.5509 |
| EJ-204 | -500 | linear | 5.7735 | 5.7735 | 7.0226 |
| EJ-204 | -500 | pol2 | -0.0855 | 0.0855 | -0.1040 |
| EJ-204 | -200 | linear | -1.2137 | 1.2137 | -1.6296 |
| EJ-204 | -200 | pol2 | 0.9364 | 0.9364 | 1.2573 |
| EJ-204 | 0 | linear | -3.9849 | 3.9849 | -5.5008 |
| EJ-204 | 0 | pol2 | 0.5570 | 0.5570 | 0.7689 |
| EJ-204 | 200 | linear | -3.5340 | 3.5340 | -4.7732 |
| EJ-204 | 200 | pol2 | -1.4879 | 1.4879 | -2.0097 |
| EJ-204 | 500 | linear | 5.9933 | 5.9933 | 7.2977 |
| EJ-204 | 500 | pol2 | 0.2698 | 0.2698 | 0.3285 |
| EJ-204 | 650 | linear | -4.4095 | 4.4095 | -5.1860 |
| EJ-204 | 650 | pol2 | -1.3143 | 1.3143 | -1.5457 |
| EJ-230 | -650 | linear | -3.9018 | 3.9018 | -4.5877 |
| EJ-230 | -650 | pol2 | -1.0768 | 1.0768 | -1.2661 |
| EJ-230 | -500 | linear | 3.7576 | 3.7576 | 4.7429 |
| EJ-230 | -500 | pol2 | -2.2865 | 2.2865 | -2.8861 |
| EJ-230 | -200 | linear | -2.8369 | 2.8369 | -4.0986 |
| EJ-230 | -200 | pol2 | -0.9767 | 0.9767 | -1.4111 |
| EJ-230 | 0 | linear | -4.2029 | 4.2029 | -6.2012 |
| EJ-230 | 0 | pol2 | -0.4783 | 0.4783 | -0.7057 |
| EJ-230 | 200 | linear | -0.5343 | 0.5343 | -0.7621 |
| EJ-230 | 200 | pol2 | 1.4459 | 1.4459 | 2.0624 |
| EJ-230 | 500 | linear | 6.4146 | 6.4146 | 7.9668 |
| EJ-230 | 500 | pol2 | 2.1772 | 2.1772 | 2.7040 |
| EJ-230 | 650 | linear | -2.0000 | 2.0000 | -2.3392 |
| EJ-230 | 650 | pol2 | 1.0954 | 1.0954 | 1.2812 |

| Material | model | LOO RMSE [ps] | LOO MAE [ps] | PRESS [ps2] | standardized RMSE | RMSE reduction vs linear |
|---|---|---|---|---|---|---|
| EJ-200 | linear | 3.5300 | 3.3282 | 87.2284 | 4.3811 | 0.0000 |
| EJ-200 | pol2 | 0.7844 | 0.6923 | 4.3074 | 0.9912 | 0.7778 |
| EJ-204 | linear | 4.2141 | 3.9017 | 124.3130 | 5.2536 | 0.0000 |
| EJ-204 | pol2 | 0.9955 | 0.8523 | 6.9374 | 1.2606 | 0.7638 |
| EJ-230 | linear | 3.7889 | 3.3783 | 100.4882 | 4.9036 | 0.0000 |
| EJ-230 | pol2 | 1.4929 | 1.3624 | 15.6013 | 1.9134 | 0.6060 |

Across all materials, LOO PRESS decreases from 312.030 to 26.846 ps2, a 91.4% reduction. Pol2 predicts better out of sample at every material and every aggregate error measure. Under the H1 decision rule this is evidence against simple overfitting: absorption of the 7.183 ps linear residual by pol2 is supported. This establishes misspecification of the linear between curve; it does not establish that the underlying physical residual is exactly zero.

The residual ordering is EJ-204 > EJ-230 > EJ-200, whereas the same-projection mixture term orders EJ-200 > EJ-204 > EJ-230. The opposite rank pattern rules out the measured mixture term as the dominant origin of the material-independent target. A single common residual gives chi2/ndf = 1.391/2 (p=0.499), so the three amplitudes are statistically consistent with one material-independent value. The mixture remains a smaller additive descriptive component, not the cause of that common structure. The fitted values and both ranks are in `material_common_fit.csv`.

The retained quadratic summaries of this descriptive curve are:

| Material | a2 [ps/m2] | paired a2 SE | diagonal chi2/ndf |
|---|---|---|---|
| EJ-200 | 3.6968 | 0.4226 | 15.2662 |
| EJ-204 | 5.3299 | 0.5247 | 21.1566 |
| EJ-230 | 5.3375 | 0.5625 | 18.0156 |

These chi2/ndf values use the historical cell-SEM weighting (7 points, 5 nominal degrees of freedom); they are lack-of-shape diagnostics, not calibrated goodness-of-fit probabilities: center subtraction and the fitted between slope correlate residuals. The high values preclude interpreting a2 as an adequate shape description.

## 5.3: profile model and binning sensitivity

All 21 original uniform 40-bin TProfile pol1 slopes and their integrated chain reproduce Step 2. The four variants form a 2×2 comparison: uniform/quantile binning × pol1/pol2. Quantiles use 40 equal-population nominal bins (duplicate boundaries removed), extrema padded by 0.5 PE, and standard ROOT TProfile bin-center abscissae and errors. For pol2, beta_i = p1 + 2*p2*mean(Npe)_i. This local-derivative convention is an explicit analysis assumption; it is not a nonlinear transport law. All coefficients, derivative errors, chi2/ndf, fit status and covariance status are in `profile_refits.csv`; the actual profiles, fitted functions and covariance matrices are saved in `profile_refits.root`.

| Variant | chi2/ndf > 5 (of 21) | median chi2/ndf | maximum chi2/ndf |
|---|---|---|---|
| uniform_pol1 | 12 | 5.5956 | 33.2485 |
| uniform_pol2 | 1 | 1.7955 | 5.3022 |
| quantile_pol1 | 17 | 6.3666 | 13.5186 |
| quantile_pol2 | 0 | 1.3257 | 1.9093 |

Pol2 reduces failures from 12/21 to 1/21 with uniform bins; quantile pol2 has 0/21 failures. The large coefficient changes therefore cannot be dismissed as numerical fit failure alone. Profile specification materially affects the within-chain remnant.

Integrate the even measured count profile trapezoidally with these local slopes; no attenuation-model refit. The residual is observed even shift minus that chain. It is distinct from the descriptive between-fit target. Changing local beta cannot alter r_between by definition, so merely showing the original 6–8 ps residual again would be a tautological robustness test.

| Material | Variant | chain residual at 500 [ps] | change vs original [ps] |
|---|---|---|---|
| EJ-200 | uniform_pol1 | 37.9821 | 0.0000 |
| EJ-200 | uniform_pol2 | 61.1193 | 23.1373 |
| EJ-200 | quantile_pol1 | 29.9748 | -8.0073 |
| EJ-200 | quantile_pol2 | 62.9319 | 24.9498 |
| EJ-200 | event_ols_total | 48.9067 | 10.9246 |
| EJ-200 | event_ols_two_counts | 48.6153 | 10.6332 |
| EJ-204 | uniform_pol1 | 62.0571 | 0.0000 |
| EJ-204 | uniform_pol2 | 84.2531 | 22.1960 |
| EJ-204 | quantile_pol1 | 49.0583 | -12.9988 |
| EJ-204 | quantile_pol2 | 85.3294 | 23.2723 |
| EJ-204 | event_ols_total | 69.2361 | 7.1790 |
| EJ-204 | event_ols_two_counts | 68.3377 | 6.2806 |
| EJ-230 | uniform_pol1 | 53.7197 | 0.0000 |
| EJ-230 | uniform_pol2 | 92.8754 | 39.1556 |
| EJ-230 | quantile_pol1 | 58.2526 | 4.5329 |
| EJ-230 | quantile_pol2 | 94.9641 | 41.2443 |
| EJ-230 | event_ols_total | 77.2709 | 23.5512 |
| EJ-230 | event_ols_two_counts | 77.0691 | 23.3494 |

The retained within-chain sensitivity uses exactly the three specifications requested in the original 5.3 contract: uniform-bin pol1 (linear nominal), uniform-bin pol2, and quantile-bin pol1. Quantile pol2 remains a documented 2×2 diagnostic but is not added to the declared three-model envelope. The envelope is quoted asymmetrically around the nominal result; it is a model-specification range, not a Gaussian standard deviation.

| Material | \|x\| [mm] | linear [ps] | pol2 [ps] | quantiles [ps] | envelope low [ps] | envelope high [ps] | nominal/max envelope deviation | zero excluded |
|---|---|---|---|---|---|---|---|---|
| EJ-200 | 200 | 5.6263 | 9.6177 | 4.9366 | 4.9366 | 9.6177 | 1.4096 | True |
| EJ-200 | 500 | 37.9821 | 61.1193 | 29.9748 | 29.9748 | 61.1193 | 1.6416 | True |
| EJ-200 | 650 | 50.8619 | 97.2349 | 37.6613 | 37.6613 | 97.2349 | 1.0968 | True |
| EJ-204 | 200 | 10.3669 | 13.0375 | 7.3753 | 7.3753 | 13.0375 | 3.4653 | True |
| EJ-204 | 500 | 62.0571 | 84.2531 | 49.0583 | 49.0583 | 84.2531 | 2.7959 | True |
| EJ-204 | 650 | 95.3351 | 139.7802 | 72.2622 | 72.2622 | 139.7802 | 2.1450 | True |
| EJ-230 | 200 | 9.0303 | 15.1625 | 9.6078 | 9.0303 | 15.1625 | 1.4726 | True |
| EJ-230 | 500 | 53.7197 | 92.8754 | 58.2526 | 53.7197 | 92.8754 | 1.3720 | True |
| EJ-230 | 650 | 82.6933 | 161.0172 | 95.5223 | 82.6933 | 161.0172 | 1.0558 | True |

This table is the G1 registered-remnant result, separate from the between-fit re-specification above. At |x|=500 the nominal residuals and conservative envelope uncertainties are 37.98 +23.14/−8.01 ps (EJ-200), 62.06 +22.20/−13.00 ps (EJ-204), and 53.72 +39.16/−0.00 ps (EJ-230). Across materials and the three declared specifications, the central values span 29.97–92.88 ps. This registered remnant is therefore not a well-defined quantity. Treating the largest one-sided excursion as a one-sigma-equivalent sensitivity scale gives only 1.64, 2.80, and 1.37 envelope units: none reaches 3. The sign is robust because all three specification values are positive, but the former fit-error significance is withdrawn. Since the envelope has no sampling distribution, these ratios are diagnostics rather than statistical z scores.

To test localized shape in the sensitivity curves, additionally define B = r(500) − [(1−w)r(200)+w*r(650)], w=[N(500)−N(200)]/[N(650)−N(200)]. This declared diagnostic removes a broad trend linear in the measured Npe profile. For r_between, B is exactly the same as for the raw mean curve: a constant beta*N term cancels. Compute the chain variant independently on each mirror (center → 200 → 500 → 650) and on the even curve. This is a shape sensitivity, not a new causal observable or an interpolated physical boundary.

| Material | Variant | mirror | target B [ps] | target B SE | chain B [ps] |
|---|---|---|---|---|---|
| EJ-200 | uniform_pol1 | even | 5.9538 | 0.6928 | 16.2791 |
| EJ-200 | uniform_pol1 | negative | 6.0585 | 0.9410 | 16.0494 |
| EJ-200 | uniform_pol1 | positive | 5.8513 | 0.9988 | 16.5125 |
| EJ-200 | uniform_pol2 | even | 5.9538 | 0.6928 | 20.3626 |
| EJ-200 | uniform_pol2 | negative | 6.0585 | 0.9410 | 20.2363 |
| EJ-200 | uniform_pol2 | positive | 5.8513 | 0.9988 | 20.4918 |
| EJ-200 | quantile_pol1 | even | 5.9538 | 0.6928 | 13.4079 |
| EJ-200 | quantile_pol1 | negative | 6.0585 | 0.9410 | 12.2930 |
| EJ-200 | quantile_pol1 | positive | 5.8513 | 0.9988 | 14.5322 |
| EJ-200 | quantile_pol2 | even | 5.9538 | 0.6928 | 20.7971 |
| EJ-200 | quantile_pol2 | negative | 6.0585 | 0.9410 | 20.7216 |
| EJ-200 | quantile_pol2 | positive | 5.8513 | 0.9988 | 20.8752 |
| EJ-200 | event_ols_total | even | 5.9538 | 0.6928 | 17.7728 |
| EJ-200 | event_ols_total | negative | 6.0585 | 0.9410 | 17.9137 |
| EJ-200 | event_ols_total | positive | 5.8513 | 0.9988 | 17.6333 |
| EJ-200 | event_ols_two_counts | even | 5.9538 | 0.6928 | 18.2031 |
| EJ-200 | event_ols_two_counts | negative | 6.0585 | 0.9410 | 18.3657 |
| EJ-200 | event_ols_two_counts | positive | 5.8513 | 0.9988 | 18.0420 |
| EJ-204 | uniform_pol1 | even | 6.8396 | 0.6962 | 21.3988 |
| EJ-204 | uniform_pol1 | negative | 5.9913 | 0.9721 | 20.4380 |
| EJ-204 | uniform_pol1 | positive | 7.6869 | 0.9854 | 22.3579 |
| EJ-204 | uniform_pol2 | even | 6.8396 | 0.6962 | 26.0315 |
| EJ-204 | uniform_pol2 | negative | 5.9913 | 0.9721 | 25.4848 |
| EJ-204 | uniform_pol2 | positive | 7.6869 | 0.9854 | 26.5772 |
| EJ-204 | quantile_pol1 | even | 6.8396 | 0.6962 | 18.5506 |
| EJ-204 | quantile_pol1 | negative | 5.9913 | 0.9721 | 17.5819 |
| EJ-204 | quantile_pol1 | positive | 7.6869 | 0.9854 | 19.5178 |
| EJ-204 | quantile_pol2 | even | 6.8396 | 0.6962 | 26.2353 |
| EJ-204 | quantile_pol2 | negative | 5.9913 | 0.9721 | 24.7896 |
| EJ-204 | quantile_pol2 | positive | 7.6869 | 0.9854 | 27.6781 |
| EJ-204 | event_ols_total | even | 6.8396 | 0.6962 | 22.9594 |
| EJ-204 | event_ols_total | negative | 5.9913 | 0.9721 | 22.0724 |
| EJ-204 | event_ols_total | positive | 7.6869 | 0.9854 | 23.8448 |
| EJ-204 | event_ols_two_counts | even | 6.8396 | 0.6962 | 23.1614 |
| EJ-204 | event_ols_two_counts | negative | 5.9913 | 0.9721 | 22.6254 |
| EJ-204 | event_ols_two_counts | positive | 7.6869 | 0.9854 | 23.6917 |
| EJ-230 | uniform_pol1 | even | 5.7362 | 0.6763 | 18.5897 |
| EJ-230 | uniform_pol1 | negative | 5.3341 | 0.9405 | 17.3620 |
| EJ-230 | uniform_pol1 | positive | 6.1387 | 0.9733 | 19.8202 |
| EJ-230 | uniform_pol2 | even | 5.7362 | 0.6763 | 26.0346 |
| EJ-230 | uniform_pol2 | negative | 5.3341 | 0.9405 | 24.3090 |
| EJ-230 | uniform_pol2 | positive | 6.1387 | 0.9733 | 27.7645 |
| EJ-230 | quantile_pol1 | even | 5.7362 | 0.6763 | 18.2041 |
| EJ-230 | quantile_pol1 | negative | 5.3341 | 0.9405 | 17.8983 |
| EJ-230 | quantile_pol1 | positive | 6.1387 | 0.9733 | 18.5099 |
| EJ-230 | quantile_pol2 | even | 5.7362 | 0.6763 | 27.3986 |
| EJ-230 | quantile_pol2 | negative | 5.3341 | 0.9405 | 26.0490 |
| EJ-230 | quantile_pol2 | positive | 6.1387 | 0.9733 | 28.7514 |
| EJ-230 | event_ols_total | even | 5.7362 | 0.6763 | 23.3212 |
| EJ-230 | event_ols_total | negative | 5.3341 | 0.9405 | 22.2518 |
| EJ-230 | event_ols_total | positive | 6.1387 | 0.9733 | 24.3929 |
| EJ-230 | event_ols_two_counts | even | 5.7362 | 0.6763 | 23.9739 |
| EJ-230 | event_ols_two_counts | negative | 5.3341 | 0.9405 | 22.5332 |
| EJ-230 | event_ols_two_counts | positive | 6.1387 | 0.9733 | 25.4185 |

The localized excess B remains positive separately on both mirrors for all four profile variants: 12.29–28.75 ps. However its magnitude changes strongly; this supports persistence of the localized shape, not invariance of the original 6–8 ps amplitude under a physically identified correction. Even-chain a2 moves from 123.10/227.05/197.23 to 238.29/336.23/384.33 ps/m2 for quantile pol2.

Profile sensitivity values are central estimates; their formal fit errors are retained, but no independent confidence claim is made from a failed profile fit or from these correlated post-fit shape differences.

## 5.4: exact source-mixture finite differences

For each END use the creator source of the actual winning photon: source_type=2 is Cherenkov, source_type=1 is scintillation. mu_C and mu_S are conditional means of the winning END timestamp, not means of source-specific minima in all events. All winner classes and Npe source counts close exactly. The interval identity used on adjacent signed positions is:

`Delta mu = f_bar Delta mu_C + (1−f_bar) Delta mu_S + (mu_C_bar−mu_S_bar) Delta f`

Dividing each term by the measured interval length gives finite-difference derivatives; integrating from x=0 and averaging the ENDs gives T0 terms with exact numerical closure. The interval-midpoint convention is declared and avoids an unreported product-rule discretization error. The three terms and derivative units ps/m, separately by END and signed interval, are in `mixture_intervals.csv`. All signed integrated terms, left/right contributions and paired bootstrap errors are in `mixture_components.csv`. No inference of a pointwise derivative inside the unmeasured intervals is possible.

| Material | \|x\| [mm] | f dmu_C [ps] | (1−f) dmu_S [ps] | mix term [ps] | mix SE [ps] | sum / observed [ps] |
|---|---|---|---|---|---|---|
| EJ-200 | 200 | -3.3150 | 2.2873 | -0.3258 | 0.3183 | -1.3535 |
| EJ-200 | 500 | -29.5895 | 23.9132 | -0.5876 | 0.3168 | -6.2639 |
| EJ-200 | 650 | -179.6839 | 159.5691 | -11.8077 | 0.4447 | -31.9225 |
| EJ-204 | 200 | -3.9196 | 3.9937 | 0.2321 | 0.2778 | 0.3062 |
| EJ-204 | 500 | -28.7536 | 31.0584 | -0.3687 | 0.2856 | 1.9362 |
| EJ-204 | 650 | -155.5207 | 151.2565 | -10.0428 | 0.3359 | -14.3070 |
| EJ-230 | 200 | -3.7862 | 5.4433 | 0.0217 | 0.2385 | 1.6788 |
| EJ-230 | 500 | -29.9157 | 36.9529 | -0.1484 | 0.2367 | 6.8889 |
| EJ-230 | 650 | -136.0443 | 143.5681 | -7.3300 | 0.2673 | 0.1939 |

At |x|=650 the raw mixing term is −11.808/−10.043/−7.330 ps, opposite in sign to the registered within-chain residual. At |x|=500 it is only −0.588/−0.369/−0.148 ps. It therefore does not explain the original positive remnant; the magnitude near −9 ps occurs at 650 mm.

| Material | \|x\| [mm] | original chain residual [ps] | raw mix [ps] | original minus raw mix [ps] |
|---|---|---|---|---|
| EJ-200 | 200 | 5.6263 | -0.3258 | 5.9521 |
| EJ-200 | 500 | 37.9821 | -0.5876 | 38.5696 |
| EJ-200 | 650 | 50.8619 | -11.8077 | 62.6696 |
| EJ-204 | 200 | 10.3669 | 0.2321 | 10.1348 |
| EJ-204 | 500 | 62.0571 | -0.3687 | 62.4258 |
| EJ-204 | 650 | 95.3351 | -10.0428 | 105.3779 |
| EJ-230 | 200 | 9.0303 | 0.0217 | 9.0086 |
| EJ-230 | 500 | 53.7197 | -0.1484 | 53.8681 |
| EJ-230 | 650 | 82.6933 | -7.3300 | 90.0232 |

A raw mixing shift and a detrended target have different definitions. For a like-for-like comparison apply the SAME operator R_N(v)=v−v_0−OLS_slope(v,N)*(N−N_0) to the integrated mixing curve. For fixed measured N this is linear, so R_N(T0)=R_N(mixing)+R_N(T0−mixing) exactly. This is a descriptive allocation of a measured curve; it does not cure circularity or identify a counterfactual effect.

| Material | target [ps] | raw mix [ps] | same-projection mix [ps] | projected mix SE [ps] | signed fraction of target | fraction SE | remaining [ps] | remaining SE [ps] |
|---|---|---|---|---|---|---|---|---|
| EJ-200 | 6.3685 | -0.5876 | 4.1569 | 0.2978 | 0.6527 | 0.0981 | 2.2115 | 0.8642 |
| EJ-204 | 7.8585 | -0.3687 | 3.7672 | 0.2626 | 0.4794 | 0.0588 | 4.0913 | 0.8774 |
| EJ-230 | 7.3385 | -0.1484 | 2.8001 | 0.2157 | 0.3816 | 0.0483 | 4.5383 | 0.8267 |

Under this identical projection the mixture receives 65.3%, 47.9%, and 38.2% of the 500 mm target (signed descriptive fractions, not independent explained fractions). The remaining effect is nonzero, especially for EJ-204/230. The large negative mixing shift at 650 changes the fitted broad trend, which is why its projected contribution at 500 becomes positive despite a small negative raw value.

The fitted cell fraction f_C is fixed as a population summary at each x; the event-level winning source is random and may correlate with photon counts. Therefore it is not justified to assert that a within-cell beta is necessarily blind to every mixture effect. The exact identity above measures the mixture term without that assertion.

## 5.5: separate scintillation and Cherenkov counts

Fit event-level T0 = intercept + beta_S*Nscint_END + beta_C*NCher_END independently in each cell, where each count sums the two ENDs. Use event OLS with HC1 covariance and retain a one-count event OLS control to separate count splitting from the original binned-profile weighting. Average coefficients/counts over mirrors then sum the two trapezoidal beta*dN chains. All 21 coefficient pairs, covariance, count correlations and conditioning are in `two_count_slopes.csv`. This remains a within-to-between transfer test; a second count does not automatically resolve position confounding.

| Material | Variant | \|x\| [mm] | residual [ps] | paired SE [ps] |
|---|---|---|---|---|
| EJ-200 | uniform_pol1 | 200 | 5.6263 | nan |
| EJ-200 | uniform_pol1 | 500 | 37.9821 | nan |
| EJ-200 | uniform_pol1 | 650 | 50.8619 | nan |
| EJ-200 | event_ols_total | 200 | 7.8198 | 0.8965 |
| EJ-200 | event_ols_total | 500 | 48.9067 | 1.1127 |
| EJ-200 | event_ols_total | 650 | 73.4198 | 1.4666 |
| EJ-200 | event_ols_two_counts | 200 | 7.7994 | 0.9080 |
| EJ-200 | event_ols_two_counts | 500 | 48.6153 | 1.3103 |
| EJ-200 | event_ols_two_counts | 650 | 71.4261 | 1.8611 |
| EJ-204 | uniform_pol1 | 200 | 10.3669 | nan |
| EJ-204 | uniform_pol1 | 500 | 62.0571 | nan |
| EJ-204 | uniform_pol1 | 650 | 95.3351 | nan |
| EJ-204 | event_ols_total | 200 | 10.8906 | 0.8577 |
| EJ-204 | event_ols_total | 500 | 69.2361 | 1.2411 |
| EJ-204 | event_ols_total | 650 | 110.1497 | 1.7503 |
| EJ-204 | event_ols_two_counts | 200 | 10.8095 | 0.8659 |
| EJ-204 | event_ols_two_counts | 500 | 68.3377 | 1.4457 |
| EJ-204 | event_ols_two_counts | 650 | 107.2092 | 2.1831 |
| EJ-230 | uniform_pol1 | 200 | 9.0303 | nan |
| EJ-230 | uniform_pol1 | 500 | 53.7197 | nan |
| EJ-230 | uniform_pol1 | 650 | 82.6933 | nan |
| EJ-230 | event_ols_total | 200 | 12.7657 | 0.8132 |
| EJ-230 | event_ols_total | 500 | 77.2709 | 1.2756 |
| EJ-230 | event_ols_total | 650 | 129.0017 | 1.8218 |
| EJ-230 | event_ols_two_counts | 200 | 12.9406 | 0.8311 |
| EJ-230 | event_ols_two_counts | 500 | 77.0691 | 1.4834 |
| EJ-230 | event_ols_two_counts | 650 | 126.2713 | 2.2418 |

Splitting the counts changes event-OLS a2 by only −4.54/−6.75/−6.15 ps/m2 relative to the matched one-count event OLS control; the residual remains +170.61/+255.81/+298.94 ps/m2. The positive localized B remains 18.20/23.16/23.97 ps. This two-count linear transfer does not close the target.

The complete quadratic comparison (residual fit, not separate fitted-curve subtraction) is:

| Material | Variant | a2 [ps/m2] | paired a2 SE | diagonal chi2/ndf | B at 500 [ps] |
|---|---|---|---|---|---|
| EJ-200 | uniform_pol1 | 123.1032 | nan | 24.8400 | 16.2791 |
| EJ-200 | uniform_pol2 | 231.2984 | nan | 5.1981 | 20.3626 |
| EJ-200 | quantile_pol1 | 90.8271 | nan | 22.8692 | 13.4079 |
| EJ-200 | quantile_pol2 | 238.2934 | nan | 5.0970 | 20.7971 |
| EJ-200 | event_ols_total | 175.1448 | 3.2632 | 11.6908 | 17.7728 |
| EJ-200 | event_ols_two_counts | 170.6059 | 4.2649 | 15.7313 | 18.2031 |
| EJ-204 | uniform_pol1 | 227.0532 | nan | 12.4813 | 21.3988 |
| EJ-204 | uniform_pol2 | 332.1716 | nan | 1.7530 | 26.0315 |
| EJ-204 | quantile_pol1 | 174.2082 | nan | 16.2557 | 18.5506 |
| EJ-204 | quantile_pol2 | 336.2260 | nan | 1.6660 | 26.2353 |
| EJ-204 | event_ols_total | 262.5611 | 3.9846 | 7.0603 | 22.9594 |
| EJ-204 | event_ols_two_counts | 255.8137 | 5.0900 | 9.8930 | 23.1614 |
| EJ-230 | uniform_pol1 | 197.2264 | nan | 10.3802 | 18.5897 |
| EJ-230 | uniform_pol2 | 379.4790 | nan | 4.0211 | 26.0346 |
| EJ-230 | quantile_pol1 | 226.4128 | nan | 2.6861 | 18.2041 |
| EJ-230 | quantile_pol2 | 384.3280 | nan | 2.3193 | 27.3986 |
| EJ-230 | event_ols_total | 305.0921 | 4.1960 | 1.9012 | 23.3212 |
| EJ-230 | event_ols_two_counts | 298.9417 | 5.2310 | 3.6534 | 23.9739 |

Quadratic chi2/ndf uses historical cell SEM weights and 5 nominal degrees of freedom for comparison; the integrated predictions also have uncertainty and shared fitted coefficients, so these are shape diagnostics. Bootstrap errors are supplied for event OLS variants; profile variants retain conditional formal errors in the sidecar. No poor quadratic or poor profile fit is used as a mechanism measurement.

## Step 5 closeout

The registered +123.10/+227.05/+197.23 ps/m2 remnant resulted from two compounded specification errors: transferring the within-cell slope to a between-position change, and forcing the between-position mean(T0)-mean(Npe) relation to be linear. The first is rejected at 16.84, 20.41, and 25.72 standard errors. For the second, the seven-position linear-versus-pol2 F tests are 140.20/113.12/46.45 with p=0.000291/0.000443/0.002423, and pol2 reduces LOO RMSE by 77.8/76.4/60.6%. Thus the 7.183 ps common residual is identified as an artifact of the linear between-curve specification.

**Step 5 conclusion:** the T0(x) nonlinearity is compatible with being entirely an Npe response once curvature is allowed and the response is estimated between positions. This does not prove that no additional physical mechanism exists; with only seven positions, this design cannot resolve such a mechanism after the Npe response is re-specified.

## Reproducibility and stopping point

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step5.py
```

Input SHA-256: `816c302c3ed37f7eb29855bcad7e734a1b6a77238f60a30361635b68377e341d`. Bootstrap seed: 460516; replicates: 2000. ROOT version: 6.40.02. Input opened read-only. Every produced PDF has `.root`, `.csv`, `.meta.json` companions. `step5_revision_tables.root` also preserves all table columns, including categorical labels. Rollback tag before edits: `pre-exec46-step5-e1-e5-20260916`.

**Step 6 has not run and still requires René’s explicit approval. No push was performed.**
