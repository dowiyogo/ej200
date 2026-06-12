# Speaker notes — EXEC_12T-B Beamer
# Intrinsic Timing and Position Reconstruction with Adjacent SiPMs in EJ-204

---

## Slide 1 — Title
This talk presents the intrinsic timing and position reconstruction for
a pair of adjacent SiPMs (IDs 28 and 29) in a 1.4 m EJ-204 scintillator bar.
All results are simulation predictions with jitter = 0 — no SPTR, no electronics.

---

## Slide 2 — Scientific questions
Introduce the seven questions that structure the analysis. Emphasise that
the results are simulation-intrinsic and will not match a real detector
without adding SPTR and electronics (see Slide 27).

---

## Slide 3 — Detector geometry
The bar is 1.4 m long, 6 cm wide, 1 cm thick. The fine scan covers 41
positions from x = -462 to -422 mm, bracketing the pair centred at -442 mm.
SiPM A is closer to the bar edge; SiPM B is slightly more central.

---

## Slide 4 — Available datasets
Three datasets: the fine scan (main analysis), the global EndTop scan (context),
and the window-dip runs (diagnostic). Emphasise that efficiency = 100% for
both k=4 and k=20 in the fine scan — no efficiency loss is observed because
NPE >> 20 throughout.

---

## Slide 5 — From hits to observables
Each photon that triggers the SiPM produces a hit with a time stamp.
We sort these by arrival time and take the k-th earliest. This is a
pure order statistic — no waveform reconstruction needed.

---

## Slide 6 — Timing estimators
Δt encodes position (the difference in photon path length to A vs B).
t+ encodes the common arrival time. These two are complementary and
their covariance is measured, not assumed.

---

## Slide 7 — 20th hit ≠ 20% CFD (IMPORTANT)
This is the most important conceptual slide. The 20th detected hit is
simply the 20th smallest arrival time — a counting statistic. The 20% CFD
of Lv et al. reconstructs a full SiPM waveform and takes the time at 20%
of its peak amplitude. These are completely different estimators.
Repeat this clearly if there are questions about the comparison.

---

## Slide 8 — Mean Δt vs x — calibrations
The 20PE slope (9.67 ps/mm) is 64% steeper than 4PE (5.88 ps/mm).
This means the position sensitivity of 20PE is much higher. However,
the calibration χ²/ndf of 43.5 vs 3.7 tells us that the linear model
fits 20PE much worse — there is significant non-linearity near the SiPMs.

---

## Slide 9 — Timing width vs position
The U-shape reflects the photon path-length variation. Near the pair
centre, photons arrive more uniformly; at the scan edges, the spread grows.
4PE is consistently narrower (~78 ps vs ~100 ps for 20PE), but this
doesn't automatically translate to better position resolution because
of the slope difference.

---

## Slide 10 — A-B correlation
For k=4, the correlation ρ ≈ 0.02 (near zero), meaning the two channels
measure mostly independent information. For k=20, ρ ≈ 0.31 — the later
hits in both channels are correlated because they tend to come from the
same photon burst (common scintillation fluctuation). This reduces the
independence of the two measurements.

---

## Slide 11 — Common-time estimator width
σ(t+) = 48 ps for 4PE and 66 ps for 20PE. Both are larger than σ(Δt)/2
due to the non-zero A-B correlation. In a real detector with SPTR of
~100-200 ps per channel, these numbers will increase substantially.

---

## Slide 12 — Δt at REF1
REF1 is at -441 mm, the position closest to the pair midpoint (where
Δt ≈ 0). The distributions are approximately Gaussian with the widths
quoted. Note that the widths differ: 4PE is narrower.

---

## Slide 13 — Δt at REF2
REF2 is at -443 mm, selected as the position with maximum NPE in channel A.
Results are very similar to REF1 (both are near the pair midpoint).

---

## Slide 14 — x_rec at REF1
Leave-one-out cross-validation: we reconstruct each event using a
calibration fitted on the other 40 positions. The 4PE distribution is
wider (lower slope) but the mean is near the true position. The 20PE
distribution is narrower (steeper slope), at the cost of larger bias.

---

## Slide 15 — x_rec at REF2
Consistent with REF1. The full bias/sigma table is in backup B4.

---

## Slide 16 — Position resolution vs threshold
The resolution curve has a complex shape. It improves as k increases from
1 to ~4, then worsens, then improves again around k=20 (because the slope
increase dominates), then worsens again at k>25 as the non-linearity takes over.
4PE: 13.3 mm mean LOO RMS68. 20PE: 10.3 mm.

---

## Slide 17 — Slope and χ²/ndf vs threshold
The slope monotonically increases with k. The χ²/ndf blows up after k≈15,
confirming that the linear calibration model becomes inadequate for high
thresholds. The 20PE calibration already has χ²/ndf = 43.5.

---

## Slide 18 — Bias vs threshold
Maximum absolute bias grows with k. For 4PE it is only 0.90 mm;
for 20PE it rises to 2.24 mm. The bias is driven by the non-linearity
at the scan edges near the SiPMs.

---

## Slide 19 — Resolution–bias trade-off
Because efficiency = 100% for all k, the only relevant trade-off is
between position resolution (RMS68) and reconstruction bias. 4PE gives
lower bias at the cost of slightly worse resolution. 20PE gives better
resolution but with higher bias. Both are on the approximate Pareto front.

---

## Slide 20 — Threshold sweep table
Selected rows from the full 30-threshold sweep. Note the RMS68(x)
minimum near k=1–2 (very low, because k=1 is just the first photon with
tiny spread) then rising, then falling at k=20, then rising again at k=30.

---

## Slide 21 — Main 4PE vs 20PE comparison
Read the table row by row. The key finding: 20PE has better mean position
resolution (10.3 vs 13.3 mm) driven by the steeper slope, but it has
worse timing spread, worse calibration linearity, and higher A-B correlation.

---

## Slide 22 — Does 20PE really improve?
The verdict is nuanced: 20PE improves position resolution via slope, but
not via timing precision. It comes with non-linearity and common-mode costs.
In a real detector, the choice would depend on the specific use case.

---

## Slide 23 — Window-dip NPE
The NPE drop at the SiPM centre (run A at x=-652, distance=0) is a
geometric effect: photons entering directly at the SiPM centre are
partially absorbed or reflect differently. Run C1 at 4 mm offset shows
the highest NPE. This is the "window dip" effect.

---

## Slide 24 — Window-dip timing
The mean arrival time t4 varies by ~14 ps across the four runs.
For t20, the variation is larger (~23 ps). More photons (higher NPE)
means earlier mean arrival time — consistent with photon statistics.

---

## Slide 25 — Global EXEC_12 context
The pair (28,29) provides the best local precision (~3 mm for ratio),
but it only covers a 40 mm window of the 1.4 m bar. The global methods
(EndTop timing, centroid) cover the full bar but with much worse precision.
No method covers the central gap. This is the fundamental limitation of
the current detector layout for X reconstruction.

---

## Slide 26 — Comparison with Lv et al.
The key message is in the last row of the table: "Direct numerical comparison
is NOT valid." The detectors differ in geometry (200×200×6 mm vs 1400×60×10 mm),
channel count (64 vs pair), method (waveform CFD vs order statistic), and
the Lv et al. simulation includes SPTR and electronics while ours does not.
Use this slide to explain the physics, not to benchmark.

---

## Slide 27 — Timing budget
This is the honest accounting. The simulation gives us the optical lower bound.
SPTR alone adds ~100 ps per channel (quadrature). Electronics add more.
Our 78 ps σ(Δt) for 4PE will become several hundred ps in a real detector.
This is not a limitation of the analysis — it is the correct scientific framing.

---

## Slide 28 — Conclusions
Read the 7 conclusions. Emphasise (7): all numbers are intrinsic simulation
predictions. The question "what timing and position precision is achievable?"
is answered for the intrinsic case; the experimental case requires all the
absent ingredients listed on Slide 27.

---

## Slide 29 — Next steps
The simulation additions (SPTR, electronics) are the priority. Once those
are in, the comparison with Lv et al. becomes more meaningful. The Y scan
is needed for 2D reconstruction.

---

## Backup slides B1–B10
These contain full data inventory, event builder validation, fit QA,
LOO methodology, window-dip metadata, and reproducibility commands.
Refer to them as needed during questions.
