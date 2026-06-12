# EXEC_12T timing and position report

## Abstract
The fourth-hit gate reproduces EXEC_11 exactly. Both 4th and 20th detected-hit samples have 100% efficiency, so native and matched comparisons coincide. The 20th hit has a broader differential timing spread (99.8 versus 78.5 ps), but a steeper position slope and better cross-validated X RMS68 (10.32 versus 13.34 mm). Its common-time spread is worse. Higher order statistics continue to improve X RMS68 through k=30 while calibration nonlinearity rises. These are intrinsic simulated hit-order estimator spreads, not experimental detector resolutions and not 20% CFD.

## Definitions
The fourth and twentieth detected hits are order statistics of accepted optical hits.
They are not the 20% constant-fraction discriminator of Lv et al., which reconstructs
a waveform, takes a leading-edge 20% crossing, and combines 64 channels.

## Provenance and 4PE gate
41 positions, 3000 events each, EJ-204/OPSC-101, Top readout, jitter zero, data
commit f431c01. EXEC_11 reproduction: PASS.

## Main comparison
metric,4th,20th
efficiency mean,1.0,1.0
slope ps/mm,5.880337826240657,9.673246069752915
mean delta-t RMS68 ps,78.47006841455072,99.83712557367672
mean CV X RMS68 mm,13.34414286488821,10.320529503460195
max abs bias mm,0.8983456623983983,2.2428394801891938
mean t+ RMS68 ps,47.99224821714844,65.67712730494983


## Window collection
The four EXEC_08b points are explanatory only. NPE rises from roughly
349.4 to 413.2 as the track moves away from
the window-centered dip, while t4/t20 spreads also change. They are not used in calibration.

## Global context
EXEC_12 shows transferable local Top precision but invalidity at the central gap and
large global-model edge biases. Y remains only a y=0 feasibility proxy.

## Comparison with Lv et al.
Lv et al.: EJ-200 200x200x6 mm plate, 64 SiPMs, waveform 20% CFD, NPE weighting,
2016 geometric pairs and CNN. This work: 1.4 m EJ-204 bar, adjacent Top pair,
hit-order statistics, empirical 1-D calibration, zero added jitter. The numerical
resolutions are not directly comparable.

## Limitations
Simulation only; no SPTR, electronics waveform, real CFD, noise, dark counts,
cross-talk, afterpulsing, synchronization uncertainty, or experimental validation.

## Conclusion
The fourth-hit gate reproduces EXEC_11 exactly. Both 4th and 20th detected-hit samples have 100% efficiency, so native and matched comparisons coincide. The 20th hit has a broader differential timing spread (99.8 versus 78.5 ps), but a steeper position slope and better cross-validated X RMS68 (10.32 versus 13.34 mm). Its common-time spread is worse. Higher order statistics continue to improve X RMS68 through k=30 while calibration nonlinearity rises. These are intrinsic simulated hit-order estimator spreads, not experimental detector resolutions and not 20% CFD.
