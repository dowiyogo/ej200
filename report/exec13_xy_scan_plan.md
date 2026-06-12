# EXEC_13 staged X-Y scan proposal

No simulation is executed by this plan.

## Minimal code changes

Add `/muon/gunY`, preserve the current default `0 mm`, set the primary vertex Y,
record `gun_y_mm` beside `gun_x_mm`, and add unit tests plus CTest coverage for
command parsing, ROOT schema, vertex values, and rejection outside the 60 mm bar.

## Stage A: observability

- X: `-600, -300, 0, +300, +600 mm`
- Y: `-25, -15, -5, 0, +5, +15, +25 mm`
- Start with 500 events per point: 35 points, 17,500 events.
- Require trajectories to remain inside `|Y| < 30 mm`; retain a 5 mm edge margin.

EXEC_07 files contain roughly 9.7–12.2 million hits and 0.62–0.77 GB per 2000
events. Stage A therefore projects roughly 90 million hits and 5.5 GB. Measure
actual CPU time on the first three points before submitting the remainder.

Evaluate held-out Y planes and held-out X columns. Proceed only if mirror-channel
profiles and centroids change monotonically with Y beyond their event spread.

## Stage B: focused scan

Run a finer Y grid only in X regions where Stage A demonstrates useful sensitivity.
Keep complete X or Y slices out of training to prevent positional leakage. Compare
calibrated centroids, covariance-aware combinations, and simple interpretable models
before considering higher-capacity regressors.
