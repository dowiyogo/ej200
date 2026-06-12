# Event-by-event position reconstruction in a 1.4 m EJ-204 scintillator bar with end and top SiPM readout

## Scope
EXEC_11 established a local result near pair (28,29). EXEC_12 evaluates global X
with 31 leave-one-position-out folds and explores, but does not calibrate, Y at y=0.
All resolutions are **simulation prediction — EJ-204 — intrinsic optical timing**.

## Data and validation
The EndTop scan contains 31 stable ROOT files and 2000 events per position. Logs
confirm OPSC-101/EJ-204 and EndTop. ROOT content confirms gun X. Jitter=0, y=0 and
perpendicular incidence are inferred from the generating script/code.

## Global X results
| method | bias_mean | bias_abs_max | sigma_core_mean | sigma_core_max | valid_mean |
|---|---|---|---|---|---|
| A_end | -0.004 | 135.958 | 14.741 | 31.329 | 1.000 |
| R_end | -0.005 | 55.626 | 21.350 | 27.048 | 1.000 |
| dt_end_pool | -0.005 | 10.310 | 20.611 | 29.925 | 1.000 |
| dt_end_weighted | -0.010 | 18.369 | 37.412 | 68.972 | 1.000 |
| local_R | 0.001 | 4.542 | 3.001 | 3.162 | 0.968 |
| x_top_centroid_raw | 0.000 | 57.886 | 4.379 | 5.625 | 1.000 |

BLUE over End weighted timing, End ratio and Top centroid has mean core spread
`8.956 mm`; its maximum absolute bias is `14.436 mm`.
The local transferable Top-pair ratio has mean core spread
`3.001 mm`, but is invalid
at x=0 because the 24 mm central gap makes pair choice ambiguous. It has position-
dependent biases up to `4.542 mm`.
Thus the ~5.2 mm local EXEC_11 result generalizes in local-pair spread, but not as a
uniformly available unbiased estimator across the full bar. Global linear Top
centroid has narrow event spread but large edge bias; End estimators degrade strongly
at extremes. The current BLUE improves individual global End methods but does not
beat the valid local-pair estimator.

## Y=0 feasibility proxy
| estimator | mean_centroid | max_abs_mean | mean_width | max_width |
|---|---|---|---|---|
| y_L_centroid | -0.256 | 0.896 | 1.395 | 2.778 |
| y_R_centroid | -0.245 | 0.887 | 1.388 | 2.841 |
| y_end_centroid | -0.194 | 0.884 | 0.631 | 0.919 |

At y_true=0, the combined end-channel light-sharing centroid is centered near zero
and has an event-by-event width shown above. This is a central-point estimator spread /
feasibility proxy, **not a Y spatial resolution**.

## Availability matrix
| Capability | Status |
|---|---|
| Local X near pair (28,29) | validated in current simulation |
| Global X along bar | evaluated with position-level LOO CV |
| Y at central incidence | exploratory proxy |
| General Y reconstruction | not yet available |
| Full X–Y reconstruction | requires Y scan |
| Experimental validation | not yet available |

## Conclusion
Current data support X reconstruction and show transferable local Top information,
but edge biases, the central gap, and model nonlinearity prevent claiming a uniform
5 mm full-bar result. Current data do not permit a validated X–Y reconstruction.
