# What is comparable to the historical 85-90 ps?

| Estimator | Cancels position? | Requires per-event amplitude? | Comparable to historical weighted mean? | Useful for ToF? | Useful for x reconstruction? |
|---|---|---|---|---|---|
| t_avg = (t_L+t_R)/2 | Yes for symmetric two-ended propagation, but fragile when one side has very low PE | No | No, unless the reference used equal endpoint weights | Yes at center; poor at edges in EXEC_27 | Only with separate t_L-t_R or calibrated propagation |
| wmean_static | Partly; it weights by per-position timing precision | No, but needs a calibration map versus x | Yes, closest if historical analysis used calibrated endpoint weights | Yes, with x-dependent calibration | Not by itself; combine with asymmetry/timing difference |
| wmean_event using Npe | Partly; it follows per-event endpoint information | Yes | Plausibly comparable if historical weighted mean used amplitude weights | Yes, but can collapse to near-end timing at edges | Npe asymmetry itself carries x information |
| GLS_residual | Yes for resolution after removing per-position mean | No per event, but needs covariance calibration | Yes for intrinsic resolution diagnostics | Yes after calibration; not a raw timestamp | No, it removes the position mean by construction |
| near-end only | No | No | No | Good local timing diagnostic, not a two-ended ToF estimator | Strong x dependence; near/far identity gives side |
| far-end only | No | No | No | Diagnostic of light transport and far-side detectability | Strong x dependence; useful for transport validation |
