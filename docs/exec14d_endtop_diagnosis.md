# EXEC_14D EndTop/End-only diagnosis

## Verdict: PHYSICAL

The apples-to-apples comparison uses the same End channels, SUM4 A/B assignment,
leading-edge implementation, absolute 4 PE threshold, and estimator
`sigma_group = std(delta_t_AB)/sqrt(2)` for both configurations. No fit window is
applied to either campaign. End-only is bootstrapped to 2000 events
with 200 replicas (seed 14014).

| x [mm] | side | ratio at 2000 events | Npe EndTop / End-only | t99 EndTop / End-only [ns] | f(t>10) EndTop / End-only |
|---:|:---:|---:|---:|---:|---:|
| -400 | L | 1.223 +/- 0.020 | 581.1 / 939.4 | 9.426 / 9.961 | 0.00697 / 0.00977 |
| -400 | R | 1.963 +/- 0.033 | 53.0 / 160.4 | 12.303 / 13.480 | 0.04341 / 0.06997 |
| +0 | L | 1.393 +/- 0.024 | 151.0 / 350.0 | 10.900 / 11.840 | 0.01756 / 0.02895 |
| +0 | R | 1.442 +/- 0.024 | 151.1 / 350.5 | 10.942 / 11.843 | 0.01822 / 0.02892 |
| +400 | L | 1.963 +/- 0.032 | 52.8 / 160.5 | 12.327 / 13.509 | 0.04256 / 0.07033 |
| +400 | R | 1.234 +/- 0.019 | 577.8 / 938.3 | 9.435 / 9.953 | 0.00703 / 0.00973 |

The lower 95% bootstrap quantile is above one in all six comparisons:
`1.186` minimum. End-only retains
`1.62` to
`3.04` times more End
photoelectrons. The corresponding pure-counting prediction,
`sqrt(Npe_End-only/Npe_EndTop)`, spans
`1.271--1.743`.
EndTop nevertheless shortens the side-specific t99 by
`0.518--1.182` ns.

## Geometry validation

- All three log pairs differ only on +Y: EndTop is instrumented with 70 Top SiPMs; End-only is wrapped. End faces, other wraps, PDE, reflectivity, and 16 End SiPMs match.

## Interpretation

The tail-drain mechanism persists, but its benefit is smaller than the timing
penalty caused by the loss of End photoelectrons. For this fast EJ-230
scintillator the EndTop configuration therefore widens intrinsic End timing.
The sign is physical within the tested configurations and is not produced by
the 10000-versus-2000 event-count difference.
