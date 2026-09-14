# EXEC42 Phase B preregistration — H3 material invariance

Published before calculation of any 21-cell EXEC42 contract result. This file
adds H3 only. It does not alter the definitions in
`EXEC40_PREREGISTRATION.md` (SHA256
`cd13d28db1811e5c74d9b17d60d1be7aa0971f02ff8768d48528126602fb5531`)
or `EXEC41_PREREGISTRATION.md` (SHA256
`933433a2882c4753301c675913df1f37db614e4afc8316ae87df6c122021d418`).
Only the completed existing ROOT files under
`/home/rrios/exec42_20260913/grid/` are analyzed; no simulation or relaunch is
authorized.

## Fixed population and uncertainty convention

V2 uses first physical BarLV encounters of scintillation photons exiting the
bar. Portable outcome 2 into AirGapYMinusPV, AirGapZPlusPV,
AirGapZMinusPV or WorldPV is escape. Sensor acceptance is not escape. Volume
IDs are the EXEC41 production dictionary. All cell metrics resample the 10,000
generated events 500 times with NumPy seed 42091401 and report `ddof=1`
standard errors. The identical bootstrap event-index weights are used in every
cell because all cells share seeds and event IDs.

For comparisons between materials at one position, compute the V2 escape
ratio in each material from the same bootstrap event resample and use the
standard deviation of the paired ratio difference. Do not combine independent
cell SEs in quadrature. A comparison is not evaluable if either cell lacks all
10,000 event rows, has a zero denominator, or lacks the required production
branches.

## H3 — invariance of escape with material

At each of x = 0, ±200, ±500 and ±650 mm, compare all three material pairs:
EJ-200 versus EJ-204, EJ-200 versus EJ-230, and EJ-204 versus EJ-230. A pair is
compatible when the absolute escape-fraction difference is no more than three
paired-bootstrap SE. H3 passes only if all 21 pairwise comparisons are
evaluable and compatible.

The position profile is descriptive. Dependence on x, especially at ±650 mm,
is expected from finite geometry and receives no PASS/FAIL without a complete
geometric first-hit calculation.

## Existing contract rules applied without modification

V1 primary is `sum(produced_scint)/(nominal yield * sum(total deposited
energy))`, with yields 10,000, 10,400 and 9,700 photons/MeV for EJ-200, EJ-204
and EJ-230. It passes when unity is inside its three-SE interval. The
non-optical-deposit denominator is reported separately and cannot replace the
primary result.

V2-H1 passes when aggregate escape minus three SE is at least 0.226172. Report
the seven frozen first-face classes and their conditional escape fractions.
The 20-bin uniform and old `p(mu)=2mu` angular comparisons remain descriptive
under EXEC41 and do not enter acceptance.

V5 uses matched boundary detections divided by independently observed
incidents. Compare it with `sum(expected_surface_pde)/sum(incident)` on the
same incident population. It passes when their difference is compatible with
zero within three event-bootstrap SE of the paired ratio difference. The
legacy independent-SD detection fraction and its orphan component are reported
separately and cannot replace or rescue matched V5.

`ready_for_acceptance` may become true only if V1 primary, V2-H1 and matched
V5 pass in all 21 cells, H3 passes all 21 comparisons, and no required result
is not evaluable. The shared simulation seeds make this a paired design; points
within a material are correlated and must not be treated as independent in any
later fit or tolerance adjustment.
