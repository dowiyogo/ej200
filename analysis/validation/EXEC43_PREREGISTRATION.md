# EXEC43 preregistration — mechanism of material-dependent first escape

Published before numerical evaluation of EXEC43. Only the 21 existing ROOT
files selected and verified in EXEC42 are used. No simulation, parameter fit,
deck edit, or golden-contract change is authorized.

## Inputs and uncertainty convention

The input inventory is
`/home/rrios/exec42_20260913/analysis/input_inventory.csv` (SHA256
`012bb6e02f3cacdfd05c418ada06c031942032882aa4f9a886bdfadc36790ea7`).
All cells contain 10,000 generated events and use seeds 26092601 and 8349041.
Uncertainties use 500 generated-event bootstrap resamples, NumPy seed 43091401,
`ddof=1`, and identical event-index weights in all cells. Comparisons at fixed
position use the bootstrap distribution of the paired difference.

The first-encounter population and portable escape outcome are unchanged from
EXEC41/42: source 1 scintillation photons at their first physical BarLV
encounter while exiting the bar; outcome 2 into AirGapYMinusPV,
AirGapZPlusPV, AirGapZMinusPV, or WorldPV is escape. Sensor acceptance is not
escape.

## Test 1 — effective refractive index (M2)

Read the exact `RINDEX` and `SCINTILLATIONCOMPONENT1` tables used by OPSC-100,
OPSC-101, and OPSC-106. Integrate each emission spectrum by piecewise-linear
trapezoids in wavelength and compute
`n_eff = integral n(lambda) S(lambda) dlambda / integral S(lambda) dlambda`,
interpolating the corresponding RINDEX table linearly. Report the emission
mean wavelength under the same convention and the geometric cone fraction
`1 - sqrt(1 - 1/n_eff^2)`.

M2 is quantitatively excluded if `max(n_eff)-min(n_eff) < 1e-4`. Otherwise,
compare the predicted ordering and magnitude with the measured aggregate V2
escape fractions without adjusting a dispersion model.

## Test 2 — path to first encounter (M1 directional prediction)

The production schema stores no creation point, encounter point, time, or
track length. It stores face and signed incidence cosine. Consequently, exact
absolute paths are evaluable only where the source-to-plane perpendicular
distance is fixed; the common dimensionless measured path factor
`s = sec(theta) = 1/cos_incidence` is evaluated for every face and outcome.
Report counts and histogram-derived median, 10th and 90th percentiles of `s`,
separately for each cell, face, and escaped/non-escaped outcome. A face-level
comparison is not evaluable when either outcome is absent.

The central directional prediction passes descriptively for an evaluable face
when median `s(non-escape) > s(escape)`. Also report exact geometric path
quantiles `d*s` for the fixed-distance faces: -Y and +Y use d=30 mm; the -X
and +X sensor patches use d=700+x and d=700-x mm. For +/-Z, photon production
spans z=-5..+5 mm, so d is unrecorded in [0,10] mm. WorldPV does not identify
the +Y/-X/+X face. Do not label an assumed midpoint path as measured.

## Test 3 — photons absent before a first encounter (M1 survival prediction)

For every cell calculate
`f_absent = 1 - N_first_encounter/N_produced_scint`, using event sums. A
material pair at fixed x supports the predicted ordering when the shorter
attenuation length has the larger absent fraction and the paired difference is
positive by more than three bootstrap SE. M1 survival ordering passes only if
EJ-200 (3800 mm) < EJ-204 (1600 mm) < EJ-230 (1200 mm) in all seven positions
and all 14 adjacent pair comparisons exceed three paired SE.

## Test 4 — quantitative attenuation closure

Use the measured EJ-200 first-face and cosine distribution as the longest-
attenuation baseline. For a fixed perpendicular distance d, reweight every
surviving encounter from lambda_0=3800 mm to lambda_1 by
`exp[-d*sec(theta)*(1/lambda_1 - 1/lambda_0)]`. For +/-Z, integrate this weight
over a uniform unrecorded production distance d in [0,10] mm. This is an
explicit geometric source model necessitated by the absent creation point. For
WorldPV, calculate an ambiguity envelope by assigning all rows in turn to
d=30 mm, d=700-|x| mm, and d=700+|x| mm; do not choose among them after seeing
the result. All other resolved distances are fixed as in Test 2.

Predict the aggregate escape fraction for EJ-204 and EJ-230 at each x from
the reweighted EJ-200 population. Propagate event-bootstrap uncertainty using
the common weights. Closure at a target cell is compatible when its observed
escape lies inside the WorldPV prediction envelope expanded by three paired
bootstrap SE at each boundary. M1 quantitative closure passes only if all 14
target cells are compatible. Report residuals and envelope widths; no fitted
parameter or empirical normalization is allowed.

Because exact path length is absent, this closure is explicitly
`GEOMETRY-CONSTRAINED ESTIMATE`, not an exact reconstruction from persisted
tracks. Even if it passes, the report must distinguish evidence supporting M1
from a direct measurement of the full path distribution.

## Test 5 — position dependence

For each of the three material pairs form the seven paired escape differences.
Report their range, mean, and fractional peak-to-peak variation. Test a
constant difference with the covariance of the seven paired-bootstrap
profiles, a pseudoinverse chi-square, its numerical rank, and p-value. Call
position dependence statistically resolved only at p<0.01. This test is
descriptive and does not alter Tests 1–4.

## Mechanism decision and proposed H3'

Call M1 confirmed within the persisted-data limits only if M2 is excluded,
the Test-2 directional comparison passes for every evaluable nonsensor face,
the Test-3 survival ordering passes, and all 14 Test-4 closures are compatible.
Otherwise report M1 as supported, disfavored, or not fully evaluable according
to the failed component; do not promote a mechanism by convenience.

If M1 is confirmed, propose but do not apply H3': first-encounter escape may
depend on material only through pre-encounter bulk-absorption filtering, and
every material-pair difference must agree with a preregistered transport
reweighting from attenuation lengths and a path distribution within a
declared tolerance. A production H3' must require persisted exact path length
or creation/encounter coordinates; the geometry-constrained proxy in this
task is insufficient to become the golden acceptance definition without
René's decision.
