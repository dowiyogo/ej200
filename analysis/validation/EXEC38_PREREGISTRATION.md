# EXEC38 preregistration — 2026-09-13

Registered before new numerical analysis. Earlier reports and schemas have been
read; their published values are not blinded and V6 is explicitly retrospective.
No transport simulation is authorized. No physics or timing parameters change.

## Inputs and estimands

Use the 21 successful EXEC34R cells (EndTop, 70 TOP sensors, three materials,
x = -650,-500,-200,0,200,500,650 mm, 10000 generated events, seeds 26092601 and
8349041, four workers, eventModulo=1, source 420addf), plus archived D0/D3
(END-only, EJ204, x=0, 2000 events, same seeds, sources 2ddf205/391fa4e).
Resolve from manifests. Hash and validate existing EXEC36 END hit caches and
EXEC36/37 same-end timestamp caches against their recorded provenance. Missing
streams must produce NOT_EVALUABLE, never a synthetic measurement or PASS.

All acceptance definitions use physical source photons, surface encounters,
destinations and timestamps. Engine-specific diagnostic enum values are not
part of the acceptance API. Source configurations and simulation commands are
preserved as provenance, not instructions to run a simulation.

Use 500 event bootstrap replicas, NumPy PCG64 seed 38091301. Resample all
generated event indices; retain the same resample across positions, ends and
materials sharing an event population, to preserve observed paired covariance.
N=2000 and N=10000 populations are separate. Report bootstrap standard errors;
all quoted errors are statistical and conditional on the frozen configuration.
At least 95% finite replicas and 200 accepted events are required. No electronics,
BLUE, test-beam numerical comparison, trimming or fitted parameter retuning.

## Predictions and tolerances

**V1:** sum(N_scintillation_produced)/(yield * sum(deposited_energy_MeV)) = 1.
Use actual recorded energy deposition, not inferred from photons. Eljen gives
10000/10400/9700 photons/MeV for EJ200/204/230. Accept |ratio-1| <= 3 standard
errors (event-paired propagation); this tests implementation of the specified
yield, not unreported manufacturer systematic precision. If dE is absent, the
production count alone is descriptive and the ratio is NOT_EVALUABLE.

**V2:** for flux-weighted isotropic incidence at a bar-air surface,
theta_c=asin(1/1.58), P(cone)=1/1.58^2; prescribed first-encounter escape band
[0.36,0.40] with Fresnel losses. Accept a measured fraction within this band
expanded by 3 event-bootstrap SE; no arbitrary percent padding. Restrict to
scintillation photons and each photon's first encounter; only bar-air surfaces
have this critical angle. Record separately coupled sensor surfaces. The
flux-weighted angular assumption is not automatically implied by a localized
isotropic source in a finite bar; applicability must be established as well as
the encounter ordering. Aggregate encounter tallies cannot pass this test.

**V3:** effective axial attenuation length < nominal datasheet light attenuation
length (3800/1600/1200 mm). The datasheet calls this light attenuation length,
not an independently certified microscopic absorption length; the source uses
these nominal values as ABSLENGTH. For each end use distance to that end,
d_left=700+x, d_right=700-x (mm), not |x| which merges near and far observations.
Fit Npe(d)=A exp(-d/lambda) with full covariance of event means, GLS in count
space. Fit A1 exp(-d/lambda1)+A2 exp(-d/lambda2) as a four-parameter alternative.
Positive amplitudes, lengths bounded [1,1e6] mm, amplitudes [1e-6,1e7]; fixed
multistart short/long length pairs (100,1000),(300,3000),(500,10000) mm and equal
amplitudes initially. Seven points; chi2 dof=5 and 3. Descriptive GOF requires
p>=0.01. Use AICc=chi2+2k+2k(k+1)/(7-k-1); require improvement >=6 and adequate
GOF to prefer two components. AICc is a model-selection heuristic, not a
calibrated likelihood-ratio test at the coincident-length boundary. Length SE
from bootstrap refits; if the model fails GOF, retain descriptive parameters
but do not claim a validated single attenuation scale. No fit-range changes.
For the single-length bound: PASS if lambda+3SE<nominal, FAIL if
lambda-3SE>nominal, otherwise INDETERMINATE. Compare the three nominal/lambda
ratios per end with paired bootstrap covariance: common-ratio GLS chi2 (2 dof),
p>=0.01 for consistency. This is a stated approximate zigzag hypothesis, not a
universal physical identity: different spectra, surface losses and scattering
can select different trapped populations. Report model, bound and ratio tests
separately; failure of either hypothesis is retained.

**V4:** c/n=299.792/1.58 mm/ns, inverse slope v_eff expected below c/n. Use the
photon-weighted mean detected arrival time relative to gun=0, sum(t)/Nphotons,
per end and position. Full bootstrap covariance, linear GLS t=a+b*d and
v_eff=1/b. Quote signed slope versus x as well. PASS if b>0 and
v_eff+3SE<c/n; FAIL if b<=0 or v_eff-3SE>c/n; otherwise INDETERMINATE.
Require linear-model p>=0.01 to interpret one speed across the profile; retain
the bound result separately if linearity fails. A slope of selected detected
photon means is not an individual photon's causal speed and need not obey a
strict causal bound without a stable path population. Diagnose this limitation
if the requested effective-speed hypothesis fails; do not change the estimator.

**V5:** Ndetected/Nincident for a specified end equals the emission-spectrum
weighted PDE, integral S(E)*PDE(hc/E)dE / integral S(E)dE, using the actual
configured emission table and the fixed Broadcom PDE file. Also disclose the
wavelength-density convention where applicable. Accept difference <=3 combined
statistical SE. No manufacturer systematic error is supplied. Incident photons
must be independent recorded arrivals at that end, including unsuccessful
detections; detected photons or inverse-PDE weights cannot reconstruct a
measured denominator. A total all-sensor count cannot replace an end count.
If transport reshapes the incident spectrum, the emission-average equality is
only a hypothesis; the exact detector closure uses the incident spectrum and
sum of its individual PDE probabilities. Do not silently conflate them.

**V6:** retrospective source audit, without recalculating the archived physical
measurements. Photon conservation expects exactly one terminal physical
destination per source photon; future tolerance is zero missing/duplicate
photons (integer identity). Reflector return probability expects 0.98; future
tolerance 3*sqrt(.98*.02/N_incident), with event-cluster SE instead if supplied.
No TIR for n2>n1: exact zero in an independently defined angular/physical
population. Mirror symmetry expects equality for mirrored configurations:
future |difference| <= 3 paired event-bootstrap SE, with no invented fixed
percentage tolerance. Worker-count reproducibility expects identical per-event
destination counts for the specified common-random-number CPU configuration,
exact tolerance zero. Different RNGs/engines use physical statistical criteria,
not CPU event identity. Existing results retain their original scope and gates;
missing or misattributed archived claims are UNVERIFIED, not retrospective PASS.

**V7:** strictly positive same-end Pearson rho(T1,T2). Primary V2 SUM4 groups
left {0,1,2,3}/{4,5,6,7}, right {8,9,10,11}/{12,13,14,15}; unchanged 0.5 ns rise,
5 ns fall, 4 PE crossing. Also report existing EndTop V1 subset as separately
labelled secondary evidence, never select the better outcome. D0/D3 V2 only.
PASS if rho-3*bootstrap_SE>0; FAIL if rho+3SE<0; otherwise INDETERMINATE
(compatible with zero). Standard deviations use ddof=1, paired accepted events,
no Gaussian fit or core cut: Pearson's variance identity is an ordinary second-
moment identity, not an identity of the EXEC35 quantile core estimator.
Report s1,s2,sDelta,sDelta/sqrt(2), Q=s1/(sDelta/sqrt(2)), F=1/sqrt(1-rho),
and Q-F with paired bootstrap SE. Equal-variance approximation compatible if
|Q-F|<=3SE. Exact unequal-variance prediction is
sqrt(2*s1^2/(s1^2+s2^2-2*rho*s1*s2)); check identity within relative 1e-10
(floating-point check only, not independent physics evidence). Report direct
gun-relative sigmas and conditional efficiency. No claim of experimental access
to a jitter-free gun time. No numerical experimental reference is used.

## Primary sources

- Eljen EJ200/EJ204: https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212
- Eljen EJ230, August 2023: https://eljentechnology.com/images/products/data_sheets/EJ-228_EJ-230.pdf
- Snell's law and flux measure sin(theta)cos(theta)dtheta (V2); exponential
  path survival exp(-L/lambda) (V3); t=L/(c/n) for constant n (V4); Bernoulli
  detection expectation (V5); conservation/reciprocity (V6); covariance identity
  Var(T1-T2)=Var(T1)+Var(T2)-2Cov(T1,T2) (V7).
- Historical commits 55766876ffa66e44e8a42461029d6276ff14d164 and
  e4110ce5f83bbf8bbb7e3e1731805e3631054c09 are diagnostic history, not predictions.

## Contract status

The golden file will contain one entry per V1–V7 with explicit definitions,
source, tolerance, measured values/errors or null, exact configurations and
source hashes. It is not certified ready for acceptance if required observables
are absent or hypotheses fail. Preserve these failures; never promote an
unavailable test to PASS or use a current simulated value as its prediction.
