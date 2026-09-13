# EXEC40: registered before the sole transport run

Base 420addf; rollback pre-exec40-20260913; all three observations are
production fields: deposited energy, first physical surface encounters, and
incident/detected counts are needed to interpret physical outputs. Existing
internal diagnostic censuses remain optional. No changes to RNG calls, track
status, materials, surfaces, geometry, or sensitive-detector decisions.

One cell only: EndTop, 70 TOP, OPSC-101/EJ-204, mu- 1 GeV vertical,
x=0 mm, 2000 events, seeds 26092601 8349041, 4 workers, eventModulo 1,
diagnostics OFF, jitter 0 ns. No smoke test may run additional beam events.

G-I.1: exactly 397.1525 pe/end and every event's L/R/T counts equal archived
EXEC29 D1. Any nonzero difference aborts physical validation, without retuning.
G-I.2: 2000 readable event rows, readable first-encounter and per-SiPM tables,
unique first photon keys, valid signed incidence cosines, independently matched
incident/detected observations. Any deficiency is reported, not inferred away.
G-I.3: report bytes and wall versus D1; also quote archived EXEC33 S2 with
matching 4-worker diagnostics-OFF configuration. A size ratio >2 requires a
proposal, not an automatic move behind the diagnostic flag.

Statistical convention: resample generated events jointly, 500 bootstrap
replicates, NumPy seed 40091301. Standard errors use ddof=1. Report point
estimates and 3-SE tolerance, fixed before looking at new results.

V1 primary: sum(produced_scint)/(10400*sum(edep_total_MeV)); require unity
within 3 bootstrap SE. Total includes every depositing particle in BarLV.
Record nonionizing and total-minus-nonionizing energy separately. Also report
the nonoptical-deposit denominator as a declared secondary diagnostic, never
as a substitute for the primary outcome. Datasheet yields: EJ200 10000,
EJ204 10400, EJ230 9700 photons/MeV.

V2 primary: scintillation photons whose first physical encounter exits BarLV;
outward refraction/transmission divided by all such encounters. Require overlap
of its 3-SE interval with [0.36,0.40]. Also report the bar-to-air subset and
all-source subset, keeping sensor absorption distinct from escape to air.
Cosine goodness of fit: 20 equal bins on [0,1], null p(mu)=2*mu; compare
normalized histogram with bin probabilities b^2-a^2, using covariance of
event-bootstrap histograms, pseudoinverse chi-square and its numerical rank;
require p>=0.01. Report conventional photon Pearson chi-square separately.
Check signed cosines in [-1e-9,1+1e-9], normal norm, and geometric versus
navigator normal agreement; never abs() or flip a normal to force this check.
Exclude nonphysical StepTooSmall/NotAtBoundary/Undefined bookkeeping steps.
The literal flux-isotropy hypothesis is retained, but first hits from a localized
source in a finite bar need not have the equilibrium isotropic-flux law.
Failure alone is not evidence of a faulty isotropic generator.

V5 primary: detected/independently observed incident photons, compare against
0.619397 (EJ204), 3 paired-bootstrap SE; report left/right/TOP/all separately.
EJ200 reference 0.606655; EJ230 0.610613. Count a unique (event,track,SiPM)
when the post-volume is a sensor, pre differs, and status is Detection,
Absorption, FresnelRefraction, Transmission, SameMaterial or coated transmission.
Terminal surface acceptance counts even though a metal-like sensitive boundary
absorbs rather than transporting the photon inside. Reflections do not count.
Detection callbacks are observed separately and matched at event end, so callback
ordering cannot create incidents. Retain unmatched counts as a closure check.
Record actual surface-efficiency expectation on the incident spectrum separately
from the emitted-spectrum reference; this diagnostic cannot rescue primary V5.

Storage: event and channel tables scale as O(events*channels), with zero-hit
channels retained. One compact first-encounter row per photon scales with
generated photons (potentially comparable to existing ~3 GB/10000-event hits).
No per-bounce history is added. Measure actual growth before proposing changes.

Full grid acceptance remains false irrespective of this one-cell outcome.
