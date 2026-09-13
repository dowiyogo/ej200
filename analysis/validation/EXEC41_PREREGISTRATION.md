# EXEC41 preregistration — corrected first-encounter hypothesis

Published before any EXEC41 numerical calculation. It supersedes the physical
interpretation of V2 in `EXEC40_PREREGISTRATION.md` (SHA256
`cd13d28db1811e5c74d9b17d60d1be7aa0971f02ff8768d48528126602fb5531`),
which remains unchanged. It does not erase the old result or alter any data.
All calculations use the already generated EJ-204 EndTop70 x=0, N=2000 cell,
seeds 26092601 and 8349041. No simulation is authorized for EXEC41.

## V2 populations and hypotheses

The analysis retains the EXEC40 primary population: first physical BarLV
surface encounter of scintillation photons exiting BarLV. Escape is portable
outcome 2 (transmission/refraction) into air or world; sensor acceptance is not
escape. Event bootstrap uses 500 paired generated-event resamples with NumPy
seed 41091301 and reports one standard error.

**H1 — solid-angle lower bound.** For isotropic emission, the emission-angle
cosine is uniform. With n=1.58, cos(theta_c)=0.773828 and the fraction in the
escape cone of a specified face is 1-cos(theta_c)=0.226172. Treat 0.226172 as a
lower bound for the aggregate first-encounter escape fraction. H1 passes if the
measured aggregate fraction minus 3 event-bootstrap SE is at least 0.226172.
Decompose both encounter share and conditional escape fraction by first face.
The two nearest ±Z faces are tested descriptively for compatibility with
0.226172: their 3-SE interval must contain that value. Attribute aggregate
excess arithmetically to escaped photons whose first encounter is on other face
classes. Do not fit geometry or change populations to obtain agreement.

The available production record names AirGapZPlus, AirGapZMinus and
AirGapYMinus directly; EndSiPMLeft, EndSiPMRight and TopSiPM identify -X, +X
and +Y sensor patches. A BarPV→WorldPV pair does not encode which uncovered
part of +Y/-X/+X was crossed, because position and normal components were not
stored. Report that irreducible class as `OPEN_+Y_OR_+/-X_UNRESOLVED`; do not
assign its rows statistically to faces. This is a data limitation, not a reason
to invent a decomposition.

**H2 — angular form.** Compare the signed incidence-cosine histogram with a
uniform density using the same 20 bins and event-bootstrap covariance method as
EXEC40. Report the deviation and the low-mu deficit/excess descriptively. H2 has
no PASS/FAIL because an analytic first-hit law requires the complete finite-bar
geometry. Also retain the old p(mu)=2mu statistic and old 0.36–0.40 FAIL beside
the corrected interpretation.

Geometry is read from code constants. At x=0, y=0 the source track spans
z=-5..+5 mm: distance to the nearest ±Z face is 0..5 mm; distance to -Y and
the nominal +Y plane is 30 mm; distances to ±X are 700 mm. Embedded SiPM
patches occupy parts of +Y and ±X and do not change those nominal face planes.

## Detection-gap diagnosis

Characterize only distinctions supported by existing records: orphan counts
by event, face and sensor are the difference between independent SD detections
and matched boundary detections. Since `sipm_hits` lacks track ID/boundary pair
and `sipm_event_counts` stores aggregates, individual orphan photon energies
and volume pairs cannot be identified after the run. Report this limitation.
Inspect the actual SD and Geant4 callback paths and enumerate any boundary
status excluded by the registered incident criterion. A proposed correction
must remain unapplied. Retain the matched result 0.61240 ± 0.00011 versus the
incident-spectrum expectation 0.61246 ± 0.00003, and quantify wavelength/PDE
differences available from existing aggregate and hit data.

## Storage reduction acceptance

Move raw toolkit `boundary_status`, `energy_eV`, `normal_norm`, `normal_valid`
and `normal_orientation_valid` behind `EJ200_ENABLE_DIAGNOSTICS`. In production,
replace `pre_volume` and `post_volume` strings with this integer dictionary:
0 unknown, 1 BarPV, 2 AirGapYMinusPV, 3 AirGapZPlusPV, 4 AirGapZMinusPV,
5 EndSiPMLeft_PV, 6 EndSiPMRight_PV, 7 TopSiPMPV, 8 WorldPV. Preserve event and
track identity, source, copy numbers, exiting direction, portable outcome and
signed incidence cosine.

Store the cosine as IEEE-754 binary32 (`G4float`, ROOT float branch): about
7 significant decimal digits and at most half an ulp rounding. On the existing
ROOT, construct a filtered representation using the dictionary and float32,
then repeat V2. Acceptance requires absolute shifts in every reported escape
fraction and each of 20 normalized cosine-bin probabilities to be <= their
original one-SE uncertainties. V1 and V5 use unchanged trees and must be
bit-identical when the filtered encounter tree is substituted. If a change
fails, do not apply that concrete reduction.

Measure the filtered encounter-tree bytes directly. Estimate a full ROOT as
filtered encounter bytes plus the measured original compressed bytes of the
unchanged sipm_hits, event_observables and sipm_event_counts trees; label this
as an estimate because different writers/compression prevent an exact full-file
comparison without regeneration. No new transport run is permitted.

## V1 investigation

Use existing per-event counts and energy fields. Check exact bookkeeping among
total, optical, nonionizing and ionizing deposits. Compare generated
scintillation photons with the Geant4 expectation based on the energy deposit
on the parent steps, including the configured `RESOLUTIONSCALE=1` fluctuation
model. Inspect the live physics configuration for Birks/quenching and
particle-dependent yield. Do not change material parameters. Distinguish an
energy-accounting denominator error from a changed effective scintillation
yield; do not select a denominator after observing agreement.

The golden reference remains `ready_for_acceptance=false`; a 21-cell grid still
requires separate authorization.
