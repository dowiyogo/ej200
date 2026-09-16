# EXEC_46 Step 3 — pure transport g(d) and Cherenkov guiding

Date: 2026-09-16

## Checkpoint verdict

Step 3 is complete. No simulation was run and all production ROOT files were read-only.
The source-separated result identifies the fast Cherenkov edge quantitatively: the
first-Cherenkov linear slopes correspond to about 148 mm/ns, within 1.2% of the
parameter-free cone-edge prediction 146.903 mm/ns.
The all-photon means are much slower because they include recirculated paths. None of
the low-order global models has acceptable absolute chi-square, so these slopes are
diagnostic summaries rather than complete models of g(d).

This is the required checkpoint. Step 4 has not been started.

## Definitions and input checks

The production macros contain `/muon/angle 0` in all 21 cells. The source maps this
to momentum `(0,0,-1)`, perpendicular to the bar x axis. The effective RINDEX read
independently for all three materials is 1.58. Thus theta_C =
50.73475 deg, theta_crit = 39.26525 deg,
c/n = 189.742062 mm/ns, and the Cherenkov cone-edge axial
velocity is 146.902903 mm/ns.

The beta=1 identities are the requested a priori approximation. A 1 GeV kinetic-energy
muon has beta=0.995424; without any fitted parameter this moves theta_C to
50.51908 deg, the axial lower edge to
39.48092 deg, and its axial velocity to
146.449825 mm/ns.

For every photon:

```text
tprop = t_detection_ns - t_creation_ns
d_direct = |x_detection - x_creation| in three dimensions
rho_detour = path_length_mm / d_direct
v_apparent = d_direct / tprop
```

The primary g(d) calculation uses microscopic d_direct. Nominal distances are used
only to pair the left and right mirror realizations at 50, 200, 500, 700, 900, 1200
and 1350 mm. At 700 mm the two faces share the same x=0 simulated events, so they are
symmetry realizations but are not statistically independent in the strict sampling
sense; the all-photon summary does not store their covariance.

## B1 — withdrawn sparse-grid boundaries

`BOUNDARY_IDENTITY_NOT_ESTABLISHED` is retained and A2 is withdrawn. The correlation
changes sign between |x|=500 and 650 mm. Both reported linear crossings are artifacts
of interpolating across a regime change. The available grid samples only END distances
{50, 200, 500, 700, 900, 1200, 1350} mm; the transition lies wholly in the 50--200 mm
gap. No new simulation is proposed here.

## Fits to mirror-combined g(d)

The mandated linear fit is through the origin, t=d/v. The quadratic and cubic models
are also through the origin. Errors are event-cluster SEMs for all photons and event
SEMs for first-by-source photons. The shaded bands in `g_d.pdf` span the two mirror
means rather than pretending that a photon-IID error describes an event cluster.

### Linear summaries

| sample | material | source | v [mm/ns] | chi2/ndf | vs cone edge | vs c/n |
|---|---|---|---:|---:|---:|---:|
| all_photons | EJ-200 | scintillation | 134.786 +/- 0.009 | 3459432.2/6 = 576572.0* | -8.25% | -28.96% |
| all_photons | EJ-200 | Cherenkov | 110.633 +/- 0.046 | 127834.9/6 = 21305.8* | -24.69% | -41.69% |
| all_photons | EJ-204 | scintillation | 138.435 +/- 0.008 | 5431304.8/6 = 905217.5* | -5.76% | -27.04% |
| all_photons | EJ-204 | Cherenkov | 118.151 +/- 0.041 | 157486.7/6 = 26247.8* | -19.57% | -37.73% |
| all_photons | EJ-230 | scintillation | 138.520 +/- 0.008 | 6025468.8/6 = 1004244.8* | -5.71% | -27.00% |
| all_photons | EJ-230 | Cherenkov | 119.639 +/- 0.040 | 184914.5/6 = 30819.1* | -18.56% | -36.95% |
| first_by_source | EJ-200 | scintillation | 186.071 +/- 0.010 | 11765.8/6 = 1961.0* | +26.66% | -1.93% |
| first_by_source | EJ-200 | Cherenkov | 148.548 +/- 0.023 | 1636.5/6 = 272.7* | +1.12% | -21.71% |
| first_by_source | EJ-204 | scintillation | 186.060 +/- 0.010 | 7826.4/6 = 1304.4* | +26.65% | -1.94% |
| first_by_source | EJ-204 | Cherenkov | 148.057 +/- 0.026 | 4879.9/6 = 813.3* | +0.79% | -21.97% |
| first_by_source | EJ-230 | scintillation | 186.155 +/- 0.010 | 5571.4/6 = 928.6* | +26.72% | -1.89% |
| first_by_source | EJ-230 | Cherenkov | 147.863 +/- 0.028 | 6938.3/6 = 1156.4* | +0.65% | -22.07% |
| first_overall | EJ-200 | first overall | 185.996 +/- 0.010 | 48058.7/6 = 8009.8* | +26.61% | -1.97% |
| first_overall | EJ-204 | first overall | 185.993 +/- 0.010 | 31988.6/6 = 5331.4* | +26.61% | -1.98% |
| first_overall | EJ-230 | first overall | 186.113 +/- 0.010 | 21371.5/6 = 3561.9* | +26.69% | -1.91% |

`*` marks an inadequate absolute fit (all entries above).

### Curvature tests

| sample | material | source | Delta chi2 linear->quadratic | |b2|/err | quadratic chi2/ndf | cubic chi2/ndf |
|---|---|---|---:|---:|---:|---:|
| all_photons | EJ-200 | Cherenkov | 94028.8 | 306.6 | 33806.2/5 = 6761.2* | 15015.9/4 = 3754.0* |
| all_photons | EJ-200 | scintillation | 2162670.8 | 1470.6 | 1296761.4/5 = 259352.3* | 515331.0/4 = 128832.7* |
| all_photons | EJ-204 | Cherenkov | 104546.0 | 323.3 | 52940.7/5 = 10588.1* | 21457.1/4 = 5364.3* |
| all_photons | EJ-204 | scintillation | 3231982.0 | 1797.8 | 2199322.7/5 = 439864.5* | 791390.1/4 = 197847.5* |
| all_photons | EJ-230 | Cherenkov | 115850.0 | 340.4 | 69064.4/5 = 13812.9* | 26066.7/4 = 6516.7* |
| all_photons | EJ-230 | scintillation | 3603775.4 | 1898.4 | 2421693.4/5 = 484338.7* | 865787.4/4 = 216446.8* |
| first_by_source | EJ-200 | Cherenkov | 1536.9 | 39.2 | 99.6/5 = 19.9* | 82.5/4 = 20.6* |
| first_by_source | EJ-200 | scintillation | 8410.5 | 91.7 | 3355.3/5 = 671.1* | 561.1/4 = 140.3* |
| first_by_source | EJ-204 | Cherenkov | 4305.5 | 65.6 | 574.5/5 = 114.9* | 223.3/4 = 55.8* |
| first_by_source | EJ-204 | scintillation | 5181.9 | 72.0 | 2644.5/5 = 528.9* | 377.7/4 = 94.4* |
| first_by_source | EJ-230 | Cherenkov | 5948.6 | 77.1 | 989.6/5 = 197.9* | 329.5/4 = 82.4* |
| first_by_source | EJ-230 | scintillation | 3457.4 | 58.8 | 2114.0/5 = 422.8* | 158.5/4 = 39.6* |
| first_overall | EJ-200 | first overall | 12377.9 | 111.3 | 35680.8/5 = 7136.2* | 25679.9/4 = 6420.0* |
| first_overall | EJ-204 | first overall | 7480.7 | 86.5 | 24507.9/5 = 4901.6* | 18230.2/4 = 4557.5* |
| first_overall | EJ-230 | first overall | 4725.4 | 68.7 | 16646.1/5 = 3329.2* | 12025.0/4 = 3006.3* |

The residual-guided cubic improves chi-square again but remains rejected.
The opening into boundary-count families explains why a single low-order curve
does not describe the conditional transport mean.

## Comparison with the historical first-PE effective velocity

The historical estimator mixes emission and transport and reported the following
finite-difference ranges. The current first-overall row removes creation time event
by event before fitting, while the source-separated rows expose the mixture.

| material | historical local range [mm/ns] | historical linear | pure first-overall | first scintillation | first Cherenkov |
|---|---:|---:|---:|---:|---:|
| EJ-200 | 172.63--183.57 | 182.216 | 185.996 | 186.071 | 148.548 |
| EJ-204 | 175.04--182.56 | 181.652 | 185.993 | 186.060 | 148.057 |
| EJ-230 | 177.26--182.76 | 181.611 | 186.113 | 186.155 | 147.863 |

The source-specific Cherenkov value is not expected to equal the historical
mixed estimator. Its agreement target is the cone-edge speed, which it meets to
better than 1.2% in all three materials. The scintillation-only first photon is
near 186 mm/ns and therefore close to c/n, as expected for the earliest isotropic
photons selected from a large population.

## B2 — mirror consistency

The full table is `mirror_consistency.csv`: 84 source-separated rows plus 21
first-overall diagnostic rows. Compact maxima are:

| sample | material | source | max |Delta mean| [ps] | max |pull| |
|---|---|---|---:|---:|
| all_photons | EJ-200 | Cherenkov | 8.012 | 1.33 |
| all_photons | EJ-200 | scintillation | 2.835 | 1.44 |
| all_photons | EJ-204 | Cherenkov | 39.654 | 2.36 |
| all_photons | EJ-204 | scintillation | 1.254 | 0.93 |
| all_photons | EJ-230 | Cherenkov | 14.693 | 1.35 |
| all_photons | EJ-230 | scintillation | 1.882 | 1.20 |
| first_by_source | EJ-200 | Cherenkov | 5.443 | 1.96 |
| first_by_source | EJ-200 | scintillation | 2.216 | 1.47 |
| first_by_source | EJ-204 | Cherenkov | 15.767 | 1.37 |
| first_by_source | EJ-204 | scintillation | 2.274 | 2.03 |
| first_by_source | EJ-230 | Cherenkov | 5.477 | 1.39 |
| first_by_source | EJ-230 | scintillation | 2.454 | 1.59 |
| first_overall | EJ-200 | first overall | 1.597 | 1.90 |
| first_overall | EJ-204 | first overall | 2.132 | 1.88 |
| first_overall | EJ-230 | first overall | 2.629 | 1.60 |

The largest pull is 2.36 and the largest absolute mirror difference is
39.654 ps. The two realizations are consistent at the declared precision.

## B3 — direct Cherenkov edge test

The primary-cone geometry applies because the gun is perpendicular to x. The
identity with theta_crit is exact only in the beta=1 limit; the configured finite
beta predicts the lower edge at 39.481 deg.
However, `source_type==2` records the creator process and not parent track identity.
It therefore includes Cherenkov photons made by secondary charged particles. The
`primary_like_proxy` requires creation x within 0.001 mm of gun x and creation y
within 0.001 mm of zero; it is a geometric proxy, not a recovered parent label.

| material | selection | N | q0.1% [deg] | median [deg] | fraction below theta_crit | modal 0.02-deg bin |
|---|---|---:|---:|---:|---:|---:|
| EJ-200 | all_source_type_2 | 20000 | 1.520 | 39.640 | 11.6600% | 39.50--39.52 |
| EJ-200 | primary_like_proxy | 11116 | 39.412 | 39.694 | 0.0270% | 39.50--39.52 |
| EJ-204 | all_source_type_2 | 20000 | 1.504 | 39.642 | 11.7850% | 39.50--39.52 |
| EJ-204 | primary_like_proxy | 11059 | 39.425 | 39.697 | 0.0000% | 39.50--39.52 |
| EJ-230 | all_source_type_2 | 20000 | 1.650 | 39.644 | 11.5400% | 39.50--39.52 |
| EJ-230 | primary_like_proxy | 11025 | 39.425 | 39.700 | 0.0181% | 39.50--39.52 |

The finite-beta edge at 39.481 deg lies inside a partly filled histogram bin.
The modal 39.50--39.52 deg bin immediately to its right therefore does not imply
a 0.019-deg discrepancy; the agreement is limited by the 0.02-deg binning and
exhibits the expected caustic pile-up. Isolated 41--43 deg teeth contain about
one photon per bin and are Poisson noise, not angular discretization. The upper
edge of the selected angular window is temporal selection, not a second geometric
cone boundary. In the beta=1 limit the cone-edge and END critical angles coincide
exactly through arccos(sin(theta_C)) = arcsin(1/n). The undifferentiated source-type-2
population fails the strict lower-bound
test: about 11.5--11.8% lies below theta_crit. The primary-like proxy satisfies
the bound for 99.97% or more of photons and displays the expected edge pile-up.
The few remaining proxy violations show that a creation-position cut cannot prove
parentage. The direct test therefore confirms the primary-cone mechanism while
also measuring the secondary contamination that prevents applying the identity to
all source-type-2 photons.

## B4 — source-separated transport

The predicted beta=1 edge velocity is 146.903
mm/ns, the finite-beta value is 146.450
mm/ns, and the group velocity is 189.742 mm/ns. The fitted first-
Cherenkov velocities are listed above: they are 0.65--1.12% above the edge
beta=1 prediction, 0.96--1.43% above the finite-beta prediction, and 21.7--22.1%
below c/n. This closes the axial guiding mechanism at
the precision allowed by a global linear summary. Its quadratic residual is still
statistically significant, so the cone edge is not the whole detected population.

All detected Cherenkov photons yield much lower linear summaries (111--120 mm/ns),
while all scintillation photons give 135--139 mm/ns. These are path-population
means, not material group velocities. `g_by_boundary.pdf` separates them into the
mandated 0, 1--2, 3--5, 6--10 and >10 encounter families; `g_by_exit_angle.pdf`
shows the corresponding final-angle stratification.

## Microscopic binning and TH2 coverage

The 10 mm microscopic bins contain 313,885,061 photons. The 0--60 ns TH2
excludes 1,998 upper-overflow photons (0.000637%) only
from the raster; all moments and fits retain them. `tprop_variance.csv` contains
the binned mean and conditional variance by material, face and source.

## Reproducibility

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/build_step3_transport.py --processes 4
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step3.py
```

Primary artifacts:

- `all_photon_cell_face.csv`: event-cluster means for 84 cell/face/source groups.
- `first_by_source.root`: first photon per event, END face and source.
- `fit_summary.csv`, `mirror_consistency.csv`, and `cherenkov_edge_summary.csv`.
- Each figure has `.pdf`, `.root`, `.csv`, and `.meta.json` sidecars.

The build wall time was 226.2 s with four read processes. No push, merge, deck
edit, or Geant4 run was performed.
