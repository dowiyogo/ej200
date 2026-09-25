# Current optical-performance metrics and Beamer v9

Date: 2026-09-14  
Host: `t0minidaq`  
Repository working tree: `/home/rrios/ej200`  
Presentation: `/home/rrios/ej200/presentations/v9/talk_v9.pdf`

## Scope and decision on simulation

The presentation answers the detector-performance question using the corrected
optical transport. It is not a repository audit or another transport-validation
deck. No Geant4 simulation was launched.

The input audit found 21 complete, readable ROOT cells under
`/home/rrios/exec42_20260913/grid/cells/`: three materials and seven positions,
with 10,000 generated events per cell. Their per-cell metadata identify:

- simulation source `b35ee84acadef12c93506e0720580c91f901fcbf`;
- Geant4 11.4.0;
- EndTop geometry, 16 END and 70 TOP SiPMs;
- vertical 1 GeV negative muons at x = 0, ±200, ±500 and ±650 mm;
- zero injected SPTR and no electronics;
- seeds 26092601 and 8349041;
- diagnostics disabled for the production binary;
- individual detected-photoelectron `time_ns`, event identifier, face and SiPM
  identifier, plus event-level energy deposition and photon/PE counts.

The source snapshot contains both corrected behaviors: reflected optical states
survive the World navigation guard, and the air-to-wrap boundary uses
`CreateMylarReflector(0.98)`. The model had already passed the production,
first-encounter escape, incident-spectrum PDE, absorption-filtering and photon
ledger checks documented by v8. Because the stored ROOTs have all individual
END and TOP hit times, every requested timing observable can be recalculated.
The Phase-0 decision was therefore `NO_NEW_SIMULATION_REQUIRED`.

The complete audit is in
`presentations/v9/sources/dataset_audit.csv`; its machine-readable decision is
in `dataset_audit.meta.json`. `input_roots.tsv` records the exact 21 ROOT paths
and their recorded SHA-256 hashes. The extraction macro opened each ROOT,
required 10,000 `event_observables` entries, cross-checked the hit counts event
by event, and observed zero mismatches in all cells.

## Observable and fitting definitions

At each END, all detected photoelectrons from its eight SiPMs are pooled and
ordered by arrival time. The baseline uses the first photoelectron:

\[
t_L=t_L^{(1)},\qquad t_R=t_R^{(1)},\qquad
T_0=\frac{t_L+t_R}{2},\qquad \Delta t=t_R-t_L.
\]

The timing model is

\[
t_L=t_0+\frac{L/2+x}{v_{\rm eff}}+\epsilon_L,
\qquad
t_R=t_0+\frac{L/2-x}{v_{\rm eff}}+\epsilon_R.
\]

Thus the first-order x dependence cancels in T0, while

\[
\Delta t\simeq-\frac{2x}{v_{\rm eff}},\qquad
v_{\rm eff}=-\frac{2}{d\langle\Delta t\rangle/dx},\qquad
\sigma_x=\frac{v_{\rm eff}}{2}\sigma_{\Delta t}.
\]

All displayed histograms and fits were built with CERN ROOT 6.40.02. A first
Gaussian pass uses the median ± twice `(q84-q16)/2`; a second ROOT `TF1`
Gaussian is fit over the resulting mean ± 2 sigma. The adopted timing
resolution is this central Gaussian sigma. The global RMS over all events is
retained separately as a tail-sensitive width. The representative fits support
that choice: EJ-230 at x=0 gives chi2/ndf = 160.4/137 = 1.17; at x=+650 mm it
gives 266.5/136 = 1.96. The latter is shown in backup to expose the increasing
non-Gaussianity near the bar end.

Npe uncertainties are event-level SEM. Timing-width uncertainties are ROOT fit
errors on sigma. The RMS uncertainty is the normal-width approximation
`RMS/sqrt(2(N-1))`. The longitudinal-resolution uncertainty propagates the
Gaussian `sigma_DeltaT` fit error and the formal linear-slope error.

## Main physical results

All values in this section were recomputed from the current corrected-transport
ROOTs unless explicitly labelled analytical.

### Photon production and light collection

At x=0 the mean deposited energy is 1.925–1.936 MeV. With the configured light
yields, the observed mean scintillation production is 19,254.5 for EJ-200,
20,070.4 for EJ-204 and 18,780.6 photons for EJ-230. This agrees with the napkin
expectation of about 20,000 photons from a muon crossing 1 cm of plastic.

| material | Npe left | Npe right | Npe/end | Npe TOP | Npe total | Npe total / produced |
|---|---:|---:|---:|---:|---:|---:|
| EJ-200 | 528.5 ± 1.3 | 528.3 ± 1.3 | 528.4 ± 1.3 | 5202.8 ± 12.7 | 6259.6 ± 15.3 | 32.542% |
| EJ-204 | 398.2 ± 1.0 | 398.1 ± 1.0 | 398.1 ± 1.0 | 4767.8 ± 11.8 | 5564.1 ± 13.8 | 27.747% |
| EJ-230 | 304.2 ± 0.8 | 304.3 ± 0.8 | 304.3 ± 0.8 | 4086.2 ± 10.0 | 4694.7 ± 11.5 | 25.021% |

The table is for x=0. Npe/end rises toward either end because the close END arm
becomes dominant: the full ranges are 528–1301 pe/end for EJ-200, 398–1273 for
EJ-204 and 304–1126 for EJ-230. This increase does not imply a better T0 because
the far arm remains the weak timing measurement.

### Propagation

For n=1.58, the analytical axial limit is c/n = 189.742 mm/ns and the minimum
center-to-END transit time is 700/(c/n) = 3.69 ns. The measured Delta-t slopes
and derived effective speeds are:

| material | slope [ns/mm] | v_eff [mm/ns] | formal fit error [mm/ns] | chi2/ndf | effective angle |
|---|---:|---:|---:|---:|---:|
| EJ-200 | -0.0109760 | 182.216 | 0.022 | 825.8/5 | 16.19 deg |
| EJ-204 | -0.0110101 | 181.652 | 0.022 | 447.9/5 | 16.79 deg |
| EJ-230 | -0.0110126 | 181.611 | 0.022 | 419.3/5 | 16.83 deg |

The speeds are 4.0–4.3% below c/n. The effective angle follows from
`cos(beta_eff)=v_eff/(c/n)` and is only a one-number path interpretation, not an
angle assigned to every photon. The large chi2 values occur because the
event-mean errors are sub-picosecond and reveal small curvature. Scaling the
formal speed errors by `sqrt(chi2/ndf)` gives about 0.28, 0.21 and 0.20 mm/ns,
respectively; the talk therefore presents v_eff as an approximate effective
slope rather than implying 0.01% model accuracy.

For EJ-230 at x=0, the fitted T0 mean is 4.0964 ± 0.0008 ns, about 0.41 ns above
the axial limit. The excess has the expected sign from transverse path
components, reflections and nonzero scintillation-emission time.

### Optical timing and longitudinal position

| material | sigma(T0), x=0 [ps] | global RMS, x=0 [ps] | sigma(T0) over x scan [ps] | sigma_x, x=0 [mm] | sigma_x over x scan [mm] |
|---|---:|---:|---:|---:|---:|
| EJ-200 | 76.49 ± 0.86 | 76.91 | 76.4–82.0 | 13.05 ± 0.15 | 13.0–14.6 |
| EJ-204 | 71.47 ± 0.80 | 72.44 | 71.5–84.3 | 12.43 ± 0.14 | 12.4–14.8 |
| EJ-230 | 67.14 ± 0.77 | 67.77 | 67.1–81.6 | 11.42 ± 0.13 | 11.4–14.8 |

EJ-230 has the best intrinsic optical timing and central position resolution,
even though it has the lowest integrated PE yield. The configured light yields
are 10,000, 10,400 and 9,700 photons/MeV; rise times are 0.9, 0.7 and 0.5 ns;
decay times are 2.1, 1.8 and 1.5 ns for EJ-200, EJ-204 and EJ-230. The
corresponding bulk attenuation lengths are 3.8, 1.6 and 1.2 m. Long attenuation
makes EJ-200 the light-collection winner, while the faster rise and decay of
EJ-230 sharpen its early photoelectron population. The measured ordering
therefore demonstrates that a single 1/sqrt(Npe) law is insufficient across
materials and positions.

### First-photon, order-statistic and mean-of-first-m estimators

For EJ-230 at x=0, the current ROOTs permit a same-event scan from k=1 to 20:

\[
T_0(k)=\frac{t_L^{(k)}+t_R^{(k)}}{2}.
\]

The central Gaussian width is 67.14 ± 0.77 ps at k=1 and the k-th-photon curve
reaches a shallow minimum, 66.00 ± 0.76 ps, at k=3. A second scan averages the
first m arrival times independently on each END before forming T0. It reaches a
clearer minimum at m=5: 57.59 ± 0.66 ps. Relative to the first-PE baseline this
is a 9.55 ps, or 14.2%, reduction in the fitted core width. Because both widths
come from the same events and no paired bootstrap was run, no uncertainty is
assigned to their difference. No position-dependent choice was used. By k=20
the k-th-photon width has grown to 86.59 ± 1.02 ps and its global RMS to 96.86
ps; for the mean-of-first-20 estimator the corresponding values are 65.18 ±
0.77 ps and 72.98 ps.

This establishes the first PE as the baseline and the arithmetic mean of the
first five PE per END as the best estimator actually tested in v9. A weighted
mean was not introduced because the optical-only data do not define validated
channel weights, and no END+TOP combination is claimed.

### Symmetry and fit-stability addendum

The apparent left-right asymmetry of the Gaussian timing curve was tested on the
same current ROOT events. For every one of the 21 cells,
`sources/t0_fit_diagnostics.root` now preserves the exact T0 histogram and ROOT
TF1 object used for the quoted resolution. Its TTree and the parallel
`t0_fit_diagnostics.csv` record fitted sigma and error, fitted mean and error,
chi2, ndf, fit status, final fit limits, global RMS and robust half-width
`q_w=(q84-q16)/2`. All 21 ROOT fit statuses are zero.

RMS and q-width uncertainties were obtained with 300 deterministic event-level
bootstrap replicas per cell. The seed series begins at 26091443. The mirrored
differences are defined as the +|x| width minus the -|x| width:

| material | |x| [mm] | Delta sigma_G [ps] | z_G | Delta RMS [ps] | z_RMS | Delta q_w [ps] | z_q |
|---|---:|---:|---:|---:|---:|---:|---:|
| EJ-200 | 200 | +3.89 ± 1.27 | +3.06 | +1.51 ± 0.78 | +1.94 | +1.77 ± 1.05 | +1.69 |
| EJ-200 | 500 | -0.64 ± 1.30 | -0.49 | +0.48 ± 0.77 | +0.63 | +0.82 ± 0.99 | +0.83 |
| EJ-200 | 650 | +2.72 ± 1.27 | +2.14 | +0.61 ± 0.79 | +0.76 | +0.33 ± 1.06 | +0.31 |
| EJ-204 | 200 | +0.08 ± 1.16 | +0.07 | -0.44 ± 0.71 | -0.61 | +0.29 ± 0.94 | +0.31 |
| EJ-204 | 500 | +1.30 ± 1.29 | +1.00 | -0.09 ± 0.82 | -0.11 | +1.27 ± 1.10 | +1.15 |
| EJ-204 | 650 | -1.18 ± 1.34 | -0.88 | +0.25 ± 0.89 | +0.28 | -0.44 ± 1.20 | -0.37 |
| EJ-230 | 200 | -0.06 ± 1.12 | -0.05 | +0.89 ± 0.71 | +1.26 | +0.90 ± 0.95 | +0.94 |
| EJ-230 | 500 | +1.91 ± 1.24 | +1.54 | +1.29 ± 0.87 | +1.48 | +2.16 ± 1.14 | +1.89 |
| EJ-230 | 650 | +0.38 ± 1.30 | +0.29 | +0.45 ± 0.92 | +0.49 | +1.31 ± 1.26 | +1.04 |

No RMS or robust-width difference reaches two standard deviations. There is no
coherent sign across material or position. Only the central Gaussian definition
produces differences above two sigma, both for EJ-200: 3.06 sigma at 200 mm and
2.14 sigma at 650 mm. At those same pairs the robust differences are only 1.69
and 0.31 sigma, and the RMS differences 1.94 and 0.76 sigma. The preregistered
decision logic therefore identifies the largest visible asymmetry primarily as
a fit-definition effect combined with finite statistics, rather than physical
left-right asymmetry.

The observable construction was checked independently. END IDs 0--7 are pooled
on the left and IDs 8--15 on the right; each list is sorted by the same arrival
time field, and T0 assigns the two arms equal weight. All 210,000 extracted
events agree exactly with the stored per-arm hit counts. Under x reflection,
left/right yields and first-time means exchange without a coherent residual;
the largest individual swapped comparison is 2.38 sigma among the 36 yield/time
checks. Normalized EJ-230 T0 distributions, centered at their own means and
overlaid with identical `[-0.30,+0.30] ns` range and 120-bin definition, agree
visually for all three mirrored pairs. No left/right construction defect was
found.

The main result remains the central Gaussian width because its representative
fits remain adequate and it has a direct fit uncertainty. The main slide now
states the robust cross-check. RMS and q-width curves are retained in backup and
do not replace the original result.

### TIR, reflector and TOP

The current corrected reflector A/B remains a valid quantitative illustration
of light recycling: for EJ-204 at x=0 and N=2000 per arm, the END-yield gain is
1.09896 ± 0.00852 in END-only and 1.01155 ± 0.00830 in EndTop. In the END-only
census the reflected geometry incurred about 21.8 times more air-to-wrap
encounters. These results support the physical picture of TIR-guided early
photons plus reflector-recovered photons with more encounters. The stored data
do not contain the complete per-track angular history needed to assign a timing
width independently to those two populations.

At x=0, TOP records 4086–5203 PE across 70 sensors. Its global first detected
photon has a fitted core width of 2.9–3.1 ps, but this is an ideal optical lower
bound selected among 70 position-local channels. It is not directly comparable
to END T0, which combines one first PE from each of two eight-sensor arms. No
validated TOP pulse/threshold estimator or END+TOP combination exists, so the
deck does not claim a detector-level gain from TOP.

## Electronics boundary

Every new timing result is optical-only. A detector result would require

\[
\sigma_{\rm detector}^2 \simeq \sigma_{\rm optical}^2+
\sigma_{\rm SPTR}^2+\sigma_{\rm frontend}^2+
\sigma_{\rm threshold}^2+\sigma_{\rm clock}^2+\cdots .
\]

There is no jointly validated SPTR, pulse shape, threshold and TDC model in the
current dataset. The presentation therefore reports neither an electronics-
convolved resolution nor an experimental total.

## Historical timing results intentionally omitted

The following v6/v7 numbers were inspected only to understand useful
observables and were not reused: 53.68 ps END, 15.20 ps TOP, 15.21 ps combined
BLUE, 177.6 mm/ns effective speed and 7.9 mm longitudinal resolution. They were
derived from transport preceding the two optical repairs. None appears as a
current result in v9. Any accidental numerical similarity would not restore
their validity.

## Metrics still unavailable

- A validated electronics-convolved timing resolution.
- A calibrated threshold/pulse estimator for TOP.
- A common, validated END+TOP time estimator at a fixed channel cost.
- A physically validated weighted-mean estimator; the only multiphoton
  combination tested here is the unweighted mean of the first m END arrivals.
- An experimental detector resolution derived from this optical result.
- A timing decomposition by full per-track TIR-guided versus
  reflector-recovered angular history.
- A fair channel-count optimisation across alternative TOP layouts.

The ideal global-first TOP observable is available and shown only as an optical
information bound. It is not promoted to a hardware performance number.

## ROOT products and reproducibility

The analysis commands were:

```bash
cd /home/rrios/ej200/presentations/v9
root -l -b -q 'macros/build_timing_dataset.C+' \
  2>&1 | tee build_timing_dataset.root.log
root -l -b -q 'macros/build_timing_summary.C+' \
  2>&1 | tee build_timing_summary.root.log
root -l -b -q 'macros/build_symmetry_diagnostics.C+' \
  2>&1 | tee build_symmetry_diagnostics.root.log
for name in geometry_current npe_vs_x delta_t_vs_x propagation_times \
  t0_fit_ej230_x0 t0_fit_ej230_x650 sigma_t0_vs_x sigma_x_vs_x \
  material_comparison order_scan npe_vs_sigma end_vs_top fit_grid_ej230 \
  mirror_overlays_ej230 width_estimators_ej230; do
  root -l -b -q "macros/${name}.C" >"figures/${name}.root.log" 2>&1
done
/home/rrios/ej200_deck_20260910/build_exec29_docs/tools/tectonic \
  talk_v9.tex --keep-logs --keep-intermediates
```

`rebuild_v9.sh` records the same sequence. The extraction produced 210,000
event rows in `sources/timing_events.root`. The fit products are in
`sources/timing_summary.root`, with numerical exports in
`timing_summary.csv`, `material_summary.csv`, `order_scan_ej230_x0.csv`,
`t0_fit_diagnostics.csv` and `symmetry_diagnostics.csv`. The exact 21 fit
histogram/function pairs are also copied into `t0_fit_diagnostics.root`.

Fifteen figure quartets were verified: each has an executable `.C` macro, a
readable `.root` sidecar, a `.pdf`, and valid `.meta.json`. The presentation was
compiled on `t0minidaq` with Tectonic 0.17.0. Compilation exited 0 and produced
24 pages: 16 main slides including the title and 8 backup slides. The log has no
overfull or underfull boxes. All pages were rasterized and visually inspected;
no missing figures, clipping, empty pages or illegible overlays were found.

Key derived-source hashes before final packaging were:

| artifact | SHA-256 |
|---|---|
| `sources/timing_events.root` | `d19a7c8ca602f19d428af1652890c8a6f6f383cf5fffa80d109bb726831340be` |
| `sources/timing_summary.root` | `e13d06c5e904de76ca3cf4678c9fe1696c82fd72d5be099c5ff71ac459d9d4fe` |
| `sources/timing_summary.csv` | `95d469ab802321f17fe1e778cfd1bfe732f05ae1b510c41ff98e1467b3a15c6e` |
| `sources/material_summary.csv` | `525b2a63508c42399dacb43e6f0f77bf551a69687a081017f2aaa1c02dc0bcf2` |
| `sources/order_scan_ej230_x0.csv` | `8658334f0e04cbf5e7928ecc05c22bba74c2b3fa9b02301a76dc22b08362a241` |
| `sources/t0_fit_diagnostics.root` | `4fb87f60265c044d48c1bd2eb25e901967b99299d1b3abac2f3ad3a7f04652bd` |
| `sources/t0_fit_diagnostics.csv` | `e1c8e1d7142731b99b23341bf7b5f65a9438c94a8145f777c288915d39e6b2ec` |
| `sources/symmetry_diagnostics.csv` | `b40daf4745482682fdd8b64d336fc09cfd9367279bbdb06eb2eae035181b9ace` |

No merge or push was performed.
