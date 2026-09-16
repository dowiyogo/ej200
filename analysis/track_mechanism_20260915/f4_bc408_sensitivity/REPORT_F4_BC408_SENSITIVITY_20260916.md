# EXEC_46 F4 — BC-408 optical-model sensitivity at EJ200 x=-650 mm

Date: 2026-09-16. Step 6 remains suspended.

Two authorized 10,000-event cells were run with the original EndTop geometry, OPSC-100 emission/timing/yield, Broadcom PDE, zero SPTR, seeds 26092601/8349041, and four workers. Only RINDEX and ABSLENGTH in isolated SSLG4 runtime copies changed. The existing constant-n cell is the paired baseline; it was not rerun.

Each new log contains the inherited Geant4 `mat031` warning that the OPSC-100 fractional masses sum to 0.998943. The identical warning is present in the existing baseline log; it was not introduced by the RINDEX/ABSLENGTH sensitivity change.

## Published input and our derivation

[Huggins, Johnson and Buckner](https://arxiv.org/abs/2608.13710) measured commercial PVT scintillators over 370–660 nm. For BC-408 their ODR model is n(lambda)=1.518+0.640 exp(-0.00423 lambda_nm). They report 28.23 +/- 2.7 mm at 372 nm and only 90% lower bounds at visible wavelengths, including >764 mm at 439 nm. They have not implemented group-velocity corrections.

Our derivative gives n(420 nm)=1.626296, n_group(420 nm)=1.818695, and v_group=164.839 mm/ns. Propagating the paper’s typical 0.0055–0.0075 weighted-index uncertainty as a fully correlated additive index uncertainty gives sigma(v_group)=0.498–0.680 mm/ns. This propagation is ours and is conditional because the paper does not publish the ODR coefficient covariance needed for a full n_group uncertainty.

For the configured 1 GeV muon (beta=0.995424), the same derivation gives a 420 nm axial cone-edge velocity of 129.629 mm/ns. The finite-beta edge moves from 37.458 deg at 370 nm to 40.174 deg at 660 nm, a 2.716 deg chromatic span. The identity between the cone edge and critical angle remains exact wavelength by wavelength in the beta=1 limit.

At 200 mm, replacing c/1.58 by the derived 420 nm group velocity changes the nominal transport time by 159.2 ps, 22.2 times the 7.183 ps common residual. This scale comparison is our derivation, not a result quoted by Huggins et al.

The MPT samples n(lambda) every 2 nm from 370 to 660 nm. It is held constant at n(370) from 200–370 nm and at n(660) from 660–800 nm. The UV clamp is an explicit hypothesis: no exponential extrapolation is used. ABSLENGTH is 28.23 mm through 372 nm, interpolates to the visible scenario at 439 nm, then is held at either 764 or 3800 mm. The 764 mm scenario represents a lower bound, not a measured central value.

## F4 results

| scenario | first Cherenkov [%] | primary caustic q05–q95 [deg] | width [deg] | first-Cher local velocity [mm/ns] | first-scint local velocity [mm/ns] | first-scint total-track effective speed [mm/ns] | scintillator-MPT group speed at detected wavelengths [mm/ns] |
|---|---:|---:|---:|---:|---:|---:|---:|
| constant_n1p58_existing | 74.10 +/- 0.44 | 39.488–41.792 | 2.305 +/- 0.105 | 149.847 +/- 0.101 | 179.618 +/- 0.134 | 189.761 +/- 0.019 | 189.742 |
| visible_lower_764mm | 82.25 +/- 0.38 | 38.370–46.569 | 8.199 +/- 0.251 | 147.455 +/- 0.083 | 168.992 +/- 0.114 | 180.605 +/- 0.078 | 166.151 |
| visible_current_3800mm | 82.02 +/- 0.38 | 38.351–46.101 | 7.749 +/- 0.212 | 147.670 +/- 0.084 | 169.294 +/- 0.111 | 180.565 +/- 0.072 | 166.198 |

The published-dispersion MPT changes every requested observable. The difference between the two absorption scenarios is the unresolved attenuation systematic; neither scenario is selected as truth. The sub-370 nm Cherenkov region is unmeasured and is the dominant model limitation because the emitted spectrum scales approximately as 1/lambda^2.

The previous 74.10% first-photon Cherenkov fraction becomes 82.25% (764 mm) and 82.02% (3800 mm). The primary-cone central-90% width changes from 2.305 deg to 8.199/7.749 deg. This empirical quantile width is distinct from both the 0.02 deg histogram bin and the 2.716 deg parameter-free edge span.

The one-cell first-Cherenkov local direct velocity moves from 149.847 to 147.455/147.670 mm/ns. The first-scintillation local direct velocity moves from 179.618 to 168.992/169.294 mm/ns.

`path_length_mm` is `G4Track::GetTrackLength()` plus the detection step (`src/SiPMSD.cc`); it includes every medium crossed between creation and detection. Consequently, path divided by propagation time is a total-track effective speed, not the scintillator group velocity. The last table column instead evaluates Geant4’s RINDEX-derived group-velocity table at each detected wavelength. Their difference is therefore expected for tracks containing non-scintillator segments and is not an MPT closure failure.

The historical 148.548 mm/ns Cherenkov and 186.071 mm/ns scintillation numbers are seven-distance, origin-constrained slopes. A single x=-650 cell cannot refit those slopes. Likewise, the [SHiP test-beam 155 mm/ns value](https://doi.org/10.1016/j.nima.2020.164398) is a global effective propagation speed from time versus position, not a microscopic group velocity or one-point d/t ratio. F4 can test whether the local transport shifts toward that scale. Numerically, the distance of the local first-scintillation proxy from 155 mm/ns decreases from 24.618 to 13.992/14.294 mm/ns, so the shift is in the requested direction. It cannot establish agreement or trigger the stated full-campaign invalidation criterion by itself. That requires at least two positions under the corrected optical model.

## Decision

The current constant-n optical model is demonstrably non-robust for transport timing at the scale of the 7 ps residual: the measured-input sensitivity shifts local photon transport by far more than 7 ps-equivalent timing and introduces chromatic broadening that the old model fixes to zero. Therefore no Step 6 mechanism attribution is defensible with the constant-n campaign. The stronger claim that the old campaign is experimentally invalidated by reproducing 155 mm/ns is **not decidable from this one-position design**.

EJ-230 remains on its existing constant-n model: BC-420 was not measured and BC-422 is not a valid substitute. It is explicitly uncorrected.

## Reproducibility

```bash
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/run_f4_bc408_sensitivity.py
env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_f4_bc408_sensitivity.py
```

The exact MPT tables, macros, logs, ROOT hashes, run times and commands are under `/home/rrios/exec46_20260915/f4_bc408_sensitivity/`. Figure sidecars are `f4_metrics.{csv,root,meta.json}`; source-selection caches are retained under `scratch/`. No Step 6 analysis, push, merge, or deck edit occurred.
