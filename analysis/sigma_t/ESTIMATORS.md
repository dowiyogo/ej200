# Independent TOP and END definitions — decision register

Paths starting `core/` below mean `upstream/analysis_core/`; `related/` means `upstream/related/`. Imported sources remain byte-identical. No choice is made here among incompatible historical definitions.

## TOP: photon order, channel selection and reporting are different axes

| Code name and exact source | Grouping / timestamp | Indices and scope | What is actually reported |
|---|---|---|---|
| `compute_tN`, `extract_group_tN`, `TOP_SUM4`; core/timing_fit_pipeline.py:80,107,190 | Rank TOP face_type=2 channels by aggregate hits in the position file; take highest 4 global IDs (subset of 16..85). Merge their hits per event, sort, take time at zero-based N−1. Events with fewer than N hits are omitted | YAML `N_VALUES` at config/exec16_config.yaml:109 = [4,20]; same N values at every position; channel choice recalculated per position | Separate Gaussian sigma curve for each fixed N, not automatic min over N. Raw event denominator is unique IDs present in hit tree, so events with no hits anywhere are absent |
| `TOP_SUM4_N1`; core/exec18_main.py:217 and core/exec16_endtop_ej204.py:541 | Four TOP channels selected by aggregate NPE (not necessarily a contiguous physical cluster); first photon of their merged stream | photon index N=1; SUM4 means **4 channels**, not the 4th photon; per-position channel ranking | Fixed-N=1 fits; EXEC_18 adds walk correction in its walk study |
| `TOP_THRk`; core/exec16_endtop_ej204.py:89,545 | one nearest TOP channel, t_k | k=(1,2,3,5); per position | Separate fixed-k curves |
| `run_C3`; core/exec17_corrections.py:313–349 | aggregate-count-selected TOP_SUM4 and TOP_SUM8, merged stream t_N | N=[2,3,4,6,8], x=0,−690,+690; group and N optimized per position | `min(rows, key=...sigma_ps)` on the same fitted rows; reports the winning N and its same-sample minimum sigma. Procedure **(b)** |
| `order_statistics`, `HOOK_THRESHOLD_SWEEP`; related/b0aaac1/analysis/exec12t_timing_threshold_analysis.py:14,23,61,76 | **two separate TOP channels (28,29)**; cache each channel's ordered arrivals; derive Δt=t_A(k)−t_B(k), t+=(t_A(k)+t_B(k))/2 | actual code caches/scans **1..30**; primary 4,20; 41 fine positions −462..−422 mm; matched20 selection also reported | per-k position widths and mean over positions; width `sigma_core=(q84−q16)/2` (RMS68), not automatically a Gaussian sigma. No same-event train/evaluate split for choosing k. Leave-one-position-out linear calibration evaluates position reconstruction, not procedure (c) for timing |
| `sig_top`, `k_opt`; related/420addf/analysis/optim/phase_sparse_top.py:42,77,96,124 | all face_type=2 channels merged; timing values evaluated on events common to left/right END and TOP. This extra END coincidence affects the TOP sample | k=1..min(50,max(1,int(percentile(cT,2)))); occupancy bound uses all TOP event counts, timing mask uses common events. Per material/layout/position | IQR/1.349 width; `nanmin(sig_top)` and `nanargmin+1` on the same sample → **(b)**. Script then combines arms using BLUE; imported for audit only |
| `k_opt`, `m_opt`, `phase_ab_optimal.csv`; related/420addf/analysis/optim/phase_ab.py:35,50,54,63,95,104,126 | **END-only**, not TOP: pool each END face separately. k: average the two k-th hits; m: average the first m hits on each END, then average left/right | k<=250, m<=100; occupancy cap exactly `min(Kmax,max(1,pct2(cL),pct2(cR)))`; candidates require both END counts >= index and >=30 valid events. Each material/position separately | IQR/1.349; minimum over the same sample → **(b)**. m is a first-m **mean**, not an order statistic. CSV has no timing error column |

The EXEC_12T title and commit `f38c9e1` say “fourth-to-thirtieth”; the executable bounds are 1..30. The primary comparison is 4th vs 20th. Earlier `tn_order_statistics.C` elsewhere computes 1..20 per-channel means/RMS; it is not this EXEC_12T sweep.

## Three distinct reporting procedures

**(a) Fixed global index.** Choose one k once, freeze it and apply it to every position. The core pipeline can report fixed N curves and EXEC_16 TOP_SUM4 uses fixed N=1. These mechanisms support fixed-index evaluation, but the supplied code does not establish an independent global selection dataset/objective. A curve evaluated at fixed k is not evidence that k was chosen without looking at that curve.

**(b) Optimum per position on the selection sample.** Implemented explicitly in EXEC_17 C3, phase_sparse_top and phase_ab. Selection and reported width use the same event sample. The selected minimum is downward-biased by selection. A fit error or bootstrap conditional on that chosen k does not remove the selection bias.

**(c) Independent selection/evaluation event samples.** No such TOP index-selection/evaluation implementation was found in the inventoried families. It must be specified after René chooses it; do not relabel the existing bootstrap or EXEC_12T leave-one-position-out calibration as (c). Independence avoids the same-sample argmin optimism; it does not guarantee every finite-sample Gaussian/core-width estimator is mathematically unbiased. Any channel ranking and walk/fit calibration trained from data must also be frozen from training data, with the event split recorded.

`1c5a4c5`'s biased-argmin caveat actually concerns an **81-point w-weight scan** on the same n=5000 sample, not a TOP k scan. It is relevant evidence about selection reporting, not an implementation of (c). The warning at `presentations/v6/talk_v6.tex:952` compares END m*=7 with m=8 while changing N_TOP; it does not isolate a TOP-hardware effect. Both references were checked as source/history, not used to choose an estimator.

## TOP fit and uncertainty

For the core Gaussian branch, core/lib/fit_engine.py:193 uses histogram mode and MAD seed, then a ROOT Gaussian fit within mode ± `FIT_WINDOW_SIGMAS`×MAD-derived sigma. The YAML defaults are factor 2, sqrt_n bins, minimum 30 events, `R Q S 0`, bootstrap 300 with seed 20260618. EXEC_16/17/19 standalone scripts have their own defaults (including 200 bootstrap replicates and FD versus sqrt_n), so there is no single universal fit configuration.

Returned `sigma_fit`, `sigma_fit_err`, `bootstrap_err` are **ns**; multiply by 1000 for ps. The ROOT parameter covariance error and bootstrap error are separate. `_bootstrap_sigma` at :124 resamples the timestamp vector, retains the original fit window/seeds and uses SciPy with bins restricted to that window; it is not an independent holdout and does not reselect k. Fewer than 10 accepted bootstrap fits returns NaN. A finite sigma alone is not proof of a successful fit: inspect `fit_status`, `chi2_ndf`, flags and uncertainty.

## END SUM4 + leading-edge (independent arm)

`related/420addf/analysis/congruent_sum4_timing.C` implements:

- `GroupIndex` / `HOOK_MAP` (:216): fixed clusters {0,1,2,3}, {4,5,6,7} on left; {8,9,10,11}, {12,13,14,15} on right. This differs from the core pipeline's four highest-count channels.
- `kSprRiseNs=0.5`, `kSprFallNs=5.0` (:45): normalized difference of exponentials, peak amplitude 1 PE for one arrival; `Pulse` (:165). No SPTR/electronic jitter, walk correction or ToT cut in this stage. LSB=24 ps is metadata only, not applied quantization.
- `kThresholdPe=4.0` / `HOOK_THR` (:49): absolute **summed waveform amplitude**, 4 PE-equivalent. It is not the time of the fourth photon and not a 4%/20% CFD fraction.
- `LeadingEdgeTime` (:170): sort arrivals, accumulate slow/fast exponential states, locate first rising interval reaching threshold, bisect 60 times. Return NaN if no crossing.
- `Earliest` / `HOOK_ENDRED` (:228): first finite cluster crossing in each end. `left=Earliest(trigger[0],trigger[1])`, `right=Earliest(trigger[2],trigger[3])` (:361); accept event only if both finite; ΔT_LR=left−right (:363).
- `FitCore` (:110): min 20 events; 100 bins over median ±8 MAD-sigma; four iterative fits within current mean ±2 current sigma. Successful fit exports Gaussian width and TF1 parameter error in ps, chi2/ndf. Failed fit falls back to sample RMS and RMS/sqrt(2(n−1)), with chi2/ndf=−1 and usedFit=false. Fewer than 20 events leaves default zero fields; these are not valid measured resolutions.
- `summary.csv` (:384): `sigma_lr_end_ps`, `sigma_lr_end_err_ps` are width/error of ΔT. `sigma_single_end_ps` and its error divide both by sqrt(2). That normalization assumes an interpretation of equal independent end timing contributions; it is not sigma of (t_L+t_R)/2. No TOP/END combination is involved.

The macro also computes a **different TOP split A/B waveform observable** and draws a historical 88 ps line. That TOP output is not the requested TOP photon-index estimator, and this full macro must not be used as the unqualified EXEC_33 pilot command.

### Test-beam mirror is a second live END fit

`related/420addf/analysis/tb_mirror_sigma_vs_x.C:10` includes the same pulse/threshold implementation. Its `FitPeakSeeded` (:105) uses 200 bins, median±8 robust width with minimum robust seed 0.02 ns; it seeds at the peak, fits within peak±2 **ns** (not ±2 sigma), seed sigma=0.5 ns and limits 0.02..5 ns. It reports raw and NPE>=p25/p50 variants (:302), plus FWHM and NPE Landau diagnostics. Its event storage is fixed at 10000 (:263). On fit failure, fields can remain zero; this requires explicit validity handling in a future pilot adapter.

`presentations/sim_status_hi_2026-06-08/talk.tex:68,441` documents SUM4 + leading-edge and sigma(ΔT_LR)/sqrt(2); :422 is the historical 88 ps comparison. The source supports the waveform methodology but the two live fit functions differ. They remain distinct; René must choose which END fit/normalization/parameters apply. No historical value is a pilot acceptance criterion.

## EXEC_18 / EXEC_19 conflict — recorded, not adjudicated

core/exec18_main.py:217 defines TOP_SUM4_N1 as minimum first arrival across selected channels. core/exec19_main.py:207–208 implements that same raw timestamp **by calling `tN_stream(...,N=1)` on the merged selected stream**. It is therefore the k=1 member of the merged TOP order-statistic family at the raw timestamp stage, not a separate min-only algorithm. The complete Chain A path then applies walk correction and its own fit; it is not numerically identical to every raw N=1 report.

core/exec19_main.py:223 implements a **separate Chain B aggregation**: per-channel first-hit timestamps → per-channel walk correction → per-channel fitted inverse-variance weights → weighted mean over at least two valid channels → fit that mean. It is not a k=1 merged-stream order statistic, though it reuses the same primitive on each channel.

core/exec19_main.py:334 states that EXEC_18's named sigma_int_B=69.2 ps used the min-stream. The archived EXEC_19_REPORT.md:10,18,21 records 80.7 ps for the weighted Chain B at x=0 and describes the old pairing as wrong. Those are historical source claims only. Their distinction is documented without selecting either chain, recalculating them or declaring either the EXEC_33 pilot estimator. No BLUE or TOP/END combination is applied.

## Electronics versions remain separate

1. `analysis_core` first-hit/order-statistic fits are labeled intrinsic; YAML `HOOK_SPTR=false` is a stub, not a waveform injection. EXEC_17/18/19 constants use SPTR=106 ps and FastIC=10 ps in post-hoc quadrature tables; Chain B scales SPTR by sqrt(N_eff). Those formulas belong to those historical chains and are not silently attached to TOP order statistics.
2. The SUM4 END response uses a **0.5/5 ns normalized difference of exponentials**, 4-PE amplitude threshold and zero injected jitter. `HOOK_*` labels are explicitly provisional.
3. `9040783` introduced a separate named waveform/dCFD registry: shortened 2/3 ns, intrinsic unknown-rise/55 ns, FastIC-measured pending at that commit. It removed the mixed 2/55 ns default and required explicit pulse/SPTR/electronics parameters.
4. The live main snapshot also contains `7be59bd`, which later fills `fastic_measured` with 2/3 ns as bank-derived and still requiring full-chain confirmation, and adds 85/92 ps bank-derived SPTR entries alongside published-other-device/OV entries. This is later than 9040783; it is not retroactively attributed to that commit.
5. `SPTR_PROVENANCE.md` records 137±4 ps FWHM intrinsic and 172 ps FWHM detector for the cited related device/OV, converting using 1/2.355 (~58/~73 ps sigma). It explicitly says average-based sqrt(k*N) propagation is not established for order statistics. These are statements in the archived provenance document, not new verification of the paper or operating-point equivalence. The later electronics inventory retains the 106 ps FWHM/sigma/device ambiguity and the OV distinction.
6. `sipm_waveform_dcfd.py:287` jitters each arrival with a Gaussian transit sigma, constructs a sampled waveform with `(1-exp(-t/tau_r))*exp(-t/tau_f)`, and obtains a fraction crossing. At :319 it adds a quantity named electronics_sigma_ns directly to waveform amplitudes as Gaussian sample noise. This implementation detail is preserved; its parameter name/unit is not reinterpreted as a proven timestamp-jitter model.

The 0.5/5 ns END response and the named 2/3 or unknown/55 ns dCFD responses are live independent versions with different mathematical pulse parameterizations. No SPTR or pulse version is unified here.
