# sigma_t source inventory — EXEC_33 gate

## Import boundary and identity

MSI `/home/reriosto/SHiP/analysis_core/` is a separate Git repository, branch `master`, source commit `295fed4e7409fc910792b188f250a1957f850d26`. Its tracked source was clean; untracked products were present. The exact remote inventory, selection, byte sizes and selected-file SHA-256 values are in [provenance/msi_inventory.json](provenance/msi_inventory.json) and [provenance/import_manifest.json](provenance/import_manifest.json).

The first import commit copies 91 files without content changes: all Python source, YAML configuration, original README, and historical out/ CSV/JSON/Markdown records. Generated ROOT/figures, LaTeX build products, caches and Git administration were inventoried but not imported. The source README remains intact under `upstream/analysis_core/README.md`; this directory's README is new documentation. No estimator was rewritten or run on a new simulation sample.

The requested estimator family extends beyond `analysis_core`. A second unchanged import preserves END SUM4/mirror code, optimization code and electronics sources from simulation main `420addf`, and EXEC_12T/EXEC_11 source from `b0aaac1` on its own historical lineage. See [provenance/related_import_manifest.json](provenance/related_import_manifest.json). `phase_ab_optimal.csv` was present as an untracked historical artifact in the t0minidaq clone; its file SHA is recorded, not assigned an invented Git commit.

## Runtime dependencies

Verified on MSI: Python **3.12.12**, ROOT/PyROOT **6.36.10**, NumPy **2.4.6**, uproot **5.7.4**, SciPy **1.17.1**, Matplotlib **3.10.9**, PyYAML **6.0.3**. Pandas is additionally imported by reporting and optimization scripts. C++ macros require ROOT/Cling with standard C++ headers. Beamer-generation scripts additionally require LaTeX; they do not estimate sigma.

`python3` on the SSH PATH resolves to Python 3.9.25 and fails to import this PyROOT build. Use **python3.12**, not an unqualified python3. This mismatch was observed, then the Python 3.12 dependency import and `timing_fit_pipeline.py --help` succeeded. No package or environment was modified.

[provenance/source_inventory.json](provenance/source_inventory.json) enumerates imports, every function and line, and all literal absolute paths in the imported Python source. Static parsing/compilation succeeded for all 32 imported Python files; this does not establish data-level correctness.

## Execution order, inputs and outputs

These are several live workflows, not one sequential pipeline. Paths below are the unmodified upstream paths. Many scripts hard-code MSI paths and historical N/geometry; directly running their main functions is not a new-pilot recipe.

| Entry / stage | Order and dependencies | Configuration / input | Output |
|---|---|---|---|
| `timing_fit_pipeline.py:804` | YAML → QA-0 branch check → ROOT/schema check → group/threshold timestamps → Gaussian fit/bootstrap → sidecars → plots/summary | `config/exec16_config.yaml`; `--config`, `--materials`, `--outdir`; historical filename builder at :60 includes literal `_300ev_`; ROOT `sipm_hits` branches `event_id,face_type,global_id,time_ns` and configured schema | per material/group/N ROOT histograms/functions/graphs, results CSV, metadata JSON, summary and plots under `outputs/` |
| `lib/gates.py:23` | Called by primary pipeline before/after fitting | branch-tip SHA, ROOT branches, fit bounds, group chi-square limits | aborts/flags; output hash manifest. The SHA gate checks a branch ref, not provenance embedded in the ROOT |
| `lib/robust_seeds.py` | Called by fit engine | finite timestamp vector; sqrt_n or FD binning | histogram peak, MAD-derived seed, bin count, diagnostics |
| `lib/fit_engine.py:193` | ROOT fit then SciPy bootstrap | fit window factor, min events, binning, bootstrap count/seed | sigma_fit/error, bootstrap error, chi2/ndf, status, histogram and TF1 |
| `lib/sidecar.py` | After per-position fits | result dictionaries and input paths/hashes | `.root`, `.csv`, metadata JSON |
| `regenerate_figures.py:362` | After primary pipeline | `outputs/*_results.csv`, ROOT histograms; `--material`; hard-coded historical labels | regenerated figures and audit summaries; no replacement fit |
| `cp5_final_figures.py:673` | After primary pipeline | same CSV/ROOT products, material constants | comparative overlays, matrix, material/fit plots, manifest |
| `exec16_endtop_ej204.py:1041` | discover positions → CP0 → L0/L1 → L2 → L3/L4 → JSON/report | hard-coded EndTop historical run directory, N_TOP=70, constants in module; fit_engine | `out/EXEC_16/L0..L4`, `results.json`, `MORNING_REPORT.md`; some stages combine arms, so do not run this main for EXEC_33 |
| `exec17_validation.py` | V1..V9 validation stages on raw ROOT and EXEC_16 products | module constants and historical paths | `out/EXEC_17/V*`, verdict/report/figures |
| `exec17_c0_patch.py:426` | C0_V5 → conditional stop → C0_V7 → C0_V2 | validation output `V7_scaling_summary.csv` plus raw ROOT | corrected validation CSV/JSON/figures, `c0_patch_summary.json` |
| `exec17_corrections.py:445` | C1 → C2 → C3 threshold sweep → C4 electronics table | EXEC_16 L2 CSV and raw positions; C3 N=[2,3,4,6,8] | C1..C4 CSV/figures and `results_corrections.json`; includes arm combination and deferred electronics, not an EXEC_33 main |
| `walk_correction.py:57,123` | fit_walk → apply_correction | timestamp and NPE vectors; median curve/model | fitted walk parameters, corrected timestamps; same-sample calibration in EXEC_18/19 |
| `exec18_main.py` | topology/label checks → walk studies → readout-model tables/report | raw historical EndTop ROOT; module constants; walk_correction and fit_engine | `out/EXEC_18/T*`, results/report/memo; historical Chain A/B pairing conflict preserved |
| `exec19_main.py:864` | integrity gate → T1 chains → T2 fraction scan → T3 topology → T4 reconciliation → canonical/report | same raw EndTop directory; first-hit, count-fraction CFD, per-channel weighting; walk/fit modules | `out/EXEC_19/T*`, canonical estimators, report/memo/results JSON |
| `exec20_main.py:632` | T1 END baseline → T3 TOP contribution → T4 reconciliation → T5 canonical | EndTop ROOT, EXEC_19 functions; optional newest END-only directory | `out/EXEC_20`; intrinsic first-hit averages and combinations, not SUM4 leading-edge |
| `exec21_diagnosis.py:660` | D1 falsification → D2 budget → D3 times → D4 theory → D6 checks | EndTop + separate END-only/dedicated raw ROOT and module constants | diagnosis, D* tables/figures, JSON |
| `exec22_analysis.py:55` | T3 ROOT analysis → optional latest T5 scan → comparisons/report | `out/EXEC_22/T3/x*mm`, N_EVENTS_T3=500, hard-coded END-only paths | first-hit END-average profiles and diagnosis outputs |
| `beamer/EXEC_*/build_beamer*.py` | After corresponding results JSON/CSV | historical result files | TeX/deck artifacts; no sigma estimator |
| `related/420addf/analysis/congruent_sum4_timing.C:323` | sequential ROOT files → per-event waveform threshold → END delta → FitCore → CSV/plots/hooks | inputDir, outputDir, nEvents; pulse/threshold/map are compile-time constants | `summary.csv`, plots, `hooks.txt`; no complete pilot ROOT/meta sidecar set |
| `related/420addf/analysis/tb_mirror_sigma_vs_x.C:249` | audit CSV → same END leading-edge → peak-seeded fit and NPE proxy cuts | audit column 0=x, 2=ROOT path, 6=OK; fixed 10000 event slots | `sigma_vs_x_End.csv`, maps/figures; distinct fit from FitCore; no full pilot sidecar set |
| `related/b0aaac1/analysis/exec12t_timing_threshold_analysis.py` | inventory/cache → reproduce EXEC_11 → analyze → window study; stage CLI in source | fixed two TOP channels, pairscan and earlier EXEC11/12 products; imports sibling exec11_pair_analysis.py | NPZ order-statistic cache, threshold/calibration/LOO-position CSVs, figures |
| `related/b0aaac1/analysis/exec12t_make_products.py` | After EXEC_12T analysis | `--out` and prior CSVs | reports/tables/deck; no new timing estimator |
| `related/420addf/analysis/optim/phase_ab.py` | per material/position raw ROOT → END k/m scan → same-sample argmin | hard-coded END-only ROOT directory; k<=250,m<=100, occupancy-limited | `phase_ab_kscan_full.csv`, `phase_ab_optimal.csv`, plots; no sigma uncertainty |
| `related/420addf/analysis/optim/phase_sparse_top.py` | raw ROOT → END m scan + TOP k scan → argmins → BLUE | sparse TOP layouts 4/8/14/20, k_TOP<=50; hard-coded paths | full/optimal CSV, combined quantities. **Do not execute:** it computes prohibited BLUE internally |
| `related/420addf/analysis/timing/sipm_waveform_dcfd.py/.cpp` | waveform/SPTR reconstruction → waveform fraction crossing → fits/output | explicit pulse-model, transit/electronics arguments; pulse_models.py registry | waveform/dCFD-specific outputs; independent live electronics path, not the SUM4 leading-edge implementation |

Raw ROOT files and historical output directories remain at their original locations. The imported records are audit evidence only. No historical sigma is a validation reference for the new pilot.
