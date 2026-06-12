# Visual QA — exec12tb_beamer.pdf — 2026-06-12

Each slide rendered at 120 DPI with pdftoppm. Checked for:
- correct figure-title mapping
- no recycled figures
- legible axes and labels
- no text cutoff
- no empty slides

| Slide | Title | Content | Status |
|-------|-------|---------|--------|
| 01 | Title | Title slide with author/affiliation | OK |
| 02 | Scientific questions | 6 numbered questions as questions (not conclusions) | OK |
| 03 | Detector geometry | fine_scan_geometry.pdf — bar+scan diagram | OK |
| 04 | Datasets | dataset_inventory table + efficiency note | OK |
| 05 | From hits to observables | order_statistic_schematic.pdf | OK |
| 06 | Timing estimators | Formulas Δt and t+ | OK |
| 07 | 20th hit ≠ 20% CFD | cfd_vs_orderstat_schematic.pdf — two panels | OK |
| 08 | Mean Δt vs x | mean_dt_vs_x_with_residuals.pdf + inline table | OK |
| 09 | Timing width | rms68_dt_vs_x_4_20.pdf | OK |
| 10 | A-B correlation | rho_ab_vs_x_4_20.pdf | OK |
| 11 | t+ width | tplus_rms68_vs_x_4_20.pdf | OK |
| 12 | Δt at REF1 | ref1_dt_overlay.pdf | OK |
| 13 | Δt at REF2 | ref2_dt_overlay.pdf | OK |
| 14 | x_rec at REF1 | ref1_xrec_overlay.pdf | OK |
| 15 | x_rec at REF2 | ref2_xrec_overlay.pdf + inline summary table | OK |
| 16 | Resolution vs threshold | sweep_rms68x_vs_k.pdf | OK |
| 17 | Slope and chi2 vs threshold | sweep_slope_vs_k.pdf + sweep_chi2_vs_k.pdf | OK |
| 18 | Bias vs threshold | sweep_bias_vs_k.pdf | OK |
| 19 | Trade-off scatter | tradeoff_resolution_vs_bias.pdf | OK |
| 20 | Threshold sweep table | threshold_sweep.tex (generated) | OK |
| 21 | 4PE vs 20PE table | threshold_4_20_comparison.tex (generated) | OK |
| 22 | Does 20PE help? | Two-column alert/example blocks with macros | OK |
| 23 | Window-dip NPE | window_dip_counts.pdf | OK |
| 24 | Window-dip timing | window_dip_t4_t20.pdf | OK |
| 25 | Global context | global_context_sigma_bias.pdf | OK |
| 26 | Lv et al. comparison | lv_comparison.tex (qualitative, NOT-VALID flag) | OK |
| 27 | Timing budget | Text — absent effects enumerated | OK |
| 28 | Conclusions | 7 numbered conclusions with macros | OK |
| 29 | Next steps | Two-column blocks | OK |
| 30 | Backup header | "Backup slides follow" | OK |
| 31 | B1 Data inventory | Text + validation checklist | OK |
| 32 | B2 Event builder | Text — test list | OK |
| 33 | B3 Fit QA | Text — procedure | OK |
| 34 | B4 Reference positions | Inline table with LOO bias/sigma | OK |
| 35 | B5 Threshold sweep | Inline numeric table with macros | OK |
| 36 | B6 LOO methodology | Text — procedure | OK |
| 37 | B7 Window-dip metadata | window_dip_summary.tex | OK |
| 38 | B8 Channel map | Text description | OK |
| 39 | B9 Limitations checklist | Two-column ✓/× list | OK |
| 40 | B10 Reproducibility | Commands + git info | OK |

## QA Gates
- [x] 40 pages rendered without errors
- [x] 20 figures, 0 SHA duplicates
- [x] Zero LaTeX errors (!), zero Overfull hbox
- [x] All \includegraphics targets exist in figures/
- [x] generated_numbers.tex: 21 macros, all verified vs CSV source
- [x] No slide uses recycled figure under wrong title
- [x] Efficiency note on Frame 4 (not plotted as flat line)
