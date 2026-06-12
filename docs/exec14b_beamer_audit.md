# EXEC_14B Beamer audit

## Initial state

- Initial HEAD: `be6743ef0328961a19dc063ad7d709fd2ab17a2d`
- Checkpoint tag: `checkpoint/pre-exec14b-beamer-fix-20260613-0041`
- Audit date: 2026-06-13 CEST
- Audited PDF: `results_ej230_analysis/report/exec13_ej230_report_full.pdf`
- PDF pages: 119
- Pages containing literal `figs/`: 104
- Raster images reported by `pdfimages -list`: 0
- Image references audited: 156
- References resolving to an existing canonical asset: 135
- Broken dispersion references: 21
- Forbidden EJ-204 data-figure references: 4 references to 3 distinct copied files

## Preliminary root causes

1. The report uses `\graphicspath{{../figs/}}` while every external image reference
   already starts with `figs/`. From `report/`, TeX therefore searches the wrong
   locations (`report/figs/...` and `../figs/figs/...`).
2. The 21 dispersion panels contain literal residual replacement backreferences
   (`\1` and `\2`) produced by `scripts/make_beamer_ej230.py`.
3. The report was compiled in draft mode, allowing missing images to become filename
   boxes rather than stopping the build.
4. `results_ej230_analysis/tables/exec14_macros.tex` is not imported by the report,
   contains duplicated macro names for different thresholds, and the visible report
   still contains hardcoded empirical results.
5. Existing EXEC_08b and EXEC_09 figures have no corresponding EJ-230 special-simulation
   ROOT/CSV products in the results tree and require provenance validation or regeneration.

## Final state

Pending reconstruction and strict PDF preflight.
