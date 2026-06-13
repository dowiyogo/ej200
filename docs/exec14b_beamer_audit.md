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

- Canonical paths: `\graphicspath{{../}}` with every external reference written
  as `figs/<asset>.png`.
- All 21 corrupted dispersion references were replaced explicitly by the proper
  `top`, `endL`, and `endR` asset names for the seven key positions.
- Asset audit: `156/156` references found, zero missing/invalid/backreference
  errors, and zero forbidden EJ-204 data-figure copies.
- Numeric audit: 119 frames import 50 generated CSV-derived tables; regeneration
  comparison reports zero differences. Macro names are alphabetic and unique,
  including distinct `Nfour` and `Ntwenty` thresholds.
- Frame parity: `119/119` frametitles in the EXEC_12 structural order. Asset
  parity: 156 references and 134 unique assets on both sides.
- Strict compile: exit code 0 with `latexmk -pdf -halt-on-error`; no missing
  files, draft setting, undefined commands, duplicate commands, missing
  characters, or overfull boxes.
- Final PDF preflight: 119 pages, 104 pages with raster images, 156 page-image
  references, zero literal `figs/`, zero filename boxes, and zero suspicious
  empty image boxes. `pdfimages -list` detects embedded images.
- Visual render: all 119 pages rendered at 130 dpi; required pages 5--8, 17--20,
  35, 37, 39, 41--49, 113, 115, 117, and 119 inspected successfully.
