# talk_v7 — corrected transport and measured light collection

This is a new English talk using the v6 Beamer class, theme, package list and inline source-note style. It presents the measured guard/reflector defects, disjoint terminal accounting, configuration-dependent light yield and the status of historical claims. It contains no corrected timing-resolution result.

Build only `talk_v7.tex`. `baseline_v6/` is the untouched copy of the v6 starting point, retained as an inactive historical scaffold. Its old figures/macros are not inputs to v7. The original `presentations/v6/` is unchanged.

All numerical text is in `results_macros_v7.tex`. `sources/value_provenance.csv` and `sources/provenance.json` record values, interpretation, configuration, N, seeds, commits, run commands and copied sidecars. Quoted unsupported historical numbers appear only in the audit section; missing run provenance remains explicitly missing. Existing untracked historical CSVs have file hashes and are not assigned an invented Git commit.

The figures are generated from archived sidecars by `scripts/build_evidence.py`; each has `.csv`, `.root`, `.meta.json` and `.pdf` files. Their ROOT/CSV contents are round-trip checked. Plot coordinates and axis ticks are generated from the registered figure data, not handwritten scientific numbers in the slides.

Canonical compilation is on t0minidaq with the existing Tectonic 0.17.0 installation previously used for EXEC_29. No LaTeX package was added. Build metrics, source/PDF hashes, page count and box warnings are recorded in `BUILD_VALIDATION.json` after compilation and in the external EXEC_32 report.

Numerical qualifications retained:

- The guard-only yields 82.304 and 395.394 pe/end both have N=500.
- 98.000104% is the V1 END-only metallic forward fraction; it is not relabeled as the EndTop B0 fraction.
- Multiplicities count encounters, not successful reflections or per-track angular histories.
- The TOP hit/scintillation-generation ratio is not unique interception.
- A zero terminal ledger residual is an accounting check, not independent validation of physics.
- All historical timing numbers remain pending corrected-physics validation. No value in v6 is replaced.
