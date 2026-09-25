#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Read-only ROOT re-analysis; no Geant4 command is invoked.
analysis=/home/rrios/ej200/analysis/timing_symmetry_20260914
(
  cd "$analysis"
  root -l -b -q 'macros/timing_fit_stability.C+' 2>&1 | tee logs/timing_fit_stability.root.log
  root -l -b -q 'macros/mirror_diagnostics.C+' 2>&1 | tee logs/mirror_diagnostics.root.log
  root -l -b -q 'macros/summarize_diagnostics.C+' 2>&1 | tee logs/summarize_diagnostics.root.log
  for name in fit_grid_EJ200 fit_grid_EJ204 fit_grid_EJ230 mirror_EJ200 mirror_EJ230 width_methods_EJ200 width_methods_EJ204 width_methods_EJ230 fit_stability_EJ200 mean5_EJ230 sigma68_vs_x; do
    root -l -b -q "macros/${name}.C" >"logs/${name}.root.log" 2>&1
  done
)
cp "$analysis"/figures/* figures/
cp "$analysis"/sources/{timing_fit_stability.csv,timing_fit_stability.root,timing_widths.csv,mirror_diagnostics.csv,mirror_diagnostics.root,symmetry_metrics.csv,symmetry_metrics.root} sources/
/home/rrios/ej200_deck_20260910/build_exec29_docs/tools/tectonic talk_v9p1.tex --keep-logs
