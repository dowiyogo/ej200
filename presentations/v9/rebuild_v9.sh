#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# The extraction macro reads only the 21 existing corrected-transport ROOT files.
# It does not launch Geant4 or generate events.
root -l -b -q 'macros/build_timing_dataset.C+' 2>&1 | tee build_timing_dataset.root.log
root -l -b -q 'macros/build_timing_summary.C+' 2>&1 | tee build_timing_summary.root.log
root -l -b -q 'macros/build_symmetry_diagnostics.C+' 2>&1 | tee build_symmetry_diagnostics.root.log

figures=(
  geometry_current npe_vs_x delta_t_vs_x propagation_times
  t0_fit_ej230_x0 t0_fit_ej230_x650 sigma_t0_vs_x sigma_x_vs_x
  material_comparison order_scan npe_vs_sigma end_vs_top
  fit_grid_ej230 mirror_overlays_ej230 width_estimators_ej230
)
for name in "${figures[@]}"; do
  root -l -b -q "macros/${name}.C" >"figures/${name}.root.log" 2>&1
done

/home/rrios/ej200_deck_20260910/build_exec29_docs/tools/tectonic \
  talk_v9.tex --keep-logs --keep-intermediates
