#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p sources figures logs
root -l -b -q 'macros/analyze_tsum_veff.C+' 2>&1 | tee logs/analyze_tsum_veff.root.log
root -l -b -q 'macros/extended_models.C+' 2>&1 | tee logs/extended_models.root.log
root -l -b -q 'macros/summarize_decision.C+' 2>&1 | tee logs/summarize_decision.root.log
root -l -b -q 'macros/collapse_goodness.C+' 2>&1 | tee logs/collapse_goodness.root.log
if [[ "${SKIP_SPECTRAL:-0}" != 1 ]]; then
  root -l -b -q 'macros/spectral_check.C+' 2>&1 | tee logs/spectral_check.root.log
elif [[ ! -r sources/spectral_check.root ]]; then
  echo 'SKIP_SPECTRAL=1 but sources/spectral_check.root is absent' >&2
  exit 1
fi
for macro in macros/*.C; do
  name="$(basename "$macro" .C)"
  case "$name" in
    analyze_tsum_veff|collapse_goodness|extended_models|figure_common|spectral_check|spectral_figures|summarize_decision) continue ;;
  esac
  root -l -b -q "$macro" >"logs/${name}.root.log" 2>&1
done
