# Tsum and distance-dependent first-arrival response

This ROOT-only analysis tests mean `T0(x)`, `Tsum(x)`, LEFT/RIGHT symmetry and
the linearity of the effective first-arrival response `g(d)` using the existing
21 corrected-transport cells. It does not run Geant4 and does not modify v9.

Rebuild the numerical artifacts and figures with:

```bash
cd /home/rrios/ej200/analysis/tsum_veff_20260914
./rebuild.sh
```

The spectral pass reads four branches from the 164 GB production dataset and
can take several minutes. Set `SKIP_SPECTRAL=1` to reuse an existing
`sources/spectral_check.root`; the script aborts if that artifact is absent.

The external report is
`/home/rrios/ej200/docs/reports/REPORT_TSUM_VEFF_POSITION_20260914.md`.
