# Current optical-performance presentation (v9)

This deck reports timing, light collection and longitudinal-position metrics
recalculated from the 21-cell corrected-transport grid. It does not run a new
simulation. The common configuration is EndTop, 16 END and 70 TOP SiPMs,
vertical 1 GeV muons, zero injected SPTR, no electronics, seven x positions
and 10,000 events per cell.

Rebuild the ROOT extraction, fits, figures and PDF on `t0minidaq` with:

```bash
cd /home/rrios/ej200/presentations/v9
./rebuild_v9.sh
```

The final timing summary is `sources/timing_summary.csv`; material-level
propagation fits are in `sources/material_summary.csv`; the k-th photon and
mean-of-first-m scans are in `sources/order_scan_ej230_x0.csv`. Every displayed figure has a `.C`,
`.root`, `.pdf` and `.meta.json` quartet.

The symmetry addendum is recorded in `sources/t0_fit_diagnostics.csv` and
`sources/symmetry_diagnostics.csv`. The companion ROOT file preserves all 21
exact T0 histograms and Gaussian fit objects. RMS and robust-width uncertainties
use 300 deterministic event-bootstrap replicas.
