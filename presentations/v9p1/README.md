# Current optical-performance presentation (v9p1)

This deck is the Gate-1 correction of v9 after the definitive timing-symmetry
diagnosis. The main END timing metric is the robust central width
`sigma68=(q84-q16)/2`; Gaussian widths remain core-shape diagnostics. No new
Geant4 simulation was run. The original `presentations/v9/` is unchanged.

Build the PDF from the already generated ROOT figures with:

```bash
cd /home/rrios/ej200/presentations/v9p1
./rebuild_v9p1.sh
```

The complete ROOT analysis and its reproducible macros live in
`/home/rrios/ej200/analysis/timing_symmetry_20260914/`. Small diagnostic ROOT
and CSV artifacts are copied into `sources/` for the deck. The event-level
210,000-event extraction remains in `presentations/v9/sources/timing_events.root`
and is not duplicated.

Every diagnostic figure has a `.C`, `.root`, `.pdf`, and `.meta.json` artifact.
The external report is `/home/rrios/ej200/docs/reports/REPORT_TIMING_SYMMETRY_20260914.md`.
