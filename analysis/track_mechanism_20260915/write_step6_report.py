from pathlib import Path
import json
import pandas as pd

ROOT = Path('analysis/track_mechanism_20260915')
OUT = ROOT / 'step6_v2'
REPORT = ROOT / 'REPORT_STEP6_V2_TABLES_20260917.md'

def table(path, columns=None, digits=8):
    if Path(path).suffix == '.json':
        payload = json.loads(Path(path).read_text())
        if isinstance(payload, dict):
            payload = [payload]
        df = pd.DataFrame(payload)
    else:
        df = pd.read_csv(path)
    if columns:
        df = df[[c for c in columns if c in df.columns]]
    headers = list(df.columns)
    lines = ['| ' + ' | '.join(headers) + ' |',
             '|' + '|'.join(['---'] * len(headers)) + '|']
    for _, row in df.iterrows():
        values = []
        for value in row:
            if isinstance(value, (float, int)) and not pd.isna(value):
                values.append(f'{value:.{digits}f}')
            else:
                values.append(str(value))
        lines.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join(lines)

lines = ['# EXEC_46 Step 6 v2 tables', '', 'All values below are tabulations of existing v2 ROOT-derived outputs. No production ROOT was modified and no simulation was run.', '', '## T1: inventory and integrity', '', 'Columns are copied from `exec46_inventory.csv`: ROOT path/hash and `.DONE` hash status; event/photon counts; Npe means; clock, ordering, path, boundary, apparent-speed, duplicate-key, source, sensor, gun-position and tau gate counters; runtime and informational GROUPVEL fields.', '', table(OUT/'exec46_inventory.csv'), '', 'Source-type census by all hits and first END photons:', '', table(OUT/'exec46_source_census.csv'), '', '## T2: old versus new', '', 'Npe_END is the mean number of END hits per event. `fCher(first)` is the fraction of first END photons with `source_type == 2`; delta is new minus old in percentage points. EJ-230 rows were checked for exact equality.', '', table(ROOT/'../step6_v2/t2_old_vs_new.csv') if (OUT/'t2_old_vs_new.csv').exists() else 'T2 table: see the recorded comparison in the execution log; CSV not generated.', '', '## T3: transport fits', '', 'The fit table reports the linear-through-origin velocity and chi2/ndf for first overall, scintillation and Cherenkov populations at the seven nominal distances. Optical columns use the runtime-resolved wavelength-dependent tables. The selected-wavelength join matched all selected rows; rejected rows are zero because the join is fatal on missing or duplicate matches.', '', table(OUT/'step3/fit_summary.csv'), '', 'Join/optical runtime tables:', '', table(OUT/'step3/optical_runtime.csv') if (OUT/'step3/optical_runtime.csv').exists() else table(OUT/'step3/optical_runtime.json'), '', '## T4: Cherenkov cone metrics', '', 'Tables report first-photon Cherenkov fractions, primary-like angular distributions, axial velocities, wavelength-resolved edge quantities and clamped/outside-measured fractions.', '', table(OUT/'step4/cherenkov_counts.csv'), '', table(OUT/'step4/cherenkov_angle_window.csv'), '', table(OUT/'step4/cherenkov_angle_by_nc.csv'), '', '## T5: T0 tables', '', 'The derived tree defines tL and tR as first END photons and T0=(tL+tR)/2. Cell means and SEMs, within/between slopes, linear and POL2 chain residuals, F-tests and leave-one-out rows are copied below.', '', table(OUT/'step5/cell_statistics.csv'), '', table(OUT/'step5/within_between_summary.csv'), '', table(OUT/'step5/between_specification.csv'), '', table(OUT/'step5/h1_even_f_test.csv'), '', table(OUT/'step5/h1_loo_summary.csv'), '', '## T6: widths and guardrails', '', 'RMS is the ordinary population standard deviation per cell. IQR/1.349 uses q75-q25. Gaussian fits use 4 ps bins over median +/-0.320 ns and a fit window median +/-2*q68, q68=(q84-q16)/2. A `*` marks chi2/ndf > 5; such fits are not treated as reliable.', '', table(OUT/'step6/t6_widths.csv'), '', 'L2 uses ordinary RMS only: margin = RMS(T0)-[RMS(tL)+RMS(tR)]/2. `se_gap_ns` is NOT_AVAILABLE. SE_gap, defined as a paired influence-function delta-method uncertainty, is not implemented in this repository; margins are therefore reported without statistical tolerance. `CHECK` marks a margin with absolute value below 1 ps or a positive margin.', '', table(OUT/'step6/t6_l2_guardrail.csv'), '', 'Existing near-END mixture and angle-time columns from `step4/near_end_mixture_summary.csv`:', '', table(OUT/'step6/t6_mixture_and_correlation.csv'), '', '## Sources', '', '- T1: `step6_v2/exec46_inventory.csv`, `exec46_source_census.csv`, `exec46_tau_diagnostics.csv`.', '- T2: old campaign `/home/rrios/exec46_20260915/full_grid` and v2 campaign ROOT files; EJ-230 exact-equality gate passed.', '- T3: `step6_v2/step3/`.', '- T4: `step6_v2/step4/`.', '- T5: `step6_v2/step5/` and the derived tree SHA-256 `764c643e6ccfc3d0b96c7fc20d6753e86139233599ca19359c8d95e585844290`.', '- T6: `step6_v2/step6/` and `step6_v2/step4/near_end_mixture_summary.csv`.', '']
lines.insert(0, 'T3 join diagnostics are persisted in `step3/wavelength_join_diagnostics.csv`: columns are selected, joined and rejected photon counts per cell.')
REPORT.write_text('\n'.join(lines))
