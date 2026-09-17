from pathlib import Path
import pandas as pd

BASE = Path('analysis/track_mechanism_20260915/step6_v2/veff')
REPORT = BASE / 'REPORT_VEFF_TIMEDIFF_20260917.md'

def render(path, columns=None, digits=8):
    df = pd.read_csv(path)
    if columns:
        df = df[[c for c in columns if c in df.columns]]
    headers = list(df.columns)
    lines = ['| ' + ' | '.join(headers) + ' |', '|' + '|'.join(['---'] * len(headers)) + '|']
    for _, row in df.iterrows():
        vals = []
        for key, value in zip(headers, row):
            if pd.isna(value): vals.append('')
            elif isinstance(value, (float, int)):
                fmt = '.8e' if 'cubic_a3' in key else f'.{digits}f'
                vals.append(format(value, fmt))
            else: vals.append(str(value))
        lines.append('| ' + ' | '.join(vals) + ' |')
    return '\n'.join(lines)

u4 = pd.read_csv(BASE/'u4_linearity.csv')
residuals = u4[u4['x_mm'].notna()]
secants = u4[u4['x_low_mm'].notna()]
lines = [
'# EXEC_46 effective propagation velocity tables', '',
'All tables use the verified derived trees. The v2 tree SHA-256 is `764c643e6ccfc3d0b96c7fc20d6753e86139233599ca19359c8d95e585844290`; the old tree SHA-256 is `816c302c3ed37f7eb29855bcad7e734a1b6a77238f60a30361635b68377e341d`. No production ROOT was modified and no simulation was run.', '',
'## U1: time difference by cell', '',
'`mean_tdiff_ns` is the event mean of `tL-tR`; `sem_tdiff_ns` is its event SEM; `n_events` is the event count.', '', render(BASE/'u1_tdiff_cells.csv'), '',
'## U2: weighted slopes', '',
'Each material/campaign fit uses the seven position means, x in mm, and weights 1/SEM^2. `chi2_flag=*` marks chi2/ndf > 5.', '', render(BASE/'u2_slope_fits.csv'), '',
'## U3: effective velocity conventions', '',
'`v_eff_two_end` is 2/abs(d<tL-tR>/dx). `v_eff_one_end` is 1/abs(d<tL>/dx) or 1/abs(d<tR>/dx), respectively. Errors are first-order propagated from slope errors. No convention is selected.', '', render(BASE/'u3_veff.csv'), '',
'## U4: linearity and local secants', '',
'For point rows, `linear_residual_ps` is observed minus weighted linear prediction. The cubic fit is y=a1*x+a3*x^3 with no intercept; its chi2/ndf is tabulated separately. Asterisk flags retain the chi2/ndf > 5 rule.', '', render(BASE/'u4_linearity.csv'), '',
'Local secant rows use adjacent positions. `secant_slope_ns_per_mm` is the difference-mean change divided by position change. The two effective-velocity columns apply the two conventions to that secant.', '', render(BASE/'u4_linearity.csv', ['campaign','material','x_low_mm','x_high_mm','secant_slope_ns_per_mm','v_eff_two_end_mm_ns','v_eff_one_end_mm_ns']), '',
'## U5: contrast table', '',
'`experimental_reference_mm_ns=155` is the reported Betancourt et al. reference value, shown without differences or compatibility calculations.', '', render(BASE/'u5_contrast.csv'), '',
'## U6: old/new control', '',
'The same U1-U5 tables contain both `old` and `v2` rows. EJ-230 numeric rows are identical between campaigns in all five CSV outputs; the maximum numeric difference is 0.', '',
'## Sources and units', '',
'- `u1_tdiff_cells.csv`: derived-tree event means and SEMs.',
'- `u2_slope_fits.csv`: weighted seven-position linear fits.',
'- `u3_veff.csv`: separate tL, tR and tL-tR fits with propagated errors.',
'- `u4_linearity.csv`: point residuals, odd cubic fits and adjacent-position secants.',
'- `u5_contrast.csv`: two-end and one-end velocities, linear chi2/ndf, local range and experimental reference.',
'- Time quantities are ns; slopes are ns/mm; velocities are mm/ns; residuals are ps.', ''
]
REPORT.write_text('\n'.join(lines))
