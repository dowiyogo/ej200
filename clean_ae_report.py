from pathlib import Path
import pandas as pd
base=Path('analysis/track_mechanism_20260915/step6_v2/veff_rank'); report=base/'REPORT_VEFF_RANK_SCAN_20260917.md'; text=report.read_text().split('\n## AE1: theta_eff tables')[0]
def md(path):
 d=pd.read_csv(path); h=list(d.columns); out=['| '+' | '.join(h)+' |','|'+'|'.join(['---']*len(h))+'|']
 for _,r in d.iterrows(): out.append('| '+' | '.join('' if pd.isna(v) else (f'{v:.8g}' if isinstance(v,(float,int)) else str(v)) for v in r)+' |')
 return '\n'.join(out)
text += '\n\n## AE1: theta_eff tables\n\n`theta_eff = arccos(d_axial/path_length_mm)` for the threshold-crossing photon. The aggregate table uses an events-weighted mean of per-cell quantiles. The per-position table is listed separately. Reference: 35.4 degrees experimental.\n\n'+md(base/'ae1_theta_cfd_by_extreme.csv')+'\n\n### By position\n\n'+md(base/'ae1_theta_cfd_by_position.csv')
text += '\n\n## AE2: geometry table\n\nThe exact event-level `<path>/<d_axial>` was not persisted; it is `NOT_AVAILABLE`. The proxy `1/cos(theta_median)` and implied velocity are tabulated without compatibility calculations.\n\n'+md(base/'ae2_geometry_cfd.csv')
text += '\n\n## AE3-AE5: configuration references\n\n`fastic_measured` and `penarodriguez_shortened` are identical `(2,3 ns)` models. The two distinct pulse forms are `(2,3 ns)` and `(0.5,2 ns)`. The MUSIC/SAMPIC chain has no model in `pulse_models.py`; comparison with 155 mm/ns is indicative only. Geometry references: simulated bar 140 cm versus test-beam 150 cm; configured muon 1 GeV/c versus test-beam 2.5 GeV/c. Blondel references are 16.1 cm/ns at center and approximately 14 cm/ns for x<20 cm.'
text += '\n\n## W6 updated reference column\n\n`difference_from_experimental_reference_mm_ns` is simulated minus 155 mm/ns; no compatibility calculation is made.\n\n'+md(base/'w6_cfd_fits.csv')+'\n'
report.write_text(text)
