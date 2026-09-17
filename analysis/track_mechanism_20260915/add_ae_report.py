from pathlib import Path
import pandas as pd
base=Path('analysis/track_mechanism_20260915/step6_v2/veff_rank'); report=base/'REPORT_VEFF_RANK_SCAN_20260917.md'
def md(path):
 d=pd.read_csv(path); h=list(d.columns); out=['| '+' | '.join(h)+' |','|'+'|'.join(['---']*len(h))+'|']
 for _,r in d.iterrows(): out.append('| '+' | '.join('' if pd.isna(v) else (f'{v:.8g}' if isinstance(v,(float,int)) else str(v)) for v in r)+' |')
 return '\n'.join(out)
text=report.read_text()
text += '\n\n## AE1: theta_eff tables\n\n`theta_eff = arccos(d_axial/path_length_mm)` for the threshold-crossing photon. The reported aggregate is a weighted mean of per-cell quantiles, weighted by events used. The reference is experimental only.\n\n'+md(base/'ae1_theta_cfd_by_extreme.csv')+'\n\n## AE2: geometry table\n\nThe exact event-level `<path_length_mm>/<d_axial_mm>` was not persisted in the partials and is `NOT_AVAILABLE`. `path_over_daxial_from_theta_median = 1/cos(theta_median)` is a proxy derived from the persisted theta median; `v_implied = v_group_MPT / proxy`.\n\n'+md(base/'ae2_geometry_cfd.csv')+'\n\n## AE3-AE5: configuration references\n\nThe two `(2,3) ns` models, `fastic_measured` and `penarodriguez_shortened`, are identical in the generated fits. The distinct pulse-shape range is the interval between `(0.5,2.0) ns` and `(2,3) ns`; for EJ-200 at CFD 14%, the endpoints are 175.163556 and 171.999485 mm/ns. No model represents MUSIC/SAMPIC, so comparison with 155 mm/ns is indicative only. The simulated bar is 140 cm versus 150 cm for the test beam; configured muon momentum is 1 GeV/c versus 2.5 GeV/c in the test beam. Blondel reference values are 16.1 cm/ns at the center and approximately 14 cm/ns for x < 20 cm.\n\n## W6 updated reference column\n\n`difference_from_experimental_reference_mm_ns` is the simulated value minus 155 mm/ns; no compatibility calculation is made.\n\n'+md(base/'w6_cfd_fits.csv')+'\n'
report.write_text(text)
