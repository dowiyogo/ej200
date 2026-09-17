from pathlib import Path
import pandas as pd
base=Path('analysis/track_mechanism_20260915/step6_v2/veff_rank'); report=base/'REPORT_VEFF_RANK_SCAN_20260917.md'
d=pd.read_csv(base/'ae1_theta_cfd_by_position.csv'); h=list(d.columns); lines=['\n## AE1: theta_eff by position\n','| '+' | '.join(h)+' |','|'+'|'.join(['---']*len(h))+'|']
for _,r in d.iterrows(): lines.append('| '+' | '.join('' if pd.isna(v) else (f'{v:.8g}' if isinstance(v,(float,int)) else str(v)) for v in r)+' |')
report.write_text(report.read_text()+'\n'.join(lines)+'\n')
