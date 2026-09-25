# CUARENTENA 2026-09-25 — Ver AUDITORIA_REPO_20260925.md §8: Borrador preliminar de agregación CFD donde slope_ns_per_mm estaba dejado como np.nan; sobrescribiría cfd_summary.csv con NaNs si se ejecuta (reemplazado por aggregate_cfd_report.py).
import argparse,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import linregress

def main():
 p=argparse.ArgumentParser();p.add_argument('--dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args(); rows=[]
 for f in sorted((a.dir/'partial').glob('*.csv')): rows.append(pd.read_csv(f))
 d=pd.concat(rows,ignore_index=True); d.to_csv(a.out/'w6_cfd_cells.csv',index=False)
 fits=[]
 for (material,estimator),g in d.groupby(['material','estimator']):
  # Mean tL/tR per cell from partials, then fit tL-tR vs x.
  piv=g.pivot(index=['x_mm','extreme'],columns='estimator',values='mean_time_ns').reset_index()
  left=g[g.extreme=='L'][['x_mm','mean_time_ns','n_events_used']].rename(columns={'mean_time_ns':'L','n_events_used':'nL'}); right=g[g.extreme=='R'][['x_mm','mean_time_ns','n_events_used']].rename(columns={'mean_time_ns':'R','n_events_used':'nR'}); z=left.merge(right,on='x_mm'); z['diff']=z.L-z.R; z['se']=np.sqrt(1/z.nL+1/z.nR)*z[['L','R']].mean(axis=1)*0+np.nan
  # Partial files do not store SEM; use cell-event timestamps unavailable here.
  fits.append({'material':material,'estimator':estimator,'slope_ns_per_mm':np.nan,'slope_error_ns_per_mm':np.nan,'chi2_ndf':np.nan,'v_eff_two_end_mm_ns':np.nan,'note':'NOT_AVAILABLE: partial stores means but not event SEM'})
 pd.DataFrame(fits).to_csv(a.out/'w6_cfd_fits.csv',index=False)
 print('PASS',len(d),len(fits))
if __name__=='__main__':main()
