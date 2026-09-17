from pathlib import Path
import json, numpy as np, pandas as pd
base=Path('analysis/track_mechanism_20260915/step6_v2/veff_rank')
parts=pd.concat([pd.read_csv(p) for p in sorted((base/'partial').glob('*.csv'))],ignore_index=True)
c=parts[parts.estimator.str.startswith('CFD_')].copy()
rows=[]
for keys,g in c.groupby(['material','pulse_model','cfd_fraction']):
 for extreme,gx in g.groupby('extreme'):
  w=gx.n_events_used.to_numpy(float); rows.append(dict(material=keys[0],pulse_model=keys[1],cfd_fraction=keys[2],extreme=extreme,theta_q05_weighted_mean_deg=np.average(gx.theta_q05_deg,weights=w),theta_median_weighted_mean_deg=np.average(gx.theta_median_deg,weights=w),theta_q95_weighted_mean_deg=np.average(gx.theta_q95_deg,weights=w),reference_experimental_theta_deg=35.4))
ae=pd.DataFrame(rows); ae.to_csv(base/'ae1_theta_cfd_by_extreme.csv',index=False)
pos=[]
for keys,g in c.groupby(['material','pulse_model','cfd_fraction','x_mm']):
 w=g.n_events_used.to_numpy(float); pos.append(dict(material=keys[0],pulse_model=keys[1],cfd_fraction=keys[2],x_mm=keys[3],theta_q05_weighted_mean_deg=np.average(g.theta_q05_deg,weights=w),theta_median_weighted_mean_deg=np.average(g.theta_median_deg,weights=w),theta_q95_weighted_mean_deg=np.average(g.theta_q95_deg,weights=w),reference_experimental_theta_deg=35.4))
pd.DataFrame(pos).to_csv(base/'ae1_theta_cfd_by_position.csv',index=False)
r=[]
for keys,g in c.groupby(['material','pulse_model','cfd_fraction']):
 for extreme,gx in g.groupby('extreme'):
  ratio=np.average(gx.mean_path_over_daxial,weights=gx.n_events_used); vgroup=189.74206202531644 if keys[0]=='EJ-230' else (166.198 if keys[0]=='EJ-200' else 171.893); r.append(dict(material=keys[0],pulse_model=keys[1],cfd_fraction=keys[2],extreme=extreme,path_over_daxial_exact=ratio,v_group_MPT_mm_ns=vgroup,v_implied_mm_ns=vgroup/ratio,reference_experimental_theta_deg=35.4))
ae2=pd.DataFrame(r); ae2.to_csv(base/'ae2_geometry_cfd.csv',index=False)
w6=pd.read_csv(base/'w6_cfd_fits.csv'); w6['difference_from_experimental_reference_mm_ns']=w6.v_eff_two_end_mm_ns-w6.experimental_reference_mm_ns; w6.to_csv(base/'w6_cfd_fits.csv',index=False,float_format='%.12g')
meta=json.loads((base/'w6_cfd_fits.meta.json').read_text()); meta.update({'ae1_theta_method':'weighted mean of per-cell theta quantiles from partials','ae2_geometry_method':'event-level path/d_axial persisted in CFD partials; weighted mean by events used','ae2_exact_ratio':'AVAILABLE','music_sampic_model_available':False}); (base/'w6_cfd_fits.meta.json').write_text(json.dumps(meta,indent=2)+'\n')
print('AE1',len(ae),'AE2',len(ae2),'W6',len(w6))
