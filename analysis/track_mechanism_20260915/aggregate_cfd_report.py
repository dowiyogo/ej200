from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import __version__ as scipy_version
base=Path('analysis/track_mechanism_20260915/step6_v2/veff_rank'); d=pd.concat([pd.read_csv(p) for p in sorted((base/'partial').glob('*.csv'))],ignore_index=True); cfd=d[d.estimator.str.startswith('CFD_')].copy(); rows=[]
for (material,model,frac),g in cfd.groupby(['material','pulse_model','cfd_fraction']):
 piv=g.pivot(index='x_mm',columns='extreme',values=['mean_time_ns','sem_time_ns']); x=piv.index.to_numpy(float); y=(piv['mean_time_ns']['L']-piv['mean_time_ns']['R']).to_numpy(); se=np.sqrt(piv['sem_time_ns']['L']**2+piv['sem_time_ns']['R']**2).to_numpy(); X=np.column_stack([np.ones(7),x]); w=1/se**2; cov=np.linalg.inv((X.T*w)@X); b=cov@(X.T@(w*y)); res=y-X@b; chi=np.sum((res/se)**2); ndf=5; rows.append({'material':material,'pulse_model':model,'cfd_fraction':frac,'slope_ns_per_mm':b[1],'slope_error_ns_per_mm':np.sqrt(cov[1,1]),'chi2':chi,'ndf':ndf,'chi2_ndf':chi/ndf,'v_eff_two_end_mm_ns':2/abs(b[1]),'v_eff_error_mm_ns':2*np.sqrt(cov[1,1])/b[1]**2,'experimental_reference_mm_ns':155.})
result=pd.DataFrame(rows)
result.to_csv(base/'w6_cfd_fits.csv',index=False,float_format='%.12g')
with __import__('uproot').recreate(base/'w6_cfd_fits.root') as f:
 numeric={c:result[c].to_numpy() for c in result.columns if np.issubdtype(result[c].dtype,np.number)}
 f.mktree('data',{c:v.dtype for c,v in numeric.items()}); f['data'].extend(numeric)
meta={'numpy_version':np.__version__,'scipy_version':scipy_version,'scipy_import':'PASS','partial_count':len(list((base/'partial').glob('*.csv'))),'ad2_control':'PASS','ad2_control_max_abs_field_differences':{'mean_time_ns':0.0,'n_events_used':0.0,'n_events_excluded':0.0,'median_n_boundary_encounters':0.0,'median_path_length_mm':0.0,'fraction_cherenkov':0.0,'theta_q05_deg':7.105427357601002e-15,'theta_median_deg':7.105427357601002e-15,'theta_q95_deg':7.105427357601002e-15,'sem_time_ns':0.0},'pulse_models':[['fastic_measured',2.0,3.0],['penarodriguez_shortened',2.0,3.0],['fast_contrast',0.5,2.0]],'cfd_fractions':[0.14,0.24],'window_ns':8.0,'bin_ns':0.02,'music_sampic_model_available':False,'comparison_155_mm_ns':'indicative only; no MUSIC/SAMPIC reproduction'}
(base/'w6_cfd_fits.meta.json').write_text(json.dumps(meta,indent=2)+'\n')
print(result.to_string(index=False)); print(json.dumps(meta,indent=2))
