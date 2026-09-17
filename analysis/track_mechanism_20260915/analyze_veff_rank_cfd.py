#!/usr/bin/env python3
"""Checkpointed rank timestamps and digital CFD for EXEC46 v2."""
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np, pandas as pd, uproot
N_VALUES=(1,2,5,20,100); MATERIALS={0:'EJ-200',1:'EJ-204',2:'EJ-230'}; POSITIONS=(-650,-500,-200,0,200,500,650)
PULSE_MODELS=(('fastic_measured',2.,3.),('penarodriguez_shortened',2.,3.),('fast_contrast',.5,2.))
CFD_FRACTIONS=(.14,.24); BIN_NS=.02; CFD_WINDOW_NS=8.
def root_path(c): return next(Path(c).glob('attempts/*/photon_hits_run000.root'))
def pulse(times,rise,fall):
 t0=float(times.min()); nb=int(CFD_WINDOW_NS/BIN_NS); b=np.floor((times-t0)/BIN_NS).astype(int); keep=(b>=0)&(b<nb); h=np.bincount(b[keep],minlength=nb).astype(float); ar=math.exp(-BIN_NS/rise); af=math.exp(-BIN_NS/fall); r=np.zeros(nb); f=np.zeros(nb)
 for i in range(1,nb): r[i]=ar*r[i-1]+h[i]; f[i]=af*f[i-1]+h[i]
 p=f-r; amp=p.max(); out=[]
 for frac in CFD_FRACTIONS:
  hit=np.flatnonzero(p>=frac*amp); out.append(t0+hit[0]*BIN_NS if len(hit) else np.nan)
 return out
def process(cell,material,x,outdir):
 target=outdir/'partial'/f'{Path(cell).name}.csv'; target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists() and len(pd.read_csv(target)): return 'SKIP'
 with uproot.open(root_path(cell)) as f: a=f['sipm_hits'].arrays(['event_id','face_type','time_ns','path_length_mm','n_boundary_encounters','source_type','x_mm','x_creation_mm','y_mm','y_creation_mm','z_mm','z_creation_mm'],library='np')
 rows=[]
 for face,tag in ((0,'L'),(1,'R')):
  m=a['face_type']==face; indices=np.flatnonzero(m); e=a['event_id'][m]; order=np.lexsort((a['time_ns'][m],e)); es=e[order]; starts=np.r_[0,np.flatnonzero(es[1:]!=es[:-1])+1]; stops=np.r_[starts[1:],len(es)]; rank=np.arange(len(es))-np.repeat(starts,np.diff(np.r_[starts,len(es)])); keep=rank<200; ev=es[keep]; rk=rank[keep]; idx=indices[order[keep]]; grouped=[a['time_ns'][indices[order[lo:hi]]] for lo,hi in zip(starts,stops)]; vals={k:np.full((10000,200),np.nan) for k in ('time_ns','path_length_mm','n_boundary_encounters','source_type','theta_eff_deg')}
  for k in ('time_ns','path_length_mm','n_boundary_encounters','source_type'): vals[k][ev,rk]=a[k][idx]
  axial=np.sqrt((a['x_mm'][idx]-a['x_creation_mm'][idx])**2+(a['y_mm'][idx]-a['y_creation_mm'][idx])**2+(a['z_mm'][idx]-a['z_creation_mm'][idx])**2); vals['theta_eff_deg'][ev,rk]=np.degrees(np.arccos(np.clip(axial/vals['path_length_mm'][ev,rk],-1,1)))
  for n in N_VALUES:
   q=np.isfinite(vals['time_ns'][:,n-1]); rows.append(dict(material=material,x_mm=x,extreme=tag,estimator=f'N{n}',mean_time_ns=float(vals['time_ns'][q,n-1].mean()),n_events_used=int(q.sum()),n_events_excluded=int((~q).sum()),median_n_boundary_encounters=float(np.median(vals['n_boundary_encounters'][q,n-1])),median_path_length_mm=float(np.median(vals['path_length_mm'][q,n-1])),fraction_cherenkov=float(np.mean(vals['source_type'][q,n-1]==2)),theta_q05_deg=float(np.quantile(vals['theta_eff_deg'][q,n-1],.05)),theta_median_deg=float(np.median(vals['theta_eff_deg'][q,n-1])),theta_q95_deg=float(np.quantile(vals['theta_eff_deg'][q,n-1],.95))))
  for model,rise,fall in PULSE_MODELS:
   for frac in CFD_FRACTIONS:
    ts=[]; paths=[]; bounds=[]; sources=[]; thetas=[]; path_over_axial=[]
    for event,times in enumerate(grouped):
     stamp=pulse(times,rise,fall)[list(CFD_FRACTIONS).index(frac)]; sel=indices[order[starts[event]:stops[event]]]; hit=sel[np.argmin(np.abs(times-stamp))] if np.isfinite(stamp) else None
     if hit is None: continue
     ts.append(stamp); paths.append(a['path_length_mm'][hit]); bounds.append(a['n_boundary_encounters'][hit]); sources.append(a['source_type'][hit]); axial=abs((700.0 if face==1 else -700.0)-a['x_creation_mm'][hit]); path_over_axial.append(a['path_length_mm'][hit]/axial); thetas.append(math.degrees(math.acos(np.clip(axial/a['path_length_mm'][hit],-1,1))))
    rows.append(dict(material=material,x_mm=x,extreme=tag,estimator=f'CFD_{model}_{int(frac*100)}',pulse_model=model,tau_rise_ns=rise,tau_fall_ns=fall,cfd_fraction=frac,mean_time_ns=float(np.mean(ts)),sem_time_ns=float(np.std(ts,ddof=1)/math.sqrt(len(ts))),n_events_used=len(ts),n_events_excluded=10000-len(ts),median_n_boundary_encounters=float(np.median(bounds)),median_path_length_mm=float(np.median(paths)),mean_path_over_daxial=float(np.mean(path_over_axial)),fraction_cherenkov=float(np.mean(np.array(sources)==2)),theta_q05_deg=float(np.quantile(thetas,.05)),theta_median_deg=float(np.median(thetas)),theta_q95_deg=float(np.quantile(thetas,.95))))
 pd.DataFrame(rows).to_csv(target,index=False); return 'DONE'
def main():
 p=argparse.ArgumentParser(); p.add_argument('--output-dir',type=Path,required=True); args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True); start=time.time(); done=[]; pending=[]; base=Path('/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells')
 for code,material in MATERIALS.items():
  for x in POSITIONS:
   cell=base / (f'{material.replace("-","")}_'+('xm' if x<0 else 'xp')+str(abs(x))); status=process(cell,material,x,args.output_dir); (done if status in ('DONE','SKIP') else pending).append(cell.name); print(cell.name,status,'elapsed',round(time.time()-start,1),flush=True)
 print(json.dumps({'completed':done,'pending':pending,'elapsed_s':time.time()-start,'resume_command':f'python3 {Path(__file__).resolve()} --output-dir {args.output_dir}'},indent=2))
if __name__=='__main__': main()
