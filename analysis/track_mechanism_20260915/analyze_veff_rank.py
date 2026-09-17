#!/usr/bin/env python3
"""Rank-timestamp and digital-CFD effective velocity tables."""
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np, pandas as pd, uproot

N_VALUES=(1,2,3,5,10,20,50,100,200)
MATERIALS={0:'EJ-200',1:'EJ-204',2:'EJ-230'}
POSITIONS=(-650,-500,-200,0,200,500,650)
EXPECTED={'v2':'764c643e6ccfc3d0b96c7fc20d6753e86139233599ca19359c8d95e585844290','old':'816c302c3ed37f7eb29855bcad7e734a1b6a77238f60a30361635b68377e341d'}
# Explicit W6 assumption: a single-photoelectron double exponential pulse.
TAU_RISE_NS=0.5; TAU_FALL_NS=2.0; CFD_FRACTIONS=(0.14,0.24); BIN_NS=0.02; CFD_WINDOW_NS=40.0

def sha256(p):
 d=hashlib.sha256();
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): d.update(b)
 return d.hexdigest()

def fit(x,y,se):
 X=np.column_stack([np.ones(len(x)),x]); w=1/se**2; cov=np.linalg.inv((X.T*w)@X); b=cov@(X.T@(w*y)); r=y-X@b; chi=float(np.sum((r/se)**2)); return b,np.sqrt(np.diag(cov)),chi,len(x)-2

def root_path(cell): return next((cell).glob('attempts/*/photon_hits_run000.root'))
def load_cell(root):
 with uproot.open(root) as f:
  return f['sipm_hits'].arrays(['event_id','face_type','time_ns','path_length_mm','n_boundary_encounters','source_type','x_mm','x_creation_mm','y_mm','y_creation_mm','z_mm','z_creation_mm'],library='np')
def ranks(a):
 out=[]
 for face in (0,1):
  m=a['face_type']==face; eid=a['event_id'][m]; order=np.lexsort((a['time_ns'][m],eid)); e=eid[order]
  starts=np.r_[0, np.flatnonzero(e[1:] != e[:-1])+1]
  ranks_in_event=np.arange(len(e))-np.repeat(starts, np.diff(np.r_[starts,len(e)]))
  keep=ranks_in_event < 200; event=e[keep]; rank=ranks_in_event[keep]; idx=order[keep]
  vals={k:np.full((10000,200),np.nan) for k in ('time_ns','path_length_mm','n_boundary_encounters','source_type','d_axial_mm')}
  for k in ('time_ns','path_length_mm','n_boundary_encounters','source_type'):
   vals[k][event,rank]=a[k][m][idx]
  vals['d_axial_mm'][event,rank]=np.sqrt((a['x_mm'][m][idx]-a['x_creation_mm'][m][idx])**2+(a['y_mm'][m][idx]-a['y_creation_mm'][m][idx])**2+(a['z_mm'][m][idx]-a['z_creation_mm'][m][idx])**2)
  out.append(vals)
 return out

def process_campaign(base,label,expected,output):
 cells=[]; start=time.time()
 for code,material in MATERIALS.items():
  for x in POSITIONS:
   root=root_path(base / (f'{material.replace("-","")}_'+('xm' if x<0 else 'xp')+str(abs(x))))
   a=load_cell(root); rr=ranks(a); cells.append((label,material,x,rr))
   print(label,material,x,'elapsed',round(time.time()-start,1),flush=True)
 return cells

def main():
 p=argparse.ArgumentParser(); p.add_argument('--v2',type=Path,required=True); p.add_argument('--old',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
 # Hashes are checked before any ROOT read.
 for label,path in [('v2',args.v2),('old',args.old)]:
  if sha256(path)!=EXPECTED[label]: raise RuntimeError(label+' tree hash mismatch')
 # Production roots are read-only; one cell at a time.
 allcells=process_campaign(Path('/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells'),'v2',EXPECTED['v2'],args.output_dir)
 allcells+=process_campaign(Path('/home/rrios/exec46_20260915/full_grid/cells'),'old',EXPECTED['old'],args.output_dir)
 rows=[]; fits=[]; diagnostics=[]
 for label,material,x,rr in allcells:
  for n in N_VALUES:
   left=rr[0]['time_ns'][:,n-1]; right=rr[1]['time_ns'][:,n-1]; valid=np.isfinite(left)&np.isfinite(right); d=left[valid]-right[valid]; sem=d.std(ddof=1)/math.sqrt(len(d)); rows.append(dict(campaign=label,material=material,x_mm=x,N=n,mean_tdiff_ns=d.mean(),sem_tdiff_ns=sem,n_events_used=len(d),n_events_excluded=10000-len(d)))
  for face,tag in ((0,'L'),(1,'R')):
   v=rr[face]; z=v['time_ns'][:,n-1]; q=np.isfinite(z)
   diagnostics.append(dict(campaign=label, material=material, x_mm=x, N=n,
     extreme=tag, n_events_used=int(q.sum()), n_events_excluded=int((~q).sum()),
     median_n_boundary_encounters=float(np.median(v['n_boundary_encounters'][q,n-1])),
     median_path_length_mm=float(np.median(v['path_length_mm'][q,n-1])),
     fraction_cherenkov=float(np.mean(v['source_type'][q,n-1]==2)),
     median_theta_eff_deg=float(np.degrees(np.arccos(np.clip(
       v['d_axial_mm'][q,n-1]/v['path_length_mm'][q,n-1], -1, 1))))))
 c=pd.DataFrame(rows); f=pd.DataFrame(fits); dg=pd.DataFrame(diagnostics); c.to_csv(args.output_dir/'w1_tdiff_rank.csv',index=False); dg.to_csv(args.output_dir/'w4_rank_diagnostics.csv',index=False); print('PASS',len(c),len(dg))
if __name__=='__main__': main()
