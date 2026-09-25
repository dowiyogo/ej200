# CUARENTENA 2026-09-25 — Este script fue identificado en
# AUDITORIA_REPO_20260925.md como potencialmente dañino si se ejecuta
# sin revisión (puede sobrescribir resultados). NO ejecutar sin antes
# leer docs/execution_logs/ correspondiente. Movido y neutralizado por
# la reorganización de 2026-09-25.
# CUARENTENA 2026-09-25 — Ver AUDITORIA_REPO_20260925.md §8: Micro-benchmark usado para medir tiempos de CPU de distintas vectorizaciones numpy/awkward sobre una celda ROOT antes de correr analyze_veff_rank_cfd.py (OBSOLETO / SCRATCH, nombrado explícitamente en la lista de cuarentena).
import numpy as np, uproot
from pathlib import Path
ROOT=next(Path('/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ200_xp0').glob('attempts/*/photon_hits_run000.root'))
B=.02; R=2.; F=3.; frac=.14; N=100
def grouped():
 with uproot.open(ROOT) as f:a=f['sipm_hits'].arrays(['event_id','face_type','time_ns'],library='np')
 out=[]
 for face in (0,1):
  m=(a['face_type']==face)&(a['event_id']<N); idx=np.flatnonzero(m); o=np.lexsort((a['time_ns'][idx],a['event_id'][idx])); e=a['event_id'][idx][o]; st=np.r_[0,np.flatnonzero(e[1:]!=e[:-1])+1]; sp=np.r_[st[1:],len(e)]
  out += [a['time_ns'][idx[o[lo:hi]]] for lo,hi in zip(st,sp)]
 return out
def cfd(t,w):
 t0=t.min(); nb=int(w/B); b=np.floor((t-t0)/B).astype(int); h=np.bincount(b[(b>=0)&(b<nb)],minlength=nb).astype(float); ar=np.exp(-B/R); af=np.exp(-B/F); r=np.zeros(nb); f=np.zeros(nb)
 for i in range(1,nb): r[i]=ar*r[i-1]+h[i]; f[i]=af*f[i-1]+h[i]
 p=f-r; peak=np.argmax(p); hit=np.flatnonzero(p>=frac*p[peak]); return t0+(hit[0] if len(hit) else np.nan)*B,t0+peak*B
G=grouped(); long=np.array([cfd(t,40.)[0] for t in G]); peaks=np.array([cfd(t,40.)[1] for t in G]); print('peak_offset_median_ns',np.median(peaks-np.array([t.min() for t in G])),'peak_offset_q95_ns',np.quantile(peaks-np.array([t.min() for t in G]),.95))
for w in (8.,10.,12.,16.):
 short=np.array([cfd(t,w)[0] for t in G]); d=(short-long)*1000; print('window_ns',w,'mean_diff_ps',np.nanmean(d),'sd_diff_ps',np.nanstd(d,ddof=1),'max_abs_ps',np.nanmax(abs(d)),'q95_abs_ps',np.nanquantile(abs(d),.95))
