#!/usr/bin/env python3
# CUARENTENA 2026-09-25 — Ver AUDITORIA_REPO_20260925.md §8: Micro-benchmark usado para medir tiempos de CPU de distintas vectorizaciones numpy/awkward sobre una celda ROOT antes de correr analyze_veff_rank_cfd.py (OBSOLETO / SCRATCH, nombrado explícitamente en la lista de cuarentena).
import time
from pathlib import Path
import numpy as np
import uproot

ROOT=next(Path('/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ200_xp0').glob('attempts/*/photon_hits_run000.root'))
BIN_NS=.02; RISE=.5; FALL=2.; WINDOW=40.; N_EVENTS=100

def cfd(times, frac):
 t0=times.min(); nb=int(WINDOW/BIN_NS); bins=np.floor((times-t0)/BIN_NS).astype(int); keep=(bins>=0)&(bins<nb); h=np.bincount(bins[keep],minlength=nb).astype(float); ar=np.exp(-BIN_NS/RISE); af=np.exp(-BIN_NS/FALL); r=np.zeros(nb); f=np.zeros(nb)
 for i in range(1,nb): r[i]=ar*r[i-1]+h[i]; f[i]=af*f[i-1]+h[i]
 p=f-r; hit=np.flatnonzero(p>=frac*p.max()); return t0+hit[0]*BIN_NS if len(hit) else np.nan

t0=time.perf_counter();
with uproot.open(ROOT) as rf: arr=rf['sipm_hits'].arrays(['event_id','face_type','time_ns'],library='np')
t1=time.perf_counter();
mask=arr['event_id']<N_EVENTS; grouped={}
for face in (0,1):
 for event in range(N_EVENTS): grouped[(event,face)]=arr['time_ns'][mask&(arr['face_type']==face)&(arr['event_id']==event)]
t2=time.perf_counter();
cfd_values=[]
for face in (0,1):
 for event in range(N_EVENTS): cfd_values.append(cfd(grouped[(event,face)],.14))
t3=time.perf_counter();
print('root_read_s',t1-t0,'per_event_filter_group_s',t2-t1,'cfd_python_s',t3-t2,'total_s',t3-t0)
print('hits_per_event_face_mean',np.mean([len(v) for v in grouped.values()]))
for q in (.01,.02,.05,.10,.14):
 diffs=[]
 for face in (0,1):
  for event in range(N_EVENTS):
   times=grouped[(event,face)]; qtime=np.quantile(times,q); explicit=cfd_values[face*N_EVENTS+event]; diffs.append((qtime-explicit)*1000)
 diffs=np.asarray(diffs); print('quantile',q,'mean_diff_ps',diffs.mean(),'sd_diff_ps',diffs.std(ddof=1),'mae_ps',np.mean(np.abs(diffs)),'p95_abs_ps',np.quantile(abs(diffs),.95))
