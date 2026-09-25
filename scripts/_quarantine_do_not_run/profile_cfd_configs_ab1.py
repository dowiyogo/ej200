# CUARENTENA 2026-09-25 — Este script fue identificado en
# AUDITORIA_REPO_20260925.md como potencialmente dañino si se ejecuta
# sin revisión (puede sobrescribir resultados). NO ejecutar sin antes
# leer docs/execution_logs/ correspondiente. Movido y neutralizado por
# la reorganización de 2026-09-25.
# CUARENTENA 2026-09-25 — Ver AUDITORIA_REPO_20260925.md §8: Micro-benchmark usado para medir tiempos de CPU de distintas vectorizaciones numpy/awkward sobre una celda ROOT antes de correr analyze_veff_rank_cfd.py (OBSOLETO / SCRATCH, nombrado explícitamente en la lista de cuarentena).
import time
from pathlib import Path
import numpy as np
import uproot

ROOT=next(Path('/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ200_xp0').glob('attempts/*/photon_hits_run000.root'))
N_EVENTS=100; BIN_NS=.02; WINDOW=40.; FRACTION=.14

def run(rise,fall):
 t0=time.perf_counter()
 with uproot.open(ROOT) as f: a=f['sipm_hits'].arrays(['event_id','face_type','time_ns'],library='np')
 t1=time.perf_counter(); groups={}
 for face in (0,1):
  m=(a['face_type']==face)&(a['event_id']<N_EVENTS); idx=np.flatnonzero(m); order=np.lexsort((a['time_ns'][idx],a['event_id'][idx])); e=a['event_id'][idx][order]; starts=np.r_[0,np.flatnonzero(e[1:]!=e[:-1])+1]; stops=np.r_[starts[1:],len(e)]
  for lo,hi,event in zip(starts,stops,e[starts]): groups[(int(event),face)]=a['time_ns'][idx[order[lo:hi]]]
 t2=time.perf_counter(); results=[]; max_times=[]
 for face in (0,1):
  for event in range(N_EVENTS):
   times=groups[(event,face)]; start=times.min(); bins=np.floor((times-start)/BIN_NS).astype(int); nb=int(WINDOW/BIN_NS); keep=(bins>=0)&(bins<nb); h=np.bincount(bins[keep],minlength=nb).astype(float); ar=np.exp(-BIN_NS/rise); af=np.exp(-BIN_NS/fall); r=np.zeros(nb); f=np.zeros(nb)
   for i in range(1,nb): r[i]=ar*r[i-1]+h[i]; f[i]=af*f[i-1]+h[i]
   pulse=f-r; peak=int(np.argmax(pulse)); cross=np.flatnonzero(pulse>=FRACTION*pulse[peak]); results.append(start+cross[0]*BIN_NS if len(cross) else np.nan); max_times.append(start+peak*BIN_NS)
 t3=time.perf_counter(); return dict(rise=rise,fall=fall,window_ns=WINDOW,bin_ns=BIN_NS,bins=nb,read_s=t1-t0,group_s=t2-t1,cfd_s=t3-t2,total_s=t3-t0,cfd_ms_per_event_face=(t3-t2)/(2*N_EVENTS)*1000,peak_offset_median_ns=float(np.median(np.asarray(max_times)-np.array([groups[(e,fa)].min() for fa in (0,1) for e in range(N_EVENTS)]))))
for config in ((.5,2.),(2.,3.)):
 print(run(*config))
