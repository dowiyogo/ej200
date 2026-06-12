#!/usr/bin/env python3
"""EXEC_12T fourth-to-thirtieth detected-hit timing analysis."""
from __future__ import annotations
import argparse,json,math,pathlib,re,subprocess
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np,pandas as pd,uproot
from scipy.stats import skew,kurtosis
import importlib.util

HOOK_FINE_SCAN_DIR=pathlib.Path("/home/reriosto/SHiP/t0minidaq/pairscan_2026-06-11")
HOOK_EXEC11_RESULTS="AUTO";HOOK_EXEC12_RESULTS="AUTO"
HOOK_WINDOW_DIP_DIR=pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec08b_window_dip")
HOOK_PAIR_IDS=(28,29);HOOK_N_EVENTS=3000;HOOK_PRIMARY_THRESHOLDS=(4,20);HOOK_THRESHOLD_SWEEP=tuple(range(1,31))
HOOK_FIT_RANGE_NSIGMA=2.5;HOOK_BOOTSTRAP_REPLICAS=2000;HOOK_RANDOM_SEED=20260612;HOOK_COMMON_SAMPLE_THRESHOLD=20
HOOK_LANGUAGE="EN";HOOK_DATA_COMMIT="f431c01";EXPECTED=np.arange(-462.,-421.)

P11=pathlib.Path(__file__).with_name("exec11_pair_analysis.py");S=importlib.util.spec_from_file_location("e11",P11);E11=importlib.util.module_from_spec(S);S.loader.exec_module(E11)
def gitsha():return subprocess.check_output(["git","rev-parse","--short","HEAD"],text=True).strip()
def xpos(p):return float(re.search(r"x([+-]?\d+(?:\.\d+)?)mm",p.name).group(1))
def latest(prefix):
 d=sorted(pathlib.Path("results").glob(prefix),key=lambda p:p.stat().st_mtime);return d[-1]
def order_statistics(event,gid,time,n=3000,kmax=30):
 out=[];counts=[]
 for ch in HOOK_PAIR_IDS:
  sel=(gid==ch)&np.isfinite(time);ev=event[sel];tm=time[sel];c=np.bincount(ev,minlength=n);a=np.full((n,kmax),np.nan)
  order=np.lexsort((tm,ev));se=ev[order];st=tm[order];starts=np.r_[0,np.cumsum(c[:-1])]
  for k in range(1,kmax+1):
   ok=c>=k;a[ok,k-1]=st[starts[ok]+k-1]
  out.append(a);counts.append(c)
 return counts[0],counts[1],out[0],out[1]
def robust(v):
 v=np.asarray(v);v=v[np.isfinite(v)]
 if not len(v):return dict(n=0,mean=np.nan,median=np.nan,sigma_core=np.nan,rms=np.nan,rms68=np.nan,mad_sigma=np.nan,skewness=np.nan,excess_kurtosis=np.nan,tail2=np.nan,tail3=np.nan)
 med=np.median(v);q16,q84=np.quantile(v,[.16,.84]);core=(q84-q16)/2
 return dict(n=len(v),mean=np.mean(v),median=med,sigma_core=core,rms=np.std(v,ddof=1),rms68=core,mad_sigma=1.4826*np.median(np.abs(v-med)),skewness=skew(v),excess_kurtosis=kurtosis(v),tail2=np.mean(np.abs(v-med)>2*core),tail3=np.mean(np.abs(v-med)>3*core))
def fit_line(x,y,e=None):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
 if e is not None:e=np.asarray(e,dtype=float)
 X=np.c_[x,np.ones(len(x))];w=np.ones(len(x)) if e is None else 1/np.square(e);cov=np.linalg.inv(X.T@(w[:,None]*X));p=cov@X.T@(w*y);res=y-X@p;return p,cov,float(np.sum(w*res*res)/(len(x)-2))
def cache_file(path,out,n_events=HOOK_N_EVENTS):
 es=[];gs=[];ts=[];guns=[];entries=bad=0;present=np.zeros(n_events,bool)
 with uproot.open(path) as f:
  tree=f["sipm_hits"];schema=list(tree.keys())
  for a in tree.iterate(["event_id","global_id","time_ns","gun_x_mm"],step_size="100 MB",library="np"):
   entries+=len(a["event_id"]);bad+=np.count_nonzero(~np.isfinite(a["time_ns"]));present[np.unique(a["event_id"])]=1;guns.extend(np.unique(a["gun_x_mm"]).tolist())
   q=np.isin(a["global_id"],HOOK_PAIR_IDS);es.append(a["event_id"][q]);gs.append(a["global_id"][q]);ts.append(a["time_ns"][q])
 ev=np.concatenate(es);gid=np.concatenate(gs);tm=np.concatenate(ts);na,nb,ta,tb=order_statistics(ev,gid,tm,n_events)
 x=xpos(path);np.savez_compressed(out/f"pair_x{x:+.1f}mm.npz",x_true_mm=x,event_id=np.arange(n_events),npe_A=na,npe_B=nb,t_A=ta,t_B=tb)
 return dict(file=path.name,x_true_mm=x,n_entries=entries,n_event_ids_present=present.sum(),event_id_min=np.flatnonzero(present)[0],event_id_max=np.flatnonzero(present)[-1],gun_x_values=";".join(map(str,sorted(set(guns)))),hits_A=na.sum(),hits_B=nb.sum(),nonfinite_time_fraction=bad/entries,schema=",".join(schema),size_bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
def load(p):
 with np.load(p) as z:return {k:z[k] for k in z.files}
def inventory_cache(out):
 files=sorted(HOOK_FINE_SCAN_DIR.glob("pairscan_x*.root"),key=xpos)
 if not np.array_equal([xpos(p) for p in files],EXPECTED):raise RuntimeError("fine-scan positions incomplete")
 rows=[];(out/"derived/order_statistics").mkdir(parents=True,exist_ok=True)
 for i,p in enumerate(files):print(f"{i+1}/41 {p.name}",flush=True);rows.append(cache_file(p,out/"derived/order_statistics"))
 pd.DataFrame(rows).to_csv(out/"analysis/data_inventory.csv",index=False)
 prov={"fine_scan":{"material":"OPSC-101/EJ-204 confirmed from pairscan.log","readout":"Top confirmed from pairscan.log","jitter":"0 ns confirmed from macros","data_commit":HOOK_DATA_COMMIT},"conceptual_distinction":"20th detected-hit timing != 20% constant-fraction timing from Lv et al.","analysis_commit":gitsha()}
 (out/"analysis/configuration_provenance.json").write_text(json.dumps(prov,indent=2)+"\n")
def events(d,k,matched=False):
 valid=(d["npe_A"]>=k)&(d["npe_B"]>=k)
 if matched:valid=(d["npe_A"]>=20)&(d["npe_B"]>=20)
 a=d["t_A"][:,k-1];b=d["t_B"][:,k-1]
 return valid,a,b,1000*(a-b),500*(a+b)
def reproduce(out):
 e11=latest("exec11_*");old=pd.read_csv(e11/"analysis/pairscan_summary_v2.csv");rows=[]
 for p in sorted((out/"derived/order_statistics").glob("*.npz"),key=lambda p:load(p)["x_true_mm"]):
  d=load(p);v,a,b,dt,tp=events(d,4);f=E11.iterative_gaussian_fit(dt[v],np.linspace(-3000,3000,301),guard_subpeak=True);o=old[old.x_true_mm==d["x_true_mm"]].iloc[0]
  rows.append(dict(x_true_mm=d["x_true_mm"],npe_A_exact=d["npe_A"].mean()==o.mean_npe_A,npe_B_exact=d["npe_B"].mean()==o.mean_npe_B,eff_exact=v.mean()==o.efficiency,mu_new=f["mean"],mu_old=o.mu_dt_ps,sigma_new=f["sigma"],sigma_old=o.sigma_dt_ps))
 r=pd.DataFrame(rows);r.to_csv(out/"analysis/exec11_reproduction_check.csv",index=False)
 p,c,ch=fit_line(r.x_true_mm.to_numpy(),r.mu_new.to_numpy());po,co,cho=fit_line(old.x_true_mm.to_numpy(),old.mu_dt_ps.to_numpy())
 ok=bool(np.allclose(r.mu_new,r.mu_old) and np.allclose(r.sigma_new,r.sigma_old) and abs(p[0]-po[0])<2*math.sqrt(c[0,0]+co[0,0]))
 (out/"logs/exec11_reproduction.log").write_text(f"pass={ok}\nslope_new={p[0]}\nslope_old={po[0]}\n")
 if not ok:raise RuntimeError("EXEC_11 4PE reproduction gate failed")
def analyze(out):
 data=[load(p) for p in sorted((out/"derived/order_statistics").glob("*.npz"),key=lambda p:load(p)["x_true_mm"])];summary=[];sweep=[];preds={4:[],20:[]}
 for k in HOOK_THRESHOLD_SWEEP:
  pos=[]
  for d in data:
   v,a,b,dt,tp=events(d,k);rd=robust(dt[v]);rt=robust(tp[v]);pos.append(dict(x_true_mm=d["x_true_mm"],threshold=k,efficiency=v.mean(),mean_dt=rd["mean"],sigma_dt_core=rd["sigma_core"],rms68_dt=rd["rms68"],tail3_dt=rd["tail3"],mean_tplus=rt["mean"],sigma_tplus_core=rt["sigma_core"],rho_ab=np.corrcoef(a[v],b[v])[0,1],**{f"dt_{q}":z for q,z in rd.items() if q in ("rms","skewness","excess_kurtosis")}))
  pf=pd.DataFrame(pos);coef,cov,chi=fit_line(pf.x_true_mm.to_numpy(),pf.mean_dt.to_numpy())
  cv=[]
  for i,d in enumerate(data):
   tr=pf.drop(i);co,_,_=fit_line(tr.x_true_mm.to_numpy(),tr.mean_dt.to_numpy());v,a,b,dt,tp=events(d,k);xr=(dt[v]-co[1])/co[0];m=robust(xr-d["x_true_mm"]);cv.append((m,xr,v))
   if k in preds:preds[k].append((d["x_true_mm"],xr,v,tp[v]))
  sweep.append(dict(threshold=k,efficiency_mean=pf.efficiency.mean(),efficiency_min=pf.efficiency.min(),slope_ps_per_mm=coef[0],slope_error=math.sqrt(cov[0,0]),calibration_chi2_ndf=chi,mean_sigma_delta_t_core=pf.sigma_dt_core.mean(),mean_rms68_delta_t=pf.rms68_dt.mean(),mean_sigma_x_core_cv=np.mean([q[0]["sigma_core"] for q in cv]),mean_rms68_x_cv=np.mean([q[0]["rms68"] for q in cv]),max_abs_bias=max(abs(q[0]["mean"]) for q in cv),tail_fraction=np.mean([q[0]["tail3"] for q in cv])))
  if k in (4,20):summary+=pos
 pd.DataFrame(sweep).to_csv(out/"analysis/threshold_sweep_summary.csv",index=False);pd.DataFrame(summary).to_csv(out/"analysis/temporal_position_summary.csv",index=False)
 for k in (4,20):
  rows=[];allx=[];alltrue=[]
  for x,xr,v,tp in preds[k]:rows.append(dict(x_true_mm=x,threshold=k,**robust(xr-x),valid_fraction=v.mean()));allx.append(xr);alltrue.append(np.full(len(xr),x))
  pd.DataFrame(rows).to_csv(out/f"analysis/position_reconstruction_{k}pe.csv",index=False);np.savez_compressed(out/f"analysis/loo_predictions_{k}pe.npz",x_rec=np.concatenate(allx),x_true=np.concatenate(alltrue))
 # matched comparison and calibrations
 comps=[]
 for k in (4,20):
  pf=pd.DataFrame(summary);pf=pf[pf.threshold==k];co,cov,chi=fit_line(pf.x_true_mm.to_numpy(),pf.mean_dt.to_numpy())
  pd.DataFrame([dict(threshold=k,slope_ps_per_mm=co[0],slope_error=math.sqrt(cov[0,0]),intercept_ps=co[1],chi2_ndf=chi,v_eff_cm_ns=2000/abs(co[0])/10)]).to_csv(out/f"analysis/calibration_{k}pe.csv",index=False)
  for matched in (False,True):
   vals=[]
   for d in data:
    v,a,b,dt,tp=events(d,k,matched);vals.append((v.mean(),robust(dt[v])["sigma_core"],robust(tp[v])["sigma_core"]))
   comps.append(dict(threshold=k,sample="matched20" if matched else "native",efficiency_mean=np.mean([q[0] for q in vals]),sigma_dt_core_mean=np.mean([q[1] for q in vals]),sigma_tplus_core_mean=np.mean([q[2] for q in vals])))
 pd.DataFrame(comps).to_csv(out/"analysis/threshold_4_20_summary.csv",index=False)
 make_figures(out,data)
def make_figures(out,data):
 sw=pd.read_csv(out/"analysis/threshold_sweep_summary.csv");sm=pd.read_csv(out/"analysis/temporal_position_summary.csv")
 for col,name,y in [("efficiency_mean","threshold_efficiency","Efficiency"),("mean_sigma_delta_t_core","threshold_sigma_dt","Mean delta-t RMS68 (ps)"),("mean_sigma_x_core_cv","threshold_sigma_x","Mean X RMS68 (mm)"),("max_abs_bias","threshold_bias","Max absolute bias (mm)"),("slope_ps_per_mm","threshold_slope","Slope (ps/mm)"),("calibration_chi2_ndf","threshold_chi2","Calibration chi2/ndf")]:
  fig,ax=plt.subplots();ax.plot(sw.threshold.to_numpy(),sw[col].to_numpy(),"o-");ax.scatter([4,20],sw.set_index("threshold").loc[[4,20],col].to_numpy(),c=["blue","red"]);ax.set(xlabel="Detected-hit order k",ylabel=y,title="Pair (28,29), EJ-204, intrinsic");ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/f"figures/{name}.pdf");plt.close(fig)
 fig,ax=plt.subplots();ax.plot(sw.efficiency_mean.to_numpy(),sw.mean_sigma_x_core_cv.to_numpy(),"o-");ax.scatter(sw.loc[sw.threshold.isin([4,20]),"efficiency_mean"].to_numpy(),sw.loc[sw.threshold.isin([4,20]),"mean_sigma_x_core_cv"].to_numpy(),c=["blue","red"]);ax.set(xlabel="Mean efficiency",ylabel="Mean X RMS68 (mm)",title="Efficiency-resolution trade-off");fig.tight_layout();fig.savefig(out/"figures/threshold_pareto.pdf");plt.close(fig)
 for col,name,y in [("efficiency","efficiency_4_20_vs_x","Efficiency"),("sigma_dt_core","timing_width_4_20_vs_x","delta-t RMS68 (ps)"),("mean_dt","mean_delta_4_20_vs_x","Mean delta-t (ps)"),("rho_ab","correlation_ab_vs_x","rho(A,B)")]:
  fig,ax=plt.subplots();[ax.plot(g.x_true_mm.to_numpy(),g[col].to_numpy(),"o-",label=f"{k}th hit") for k,g in sm.groupby("threshold")];ax.set(xlabel="X (mm)",ylabel=y);ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/f"figures/{name}.pdf");plt.close(fig)
 refs=pd.read_csv(latest("exec11_*")/"analysis/reference_positions.csv")
 for _,r in refs.iterrows():
  d=min(data,key=lambda q:abs(q["x_true_mm"]-r.x_true_mm));fig,ax=plt.subplots()
  for k in (4,20):
   v,a,b,dt,tp=events(d,k);ax.hist(dt[v],bins=80,histtype="step",density=True,label=f"{k}th detected hit")
  ax.set(xlabel="delta-t (ps)",ylabel="density",title=f"{r.reference}: x={r.x_true_mm:.0f} mm");ax.legend();fig.tight_layout();fig.savefig(out/f"figures/{r.reference.lower()}_delta_t_4_20.pdf");plt.close(fig)
 # conceptual diagrams
 fig,ax=plt.subplots(figsize=(8,3));hits=np.array([.2,.35,.5,.7,.9,1.2,1.5,1.8,2.1,2.4,2.8,3.1,3.4,3.8,4.1,4.5,4.9,5.3,5.7,6.1]);ax.eventplot(hits);ax.axvline(hits[3],c="blue",label="4th detected hit");ax.axvline(hits[19],c="red",label="20th detected hit");ax.set_yticks([]);ax.set_xlabel("arrival time");ax.legend();fig.tight_layout();fig.savefig(out/"figures/order_statistic_schematic.pdf");plt.close(fig)
def window(out):
 rows=[]
 for p in sorted(HOOK_WINDOW_DIP_DIR.glob("*.root"),key=xpos):
  x=xpos(p);n=2000;es=[];gs=[];ts=[]
  with uproot.open(p) as f:
   for a in f["sipm_hits"].iterate(["event_id","global_id","time_ns"],step_size="100 MB",library="np"):
    es.append(a["event_id"]);gs.append(a["global_id"]);ts.append(a["time_ns"])
  ev=np.concatenate(es);gid=np.concatenate(gs);tm=np.concatenate(ts);local=int(np.argmin(abs(E11.np.r_[(-692+20*np.arange(35)),(12+20*np.arange(35))]-x)))+16
  na,nb,ta,tb=order_statistics(ev,np.where(gid==local,28,29),tm,n) # A=nearest, B=dummy
  for k in (4,20):v=na>=k;rows.append(dict(run=p.stem,x_true_mm=x,nearest_top_id=local,distance_to_center_mm=abs(((-692+20*(local-16)) if local-16<35 else (12+20*(local-51)))-x),threshold=k,mean_npe=na.mean(),efficiency=v.mean(),mean_time_ns=np.nanmean(ta[:,k-1]),sigma_time_ns=robust(ta[v,k-1])["sigma_core"]))
 pd.DataFrame(rows).to_csv(out/"analysis/window_dip_summary.csv",index=False)
def main():
 p=argparse.ArgumentParser();p.add_argument("stage",choices=["cache","reproduce","analyze","window"]);p.add_argument("--output-dir",type=pathlib.Path,required=True);a=p.parse_args()
 for d in ["analysis","derived/order_statistics","figures","tables","report","beamer","logs"]:(a.output_dir/d).mkdir(parents=True,exist_ok=True)
 {"cache":inventory_cache,"reproduce":reproduce,"analyze":analyze,"window":window}[a.stage](a.output_dir)
if __name__=="__main__":main()
