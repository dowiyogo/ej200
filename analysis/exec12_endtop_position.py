#!/usr/bin/env python3
"""EXEC_12 global EndTop position analysis from per-hit ROOT trees."""
from __future__ import annotations
import argparse, json, math, pathlib, re, subprocess
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, uproot
from scipy.interpolate import PchipInterpolator

DATA=pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000")
EXPECTED=np.array([-690,-670,-650,-600,-550,-500,-450,-400,-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,600,650,670,690],float)
N_EVENTS=2000; N_CHANNELS=86; THRESH=4
END_Y=(np.arange(8)-3.5)*7.5
TOP_X=np.r_[(-692+20*np.arange(35)),(12+20*np.arange(35))]
CONTEXT="simulation prediction — EJ-204 — intrinsic optical timing"

def xpos(path): return float(re.search(r"x(-?\d+)mm",path.name).group(1))
def top_x(local): return float(TOP_X[local])
def end_y(local): return float(END_Y[local])
def pooled_fourth(event, time, n=N_EVENTS):
    out=np.full(n,np.nan); order=np.lexsort((time,event)); ev=event[order]; t=time[order]
    counts=np.bincount(ev,minlength=n); starts=np.r_[0,np.cumsum(counts[:-1])]; ok=counts>=4
    out[ok]=t[starts[ok]+3]; return out
def weighted_time(times,counts):
    valid=np.isfinite(times)&(counts>=4); den=np.sum(np.where(valid,counts,0),axis=1)
    return np.divide(np.nansum(np.where(valid,times*counts,np.nan),axis=1),den,out=np.full(len(den),np.nan),where=den>0)
def centroid(counts,pos):
    den=counts.sum(1); return np.divide(counts@pos,den,out=np.full(len(den),np.nan),where=den>0)
def blue(cov):
    if np.linalg.cond(cov)>1e10: raise np.linalg.LinAlgError("ill-conditioned covariance")
    one=np.ones(len(cov)); inv=np.linalg.inv(cov); return inv@one/(one@inv@one)
def fit_poly(raw,x,deg=1):
    return np.polyfit(raw,x,deg)
def metrics(values,true):
    v=values[np.isfinite(values)]
    if len(v)==0:return dict(mean_x_rec=np.nan,bias=np.nan,sigma_core=np.nan,rms=np.nan,rms68=np.nan,median=np.nan,mad_sigma=np.nan,fraction_outside_2sigma=np.nan,fraction_outside_3sigma=np.nan,valid_fraction=0.)
    r=v-true; med=np.median(v); mad=1.4826*np.median(np.abs(v-med))
    q16,q84=np.quantile(r,[.16,.84]); core=(q84-q16)/2
    return dict(mean_x_rec=np.mean(v),bias=np.mean(r),sigma_core=core,rms=np.sqrt(np.mean(r*r)),
      rms68=core,median=med,mad_sigma=mad,fraction_outside_2sigma=np.mean(np.abs(r-np.median(r))>2*core),
      fraction_outside_3sigma=np.mean(np.abs(r-np.median(r))>3*core),valid_fraction=len(v)/len(values))
def local_pair(counts,t4):
    n=len(counts); raw=np.full(n,np.nan); center=np.full(n,np.nan); status=np.full(n,"ok",object); ids=np.full((n,2),-1)
    for e in range(n):
        m=int(np.argmax(counts[e])); candidates=[j for j in (m-1,m+1) if 0<=j<70 and abs(TOP_X[j]-TOP_X[m])==20]
        if not candidates: status[e]="edge_or_gap"; continue
        excluded=[j for j in (m-1,m+1) if 0<=j<70 and abs(TOP_X[j]-TOP_X[m])!=20]
        if excluded and max(counts[e,j] for j in excluded)>max(counts[e,j] for j in candidates):
            status[e]="gap_ambiguous"; continue
        j=max(candidates,key=lambda q: counts[e,q]); a,b=sorted((m,j)); ids[e]=[a+16,b+16]; center[e]=(TOP_X[a]+TOP_X[b])/2
        if counts[e,a]>0 and counts[e,b]>0: raw[e]=math.log(counts[e,a]/counts[e,b])
        else: status[e]="zero_count"
    return raw,center,status,ids

def derive_file(path,out):
    x=xpos(path); counts=np.zeros((N_EVENTS,N_CHANNELS),np.int32); chunks=[]; present=np.zeros(N_EVENTS,bool)
    face=np.zeros(3,np.int64); entries=bad=0; guns=[]
    with uproot.open(path) as f:
      tree=f["sipm_hits"]; schema=",".join(tree.keys())
      for a in tree.iterate(["event_id","global_id","face_type","time_ns","gun_x_mm"],step_size="120 MB",library="np"):
        ev=a["event_id"].astype(int); gid=a["global_id"].astype(int); tm=a["time_ns"]; entries+=len(ev); bad+=np.count_nonzero(~np.isfinite(tm))
        present[np.unique(ev)]=1; np.add.at(counts,(ev,gid),1); face+=np.bincount(a["face_type"],minlength=3); guns.extend(np.unique(a["gun_x_mm"]).tolist())
        chunks.append((ev.copy(),gid.copy(),tm.copy()))
    ev=np.concatenate([c[0] for c in chunks]); gid=np.concatenate([c[1] for c in chunks]); tm=np.concatenate([c[2] for c in chunks])
    key=ev*N_CHANNELS+gid; order=np.lexsort((tm,key)); sk=key[order]; st=tm[order]; bc=np.bincount(sk,minlength=N_EVENTS*N_CHANNELS)
    starts=np.r_[0,np.cumsum(bc[:-1])]; t4=np.full(N_EVENTS*N_CHANNELS,np.nan); ok=bc>=4; t4[ok]=st[starts[ok]+3]; t4=t4.reshape(N_EVENTS,N_CHANNELS)
    lp=pooled_fourth(ev[gid<8],tm[gid<8]); rp=pooled_fourth(ev[(gid>=8)&(gid<16)],tm[(gid>=8)&(gid<16)])
    L,R,T=counts[:,:8],counts[:,8:16],counts[:,16:]; lt,rt,tt=t4[:,:8],t4[:,8:16],t4[:,16:]
    raw,pc,ps,pids=local_pair(T,tt)
    arrays=dict(x_true_mm=np.full(N_EVENTS,x),event_id=np.arange(N_EVENTS),counts=counts,t4=t4,
      t_L_pool=lp,t_R_pool=rp,dt_end_pool=lp-rp,dt_end_weighted=weighted_time(lt,L)-weighted_time(rt,R),
      npe_L_total=L.sum(1),npe_R_total=R.sum(1),npe_T_total=T.sum(1),R_end=np.log(L.sum(1)/R.sum(1)),
      A_end=(L.sum(1)-R.sum(1))/(L.sum(1)+R.sum(1)),x_top_centroid_raw=centroid(T,TOP_X),
      y_L_centroid=centroid(L,END_Y),y_R_centroid=centroid(R,END_Y),y_end_centroid=centroid(L+R,END_Y),
      local_R=raw,local_pair_center=pc,local_pair_status=ps,local_pair_ids=pids)
    np.savez_compressed(out/f"events_x{x:+.0f}mm.npz",**arrays)
    return dict(file=path.name,x_true_mm=x,size_bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,n_entries=entries,
      tree="sipm_hits",schema=schema,n_event_ids_present=present.sum(),event_id_min=np.flatnonzero(present)[0],event_id_max=np.flatnonzero(present)[-1],
      gun_x_values=";".join(map(str,sorted(set(guns)))),nonfinite_time_fraction=bad/entries,hits_end_left=face[0],hits_end_right=face[1],hits_top=face[2])
def load(path): 
    with np.load(path,allow_pickle=True) as z:return {k:z[k] for k in z.files}
def derive(out):
    files=sorted(DATA.glob("photon_hits_x*mm.root"),key=xpos)
    if not np.array_equal([xpos(p) for p in files],EXPECTED): raise RuntimeError("position inventory mismatch")
    (out/"derived/events").mkdir(parents=True,exist_ok=True); rows=[]
    for i,p in enumerate(files): print(f"{i+1}/31 {p.name}",flush=True); rows.append(derive_file(p,out/"derived/events"))
    pd.DataFrame(rows).to_csv(out/"analysis/data_inventory.csv",index=False)

def cv(out):
    files=sorted((out/"derived/events").glob("*.npz"),key=lambda p:load(p)["x_true_mm"][0])
    data=[load(p) for p in files]; methods=["dt_end_pool","dt_end_weighted","R_end","A_end","x_top_centroid_raw","local_R"]
    rows=[]; cal=[]; predictions=[]
    for test,d in enumerate(data):
      train=[q for i,q in enumerate(data) if i!=test]; xt=np.array([q["x_true_mm"][0] for q in train])
      for method in methods:
        if method=="local_R":
          # transferable local coordinate: train true offset relative to selected pair center
          rr=np.concatenate([q[method] for q in train]); yy=np.concatenate([q["x_true_mm"]-q["local_pair_center"] for q in train]); ok=np.isfinite(rr)&np.isfinite(yy)
          coef=np.polyfit(rr[ok],yy[ok],1); pred=d["local_pair_center"]+np.polyval(coef,d[method])
        else:
          means=np.array([np.nanmean(q[method]) for q in train])
          coef=fit_poly(means,xt,1); pred=np.polyval(coef,d[method])
        m=metrics(pred,d["x_true_mm"][0]); rows.append(dict(x_true_mm=d["x_true_mm"][0],method=method,**m))
        cal.append(dict(test_x=d["x_true_mm"][0],method=method,coefficients=json.dumps(coef.tolist()),train_positions=30))
        predictions.append(pd.DataFrame({"x_true_mm":d["x_true_mm"],"event_id":d["event_id"],"method":method,"x_rec":pred}))
    summary=pd.DataFrame(rows); summary.to_csv(out/"analysis/x_reconstruction_summary.csv",index=False)
    pd.DataFrame(cal).to_csv(out/"analysis/cv_calibrations.csv",index=False); pd.concat(predictions).to_csv(out/"analysis/cv_predictions.csv.gz",index=False,compression="gzip")
    pd.DataFrame([{"models":"all predefined linear primary models","selection":"none on outer test fold","folds":31}]).to_csv(out/"analysis/cv_model_selection.csv",index=False)
    return data,summary
def combine(out,data):
    preds=pd.read_csv(out/"analysis/cv_predictions.csv.gz"); base=["dt_end_weighted","R_end","x_top_centroid_raw"]; rows=[]; weights=[]
    for d in data:
      x=d["x_true_mm"][0]; test=preds[preds.x_true_mm==x].pivot(index="event_id",columns="method",values="x_rec")
      train=preds[preds.x_true_mm!=x].pivot_table(index=["x_true_mm","event_id"],columns="method",values="x_rec")
      complete=np.all(np.isfinite(train[base]),axis=1); train=train.loc[complete]
      truth=train.index.get_level_values(0).to_numpy(); resid=train[base].to_numpy()-truth[:,None]; cov=np.cov(resid,rowvar=False); cond=np.linalg.cond(cov)
      try:w=blue(cov); used=base
      except np.linalg.LinAlgError:
        used=base[:-1]; cov=np.cov(resid[:,:-1],rowvar=False); w=blue(cov)
      valid=np.all(np.isfinite(test[used]),axis=1); val=test.loc[valid,used].to_numpy()@w; rows.append(dict(x_true_mm=x,method="BLUE",**metrics(val,x)))
      weights.append(dict(test_x=x,estimators=";".join(used),weights=json.dumps(w.tolist()),condition_number=cond,covariance=json.dumps(cov.tolist())))
    pd.DataFrame(rows).to_csv(out/"analysis/blue_summary.csv",index=False); pd.DataFrame(weights).to_csv(out/"analysis/blue_weights.csv",index=False)
def ystudy(out,data):
    rows=[]
    for d in data:
      x=d["x_true_mm"][0]
      for name in ["y_L_centroid","y_R_centroid","y_end_centroid"]:
        v=d[name]; rows.append(dict(x_true_mm=x,estimator=name,mean=np.nanmean(v),width=np.nanstd(v,ddof=1),median=np.nanmedian(v)))
    pd.DataFrame(rows).to_csv(out/"analysis/y0_feasibility_summary.csv",index=False)
def plots(out):
    s=pd.read_csv(out/"analysis/x_reconstruction_summary.csv"); b=pd.read_csv(out/"analysis/blue_summary.csv"); s=pd.concat([s,b])
    for col,name,ylabel in [("bias","bias_vs_x","Bias (mm)"),("sigma_core","sigma_core_vs_x","Gaussian-core proxy RMS68 (mm)"),("rms68","rms68_vs_x","RMS68 (mm)"),("valid_fraction","valid_fraction_vs_x","Valid fraction")]:
      fig,ax=plt.subplots(figsize=(9,5)); [ax.plot(g.x_true_mm.to_numpy(),g[col].to_numpy(),"o-",ms=3,label=m) for m,g in s.groupby("method")]; ax.set(xlabel="True X (mm)",ylabel=ylabel,title=CONTEXT); ax.grid(alpha=.2); ax.legend(fontsize=7,ncol=2); fig.tight_layout(); fig.savefig(out/f"figures/{name}.pdf"); plt.close(fig)
    y=pd.read_csv(out/"analysis/y0_feasibility_summary.csv")
    for col,name in [("mean","y_centroid_mean_vs_x"),("width","y_centroid_width_vs_x")]:
      fig,ax=plt.subplots(figsize=(9,5)); [ax.plot(g.x_true_mm.to_numpy(),g[col].to_numpy(),"o-",label=m) for m,g in y.groupby("estimator")]; ax.set(xlabel="X (mm)",ylabel=f"Y centroid {col} (mm)",title="y_true=0 feasibility proxy"); ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/f"figures/{name}.pdf");plt.close(fig)
    pred=pd.read_csv(out/"analysis/cv_predictions.csv.gz")
    for method,name in [("dt_end_pool","calibration_end_timing"),("R_end","calibration_end_ratio"),("x_top_centroid_raw","calibration_top_centroid")]:
      g=s[s.method==method];fig,ax=plt.subplots(figsize=(9,5));ax.errorbar(g.x_true_mm.to_numpy(),g.mean_x_rec.to_numpy(),yerr=g.sigma_core.to_numpy(),fmt="o");ax.plot(EXPECTED,EXPECTED,"--");ax.set(xlabel="True X (mm)",ylabel="Reconstructed X (mm)",title=f"{method} LOO prediction — {CONTEXT}");ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/f"figures/{name}.pdf");plt.close(fig)
    chosen=[-690,-450,0,450,690];fig,axs=plt.subplots(1,5,figsize=(16,3.5),sharey=True)
    for ax,x in zip(axs,chosen):
      q=pred[(pred.x_true_mm==x)&(pred.method=="local_R")];v=(q.x_rec-x).dropna().to_numpy()
      if len(v):ax.hist(v,bins=40)
      else:ax.text(.5,.5,"no valid local pair",ha="center",transform=ax.transAxes)
      ax.set_title(f"x={x}");ax.set_xlabel("residual (mm)")
    axs[0].set_ylabel("events");fig.suptitle("Local Top pair residuals; gap position has no valid events");fig.tight_layout();fig.savefig(out/"figures/residual_distributions_selected_positions.pdf");plt.close(fig)
    data=[load(p) for p in sorted((out/"derived/events").glob("*.npz"),key=lambda p:load(p)["x_true_mm"][0])]
    prof=[]
    for d in data:
      c=d["counts"];prof.append([d["x_true_mm"][0],*c[:,:16].mean(0)])
    prof=np.array(prof);fig,ax=plt.subplots(figsize=(9,5))
    for j in range(8):ax.plot(prof[:,0],prof[:,1+j]-prof[:,1+7-j],label=f"L {j}-(7-{j})")
    ax.set(xlabel="X (mm)",ylabel="mean hit-count mirror difference",title="End-channel symmetry at y_true=0");ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/"figures/end_channel_symmetry.pdf");plt.close(fig)
    yl=np.concatenate([d["y_L_centroid"] for d in data]);yr=np.concatenate([d["y_R_centroid"] for d in data]);fig,ax=plt.subplots(figsize=(6,6));ax.hexbin(yl,yr,gridsize=50,mincnt=1);ax.set(xlabel="y_L centroid (mm)",ylabel="y_R centroid (mm)",title="y_true=0 feasibility proxy");fig.tight_layout();fig.savefig(out/"figures/y_left_vs_y_right.pdf");plt.close(fig)
def main():
 p=argparse.ArgumentParser();p.add_argument("stage",choices=["derive","analyze"]);p.add_argument("--output-dir",type=pathlib.Path,required=True);a=p.parse_args()
 for d in ["analysis","figures","tables","derived/events","logs","report"]:(a.output_dir/d).mkdir(parents=True,exist_ok=True)
 if a.stage=="derive":derive(a.output_dir)
 else:
  data,s=cv(a.output_dir);combine(a.output_dir,data);ystudy(a.output_dir,data);plots(a.output_dir)
if __name__=="__main__":main()
