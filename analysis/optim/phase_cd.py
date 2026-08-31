#!/usr/bin/env python3
"""Phase C: TOP SiPM ablation (scan_resolution). Phase D: END SiPM ablation (scan_end_vikuiti)."""
import uproot, numpy as np, matplotlib, pandas as pd, time, warnings
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')
rng = np.random.default_rng(42)

REPO    = Path("/home/rrios/ej200")
SCANRES = REPO / "results/scan_resolution"
ENDVIK  = REPO / "results/scan_end_vikuiti"
OUT     = REPO / "analysis/optim"
OUT.mkdir(exist_ok=True)

COLORS = {'ej204_tir':'#FF9800','ej204_vik':'#FF5722',
          'ej230_tir':'#8BC34A','ej230_vik':'#4CAF50'}
LABELS = {'ej204_tir':'EJ-204 TIR','ej204_vik':'EJ-204 Vik',
          'ej230_tir':'EJ-230 TIR','ej230_vik':'EJ-230 Vik'}
SCINTS_END  = ['ej204','ej230']
COLORS_END  = {'ej200':'#2196F3','ej204':'#FF5722','ej230':'#4CAF50'}
LABELS_END  = {'ej200':'EJ-200','ej204':'EJ-204','ej230':'EJ-230'}

def rsigma(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if len(a) < 20: return np.nan
    q1,q3 = np.percentile(a,[25,75])
    return (q3-q1)/1.349

# ── Phase C: TOP ablation ──────────────────────────────────────────────────────
N_LEVELS_TOP = [70,50,35,25,20,14,10,7,5,4,3,2,1]
N_SEEDS      = 8
KMAX_TOP     = 25

def analyze_top_ablation(path):
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','local_id','time_ns'], library='np')
    ev, ft, lid, tn = a['event_id'], a['face_type'], a['local_id'], a['time_ns']
    n_ev = np.unique(ev).size

    mask_top = ft == 2
    ev_t, lid_t, tn_t = ev[mask_top], lid[mask_top], tn[mask_top]
    idx = np.lexsort([tn_t, ev_t])
    ev_t, lid_t, tn_t = ev_t[idx], lid_t[idx], tn_t[idx]
    all_lids = np.unique(lid_t)
    N_full   = len(all_lids)

    results = []
    for N in N_LEVELS_TOP:
        N = min(N, N_full)
        sigs, kpts = [], []
        for seed in range(N_SEEDS):
            if N == N_full:
                chosen = all_lids
            else:
                chosen = rng.choice(all_lids, size=N, replace=False)
            sub = np.isin(lid_t, chosen)
            ev_s, tn_s = ev_t[sub], tn_t[sub]
            uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
            Ke = min(KMAX_TOP, max(1, int(np.percentile(cnts,5)))) if len(cnts) else 1
            best_s, best_k = np.inf, 1
            for k in range(1, Ke+1):
                mk = cnts >= k
                if mk.sum() < 30: continue
                T0 = tn_s[starts[mk]+(k-1)]
                s  = rsigma(T0-T0.mean())*1e3
                if np.isfinite(s) and s < best_s: best_s,best_k = s,k
            sigs.append(best_s); kpts.append(best_k)
        results.append(dict(N_sipm=N, sig_mean=np.nanmean(sigs), sig_std=np.nanstd(sigs),
                            k_opt=float(np.nanmean(kpts)), npe_mean=ev_t[np.isin(lid_t,chosen)].size/n_ev))
    return results

# ── Phase D: END ablation ──────────────────────────────────────────────────────
N_LEVELS_END = [8,7,6,5,4,3,2,1]
N_SEEDS_END  = 8
KMAX_END     = 200

def analyze_end_ablation(path):
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','local_id','time_ns'], library='np')
    ev, ft, lid, tn = a['event_id'], a['face_type'], a['local_id'], a['time_ns']
    n_ev = np.unique(ev).size
    all_lids = np.arange(8)

    results = []
    for N in N_LEVELS_END:
        sigs, kpts = [], []
        for seed in range(N_SEEDS_END):
            chosen = all_lids if N==8 else rng.choice(all_lids, size=N, replace=False)
            sides = {}
            for ft_id in [0,1]:
                m = (ft==ft_id) & np.isin(lid, chosen)
                ev_m, tn_m = ev[m], tn[m]
                if len(ev_m)==0: sides[ft_id]=None; continue
                idx = np.lexsort([tn_m, ev_m])
                ev_s, tn_s = ev_m[idx], tn_m[idx]
                uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
                sides[ft_id] = dict(uniq=uniq, starts=starts, cnts=cnts, t=tn_s)
            if sides[0] is None or sides[1] is None: continue
            com = np.intersect1d(sides[0]['uniq'], sides[1]['uniq'])
            iL  = np.searchsorted(sides[0]['uniq'], com)
            iR  = np.searchsorted(sides[1]['uniq'], com)
            cLc = sides[0]['cnts'][iL]; cRc = sides[1]['cnts'][iR]
            sLc = sides[0]['starts'][iL]; sRc = sides[1]['starts'][iR]
            tL  = sides[0]['t'];        tR  = sides[1]['t']
            Ke  = min(KMAX_END, max(1, int(np.percentile(cLc,2)), int(np.percentile(cRc,2))))
            best_s, best_k = np.inf, 1
            for k in range(1, Ke+1):
                mk = (cLc>=k)&(cRc>=k)
                if mk.sum()<30: continue
                T0 = (tL[sLc[mk]+(k-1)] + tR[sRc[mk]+(k-1)])/2
                s  = rsigma(T0-T0.mean())*1e3
                if np.isfinite(s) and s < best_s: best_s,best_k = s,k
            sigs.append(best_s); kpts.append(best_k)
        results.append(dict(N_sipm_per_end=N, N_total=N*2,
                            sig_mean=np.nanmean(sigs), sig_std=np.nanstd(sigs),
                            k_opt=float(np.nanmean(kpts))))
    return results

# ── run Phase C ───────────────────────────────────────────────────────────────
configs_c = [
    ('ej204_tir','ej204','tir'),('ej204_vik','ej204','vik'),
    ('ej230_tir','ej230','tir'),('ej230_vik','ej230','vik'),
]
rows_c = []
print("="*60); print("PHASE C: TOP SiPM ablation"); print("="*60)
for key, scint, refl in configs_c:
    fname = SCANRES/f"photon_hits_{scint}_{refl}_x0.root"
    if not fname.exists(): print(f"  MISSING {fname.name}"); continue
    t0 = time.time()
    res = analyze_top_ablation(fname)
    dt = time.time()-t0
    print(f"  {key:14s}  ({dt:.1f}s)")
    for r in res:
        r.update(dict(config=key,scint=scint,reflector=refl,x=0))
        rows_c.append(r)
        print(f"    N={r['N_sipm']:3d}  σ={r['sig_mean']:6.1f}±{r['sig_std']:.1f}ps  k*={r['k_opt']:.1f}")

df_c = pd.DataFrame(rows_c)
df_c.to_csv(OUT/'phase_c_top_ablation.csv', index=False)
print(f"\nSaved: {OUT/'phase_c_top_ablation.csv'}")

# ── run Phase D ───────────────────────────────────────────────────────────────
rows_d = []
print(); print("="*60); print("PHASE D: END SiPM ablation"); print("="*60)
for scint in ['ej204','ej230']:
    for x in [0, -600, 600]:
        fname = ENDVIK/f"photon_hits_{scint}_x{x}.root"
        if not fname.exists(): print(f"  MISSING {fname.name}"); continue
        t0 = time.time()
        res = analyze_end_ablation(fname)
        dt = time.time()-t0
        print(f"  {scint:8s} x={x:+4d}  ({dt:.1f}s)")
        for r in res:
            r.update(dict(scint=scint, x=x))
            rows_d.append(r)
            print(f"    N={r['N_sipm_per_end']:2d}/end  σ={r['sig_mean']:6.1f}±{r['sig_std']:.1f}ps  k*={r['k_opt']:.1f}")

df_d = pd.DataFrame(rows_d)
df_d.to_csv(OUT/'phase_d_end_ablation.csv', index=False)
print(f"\nSaved: {OUT/'phase_d_end_ablation.csv'}")

# ── plots ─────────────────────────────────────────────────────────────────────

# Fig 4: TOP ablation σ vs N_sipm
fig, axes = plt.subplots(1,2,figsize=(12,5))
for key,scint,refl in configs_c:
    sub = df_c[df_c.config==key].sort_values('N_sipm')
    c = COLORS[key]
    axes[0].errorbar(sub.N_sipm, sub.sig_mean, yerr=sub.sig_std,
                     fmt='o-', color=c, label=LABELS[key], lw=1.5, ms=5)
    axes[1].plot(sub.N_sipm, sub.k_opt, 'o-', color=c, label=LABELS[key], lw=1.5)
for ax in axes:
    ax.set_xlabel('N TOP SiPMs'); ax.legend(); ax.grid(alpha=0.3)
axes[0].set_ylabel('σ(T0, k*)  [ps]'); axes[0].set_title('TOP SiPM ablation: timing resolution')
axes[1].set_ylabel('k*');              axes[1].set_title('Optimal k vs N SiPMs')
fig.suptitle('Phase C: TOP SiPM ablation (scan_resolution x=0)', fontsize=12)
fig.tight_layout(); fig.savefig(OUT/'fig4_top_ablation.png', dpi=150)
print("Saved fig4_top_ablation.png")

# Fig 5: END ablation σ vs N_sipm/end (x=0)
fig, ax = plt.subplots(figsize=(8,6))
for scint in ['ej204','ej230']:
    sub = df_d[(df_d.scint==scint)&(df_d.x==0)].sort_values('N_sipm_per_end')
    c = COLORS_END[scint]
    ax.errorbar(sub.N_sipm_per_end, sub.sig_mean, yerr=sub.sig_std,
                fmt='o-', color=c, label=LABELS_END[scint], lw=1.5, ms=5)
ax.set_xlabel('SiPMs per end face'); ax.set_ylabel('σ(T0, k*)  [ps]')
ax.set_title('END+Vikuiti: ablation (software masking, x=0)')
ax.legend(); ax.grid(alpha=0.3)
ax.set_xticks(N_LEVELS_END)
fig.tight_layout(); fig.savefig(OUT/'fig5_end_ablation.png', dpi=150)
print("Saved fig5_end_ablation.png")

plt.close('all')
print("\nPhase C+D complete.")
