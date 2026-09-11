#!/usr/bin/env python3
"""Phase A+B: k-scan and m-average for END+Vikuiti (all 3 scintillators, 13 positions)."""
import uproot, numpy as np, matplotlib, pandas as pd, time, warnings
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')

REPO   = Path("/home/rrios/ej200")
ENDVIK = REPO / "results/scan_end_vikuiti"
OUT    = REPO / "analysis/optim"
OUT.mkdir(exist_ok=True)

POSITIONS = [-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600]
SCINTS    = ['ej200','ej204','ej230']
COLORS    = {'ej200':'#2196F3','ej204':'#FF5722','ej230':'#4CAF50'}
LABELS    = {'ej200':'EJ-200','ej204':'EJ-204','ej230':'EJ-230'}

def rsigma(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if len(a) < 20: return np.nan
    q1,q3 = np.percentile(a,[25,75])
    return (q3-q1)/1.349

def prep_side(ev, ft, tn, ft_id):
    m = ft == ft_id
    ev_m, tn_m = ev[m], tn[m]
    idx = np.lexsort([tn_m, ev_m])
    ev_s, tn_s = ev_m[idx], tn_m[idx]
    uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
    return uniq, starts, cnts, tn_s

def analyze_file(path, Kmax=250, Mmax=100):
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
    ev, ft, tn = a['event_id'], a['face_type'], a['time_ns']
    n_ev = np.unique(ev).size

    uL,sL,cL,tL = prep_side(ev,ft,tn,0)
    uR,sR,cR,tR = prep_side(ev,ft,tn,1)

    com = np.intersect1d(uL, uR)
    iL  = np.searchsorted(uL, com); iR = np.searchsorted(uR, com)
    cLc, cRc = cL[iL], cR[iR]
    sLc, sRc = sL[iL], sR[iR]

    npe_L, npe_R = float(cLc.mean()), float(cRc.mean())
    Kmax_eff = min(Kmax, max(1, int(np.percentile(cLc,2)), int(np.percentile(cRc,2))))
    Mmax_eff = min(Mmax, Kmax_eff)

    sig_k = np.full(Kmax_eff, np.nan)
    eff_k = np.full(Kmax_eff, np.nan)
    for k in range(1, Kmax_eff+1):
        mk = (cLc>=k)&(cRc>=k)
        eff_k[k-1] = mk.sum()/n_ev
        if mk.sum()<30: continue
        T0 = (tL[sLc[mk]+(k-1)] + tR[sRc[mk]+(k-1)])/2
        sig_k[k-1] = rsigma(T0-T0.mean())*1e3

    sig_m = np.full(Mmax_eff, np.nan)
    eff_m = np.full(Mmax_eff, np.nan)
    for m in range(1, Mmax_eff+1):
        mm = (cLc>=m)&(cRc>=m)
        eff_m[m-1] = mm.sum()/n_ev
        if mm.sum()<30: continue
        iL2 = sLc[mm,None]+np.arange(m); iR2 = sRc[mm,None]+np.arange(m)
        T0 = (tL[iL2].mean(1)+tR[iR2].mean(1))/2
        sig_m[m-1] = rsigma(T0-T0.mean())*1e3

    return dict(sig_k=sig_k, eff_k=eff_k, sig_m=sig_m, eff_m=eff_m,
                npe_L=npe_L, npe_R=npe_R, n_ev=n_ev,
                Kmax_eff=Kmax_eff, Mmax_eff=Mmax_eff)

# ── run all files ─────────────────────────────────────────────────────────────
rows_k = []   # per (scint, x, k) row
rows_opt = [] # per (scint, x) summary

print(f"{'scint':8s} {'x':>6s}  {'k*':>5s}  {'σ(k*) ps':>9s}  {'m*':>5s}  {'σ(m*) ps':>9s}  {'Npe_L':>6s}  {'time':>5s}")
print("-"*75)

for scint in SCINTS:
    for x in POSITIONS:
        fname = ENDVIK/f"photon_hits_{scint}_x{x}.root"
        if not fname.exists(): print(f"  MISSING {fname.name}"); continue
        t0 = time.time()
        r = analyze_file(fname)
        dt = time.time()-t0

        # optimal k
        valid_k = np.isfinite(r['sig_k'])
        k_arr = np.arange(1, r['Kmax_eff']+1)
        if valid_k.any():
            k_opt  = int(k_arr[valid_k][np.argmin(r['sig_k'][valid_k])])
            sig_opt= float(np.nanmin(r['sig_k']))
        else:
            k_opt, sig_opt = 0, np.nan

        # optimal m
        valid_m = np.isfinite(r['sig_m'])
        m_arr = np.arange(1, r['Mmax_eff']+1)
        if valid_m.any():
            m_opt  = int(m_arr[valid_m][np.argmin(r['sig_m'][valid_m])])
            sig_m_opt = float(np.nanmin(r['sig_m']))
        else:
            m_opt, sig_m_opt = 0, np.nan

        print(f"{scint:8s} {x:+6d}  {k_opt:5d}  {sig_opt:9.1f}  {m_opt:5d}  {sig_m_opt:9.1f}  {r['npe_L']:6.0f}  {dt:5.1f}s")

        # collect per-k rows
        for i,k in enumerate(k_arr):
            rows_k.append(dict(scint=scint,x=x,k=int(k),
                               sig_k=float(r['sig_k'][i]),eff_k=float(r['eff_k'][i])))
        for i,m in enumerate(m_arr):
            rows_k[-len(m_arr)+i].update(
                dict(m=int(m),sig_m=float(r['sig_m'][i]),eff_m=float(r['eff_m'][i]))
            ) if i < len(r['sig_m']) else None

        rows_opt.append(dict(scint=scint,x=x,k_opt=k_opt,sig_k_opt=sig_opt,
                             m_opt=m_opt,sig_m_opt=sig_m_opt,
                             npe_L=r['npe_L'],npe_R=r['npe_R'],n_ev=r['n_ev']))

# save csvs
pd.DataFrame(rows_k).to_csv(OUT/'phase_ab_kscan_full.csv', index=False)
df_opt = pd.DataFrame(rows_opt)
df_opt.to_csv(OUT/'phase_ab_optimal.csv', index=False)
print(f"\nSaved: {OUT/'phase_ab_kscan_full.csv'}")
print(f"Saved: {OUT/'phase_ab_optimal.csv'}")

# ── plots ─────────────────────────────────────────────────────────────────────

# Fig 1: σ(k) vs k at x=0 for all scintillators (Phase A)
fig, axes = plt.subplots(1,2,figsize=(12,5))
for scint in SCINTS:
    fname = ENDVIK/f"photon_hits_{scint}_x0.root"
    if not fname.exists(): continue
    r = analyze_file(fname)
    k_arr = np.arange(1, r['Kmax_eff']+1)
    m_arr = np.arange(1, r['Mmax_eff']+1)
    axes[0].plot(k_arr, r['sig_k'], color=COLORS[scint], lw=1.5, label=LABELS[scint])
    axes[1].plot(m_arr, r['sig_m'], color=COLORS[scint], lw=1.5, label=LABELS[scint])
    # mark optimal
    vk = np.isfinite(r['sig_k'])
    if vk.any():
        k_opt = k_arr[vk][np.argmin(r['sig_k'][vk])]
        axes[0].axvline(k_opt, color=COLORS[scint], ls='--', alpha=0.5)
    vm = np.isfinite(r['sig_m'])
    if vm.any():
        m_opt = m_arr[vm][np.argmin(r['sig_m'][vm])]
        axes[1].axvline(m_opt, color=COLORS[scint], ls='--', alpha=0.5)

for ax, lbl in zip(axes, ['k-th photon (T0 = [tL^k+tR^k]/2)','mean of first m photons']):
    ax.set_xlabel(lbl.split()[0]); ax.set_ylabel('σ(T0)  [ps]')
    ax.set_title(lbl); ax.legend(); ax.grid(alpha=0.3); ax.set_ylim(bottom=0)
axes[0].set_xlabel('k'); axes[1].set_xlabel('m')
fig.suptitle('END+Vikuiti timing resolution vs estimator order (x=0)', fontsize=12)
fig.tight_layout()
fig.savefig(OUT/'fig1_sigma_vs_k_x0.png', dpi=150)
print("Saved fig1_sigma_vs_k_x0.png")

# Fig 2: optimal σ and k* vs position
fig, axes = plt.subplots(1,2,figsize=(12,5))
for scint in SCINTS:
    sub = df_opt[df_opt.scint==scint].sort_values('x')
    axes[0].plot(sub.x, sub.sig_k_opt, 'o-', color=COLORS[scint], label=LABELS[scint])
    axes[1].plot(sub.x, sub.k_opt,     'o-', color=COLORS[scint], label=LABELS[scint])
axes[0].set_xlabel('x  [mm]'); axes[0].set_ylabel('σ(T0, k*)  [ps]')
axes[0].set_title('Optimal resolution vs position'); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].set_xlabel('x  [mm]'); axes[1].set_ylabel('k*')
axes[1].set_title('Optimal k* vs position'); axes[1].legend(); axes[1].grid(alpha=0.3)
fig.suptitle('END+Vikuiti: position scan at optimal k', fontsize=12)
fig.tight_layout()
fig.savefig(OUT/'fig2_optimal_vs_x.png', dpi=150)
print("Saved fig2_optimal_vs_x.png")

# Fig 3: σ vs efficiency Pareto at x=0
fig, ax = plt.subplots(figsize=(8,6))
for scint in SCINTS:
    fname = ENDVIK/f"photon_hits_{scint}_x0.root"
    if not fname.exists(): continue
    r = analyze_file(fname)
    k_arr = np.arange(1, r['Kmax_eff']+1)
    mask = np.isfinite(r['sig_k']) & np.isfinite(r['eff_k'])
    ax.plot(r['eff_k'][mask]*100, r['sig_k'][mask], color=COLORS[scint],
            lw=1.5, label=LABELS[scint])
    # mark k=1,2,5,10
    for km in [1,2,5,10,20,50]:
        if km <= r['Kmax_eff']:
            ax.plot(r['eff_k'][km-1]*100, r['sig_k'][km-1], 'o',
                    color=COLORS[scint], ms=5)
            if km in [1,5,20] and np.isfinite(r['sig_k'][km-1]):
                ax.annotate(f'k={km}', (r['eff_k'][km-1]*100, r['sig_k'][km-1]),
                            textcoords='offset points', xytext=(5,3), fontsize=7)
ax.set_xlabel('Efficiency  [%]'); ax.set_ylabel('σ(T0)  [ps]')
ax.set_title('σ vs efficiency trade-off at x=0 (END+Vikuiti)')
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT/'fig3_sigma_vs_efficiency.png', dpi=150)
print("Saved fig3_sigma_vs_efficiency.png")

plt.close('all')
print("\nPhase A+B complete.")
