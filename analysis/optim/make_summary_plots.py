#!/usr/bin/env python3
"""Generate polished summary figures for the Beamer."""
import uproot, numpy as np, matplotlib, pandas as pd, warnings
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

warnings.filterwarnings('ignore')

REPO   = Path("/home/rrios/ej200")
ENDVIK = REPO / "results/scan_end_vikuiti"
OUT    = REPO / "analysis/optim"

COLORS = {'ej200':'#1565C0','ej204':'#BF360C','ej230':'#2E7D32'}
LABELS = {'ej200':'EJ-200  (τ=2.1 ns, 10000 ph/MeV)',
          'ej204':'EJ-204  (τ=1.8 ns, 10400 ph/MeV)',
          'ej230':'EJ-230  (τ=1.5 ns,  9700 ph/MeV)'}
LS     = {'ej200':'-','ej204':'--','ej230':'-.'}

df_opt = pd.read_csv(OUT/'phase_ab_optimal.csv')
df_ks  = pd.read_csv(OUT/'phase_ab_kscan_full.csv')
df_c   = pd.read_csv(OUT/'phase_c_top_ablation.csv')
df_d   = pd.read_csv(OUT/'phase_d_end_ablation.csv')

def rsigma(a):
    a=np.asarray(a,float);a=a[np.isfinite(a)]
    if len(a)<20: return np.nan
    q1,q3=np.percentile(a,[25,75]); return (q3-q1)/1.349

def prep_side(ev,ft,tn,ft_id):
    m=ft==ft_id; ev_m,tn_m=ev[m],tn[m]
    idx=np.lexsort([tn_m,ev_m]); ev_s,tn_s=ev_m[idx],tn_m[idx]
    uniq,starts,cnts=np.unique(ev_s,return_index=True,return_counts=True)
    return uniq,starts,cnts,tn_s

def get_km_curves(path, Kmax=120, Mmax=80):
    with uproot.open(path) as f:
        a=f['sipm_hits'].arrays(['event_id','face_type','time_ns'],library='np')
    ev,ft,tn=a['event_id'],a['face_type'],a['time_ns']
    n_ev=np.unique(ev).size
    uL,sL,cL,tL=prep_side(ev,ft,tn,0); uR,sR,cR,tR=prep_side(ev,ft,tn,1)
    com=np.intersect1d(uL,uR); iL=np.searchsorted(uL,com); iR=np.searchsorted(uR,com)
    cLc,cRc=cL[iL],cR[iR]; sLc,sRc=sL[iL],sR[iR]
    Ke=min(Kmax,max(1,int(np.percentile(cLc,2)),int(np.percentile(cRc,2))))
    Me=min(Mmax,Ke)
    sig_k=np.full(Ke,np.nan); sig_m=np.full(Me,np.nan)
    for k in range(1,Ke+1):
        mk=(cLc>=k)&(cRc>=k)
        if mk.sum()<30: continue
        T0=(tL[sLc[mk]+(k-1)]+tR[sRc[mk]+(k-1)])/2
        sig_k[k-1]=rsigma(T0-T0.mean())*1e3
    for m in range(1,Me+1):
        mm=(cLc>=m)&(cRc>=m)
        if mm.sum()<30: continue
        idxL=sLc[mm,None]+np.arange(m); idxR=sRc[mm,None]+np.arange(m)
        T0=(tL[idxL].mean(1)+tR[idxR].mean(1))/2
        sig_m[m-1]=rsigma(T0-T0.mean())*1e3
    return np.arange(1,Ke+1),sig_k,np.arange(1,Me+1),sig_m

# ──────────────────────────────────────────────────────────────────────────────
# FIG A: σ(k) and σ(m) at x=0 for all scintillators — side by side
# ──────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1,2,figsize=(13,5))
for scint in ['ej200','ej204','ej230']:
    path = ENDVIK/f"photon_hits_{scint}_x0.root"
    k_arr, sk, m_arr, sm = get_km_curves(path)
    c  = COLORS[scint]; ls = LS[scint]
    # k-scan
    vk = np.isfinite(sk)
    axes[0].plot(k_arr[vk], sk[vk], color=c, ls=ls, lw=2, label=LABELS[scint])
    if vk.any():
        k_opt = k_arr[vk][np.argmin(sk[vk])]
        axes[0].axvline(k_opt, color=c, lw=0.8, ls=':', alpha=0.7)
        axes[0].annotate(f'k*={k_opt}', xy=(k_opt, np.nanmin(sk)),
                         xytext=(k_opt+2, np.nanmin(sk)+3), fontsize=8, color=c)
    # m-average
    vm = np.isfinite(sm)
    axes[1].plot(m_arr[vm], sm[vm], color=c, ls=ls, lw=2, label=LABELS[scint])
    if vm.any():
        m_opt = m_arr[vm][np.argmin(sm[vm])]
        axes[1].axvline(m_opt, color=c, lw=0.8, ls=':', alpha=0.7)
        axes[1].annotate(f'm*={m_opt}', xy=(m_opt, np.nanmin(sm)),
                         xytext=(m_opt+1, np.nanmin(sm)+3), fontsize=8, color=c)

for ax, xl, tl in zip(axes,
        ['k  (order of photon)', 'm  (mean of first m photons)'],
        ['Phase A: k-th photon estimator\nT₀ = [t_L^(k) + t_R^(k)] / 2',
         'Phase B: m-average estimator\nT₀ = [mean(t_L^1..m) + mean(t_R^1..m)] / 2']):
    ax.set_xlabel(xl, fontsize=12)
    ax.set_ylabel('σ(T₀)  [ps]', fontsize=12)
    ax.set_title(tl, fontsize=10.5)
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    ax.set_xlim(left=1); ax.set_ylim(40,110)
fig.suptitle('END+Vikuiti timing resolution vs estimator — x = 0 mm', fontsize=13)
fig.tight_layout()
fig.savefig(OUT/'figA_km_curves_x0.png', dpi=150, bbox_inches='tight')
print("figA_km_curves_x0.png")

# ──────────────────────────────────────────────────────────────────────────────
# FIG B: optimal σ vs position — k* and m*
# ──────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1,2,figsize=(13,5))
for scint in ['ej200','ej204','ej230']:
    sub = df_opt[df_opt.scint==scint].sort_values('x')
    c = COLORS[scint]; ls = LS[scint]
    axes[0].plot(sub.x, sub.sig_k_opt, 'o', color=c, ls=ls, ms=6, lw=1.8, label=LABELS[scint].split()[0])
    axes[1].plot(sub.x, sub.sig_m_opt, 'o', color=c, ls=ls, ms=6, lw=1.8, label=LABELS[scint].split()[0])
for ax, tl in zip(axes,['Best k* estimator: σ(T₀, k*)','Best m* estimator: σ(T₀, m*)']):
    ax.set_xlabel('x position  [mm]', fontsize=12)
    ax.set_ylabel('σ(T₀)  [ps]', fontsize=12)
    ax.set_title(tl); ax.legend(fontsize=10); ax.grid(alpha=0.25)
    ax.set_ylim(40,75)
    ax.axhline(59, color='grey', ls=':', lw=1, alpha=0.5, label='59 ps ref')
fig.suptitle('END+Vikuiti: timing resolution vs muon position', fontsize=13)
fig.tight_layout()
fig.savefig(OUT/'figB_sigma_vs_x.png', dpi=150, bbox_inches='tight')
print("figB_sigma_vs_x.png")

# ──────────────────────────────────────────────────────────────────────────────
# FIG C: END ablation — σ vs N SiPMs/end
# ──────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8,5.5))
for scint, c in [('ej204','#BF360C'),('ej230','#2E7D32')]:
    sub0 = df_d[(df_d.scint==scint)&(df_d.x==0)].sort_values('N_sipm_per_end')
    ax.errorbar(sub0.N_sipm_per_end, sub0.sig_mean, yerr=sub0.sig_std,
                fmt='o-', color=c, lw=2, ms=7, capsize=4,
                label=f'{scint.upper()} x=0')
    subm = df_d[(df_d.scint==scint)&(df_d.x==600)].sort_values('N_sipm_per_end')
    ax.errorbar(subm.N_sipm_per_end, subm.sig_mean, yerr=subm.sig_std,
                fmt='s--', color=c, lw=1.5, ms=6, capsize=3, alpha=0.6,
                label=f'{scint.upper()} x=+600 mm')
# add 1/sqrt(N) reference
n_ref = np.linspace(1,8,100)
ax.plot(n_ref, 56*np.sqrt(8/n_ref*2)/np.sqrt(2), 'k:', lw=1.2, alpha=0.5,
        label='1/√Nₚₕ scaling')
ax.set_xlabel('SiPMs per end face  (both ends same)', fontsize=12)
ax.set_ylabel('σ(T₀, k*)  [ps]', fontsize=12)
ax.set_title('Phase D: END SiPM ablation — software masking of local_id\n'
             '(Vikuiti reflector, optimal k scanned independently)', fontsize=10.5)
ax.set_xticks(range(1,9)); ax.grid(alpha=0.25); ax.legend(fontsize=9)
ax.set_ylim(40,165)
fig.tight_layout()
fig.savefig(OUT/'figC_end_ablation.png', dpi=150, bbox_inches='tight')
print("figC_end_ablation.png")

# ──────────────────────────────────────────────────────────────────────────────
# FIG D: TOP ablation — σ vs N top SiPMs
# ──────────────────────────────────────────────────────────────────────────────
CTOP = {'ej204_tir':'#FF9800','ej204_vik':'#BF360C',
        'ej230_tir':'#8BC34A','ej230_vik':'#2E7D32'}
LTOP = {'ej204_tir':'EJ-204 TIR','ej204_vik':'EJ-204 Vikuiti',
        'ej230_tir':'EJ-230 TIR','ej230_vik':'EJ-230 Vikuiti'}
fig, ax = plt.subplots(figsize=(8,5.5))
for key in ['ej230_vik','ej230_tir','ej204_vik','ej204_tir']:
    sub = df_c[df_c.config==key].sort_values('N_sipm')
    sub_ok = sub[sub.sig_mean < 2000]
    ax.errorbar(sub_ok.N_sipm, sub_ok.sig_mean, yerr=sub_ok.sig_std,
                fmt='o-', color=CTOP[key], lw=1.8, ms=6, capsize=3,
                label=LTOP[key])
# reference lines from END+Vik at full SiPMs
ax.axhline(56, color='#2E7D32', ls=':', lw=1.5, alpha=0.7, label='END+Vik EJ-230 (k*=2)')
ax.axhline(60, color='#BF360C', ls=':', lw=1.5, alpha=0.7, label='END+Vik EJ-204 (k*=3)')
ax.set_xlabel('Number of TOP SiPMs', fontsize=12)
ax.set_ylabel('σ(T₀, k*)  [ps]', fontsize=12)
ax.set_title('Phase C: TOP SiPM ablation (scan_resolution x=0)\n'
             'k* optimized for each N; 8 random seeds; dotted = END+Vik reference', fontsize=10.5)
ax.set_xscale('log'); ax.set_ylim(0,700); ax.grid(alpha=0.25); ax.legend(fontsize=9)
ax.set_xticks([1,2,3,5,10,20,35,50,70]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
fig.tight_layout()
fig.savefig(OUT/'figD_top_ablation.png', dpi=150, bbox_inches='tight')
print("figD_top_ablation.png")

# ──────────────────────────────────────────────────────────────────────────────
# FIG E: Pareto frontier — σ vs N_total SiPMs (END and TOP separately)
# ──────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9,6))

# END data: Phase D at x=0 for ej204 and ej230
for scint, c, lbl in [('ej230','#2E7D32','EJ-230 END+Vik'),
                       ('ej204','#BF360C','EJ-204 END+Vik')]:
    sub = df_d[(df_d.scint==scint)&(df_d.x==0)].sort_values('N_total')
    ax.errorbar(sub.N_total, sub.sig_mean, yerr=sub.sig_std,
                fmt='o-', color=c, lw=2, ms=8, capsize=4, label=lbl, zorder=5)

# TOP data: Phase C for ej230_vik and ej204_vik (best reflector per scint)
for key, c, lbl in [('ej230_vik','#8BC34A','EJ-230 TOP+Vik'),
                     ('ej204_vik','#FF9800','EJ-204 TOP+Vik')]:
    sub = df_c[(df_c.config==key)].sort_values('N_sipm')
    sub_ok = sub[sub.sig_mean < 2000]
    ax.errorbar(sub_ok.N_sipm, sub_ok.sig_mean, yerr=sub_ok.sig_std,
                fmt='s--', color=c, lw=1.5, ms=6, capsize=3, label=lbl,
                alpha=0.85, zorder=4)

ax.set_xlabel('Total SiPM count', fontsize=12)
ax.set_ylabel('σ(T₀, k*)  [ps]', fontsize=12)
ax.set_title('Phase E: Pareto frontier — timing resolution vs SiPM count  (x = 0 mm)\n'
             'Filled markers = END+Vikuiti readout   Open markers = TOP readout', fontsize=10.5)
ax.legend(fontsize=10); ax.grid(alpha=0.25)
ax.set_xscale('log')
ax.set_ylim(0,400); ax.set_xlim(1,100)
ax.set_xticks([2,4,8,16,25,50,70]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.axvline(16, color='k', lw=0.8, ls=':', alpha=0.4)
ax.annotate('16 END SiPMs\n(full config)', xy=(16,180), fontsize=8, ha='center', color='k')
fig.tight_layout()
fig.savefig(OUT/'figE_pareto.png', dpi=150, bbox_inches='tight')
print("figE_pareto.png")

plt.close('all')
print("\nAll summary figures done.")
