#!/usr/bin/env python3
"""Phase E: Pareto frontier — σ(T0) vs SiPM count, combining all phases."""
import numpy as np, matplotlib, pandas as pd, warnings
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')

OUT = Path("/home/rrios/ej200/analysis/optim")

# Load results
df_opt = pd.read_csv(OUT/'phase_ab_optimal.csv')
df_c   = pd.read_csv(OUT/'phase_c_top_ablation.csv')
df_d   = pd.read_csv(OUT/'phase_d_end_ablation.csv')

# Build unified Pareto dataframe
points = []

# END+Vikuiti ablation at x=0 (N_total = N_per_end * 2)
for _, row in df_d[df_d.x==0].iterrows():
    points.append(dict(
        label=f"END+Vik {row.scint.upper()} {int(row.N_total)} SiPMs",
        scint=row.scint, config='END+Vik', readout='END',
        N_sipm=int(row.N_total),
        sig=row.sig_mean, sig_err=row.sig_std,
    ))

# TOP configurations at x=0 (N_total = N_sipm top SiPMs, 0 end SiPMs)
for _, row in df_c[df_c.x==0].iterrows():
    points.append(dict(
        label=f"TOP {row.config.upper()} {int(row.N_sipm)} SiPMs",
        scint=row.scint, config=row.config, readout='TOP',
        N_sipm=int(row.N_sipm),
        sig=row.sig_mean, sig_err=row.sig_std,
    ))

# END+Vik full (16 SiPMs) from phase_ab at x=0
for scint in ['ej200','ej204','ej230']:
    sub = df_opt[(df_opt.scint==scint)&(df_opt.x==0)]
    if len(sub)==0: continue
    row = sub.iloc[0]
    points.append(dict(
        label=f"END+Vik {scint.upper()} k*={int(row.k_opt)}",
        scint=scint, config='END+Vik_full', readout='END',
        N_sipm=16, sig=row.sig_k_opt, sig_err=0,
    ))

df_pareto = pd.DataFrame(points)
df_pareto.to_csv(OUT/'phase_e_pareto.csv', index=False)
print("Pareto points:")
print(df_pareto[['label','N_sipm','sig']].sort_values('N_sipm').to_string(index=False))

# Compute Pareto frontier: no other config with fewer SiPMs AND better σ
df_sorted = df_pareto.sort_values('N_sipm')
pareto_mask = []
best_sig = np.inf
for _, row in df_sorted.iterrows():
    if row.sig <= best_sig:
        best_sig = row.sig
        pareto_mask.append(True)
    else:
        pareto_mask.append(False)
df_pareto['pareto'] = pareto_mask

# ── plot ───────────────────────────────────────────────────────────────────────
COLORS_R = {'END':'#E91E63','TOP':'#2196F3'}
MARKERS  = {'ej200':'o','ej204':'s','ej230':'^','ej204_tir':'D','ej204_vik':'s',
            'ej230_tir':'D','ej230_vik':'^','END+Vik_full':'*'}

fig, ax = plt.subplots(figsize=(11,7))

for readout, grp in df_pareto.groupby('readout'):
    c = COLORS_R[readout]
    for _, row in grp.iterrows():
        mk = MARKERS.get(row.scint, 'o') if readout=='END' else MARKERS.get(row.config,'o')
        kws = dict(ms=9 if row.pareto else 6, zorder=5 if row.pareto else 3,
                   mec='k' if row.pareto else c, mew=1.5 if row.pareto else 0.5)
        ax.errorbar(row.N_sipm, row.sig, yerr=row.sig_err,
                    fmt=mk, color=c, **kws)

# Draw Pareto frontier line
pf = df_pareto[df_pareto.pareto].sort_values('N_sipm')
ax.step(pf.N_sipm, pf.sig, where='post', color='k', lw=1, ls='--', alpha=0.5, label='Pareto frontier')

# Custom legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
legend_els = [
    Patch(color='#E91E63', label='END readout'),
    Patch(color='#2196F3', label='TOP readout'),
    Line2D([0],[0],color='k',ls='--',lw=1,label='Pareto frontier'),
    Line2D([0],[0],marker='o',color='w',mec='k',mfc='#E91E63',ms=9,lw=0,label='EJ-200'),
    Line2D([0],[0],marker='s',color='w',mec='k',mfc='#E91E63',ms=9,lw=0,label='EJ-204'),
    Line2D([0],[0],marker='^',color='w',mec='k',mfc='#E91E63',ms=9,lw=0,label='EJ-230'),
    Line2D([0],[0],marker='D',color='w',mec='k',mfc='#2196F3',ms=9,lw=0,label='TIR'),
    Line2D([0],[0],marker='s',color='w',mec='k',mfc='#2196F3',ms=9,lw=0,label='Vikuiti'),
]
ax.legend(handles=legend_els, loc='upper right', fontsize=9, ncol=2)

ax.set_xlabel('Total SiPM count', fontsize=13)
ax.set_ylabel('σ(T0)  [ps]', fontsize=13)
ax.set_title('Pareto frontier: timing resolution vs SiPM count  (x=0)', fontsize=13)
ax.grid(alpha=0.3)
ax.set_xscale('log')
ax.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(OUT/'fig6_pareto_frontier.png', dpi=150)
print("\nSaved fig6_pareto_frontier.png")

# Summary table
print("\n=== PARETO-OPTIMAL CONFIGURATIONS ===")
print(f"{'Config':45s}  {'N_SiPM':>7s}  {'σ [ps]':>8s}")
print("-"*65)
for _, row in pf.iterrows():
    print(f"{row.label:45s}  {row.N_sipm:7d}  {row.sig:8.1f}")
