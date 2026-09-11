#!/usr/bin/env python3
"""
EndSparseTop analysis: END+Vikuiti + N sparse top SiPMs.
Computes:
  - T0_end   : m-average of end SiPMs (face_type 0,1)
  - T0_top   : k-th order stat of top SiPMs (face_type 2)
  - T0_comb  : BLUE-weighted combination of T0_end and T0_top
Outputs:
  phase_sparse_top_full.csv   — σ vs m/k for every config+position
  phase_sparse_top_optimal.csv — best σ per config+position
"""
import uproot, numpy as np, pandas as pd
from pathlib import Path
from itertools import product

REPO   = Path("/home/rrios/ej200")
INDIR  = REPO / "results/scan_end_vik_sparse_top_v2"
ENDVIK = REPO / "results/scan_end_vik_sparse_top_v2"   # same dir, has all faces
ENDVIK_REF = REPO / "results/scan_end_vikuiti"          # End-only reference
OUT    = REPO / "analysis/optim"

SCINTS = ["ej230", "ej204"]
NTOPS  = [4, 8, 14, 20]
XPOS   = list(range(-600, 601, 100))

def rsigma(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 20:
        return np.nan
    q1, q3 = np.percentile(a, [25, 75])
    return (q3 - q1) / 1.349

def prep_side(ev, ft, tn, ft_id):
    m = ft == ft_id
    ev_m, tn_m = ev[m], tn[m]
    idx = np.lexsort([tn_m, ev_m])
    ev_s, tn_s = ev_m[idx], tn_m[idx]
    uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
    return uniq, starts, cnts, tn_s

def analyse_file(path, Mmax=80, Kmax_end=80, Kmax_top=50):
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id', 'face_type', 'time_ns'], library='np')
    ev, ft, tn = a['event_id'], a['face_type'], a['time_ns']

    uL, sL, cL, tL = prep_side(ev, ft, tn, 0)
    uR, sR, cR, tR = prep_side(ev, ft, tn, 1)
    uT, sT, cT, tTop = prep_side(ev, ft, tn, 2)

    # Common end events
    com_end = np.intersect1d(uL, uR)
    iL = np.searchsorted(uL, com_end); iR = np.searchsorted(uR, com_end)
    cLc, cRc = cL[iL], cR[iR]
    sLc, sRc = sL[iL], sR[iR]

    # Common end+top events (events that have both end and top hits)
    com_all = np.intersect1d(com_end, uT)
    ia = np.searchsorted(com_end, com_all)
    iT = np.searchsorted(uT, com_all)
    cLa, cRa = cLc[ia], cRc[ia]
    sLa, sRa = sLc[ia], sRc[ia]
    cTa = cT[iT]; sTa = sT[iT]

    # ── T0_end: m-average on end SiPMs ──
    Me = min(Mmax, max(1, int(np.percentile(cLc, 2)), int(np.percentile(cRc, 2))))
    sig_end = np.full(Me, np.nan)
    for m in range(1, Me + 1):
        mask = (cLc >= m) & (cRc >= m)
        if mask.sum() < 30:
            continue
        iLi = sLc[mask, None] + np.arange(m)
        iRi = sRc[mask, None] + np.arange(m)
        T0 = (tL[iLi].mean(1) + tR[iRi].mean(1)) / 2.0
        sig_end[m - 1] = rsigma(T0 - T0.mean()) * 1e3

    # ── T0_top: k-th order stat on top SiPMs ──
    if len(uT) > 30:
        Kt = min(Kmax_top, max(1, int(np.percentile(cT, 2))))
        sig_top = np.full(Kt, np.nan)
        for k in range(1, Kt + 1):
            mask = cTa >= k
            if mask.sum() < 30:
                continue
            T0t = tTop[sTa[mask] + (k - 1)]
            sig_top[k - 1] = rsigma(T0t - T0t.mean()) * 1e3
    else:
        Kt = 0
        sig_top = np.array([])

    # ── T0_comb: BLUE combination of T0_end(m*) and T0_top(k*) ──
    # Use optimal m* and k*, then BLUE-combine
    sig_comb_best = np.nan
    m_opt_comb = np.nan
    k_opt_comb = np.nan
    if np.any(np.isfinite(sig_end)) and len(sig_top) > 0 and np.any(np.isfinite(sig_top)):
        m_opt = int(np.nanargmin(sig_end)) + 1
        k_opt = int(np.nanargmin(sig_top)) + 1
        # Compute T0_end and T0_top vectors for BLUE
        mask_e = (cLa >= m_opt) & (cRa >= m_opt) & (cTa >= k_opt)
        if mask_e.sum() >= 30:
            iLi = sLa[mask_e, None] + np.arange(m_opt)
            iRi = sRa[mask_e, None] + np.arange(m_opt)
            T0e = (tL[iLi].mean(1) + tR[iRi].mean(1)) / 2.0
            T0t = tTop[sTa[mask_e] + (k_opt - 1)]
            # BLUE: T0_comb = w*T0e + (1-w)*T0t, w = sigma_t^2 / (sigma_e^2 + sigma_t^2)
            se2 = rsigma(T0e - T0e.mean()) ** 2
            st2 = rsigma(T0t - T0t.mean()) ** 2
            cov = np.cov(T0e, T0t)[0, 1]
            denom = se2 + st2 - 2 * cov
            if denom > 0:
                w = (st2 - cov) / denom
                w = np.clip(w, 0.0, 1.0)
            else:
                w = 0.5
            T0c = w * T0e + (1 - w) * T0t
            sig_comb_best = rsigma(T0c - T0c.mean()) * 1e3
            m_opt_comb = m_opt
            k_opt_comb = k_opt

    return {
        'sig_end': sig_end,
        'sig_top': sig_top,
        'sig_end_opt': np.nanmin(sig_end) if np.any(np.isfinite(sig_end)) else np.nan,
        'm_opt': int(np.nanargmin(sig_end)) + 1 if np.any(np.isfinite(sig_end)) else np.nan,
        'sig_top_opt': np.nanmin(sig_top) if (len(sig_top) > 0 and np.any(np.isfinite(sig_top))) else np.nan,
        'k_opt': int(np.nanargmin(sig_top)) + 1 if (len(sig_top) > 0 and np.any(np.isfinite(sig_top))) else np.nan,
        'sig_comb': sig_comb_best,
        'm_comb': m_opt_comb,
        'k_comb': k_opt_comb,
        'npe_L': cLc.mean() if len(cLc) > 0 else np.nan,
        'npe_R': cRc.mean() if len(cRc) > 0 else np.nan,
        'npe_T': cT.mean() if len(cT) > 0 else np.nan,
    }

# ── End-only reference at x=0 ────────────────────────────────────────────────
print("Loading End+Vik reference (scan_end_vikuiti)...")
endvik_ref = {}
for scint in SCINTS:
    path = ENDVIK_REF / f"photon_hits_{scint}_x0.root"
    if not path.exists():
        continue
    res = analyse_file(path)
    endvik_ref[scint] = res['sig_end_opt']
    print(f"  {scint}: σ_end_ref = {res['sig_end_opt']:.1f} ps (m*={res['m_opt']})")

# ── Main loop ─────────────────────────────────────────────────────────────────
rows_opt = []
rows_full = []

total = len(SCINTS) * len(NTOPS) * len(XPOS)
done = 0
for scint, ntop, x in product(SCINTS, NTOPS, XPOS):
    tag = f"{scint}_ntop{ntop}"
    path = INDIR / f"photon_hits_{tag}_x{x}.root"
    if not path.exists():
        done += 1
        continue
    res = analyse_file(path)
    done += 1
    if done % 10 == 0:
        print(f"  {done}/{total}: {scint} N={ntop} x={x} → "
              f"σ_end={res['sig_end_opt']:.1f} σ_top={res['sig_top_opt']:.1f} "
              f"σ_comb={res['sig_comb']:.1f} ps")

    rows_opt.append({
        'scint': scint, 'ntop': ntop, 'x': x,
        'n_total_sipm': 16 + ntop,
        'sig_end_opt': res['sig_end_opt'], 'm_opt': res['m_opt'],
        'sig_top_opt': res['sig_top_opt'], 'k_opt': res['k_opt'],
        'sig_comb': res['sig_comb'],
        'npe_L': res['npe_L'], 'npe_R': res['npe_R'], 'npe_T': res['npe_T'],
    })

    for m, s in enumerate(res['sig_end'], 1):
        rows_full.append({'scint': scint, 'ntop': ntop, 'x': x,
                          'estimator': 'end_m', 'param': m, 'sigma': s})
    for k, s in enumerate(res['sig_top'], 1):
        rows_full.append({'scint': scint, 'ntop': ntop, 'x': x,
                          'estimator': 'top_k', 'param': k, 'sigma': s})

df_opt  = pd.DataFrame(rows_opt)
df_full = pd.DataFrame(rows_full)
df_opt.to_csv(OUT / 'phase_sparse_top_optimal.csv', index=False)
df_full.to_csv(OUT / 'phase_sparse_top_full.csv', index=False)
print(f"\nWrote phase_sparse_top_optimal.csv ({len(df_opt)} rows)")
print(f"Wrote phase_sparse_top_full.csv ({len(df_full)} rows)")

# ── Quick summary table ───────────────────────────────────────────────────────
print("\n── x=0 summary ──────────────────────────────────────────────────────")
print(f"{'Config':25s}  {'N_tot':5s}  {'σ_end':7s}  {'σ_top':7s}  {'σ_comb':7s}")
print("-" * 60)
for scint in SCINTS:
    ref = endvik_ref.get(scint, np.nan)
    print(f"{'End+Vik ref':25s}  {'16':5s}  {ref:7.1f}  {'—':>7s}  {'—':>7s}")
    sub = df_opt[(df_opt.scint == scint) & (df_opt.x == 0)].sort_values('ntop')
    for _, row in sub.iterrows():
        print(f"  EST N={int(row.ntop):<3d} ({scint}):    "
              f"  {int(row.n_total_sipm):5d}  "
              f"{row.sig_end_opt:7.1f}  "
              f"{row.sig_top_opt:7.1f}  "
              f"{row.sig_comb:7.1f}")
    print()
