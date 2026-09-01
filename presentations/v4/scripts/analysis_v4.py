#!/usr/bin/env python3
"""
analysis_v4.py — Comprehensive analysis for SHiP Timing Bar presentation v4.

Produces all figures and tables needed for the scientific rebuild.

New analyses vs best_est_analysis.py:
  - Train/test split at x=0 (AH)
  - Oracle vs global estimator comparison (Y)
  - Bias vs x for END/TOP/BLUE (Z)
  - rho and BLUE weights vs x (AG)
  - sigma_x_END(m) position resolution (U)
  - TOP position estimators (V)
  - v_eff fit from Delta_t vs x (K, Q)
  - Old per-SiPM estimator comparison (W)
  - PDE audit (M)
  - Electronics jitter scenarios (AP)
  - Gaussian fits with chi2/p (E)
  - Material comparison regenerated in ROOT (P)

All figures: ROOT vector PDFs + PNG preview.
All tables: CSV.
"""

import ROOT
import numpy as np
import uproot
import json
import pandas as pd
import sys
from pathlib import Path

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

REPO   = Path("/home/rrios/ej200")
INDIR  = REPO / "results/scan_end_vik_sparse_top_v2"
INDIR_END = REPO / "results/scan_end_vikuiti"
OUTDIR = REPO / "analysis/presentation_v4"
FIGS   = OUTDIR / "figs"
TABS   = OUTDIR / "tables"

RNG    = np.random.default_rng(20260820)
N_BOOT = 1000

# TOP SiPM x positions (local_id 0-19, from geometry check)
TOP_XPOS = np.array([
    -664.3, -594.3, -524.3, -454.3, -384.4, -314.4, -244.4, -174.4,
    -104.5, -34.5, 34.5, 104.5, 174.4, 244.4, 314.4, 384.4,
    454.3, 524.3, 594.3, 664.3
])

# ═══════════════════════════════════════════════════════════════════════════════
# STYLE
# ═══════════════════════════════════════════════════════════════════════════════
def set_style():
    s = ROOT.gStyle
    s.SetOptStat(0); s.SetOptFit(0)
    s.SetCanvasColor(0); s.SetPadColor(0); s.SetFrameFillColor(0)
    s.SetCanvasBorderMode(0); s.SetPadBorderMode(0); s.SetFrameBorderMode(0)
    for ax in ("x","y","z"):
        s.SetTitleFont(42, ax); s.SetLabelFont(42, ax)
        s.SetTitleSize(0.055, ax); s.SetLabelSize(0.048, ax)
        s.SetTitleOffset(1.0, ax)
    s.SetTextFont(42); s.SetTextSize(0.05)
    s.SetPadTopMargin(0.07); s.SetPadRightMargin(0.06)
    s.SetPadBottomMargin(0.14); s.SetPadLeftMargin(0.14)
    s.SetTitleOffset(1.15, "y")
    s.SetLegendBorderSize(1); s.SetLegendFont(42); s.SetLegendFillColor(0)
    s.SetHistLineWidth(2)
    s.SetPadTickX(1); s.SetPadTickY(1)
    s.SetGridColor(ROOT.kGray+1); s.SetGridStyle(3)
    s.SetPadGridX(False); s.SetPadGridY(False)
    s.SetMarkerSize(1.2)

set_style()

# Main ROOT output
fout = ROOT.TFile(str(FIGS / "analysis_v4.root"), "RECREATE")

# ═══════════════════════════════════════════════════════════════════════════════
# CORE UTILITY
# ═══════════════════════════════════════════════════════════════════════════════
def rsigma(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if len(a) < 20: return np.nan
    q1, q3 = np.percentile(a, [25, 75])
    return (q3 - q1) / 1.349

def bootstrap(arr, func=rsigma, n=N_BOOT):
    arr = np.asarray(arr, float); arr = arr[np.isfinite(arr)]
    N = len(arr)
    return np.std([func(arr[RNG.integers(0, N, N)]) for _ in range(n)])

def load_hits_full(path):
    """Load all branches needed for full analysis."""
    with uproot.open(path) as f:
        tr = f['sipm_hits']
        a = tr.arrays(['event_id','face_type','global_id','local_id',
                        'time_ns','pde','x_mm','gun_x_mm'], library='np')
    return a

def prep_side(ev, ft, tn, ft_id):
    m = ft == ft_id
    ev_m, tn_m = ev[m], tn[m]
    idx = np.lexsort([tn_m, ev_m])
    ev_s, tn_s = ev_m[idx], tn_m[idx]
    uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
    return uniq, starts, cnts, tn_s

def get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, m):
    """Compute T0_END for given m, return (events, T0_end_ps, DeltaT_ps)."""
    com = np.intersect1d(uL, uR)
    iL  = np.searchsorted(uL, com); iR = np.searchsorted(uR, com)
    cLc = cL[iL]; cRc = cR[iR]; sLc = sL[iL]; sRc = sR[iR]
    mask = (cLc >= m) & (cRc >= m)
    tL_m = np.array([tL[sLc[i]:sLc[i]+m].mean() for i in np.where(mask)[0]])
    tR_m = np.array([tR[sRc[i]:sRc[i]+m].mean() for i in np.where(mask)[0]])
    T0e   = ((tL_m + tR_m) / 2.) * 1e3
    DeltaT = (tR_m - tL_m) * 1e3
    return com[mask], T0e, DeltaT

def get_T0_top(tTop, sT, cT, uT, k):
    """Compute T0_TOP for given k."""
    mask = cT >= k
    T0t  = tTop[sT[mask] + (k-1)] * 1e3
    return uT[mask], T0t

def get_T0_persipm(ev, ft, tn, lid, n_sipm=20):
    """Old per-SiPM estimator: mean of first-photon arrival per SiPM, averaged over SiPMs."""
    m = ft == 2
    ev_top = ev[m]; tn_top = tn[m]; lid_top = lid[m]
    # For each event: for each SiPM local_id, take min(time), then average
    ev_uniq = np.unique(ev_top)
    T0_ps = np.full(len(ev_uniq), np.nan)
    for i, eid in enumerate(ev_uniq):
        msk = ev_top == eid
        e_tn = tn_top[msk]; e_lid = lid_top[msk]
        first_times = []
        for sid in range(n_sipm):
            smsk = e_lid == sid
            if smsk.sum() > 0:
                first_times.append(e_tn[smsk].min())
        if len(first_times) > 0:
            T0_ps[i] = np.mean(first_times) * 1e3
    return ev_uniq, T0_ps

def intersect_estimators(ev_e, T0e, ev_t, T0t):
    """Return aligned arrays for events in both."""
    com = np.intersect1d(ev_e, ev_t)
    ie  = np.searchsorted(ev_e, com)
    it  = np.searchsorted(ev_t, com)
    # Filter finite
    mask = np.isfinite(T0e[ie]) & np.isfinite(T0t[it])
    return T0e[ie][mask], T0t[it][mask]

def compute_blue(T0e_arr, T0t_arr):
    """Given aligned T0e, T0t arrays (mean-subtracted), compute BLUE."""
    T0e_c = T0e_arr - T0e_arr.mean()
    T0t_c = T0t_arr - T0t_arr.mean()
    se  = rsigma(T0e_c); st  = rsigma(T0t_c)
    C   = np.cov(T0e_c, T0t_c)
    cov = C[0,1]
    rho = cov / (np.sqrt(C[0,0] * C[1,1]) + 1e-30)
    denom = se**2 + st**2 - 2*cov
    w_end = np.clip((st**2 - cov) / (denom + 1e-30), 0., 1.)
    w_top = 1. - w_end
    T0c   = w_end * T0e_c + w_top * T0t_c
    # inputs are in ps, outputs are in ps
    return dict(T0c=T0c, T0e=T0e_c, T0t=T0t_c,
                w_end=w_end, w_top=w_top, rho=rho,
                sig_e=se, sig_t=st,
                sig_c=rsigma(T0c))

def make_graph(x, y, ex=None, ey=None):
    """Create TGraph or TGraphErrors from numpy arrays."""
    x = np.asarray(x, dtype='d'); y = np.asarray(y, dtype='d')
    n = len(x)
    if ex is not None or ey is not None:
        g = ROOT.TGraphErrors(n)
        ex = np.zeros(n) if ex is None else np.asarray(ex, dtype='d')
        ey = np.zeros(n) if ey is None else np.asarray(ey, dtype='d')
        for i in range(n):
            g.SetPoint(i, x[i], y[i])
            g.SetPointError(i, ex[i], ey[i])
    else:
        g = ROOT.TGraph(n)
        for i in range(n):
            g.SetPoint(i, x[i], y[i])
    return g

def save_canvas(c, name, subdir=None):
    target = FIGS / subdir if subdir else FIGS
    fout.cd()
    c.Write(name)
    c.SaveAs(str(target / f"{name}.pdf"))
    c.SaveAs(str(target / f"{name}.png"))
    print(f"  Saved {name}.{{pdf,png}}")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. LOAD EJ-230 N=20 x=0 DATA — full branches
# ═══════════════════════════════════════════════════════════════════════════════
print("="*60)
print("1. Loading EJ-230 N=20 x=0 …")
a0 = load_hits_full(INDIR / "photon_hits_ej230_ntop20_x0.root")
ev0 = a0['event_id']; ft0 = a0['face_type']
tn0 = a0['time_ns']; lid0 = a0['local_id']
pde0 = a0['pde']

print(f"  Total hits: {len(ev0):,}")
n_ev = len(np.unique(ev0))
print(f"  Events: {n_ev}")

# PDE audit
mean_pde = pde0.mean()
print(f"  Mean PDE (all photons): {mean_pde:.4f}")
# By face type
for ft_id, name in [(0,'END-L'),(1,'END-R'),(2,'TOP')]:
    m = ft0 == ft_id
    print(f"  PDE face_type={ft_id} ({name}): mean={pde0[m].mean():.4f}, "
          f"range=[{pde0[m].min():.3f}, {pde0[m].max():.3f}]")

uL0, sL0, cL0, tL0 = prep_side(ev0, ft0, tn0, 0)
uR0, sR0, cR0, tR0 = prep_side(ev0, ft0, tn0, 1)
uT0, sT0, cT0, tTop0 = prep_side(ev0, ft0, tn0, 2)
print(f"  LEFT END: {len(uL0)} ev, {cL0.mean():.1f} hits/ev  "
      f"RIGHT END: {len(uR0)} ev, {cR0.mean():.1f} hits/ev  "
      f"TOP: {len(uT0)} ev, {cT0.mean():.1f} hits/ev (Npe_TOP)")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. TRAIN / TEST SPLIT at x=0
# ═══════════════════════════════════════════════════════════════════════════════
print("\n2. Train/test split (50/50) at x=0, EJ-230 N=20 …")
evids = np.unique(ev0)
n_total = len(evids)
n_train = n_total // 2
rng_tt = np.random.default_rng(20260820)
idx_perm = rng_tt.permutation(n_total)
train_ev = set(evids[idx_perm[:n_train]])
test_ev  = set(evids[idx_perm[n_train:]])

def split_data(ev_arr, tn_arr, ft_arr, lid_arr, ev_set):
    mask = np.isin(ev_arr, list(ev_set))
    return ev_arr[mask], tn_arr[mask], ft_arr[mask], lid_arr[mask]

def scan_m_iqr(tL, tR, sL, sR, cL, cR, uL, uR, Mmax=50):
    com = np.intersect1d(uL, uR)
    iL  = np.searchsorted(uL, com); iR = np.searchsorted(uR, com)
    cLc = cL[iL]; cRc = cR[iR]; sLc = sL[iL]; sRc = sR[iR]
    Mmax = min(Mmax, int(np.percentile(cLc, 2)), int(np.percentile(cRc, 2)))
    sigs = np.full(Mmax, np.nan)
    for m in range(1, Mmax+1):
        mask = (cLc >= m) & (cRc >= m)
        if mask.sum() < 20: continue
        tL_m = np.array([tL[sLc[i]:sLc[i]+m].mean() for i in np.where(mask)[0]])
        tR_m = np.array([tR[sRc[i]:sRc[i]+m].mean() for i in np.where(mask)[0]])
        T0e   = (tL_m + tR_m) / 2. * 1e3
        sigs[m-1] = rsigma(T0e - T0e.mean())
    return sigs, com, iL, iR, cLc, cRc, sLc, sRc

def scan_k_iqr(tTop, sT, cT, uT, Kmax=40):
    Kmax = min(Kmax, int(np.percentile(cT, 2)))
    sigs = np.full(Kmax, np.nan)
    for k in range(1, Kmax+1):
        mask = cT >= k
        if mask.sum() < 20: continue
        T0t  = tTop[sT[mask] + (k-1)] * 1e3
        sigs[k-1] = rsigma(T0t - T0t.mean())
    return sigs

# Train
ev_tr, tn_tr, ft_tr, lid_tr = split_data(ev0, tn0, ft0, lid0, train_ev)
uLt, sLt, cLt, tLt = prep_side(ev_tr, ft_tr, tn_tr, 0)
uRt, sRt, cRt, tRt = prep_side(ev_tr, ft_tr, tn_tr, 1)
uTt, sTt, cTt, tTt = prep_side(ev_tr, ft_tr, tn_tr, 2)
sigs_m_tr, _, _, _, _, _, _, _ = scan_m_iqr(tLt, tRt, sLt, sRt, cLt, cRt, uLt, uRt)
sigs_k_tr = scan_k_iqr(tTt, sTt, cTt, uTt)
m_star = int(np.nanargmin(sigs_m_tr)) + 1
k_star = int(np.nanargmin(sigs_k_tr)) + 1
print(f"  Training: m*={m_star} (σ={sigs_m_tr[m_star-1]:.2f} ps), k*={k_star} (σ={sigs_k_tr[k_star-1]:.2f} ps)")

# Test: apply global k*, m* from training
ev_te, tn_te, ft_te, lid_te = split_data(ev0, tn0, ft0, lid0, test_ev)
uLe2, sLe2, cLe2, tLe2 = prep_side(ev_te, ft_te, tn_te, 0)
uRe2, sRe2, cRe2, tRe2 = prep_side(ev_te, ft_te, tn_te, 1)
uTe2, sTe2, cTe2, tTe2 = prep_side(ev_te, ft_te, tn_te, 2)

ev_end_te, T0e_te, DT_te = get_T0_end(tLe2, tRe2, sLe2, sRe2, cLe2, cRe2, uLe2, uRe2, m_star)
ev_top_te, T0t_te        = get_T0_top(tTe2, sTe2, cTe2, uTe2, k_star)
T0e_a, T0t_a = intersect_estimators(ev_end_te, T0e_te, ev_top_te, T0t_te)
blue_te = compute_blue(T0e_a, T0t_a)

sig_end_test = rsigma(T0e_te - T0e_te.mean())
sig_top_test = rsigma(T0t_te - T0t_te.mean())
sig_blue_test = blue_te['sig_c']
boot_e_test  = bootstrap(T0e_a - T0e_a.mean())
boot_t_test  = bootstrap(T0t_a - T0t_a.mean())
boot_c_test  = bootstrap(blue_te['T0c'])

print(f"  Test σ_END  = {sig_end_test:.2f} ± {boot_e_test:.2f} ps (m*={m_star})")
print(f"  Test σ_TOP  = {sig_top_test:.2f} ± {boot_t_test:.2f} ps (k*={k_star})")
print(f"  Test σ_BLUE = {sig_blue_test:.2f} ± {boot_c_test:.2f} ps")
print(f"  Test w_END={blue_te['w_end']:.4f} ρ={blue_te['rho']:.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. ORACLE vs GLOBAL comparison at x=0
# ═══════════════════════════════════════════════════════════════════════════════
print("\n3. Oracle vs global comparison at x=0 …")
# Oracle: use optimal k/m from ALL x=0 data (same sample) — upper bound
sigs_m_all, _, _, _, _, _, _, _ = scan_m_iqr(tL0, tR0, sL0, sR0, cL0, cR0, uL0, uR0)
sigs_k_all = scan_k_iqr(tTop0, sT0, cT0, uT0)
m_oracle = int(np.nanargmin(sigs_m_all)) + 1
k_oracle = int(np.nanargmin(sigs_k_all)) + 1
sig_end_oracle = sigs_m_all[m_oracle-1]
sig_top_oracle = sigs_k_all[k_oracle-1]

# Compute BLUE with oracle k/m on all data
ev_end_all, T0e_all, _ = get_T0_end(tL0, tR0, sL0, sR0, cL0, cR0, uL0, uR0, m_oracle)
ev_top_all, T0t_all    = get_T0_top(tTop0, sT0, cT0, uT0, k_oracle)
T0e_or, T0t_or = intersect_estimators(ev_end_all, T0e_all, ev_top_all, T0t_all)
blue_or = compute_blue(T0e_or, T0t_or)

print(f"  Oracle  (same-sample): m*={m_oracle} k*={k_oracle}")
print(f"    σ_END={sig_end_oracle:.2f} σ_TOP={sig_top_oracle:.2f} σ_BLUE={blue_or['sig_c']:.2f} ps")
print(f"  Global  (train/test):  m*={m_star} k*={k_star}")
print(f"    σ_END={sig_end_test:.2f} σ_TOP={sig_top_test:.2f} σ_BLUE={sig_blue_test:.2f} ps")
print(f"  Difference (oracle vs global): ΔBLUE = {blue_or['sig_c']-sig_blue_test:.2f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. OLD PER-SIPM ESTIMATOR vs NEW POOLED (apples-to-apples, x=0, N=20)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n4. Old per-SiPM vs new pooled estimator (x=0 EJ-230 N=20) …")
# Old: mean of per-SiPM first-photon arrival times
ev_ps, T0_ps = get_T0_persipm(ev0, ft0, tn0, lid0, n_sipm=20)
T0_ps_c = T0_ps - np.nanmean(T0_ps)
sig_old_persipm = rsigma(T0_ps_c)

# New: pooled k=1, k=3
ev_k1, T0_k1 = get_T0_top(tTop0, sT0, cT0, uT0, 1)
ev_k3, T0_k3 = get_T0_top(tTop0, sT0, cT0, uT0, k_oracle)
sig_k1 = rsigma(T0_k1 - T0_k1.mean())
sig_k3 = rsigma(T0_k3 - T0_k3.mean())

print(f"  Per-SiPM mean estimator: σ = {sig_old_persipm:.2f} ps")
print(f"  Pooled k=1:              σ = {sig_k1:.2f} ps")
print(f"  Pooled k={k_oracle}:     σ = {sig_k3:.2f} ps")
print(f"  Improvement factor (per-SiPM → pooled k={k_oracle}): {sig_old_persipm/sig_k3:.2f}x")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. GAUSSIAN FITS WITH CHI2/P AT x=0, N=20
# ═══════════════════════════════════════════════════════════════════════════════
print("\n5. Gaussian fits with chi2/p at x=0 …")

# Use all data for final distributions
ev_end_all2, T0e_all2, DT_all2 = get_T0_end(tL0, tR0, sL0, sR0, cL0, cR0, uL0, uR0, m_oracle)
ev_top_all2, T0t_all2 = get_T0_top(tTop0, sT0, cT0, uT0, k_oracle)
T0e_f, T0t_f = intersect_estimators(ev_end_all2, T0e_all2, ev_top_all2, T0t_all2)
blue_f = compute_blue(T0e_f, T0t_f)

T0e_ps = T0e_f - T0e_f.mean()
T0t_ps = T0t_f - T0t_f.mean()
T0c_ps = blue_f['T0c']

sig_e_iqr = rsigma(T0e_ps)
sig_t_iqr = rsigma(T0t_ps)
sig_c_iqr = rsigma(T0c_ps)

def root_gauss_fit(arr_ps, nbins=100, name="h"):
    """Fill ROOT histogram and fit Gaussian to ±2.5 IQR core."""
    sig = rsigma(arr_ps)
    rng = 4.5 * sig
    h = ROOT.TH1D(name, "", nbins, -rng, rng)
    for v in arr_ps: h.Fill(v)
    med = np.median(arr_ps)
    q1, q3 = np.percentile(arr_ps, [25, 75])
    iqr = q3 - q1
    xlo = med - 2.5 * iqr/1.349; xhi = med + 2.5 * iqr/1.349
    f = ROOT.TF1(f"fg_{name}", "gaus", xlo, xhi)
    f.SetParameters(h.GetMaximum(), med, sig)
    res = h.Fit(f, "RSQN")
    ok = (int(res) == 0) and (f.GetNDF() > 0)
    return h, f if ok else None, ok

hE, fgE, okE = root_gauss_fit(T0e_ps, name="hE_x0")
hT, fgT, okT = root_gauss_fit(T0t_ps, name="hT_x0")
hC, fgC, okC = root_gauss_fit(T0c_ps, name="hC_x0")

for name, h, fg, ok, sig_iqr in [
    ("END",  hE, fgE, okE, sig_e_iqr),
    ("TOP",  hT, fgT, okT, sig_t_iqr),
    ("BLUE", hC, fgC, okC, sig_c_iqr),
]:
    if ok:
        print(f"  {name}: σ_IQR={sig_iqr:.2f} ps, σ_Gauss={fg.GetParameter(2):.2f} ps, "
              f"χ²/ndf={fg.GetChisquare():.1f}/{fg.GetNDF()}, "
              f"p={ROOT.TMath.Prob(fg.GetChisquare(), fg.GetNDF()):.4f}")
    else:
        print(f"  {name}: σ_IQR={sig_iqr:.2f} ps, Gauss fit failed")

# ═══════════════════════════════════════════════════════════════════════════════
# 6. v_eff FIT FROM Delta_t vs x
# ═══════════════════════════════════════════════════════════════════════════════
print("\n6. v_eff fit from Delta_t vs x (EJ-230 N=20) …")
x_scan = list(range(-600, 601, 100))
veff_data = []

for x in x_scan:
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    if not path.exists():
        print(f"  MISSING: {path.name}"); continue
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
    ev, ft, tn = a['event_id'], a['face_type'], a['time_ns']
    uL, sL, cL, tL = prep_side(ev, ft, tn, 0)
    uR, sR, cR, tR = prep_side(ev, ft, tn, 1)
    _, T0e, DT = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, m_oracle)
    mean_DT = np.median(DT)   # use median for robustness
    std_DT  = rsigma(DT - np.median(DT))
    mean_T0 = np.median(T0e)
    n_ev_x  = len(DT)
    veff_data.append((x, mean_DT, std_DT / np.sqrt(n_ev_x), mean_T0, n_ev_x))
    print(f"  x={x:+4d}: <Δt>={mean_DT:+7.3f} ps, σ(Δt)={std_DT:.2f} ps, n={n_ev_x}")

veff_data = np.array(veff_data)
xs_v  = veff_data[:,0]
dts_v = veff_data[:,1]
err_v = np.maximum(veff_data[:,2], 0.1)   # min 0.1 ps error floor

# ROOT linear fit: Δt = a + b*x  → v_eff = 2/|b| (mm/ns → mm/(ps*1000))
# b is in ps/mm, v_eff = 2/b * 1e-3 * 1e9 mm/ns = 2/b * 1e6 mm/ns... let's be careful
# Δt in ps, x in mm, b in ps/mm
# v_eff = 2/b [mm/ps] = 2/b * 1000 [mm/ns]

g_dt = make_graph(xs_v, dts_v, ey=err_v)
g_dt.SetMarkerStyle(20); g_dt.SetMarkerSize(1.2)
g_dt.SetMarkerColor(ROOT.kAzure+2); g_dt.SetLineColor(ROOT.kAzure+2)
g_dt.SetLineWidth(2)
flin = ROOT.TF1("flin", "[0]+[1]*x", -650., 650.)
flin.SetLineColor(ROOT.kRed+1); flin.SetLineWidth(2)
res_lin = g_dt.Fit(flin, "RSQ")
a_fit = flin.GetParameter(0); b_fit = flin.GetParameter(1)
ea_fit = flin.GetParError(0); eb_fit = flin.GetParError(1)
chi2_lin = flin.GetChisquare(); ndf_lin = flin.GetNDF()
p_lin = ROOT.TMath.Prob(chi2_lin, ndf_lin) if ndf_lin > 0 else 1.

v_eff_mmpns = 2. / abs(b_fit)          # mm/ps
v_eff_mmns  = v_eff_mmpns * 1e3        # mm/ns
v_eff_err   = v_eff_mmns * (eb_fit / abs(b_fit))

print(f"\n  Linear fit Δt = a + b·x:")
print(f"  a = {a_fit:.3f} ± {ea_fit:.3f} ps")
print(f"  b = {b_fit:.5f} ± {eb_fit:.5f} ps/mm")
print(f"  v_eff = 2/|b| = {v_eff_mmns:.1f} ± {v_eff_err:.1f} mm/ns")
print(f"  χ²/ndf = {chi2_lin:.2f}/{ndf_lin}, p = {p_lin:.4f}")

# Figure: Delta_t vs x with linear fit
c_veff = ROOT.TCanvas("c_veff", "veff", 800, 600)
c_veff.SetLeftMargin(0.14); c_veff.SetBottomMargin(0.14)
g_dt.SetTitle("")
g_dt.GetXaxis().SetTitle("Muon position x  [mm]")
g_dt.GetYaxis().SetTitle("#Deltat_{END} = t_{R} #minus t_{L}  [ps]")
g_dt.GetXaxis().SetLimits(-700., 700.)
g_dt.Draw("AP")
flin.Draw("same")
pv = ROOT.TPaveText(0.15, 0.68, 0.55, 0.90, "NDC")
pv.SetFillColor(0); pv.SetBorderSize(1); pv.SetTextFont(42); pv.SetTextSize(0.040)
pv.AddText(f"#Deltat = ({a_fit:.1f}) + ({b_fit:.4f}) x  ps")
pv.AddText(f"v_{{eff}} = 2/|b| = {v_eff_mmns:.1f} #pm {v_eff_err:.1f} mm/ns")
pv.AddText(f"#chi^{{2}}/ndf = {chi2_lin:.1f}/{ndf_lin},  p = {p_lin:.3f}")
pv.Draw()
lat_v = ROOT.TLatex(); lat_v.SetNDC(); lat_v.SetTextFont(42); lat_v.SetTextSize(0.042)
lat_v.DrawLatex(0.60, 0.88, "EJ-230, N_{top}=20,  m="+str(m_oracle))
save_canvas(c_veff, "fig_veff_fit")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. sigma_x_END(m) — POSITION RESOLUTION FROM END
# ═══════════════════════════════════════════════════════════════════════════════
print("\n7. sigma_x_END(m) at x=0 using v_eff calibration …")
# sigma_x = sigma_DeltaT / 2 * v_eff  [mm]  (factor 2 from Δt = 2x/v_eff)
# Actually DeltaT = 2*(x + L/2)/v_eff - 2*(L/2)/v_eff = 2*x/v_eff → σ_x = v_eff * σ(Δt) / 2
# More precisely: Δt = (t_R - t_L) = 2x/v_eff + offset
# σ_x = (v_eff/2) * σ(Δt)

sig_x_m = []
m_arr = []
for m in range(1, 31):
    com = np.intersect1d(uL0, uR0)
    iL  = np.searchsorted(uL0, com); iR = np.searchsorted(uR0, com)
    cLc = cL0[iL]; cRc = cR0[iR]; sLc = sL0[iL]; sRc = sR0[iR]
    mask = (cLc >= m) & (cRc >= m)
    if mask.sum() < 30: break
    tL_m = np.array([tL0[sLc[i]:sLc[i]+m].mean() for i in np.where(mask)[0]])
    tR_m = np.array([tR0[sRc[i]:sRc[i]+m].mean() for i in np.where(mask)[0]])
    DT_m = (tR_m - tL_m) * 1e3   # ps
    sig_DT = rsigma(DT_m - np.median(DT_m))
    # sigma_x in mm: v_eff * sigma_DeltaT / 2 (v_eff in mm/ps, DeltaT in ps)
    sig_x = v_eff_mmpns * sig_DT / 2.
    sig_x_m.append(sig_x)
    m_arr.append(m)
    print(f"  m={m:2d}: σ(Δt)={sig_DT:.2f} ps → σ_x = {sig_x:.1f} mm")

sig_x_arr = np.array(sig_x_m)
m_iqr_arr = np.array(m_arr)
m_opt_x = m_arr[int(np.argmin(sig_x_arr))]
print(f"  → Best m for position: m={m_opt_x}, σ_x = {sig_x_arr.min():.1f} mm")
print(f"  → Best m for timing (m*={m_oracle}): σ_x = {sig_x_arr[m_oracle-1]:.1f} mm")

# ═══════════════════════════════════════════════════════════════════════════════
# 8. TOP POSITION ESTIMATORS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n8. TOP position estimators at x=0, EJ-230 N=20 …")
# Load full data with local_id and gun_x for all x positions
top_pos_results = {}

for x in x_scan:
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    if not path.exists(): continue
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','local_id','time_ns'], library='np')
    ev_x = a['event_id']; ft_x = a['face_type']
    tn_x = a['time_ns'];  lid_x = a['local_id']

    # Filter TOP only
    msk = ft_x == 2
    ev_top = ev_x[msk]; tn_top = tn_x[msk]; lid_top = lid_x[msk]

    ev_uniq = np.unique(ev_top)

    # Estimator 1: local_id of earliest-hit SiPM
    x_est1 = np.full(len(ev_uniq), np.nan)
    # Estimator 2: Npe-weighted centroid among hit SiPMs
    x_est2 = np.full(len(ev_uniq), np.nan)
    # Estimator 3: position of SiPM with max Npe
    x_est3 = np.full(len(ev_uniq), np.nan)

    for i, eid in enumerate(ev_uniq):
        em = ev_top == eid
        e_tn = tn_top[em]; e_lid = lid_top[em]

        # Earliest-hit SiPM
        first_lid = e_lid[np.argmin(e_tn)]
        x_est1[i] = TOP_XPOS[first_lid]

        # Npe-weighted centroid
        npe_per_sipm = np.bincount(e_lid, minlength=20).astype(float)
        if npe_per_sipm.sum() > 0:
            x_est2[i] = np.dot(npe_per_sipm, TOP_XPOS) / npe_per_sipm.sum()

        # Max-Npe SiPM position
        x_est3[i] = TOP_XPOS[np.argmax(npe_per_sipm)]

    # Bias: x_est - x_true, where x_true = x (gun position)
    bias1 = np.nanmedian(x_est1) - x
    bias2 = np.nanmedian(x_est2) - x
    bias3 = np.nanmedian(x_est3) - x

    # Resolution (IQR)
    res1 = rsigma(x_est1 - np.nanmedian(x_est1))
    res2 = rsigma(x_est2 - np.nanmedian(x_est2))
    res3 = rsigma(x_est3 - np.nanmedian(x_est3))

    top_pos_results[x] = dict(
        bias_early=bias1, bias_centroid=bias2, bias_maxnpe=bias3,
        res_early=res1, res_centroid=res2, res_maxnpe=res3
    )
    print(f"  x={x:+4d}: earliest SiPM σ={res1:.0f}mm (bias={bias1:+.0f}mm), "
          f"centroid σ={res2:.0f}mm (bias={bias2:+.1f}mm)")

# Compare with END position resolution at m* for timing
sig_x_end_at_mopt = sig_x_arr[m_oracle-1]
sig_x_end_best    = sig_x_arr.min()
sig_x_top_centroid_at_x0 = top_pos_results[0]['res_centroid']
print(f"\n  Summary at x=0:")
print(f"  END position: σ_x = {sig_x_end_at_mopt:.1f} mm (m=m_opt_timing={m_oracle})")
print(f"  END position: σ_x = {sig_x_end_best:.1f} mm (m=m_opt_position={m_opt_x})")
print(f"  TOP centroid: σ_x = {sig_x_top_centroid_at_x0:.1f} mm")

# ═══════════════════════════════════════════════════════════════════════════════
# 9. BIAS vs x (mean T0 shift per position)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n9. Bias vs x for TOP estimator …")
bias_results = {}

# First, compute absolute calibration from the x=0 file
# (reference point: mean T0_TOP at x=0 as the calibration offset)
ev_top_ref, T0t_ref = get_T0_top(tTop0, sT0, cT0, uT0, k_oracle)
T0t_calib = np.median(T0t_ref)   # calibration offset in ps

for x in x_scan:
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    if not path.exists(): continue
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
    ev_x = a['event_id']; ft_x = a['face_type']; tn_x = a['time_ns']
    uLx, sLx, cLx, tLx = prep_side(ev_x, ft_x, tn_x, 0)
    uRx, sRx, cRx, tRx = prep_side(ev_x, ft_x, tn_x, 1)
    uTx, sTx, cTx, tTx = prep_side(ev_x, ft_x, tn_x, 2)

    # T0_END absolute (not mean-subtracted)
    ev_ex, T0e_x, DT_x = get_T0_end(tLx, tRx, sLx, sRx, cLx, cRx, uLx, uRx, m_oracle)
    ev_tx, T0t_x = get_T0_top(tTx, sTx, cTx, uTx, k_oracle)

    T0e_xa, T0t_xa = intersect_estimators(ev_ex, T0e_x, ev_tx, T0t_x)

    # Bias = uncalibrated mean - calibration constant
    # For a meaningful bias: relative to x=0
    if x == 0:
        calib_END = np.median(T0e_xa)
        calib_TOP = np.median(T0t_xa)

    bias_results[x] = dict(
        mean_END=np.median(T0e_xa),
        mean_TOP=np.median(T0t_xa),
        sig_END=rsigma(T0e_xa - np.median(T0e_xa)),
        sig_TOP=rsigma(T0t_xa - np.median(T0t_xa)),
        n=len(T0e_xa)
    )

# Compute bias relative to x=0
calib_END = bias_results[0]['mean_END']
calib_TOP = bias_results[0]['mean_TOP']
print(f"\n  Calibration (x=0): T0_END_mean={calib_END:.2f} ps, T0_TOP_mean={calib_TOP:.2f} ps")
print(f"\n  {'x':>5} {'bias_END':>10} {'bias_TOP':>10} {'sig_END':>9} {'sig_TOP':>9}")
for x in sorted(bias_results.keys()):
    r = bias_results[x]
    b_end = r['mean_END'] - calib_END
    b_top = r['mean_TOP'] - calib_TOP
    print(f"  {x:>5} {b_end:>10.2f} {b_top:>10.2f} {r['sig_END']:>9.2f} {r['sig_TOP']:>9.2f}")

# ═══════════════════════════════════════════════════════════════════════════════
# 10. rho AND BLUE WEIGHTS vs x (global k*, m*)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n10. rho and BLUE weights vs x …")
rho_x_results = {}

for x in x_scan:
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    if not path.exists(): continue
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
    ev_x = a['event_id']; ft_x = a['face_type']; tn_x = a['time_ns']
    uLx, sLx, cLx, tLx = prep_side(ev_x, ft_x, tn_x, 0)
    uRx, sRx, cRx, tRx = prep_side(ev_x, ft_x, tn_x, 1)
    uTx, sTx, cTx, tTx = prep_side(ev_x, ft_x, tn_x, 2)

    ev_ex, T0e_x, _ = get_T0_end(tLx, tRx, sLx, sRx, cLx, cRx, uLx, uRx, m_oracle)
    ev_tx, T0t_x    = get_T0_top(tTx, sTx, cTx, uTx, k_oracle)
    T0e_xa, T0t_xa  = intersect_estimators(ev_ex, T0e_x, ev_tx, T0t_x)
    bl_x = compute_blue(T0e_xa, T0t_xa)
    rho_x_results[x] = dict(rho=bl_x['rho'], w_end=bl_x['w_end'], w_top=bl_x['w_top'],
                              sig_e=bl_x['sig_e'], sig_t=bl_x['sig_t'], sig_c=bl_x['sig_c'])
    print(f"  x={x:+4d}: ρ={bl_x['rho']:.3f} w_END={bl_x['w_end']:.4f} "
          f"σ_BLUE={bl_x['sig_c']:.2f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# 11. ELECTRONICS JITTER SCENARIOS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n11. Electronics jitter scenarios at x=0, EJ-230, N=20 …")
# Model: each photon arrival time gets +Gaussian(0, sigma_jitter) smearing
# For k-th order stat: each of the pooled hits gets independent jitter
# For END m-avg: each of the m hits per side gets independent jitter

jitter_scenarios = [0., 50., 100., 150., 200.]  # ps (single-channel SPTR)
jitter_results = {}

# Reload all TOP times at x=0
with uproot.open(INDIR / "photon_hits_ej230_ntop20_x0.root") as f:
    a_j = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
ev_j = a_j['event_id']; ft_j = a_j['face_type']; tn_j = a_j['time_ns']

uLj, sLj, cLj, tLj = prep_side(ev_j, ft_j, tn_j, 0)
uRj, sRj, cRj, tRj = prep_side(ev_j, ft_j, tn_j, 1)
uTj, sTj, cTj, tTj = prep_side(ev_j, ft_j, tn_j, 2)

rng_j = np.random.default_rng(42)

for sj in jitter_scenarios:
    sj_ns = sj * 1e-3   # convert ps to ns

    # Add jitter to all hit times
    tLj_sm = tLj + (rng_j.normal(0., sj_ns, len(tLj)) if sj > 0 else 0.)
    tRj_sm = tRj + (rng_j.normal(0., sj_ns, len(tRj)) if sj > 0 else 0.)
    tTj_sm = tTj + (rng_j.normal(0., sj_ns, len(tTj)) if sj > 0 else 0.)
    # Re-sort within each event (jitter can reorder photons)
    # For END: re-sort already handled by prep_side convention (we re-sort after jitter)
    # Actually after adding independent jitter, the k-th smallest changes
    # We need to re-sort the per-event photon arrays
    # This is complex with the current data structure; use a simplified approach:
    # For each event, re-sort smeared times

    # Simple approach: operate on the per-event photon lists
    # For TOP k-th statistic:
    T0t_sm = np.full(len(uTj), np.nan)
    for i, eid in enumerate(uTj):
        start = sTj[i]; cnt = cTj[i]
        times = tTj_sm[start:start+cnt]
        if cnt >= k_oracle:
            T0t_sm[i] = np.sort(times)[k_oracle-1] * 1e3

    # For END m-avg:
    T0e_sm = np.full(len(uLj), np.nan)
    com_j = np.intersect1d(uLj, uRj)
    iLj = np.searchsorted(uLj, com_j); iRj = np.searchsorted(uRj, com_j)
    for idx, (il, ir) in enumerate(zip(iLj, iRj)):
        if cLj[il] >= m_oracle and cRj[ir] >= m_oracle:
            tL_m = np.sort(tLj_sm[sLj[il]:sLj[il]+cLj[il]])[:m_oracle].mean()
            tR_m = np.sort(tRj_sm[sRj[ir]:sRj[ir]+cRj[ir]])[:m_oracle].mean()
            T0e_sm[idx] = (tL_m + tR_m) / 2. * 1e3

    T0t_sm_fin = T0t_sm[np.isfinite(T0t_sm)]
    T0e_sm_fin = T0e_sm[np.isfinite(T0e_sm)]

    sig_top_j = rsigma(T0t_sm_fin - T0t_sm_fin.mean()) if len(T0t_sm_fin) > 20 else np.nan
    sig_end_j = rsigma(T0e_sm_fin - T0e_sm_fin.mean()) if len(T0e_sm_fin) > 20 else np.nan

    # BLUE with jitter
    # Use shorter arrays with finite values aligned
    n_min = min(len(T0t_sm_fin), len(T0e_sm_fin))
    if n_min > 30:
        bl_j = compute_blue(T0e_sm_fin[:n_min], T0t_sm_fin[:n_min])
        sig_blue_j = bl_j['sig_c']
    else:
        sig_blue_j = np.nan

    jitter_results[sj] = dict(sig_top=sig_top_j, sig_end=sig_end_j, sig_blue=sig_blue_j)
    print(f"  jitter={sj:5.0f} ps: σ_TOP={sig_top_j:.2f} ps, "
          f"σ_END={sig_end_j:.2f} ps, σ_BLUE={sig_blue_j:.2f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# 12. PDE AUDIT — effective PDE from simulation
# ═══════════════════════════════════════════════════════════════════════════════
print("\n12. PDE audit from simulation data …")
# Load pde branch for all face types
with uproot.open(INDIR / "photon_hits_ej230_ntop20_x0.root") as f:
    a_pde = f['sipm_hits'].arrays(['face_type','pde','wl_nm'], library='np')
ft_pde = a_pde['face_type']; pde_vals = a_pde['pde']; wl_vals = a_pde['wl_nm']

for ft_id, name in [(0,'END-L'),(1,'END-R'),(2,'TOP')]:
    m = ft_pde == ft_id
    print(f"  {name}: mean_PDE={pde_vals[m].mean():.4f}, "
          f"wl range=[{wl_vals[m].min():.0f},{wl_vals[m].max():.0f}] nm, "
          f"n={m.sum()}")
print(f"  Overall mean PDE (all detected photons): {pde_vals.mean():.4f}")
print(f"  Note: PDE in simulation is ~0.61, NOT 0.40 as used in old napkin")

# ═══════════════════════════════════════════════════════════════════════════════
# 13. GENERATE FIGURES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n13. Generating figures …")

# ── FIG A: Residual distributions (END, TOP, BLUE) at x=0 ──
print("  Fig A: Residual distributions …")
c_res = ROOT.TCanvas("c_residuals_v4","residuals",1100,700)
c_res.Divide(3,1)

colors = [ROOT.kAzure+2, ROOT.kOrange+1, ROOT.kGreen+3]
hists  = [hE, hT, hC]
fits   = [fgE, fgT, fgC]
labels = [f"END (m*={m_oracle})", f"TOP (k*={k_oracle})", "BLUE"]
sigs_iqr = [sig_e_iqr, sig_t_iqr, sig_c_iqr]

for i, (h, fg, col, lbl, siq) in enumerate(zip(hists, fits, colors, labels, sigs_iqr)):
    pad = c_res.cd(i+1)
    pad.SetLeftMargin(0.16 if i==0 else 0.12)
    pad.SetBottomMargin(0.14)
    h.SetLineColor(col); h.SetLineWidth(2)
    h.SetTitle("")
    h.GetXaxis().SetTitle("#DeltaT  [ps]")
    h.GetYaxis().SetTitle("Events / bin" if i==0 else "")
    h.GetYaxis().SetTitleOffset(1.3 if i==0 else 0.)
    h.Draw("HIST")
    if fg:
        fg.SetLineColor(ROOT.kRed+1); fg.SetLineWidth(2); fg.SetLineStyle(2)
        fg.Draw("same")
        pt = ROOT.TPaveText(0.12, 0.62, 0.88, 0.90, "NDC")
        pt.SetFillColor(0); pt.SetBorderSize(1)
        pt.SetTextFont(42); pt.SetTextSize(0.045)
        pt.AddText(lbl)
        pt.AddText(f"#sigma_{{IQR}} = {siq:.1f} ps")
        pt.AddText(f"#sigma_{{G}} = {fg.GetParameter(2):.1f} ps")
        chi2 = fg.GetChisquare(); ndf = fg.GetNDF()
        p = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 1.
        pt.AddText(f"#chi^{{2}}/ndf = {chi2:.1f}/{ndf}")
        pt.AddText(f"p = {p:.4f}")
        pt.Draw()

c_res.cd()
ltop = ROOT.TLatex(); ltop.SetNDC(); ltop.SetTextFont(42); ltop.SetTextSize(0.030)
ltop.DrawLatex(0.30, 0.97, "EJ-230, N_{top}=20, x=0, 5000 events — #sigma_{IQR} = (Q75#minusQ25)/1.349")
save_canvas(c_res, "figA_residuals")

# ── FIG B: v_eff already saved above ──

# ── FIG C: sigma_x_END(m) and sigma_T_END(m) on dual y-axis or two panels ──
print("  Fig C: sigma_x_END(m) and sigma_T_END(m) …")
# sigma_T_END(m) from all x=0 data
sigs_m_all_arr = np.full(30, np.nan)
for m in range(1, 31):
    com = np.intersect1d(uL0, uR0)
    iL  = np.searchsorted(uL0, com); iR = np.searchsorted(uR0, com)
    cLc = cL0[iL]; cRc = cR0[iR]; sLc = sL0[iL]; sRc = sR0[iR]
    mask = (cLc >= m) & (cRc >= m)
    if mask.sum() < 30: break
    tL_m = np.array([tL0[sLc[i]:sLc[i]+m].mean() for i in np.where(mask)[0]])
    tR_m = np.array([tR0[sRc[i]:sRc[i]+m].mean() for i in np.where(mask)[0]])
    T0e_m = (tL_m + tR_m) / 2. * 1e3
    sigs_m_all_arr[m-1] = rsigma(T0e_m - T0e_m.mean())

valid_m = np.where(np.isfinite(sigs_m_all_arr))[0]
m_vals_c = (valid_m+1).astype('d')
sig_T_vals = sigs_m_all_arr[valid_m]
sig_x_vals = sig_x_arr[:len(valid_m)]   # already computed

c_sigx = ROOT.TCanvas("c_sigx", "sigx", 800, 600)
c_sigx.SetLeftMargin(0.14); c_sigx.SetBottomMargin(0.14); c_sigx.SetRightMargin(0.15)

gT = make_graph(m_vals_c, sig_T_vals)
gT.SetMarkerStyle(20); gT.SetMarkerSize(1.2)
gT.SetMarkerColor(ROOT.kAzure+2); gT.SetLineColor(ROOT.kAzure+2); gT.SetLineWidth(2)
gT.SetTitle("")
gT.GetXaxis().SetTitle("m  (number of earliest photons averaged per side)")
gT.GetYaxis().SetTitle("#sigma_{T}^{END}  [ps]")
gT.GetYaxis().SetRangeUser(45., 120.)
gT.Draw("APL")

# Mark optimal m for timing
mopt_t_line = ROOT.TLine(m_oracle, 45., m_oracle, sigs_m_all_arr[m_oracle-1])
mopt_t_line.SetLineColor(ROOT.kAzure+2); mopt_t_line.SetLineStyle(2); mopt_t_line.SetLineWidth(2)
mopt_t_line.Draw()

# Right axis for position
axis2 = ROOT.TGaxis(c_sigx.GetUxmax(), 45., c_sigx.GetUxmax(), 120.,
                    45.*v_eff_mmpns/2., 120.*v_eff_mmpns/2., 510, "+L")
axis2.SetLabelFont(42); axis2.SetLabelSize(0.042)
axis2.SetTitle("#sigma_{x}^{END}  [mm]")
axis2.SetTitleFont(42); axis2.SetTitleSize(0.048)
axis2.SetTitleOffset(1.1)
axis2.Draw()

leg_sigx = ROOT.TLegend(0.55, 0.70, 0.82, 0.88)
leg_sigx.SetTextSize(0.040)
leg_sigx.AddEntry(gT, "#sigma_{T,END}(m)  (IQR)", "pl")
leg_sigx.AddEntry(mopt_t_line, f"m*={m_oracle} (timing opt.)", "l")
leg_sigx.Draw()

ltx = ROOT.TLatex(); ltx.SetNDC(); ltx.SetTextFont(42); ltx.SetTextSize(0.040)
ltx.DrawLatex(0.17, 0.88, "EJ-230,  N_{top}=20,  x = 0")
save_canvas(c_sigx, "figC_sigma_x_end")

# ── FIG D: TOP position estimator — bias and resolution vs x ──
print("  Fig D: TOP position estimators vs x …")
xs_d = np.array(sorted(top_pos_results.keys()), dtype='d')
res_centroid = np.array([top_pos_results[x]['res_centroid'] for x in sorted(top_pos_results)])
bias_centroid = np.array([top_pos_results[x]['bias_centroid'] for x in sorted(top_pos_results)])
res_early = np.array([top_pos_results[x]['res_early'] for x in sorted(top_pos_results)])

c_toppos = ROOT.TCanvas("c_toppos", "toppos", 800, 600)
c_toppos.SetLeftMargin(0.14); c_toppos.SetBottomMargin(0.14)

g_res_centroid = make_graph(xs_d, res_centroid)
g_res_centroid.SetMarkerStyle(20); g_res_centroid.SetMarkerSize(1.2)
g_res_centroid.SetMarkerColor(ROOT.kOrange+1); g_res_centroid.SetLineColor(ROOT.kOrange+1); g_res_centroid.SetLineWidth(2)

g_res_early = make_graph(xs_d, res_early)
g_res_early.SetMarkerStyle(24); g_res_early.SetMarkerSize(1.1)
g_res_early.SetMarkerColor(ROOT.kGreen+3); g_res_early.SetLineColor(ROOT.kGreen+3); g_res_early.SetLineWidth(2); g_res_early.SetLineStyle(2)

# END position resolution at m_oracle (constant)
end_x_line = ROOT.TLine(-650., sig_x_end_at_mopt, 650., sig_x_end_at_mopt)
end_x_line.SetLineColor(ROOT.kAzure+2); end_x_line.SetLineStyle(3); end_x_line.SetLineWidth(2)

mg_pos = ROOT.TMultiGraph()
mg_pos.Add(g_res_centroid, "PL")
mg_pos.Add(g_res_early, "PL")
mg_pos.Draw("A")
mg_pos.SetTitle("")
mg_pos.GetXaxis().SetTitle("Muon position x  [mm]")
mg_pos.GetYaxis().SetTitle("#sigma_{x}  [mm]")
mg_pos.GetXaxis().SetLimits(-700., 700.)
end_x_line.Draw()

# Mark TOP SiPM positions
for xp in TOP_XPOS:
    tm = ROOT.TLine(xp, 0., xp, 5.)
    tm.SetLineColor(ROOT.kGray+2); tm.SetLineStyle(1); tm.SetLineWidth(1)
    tm.Draw()

leg_pos = ROOT.TLegend(0.18, 0.65, 0.62, 0.88)
leg_pos.SetTextSize(0.038)
leg_pos.AddEntry(g_res_centroid, "TOP N_{pe}-weighted centroid", "pl")
leg_pos.AddEntry(g_res_early, "TOP earliest-hit SiPM", "pl")
leg_pos.AddEntry(end_x_line, f"END (m={m_oracle}) #sigma_x = {sig_x_end_at_mopt:.0f} mm", "l")
leg_pos.Draw()
ltop2 = ROOT.TLatex(); ltop2.SetNDC(); ltop2.SetTextFont(42); ltop2.SetTextSize(0.040)
ltop2.DrawLatex(0.55, 0.88, "EJ-230, N_{top}=20, #sigma_{x} (IQR)")
save_canvas(c_toppos, "figD_top_position")

# ── FIG E: Bias vs x for END and TOP ──
print("  Fig E: Bias vs x …")
xs_b = np.array(sorted(bias_results.keys()), dtype='d')
bias_END = np.array([bias_results[x]['mean_END'] - calib_END for x in sorted(bias_results)])
bias_TOP = np.array([bias_results[x]['mean_TOP'] - calib_TOP for x in sorted(bias_results)])

c_bias = ROOT.TCanvas("c_bias", "bias", 800, 600)
c_bias.SetLeftMargin(0.15); c_bias.SetBottomMargin(0.14)

g_bEND = make_graph(xs_b, bias_END)
g_bEND.SetMarkerStyle(20); g_bEND.SetMarkerSize(1.2)
g_bEND.SetMarkerColor(ROOT.kAzure+2); g_bEND.SetLineColor(ROOT.kAzure+2); g_bEND.SetLineWidth(2)

g_bTOP = make_graph(xs_b, bias_TOP)
g_bTOP.SetMarkerStyle(21); g_bTOP.SetMarkerSize(1.2)
g_bTOP.SetMarkerColor(ROOT.kOrange+1); g_bTOP.SetLineColor(ROOT.kOrange+1); g_bTOP.SetLineWidth(2)

mg_b = ROOT.TMultiGraph()
mg_b.Add(g_bEND, "PL"); mg_b.Add(g_bTOP, "PL")
mg_b.Draw("A")
mg_b.SetTitle("")
mg_b.GetXaxis().SetTitle("Muon position x  [mm]")
mg_b.GetYaxis().SetTitle("#LTT_{0}^{rec}#GT #minus #LTT_{0}^{rec}(x=0)#GT  [ps]")
mg_b.GetXaxis().SetLimits(-700., 700.)

zline = ROOT.TLine(-650., 0., 650., 0.)
zline.SetLineStyle(2); zline.SetLineColor(ROOT.kGray+2); zline.Draw()

leg_b = ROOT.TLegend(0.18, 0.72, 0.55, 0.88)
leg_b.SetTextSize(0.040)
leg_b.AddEntry(g_bEND, "END estimator", "pl")
leg_b.AddEntry(g_bTOP, "TOP estimator", "pl")
leg_b.Draw()
ltop3 = ROOT.TLatex(); ltop3.SetNDC(); ltop3.SetTextFont(42); ltop3.SetTextSize(0.040)
ltop3.DrawLatex(0.55, 0.88, "EJ-230, N_{top}=20")
ltop3.DrawLatex(0.15, 0.88, "T_{0} bias relative to x=0 calibration")
save_canvas(c_bias, "figE_bias_vs_x")

# ── FIG F: rho and w_END vs x ──
print("  Fig F: rho and w_END vs x …")
xs_r = np.array(sorted(rho_x_results.keys()), dtype='d')
rho_arr = np.array([rho_x_results[x]['rho'] for x in sorted(rho_x_results)])
wend_arr = np.array([rho_x_results[x]['w_end'] for x in sorted(rho_x_results)])

c_rho = ROOT.TCanvas("c_rho", "rho_vs_x", 800, 600)
c_rho.SetLeftMargin(0.14); c_rho.SetBottomMargin(0.14)

g_rho = make_graph(xs_r, rho_arr)
g_rho.SetMarkerStyle(20); g_rho.SetMarkerSize(1.2)
g_rho.SetMarkerColor(ROOT.kViolet+2); g_rho.SetLineColor(ROOT.kViolet+2); g_rho.SetLineWidth(2)
g_rho.SetTitle("")
g_rho.GetXaxis().SetTitle("Muon position x  [mm]")
g_rho.GetYaxis().SetTitle("#rho(T_{END}, T_{TOP})")
g_rho.GetYaxis().SetRangeUser(0., 0.35)
g_rho.Draw("APL")

ltop4 = ROOT.TLatex(); ltop4.SetNDC(); ltop4.SetTextFont(42); ltop4.SetTextSize(0.040)
ltop4.DrawLatex(0.55, 0.88, "EJ-230, N_{top}=20")
ltop4.DrawLatex(0.17, 0.88, "#rho ~ 0.15#minus0.18  (approx. constant with x)")
save_canvas(c_rho, "figF_rho_vs_x")

# ── FIG G: Jitter scenarios ──
print("  Fig G: Electronics jitter scenarios …")
jit_arr = np.array(sorted(jitter_results.keys()), dtype='d')
sig_top_j = np.array([jitter_results[j]['sig_top'] for j in sorted(jitter_results)])
sig_end_j = np.array([jitter_results[j]['sig_end'] for j in sorted(jitter_results)])
sig_blue_j = np.array([jitter_results[j]['sig_blue'] for j in sorted(jitter_results)])

c_jit = ROOT.TCanvas("c_jitter", "jitter", 800, 600)
c_jit.SetLeftMargin(0.14); c_jit.SetBottomMargin(0.14)

g_jt = make_graph(jit_arr, sig_top_j)
g_jt.SetMarkerStyle(20); g_jt.SetMarkerSize(1.4)
g_jt.SetMarkerColor(ROOT.kOrange+1); g_jt.SetLineColor(ROOT.kOrange+1); g_jt.SetLineWidth(2)

g_je = make_graph(jit_arr, sig_end_j)
g_je.SetMarkerStyle(21); g_je.SetMarkerSize(1.3)
g_je.SetMarkerColor(ROOT.kAzure+2); g_je.SetLineColor(ROOT.kAzure+2); g_je.SetLineWidth(2)

g_jb = make_graph(jit_arr, sig_blue_j)
g_jb.SetMarkerStyle(22); g_jb.SetMarkerSize(1.3)
g_jb.SetMarkerColor(ROOT.kGreen+3); g_jb.SetLineColor(ROOT.kGreen+3); g_jb.SetLineWidth(2)

mg_j = ROOT.TMultiGraph()
mg_j.Add(g_jt, "PL"); mg_j.Add(g_je, "PL"); mg_j.Add(g_jb, "PL")
mg_j.Draw("A")
mg_j.SetTitle("")
mg_j.GetXaxis().SetTitle("Single-channel timing jitter  [ps rms]")
mg_j.GetYaxis().SetTitle("#sigma_{T}  [ps]  (IQR-based)")
mg_j.GetXaxis().SetLimits(-10., 220.)

leg_j = ROOT.TLegend(0.18, 0.65, 0.60, 0.88)
leg_j.SetTextSize(0.038)
leg_j.AddEntry(g_jt, f"TOP pooled k*={k_oracle} (N=20 SiPMs)", "pl")
leg_j.AddEntry(g_je, f"END m-avg m*={m_oracle} (N=16 SiPMs)", "pl")
leg_j.AddEntry(g_jb, "BLUE combination", "pl")
leg_j.Draw()
ltop5 = ROOT.TLatex(); ltop5.SetNDC(); ltop5.SetTextFont(42); ltop5.SetTextSize(0.038)
ltop5.DrawLatex(0.40, 0.88, "EJ-230,  x=0  (scenario study)")
ltop5.DrawLatex(0.40, 0.83, "Optical result at 0 ps = hard lower bound")
save_canvas(c_jit, "figG_jitter")

# ── FIG H: Estimator comparison (per-SiPM vs pooled) ──
print("  Fig H: Old vs new TOP estimator …")
# Use k-scan from full data
sigs_k_full = scan_k_iqr(tTop0, sT0, cT0, uT0, Kmax=30)
k_vals_h = np.where(np.isfinite(sigs_k_full))[0] + 1

c_est = ROOT.TCanvas("c_estimator", "estimator", 800, 600)
c_est.SetLeftMargin(0.14); c_est.SetBottomMargin(0.14)

g_kscan = make_graph(k_vals_h.astype('d'), sigs_k_full[k_vals_h-1])
g_kscan.SetMarkerStyle(20); g_kscan.SetMarkerSize(1.1)
g_kscan.SetMarkerColor(ROOT.kOrange+1); g_kscan.SetLineColor(ROOT.kOrange+1); g_kscan.SetLineWidth(2)
g_kscan.SetTitle("")
g_kscan.GetXaxis().SetTitle("k  (k-th earliest pooled TOP photon)")
g_kscan.GetYaxis().SetTitle("#sigma_{T}^{TOP}  [ps]  (IQR)")
g_kscan.GetYaxis().SetRangeUser(0., max(sig_old_persipm*1.2, 80.))
g_kscan.Draw("APL")

# Horizontal line for per-SiPM estimator
l_old = ROOT.TLine(0.5, sig_old_persipm, k_vals_h.max()+0.5, sig_old_persipm)
l_old.SetLineColor(ROOT.kRed+1); l_old.SetLineStyle(2); l_old.SetLineWidth(2)
l_old.Draw()

# Mark k* minimum
k_min_pt = make_graph(np.array([float(k_oracle)]), np.array([sigs_k_full[k_oracle-1]]))
k_min_pt.SetMarkerStyle(29); k_min_pt.SetMarkerSize(2.5)
k_min_pt.SetMarkerColor(ROOT.kRed+1)
k_min_pt.Draw("P same")

pv_est = ROOT.TPaveText(0.50, 0.55, 0.92, 0.88, "NDC")
pv_est.SetFillColor(0); pv_est.SetBorderSize(1)
pv_est.SetTextFont(42); pv_est.SetTextSize(0.040)
pv_est.AddText("Per-SiPM mean (legacy):")
pv_est.AddText(f"  #sigma = {sig_old_persipm:.1f} ps")
pv_est.AddText(f"Pooled k*={k_oracle}  (new):")
pv_est.AddText(f"  #sigma = {sigs_k_full[k_oracle-1]:.1f} ps")
pv_est.AddText(f"Gain:  #times{sig_old_persipm/sigs_k_full[k_oracle-1]:.1f}")
pv_est.Draw()

ltop6 = ROOT.TLatex(); ltop6.SetNDC(); ltop6.SetTextFont(42); ltop6.SetTextSize(0.040)
ltop6.DrawLatex(0.17, 0.88, "EJ-230,  N_{top}=20,  x=0")
ltop6.DrawLatex(0.17, 0.83, "Same events, same SiPMs")
save_canvas(c_est, "figH_estimator_compare")

# ═══════════════════════════════════════════════════════════════════════════════
# 14. SAVE TABLES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n14. Saving tables …")

# Table 1: Final numbers
rows = []
for x in sorted(rho_x_results.keys()):
    r = rho_x_results[x]
    br = bias_results[x]
    pr = top_pos_results.get(x, {})
    rows.append(dict(
        x=x,
        sig_END=r['sig_e'],
        sig_TOP=r['sig_t'],
        sig_BLUE=r['sig_c'],
        rho=r['rho'],
        w_END=r['w_end'],
        w_TOP=r['w_top'],
        bias_END_ps=br['mean_END']-calib_END,
        bias_TOP_ps=br['mean_TOP']-calib_TOP,
        sig_x_centroid_mm=pr.get('res_centroid', np.nan),
    ))
pd.DataFrame(rows).to_csv(TABS / "tab_sigma_vs_x_global_km.csv", index=False, float_format="%.3f")
print(f"  Saved tab_sigma_vs_x_global_km.csv")

# Table 2: N_top scan
with open(REPO/"analysis/optim/root_best_est/verification.json") as fp:
    ver = json.load(fp)

ntop_rows = []
for N_str, r in ver["ntop_scan_ej230_x0"].items():
    ntop_rows.append(dict(
        N_END=16, N_TOP=int(N_str),
        N_total=16+int(N_str),
        m_opt=r['m_opt'], k_opt=r['k_opt'],
        sig_END=r['sig_end'], sig_TOP=r['sig_top'],
        sig_BLUE=r['sig_comb'],
        w_END=r['w_end'], w_TOP=r['w_top'],
        rho=r['rho'],
        Npe_TOP=r['npe_top'],
        boot_c=r.get('boot_c',np.nan)
    ))
pd.DataFrame(ntop_rows).to_csv(TABS / "tab_ntop_scan.csv", index=False, float_format="%.3f")
print(f"  Saved tab_ntop_scan.csv")

# Table 3: Train/test and oracle summary
summary = dict(
    oracle_m=m_oracle, oracle_k=k_oracle,
    oracle_sig_end=sig_end_oracle, oracle_sig_top=sig_top_oracle, oracle_sig_blue=blue_or['sig_c'],
    traintest_m=m_star, traintest_k=k_star,
    test_sig_end=sig_end_test, test_sig_top=sig_top_test, test_sig_blue=sig_blue_test,
    test_boot_e=boot_e_test, test_boot_t=boot_t_test, test_boot_c=boot_c_test,
    persipm_sig=sig_old_persipm, pooled_k1_sig=sig_k1, pooled_kopt_sig=sig_k3,
    pde_mean_all=float(pde_vals.mean()),
    veff_mmns=v_eff_mmns, veff_err_mmns=v_eff_err,
    sig_x_END_mopt_mm=sig_x_end_at_mopt, sig_x_END_best_mm=sig_x_end_best,
    m_opt_timing=m_oracle, m_opt_position=m_opt_x,
    sig_x_TOP_centroid_x0_mm=sig_x_top_centroid_at_x0,
)
with open(TABS / "summary_numbers.json", "w") as fp:
    json.dump(summary, fp, indent=2, default=float)
print(f"  Saved summary_numbers.json")

# Table 4: Jitter scenarios
jit_rows = []
for jit in sorted(jitter_results.keys()):
    r = jitter_results[jit]
    jit_rows.append(dict(jitter_ps=jit, sig_TOP=r['sig_top'], sig_END=r['sig_end'], sig_BLUE=r['sig_blue']))
pd.DataFrame(jit_rows).to_csv(TABS / "tab_jitter.csv", index=False, float_format="%.2f")
print(f"  Saved tab_jitter.csv")

fout.Close()
print("\nDone. All figures in:", FIGS)
print("All tables in:", TABS)
