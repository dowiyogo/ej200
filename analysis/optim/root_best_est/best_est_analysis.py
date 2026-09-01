#!/usr/bin/env python3
"""
best_est_analysis.py — Independent reanalysis of EndSparseTop optimal configuration.

Reads ROOT files from results/scan_end_vik_sparse_top_v2/
Writes ROOT figures to analysis/optim/root_best_est/

Produces:
  fig1_kscan.{pdf,png}       — sigma_TOP(k) for N=20 EJ-230 x=0
  fig2_residuals.{pdf,png}   — T residual histograms (END / TOP / BLUE)
  fig3_sigma_vs_x.{pdf,png}  — sigma vs position for N=20
  fig4_sigma_vs_ntop.{pdf,png} — sigma_TOP vs N_top with fit
  figures.root               — all ROOT objects
  verification.json          — all numbers for presentation cross-check
"""

import ROOT
import numpy as np
import uproot
import json
import sys
from pathlib import Path

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

REPO   = Path("/home/rrios/ej200")
INDIR  = REPO / "results/scan_end_vik_sparse_top_v2"
OUTDIR = REPO / "analysis/optim/root_best_est"
OUTDIR.mkdir(parents=True, exist_ok=True)

N_BOOT = 400   # bootstrap iterations
RNG    = np.random.default_rng(2026)

# ═══════════════════════════════════════════════════════════════════════════════
# ROOT STYLE
# ═══════════════════════════════════════════════════════════════════════════════

def set_style():
    s = ROOT.gStyle
    s.SetOptStat(0); s.SetOptFit(0)
    s.SetCanvasColor(0); s.SetPadColor(0); s.SetFrameFillColor(0)
    s.SetCanvasBorderMode(0); s.SetPadBorderMode(0); s.SetFrameBorderMode(0)
    for ax in ("x","y","z"):
        s.SetTitleFont(42, ax); s.SetLabelFont(42, ax)
        s.SetTitleSize(0.052, ax); s.SetLabelSize(0.047, ax)
        s.SetTitleOffset(1.0, ax)
    s.SetTextFont(42); s.SetTextSize(0.05)
    s.SetPadTopMargin(0.07); s.SetPadRightMargin(0.05)
    s.SetPadBottomMargin(0.13); s.SetPadLeftMargin(0.15)
    s.SetTitleOffset(1.25, "y")
    s.SetLegendBorderSize(0); s.SetLegendFont(42); s.SetLegendFillColor(0)
    s.SetHistLineWidth(2)
    s.SetPadTickX(1); s.SetPadTickY(1)
    s.SetGridColor(ROOT.kGray+1); s.SetGridStyle(3)
    s.SetPadGridX(True); s.SetPadGridY(True)

set_style()
fout = ROOT.TFile(str(OUTDIR / "figures.root"), "RECREATE")

# ═══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def rsigma(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 20: return np.nan
    q1, q3 = np.percentile(a, [25, 75])
    return (q3 - q1) / 1.349

def bootstrap_rsigma(arr, n_boot=N_BOOT):
    arr = np.asarray(arr, float); arr = arr[np.isfinite(arr)]
    n = len(arr)
    return np.std([rsigma(arr[RNG.integers(0,n,n)]) for _ in range(n_boot)])

def load_hits(path):
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id','face_type','time_ns'], library='np')
    return a['event_id'], a['face_type'], a['time_ns']

def prep_side(ev, ft, tn, ft_id):
    m = ft == ft_id
    ev_m, tn_m = ev[m], tn[m]
    idx = np.lexsort([tn_m, ev_m])
    ev_s, tn_s = ev_m[idx], tn_m[idx]
    uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
    return uniq, starts, cnts, tn_s

def scan_m_end(tL, tR, sL, sR, cL, cR, uL, uR, Mmax=120):
    """Scan m for T0_END m-average. Returns arrays of sigma vs m."""
    com = np.intersect1d(uL, uR)
    iL = np.searchsorted(uL, com); iR = np.searchsorted(uR, com)
    cLc = cL[iL]; cRc = cR[iR]; sLc = sL[iL]; sRc = sR[iR]
    Mmax = min(Mmax, int(np.percentile(cLc,2)), int(np.percentile(cRc,2)))
    sigs = np.full(Mmax, np.nan)
    for m in range(1, Mmax+1):
        mask = (cLc>=m)&(cRc>=m)
        if mask.sum()<30: continue
        T0e = (tL[sLc[mask,None]+np.arange(m)].mean(1)+
               tR[sRc[mask,None]+np.arange(m)].mean(1))/2.
        sigs[m-1] = rsigma(T0e - T0e.mean())*1e3
    return sigs, com, iL, iR, cLc, cRc, sLc, sRc

def scan_k_top(tTop, sT, cT, uT, Kmax=80):
    """Scan k for T0_TOP k-th order stat. Returns arrays of sigma vs k."""
    Kmax = min(Kmax, int(np.percentile(cT,2)))
    sigs = np.full(Kmax, np.nan)
    for k in range(1, Kmax+1):
        mask = cT>=k
        if mask.sum()<30: continue
        T0t = tTop[sT[mask]+(k-1)]
        sigs[k-1] = rsigma(T0t - T0t.mean())*1e3
    return sigs

def compute_blue(tL, tR, tTop, sLc, sRc, sTa, cLa, cRa, cTa, m_opt, k_opt):
    """
    BLUE combination at fixed m*, k*.
    Returns dict with T0e, T0t, T0c arrays and weights/correlations.
    """
    mask = (cLa>=m_opt)&(cRa>=m_opt)&(cTa>=k_opt)
    if mask.sum()<30:
        return None
    T0e = (tL[sLc[mask,None]+np.arange(m_opt)].mean(1)+
           tR[sRc[mask,None]+np.arange(m_opt)].mean(1))/2.
    T0t = tTop[sTa[mask]+(k_opt-1)]

    # Mean-subtract so residuals are centred
    T0e_c = T0e - T0e.mean()
    T0t_c = T0t - T0t.mean()

    # IQR-based sigmas
    se  = rsigma(T0e_c)
    st  = rsigma(T0t_c)
    se2 = se**2; st2 = st**2

    # Full empirical covariance
    C   = np.cov(T0e_c, T0t_c)
    cov = C[0,1]
    var_e_emp = C[0,0]; var_t_emp = C[1,1]
    rho = cov / np.sqrt(var_e_emp * var_t_emp)

    # BLUE weight (from IQR sigmas + empirical cov)
    denom = se2 + st2 - 2*cov
    if denom > 1e-30:
        w_end = np.clip((st2-cov)/denom, 0., 1.)
    else:
        w_end = 0.5
    w_top = 1.-w_end

    T0c = w_end*T0e_c + w_top*T0t_c
    sig_comb_ev = rsigma(T0c - T0c.mean())*1e3

    # Verify via covariance matrix
    sig_comb_cov = np.sqrt(max(0, w_end**2*se2 + w_top**2*st2 + 2*w_end*w_top*cov))*1e3

    return dict(
        T0e=T0e_c, T0t=T0t_c, T0c=T0c,
        w_end=w_end, w_top=w_top, rho=rho,
        sig_e=se*1e3, sig_t=st*1e3,
        sig_comb_ev=sig_comb_ev, sig_comb_cov=sig_comb_cov,
        n_ev=mask.sum(), mask=mask
    )

def gauss_fit(h):
    """Fit Gaussian to core (±2σ_IQR around median). Returns TF1 or None."""
    n = h.GetEntries()
    if n < 30: return None
    # estimate core range from quantiles
    q = np.array([0.25, 0.5, 0.75])
    xq = np.zeros(3)
    h.GetQuantiles(3, xq, q)
    iqr = (xq[2]-xq[0])/1.349
    med = xq[1]
    xlo = med - 2.5*iqr; xhi = med + 2.5*iqr
    f = ROOT.TF1("fg","gaus", xlo, xhi)
    f.SetLineColor(ROOT.kRed+1); f.SetLineWidth(2)
    res = h.Fit(f, "RSQN")   # fit in range, quiet, no draw
    if int(res) != 0: return None
    return f

def save_canvas(c, name):
    fout.cd()
    c.Write(name)
    c.SaveAs(str(OUTDIR / f"{name}.pdf"))
    c.SaveAs(str(OUTDIR / f"{name}.png"))
    print(f"  Saved {name}.{{pdf,png}}")

# ═══════════════════════════════════════════════════════════════════════════════
# LOAD MAIN FILE: EJ-230 N=20 x=0
# ═══════════════════════════════════════════════════════════════════════════════
print("="*60)
print("Loading EJ-230 N=20 x=0 …")
ev, ft, tn = load_hits(INDIR/"photon_hits_ej230_ntop20_x0.root")
print(f"  Total hits: {len(ev):,}")

uL,sL,cL,tL = prep_side(ev,ft,tn,0)
uR,sR,cR,tR = prep_side(ev,ft,tn,1)
uT,sT,cT,tTop = prep_side(ev,ft,tn,2)
print(f"  face_type=0: {len(uL)} events, {cL.mean():.1f} hits/ev")
print(f"  face_type=1: {len(uR)} events, {cR.mean():.1f} hits/ev")
print(f"  face_type=2: {len(uT)} events, {cT.mean():.1f} hits/ev")

# Scan m (END)
print("  Scanning m (END m-average) …")
sigs_m, com_end, iL_ce, iR_ce, cLc, cRc, sLc, sRc = scan_m_end(tL,tR,sL,sR,cL,cR,uL,uR)
m_opt = int(np.nanargmin(sigs_m))+1
sig_end_opt = sigs_m[m_opt-1]
print(f"  m* = {m_opt}, sigma_END = {sig_end_opt:.2f} ps")

# Scan k (TOP)
print("  Scanning k (TOP k-th order stat) …")
sigs_k = scan_k_top(tTop,sT,cT,uT)
k_opt = int(np.nanargmin(sigs_k))+1
sig_top_opt = sigs_k[k_opt-1]
print(f"  k* = {k_opt}, sigma_TOP = {sig_top_opt:.2f} ps")

# Intersection of end+top events for BLUE
com_all = np.intersect1d(com_end, uT)
ia  = np.searchsorted(com_end, com_all)
iT  = np.searchsorted(uT, com_all)
cLa = cLc[ia]; cRa = cRc[ia]; sLa = sLc[ia]; sRa = sRc[ia]
cTa = cT[iT];  sTa = sT[iT]

# BLUE
print("  Computing BLUE …")
blue = compute_blue(tL,tR,tTop,sLa,sRa,sTa,cLa,cRa,cTa,m_opt,k_opt)
print(f"  w_END={blue['w_end']:.4f}, w_TOP={blue['w_top']:.4f}")
print(f"  rho(T_END,T_TOP) = {blue['rho']:.4f}")
print(f"  sigma_BLUE (event-by-event): {blue['sig_comb_ev']:.2f} ps")
print(f"  sigma_BLUE (from cov matrix): {blue['sig_comb_cov']:.2f} ps")

# Bootstrap uncertainties
print("  Bootstrapping uncertainties …")
T0e_c = blue['T0e']; T0t_c = blue['T0t']; T0c = blue['T0c']
boot_e   = bootstrap_rsigma(T0e_c)*1e3
boot_t   = bootstrap_rsigma(T0t_c)*1e3
boot_c   = bootstrap_rsigma(T0c)*1e3
print(f"  δ(sigma_END) = {boot_e:.2f} ps, δ(sigma_TOP) = {boot_t:.2f} ps, δ(sigma_BLUE) = {boot_c:.2f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# EJ-204 N=20 x=0
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading EJ-204 N=20 x=0 …")
ev4,ft4,tn4 = load_hits(INDIR/"photon_hits_ej204_ntop20_x0.root")
uL4,sL4,cL4,tL4 = prep_side(ev4,ft4,tn4,0)
uR4,sR4,cR4,tR4 = prep_side(ev4,ft4,tn4,1)
uT4,sT4,cT4,tTop4 = prep_side(ev4,ft4,tn4,2)
print(f"  face_type=2: {len(uT4)} events, {cT4.mean():.1f} hits/ev")

sigs_m4,com4,iL4c,iR4c,cL4c,cR4c,sL4c,sR4c = scan_m_end(tL4,tR4,sL4,sR4,cL4,cR4,uL4,uR4)
sigs_k4 = scan_k_top(tTop4,sT4,cT4,uT4)
m_opt4 = int(np.nanargmin(sigs_m4))+1
k_opt4 = int(np.nanargmin(sigs_k4))+1
sig_end4 = sigs_m4[m_opt4-1]; sig_top4 = sigs_k4[k_opt4-1]

com_all4 = np.intersect1d(com4, uT4)
ia4=np.searchsorted(com4,com_all4); iT4=np.searchsorted(uT4,com_all4)
cL4a=cL4c[ia4]; cR4a=cR4c[ia4]; sL4a=sL4c[ia4]; sR4a=sR4c[ia4]
cT4a=cT4[iT4];  sT4a=sT4[iT4]

blue4 = compute_blue(tL4,tR4,tTop4,sL4a,sR4a,sT4a,cL4a,cR4a,cT4a,m_opt4,k_opt4)
boot_t4 = bootstrap_rsigma(blue4['T0t'])*1e3
boot_c4 = bootstrap_rsigma(blue4['T0c'])*1e3
print(f"  EJ-204: m*={m_opt4}, k*={k_opt4}")
print(f"  sigma_END={sig_end4:.2f} sigma_TOP={sig_top4:.2f} sigma_BLUE={blue4['sig_comb_ev']:.2f} ps")
print(f"  δ(sigma_TOP)={boot_t4:.2f} δ(sigma_BLUE)={boot_c4:.2f} ps")
print(f"  w_END={blue4['w_end']:.4f} rho={blue4['rho']:.4f}")

# EJ-204 vs EJ-230 statistical significance
diff = blue4['sig_comb_ev'] - blue['sig_comb_ev']
diff_unc = np.sqrt(boot_c4**2 + boot_c**2)
print(f"\n  EJ-204 - EJ-230 BLUE = {diff:.2f} ± {diff_unc:.2f} ps ({abs(diff)/diff_unc:.1f} sigma)")

# ═══════════════════════════════════════════════════════════════════════════════
# N_top SCAN: N=4,8,14,20 at x=0, EJ-230
# ═══════════════════════════════════════════════════════════════════════════════
print("\nN_top scan (EJ-230, x=0) …")
ntop_vals = [4, 8, 14, 20]
ntop_results = {}
for N in ntop_vals:
    path = INDIR / f"photon_hits_ej230_ntop{N}_x0.root"
    evN,ftN,tnN = load_hits(path)
    uLN,sLN,cLN,tLN = prep_side(evN,ftN,tnN,0)
    uRN,sRN,cRN,tRN = prep_side(evN,ftN,tnN,1)
    uTN,sTN,cTN,tTN = prep_side(evN,ftN,tnN,2)
    sigsN_m,comN,iLNc,iRNc,cLNc,cRNc,sLNc,sRNc = scan_m_end(tLN,tRN,sLN,sRN,cLN,cRN,uLN,uRN)
    sigsN_k = scan_k_top(tTN,sTN,cTN,uTN)
    mN = int(np.nanargmin(sigsN_m))+1
    kN = int(np.nanargmin(sigsN_k))+1

    com_allN = np.intersect1d(comN, uTN)
    iaN=np.searchsorted(comN,com_allN); iTN=np.searchsorted(uTN,com_allN)
    cLNa=cLNc[iaN]; cRNa=cRNc[iaN]; sLNa=sLNc[iaN]; sRNa=sRNc[iaN]
    cTNa=cTN[iTN];  sTNa=sTN[iTN]

    blueN = compute_blue(tLN,tRN,tTN,sLNa,sRNa,sTNa,cLNa,cRNa,cTNa,mN,kN)
    boot_tN = bootstrap_rsigma(blueN['T0t'])*1e3
    boot_cN = bootstrap_rsigma(blueN['T0c'])*1e3
    ntop_results[N] = dict(
        m_opt=mN, k_opt=kN,
        sig_end=sigsN_m[mN-1], sig_top=sigsN_k[kN-1],
        sig_comb=blueN['sig_comb_ev'], sig_comb_cov=blueN['sig_comb_cov'],
        w_end=blueN['w_end'], w_top=blueN['w_top'], rho=blueN['rho'],
        npe_top=cTN.mean(),
        boot_t=boot_tN, boot_c=boot_cN,
        sigs_k=sigsN_k   # keep full k-scan for each N
    )
    print(f"  N={N:2d}: m*={mN} k*={kN} "
          f"σ_end={sigsN_m[mN-1]:.1f} σ_top={sigsN_k[kN-1]:.1f} "
          f"σ_BLUE={blueN['sig_comb_ev']:.1f}±{boot_cN:.1f} ps  "
          f"w_TOP={blueN['w_top']:.3f} ρ={blueN['rho']:.3f}")

# ═══════════════════════════════════════════════════════════════════════════════
# POSITION SCAN: N=20 EJ-230, all x
# ═══════════════════════════════════════════════════════════════════════════════
print("\nPosition scan (EJ-230 N=20, all x) …")
xpos_vals = list(range(-600, 601, 100))
pos_results = {}
for x in xpos_vals:
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    if not path.exists():
        print(f"  MISSING: {path.name}"); continue
    evX,ftX,tnX = load_hits(path)
    uLX,sLX,cLX,tLX = prep_side(evX,ftX,tnX,0)
    uRX,sRX,cRX,tRX = prep_side(evX,ftX,tnX,1)
    uTX,sTX,cTX,tTX = prep_side(evX,ftX,tnX,2)
    sigsX_m,comX,_,_,cLXc,cRXc,sLXc,sRXc = scan_m_end(tLX,tRX,sLX,sRX,cLX,cRX,uLX,uRX)
    sigsX_k = scan_k_top(tTX,sTX,cTX,uTX)
    mX = int(np.nanargmin(sigsX_m))+1
    kX = int(np.nanargmin(sigsX_k))+1

    com_allX = np.intersect1d(comX, uTX)
    iaX=np.searchsorted(comX,com_allX); iTX2=np.searchsorted(uTX,com_allX)
    cLXa=cLXc[iaX]; cRXa=cRXc[iaX]; sLXa=sLXc[iaX]; sRXa=sRXc[iaX]
    cTXa=cTX[iTX2]; sTXa=sTX[iTX2]

    blueX = compute_blue(tLX,tRX,tTX,sLXa,sRXa,sTXa,cLXa,cRXa,cTXa,mX,kX)
    boot_cX = bootstrap_rsigma(blueX['T0c'])*1e3
    pos_results[x] = dict(
        m_opt=mX, k_opt=kX,
        sig_end=sigsX_m[mX-1], sig_top=sigsX_k[kX-1],
        sig_comb=blueX['sig_comb_ev'], boot_c=boot_cX,
        npe_top=cTX.mean()
    )
    print(f"  x={x:+4d}: k*={kX} σ_TOP={sigsX_k[kX-1]:.1f} σ_BLUE={blueX['sig_comb_ev']:.1f}±{boot_cX:.1f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — σ_TOP(k) for N=20 EJ-230 x=0
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFigure 1: k-scan …")
c1 = ROOT.TCanvas("c_kscan","k-scan",800,600)
c1.SetLeftMargin(0.15); c1.SetBottomMargin(0.14)

# Main graph: N=20 k-scan with bootstrap error bars
k_vals = np.where(np.isfinite(sigs_k))[0]+1
s_vals = sigs_k[np.isfinite(sigs_k)]
# Quick boot per k (just a few points for error bars — use Poisson approx for speed)
# Analytical approx: delta_sigma ~ sigma / sqrt(2*Nev) for IQR
n_ev_top = int((cT>=k_opt).sum())
delta_k = np.array([sigs_k[k-1]/np.sqrt(2*(cT>=k).sum()) if np.isfinite(sigs_k[k-1]) else 0
                    for k in k_vals])

g20 = ROOT.TGraphErrors(len(k_vals),
      k_vals.astype('d'), s_vals,
      np.zeros(len(k_vals)), delta_k)
g20.SetMarkerStyle(20); g20.SetMarkerSize(0.9)
g20.SetMarkerColor(ROOT.kAzure+2); g20.SetLineColor(ROOT.kAzure+2)
g20.SetLineWidth(2)
g20.GetXaxis().SetTitle("k  (k-th earliest photon on TOP face)")
g20.GetYaxis().SetTitle("#sigma(T_{0}^{TOP})  [ps]")
g20.GetYaxis().SetRangeUser(0., 80.)
g20.SetTitle("")
g20.Draw("APL")

# Mark minimum
kstar = ROOT.TGraph(1, np.array([float(k_opt)]), np.array([sig_top_opt]))
kstar.SetMarkerStyle(29); kstar.SetMarkerSize(2.2)
kstar.SetMarkerColor(ROOT.kRed+1)
kstar.Draw("P same")

# Line at k*
lk = ROOT.TLine(k_opt, 0., k_opt, sig_top_opt)
lk.SetLineColor(ROOT.kRed+1); lk.SetLineStyle(2); lk.SetLineWidth(2)
lk.Draw()

# Also overlay N=8 in grey for comparison
sigs_k8 = ntop_results[8]['sigs_k']
k8_vals = np.where(np.isfinite(sigs_k8))[0]+1
s8_vals = sigs_k8[np.isfinite(sigs_k8)]
g8 = ROOT.TGraph(len(k8_vals), k8_vals.astype('d'), s8_vals)
g8.SetMarkerStyle(24); g8.SetMarkerSize(0.8)
g8.SetMarkerColor(ROOT.kGray+2); g8.SetLineColor(ROOT.kGray+2)
g8.SetLineWidth(1); g8.SetLineStyle(3)
g8.Draw("PL same")

leg1 = ROOT.TLegend(0.55, 0.62, 0.93, 0.90)
leg1.SetTextSize(0.043)
leg1.AddEntry(g20, f"N_{{top}}=20,  k*={k_opt},  #sigma*={sig_top_opt:.1f} ps","pl")
leg1.AddEntry(kstar, "Minimum  k*","p")
leg1.AddEntry(g8,  f"N_{{top}}=8   (ref)","pl")
leg1.Draw()

lat1 = ROOT.TLatex()
lat1.SetNDC(); lat1.SetTextSize(0.040); lat1.SetTextFont(42)
lat1.DrawLatex(0.17, 0.88, "EJ-230,  x = 0")
lat1.DrawLatex(0.17, 0.83, "IQR-based #sigma, 5000 events")

save_canvas(c1,"fig1_kscan")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Residual histograms: T_END, T_TOP, T_BLUE
# ═══════════════════════════════════════════════════════════════════════════════
print("Figure 2: residual histograms …")

# Convert to ps for histogram
T0e_ps = blue['T0e']*1e3
T0t_ps = blue['T0t']*1e3
T0c_ps = blue['T0c']*1e3
# Align means (already mean-subtracted)

# Range: ±4 * largest sigma
rng_hist = 4.*max(blue['sig_e'],blue['sig_t'])
nbins = 100

hE = ROOT.TH1D("hE","",nbins,-rng_hist,rng_hist)
hT = ROOT.TH1D("hT","",nbins,-rng_hist,rng_hist)
hC = ROOT.TH1D("hC","",nbins,-rng_hist,rng_hist)
for name,(h,arr) in zip(["hE","hT","hC"],
                         [("hE",T0e_ps),("hT",T0t_ps),("hC",T0c_ps)]):
    hh = eval(name)
    for v in arr: hh.Fill(v)

hE.SetLineColor(ROOT.kAzure+3);   hE.SetLineWidth(2)
hT.SetLineColor(ROOT.kOrange+1);  hT.SetLineWidth(2)
hC.SetLineColor(ROOT.kGreen+2);   hC.SetLineWidth(2); hC.SetLineStyle(1)

# Gaussian fits
fgE = gauss_fit(hE)
fgT = gauss_fit(hT)
fgC = gauss_fit(hC)

c2 = ROOT.TCanvas("c_residuals","residuals",900,620)
c2.SetLeftMargin(0.15); c2.SetBottomMargin(0.14)

ymax = max(hE.GetMaximum(), hT.GetMaximum(), hC.GetMaximum())*1.35
hE.GetXaxis().SetTitle("#DeltaT  =  T_{0}^{rec} #minus #LTT_{0}^{rec}#GT  [ps]")
hE.GetYaxis().SetTitle("Events / bin")
hE.SetMaximum(ymax); hE.SetMinimum(0.)
hE.SetTitle("")
hE.Draw("HIST")
hT.Draw("HIST same")
hC.Draw("HIST same")

if fgE: fgE.SetLineColor(ROOT.kAzure+3);   fgE.Draw("same")
if fgT: fgT.SetLineColor(ROOT.kOrange+1);  fgT.Draw("same")
if fgC: fgC.SetLineColor(ROOT.kGreen+2);   fgC.Draw("same")

def fmt_stat(h,fg,lbl,col):
    s_iqr = rsigma(np.array([h.GetBinCenter(i+1) for i in range(nbins)
                               for _ in range(int(h.GetBinContent(i+1)))]))
    s_str = f"#sigma_{{IQR}}={s_iqr:.1f} ps"
    if fg and fg.GetNDF()>0:
        s_str += f",  #sigma_{{G}}={fg.GetParameter(2):.1f} ps"
    return f"{lbl}: {s_str}"

leg2 = ROOT.TLegend(0.38, 0.60, 0.93, 0.90)
leg2.SetTextSize(0.037)
leg2.AddEntry(hE, fmt_stat(hE,fgE,"END  (m*=%d)"%m_opt,ROOT.kAzure+3),"l")
leg2.AddEntry(hT, fmt_stat(hT,fgT,"TOP  (k*=%d)"%k_opt,ROOT.kOrange+1),"l")
leg2.AddEntry(hC, fmt_stat(hC,fgC,"BLUE (w_{TOP}=%.3f)"%blue['w_top'],ROOT.kGreen+2),"l")
leg2.Draw()

lat2 = ROOT.TLatex(); lat2.SetNDC(); lat2.SetTextSize(0.038); lat2.SetTextFont(42)
lat2.DrawLatex(0.17,0.88,"EJ-230,  N_{top}=20,  x = 0,  5000 events")
lat2.DrawLatex(0.17,0.83,"Solid: histogram   Dashed: Gaussian core fit")

save_canvas(c2,"fig2_residuals")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — σ vs x for N=20 EJ-230
# ═══════════════════════════════════════════════════════════════════════════════
print("Figure 3: sigma vs x …")
xs  = sorted(pos_results.keys())
sig_comb_arr = np.array([pos_results[x]['sig_comb'] for x in xs])
sig_err_arr  = np.array([pos_results[x]['boot_c']   for x in xs])
sig_end_arr  = np.array([pos_results[x]['sig_end']  for x in xs])

xs_d = np.array(xs, dtype='d')
g_blue_x = ROOT.TGraphErrors(len(xs), xs_d, sig_comb_arr, np.zeros(len(xs)), sig_err_arr)
g_blue_x.SetMarkerStyle(20); g_blue_x.SetMarkerSize(1.1)
g_blue_x.SetMarkerColor(ROOT.kGreen+2); g_blue_x.SetLineColor(ROOT.kGreen+2); g_blue_x.SetLineWidth(2)

g_end_x = ROOT.TGraph(len(xs), xs_d, sig_end_arr)
g_end_x.SetMarkerStyle(24); g_end_x.SetMarkerSize(1.0)
g_end_x.SetMarkerColor(ROOT.kAzure+3); g_end_x.SetLineColor(ROOT.kAzure+3)
g_end_x.SetLineWidth(2); g_end_x.SetLineStyle(2)

c3 = ROOT.TCanvas("c_vs_x","sigma vs x",880,620)
c3.SetLeftMargin(0.15); c3.SetBottomMargin(0.14)

mg3 = ROOT.TMultiGraph()
mg3.Add(g_end_x,"PL"); mg3.Add(g_blue_x,"PL")
mg3.Draw("A")
mg3.GetXaxis().SetTitle("Muon position x  [mm]")
mg3.GetYaxis().SetTitle("#sigma(T_{0})  [ps]")
mg3.GetYaxis().SetRangeUser(0., 80.)
mg3.SetTitle("")
mg3.GetXaxis().SetLimits(-700.,700.)

# Constant fit to BLUE
fc = ROOT.TF1("fc","[0]",-650.,650.)
fc.SetLineColor(ROOT.kRed+1); fc.SetLineWidth(2); fc.SetLineStyle(2)
res3 = g_blue_x.Fit(fc,"RSQ")
C_fit  = fc.GetParameter(0)
eC_fit = fc.GetParError(0)
chi2   = fc.GetChisquare()
ndf    = fc.GetNDF()
pval   = ROOT.TMath.Prob(chi2, ndf) if ndf>0 else 1.

leg3 = ROOT.TLegend(0.18, 0.65, 0.60, 0.90)
leg3.SetTextSize(0.040)
leg3.AddEntry(g_blue_x,"BLUE: END+TOP  (N_{top}=20)","pl")
leg3.AddEntry(fc,f"Const fit: C = {C_fit:.1f}#pm{eC_fit:.1f} ps","l")
leg3.AddEntry(g_end_x,"END only (m*)","pl")
leg3.Draw()

lat3 = ROOT.TLatex(); lat3.SetNDC(); lat3.SetTextSize(0.038); lat3.SetTextFont(42)
lat3.DrawLatex(0.18,0.57,f"#chi^{{2}}/ndf = {chi2:.1f}/{ndf}")
lat3.DrawLatex(0.18,0.52,f"p-value = {pval:.3f}")
lat3.DrawLatex(0.18,0.47,f"Range: {sig_comb_arr.min():.1f}#{{}}{sig_comb_arr.max():.1f} ps")
lat3.DrawLatex(0.55,0.88,"EJ-230,  N_{top}=20")

save_canvas(c3,"fig3_sigma_vs_x")
print(f"  Const fit: C = {C_fit:.2f} ± {eC_fit:.2f} ps, chi2/ndf={chi2:.1f}/{ndf}, p={pval:.3f}")
print(f"  Range: {sig_comb_arr.min():.2f}–{sig_comb_arr.max():.2f} ps")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — σ_TOP vs N_top with √N fit
# ═══════════════════════════════════════════════════════════════════════════════
print("Figure 4: sigma_TOP vs N_top …")
ntop_arr = np.array(ntop_vals, dtype='d')
stop_arr = np.array([ntop_results[N]['sig_top'] for N in ntop_vals])
etop_arr = np.array([ntop_results[N]['boot_t']  for N in ntop_vals])
npe_arr  = np.array([ntop_results[N]['npe_top'] for N in ntop_vals])

# Fit: sigma = A / sqrt(Npe)
gNpe = ROOT.TGraphErrors(4, npe_arr, stop_arr, np.zeros(4), etop_arr)
gNpe.SetMarkerStyle(20); gNpe.SetMarkerSize(1.3)
gNpe.SetMarkerColor(ROOT.kOrange+1); gNpe.SetLineColor(ROOT.kOrange+1); gNpe.SetLineWidth(2)

npe_lo = npe_arr.min()*0.6; npe_hi = npe_arr.max()*1.1
ffit = ROOT.TF1("ffit","[0]/sqrt(x)", npe_lo, npe_hi)
ffit.SetParameter(0, 15.*np.sqrt(500.))
ffit.SetLineColor(ROOT.kRed+1); ffit.SetLineWidth(2); ffit.SetLineStyle(2)
resFit = gNpe.Fit(ffit,"RSQ")
A_fit  = ffit.GetParameter(0); eA_fit = ffit.GetParError(0)
chi2f  = ffit.GetChisquare(); ndff = ffit.GetNDF()
pvalf  = ROOT.TMath.Prob(chi2f,ndff) if ndff>0 else 1.

c4 = ROOT.TCanvas("c_vs_ntop","sigma vs ntop",880,620)
c4.SetLeftMargin(0.15); c4.SetBottomMargin(0.14)
c4.SetLogx(); c4.SetLogy()

# Also: sigma vs N_top (integer)
gNtop = ROOT.TGraphErrors(4, ntop_arr, stop_arr, np.zeros(4), etop_arr)
gNtop.SetMarkerStyle(21); gNtop.SetMarkerSize(1.2)
gNtop.SetMarkerColor(ROOT.kAzure+2); gNtop.SetLineColor(ROOT.kAzure+2); gNtop.SetLineWidth(2)

gNpe.GetXaxis().SetTitle("N_{pe}^{TOP} (mean photons detected)")
gNpe.GetYaxis().SetTitle("#sigma(T_{0}^{TOP})  [ps]")
gNpe.SetTitle("")
gNpe.GetXaxis().SetLimits(npe_lo, npe_hi)
gNpe.GetYaxis().SetRangeUser(10., 120.)
gNpe.Draw("APL")
ffit.Draw("same")

# Annotate data points with N label
for i,N in enumerate(ntop_vals):
    t4 = ROOT.TLatex()
    t4.SetTextSize(0.038); t4.SetTextFont(42); t4.SetTextColor(ROOT.kOrange+1)
    t4.DrawLatex(npe_arr[i]*1.04, stop_arr[i], f"N={N}")

leg4 = ROOT.TLegend(0.50, 0.65, 0.93, 0.88)
leg4.SetTextSize(0.040)
leg4.AddEntry(gNpe,"#sigma_{TOP}^{IQR}","pl")
leg4.AddEntry(ffit,f"A/#sqrt{{N_{{pe}}}}, A={A_fit:.0f} ps^{{1/2}}","l")
leg4.Draw()

lat4 = ROOT.TLatex(); lat4.SetNDC(); lat4.SetTextSize(0.038); lat4.SetTextFont(42)
lat4.DrawLatex(0.20,0.25,f"#chi^{{2}}/ndf = {chi2f:.2f}/{ndff},  p = {pvalf:.3f}")
lat4.DrawLatex(0.55,0.88,"EJ-230,  x = 0")

save_canvas(c4,"fig4_sigma_vs_ntop")
print(f"  Fit A/sqrt(Npe): A = {A_fit:.1f}±{eA_fit:.1f}, chi2/ndf={chi2f:.2f}/{ndff}, p={pvalf:.3f}")

# ═══════════════════════════════════════════════════════════════════════════════
# SAVE ROOT OBJECTS
# ═══════════════════════════════════════════════════════════════════════════════
fout.cd()
for obj in [g20, kstar, g8, hE, hT, hC, g_blue_x, g_end_x, gNpe, ffit]:
    if obj: obj.Write()
if fgE: fgE.Write("fgE")
if fgT: fgT.Write("fgT")
if fgC: fgC.Write("fgC")
fout.Close()
print(f"\nSaved figures.root")

# ═══════════════════════════════════════════════════════════════════════════════
# VERIFICATION TABLE (JSON)
# ═══════════════════════════════════════════════════════════════════════════════
# EJ-230 END+Vik reference (from existing CSV)
import pandas as pd
csv = pd.read_csv(REPO/"analysis/optim/phase_sparse_top_optimal.csv")

def get_csv(sc,N,x):
    row = csv[(csv.scint==sc)&(csv.ntop==N)&(csv.x==x)]
    if row.empty: return {}
    r = row.iloc[0]
    return dict(sig_end=r.sig_end_opt, m_opt=int(r.m_opt),
                sig_top=r.sig_top_opt, k_opt=int(r.k_opt),
                sig_comb=r.sig_comb)

ver = {
    "ej230_n20_x0": {
        "from_root_files": {
            "sig_end": sig_end_opt, "m_opt": m_opt,
            "sig_top": sig_top_opt, "k_opt": k_opt,
            "sig_comb": blue['sig_comb_ev'],
            "sig_comb_fromcov": blue['sig_comb_cov'],
            "w_end": blue['w_end'], "w_top": blue['w_top'],
            "rho": blue['rho'], "npe_top": float(cT.mean()),
            "boot_e": boot_e, "boot_t": boot_t, "boot_c": boot_c
        },
        "from_csv": get_csv("ej230",20,0)
    },
    "ej204_n20_x0": {
        "from_root_files": {
            "sig_end": sig_end4, "m_opt": m_opt4,
            "sig_top": sig_top4, "k_opt": k_opt4,
            "sig_comb": blue4['sig_comb_ev'],
            "w_end": blue4['w_end'], "w_top": blue4['w_top'],
            "rho": blue4['rho'], "npe_top": float(cT4.mean()),
            "boot_t": boot_t4, "boot_c": boot_c4
        },
        "from_csv": get_csv("ej204",20,0)
    },
    "ej230_ej204_comparison": {
        "sig_comb_ej230": blue['sig_comb_ev'], "err_ej230": boot_c,
        "sig_comb_ej204": blue4['sig_comb_ev'], "err_ej204": boot_c4,
        "diff_ps": float(blue4['sig_comb_ev'] - blue['sig_comb_ev']),
        "diff_sigma": float(abs(blue4['sig_comb_ev'] - blue['sig_comb_ev']) /
                           np.sqrt(boot_c**2 + boot_c4**2)),
        "statistically_equivalent": bool(
            abs(blue4['sig_comb_ev'] - blue['sig_comb_ev']) <
            2*np.sqrt(boot_c**2 + boot_c4**2))
    },
    "ntop_scan_ej230_x0": {
        str(N): {k: (float(v) if not isinstance(v,np.ndarray) else list(v[:20]))
                 for k,v in ntop_results[N].items() if k!='sigs_k'}
        for N in ntop_vals
    },
    "pos_scan_ej230_n20": {
        str(x): {k: float(v) for k,v in pos_results[x].items()}
        for x in pos_results
    },
    "position_scan_summary": {
        "min_sig_comb": float(sig_comb_arr.min()),
        "max_sig_comb": float(sig_comb_arr.max()),
        "peak_to_peak_ps": float(sig_comb_arr.max()-sig_comb_arr.min()),
        "const_fit_C": float(C_fit), "const_fit_err": float(eC_fit),
        "chi2_ndf": f"{chi2:.1f}/{ndf}", "pvalue": float(pval)
    },
    "fit_sigma_vs_npe": {
        "A": float(A_fit), "A_err": float(eA_fit),
        "chi2_ndf": f"{chi2f:.2f}/{ndff}", "pvalue": float(pvalf),
        "formula": "sigma = A/sqrt(Npe)"
    }
}

with open(OUTDIR/"verification.json","w") as fp:
    json.dump(ver, fp, indent=2)
print(f"Saved verification.json")

# ═══════════════════════════════════════════════════════════════════════════════
# PRINT SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("VERIFICATION SUMMARY")
print("="*60)
print("\n── EJ-230 N=20 x=0 (from ROOT files, independent reanalysis) ──")
print(f"  sigma_END  = {sig_end_opt:.2f} ± {boot_e:.2f} ps   (m*={m_opt})")
print(f"  sigma_TOP  = {sig_top_opt:.2f} ± {boot_t:.2f} ps   (k*={k_opt})")
print(f"  sigma_BLUE = {blue['sig_comb_ev']:.2f} ± {boot_c:.2f} ps")
print(f"  sigma_BLUE(from cov) = {blue['sig_comb_cov']:.2f} ps  (cross-check)")
print(f"  w_END  = {blue['w_end']:.4f}")
print(f"  w_TOP  = {blue['w_top']:.4f}")
print(f"  rho(T_END, T_TOP) = {blue['rho']:.4f}")
print(f"\n── EJ-204 N=20 x=0 ──")
print(f"  sigma_BLUE = {blue4['sig_comb_ev']:.2f} ± {boot_c4:.2f} ps")
print(f"  w_END = {blue4['w_end']:.4f}  w_TOP = {blue4['w_top']:.4f}")
d = blue4['sig_comb_ev'] - blue['sig_comb_ev']
unc = np.sqrt(boot_c**2+boot_c4**2)
eq = abs(d)/unc < 2.
print(f"  EJ-204 vs EJ-230: diff = {d:+.2f} ps ({abs(d)/unc:.1f}σ)")
print(f"  Statistically equivalent: {'YES (<2σ)' if eq else 'NO'}")
print(f"\n── BLUE weights vs N_top (EJ-230, x=0) ──")
print(f"  {'N':>4} {'sig_end':>8} {'sig_top':>8} {'rho':>7} {'w_END':>7} {'w_TOP':>7} {'sig_BLUE':>9}")
for N in ntop_vals:
    r = ntop_results[N]
    print(f"  {N:>4} {r['sig_end']:>8.1f} {r['sig_top']:>8.1f} "
          f"{r['rho']:>7.4f} {r['w_end']:>7.4f} {r['w_top']:>7.4f} "
          f"{r['sig_comb']:>8.1f}±{r['boot_c']:.1f}")
print(f"\n── Position scan (EJ-230 N=20) ──")
print(f"  Range: {sig_comb_arr.min():.2f}–{sig_comb_arr.max():.2f} ps")
print(f"  Const fit: C = {C_fit:.2f}±{eC_fit:.2f} ps, chi2/ndf={chi2:.1f}/{ndf}, p={pval:.3f}")
print(f"\n── σ vs N_pe fit: A/sqrt(Npe) ──")
print(f"  A = {A_fit:.1f}±{eA_fit:.1f}, chi2/ndf={chi2f:.2f}/{ndff}, p={pvalf:.3f}")
print(f"  Sanity: tau_rise/sqrt(1341)=500ps/sqrt(1341) = {500./np.sqrt(1341.):.1f} ps")
print("="*60)
print("Done.")
