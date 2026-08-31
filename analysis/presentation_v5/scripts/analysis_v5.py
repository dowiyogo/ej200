#!/usr/bin/env python3
"""
analysis_v5.py — Scientific corrections and new figures for talk_v5.

New vs v4:
  1. Global-fixed sigma(x): sigma_BLUE vs x using global k=3, m=8 across all x
  2. Oracle sigma(x): locally optimal k(x), m(x) — diagnostic only
  3. Δt nonlinearity: residual r(x) = <Δt(x)> - (a + b*x)
  4. TOP position leave-one-out: train on N-1 x positions, predict the held-out one
  5. PDE clarification: correct definition of the 0.607 number
  6. Corrected material table: EJ-230 att = 1200 mm (from code)

Verified material properties from src/Materials.cc:
  EJ-200: att=3800mm, tau=2.1ns, Y=10000/MeV, n=1.58
  EJ-204: att=1600mm, tau=1.8ns, Y=10400/MeV, n=1.58
  EJ-230: att=1200mm, tau=1.5ns, Y=9700/MeV,  n=1.58
  Reflector: dielectric_dielectric, R=0.95 constant (CreateBarSkinReflector)
"""

import ROOT
import numpy as np
import uproot
import json
from pathlib import Path

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

REPO   = Path("/home/rrios/ej200")
INDIR  = REPO / "results/scan_end_vik_sparse_top_v2"
OUTDIR = REPO / "analysis/presentation_v5"
FIGS   = OUTDIR / "figs"
TABS   = OUTDIR / "tables"

RNG    = np.random.default_rng(20260820)
N_BOOT = 500   # reduced for speed; 1000 for final

# Global estimator parameters (from v4 train/test at x=0)
K_GLOBAL = 3
M_GLOBAL = 8

TOP_XPOS = np.array([
    -664.3, -594.3, -524.3, -454.3, -384.4, -314.4, -244.4, -174.4,
    -104.5,  -34.5,  34.5,  104.5,  174.4,  244.4,  314.4,  384.4,
     454.3,  524.3,  594.3,  664.3
])

# Colour scheme — must be consistent across all plots
COL_END  = ROOT.kAzure+2       # blue
COL_TOP  = ROOT.kOrange+1      # orange
COL_BLUE = ROOT.kGreen+3       # green
COL_ORA  = ROOT.kGray+2        # oracle (grey/dashed)
COL_REF  = ROOT.kRed+1         # reference/fit


# ═══════════════════════════════════════════════════════════════════════════════
# STYLE
# ═══════════════════════════════════════════════════════════════════════════════
def set_style():
    s = ROOT.gStyle
    s.SetOptStat(0); s.SetOptFit(0)
    s.SetCanvasColor(0); s.SetPadColor(0); s.SetFrameFillColor(0)
    s.SetCanvasBorderMode(0); s.SetPadBorderMode(0); s.SetFrameBorderMode(0)
    for ax in ("x", "y", "z"):
        s.SetTitleFont(42, ax); s.SetLabelFont(42, ax)
        s.SetTitleSize(0.060, ax); s.SetLabelSize(0.052, ax)
        s.SetTitleOffset(1.0, ax)
    s.SetTextFont(42); s.SetTextSize(0.05)
    s.SetPadTopMargin(0.08); s.SetPadRightMargin(0.06)
    s.SetPadBottomMargin(0.15); s.SetPadLeftMargin(0.15)
    s.SetTitleOffset(1.1, "y")
    s.SetLegendBorderSize(1); s.SetLegendFont(42); s.SetLegendFillColor(0)
    s.SetHistLineWidth(2)
    s.SetPadTickX(1); s.SetPadTickY(1)
    s.SetMarkerSize(1.3)

set_style()

fout = ROOT.TFile(str(FIGS / "analysis_v5.root"), "RECREATE")


# ═══════════════════════════════════════════════════════════════════════════════
# CORE UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════
def rsigma(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 20:
        return np.nan
    q1, q3 = np.percentile(a, [25, 75])
    return (q3 - q1) / 1.349


def bootstrap(arr, func=rsigma, n=N_BOOT):
    arr = np.asarray(arr, float)
    arr = arr[np.isfinite(arr)]
    N = len(arr)
    return np.std([func(arr[RNG.integers(0, N, N)]) for _ in range(n)])


def make_graph(x, y, ex=None, ey=None):
    x = np.asarray(x, dtype='d')
    y = np.asarray(y, dtype='d')
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


def prep_side(ev, ft, tn, ft_id):
    m = ft == ft_id
    ev_m, tn_m = ev[m], tn[m]
    idx = np.lexsort([tn_m, ev_m])
    ev_s, tn_s = ev_m[idx], tn_m[idx]
    uniq, starts, cnts = np.unique(ev_s, return_index=True, return_counts=True)
    return uniq, starts, cnts, tn_s


def get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, m):
    """m-average T0_END estimator. Returns (events, T0_ps, DeltaT_ps)."""
    com = np.intersect1d(uL, uR)
    iL = np.searchsorted(uL, com)
    iR = np.searchsorted(uR, com)
    cLc = cL[iL]; cRc = cR[iR]; sLc = sL[iL]; sRc = sR[iR]
    mask = (cLc >= m) & (cRc >= m)
    idx_ok = np.where(mask)[0]
    tL_m = np.array([tL[sLc[i]:sLc[i]+m].mean() for i in idx_ok])
    tR_m = np.array([tR[sRc[i]:sRc[i]+m].mean() for i in idx_ok])
    T0e = (tL_m + tR_m) / 2. * 1e3   # ps
    DeltaT = (tR_m - tL_m) * 1e3      # ps
    return com[mask], T0e, DeltaT


def get_T0_top(tTop, sT, cT, uT, k):
    """k-th order statistic T0_TOP. Returns (events, T0_ps)."""
    mask = cT >= k
    T0t = tTop[sT[mask] + (k - 1)] * 1e3   # ps
    return uT[mask], T0t


def intersect_estimators(ev_e, T0e, ev_t, T0t):
    com = np.intersect1d(ev_e, ev_t)
    ie = np.searchsorted(ev_e, com)
    it = np.searchsorted(ev_t, com)
    mask = np.isfinite(T0e[ie]) & np.isfinite(T0t[it])
    return T0e[ie][mask], T0t[it][mask]


def compute_blue(T0e_arr, T0t_arr):
    """BLUE combination. Covariance matrix uses ordinary Var/Cov; σIQR reported separately."""
    T0e_c = T0e_arr - T0e_arr.mean()
    T0t_c = T0t_arr - T0t_arr.mean()
    # Ordinary covariance matrix (as BLUE requires)
    C = np.cov(T0e_c, T0t_c)
    cov12 = C[0, 1]
    var_e = C[0, 0]; var_t = C[1, 1]
    # BLUE weights
    denom = var_e + var_t - 2 * cov12
    w_end = np.clip((var_t - cov12) / (denom + 1e-30), 0., 1.)
    w_top = 1. - w_end
    T0c = w_end * T0e_c + w_top * T0t_c
    rho = cov12 / (np.sqrt(var_e * var_t) + 1e-30)
    return dict(
        T0c=T0c, T0e=T0e_c, T0t=T0t_c,
        w_end=w_end, w_top=w_top, rho=rho,
        var_e=var_e, var_t=var_t, cov12=cov12,
        sig_e=rsigma(T0e_c),   # σIQR of END residuals
        sig_t=rsigma(T0t_c),   # σIQR of TOP residuals
        sig_c=rsigma(T0c)      # σIQR of BLUE residuals
    )


def load_x(x, ntop=20):
    """Load and prep all sides for given x position and N_TOP configuration."""
    path = INDIR / f"photon_hits_ej230_ntop{ntop}_x{x}.root"
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id', 'face_type', 'time_ns'], library='np')
    ev, ft, tn = a['event_id'], a['face_type'], a['time_ns']
    uL, sL, cL, tL = prep_side(ev, ft, tn, 0)
    uR, sR, cR, tR = prep_side(ev, ft, tn, 1)
    uT, sT, cT, tTop = prep_side(ev, ft, tn, 2)
    return uL, sL, cL, tL, uR, sR, cR, tR, uT, sT, cT, tTop


def save_canvas(c, name):
    fout.cd()
    c.Write(name)
    c.SaveAs(str(FIGS / f"{name}.pdf"))
    c.SaveAs(str(FIGS / f"{name}.png"))
    print(f"  Saved {name}.{{pdf,png}}")


x_scan = list(range(-600, 601, 100))   # 13 positions
NTOP_VALS = [4, 8, 14, 20]
print(f"x scan: {x_scan}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: v_eff fit + residual nonlinearity
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("1. v_eff fit and nonlinearity residual")

veff_data = []
for x in x_scan:
    uL, sL, cL, tL, uR, sR, cR, tR, uT, sT, cT, tTop = load_x(x)
    _, _, DT = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, M_GLOBAL)
    n_ev = len(DT)
    mean_DT = np.median(DT)
    std_DT  = rsigma(DT - np.median(DT))
    err_DT  = std_DT / np.sqrt(n_ev) if n_ev > 0 else 0.1
    veff_data.append((x, mean_DT, max(err_DT, 0.1)))
    print(f"  x={x:+4d}: <Δt>={mean_DT:+8.3f} ps  σ(Δt)={std_DT:.2f} ps  n={n_ev}")

veff_data = np.array(veff_data)
xs_v  = veff_data[:, 0]
dts_v = veff_data[:, 1]
err_v = veff_data[:, 2]

g_dt = make_graph(xs_v, dts_v, ey=err_v)
g_dt.SetMarkerStyle(20); g_dt.SetMarkerSize(1.3)
g_dt.SetMarkerColor(COL_END); g_dt.SetLineColor(COL_END); g_dt.SetLineWidth(2)

flin = ROOT.TF1("flin", "[0]+[1]*x", -650., 650.)
flin.SetLineColor(COL_REF); flin.SetLineWidth(2)
g_dt.Fit(flin, "RSQ")
a_fit = flin.GetParameter(0); b_fit = flin.GetParameter(1)
ea_fit = flin.GetParError(0); eb_fit = flin.GetParError(1)
chi2_lin = flin.GetChisquare(); ndf_lin = flin.GetNDF()
p_lin = ROOT.TMath.Prob(chi2_lin, ndf_lin) if ndf_lin > 0 else 1.

v_eff_mmpns = 2. / abs(b_fit)
v_eff_mmns  = v_eff_mmpns * 1e3
v_eff_err   = v_eff_mmns * (eb_fit / abs(b_fit))

print(f"\n  Linear fit: a={a_fit:.2f}±{ea_fit:.2f} ps, b={b_fit:.5f}±{eb_fit:.5f} ps/mm")
print(f"  v_eff = {v_eff_mmns:.1f} ± {v_eff_err:.1f} mm/ns  (stat.)")
print(f"  χ²/ndf = {chi2_lin:.1f}/{ndf_lin}  (large → significant nonlinearity)")

# Residuals r(x) = Δt(x) - (a + bx)
res_v = dts_v - (a_fit + b_fit * xs_v)
print(f"  Residuals r(x) [ps]:")
for xi, ri in zip(xs_v, res_v):
    print(f"    x={xi:+4.0f}: r={ri:+.2f} ps")
rms_res = np.sqrt(np.mean(res_v**2))
print(f"  Residual RMS = {rms_res:.2f} ps")

# --- FIG v5_veff: Δt vs x with linear fit ---
c_veff = ROOT.TCanvas("c_veff_v5", "veff_v5", 900, 640)
c_veff.SetLeftMargin(0.15); c_veff.SetBottomMargin(0.15)
g_dt.SetTitle("")
g_dt.GetXaxis().SetTitle("Muon position  x  [mm]")
g_dt.GetYaxis().SetTitle("#Deltat_{END}  [ps]  (t_{R} #minus t_{L})")
g_dt.GetXaxis().SetLimits(-720., 720.)
g_dt.Draw("AP")
flin.Draw("same")

pv = ROOT.TPaveText(0.17, 0.65, 0.60, 0.90, "NDC")
pv.SetFillColor(0); pv.SetBorderSize(1); pv.SetTextFont(42); pv.SetTextSize(0.042)
pv.AddText(f"#Deltat = {a_fit:.1f} + ({b_fit:.4f}) x  [ps/mm]")
pv.AddText(f"v_{{eff}} = 2/|b| = {v_eff_mmns:.0f} #pm {v_eff_err:.0f} mm/ns")
pv.AddText(f"#chi^{{2}}/ndf = {chi2_lin:.0f}/{ndf_lin} (p = {p_lin:.2e})")
pv.AddText(f"Residual RMS = {rms_res:.1f} ps")
pv.Draw()
lat_v = ROOT.TLatex(); lat_v.SetNDC(); lat_v.SetTextFont(42); lat_v.SetTextSize(0.042)
lat_v.DrawLatex(0.60, 0.88, "EJ-230,  N_{top}=20,  m=" + str(M_GLOBAL))
save_canvas(c_veff, "v5_veff_fit")

# --- FIG v5_veff_residual: nonlinearity residual ---
g_res = make_graph(xs_v, res_v, ey=err_v)
g_res.SetMarkerStyle(20); g_res.SetMarkerSize(1.3)
g_res.SetMarkerColor(COL_END); g_res.SetLineColor(COL_END); g_res.SetLineWidth(2)
g_res.SetTitle("")

c_resv = ROOT.TCanvas("c_resv", "resv", 900, 500)
c_resv.SetLeftMargin(0.15); c_resv.SetBottomMargin(0.16)
g_res.GetXaxis().SetTitle("Muon position  x  [mm]")
g_res.GetYaxis().SetTitle("Residual  r(x)  [ps]")
g_res.GetXaxis().SetLimits(-720., 720.)
g_res.Draw("AP")
zero = ROOT.TLine(-700., 0., 700., 0.)
zero.SetLineColor(ROOT.kGray+1); zero.SetLineStyle(2); zero.SetLineWidth(1)
zero.Draw()
lat_r = ROOT.TLatex(); lat_r.SetNDC(); lat_r.SetTextFont(42); lat_r.SetTextSize(0.042)
lat_r.DrawLatex(0.55, 0.88, f"RMS = {rms_res:.1f} ps")
lat_r.DrawLatex(0.17, 0.88, "Nonlinearity: r(x) = #Deltat(x) #minus (a+bx)")
save_canvas(c_resv, "v5_veff_residual")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: GLOBAL-FIXED sigma(x) — main deployable result
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print(f"2. Global-fixed sigma(x): k={K_GLOBAL}, m={M_GLOBAL}")

global_results = {}

for x in x_scan:
    uL, sL, cL, tL, uR, sR, cR, tR, uT, sT, cT, tTop = load_x(x)
    ev_e, T0e, _ = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, M_GLOBAL)
    ev_t, T0t    = get_T0_top(tTop, sT, cT, uT, K_GLOBAL)
    T0e_a, T0t_a = intersect_estimators(ev_e, T0e, ev_t, T0t)
    bl = compute_blue(T0e_a, T0t_a)

    sig_e = rsigma(T0e_a - T0e_a.mean())
    sig_t = rsigma(T0t_a - T0t_a.mean())
    sig_c = bl['sig_c']
    boot_c = bootstrap(bl['T0c'])

    global_results[x] = dict(
        sig_e=sig_e, sig_t=sig_t, sig_c=sig_c, boot_c=boot_c,
        w_end=bl['w_end'], rho=bl['rho'],
        var_e=bl['var_e'], var_t=bl['var_t'], cov12=bl['cov12']
    )
    print(f"  x={x:+4d}: σ_END={sig_e:.1f} σ_TOP={sig_t:.1f} σ_BLUE={sig_c:.1f}±{boot_c:.1f} ps"
          f"  w_END={bl['w_end']:.3f} ρ={bl['rho']:.3f}")

xs_g   = np.array(sorted(global_results.keys()), dtype='d')
sig_c_g  = np.array([global_results[x]['sig_c']  for x in sorted(global_results)], dtype='d')
boot_c_g = np.array([global_results[x]['boot_c'] for x in sorted(global_results)], dtype='d')
sig_e_g  = np.array([global_results[x]['sig_e']  for x in sorted(global_results)], dtype='d')
sig_t_g  = np.array([global_results[x]['sig_t']  for x in sorted(global_results)], dtype='d')

print(f"\n  GLOBAL σ_BLUE: mean={sig_c_g.mean():.1f} min={sig_c_g.min():.1f} max={sig_c_g.max():.1f} ps")
print(f"  Peak-to-peak variation: {sig_c_g.max()-sig_c_g.min():.1f} ps")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: ORACLE sigma(x) — local optimal k,m per position
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("3. Oracle sigma(x): locally optimal k(x), m(x)")

oracle_results = {}

for x in x_scan:
    uL, sL, cL, tL, uR, sR, cR, tR, uT, sT, cT, tTop = load_x(x)

    # Scan m for END
    best_sig_e = 1e9; best_m = 1
    for m in range(1, 25):
        ev_e, T0e, _ = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, m)
        if len(T0e) < 30: break
        s = rsigma(T0e - T0e.mean())
        if np.isfinite(s) and s < best_sig_e:
            best_sig_e = s; best_m = m

    # Scan k for TOP
    best_sig_t = 1e9; best_k = 1
    for k in range(1, 25):
        ev_t, T0t = get_T0_top(tTop, sT, cT, uT, k)
        if len(T0t) < 30: break
        s = rsigma(T0t - T0t.mean())
        if np.isfinite(s) and s < best_sig_t:
            best_sig_t = s; best_k = k

    ev_e_opt, T0e_opt, _ = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, best_m)
    ev_t_opt, T0t_opt    = get_T0_top(tTop, sT, cT, uT, best_k)
    T0e_a, T0t_a = intersect_estimators(ev_e_opt, T0e_opt, ev_t_opt, T0t_opt)
    bl = compute_blue(T0e_a, T0t_a)

    oracle_results[x] = dict(
        k=best_k, m=best_m,
        sig_e=best_sig_e, sig_t=best_sig_t,
        sig_c=bl['sig_c'], w_end=bl['w_end']
    )
    print(f"  x={x:+4d}: k*={best_k} m*={best_m} σ_END={best_sig_e:.1f}"
          f" σ_TOP={best_sig_t:.1f} σ_BLUE={bl['sig_c']:.1f} ps")

xs_o   = np.array(sorted(oracle_results.keys()), dtype='d')
sig_c_o = np.array([oracle_results[x]['sig_c'] for x in sorted(oracle_results)], dtype='d')
print(f"\n  ORACLE σ_BLUE: mean={sig_c_o.mean():.1f} min={sig_c_o.min():.1f} max={sig_c_o.max():.1f} ps")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG v5_sigma_vs_x: global (main) + oracle (reference)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n  Plotting sigma vs x …")

# TOP sensor x positions as thin vertical lines
c_svx = ROOT.TCanvas("c_svx", "svx", 900, 640)
c_svx.SetLeftMargin(0.15); c_svx.SetBottomMargin(0.15); c_svx.SetRightMargin(0.06)

# Build multigraph
g_blue_g = make_graph(xs_g, sig_c_g, ey=boot_c_g)
g_blue_g.SetMarkerStyle(20); g_blue_g.SetMarkerSize(1.3)
g_blue_g.SetMarkerColor(COL_BLUE); g_blue_g.SetLineColor(COL_BLUE); g_blue_g.SetLineWidth(3)

g_end_g = make_graph(xs_g, sig_e_g)
g_end_g.SetMarkerStyle(21); g_end_g.SetMarkerSize(1.1)
g_end_g.SetMarkerColor(COL_END); g_end_g.SetLineColor(COL_END)
g_end_g.SetLineWidth(2); g_end_g.SetLineStyle(2)

g_top_g = make_graph(xs_g, sig_t_g)
g_top_g.SetMarkerStyle(22); g_top_g.SetMarkerSize(1.1)
g_top_g.SetMarkerColor(COL_TOP); g_top_g.SetLineColor(COL_TOP)
g_top_g.SetLineWidth(2); g_top_g.SetLineStyle(3)

g_ora_g = make_graph(xs_o, sig_c_o)
g_ora_g.SetMarkerStyle(24); g_ora_g.SetMarkerSize(1.1)
g_ora_g.SetMarkerColor(COL_ORA); g_ora_g.SetLineColor(COL_ORA)
g_ora_g.SetLineWidth(2); g_ora_g.SetLineStyle(4)

mg = ROOT.TMultiGraph()
mg.Add(g_end_g, "PL")
mg.Add(g_top_g, "PL")
mg.Add(g_blue_g, "PL")
mg.Add(g_ora_g, "PL")
mg.Draw("A")
mg.SetTitle("")
mg.GetXaxis().SetTitle("Muon position  x  [mm]")
mg.GetYaxis().SetTitle("#sigma_{IQR}  [ps]")
mg.GetXaxis().SetLimits(-720., 720.)
mg.GetYaxis().SetRangeUser(0., 80.)

# Mark TOP SiPM positions
for xp in TOP_XPOS:
    tick = ROOT.TLine(xp, 0., xp, 3.)
    tick.SetLineColor(COL_TOP); tick.SetLineWidth(1)
    tick.Draw()

leg_svx = ROOT.TLegend(0.55, 0.56, 0.93, 0.90)
leg_svx.SetTextSize(0.040)
leg_svx.AddEntry(g_blue_g, f"BLUE  (global k={K_GLOBAL}, m={M_GLOBAL})", "pl")
leg_svx.AddEntry(g_top_g,  f"TOP   (global k={K_GLOBAL})", "pl")
leg_svx.AddEntry(g_end_g,  f"END   (global m={M_GLOBAL})", "pl")
leg_svx.AddEntry(g_ora_g,  "BLUE oracle k(x),m(x)", "pl")
leg_svx.Draw()

lat_s = ROOT.TLatex(); lat_s.SetNDC(); lat_s.SetTextFont(42); lat_s.SetTextSize(0.040)
lat_s.DrawLatex(0.17, 0.88, "EJ-230,  N_{top}=20")
lat_s.DrawLatex(0.17, 0.83, f"Mean #sigma_{{BLUE}} = {sig_c_g.mean():.1f} ps,  peak-to-peak = {sig_c_g.max()-sig_c_g.min():.1f} ps")

save_canvas(c_svx, "v5_sigma_vs_x")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: TOP position leave-one-out validation
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("4. TOP position leave-one-out validation")

# For each x: train calibration x=f(centroid) using all OTHER x positions,
# then predict the held-out x centroid.

def get_centroid(uT, sT, cT, tTop):
    """N_pe-weighted centroid of TOP SiPMs for each event."""
    # We need local_id — reload with local_id
    pass

# We need to reload with local_id for this analysis
def load_x_full(x):
    path = INDIR / f"photon_hits_ej230_ntop20_x{x}.root"
    with uproot.open(path) as f:
        a = f['sipm_hits'].arrays(['event_id', 'face_type', 'local_id', 'time_ns'], library='np')
    ev, ft, lid, tn = a['event_id'], a['face_type'], a['local_id'], a['time_ns']
    return ev, ft, lid, tn


def compute_centroids(ev, ft, lid, tn, x_true):
    """Compute N_pe-weighted centroid for each event. Returns array of centroids."""
    mask = ft == 2
    ev_top = ev[mask]; lid_top = lid[mask]
    ev_uniq = np.unique(ev_top)
    centroids = np.full(len(ev_uniq), np.nan)
    for i, eid in enumerate(ev_uniq):
        m = ev_top == eid
        e_lid = lid_top[m]
        npe = np.bincount(e_lid, minlength=20).astype(float)
        if npe.sum() > 0:
            centroids[i] = np.dot(npe, TOP_XPOS) / npe.sum()
    return ev_uniq, centroids


# Load centroids for all x positions
print("  Loading centroids for all x positions …")
centroids_all = {}
for x in x_scan:
    ev, ft, lid, tn = load_x_full(x)
    ev_uniq, cents = compute_centroids(ev, ft, lid, tn, x)
    centroids_all[x] = (ev_uniq, cents)
    med_c = np.nanmedian(cents)
    sig_c_pos = rsigma(cents[np.isfinite(cents)])
    print(f"  x={x:+4d}: centroid median={med_c:+.1f} mm, σ_centroid={sig_c_pos:.1f} mm")

# Leave-one-out calibration
print("\n  Leave-one-out calibration …")
loo_results = []

for x_hold in x_scan:
    # Training data: all other x positions
    train_x   = [xi for xi in x_scan if xi != x_hold]
    train_med = [np.nanmedian(centroids_all[xi][1]) for xi in train_x]

    # Fit linear calibration: x = a + b * centroid_median
    train_x_arr = np.array(train_x, float)
    train_c_arr = np.array(train_med, float)

    # Polynomial degree 1 calibration (linear)
    ok = np.isfinite(train_c_arr)
    poly1 = np.polyfit(train_c_arr[ok], train_x_arr[ok], 1)
    cal_linear = np.poly1d(poly1)

    # Also try degree 3
    poly3 = np.polyfit(train_c_arr[ok], train_x_arr[ok], 3)
    cal_cubic = np.poly1d(poly3)

    # Predict held-out x
    cents_hold = centroids_all[x_hold][1]
    ok_h = np.isfinite(cents_hold)
    x_pred_lin = cal_linear(cents_hold[ok_h])
    x_pred_cub = cal_cubic(cents_hold[ok_h])

    bias_lin = np.median(x_pred_lin) - x_hold
    bias_cub = np.median(x_pred_cub) - x_hold
    res_lin  = rsigma(x_pred_lin - np.median(x_pred_lin))
    res_cub  = rsigma(x_pred_cub - np.median(x_pred_cub))

    loo_results.append(dict(
        x=x_hold, bias_lin=bias_lin, res_lin=res_lin,
        bias_cub=bias_cub, res_cub=res_cub
    ))
    print(f"  x={x_hold:+4d}: linear bias={bias_lin:+.1f} mm, σ={res_lin:.1f} mm  |"
          f"  cubic bias={bias_cub:+.1f} mm, σ={res_cub:.1f} mm")

xs_loo  = np.array([r['x'] for r in loo_results], dtype='d')
res_lin = np.array([r['res_lin'] for r in loo_results], dtype='d')
res_cub = np.array([r['res_cub'] for r in loo_results], dtype='d')
bia_lin = np.array([r['bias_lin'] for r in loo_results], dtype='d')
print(f"\n  LOO linear: mean σ_x = {res_lin.mean():.1f} mm, max |bias| = {np.abs(bia_lin).max():.1f} mm")
print(f"  LOO cubic:  mean σ_x = {res_cub.mean():.1f} mm")

# FIG v5_top_position_loo
g_lin_loo = make_graph(xs_loo, res_lin)
g_lin_loo.SetMarkerStyle(20); g_lin_loo.SetMarkerSize(1.3)
g_lin_loo.SetMarkerColor(COL_TOP); g_lin_loo.SetLineColor(COL_TOP); g_lin_loo.SetLineWidth(2)

g_cub_loo = make_graph(xs_loo, res_cub)
g_cub_loo.SetMarkerStyle(22); g_cub_loo.SetMarkerSize(1.1)
g_cub_loo.SetMarkerColor(COL_ORA); g_cub_loo.SetLineColor(COL_ORA)
g_cub_loo.SetLineWidth(2); g_cub_loo.SetLineStyle(2)

# END position reference from global results at x=0
sig_x_end_ref = v_eff_mmpns * rsigma(
    (lambda uL,sL,cL,tL,uR,sR,cR,tR,*_: get_T0_end(tL,tR,sL,sR,cL,cR,uL,uR,M_GLOBAL)[2])
    (*load_x(0))
) / 2.
print(f"  END position resolution (m={M_GLOBAL}): σ_x = {sig_x_end_ref:.1f} mm")

c_loo = ROOT.TCanvas("c_loo", "loo", 900, 640)
c_loo.SetLeftMargin(0.15); c_loo.SetBottomMargin(0.15)

mg_loo = ROOT.TMultiGraph()
mg_loo.Add(g_lin_loo, "PL")
mg_loo.Add(g_cub_loo, "PL")
mg_loo.Draw("A")
mg_loo.SetTitle("")
mg_loo.GetXaxis().SetTitle("Muon position  x  [mm]")
mg_loo.GetYaxis().SetTitle("#sigma_{x}  [mm]  (IQR)")
mg_loo.GetXaxis().SetLimits(-720., 720.)
mg_loo.GetYaxis().SetRangeUser(0., 30.)

# END position reference line
end_line = ROOT.TLine(-700., sig_x_end_ref, 700., sig_x_end_ref)
end_line.SetLineColor(COL_END); end_line.SetLineStyle(2); end_line.SetLineWidth(2)
end_line.Draw()

leg_loo = ROOT.TLegend(0.55, 0.60, 0.93, 0.90)
leg_loo.SetTextSize(0.040)
leg_loo.AddEntry(g_lin_loo, f"TOP LOO linear  (mean={res_lin.mean():.1f} mm)", "pl")
leg_loo.AddEntry(g_cub_loo, f"TOP LOO cubic   (mean={res_cub.mean():.1f} mm)", "pl")
leg_loo.AddEntry(end_line,  f"END #sigma_{{x}} = {sig_x_end_ref:.1f} mm (m={M_GLOBAL})", "l")
leg_loo.Draw()

lat_loo = ROOT.TLatex(); lat_loo.SetNDC(); lat_loo.SetTextFont(42); lat_loo.SetTextSize(0.042)
lat_loo.DrawLatex(0.17, 0.88, "TOP N_{pe} centroid — leave-one-out spatial test")
lat_loo.DrawLatex(0.17, 0.83, "EJ-230,  N_{top}=20")
save_canvas(c_loo, "v5_top_position_loo")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: PDE audit — define what 0.607 means
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("5. PDE audit")

path0 = INDIR / "photon_hits_ej230_ntop20_x0.root"
with uproot.open(path0) as f:
    a_pde = f['sipm_hits'].arrays(['face_type', 'pde'], library='np')
ft_pde = a_pde['face_type']
pde_v  = a_pde['pde']

mean_all = pde_v.mean()
n_detected = len(pde_v)
print(f"  Detected photons: {n_detected:,}")
print(f"  Mean PDE (over detected photons): {mean_all:.4f}")
print(f"  NOTE: This is selection-biased — photons are detected BECAUSE their PDE > 0")
print(f"  The correct interpretation: the Geant4 SiPM model applies wavelength-dependent")
print(f"  efficiency; each detected photon has EFFICIENCY (PDE) from the SiPM curve.")
print(f"  The 'mean PDE = {mean_all:.3f}' is the mean efficiency of the accepted population,")
print(f"  NOT the overall fraction of incident photons detected.")
print(f"  For the photon budget, use N_detected / N_incident rather than this mean.")
for ft_id, name in [(0, 'END-L'), (1, 'END-R'), (2, 'TOP')]:
    m = ft_pde == ft_id
    if m.sum() > 0:
        print(f"  {name}: mean_PDE_detected = {pde_v[m].mean():.4f}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: N_TOP scan — use actual physical geometry files
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("6. N_TOP scan at x=0 — physical geometry files")

ntop_results = {}

for N in NTOP_VALS:
    uL, sL, cL, tL, uR, sR, cR, tR, uT, sT, cT, tTop = load_x(0, ntop=N)

    # Scan k (oracle for this physical geometry)
    best_st = 1e9; best_k_n = 1
    for k in range(1, 30):
        ev_t, T0t = get_T0_top(tTop, sT, cT, uT, k)
        if len(T0t) < 30: break
        s = rsigma(T0t - T0t.mean())
        if np.isfinite(s) and s < best_st:
            best_st = s; best_k_n = k

    # BLUE with global m and oracle k for this N
    ev_e, T0e_n, _ = get_T0_end(tL, tR, sL, sR, cL, cR, uL, uR, M_GLOBAL)
    ev_t_n, T0t_n  = get_T0_top(tTop, sT, cT, uT, best_k_n)
    T0e_a, T0t_a   = intersect_estimators(ev_e, T0e_n, ev_t_n, T0t_n)
    bl_n = compute_blue(T0e_a, T0t_a)

    ntop_results[N] = dict(
        k_opt=best_k_n, sig_t=best_st,
        sig_e=rsigma(T0e_a - T0e_a.mean()),
        sig_c=bl_n['sig_c'], w_end=bl_n['w_end'], rho=bl_n['rho'],
        npe_mean=cT.mean() if len(cT) > 0 else 0
    )
    print(f"  N={N:2d}: k*={best_k_n} σ_TOP={best_st:.1f} σ_END={ntop_results[N]['sig_e']:.1f}"
          f" σ_BLUE={bl_n['sig_c']:.1f} ps  w_END={bl_n['w_end']:.3f} ρ={bl_n['rho']:.3f}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG v5_pareto: channel count vs timing — physical configurations only
# ═══════════════════════════════════════════════════════════════════════════════
print("\n  Plotting Pareto frontier …")

# Physical configurations:
# 1. END-only (actual END-only geometry): N=16, sigma ~ 50.5 ps (from best_est results)
# 2. Hybrid: N=16+N_TOP, use results above
# Note: hybrid END sigma is degraded because TOP SiPMs intercept photons

N_END = 16
# Read END-only result from stored CSV if available
end_only_file = REPO / "analysis/optim/root_best_est/output_summary.json"
try:
    import json as json_mod
    with open(end_only_file) as fp:
        eo = json_mod.load(fp)
    sig_end_only = eo.get('sig_end_opt', 50.5)
    print(f"  END-only σ from file: {sig_end_only:.1f} ps")
except:
    sig_end_only = 50.5   # from best_est output
    print(f"  END-only σ (default): {sig_end_only:.1f} ps")

xs_pareto_phys = []
ys_pareto_phys = []
labels_pareto  = []

# END-only point
xs_pareto_phys.append(N_END)
ys_pareto_phys.append(sig_end_only)
labels_pareto.append("END-only\n(N=16, physical)")

for N in NTOP_VALS:
    xs_pareto_phys.append(N_END + N)
    ys_pareto_phys.append(ntop_results[N]['sig_c'])
    labels_pareto.append(f"BLUE N_{{TOP}}={N}")

xs_p = np.array(xs_pareto_phys, dtype='d')
ys_p = np.array(ys_pareto_phys, dtype='d')

g_pareto_phys = make_graph(xs_p, ys_p)
g_pareto_phys.SetMarkerStyle(20); g_pareto_phys.SetMarkerSize(1.5)
g_pareto_phys.SetMarkerColor(COL_BLUE); g_pareto_phys.SetLineColor(COL_BLUE); g_pareto_phys.SetLineWidth(2)

c_par = ROOT.TCanvas("c_pareto_v5", "pareto_v5", 900, 640)
c_par.SetLeftMargin(0.15); c_par.SetBottomMargin(0.15)

g_pareto_phys.SetTitle("")
g_pareto_phys.GetXaxis().SetTitle("Total SiPM channels per bar  (N_{END} + N_{TOP})")
g_pareto_phys.GetYaxis().SetTitle("#sigma_{IQR}  [ps]")
g_pareto_phys.GetXaxis().SetLimits(10., 42.)
g_pareto_phys.GetYaxis().SetRangeUser(0., 65.)
g_pareto_phys.Draw("APL")

# Add point labels
lat_p = ROOT.TLatex(); lat_p.SetTextFont(42); lat_p.SetTextSize(0.038)
for xi, yi, lbl in zip(xs_p, ys_p, labels_pareto):
    first_line = lbl.split('\n')[0]
    offset_x = 0.3; offset_y = 1.5
    lat_p.DrawLatex(xi + offset_x, yi + offset_y, first_line)

# Mark END-only with a different marker
g_end_only_pt = ROOT.TGraph(1, np.array([float(N_END)]), np.array([sig_end_only]))
g_end_only_pt.SetMarkerStyle(22); g_end_only_pt.SetMarkerSize(1.8)
g_end_only_pt.SetMarkerColor(COL_END)
g_end_only_pt.Draw("P same")

leg_par = ROOT.TLegend(0.55, 0.60, 0.93, 0.88)
leg_par.SetTextSize(0.038)
leg_par.AddEntry(g_end_only_pt, "END-only geometry (physical)", "p")
leg_par.AddEntry(g_pareto_phys, "Hybrid BLUE (END+TOP physical)", "pl")
leg_par.Draw()

lat_par = ROOT.TLatex(); lat_par.SetNDC(); lat_par.SetTextFont(42); lat_par.SetTextSize(0.040)
lat_par.DrawLatex(0.17, 0.88, "EJ-230,  x=0,  global k*/m*")
lat_par.DrawLatex(0.17, 0.83, "Global k(N_{TOP}) from oracle scan,  m=" + str(M_GLOBAL))
save_canvas(c_par, "v5_pareto")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("Summary numbers for FINAL_NUMBERS.md")

g_x0 = global_results[0]
o_x0 = oracle_results[0]

print(f"\nGlobal fixed (k={K_GLOBAL}, m={M_GLOBAL}) at x=0:")
print(f"  σ_END   = {g_x0['sig_e']:.2f} ps")
print(f"  σ_TOP   = {g_x0['sig_t']:.2f} ps")
print(f"  σ_BLUE  = {g_x0['sig_c']:.2f} ± {g_x0['boot_c']:.2f} ps")
print(f"  w_END   = {g_x0['w_end']:.4f}")
print(f"  ρ(E,T)  = {g_x0['rho']:.4f}")

print(f"\nOracle (local opt k,m) at x=0:")
print(f"  k* = {o_x0['k']}, m* = {o_x0['m']}")
print(f"  σ_BLUE_oracle = {o_x0['sig_c']:.2f} ps")
print(f"  Oracle - global bias = {o_x0['sig_c'] - g_x0['sig_c']:+.2f} ps")

print(f"\nGlobal σ_BLUE across x:")
print(f"  Mean     = {sig_c_g.mean():.1f} ps")
print(f"  Min      = {sig_c_g.min():.1f} ps")
print(f"  Max      = {sig_c_g.max():.1f} ps")
print(f"  Pk-to-pk = {sig_c_g.max()-sig_c_g.min():.1f} ps")

print(f"\nOracle σ_BLUE across x:")
print(f"  Mean     = {sig_c_o.mean():.1f} ps")
print(f"  Min      = {sig_c_o.min():.1f} ps")
print(f"  Max      = {sig_c_o.max():.1f} ps")

print(f"\nv_eff = {v_eff_mmns:.1f} ± {v_eff_err:.1f} mm/ns (stat. only)")
print(f"  χ²/ndf = {chi2_lin:.0f}/{ndf_lin}  Nonlinearity residual RMS = {rms_res:.1f} ps")

print(f"\nN_TOP scan (σ_BLUE at x=0):")
for N in NTOP_VALS:
    r = ntop_results[N]
    print(f"  N_TOP={N:2d}: k*={r['k_opt']} σ_TOP={r['sig_t']:.1f} σ_BLUE={r['sig_c']:.1f} ps  w_END={r['w_end']:.3f}")

print(f"\nEND-only (physical geometry): σ = {sig_end_only:.1f} ps")

print(f"\nTOP position leave-one-out (N_pe centroid):")
print(f"  Linear cal:  mean σ_x = {res_lin.mean():.1f} mm, max|bias| = {np.abs(bia_lin).max():.1f} mm")
print(f"  Cubic cal:   mean σ_x = {res_cub.mean():.1f} mm")

fout.Close()
print("\nDone. All figures saved to", FIGS)
