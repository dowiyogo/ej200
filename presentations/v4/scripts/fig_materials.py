#!/usr/bin/env python3
"""
fig_materials.py — Material comparison (EJ-200, EJ-204, EJ-230) from END-only scan.
Reads phase_ab_optimal.csv and generates ROOT figures for the presentation.
"""
import ROOT, numpy as np, pandas as pd
from pathlib import Path

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

REPO  = Path("/home/rrios/ej200")
FIGS  = REPO / "analysis/presentation_v4/figs"
TABS  = REPO / "analysis/presentation_v4/tables"

def set_style():
    s = ROOT.gStyle
    s.SetOptStat(0); s.SetOptFit(0)
    s.SetCanvasColor(0); s.SetPadColor(0); s.SetFrameFillColor(0)
    s.SetCanvasBorderMode(0); s.SetPadBorderMode(0); s.SetFrameBorderMode(0)
    for ax in ("x","y","z"):
        s.SetTitleFont(42,ax); s.SetLabelFont(42,ax)
        s.SetTitleSize(0.055,ax); s.SetLabelSize(0.048,ax)
        s.SetTitleOffset(1.0,ax)
    s.SetTextFont(42); s.SetTextSize(0.05)
    s.SetPadTopMargin(0.07); s.SetPadRightMargin(0.22)
    s.SetPadBottomMargin(0.14); s.SetPadLeftMargin(0.14)
    s.SetTitleOffset(1.15,"y")
    s.SetLegendBorderSize(1); s.SetLegendFont(42); s.SetLegendFillColor(0)
    s.SetHistLineWidth(2); s.SetPadTickX(1); s.SetPadTickY(1)
    s.SetMarkerSize(1.3)

set_style()

def make_graph(x, y, ey=None):
    x = np.asarray(x,'d'); y = np.asarray(y,'d')
    n = len(x)
    if ey is not None:
        g = ROOT.TGraphErrors(n)
        ey = np.asarray(ey,'d')
        for i in range(n): g.SetPoint(i,x[i],y[i]); g.SetPointError(i,0.,ey[i])
    else:
        g = ROOT.TGraph(n)
        for i in range(n): g.SetPoint(i,x[i],y[i])
    return g

df = pd.read_csv(REPO/"analysis/optim/phase_ab_optimal.csv")

mats = ['ej200','ej204','ej230']
labels  = {'ej200':'EJ-200 (Y=10000/MeV, #tau=2.1 ns)',
           'ej204':'EJ-204 (Y=10400/MeV, #tau=1.8 ns)',
           'ej230':'EJ-230 (Y=9700/MeV, #tau=1.5 ns)'}
colors  = {'ej200':ROOT.kAzure+2, 'ej204':ROOT.kOrange+1, 'ej230':ROOT.kGreen+3}
markers = {'ej200':20, 'ej204':21, 'ej230':22}
lstyle  = {'ej200':1, 'ej204':2, 'ej230':3}

# ── FIGURE M1: sigma_T_END (m-opt) vs x for three materials ──
c1 = ROOT.TCanvas("c_mat_end","mat end",900,620)
c1.SetLeftMargin(0.14); c1.SetBottomMargin(0.14); c1.SetRightMargin(0.23)

mg = ROOT.TMultiGraph()
graphs = {}
for sc in mats:
    sub = df[df.scint==sc].sort_values('x')
    xs = sub.x.values.astype('d')
    ys = sub.sig_m_opt.values.astype('d')
    g = make_graph(xs, ys)
    g.SetMarkerStyle(markers[sc]); g.SetMarkerSize(1.3)
    g.SetMarkerColor(colors[sc]); g.SetLineColor(colors[sc])
    g.SetLineWidth(2); g.SetLineStyle(lstyle[sc])
    mg.Add(g,"PL")
    graphs[sc] = g

mg.Draw("A")
mg.SetTitle("")
mg.GetXaxis().SetTitle("Muon position x  [mm]")
mg.GetYaxis().SetTitle("#sigma_{T}^{END,m*}  [ps]  (IQR)")
mg.GetXaxis().SetLimits(-700.,700.)
mg.GetYaxis().SetRangeUser(40.,75.)

# Legend placed in right margin
leg1 = ROOT.TLegend(0.78, 0.45, 0.99, 0.88)
leg1.SetTextSize(0.036); leg1.SetBorderSize(1)
for sc in mats:
    # get x=0 sigma
    s0 = df[(df.scint==sc)&(df.x==0)].sig_m_opt.values[0]
    m0 = int(df[(df.scint==sc)&(df.x==0)].m_opt.values[0])
    leg1.AddEntry(graphs[sc], f"{labels[sc].split('(')[0].strip()}", "pl")
    leg1.AddEntry(ROOT.TObject(),"", "")
leg1.SetNColumns(1)
# Simpler legend
leg1.Clear()
for sc in mats:
    s0 = df[(df.scint==sc)&(df.x==0)].sig_m_opt.values[0]
    m0 = int(df[(df.scint==sc)&(df.x==0)].m_opt.values[0])
    leg1.AddEntry(graphs[sc], f"{sc.upper()}: {s0:.1f} ps (x=0)", "pl")
leg1.Draw()

lat1 = ROOT.TLatex(); lat1.SetNDC(); lat1.SetTextFont(42); lat1.SetTextSize(0.040)
lat1.DrawLatex(0.16, 0.88, "END + Vikuiti,  m-optimal estimator")
lat1.DrawLatex(0.16, 0.83, "END-only geometry  (no TOP SiPMs)")

c1.SaveAs(str(FIGS/"figM1_mat_sigma_end.pdf"))
c1.SaveAs(str(FIGS/"figM1_mat_sigma_end.png"))
print("Saved figM1_mat_sigma_end.{pdf,png}")

# ── FIGURE M2: Npe vs x for three materials ──
c2 = ROOT.TCanvas("c_mat_npe","mat npe",900,620)
c2.SetLeftMargin(0.14); c2.SetBottomMargin(0.14); c2.SetRightMargin(0.23)

mg2 = ROOT.TMultiGraph()
graphs2 = {}
for sc in mats:
    sub = df[df.scint==sc].sort_values('x')
    xs = sub.x.values.astype('d')
    # Average L and R Npe
    npe_avg = ((sub.npe_L + sub.npe_R)/2.).values.astype('d')
    g = make_graph(xs, npe_avg)
    g.SetMarkerStyle(markers[sc]); g.SetMarkerSize(1.3)
    g.SetMarkerColor(colors[sc]); g.SetLineColor(colors[sc])
    g.SetLineWidth(2); g.SetLineStyle(lstyle[sc])
    mg2.Add(g,"PL")
    graphs2[sc] = g

mg2.Draw("A")
mg2.SetTitle("")
mg2.GetXaxis().SetTitle("Muon position x  [mm]")
mg2.GetYaxis().SetTitle("N_{pe}/end  (mean detected photons per face)")
mg2.GetXaxis().SetLimits(-700.,700.)

leg2 = ROOT.TLegend(0.78, 0.45, 0.99, 0.88)
leg2.SetTextSize(0.038); leg2.SetBorderSize(1)
for sc in mats:
    n0 = df[(df.scint==sc)&(df.x==0)][['npe_L','npe_R']].values.mean()
    leg2.AddEntry(graphs2[sc], f"{sc.upper()}: {n0:.0f} pe (x=0)", "pl")
leg2.Draw()

lat2 = ROOT.TLatex(); lat2.SetNDC(); lat2.SetTextFont(42); lat2.SetTextSize(0.040)
lat2.DrawLatex(0.16, 0.88, "END + Vikuiti")
lat2.DrawLatex(0.16, 0.83, "More photons #rightarrow better SN, but slower scintillation")

c2.SaveAs(str(FIGS/"figM2_mat_npe.pdf"))
c2.SaveAs(str(FIGS/"figM2_mat_npe.png"))
print("Saved figM2_mat_npe.{pdf,png}")

# ── FIGURE M3: Two-panel: sigma and Npe vs x for EJ-230 only ──
c3 = ROOT.TCanvas("c_ej230_endonly","ej230 end",900,620)
c3.SetLeftMargin(0.14); c3.SetBottomMargin(0.14); c3.SetRightMargin(0.06)

sub230 = df[df.scint=='ej230'].sort_values('x')
xs230 = sub230.x.values.astype('d')
ys230 = sub230.sig_m_opt.values.astype('d')
npe230 = ((sub230.npe_L + sub230.npe_R)/2.).values.astype('d')

g230s = make_graph(xs230, ys230)
g230s.SetMarkerStyle(22); g230s.SetMarkerSize(1.3)
g230s.SetMarkerColor(ROOT.kGreen+3); g230s.SetLineColor(ROOT.kGreen+3); g230s.SetLineWidth(2)
g230s.SetTitle("");
g230s.GetXaxis().SetTitle("Muon position x  [mm]")
g230s.GetYaxis().SetTitle("#sigma_{T}^{END}  [ps]")
g230s.GetXaxis().SetLimits(-700.,700.)
g230s.Draw("APL")

lat3 = ROOT.TLatex(); lat3.SetNDC(); lat3.SetTextFont(42); lat3.SetTextSize(0.042)
lat3.DrawLatex(0.55, 0.88, "EJ-230 END-only")

c3.SaveAs(str(FIGS/"figM3_ej230_endonly.pdf"))
c3.SaveAs(str(FIGS/"figM3_ej230_endonly.png"))
print("Saved figM3_ej230_endonly.{pdf,png}")

print("\nDone. Material comparison figures saved.")
