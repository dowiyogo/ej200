// new_analysis_plots.C — supplemental analysis figures for napkin presentation
// Generates from hardcoded verified data (verification.json + phase_ab_kscan_full.csv)
// Run: root -l -b -q new_analysis_plots.C
// Output: fig_end_mscan.pdf, fig_ntop_scan_full.pdf, fig_blue_weights.pdf,
//         fig_sigma_vs_x_new.pdf, fig_pareto.pdf

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TArrow.h"
#include <cmath>

void SetStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetTitleSize(0.055,"xyz");
    gStyle->SetTitleOffset(1.15,"y");
    gStyle->SetTitleOffset(1.0,"x");
    gStyle->SetTickLength(0.03,"xy");
}

// m-scan: σ_END(m) for EJ-230, x=0, END-only Vikuiti geometry (from phase_ab_kscan_full.csv)
// n_events = 20000 per point
static const int N_MSCAN = 50;
static const double mscan_m[N_MSCAN] = {
     1, 2, 3, 4, 5, 6, 7, 8, 9,10,
    11,12,13,14,15,16,17,18,19,20,
    21,22,23,24,25,26,27,28,29,30,
    31,32,33,34,35,36,37,38,39,40,
    41,42,43,44,45,46,47,48,49,50};
static const double mscan_s[N_MSCAN] = {
    59.122,53.566,52.195,51.220,50.713,50.848,50.494,50.822,51.196,51.214,
    51.546,51.897,52.290,52.443,52.648,53.185,53.353,53.370,53.557,54.044,
    54.259,54.131,54.279,54.532,54.947,55.649,56.430,57.003,57.396,58.177,
    58.451,58.872,59.404,59.959,59.955,60.221,60.926,61.476,61.886,62.308,
    62.497,63.065,63.320,63.552,63.934,64.424,64.731,64.938,65.135,65.314};

// ─── Figure 1: σ_END(m) vs m ─────────────────────────────────────────────────
void FigEndMscan() {
    double n_ev = 20000.;
    TGraphErrors* g = new TGraphErrors(N_MSCAN);
    for (int i = 0; i < N_MSCAN; ++i) {
        g->SetPoint(i, mscan_m[i], mscan_s[i]);
        g->SetPointError(i, 0., mscan_s[i]/std::sqrt(2.*n_ev));
    }
    g->SetLineColor(kRed+1); g->SetMarkerColor(kRed+1);
    g->SetMarkerStyle(20); g->SetMarkerSize(0.8);
    g->SetLineWidth(2);

    TCanvas* c = new TCanvas("cmscan","",800,580);
    g->Draw("APL");
    g->GetXaxis()->SetTitle("m  (number of END photons averaged per face)");
    g->GetYaxis()->SetTitle("#sigma_{END}(m)  [ps]");
    g->GetXaxis()->SetLimits(0, 51);
    g->GetXaxis()->SetNdivisions(510);
    g->GetYaxis()->SetRangeUser(45, 75);
    g->GetYaxis()->SetNdivisions(506);

    // mark m* = 7 (minimum in END-only geometry)
    double sig_star = 50.494;
    TLine* lv = new TLine(7, 45, 7, sig_star);
    lv->SetLineStyle(2); lv->SetLineColor(kBlue+1); lv->SetLineWidth(2); lv->Draw();
    TLine* lh = new TLine(0, sig_star, 7, sig_star);
    lh->SetLineStyle(2); lh->SetLineColor(kBlue+1); lh->SetLineWidth(2); lh->Draw();

    // mark m*=8 in combined geometry (slightly worse: 53.7 ps, due to photon loss to TOP)
    double sig_comb = 53.68;
    TLine* lv2 = new TLine(8, 45, 8, sig_comb);
    lv2->SetLineStyle(3); lv2->SetLineColor(kGreen+2); lv2->SetLineWidth(2); lv2->Draw();
    TLine* lh2 = new TLine(0, sig_comb, 8, sig_comb);
    lh2->SetLineStyle(3); lh2->SetLineColor(kGreen+2); lh2->SetLineWidth(2); lh2->Draw();

    TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.043);
    lat.SetTextColor(kBlue+1);
    lat.DrawLatex(8.5, sig_star+0.5, Form("m*=7,  #sigma=%.1f ps  (END-only geometry)", sig_star));
    lat.SetTextColor(kGreen+2);
    lat.DrawLatex(9.5, sig_comb+0.5, Form("m*=8,  #sigma=%.1f ps  (combined geometry, photon loss to TOP)", sig_comb));

    lat.SetTextColor(kGray+1); lat.SetTextSize(0.038);
    lat.DrawLatex(32, 72.5, "EJ-230,  x = 0,  Vikuiti reflector");

    c->Print("fig_end_mscan.pdf");
    delete c;
    printf("  fig_end_mscan.pdf\n");
}

// ─── Figure 2: σ_END/TOP/BLUE vs N_top with BLUE weights panel ───────────────
void FigNtopScanFull() {
    // verified data: ntop_scan_ej230_x0
    const int NN = 4;
    double ntop_v[NN]    = {4.,    8.,    14.,   20.};
    double sig_end_v[NN] = {50.96, 50.21, 51.04, 53.68};
    double sig_top_v[NN] = {78.20, 32.00, 18.67, 15.20};
    double sig_comb_v[NN]= {45.47, 28.92, 18.16, 15.21};
    double boot_t_v[NN]  = {1.80,  0.56,  0.34,  0.26};
    double boot_c_v[NN]  = {0.82,  0.48,  0.33,  0.27};
    double w_end_v[NN]   = {0.741, 0.212, 0.024, 0.013};
    double rho_v[NN]     = {0.149, 0.180, 0.182, 0.159};
    double npe_top_v[NN] = {287.,  565.,  962.,  1341.};

    TGraphErrors *g_top  = new TGraphErrors(NN);
    TGraphErrors *g_end  = new TGraphErrors(NN);
    TGraphErrors *g_comb = new TGraphErrors(NN);
    TGraph       *g_wend = new TGraph(NN);
    TGraph       *g_rho  = new TGraph(NN);

    for (int i = 0; i < NN; ++i) {
        g_end ->SetPoint(i, ntop_v[i], sig_end_v[i]);  g_end ->SetPointError(i,0,0);
        g_top ->SetPoint(i, ntop_v[i], sig_top_v[i]);  g_top ->SetPointError(i,0,boot_t_v[i]);
        g_comb->SetPoint(i, ntop_v[i], sig_comb_v[i]); g_comb->SetPointError(i,0,boot_c_v[i]);
        g_wend->SetPoint(i, ntop_v[i], w_end_v[i]);
        g_rho ->SetPoint(i, ntop_v[i], rho_v[i]);
    }

    auto StyleG = [](TGraph* g, Color_t col, Style_t ms, int lw=2) {
        g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetMarkerStyle(ms); g->SetMarkerSize(1.2); g->SetLineWidth(lw);
    };
    StyleG(g_end,  kBlue+1,  21);
    StyleG(g_top,  kRed+1,   22);
    StyleG(g_comb, kGreen+2, 20);
    StyleG(g_wend, kViolet+1, 24);
    StyleG(g_rho,  kOrange+1, 25);

    // Two-panel canvas
    TCanvas* c = new TCanvas("cntop","",900,700);
    c->Divide(1,2,0,0);

    // Top pad: sigma vs N_top
    c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.01);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(g_end,  "PL"); mg->Add(g_top,  "PL"); mg->Add(g_comb, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("");
    mg->GetYaxis()->SetTitle("#sigma(T_{0})  [ps]");
    mg->GetXaxis()->SetLimits(2, 22);
    mg->GetYaxis()->SetRangeUser(0, 95);
    mg->GetXaxis()->SetNdivisions(504);
    mg->GetYaxis()->SetNdivisions(506);
    gPad->Modified(); gPad->Update();

    TLegend* leg = new TLegend(0.55,0.55,0.93,0.88);
    leg->AddEntry(g_end,  "#sigma_{END}   (m*-average, 16 SiPMs)",  "pl");
    leg->AddEntry(g_top,  "#sigma_{TOP}   (k*-th photon, N SiPMs)", "pl");
    leg->AddEntry(g_comb, "#sigma_{BLUE}  (BLUE combination)",       "pl");
    leg->SetTextSize(0.052); leg->Draw();

    TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.05);
    lat.SetNDC(); lat.DrawLatex(0.15, 0.88, "EJ-230,  x = 0,  sparse TOP + END Vikuiti");

    // Bottom pad: w_end and rho vs N_top
    c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.01); gPad->SetBottomMargin(0.18);

    TMultiGraph* mg2 = new TMultiGraph();
    mg2->Add(g_wend, "PL"); mg2->Add(g_rho, "PL");
    mg2->Draw("A");
    mg2->GetXaxis()->SetTitle("N_{top}  (number of TOP SiPMs)");
    mg2->GetYaxis()->SetTitle("weight / correlation");
    mg2->GetXaxis()->SetLimits(2, 22);
    mg2->GetYaxis()->SetRangeUser(0, 1.0);
    mg2->GetXaxis()->SetNdivisions(504);
    mg2->GetYaxis()->SetNdivisions(505);
    gPad->Modified(); gPad->Update();

    TLegend* leg2 = new TLegend(0.40,0.55,0.93,0.88);
    leg2->AddEntry(g_wend, "w_{END}  (BLUE weight on END)",      "pl");
    leg2->AddEntry(g_rho,  "#rho(END,TOP)  (empirical correlation)", "pl");
    leg2->SetTextSize(0.06); leg2->Draw();

    c->Print("fig_ntop_scan_full.pdf");
    delete c;
    printf("  fig_ntop_scan_full.pdf\n");
}

// ─── Figure 3: σ_END/TOP/BLUE vs x (position scan, N=20) ────────────────────
void FigSigmaVsX() {
    // verified data: pos_scan_ej230_n20
    const int NX = 13;
    double x_v[NX]   = {-600,-500,-400,-300,-200,-100,  0, 100, 200, 300, 400, 500, 600};
    double se_v[NX]  = {59.02,60.09,55.42,56.27,53.56,52.83,53.68,53.05,54.53,54.95,57.36,58.81,60.87};
    double st_v[NX]  = {15.42,17.80,15.92,15.95,17.64,14.91,15.20,15.62,18.20,15.85,15.78,17.46,15.28};
    double sc_v[NX]  = {15.06,17.41,15.36,15.72,17.24,14.64,15.21,14.97,17.79,15.57,15.33,16.93,15.05};
    double bc_v[NX]  = {0.30, 0.29, 0.26, 0.29, 0.31, 0.27, 0.26, 0.23, 0.35, 0.31, 0.24, 0.26, 0.27};

    TGraph       *g_end  = new TGraph(NX, x_v, se_v);
    TGraphErrors *g_top  = new TGraphErrors(NX);
    TGraphErrors *g_comb = new TGraphErrors(NX);
    for (int i = 0; i < NX; ++i) {
        g_top ->SetPoint(i, x_v[i], st_v[i]); g_top ->SetPointError(i,0,0);
        g_comb->SetPoint(i, x_v[i], sc_v[i]); g_comb->SetPointError(i,0,bc_v[i]);
    }

    auto StyleG = [](TGraph* g, Color_t col, Style_t ms, int lw=2) {
        g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetMarkerStyle(ms); g->SetMarkerSize(1.0); g->SetLineWidth(lw);
    };
    StyleG(g_end,  kBlue+1, 21);
    StyleG(g_top,  kRed+1,  22);
    StyleG(g_comb, kGreen+2, 20);

    TCanvas* c = new TCanvas("csigx","",800,580);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(g_end,  "PL"); mg->Add(g_top, "PL"); mg->Add(g_comb, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("Muon position  x  [mm]");
    mg->GetYaxis()->SetTitle("#sigma(T_{0})  [ps]");
    mg->GetXaxis()->SetLimits(-650, 650);
    mg->GetYaxis()->SetRangeUser(0, 75);
    mg->GetXaxis()->SetNdivisions(506);
    mg->GetYaxis()->SetNdivisions(507);
    gPad->Modified(); gPad->Update();

    TLegend* leg = new TLegend(0.18, 0.65, 0.60, 0.90);
    leg->AddEntry(g_end,  "#sigma_{END}   (m*-average, 16 SiPMs)", "pl");
    leg->AddEntry(g_top,  "#sigma_{TOP}   (k*-th photon, 20 SiPMs)","pl");
    leg->AddEntry(g_comb, "#sigma_{BLUE}  (BLUE combination)","pl");
    leg->SetTextSize(0.042); leg->Draw();

    TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.040); lat.SetNDC();
    lat.DrawLatex(0.18, 0.59, "EJ-230,  N_{top}=20,  Vikuiti reflector");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.55, 0.20, Form("#sigma_{BLUE}: %.1f#font[122]{-}%.1f ps  (%.0f%% peak-to-peak)",
        14.64, 17.79, (17.79-14.64)/14.64*100.));

    c->Print("fig_sigma_vs_x_new.pdf");
    delete c;
    printf("  fig_sigma_vs_x_new.pdf\n");
}

// ─── Figure 4: Pareto frontier ────────────────────────────────────────────────
// Channel accounting (AUDITED):
//   phase_e_pareto N_sipm for END configs = TOTAL SiPMs (both ends combined):
//     N_sipm=2  →  1 per end ×2 = 2 total
//     N_sipm=16 →  8 per end ×2 = 16 total   ← standard design
//   BLUE series: N_total = N_top + 16 (END fixed at 8/end × 2)
//   TOP TIR/VIK: N_sipm = total TOP SiPMs (phase_e_pareto)
void FigPareto() {
    // END+Vik EJ-230 ablation (k*-th photon estimator from END face, both ends)
    // N_sipm values are TOTAL (directly from phase_e_pareto, already both ends combined)
    const int NE = 6;
    double ne_x[NE] = { 2,  4,  8, 12, 14, 16};
    double ne_s[NE] = {131.24, 97.90, 73.40, 63.20, 59.22, 56.04};
    double ne_e[NE] = {  0.99,  1.56,  0.84,  0.29,  0.47,  0.00};

    // m-average result for 16 SiPMs total (8/end) in combined geometry
    // END-only: 50.49 ps (m*=7); combined (with TOP): 53.68±0.93 ps (m*=8)
    const int NM = 1;
    double nm_x[NM] = {16.};
    double nm_s[NM] = {53.68};
    double nm_e[NM] = { 0.93};

    // TOP + TIR reference: N_sipm = total TOP SiPMs, per-SiPM estimator
    const int NT = 2;
    double nt_x[NT] = {50., 70.};
    double nt_s[NT] = {81.97, 70.22};
    double nt_e[NT] = {14.78,  0.00};

    // BLUE(END+TOP) sparse Vikuiti: N_total = 16 END + N_top
    const int NB = 4;
    double nb_x[NB] = {20., 24., 30., 36.}; // 16+4, 16+8, 16+14, 16+20
    double nb_s[NB] = {45.47, 28.92, 18.16, 15.21};
    double nb_e[NB] = { 0.82,  0.48,  0.33,  0.26};

    auto MkGE = [](int n, double* x, double* s, double* e) {
        TGraphErrors* g = new TGraphErrors(n, x, s, nullptr, e);
        return g;
    };
    TGraphErrors *gEnd  = MkGE(NE, ne_x, ne_s, ne_e);
    TGraphErrors *gMAvg = MkGE(NM, nm_x, nm_s, nm_e);
    TGraphErrors *gTIR  = MkGE(NT, nt_x, nt_s, nt_e);
    TGraphErrors *gBLUE = MkGE(NB, nb_x, nb_s, nb_e);

    auto StyleG = [](TGraph* g, Color_t col, Style_t ms, int lw=2) {
        g->SetLineColor(col); g->SetMarkerColor(col);
        g->SetMarkerStyle(ms); g->SetMarkerSize(1.2); g->SetLineWidth(lw);
    };
    StyleG(gEnd,  kBlue+1, 20);
    StyleG(gMAvg, kBlue+3, 29, 0); // open star, no line (single point)
    StyleG(gTIR,  kRed+1,  22);
    StyleG(gBLUE, kGreen+2,21);

    gEnd ->Sort(); gTIR ->Sort(); gBLUE->Sort();

    TCanvas* c = new TCanvas("cpareto","",900,620);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.06);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gEnd,  "PL");
    mg->Add(gMAvg, "P");
    mg->Add(gTIR,  "PL");
    mg->Add(gBLUE, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("Total SiPM count per bar  N_{SiPM}");
    mg->GetYaxis()->SetTitle("#sigma(T_{0})  [ps]");
    mg->GetXaxis()->SetLimits(0, 80);
    mg->GetXaxis()->SetNdivisions(508);
    gPad->Modified(); gPad->Update();

    // Vertical guide at N_END = 16 (standard END design)
    TLine* lv = new TLine(16, mg->GetYaxis()->GetXmin(), 16, 115);
    lv->SetLineStyle(3); lv->SetLineColor(kGray+1); lv->SetLineWidth(1); lv->Draw();
    TLatex annot; annot.SetTextFont(42); annot.SetTextSize(0.032); annot.SetTextAngle(90);
    annot.SetTextColor(kGray+1);
    annot.DrawLatex(17., 25., "16 SiPMs (standard END)");

    TLegend* leg = new TLegend(0.40, 0.58, 0.93, 0.90);
    leg->AddEntry(gEnd,  "END+Vik EJ-230  (k*-th order stat, ablation)",  "pl");
    leg->AddEntry(gMAvg, "END+Vik EJ-230  (m*-average, combined geom.)", "p");
    leg->AddEntry(gTIR,  "TOP TIR EJ-230  (per-SiPM estimator, ref.)",   "pl");
    leg->AddEntry(gBLUE, "BLUE(END+TOP)   EJ-230+Vik  (sparse TOP)",     "pl");
    leg->SetTextSize(0.036); leg->Draw();

    TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.038);
    lat.SetTextColor(kGreen+2);
    lat.DrawLatex(37., 12.8, "15.2 ps  (N_{top}=20,  36 total)");

    c->Print("fig_pareto.pdf");
    delete c;
    printf("  fig_pareto.pdf\n");
}

// ─── Main ─────────────────────────────────────────────────────────────────────
void new_analysis_plots() {
    SetStyle();
    gROOT->SetBatch(kTRUE);
    FigEndMscan();
    FigNtopScanFull();
    FigSigmaVsX();
    FigPareto();
    printf("\nDone. Four new figures written.\n");
}
