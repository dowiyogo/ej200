// fig_gen.C  —  Generate figures for napkin_first_principles deck using ROOT
// Run: root -l -b -q fig_gen.C
// Outputs: figs/fig_npe_x.pdf  fig_sigma_t_x.pdf  fig_survival_RN.pdf  fig_bulk_survival.pdf

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLatex.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

// ── Style ──────────────────────────────────────────────────────────────────
void SetStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetHistLineWidth(2);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetTitleSize(0.055,"xyz");
    gStyle->SetTitleOffset(1.1,"y");
    gStyle->SetTickLength(0.03,"xy");
}

// ── Load CSV ───────────────────────────────────────────────────────────────
struct Row { std::string scint; double x,k_opt,sig_k_opt,m_opt,sig_m_opt,npe_L,npe_R; int n_ev; };

std::vector<Row> LoadCSV(const char* path) {
    std::vector<Row> rows;
    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        Row r;
        std::string tok;
        std::getline(ss,tok,','); r.scint = tok;
        std::getline(ss,tok,','); r.x          = std::stod(tok);
        std::getline(ss,tok,','); r.k_opt      = std::stod(tok);
        std::getline(ss,tok,','); r.sig_k_opt  = std::stod(tok);
        std::getline(ss,tok,','); r.m_opt      = std::stod(tok);
        std::getline(ss,tok,','); r.sig_m_opt  = std::stod(tok);
        std::getline(ss,tok,','); r.npe_L      = std::stod(tok);
        std::getline(ss,tok,','); r.npe_R      = std::stod(tok);
        std::getline(ss,tok,','); r.n_ev       = std::stoi(tok);
        rows.push_back(r);
    }
    return rows;
}

// ── Build per-material TGraphErrors ───────────────────────────────────────
// npe_err = sqrt(mean_npe / n_ev)  (Poisson approx on mean)
// sig_err = sigma / sqrt(2*n_ev)   (chi distribution approx)
void BuildGraphs(const std::vector<Row>& rows,
                 const std::string& mat,
                 TGraphErrors*& g_npe,
                 TGraphErrors*& g_sig)
{
    g_npe = new TGraphErrors(); g_npe->SetName(("npe_"+mat).c_str());
    g_sig = new TGraphErrors(); g_sig->SetName(("sig_"+mat).c_str());
    int i = 0;
    for (const auto& r : rows) {
        if (r.scint != mat) continue;
        double npe  = 0.5*(r.npe_L + r.npe_R);
        double enpe = std::sqrt(npe / r.n_ev);
        double sig  = r.sig_m_opt;
        double esig = sig / std::sqrt(2.0 * r.n_ev);
        g_npe->SetPoint(i, r.x, npe);  g_npe->SetPointError(i, 0, enpe);
        g_sig->SetPoint(i, r.x, sig);  g_sig->SetPointError(i, 0, esig);
        ++i;
    }
}

void StyleGraph(TGraphErrors* g, Color_t col, Style_t mstyle=20) {
    g->SetLineColor(col);  g->SetMarkerColor(col);
    g->SetMarkerStyle(mstyle); g->SetMarkerSize(1.1);
    g->SetLineWidth(2);
}

// ══════════════════════════════════════════════════════════════════════════
void fig_gen() {
    SetStyle();
    gROOT->SetBatch(kTRUE);

    const char* CSV = "/home/rrios/ej200/analysis/optim/phase_ab_optimal.csv";
    auto rows = LoadCSV(CSV);

    // material labels and colours (ROOT palette)
    const char* mats[3] = {"ej200","ej204","ej230"};
    const char* labs[3] = {"EJ-200  (#lambda = 380 cm)","EJ-204  (#lambda = 160 cm)","EJ-230  (#lambda = 120 cm)"};
    Color_t     cols[3] = {kBlue+1, kGreen+2, kRed+1};
    Style_t     mark[3] = {20, 21, 22};

    TGraphErrors *gNpe[3], *gSig[3];
    for (int m=0; m<3; ++m) BuildGraphs(rows, mats[m], gNpe[m], gSig[m]);
    for (int m=0; m<3; ++m) { StyleGraph(gNpe[m],cols[m],mark[m]); StyleGraph(gSig[m],cols[m],mark[m]); }

    // ── Figure 1: Npe/end vs x ────────────────────────────────────────────
    {
        TCanvas* c = new TCanvas("c1","",700,500);
        TMultiGraph* mg = new TMultiGraph();
        for (int m=0; m<3; ++m) mg->Add(gNpe[m],"PL");
        mg->Draw("A");
        mg->GetXaxis()->SetTitle("Muon position  x  [mm]");
        mg->GetYaxis()->SetTitle("N_{pe} per END face");
        mg->GetXaxis()->SetLimits(-650, 650);
        mg->GetYaxis()->SetRangeUser(0, 1800);
        mg->GetXaxis()->SetNdivisions(506);

        TLegend* leg = new TLegend(0.53,0.62,0.93,0.82);
        for (int m=0; m<3; ++m) leg->AddEntry(gNpe[m], labs[m], "pl");
        leg->SetTextSize(0.042); leg->Draw();

        TLine* l0 = new TLine(0,0,0,1800);
        l0->SetLineStyle(2); l0->SetLineColor(kGray+1); l0->Draw();

        c->Print("figs/fig_npe_x.pdf");
        delete c; delete mg;
    }

    // ── Figure 2: sigma(T0) vs x ──────────────────────────────────────────
    {
        TCanvas* c = new TCanvas("c2","",700,500);
        TMultiGraph* mg = new TMultiGraph();
        for (int m=0; m<3; ++m) mg->Add(gSig[m],"PL");
        mg->Draw("A");
        mg->GetXaxis()->SetTitle("Muon position  x  [mm]");
        mg->GetYaxis()->SetTitle("#sigma(T_{0})  [ps]");
        mg->GetXaxis()->SetLimits(-650, 650);
        mg->GetYaxis()->SetRangeUser(0, 80);   // start from 0, show scale
        mg->GetXaxis()->SetNdivisions(506);

        TLegend* leg = new TLegend(0.53,0.15,0.93,0.35);
        for (int m=0; m<3; ++m) leg->AddEntry(gSig[m], labs[m], "pl");
        leg->SetTextSize(0.042); leg->Draw();

        TLine* l0 = new TLine(0,0,0,80);
        l0->SetLineStyle(2); l0->SetLineColor(kGray+1); l0->Draw();

        c->Print("figs/fig_sigma_t_x.pdf");
        delete c; delete mg;
    }

    // ── Figure 3: R^N vs N  (analytic) ────────────────────────────────────
    {
        TCanvas* c = new TCanvas("c3","",600,480);
        const int NP = 101;
        double xN[NP], yR[NP];
        for (int i=0; i<NP; ++i) { xN[i]=i; yR[i]=std::pow(0.98, i); }
        TGraph* g = new TGraph(NP, xN, yR);
        g->SetLineColor(kBlue+1); g->SetLineWidth(2);
        g->Draw("AL");
        g->GetXaxis()->SetTitle("Number of reflections  N");
        g->GetYaxis()->SetTitle("Survival  R^{N}");
        g->GetXaxis()->SetLimits(0,100);
        g->GetYaxis()->SetRangeUser(0,1.05);
        g->GetXaxis()->SetNdivisions(505);
        g->GetYaxis()->SetNdivisions(506);

        // annotations
        TLatex lat; lat.SetTextSize(0.045); lat.SetTextFont(42);
        struct Ann { int N; const char* lbl; };
        Ann anns[] = {{10,"0.98^{10}=0.82"},{50,"0.98^{50}=0.36"},{100,"0.98^{100}=0.13"}};
        for (auto& a : anns) {
            double yv = std::pow(0.98, a.N);
            TLine* lv = new TLine(a.N, 0, a.N, yv);
            lv->SetLineStyle(3); lv->SetLineColor(kGray+1); lv->Draw();
            lat.DrawLatex(a.N+2, yv+0.03, a.lbl);
        }
        TLatex t2; t2.SetTextSize(0.042); t2.SetNDC();
        t2.DrawLatex(0.55,0.85,"R = 0.98 (Vikuiti 3M ESR)");

        c->Print("figs/fig_survival_RN.pdf");
        delete c;
    }

    // ── Figure 4: Bulk survival bar chart ────────────────────────────────
    {
        TCanvas* c = new TCanvas("c4","",550,480);
        const char* mlabs[3] = {"EJ-200","EJ-204","EJ-230"};
        double latts[3] = {380,160,120}; // cm
        double survs[3];
        for (int m=0; m<3; ++m) survs[m] = std::exp(-70.0/latts[m]);

        TH1F* h = new TH1F("hbulk","",3,0,3);
        for (int m=0; m<3; ++m) {
            h->SetBinContent(m+1, survs[m]);
            h->GetXaxis()->SetBinLabel(m+1, mlabs[m]);
        }
        h->SetFillColor(kBlue+1);
        h->GetYaxis()->SetTitle("Bulk survival  exp(-70 cm / #lambda_{att})");
        h->GetYaxis()->SetRangeUser(0, 1.1);
        h->GetYaxis()->SetNdivisions(506);
        h->SetBarWidth(0.6); h->SetBarOffset(0.2);
        h->Draw("B");

        // different colours per bar — need separate single-bin histos
        delete h;
        Color_t bcols[3] = {kBlue+1, kGreen+2, kRed+1};
        TH1F* hb[3];
        for (int m=0; m<3; ++m) {
            hb[m] = new TH1F(Form("hb%d",m),"",3,0,3);
            hb[m]->SetBinContent(m+1, survs[m]);
            hb[m]->GetXaxis()->SetBinLabel(m+1, mlabs[m]);
            hb[m]->SetFillColor(bcols[m]);
            hb[m]->SetBarWidth(0.6); hb[m]->SetBarOffset(0.2);
            hb[m]->GetYaxis()->SetTitle("Bulk survival  exp(-70 cm / #lambda_{att})");
            hb[m]->GetYaxis()->SetRangeUser(0, 1.1);
            hb[m]->GetYaxis()->SetNdivisions(506);
            if (m==0) hb[m]->Draw("B");
            else      hb[m]->Draw("B SAME");
        }
        TLatex lat; lat.SetTextSize(0.055); lat.SetTextFont(42); lat.SetTextAlign(21);
        for (int m=0; m<3; ++m)
            lat.DrawLatex(m+0.5, survs[m]+0.025, Form("%.2f", survs[m]));

        TLatex tit; tit.SetNDC(); tit.SetTextSize(0.042); tit.SetTextFont(42);
        tit.DrawLatex(0.15, 0.92, "Axial best case:  centre #rightarrow one end");

        c->Print("figs/fig_bulk_survival.pdf");
        delete c;
    }

    printf("\nAll figures written to figs/\n");
}
