// compare_vikuiti_tir.cxx
// Compares EJ-228 cylinder light yield: Vikuiti R=0.98 vs TIR-only
// Run with: root -l -q compare_vikuiti_tir.cxx
//
// Expects (relative to working dir):
//   photon_hits_vikuiti.root
//   photon_hits_tir_only.root

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"

#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <fstream>

// ── Per-event accumulator ────────────────────────────────────────────────────

struct PerEvt {
    int total = 0;
    int top   = 0;  // face_type == 0
    int bot   = 0;  // face_type == 1
    int sipm[8] = {};
};

std::map<int,PerEvt> BuildEventMap(TTree* tree) {
    int event_id, face_type, global_id;
    tree->SetBranchAddress("event_id",  &event_id);
    tree->SetBranchAddress("face_type", &face_type);
    tree->SetBranchAddress("global_id", &global_id);

    std::map<int,PerEvt> ev;
    Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        auto& e = ev[event_id];
        e.total++;
        if (face_type == 0) e.top++;
        else                 e.bot++;
        if (global_id >= 0 && global_id < 8) e.sipm[global_id]++;
    }
    return ev;
}

// ── Stats helper ─────────────────────────────────────────────────────────────

struct Stats { double mean, sigma, median, p25, p75; };

Stats CalcStats(const std::vector<int>& v) {
    double s = 0, s2 = 0;
    for (int x : v) { s += x; s2 += (double)x*x; }
    double n = v.size();
    double mean  = s / n;
    double sigma = std::sqrt(s2/n - mean*mean);
    // median / quartiles via sorted copy
    auto sv = v;
    std::sort(sv.begin(), sv.end());
    auto at = [&](double q){ return sv[(int)(q*(n-1))]; };
    return {mean, sigma, (double)at(0.5), (double)at(0.25), (double)at(0.75)};
}

// ── Main ─────────────────────────────────────────────────────────────────────

void compare_vikuiti_tir() {

    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetHistLineWidth(2);

    // ── Open files ────────────────────────────────────────────────────────────
    TFile* fV = TFile::Open("photon_hits_vikuiti.root",  "READ");
    TFile* fT = TFile::Open("photon_hits_tir_only.root", "READ");
    if (!fV || !fT) { printf("ERROR: cannot open ROOT files\n"); return; }

    TTree* tV = (TTree*)fV->Get("sipm_hits");
    TTree* tT = (TTree*)fT->Get("sipm_hits");

    // ── Build per-event maps ──────────────────────────────────────────────────
    auto evV = BuildEventMap(tV);
    auto evT = BuildEventMap(tT);

    // Convert to vectors for stats
    std::vector<int> npe_V, npe_T, npeTop_V, npeTop_T, npeBot_V, npeBot_T;
    for (auto& [id,e] : evV) { npe_V.push_back(e.total); npeTop_V.push_back(e.top); npeBot_V.push_back(e.bot); }
    for (auto& [id,e] : evT) { npe_T.push_back(e.total); npeTop_T.push_back(e.top); npeBot_T.push_back(e.bot); }

    Stats sV    = CalcStats(npe_V);
    Stats sT    = CalcStats(npe_T);
    Stats sTopV = CalcStats(npeTop_V);
    Stats sTopT = CalcStats(npeTop_T);
    Stats sBotV = CalcStats(npeBot_V);
    Stats sBotT = CalcStats(npeBot_T);

    printf("\n");
    printf("=============================================================\n");
    printf("  EJ-228 Cylinder  —  Vikuiti vs TIR-only comparison\n");
    printf("  5000 events, 1 GeV mu-, horizontal through cylinder centre\n");
    printf("=============================================================\n");
    printf("                    Vikuiti R=0.98    TIR-only    Ratio\n");
    printf("  Total hits       :  %7.0f          %7.0f     %.3f×\n",
           (double)tV->GetEntries(), (double)tT->GetEntries(),
           (double)tV->GetEntries()/(double)tT->GetEntries());
    printf("  Npe/event (mean) :  %7.1f          %7.1f     %.3f×\n",
           sV.mean, sT.mean, sV.mean/sT.mean);
    printf("  Npe/event (σ)    :  %7.1f          %7.1f\n", sV.sigma, sT.sigma);
    printf("  Npe/event (med)  :  %7.1f          %7.1f\n", sV.median, sT.median);
    printf("  Top cap (mean)   :  %7.1f          %7.1f     %.3f×\n",
           sTopV.mean, sTopT.mean, sTopV.mean/sTopT.mean);
    printf("  Bot cap (mean)   :  %7.1f          %7.1f     %.3f×\n",
           sBotV.mean, sBotT.mean, sBotV.mean/sBotT.mean);
    printf("  Efficiency       :  %6.2f %%         %6.2f %%\n",
           100.*tV->GetEntries()/(5000.*47700.), 100.*tT->GetEntries()/(5000.*47700.));
    printf("=============================================================\n\n");

    // Write TSV for Beamer
    std::ofstream tsv("compare_summary.tsv");
    tsv << "config\tnpe_mean\tnpe_sigma\tnpe_median\tnpe_top_mean\tnpe_bot_mean\ttotal_hits\n";
    tsv << "Vikuiti\t" << sV.mean << "\t" << sV.sigma << "\t" << sV.median
        << "\t" << sTopV.mean << "\t" << sBotV.mean << "\t" << tV->GetEntries() << "\n";
    tsv << "TIR-only\t" << sT.mean << "\t" << sT.sigma << "\t" << sT.median
        << "\t" << sTopT.mean << "\t" << sBotT.mean << "\t" << tT->GetEntries() << "\n";
    tsv.close();

    // Per-SiPM means
    printf("  Per-SiPM mean Npe/event:\n");
    printf("  ID   Face    Vikuiti    TIR-only   Ratio\n");
    for (int g = 0; g < 8; ++g) {
        std::vector<int> vsipm, tsipm;
        for (auto& [id,e] : evV) vsipm.push_back(e.sipm[g]);
        for (auto& [id,e] : evT) tsipm.push_back(e.sipm[g]);
        Stats sv2 = CalcStats(vsipm);
        Stats st2 = CalcStats(tsipm);
        printf("  %-4d %-8s %7.1f    %7.1f    %.3f×\n",
               g, g<4?"Top":"Bot", sv2.mean, st2.mean, sv2.mean/st2.mean);
    }
    printf("\n");

    // ── Histograms ────────────────────────────────────────────────────────────

    int xmax = (int)(sV.mean * 1.6);
    xmax = ((xmax / 500) + 1) * 500;

    TH1D* hV    = new TH1D("hV",    ";N_{pe} per event;Counts/bin", 100, 0, xmax);
    TH1D* hT    = new TH1D("hT",    ";N_{pe} per event;Counts/bin", 100, 0, xmax);
    TH1D* hTopV = new TH1D("hTopV", ";N_{pe} per event (top cap);Counts/bin", 100, 0, xmax/2);
    TH1D* hTopT = new TH1D("hTopT", ";N_{pe} per event (top cap);Counts/bin", 100, 0, xmax/2);
    TH1D* hBotV = new TH1D("hBotV", ";N_{pe} per event (bottom cap);Counts/bin", 100, 0, xmax/2);
    TH1D* hBotT = new TH1D("hBotT", ";N_{pe} per event (bottom cap);Counts/bin", 100, 0, xmax/2);

    for (auto& [id,e] : evV) { hV->Fill(e.total); hTopV->Fill(e.top); hBotV->Fill(e.bot); }
    for (auto& [id,e] : evT) { hT->Fill(e.total); hTopT->Fill(e.top); hBotT->Fill(e.bot); }

    // Colors
    hV->SetLineColor(kBlue+1);    hV->SetFillColorAlpha(kBlue+1,  0.25);
    hT->SetLineColor(kOrange+7);  hT->SetFillColorAlpha(kOrange+7,0.25);
    hTopV->SetLineColor(kBlue+1); hTopV->SetFillColorAlpha(kBlue+1, 0.25);
    hTopT->SetLineColor(kOrange+7); hTopT->SetFillColorAlpha(kOrange+7,0.25);
    hBotV->SetLineColor(kBlue+1); hBotV->SetFillColorAlpha(kBlue+1, 0.25);
    hBotT->SetLineColor(kOrange+7); hBotT->SetFillColorAlpha(kOrange+7,0.25);

    // ── Canvas 1: Total Npe comparison ───────────────────────────────────────
    {
        TCanvas* c1 = new TCanvas("c1","Npe total comparison",900,600);
        c1->SetLeftMargin(0.12);

        double ymax = std::max(hV->GetMaximum(), hT->GetMaximum()) * 1.15;
        hV->SetMaximum(ymax);
        hT->SetMaximum(ymax);

        hV->Draw("HIST");
        hT->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(hV, Form("Vikuiti R=0.98   #mu=%.0f", sV.mean), "f");
        leg->AddEntry(hT, Form("TIR-only         #mu=%.0f", sT.mean), "f");
        leg->Draw();

        TLatex lat;
        lat.SetNDC(); lat.SetTextSize(0.038);
        lat.DrawLatex(0.15, 0.92, "EJ-228 cylinder  25#times25 mm^{2}  |  1 GeV #mu^{-}  |  5000 events");
        lat.SetTextSize(0.033);
        lat.DrawLatex(0.55, 0.67, Form("Ratio V/TIR = %.2f#times", sV.mean/sT.mean));

        TLine* lV = new TLine(sV.mean,0,sV.mean,ymax*0.6);
        TLine* lT = new TLine(sT.mean,0,sT.mean,ymax*0.6);
        lV->SetLineColor(kBlue+1);    lV->SetLineStyle(2); lV->SetLineWidth(2); lV->Draw();
        lT->SetLineColor(kOrange+7);  lT->SetLineStyle(2); lT->SetLineWidth(2); lT->Draw();

        c1->SaveAs("npe_total_comparison.pdf");
        c1->SaveAs("npe_total_comparison.png");
        printf("  Saved: npe_total_comparison.pdf/png\n");
    }

    // ── Canvas 2: Per-cap comparison ─────────────────────────────────────────
    {
        TCanvas* c2 = new TCanvas("c2","Per-cap comparison",1200,500);
        c2->Divide(2,1);

        c2->cd(1); gPad->SetLeftMargin(0.14);
        double ym = std::max(hTopV->GetMaximum(), hTopT->GetMaximum()) * 1.15;
        hTopV->SetMaximum(ym); hTopT->SetMaximum(ym);
        hTopV->SetTitle("Top cap (IDs 0-3)");
        hTopV->Draw("HIST"); hTopT->Draw("HIST SAME");
        TLegend* l1 = new TLegend(0.45,0.72,0.89,0.88);
        l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.038);
        l1->AddEntry(hTopV, Form("Vikuiti #mu=%.0f", sTopV.mean), "f");
        l1->AddEntry(hTopT, Form("TIR-only #mu=%.0f", sTopT.mean), "f");
        l1->Draw();

        c2->cd(2); gPad->SetLeftMargin(0.14);
        ym = std::max(hBotV->GetMaximum(), hBotT->GetMaximum()) * 1.15;
        hBotV->SetMaximum(ym); hBotT->SetMaximum(ym);
        hBotV->SetTitle("Bottom cap (IDs 4-7)");
        hBotV->Draw("HIST"); hBotT->Draw("HIST SAME");
        TLegend* l2 = new TLegend(0.45,0.72,0.89,0.88);
        l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.038);
        l2->AddEntry(hBotV, Form("Vikuiti #mu=%.0f", sBotV.mean), "f");
        l2->AddEntry(hBotT, Form("TIR-only #mu=%.0f", sBotT.mean), "f");
        l2->Draw();

        c2->SaveAs("npe_percap_comparison.pdf");
        c2->SaveAs("npe_percap_comparison.png");
        printf("  Saved: npe_percap_comparison.pdf/png\n");
    }

    // ── Canvas 3: Per-SiPM — side-by-side via TGraph filled boxes + log ─────
    // Uses two TGraph with marker style 22 (triangle) to avoid the
    // SetBarOffset+SetBinLabel axis-collision issue of the BAR draw option.
    // Each SiPM gets two markers separated by ±0.22 on the x-axis.
    {
        TCanvas* c3 = new TCanvas("c3","Per-SiPM",900,520);
        c3->SetLeftMargin(0.14);
        c3->SetBottomMargin(0.13);
        c3->SetLogy();

        // Collect per-SiPM means
        double meanV[8], meanT[8];
        for (int g = 0; g < 8; ++g) {
            std::vector<int> vv, tv;
            for (auto& [id,e] : evV) vv.push_back(e.sipm[g]);
            for (auto& [id,e] : evT) tv.push_back(e.sipm[g]);
            meanV[g] = CalcStats(vv).mean;
            meanT[g] = CalcStats(tv).mean;
        }

        // Draw frame with custom axis
        double yLo = 80, yHi = 2600;
        TH1D* hFrame = new TH1D("hFrame",
            "Mean N_{pe}/event per SiPM;Global SiPM ID;Mean N_{pe} / event",
            8, -0.5, 7.5);
        hFrame->SetMinimum(yLo);
        hFrame->SetMaximum(yHi);
        hFrame->GetXaxis()->SetNdivisions(8);
        hFrame->GetXaxis()->SetLabelSize(0.044);
        hFrame->GetYaxis()->SetLabelSize(0.044);
        hFrame->GetXaxis()->SetTitleSize(0.047);
        hFrame->GetYaxis()->SetTitleSize(0.047);
        hFrame->GetXaxis()->SetTitleOffset(1.0);
        hFrame->SetLineColor(0); // invisible
        hFrame->Draw("HIST");    // just the frame

        // Filled rectangles drawn as TBox pairs for each SiPM
        const double hw = 0.22;   // half-width of each bar in x-axis units
        const double gap = 0.00;  // gap between Vikuiti and TIR-only bars
        for (int g = 0; g < 8; ++g) {
            double cx = g;
            // Vikuiti bar (left)
            TBox* bV = new TBox(cx - hw - gap, yLo, cx - gap, meanV[g]);
            bV->SetFillColor(kBlue+1); bV->SetLineColor(kBlue+3); bV->SetLineWidth(1);
            bV->Draw("l");
            // TIR-only bar (right)
            TBox* bT = new TBox(cx + gap, yLo, cx + hw + gap, meanT[g]);
            bT->SetFillColor(kOrange+7); bT->SetLineColor(kOrange+4); bT->SetLineWidth(1);
            bT->Draw("l");
        }

        // Redraw axes on top so bars don't cover them
        gPad->RedrawAxis();

        // Divider line between top and bottom cap groups
        TLine* divL = new TLine(3.5, yLo, 3.5, yHi*0.35);
        divL->SetLineStyle(2); divL->SetLineColor(kGray+1); divL->SetLineWidth(2);
        divL->Draw();

        // x-axis integer labels
        TLatex xLab; xLab.SetTextSize(0.042); xLab.SetTextAlign(22);
        xLab.SetTextColor(kGray+2);
        for (int g = 0; g < 8; ++g)
            xLab.DrawLatex(g, yLo*0.55, std::to_string(g).c_str());

        // Legend
        TLegend* l3 = new TLegend(0.55, 0.76, 0.88, 0.90);
        l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.040);
        // dummy entries using colored squares
        TBox* dV = new TBox(0,0,1,1); dV->SetFillColor(kBlue+1);
        TBox* dT = new TBox(0,0,1,1); dT->SetFillColor(kOrange+7);
        l3->AddEntry(dV, Form("Vikuiti   #mu/SiPM=%.0f", sV.mean/8.), "f");
        l3->AddEntry(dT, Form("TIR-only  #mu/SiPM=%.0f", sT.mean/8.), "f");
        l3->Draw();

        TLatex lat3; lat3.SetNDC(); lat3.SetTextSize(0.032);
        lat3.DrawLatex(0.15, 0.93,
            "EJ-228 cylinder | IDs 0-3: top cap | IDs 4-7: bottom cap | log scale");
        lat3.SetTextColor(kGray+2); lat3.SetTextSize(0.036);
        lat3.DrawLatex(0.25, 0.16, "top cap");
        lat3.DrawLatex(0.63, 0.16, "bottom cap");

        c3->SaveAs("npe_persipm_comparison.pdf");
        c3->SaveAs("npe_persipm_comparison.png");
        printf("  Saved: npe_persipm_comparison.pdf/png  (TBox side-by-side, log)\n");
    }

    fV->Close(); fT->Close();
    printf("\nDone. Check: npe_total_comparison.pdf, npe_percap_comparison.pdf, npe_persipm_comparison.pdf\n\n");
}
