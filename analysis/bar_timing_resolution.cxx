// bar_timing_resolution.cxx
// Timing resolution sigma(k) for the EJ-204 vs EJ-230 bar (TIR-only)
// Method identical to cylinder analysis (analysis/timing_resolution.cxx in ej204 repo):
//   sort photon times on one end face, fit Gaussian to k-th photon distribution.
//
// Run: root -l -q 'bar_timing_resolution.cxx'
//      (from build_t0minidaq/, where the ROOT files live)
//
// Input: photon_hits_ej204_tir.root  photon_hits_ej230_tir.root
// Face filter: 0 = left end (face_type=0, IDs 0-7), 8 SiPMs

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLine.h"
#include "TPad.h"
#include "TArrow.h"

#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <cmath>
#include <cstdio>
#include <fstream>

static const std::vector<int> kShowK = {1, 2, 5, 10, 20, 50};
static int kMaxK = 50;

struct EventData {
    int totalCount = 0;
    std::vector<double> sortedTimes;  // first kMaxK times (ns), sorted
};

// face_type filter: 0=left-end, 1=right-end, -1=both ends, 2=top, -2=all
std::map<int, EventData>
LoadTimes(TTree* tree, int faceFilter) {
    int    event_id, face_type;
    double time_ns;
    tree->SetBranchAddress("event_id",  &event_id);
    tree->SetBranchAddress("face_type", &face_type);
    tree->SetBranchAddress("time_ns",   &time_ns);

    std::map<int, std::vector<double>> raw;
    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        bool accept = (faceFilter == -2) ||
                      (faceFilter == -1 && (face_type == 0 || face_type == 1)) ||
                      (face_type == faceFilter);
        if (!accept) continue;
        raw[event_id].push_back(time_ns);
    }

    std::map<int, EventData> ev;
    for (auto& [id, v] : raw) {
        EventData& ed = ev[id];
        ed.totalCount = (int)v.size();
        if ((int)v.size() > kMaxK)
            std::partial_sort(v.begin(), v.begin() + kMaxK, v.end());
        else
            std::sort(v.begin(), v.end());
        int take = std::min((int)v.size(), kMaxK);
        ed.sortedTimes.assign(v.begin(), v.begin() + take);
    }
    return ev;
}

std::vector<double>
GetKthTimes_ps(const std::map<int, EventData>& ev, int k) {
    std::vector<double> vals;
    vals.reserve(ev.size());
    for (auto& [id, e] : ev)
        if ((int)e.sortedTimes.size() >= k)
            vals.push_back(e.sortedTimes[k-1] * 1e3);  // ns -> ps

    if (vals.empty()) return vals;
    auto sv = vals;
    std::sort(sv.begin(), sv.end());
    double med = sv[sv.size()/2];
    for (double& t : vals) t -= med;
    return vals;
}

struct FitResult { double sigGauss, sigErr, sigIQR; };

FitResult FitGaussIQR(const std::vector<double>& vals,
                       TH1D*& hOut,
                       const std::string& name, int color) {
    if (vals.size() < 20) return {0,0,0};

    auto sv = vals;
    std::sort(sv.begin(), sv.end());
    int n = (int)sv.size();
    auto at = [&](double q){ return sv[(int)(q*(n-1))]; };

    double q25 = at(0.25), q75 = at(0.75);
    double iqr   = q75 - q25;
    double sigIQR = iqr / 1.349;
    double q16 = at(0.16), q84 = at(0.84);
    double sigApprox = (q84 - q16) / 2.0;

    double hLo = -5.0*sigIQR, hHi = 5.0*sigIQR;
    int nbins = std::max(40, std::min(150, n/20));
    hOut = new TH1D(name.c_str(), ";#Deltat_{k} [ps];Events / bin",
                    nbins, hLo, hHi);
    hOut->SetLineColor(color); hOut->SetLineWidth(2);
    hOut->SetFillColorAlpha(color, 0.28);
    for (double v : vals) hOut->Fill(v);

    double fitLo = -2.0*sigApprox, fitHi = 2.0*sigApprox;
    TF1* f = new TF1((name+"_gaus").c_str(), "gaus", fitLo, fitHi);
    f->SetParameter(0, hOut->GetMaximum());
    f->SetParameter(1, 0.0);
    f->SetParameter(2, sigApprox);
    int fcolor = (color == kRed+1) ? kRed+3 : kBlue+3;
    f->SetLineColor(fcolor); f->SetLineWidth(2);
    hOut->Fit(f, "RQL0");

    return {std::abs(f->GetParameter(2)), f->GetParError(2), sigIQR};
}

// ─────────────────────────────────────────────────────────────────────────────
void bar_timing_resolution() {

    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetHistLineWidth(2);

    const int col4  = kBlue+1;   // EJ-204: blue
    const int col30 = kRed+1;    // EJ-230: red

    // ── Open files ────────────────────────────────────────────────────────────
    TFile* f4  = TFile::Open("photon_hits_ej204_tir.root",  "READ");
    TFile* f30 = TFile::Open("photon_hits_ej230_tir.root",  "READ");
    if (!f4 || !f30) {
        printf("ERROR: cannot open ROOT files.\n"
               "  Expected: photon_hits_ej204_tir.root  photon_hits_ej230_tir.root\n");
        return;
    }

    TTree* t4  = (TTree*)f4->Get("sipm_hits");
    TTree* t30 = (TTree*)f30->Get("sipm_hits");
    if (!t4 || !t30) { printf("ERROR: 'sipm_hits' tree not found\n"); return; }

    // face_type: 0=left-end, 1=right-end, -1=both-ends, 2=top-SiPMs
    // With TIR-only, end SiPMs get only ~0.4 photons/event (too few for k analysis).
    // Top SiPMs (face_type=2) get ~100-120 photons/event — used for sigma(k) study.
    const int faceFilter = 2;
    const char* faceName = "top SiPMs (face_type=2, 70 SiPMs)";

    printf("\n[bar_timing_resolution]  face: %s\n", faceName);
    printf("  Loading EJ-204 TIR-only times ...\n");
    auto ev4  = LoadTimes(t4,  faceFilter);
    printf("  Loading EJ-230 TIR-only times ...\n");
    auto ev30 = LoadTimes(t30, faceFilter);
    printf("  Events: EJ-204=%zu  EJ-230=%zu\n", ev4.size(), ev30.size());

    // True mean Npe per event at face
    long long n4=0, n30=0;
    for (auto& [id,e] : ev4 ) n4  += e.totalCount;
    for (auto& [id,e] : ev30) n30 += e.totalCount;
    double mean4  = ev4.empty()  ? 0.0 : (double)n4  / ev4.size();
    double mean30 = ev30.empty() ? 0.0 : (double)n30 / ev30.size();
    printf("  Mean Npe/event (%s): EJ-204=%.1f  EJ-230=%.1f  ratio=%.2f\n",
           faceName, mean4, mean30,
           mean30 > 0 ? mean4/mean30 : 0.0);

    // Adjust kMaxK to the actual data range
    int maxK4  = 0, maxK30 = 0;
    for (auto& [id,e] : ev4 ) maxK4  = std::max(maxK4,  (int)e.sortedTimes.size());
    for (auto& [id,e] : ev30) maxK30 = std::max(maxK30, (int)e.sortedTimes.size());
    int kEffMax = std::min(kMaxK, std::min(maxK4, maxK30));
    if (kEffMax < 2) { printf("Too few photons (kEffMax=%d) — check face filter.\n", kEffMax); return; }
    printf("  Using kMaxK = %d  (data has up to %d / %d photons per event)\n",
           kEffMax, maxK4, maxK30);

    // ── σ(k) for k=1..kEffMax ────────────────────────────────────────────────
    std::vector<double> kVals(kEffMax);
    std::vector<double> sigG4(kEffMax), sigGe4(kEffMax), sigI4(kEffMax);
    std::vector<double> sigG30(kEffMax), sigGe30(kEffMax), sigI30(kEffMax);
    std::iota(kVals.begin(), kVals.end(), 1.0);

    std::map<int,TH1D*> hShow4, hShow30;
    const std::vector<int> showKList = {1, 2, 5, 10, 20, std::min(50, kEffMax)};

    printf("  Computing sigma_gauss(k) for k=1..%d ...\n", kEffMax);
    for (int k = 1; k <= kEffMax; ++k) {
        auto v4  = GetKthTimes_ps(ev4,  k);
        auto v30 = GetKthTimes_ps(ev30, k);

        TH1D *h4=nullptr, *h30=nullptr;
        auto r4  = FitGaussIQR(v4,  h4,  "h4_k"  + std::to_string(k), col4);
        auto r30 = FitGaussIQR(v30, h30, "h30_k" + std::to_string(k), col30);

        sigG4[k-1]=r4.sigGauss;  sigGe4[k-1]=r4.sigErr;  sigI4[k-1]=r4.sigIQR;
        sigG30[k-1]=r30.sigGauss; sigGe30[k-1]=r30.sigErr; sigI30[k-1]=r30.sigIQR;

        bool show = std::find(showKList.begin(), showKList.end(), k) != showKList.end();
        if (show) { hShow4[k]=h4; hShow30[k]=h30; }
        else       { delete h4;   delete h30; }
    }

    printf("\n  k   sigma_G4    sigma_IQR4   sigma_G30   sigma_IQR30  ratio(4/30)\n");
    for (int k : showKList)
        if (k <= kEffMax)
            printf("  %-3d  %7.2f ps  %7.2f ps   %7.2f ps  %7.2f ps    %.3f\n",
                   k, sigG4[k-1], sigI4[k-1], sigG30[k-1], sigI30[k-1],
                   sigG30[k-1] > 0 ? sigG4[k-1]/sigG30[k-1] : 0.0);
    printf("\n");

    // ── Find k_opt via parabola fit ───────────────────────────────────────────
    int kSearchLo = 2, kSearchHi = std::min(30, kEffMax);
    auto findMin = [&](const std::vector<double>& sig) {
        int kBest = kSearchLo;
        double sMin = sig[kSearchLo-1];
        for (int i = kSearchLo; i <= kSearchHi; ++i)
            if (sig[i-1] < sMin && sig[i-1] > 0) { sMin = sig[i-1]; kBest = i; }
        return std::make_pair(kBest, sMin);
    };

    auto refineMin = [&](int kOpt, const std::vector<double>& sig) -> std::pair<double,double> {
        int lo2 = std::max(kSearchLo, kOpt-4);
        int hi2 = std::min(kSearchHi, kOpt+4);
        TGraph* gp = new TGraph();
        int npts = 0;
        for (int i = lo2; i <= hi2; ++i)
            if (sig[i-1] > 0) { gp->SetPoint(npts++, (double)i, sig[i-1]); }
        if (npts < 3) { delete gp; return {(double)kOpt, sig[kOpt-1]}; }
        TF1* fp = new TF1("fp_tmp","pol2", lo2, hi2);
        gp->Fit(fp,"QR0");
        double a = fp->GetParameter(2), b = fp->GetParameter(1);
        double kFine = (a > 0) ? -b/(2.*a) : (double)kOpt;
        double sFine = fp->Eval(kFine);
        delete gp;
        return {kFine, sFine};
    };

    auto [kOpt4_i,  sMin4]  = findMin(sigG4);
    auto [kOpt30_i, sMin30] = findMin(sigG30);
    auto [kOpt4_f,  sOpt4]  = refineMin(kOpt4_i,  sigG4);
    auto [kOpt30_f, sOpt30] = refineMin(kOpt30_i, sigG30);

    printf("  Optimal k (parabola fit, search k=%d..%d):\n", kSearchLo, kSearchHi);
    printf("    EJ-204: k_opt = %.1f  ->  sigma_opt = %.2f ps\n", kOpt4_f, sOpt4);
    printf("    EJ-230: k_opt = %.1f  ->  sigma_opt = %.2f ps\n", kOpt30_f, sOpt30);
    printf("\n");

    // TSV output
    {
        std::ofstream tsv("bar_timing_resolution.tsv");
        tsv << "k\tsigG_EJ204\tsigIQR_EJ204\tsigG_EJ230\tsigIQR_EJ230\n";
        for (int i = 0; i < kEffMax; ++i)
            tsv << (i+1) << "\t" << sigG4[i] << "\t" << sigI4[i]
                << "\t" << sigG30[i] << "\t" << sigI30[i] << "\n";
    }

    // ── Canvas 1: sigma(k) vs k ───────────────────────────────────────────────
    {
        TGraphErrors* g4  = new TGraphErrors(kEffMax);
        TGraphErrors* g30 = new TGraphErrors(kEffMax);
        TGraph* g4I  = new TGraph(kEffMax);
        TGraph* g30I = new TGraph(kEffMax);
        for (int i = 0; i < kEffMax; ++i) {
            g4->SetPoint(i,  kVals[i], sigG4[i]);  g4->SetPointError(i, 0, sigGe4[i]);
            g30->SetPoint(i, kVals[i], sigG30[i]); g30->SetPointError(i, 0, sigGe30[i]);
            g4I->SetPoint(i,  kVals[i], sigI4[i]);
            g30I->SetPoint(i, kVals[i], sigI30[i]);
        }
        g4->SetMarkerColor(col4);   g4->SetLineColor(col4);  g4->SetMarkerStyle(20); g4->SetMarkerSize(0.7);
        g30->SetMarkerColor(col30); g30->SetLineColor(col30); g30->SetMarkerStyle(21); g30->SetMarkerSize(0.7);
        g4I->SetLineColor(col4);    g4I->SetLineStyle(2);    g4I->SetLineWidth(1);
        g30I->SetLineColor(col30);  g30I->SetLineStyle(2);   g30I->SetLineWidth(1);

        TCanvas* c1 = new TCanvas("c_sigma","sigma_vs_k",900,600);
        c1->SetLeftMargin(0.12); c1->SetBottomMargin(0.12);

        // y range
        std::vector<double> yAll;
        for (int i=0;i<kEffMax;i++){
            if(sigG4[i]>0)  yAll.push_back(sigG4[i]);
            if(sigG30[i]>0) yAll.push_back(sigG30[i]);
            if(sigI4[i]>0)  yAll.push_back(sigI4[i]);
            if(sigI30[i]>0) yAll.push_back(sigI30[i]);
        }
        double ymax = yAll.empty() ? 100 : *std::max_element(yAll.begin(),yAll.end())*1.30;

        g4->SetTitle(Form(";k-th photon rank (1 = first);#sigma(t_{k}) [ps] (%s)", faceName));
        g4->GetYaxis()->SetRangeUser(0, ymax);
        g4->GetXaxis()->SetTitleSize(0.047); g4->GetYaxis()->SetTitleSize(0.047);
        g4->Draw("APE");
        g30->Draw("PE SAME");
        g4I->Draw("L SAME");
        g30I->Draw("L SAME");

        TLegend* leg = new TLegend(0.58,0.60,0.93,0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(g4,   "EJ-204 (TIR-only): #sigma_{Gauss}","pe");
        leg->AddEntry(g4I,  "EJ-204 (TIR-only): #sigma_{IQR}","l");
        leg->AddEntry(g30,  "EJ-230 (TIR-only): #sigma_{Gauss}","pe");
        leg->AddEntry(g30I, "EJ-230 (TIR-only): #sigma_{IQR}","l");
        leg->Draw();

        // k_opt markers
        auto drawOpt = [&](double kf, double sf, int col, const char* lbl, double offY) {
            TLine* ll = new TLine(kf, 0, kf, sf+0.5);
            ll->SetLineColor(col); ll->SetLineStyle(2); ll->SetLineWidth(2); ll->Draw();
            TArrow* ar = new TArrow(kf, sf+1.8+offY, kf, sf+0.2, 0.015, "|>");
            ar->SetFillColor(col); ar->SetLineColor(col); ar->SetLineWidth(2); ar->Draw();
            TLatex lx; lx.SetTextSize(0.031); lx.SetTextColor(col); lx.SetTextAlign(22);
            lx.DrawLatex(kf, sf+2.4+offY, Form("%s k_{opt}=%.0f, #sigma=%.1fps", lbl, kf, sf));
        };
        double offV = (std::abs(kOpt4_f - kOpt30_f) < 3) ? 3.0 : 0.0;
        drawOpt(kOpt4_f,  sOpt4,  col4,  "EJ-204", offV);
        drawOpt(kOpt30_f, sOpt30, col30, "EJ-230", 0.0);

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.031);
        lat.DrawLatex(0.13,0.93,
            "EJ-204 vs EJ-230 bar (TIR-only) | 1.4mx60mmx10mm | 1 GeV #mu^{-} at x=0 | 5000 ev");

        c1->SaveAs("bar_timing_sigma_vs_k.pdf");
        c1->SaveAs("bar_timing_sigma_vs_k.png");
        printf("  Saved: bar_timing_sigma_vs_k.pdf/png\n");
    }

    // ── Canvas 2: 6-panel distributions ──────────────────────────────────────
    {
        // Only show k values that we actually have
        std::vector<int> panelK;
        for (int k : showKList)
            if (k <= kEffMax && hShow4.count(k) && hShow30.count(k))
                panelK.push_back(k);

        int ncols = 3, nrows = std::max(1, (int)(panelK.size()+2)/3);
        TCanvas* c2 = new TCanvas("c_fits","Timing distributions",
                                   1200, 400*nrows);
        c2->Divide(ncols, nrows, 0.004, 0.006);

        for (int ip = 0; ip < (int)panelK.size(); ++ip) {
            int k = panelK[ip];
            c2->cd(ip+1);
            gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.10);

            TH1D* h4  = hShow4[k];
            TH1D* h30 = hShow30[k];
            if (!h4 || !h30) continue;

            double ymax2 = std::max(h4->GetMaximum(), h30->GetMaximum()) * 1.35;
            h4->SetMaximum(ymax2);
            h4->SetTitle(Form(";#Deltat_{k=%d} [ps];Events/bin", k));
            h4->GetXaxis()->SetTitleSize(0.065); h4->GetXaxis()->SetLabelSize(0.055);
            h4->GetYaxis()->SetTitleSize(0.065); h4->GetYaxis()->SetLabelSize(0.055);
            h4->GetXaxis()->SetTitleOffset(1.05); h4->GetYaxis()->SetTitleOffset(1.10);
            h4->Draw("HIST");
            h30->Draw("HIST SAME");

            // Redraw fits
            TF1* f4r  = h4->GetFunction((std::string("h4_k") +std::to_string(k)+"_gaus").c_str());
            TF1* f30r = h30->GetFunction((std::string("h30_k")+std::to_string(k)+"_gaus").c_str());
            if (f4r)  { f4r->SetRange(f4r->GetXmin(), f4r->GetXmax()); f4r->Draw("SAME"); }
            if (f30r) { f30r->SetRange(f30r->GetXmin(),f30r->GetXmax()); f30r->Draw("SAME"); }

            TLegend* lp = new TLegend(0.35,0.68,0.97,0.90);
            lp->SetBorderSize(0); lp->SetFillStyle(0); lp->SetTextSize(0.052);
            lp->AddEntry(h4,  Form("EJ-204 #sigma=%.1fps",sigG4[k-1]), "lf");
            lp->AddEntry(h30, Form("EJ-230 #sigma=%.1fps",sigG30[k-1]),"lf");
            lp->Draw();

            TLatex ptxt; ptxt.SetNDC(); ptxt.SetTextSize(0.058); ptxt.SetTextAlign(22);
            ptxt.DrawLatex(0.5, 0.97, Form("k = %d", k));
        }

        c2->SaveAs("bar_timing_fits_panels.pdf");
        c2->SaveAs("bar_timing_fits_panels.png");
        printf("  Saved: bar_timing_fits_panels.pdf/png\n");
        printf("\nDone.\n");
    }
}
