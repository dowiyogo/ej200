// timing_resolution.cxx
// Timing resolution σ(k) as a function of k-th photon rank
// for EJ-228 cylinder: Vikuiti R=0.98 vs TIR-only
//
// Method:
//   For each event, sort the photon arrival times on the top cap.
//   Collect t_k (k-th earliest) across all events → distribution.
//   Fit a Gaussian to the CORE using IQR-based fit range (robust against
//   the right-skew intrinsic to order-statistic distributions).
//   Plot σ_gauss(k) vs k for both configurations, plus 6 example panels.
//
// Run: root -l -q 'timing_resolution.cxx("top")'
//   or root -l -q 'timing_resolution.cxx("all")'  (combine both caps)
//
// Input (same dir): photon_hits_vikuiti.root  photon_hits_tir_only.root

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphErrors.h"
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
static const int kMaxK = 50;   // maximum k to compute
static int gFaceFilter = 0;    // 0=top, 1=bot, -1=all

// ── Per-event data ────────────────────────────────────────────────────────────
struct EventData {
    int totalCount = 0;                // total hits BEFORE truncation
    std::vector<double> sortedTimes;   // first kMaxK times (ns), sorted
};

// ── Load and sort photon times ────────────────────────────────────────────────
std::map<int, EventData>
LoadTimes(TTree* tree, int faceFilter) {
    int    event_id, face_type;
    double time_ns;
    tree->SetBranchAddress("event_id",  &event_id);
    tree->SetBranchAddress("face_type", &face_type);
    tree->SetBranchAddress("time_ns",   &time_ns);

    // Accumulate all times (will be large but manageable for 5000 events × ~4800 photons = 24M)
    std::map<int, std::vector<double>> raw;
    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (faceFilter >= 0 && face_type != faceFilter) continue;
        raw[event_id].push_back(time_ns);
    }

    std::map<int, EventData> ev;
    for (auto& [id, v] : raw) {
        EventData& ed = ev[id];
        ed.totalCount = (int)v.size();
        // Partial sort: only need the first kMaxK
        if ((int)v.size() > kMaxK)
            std::partial_sort(v.begin(), v.begin() + kMaxK, v.end());
        else
            std::sort(v.begin(), v.end());
        int take = std::min((int)v.size(), kMaxK);
        ed.sortedTimes.assign(v.begin(), v.begin() + take);
    }
    return ev;
}

// ── Extract k-th times (ps) relative to median across events ─────────────────
// Using median (not mean) as reference: more robust for skewed distributions.
std::vector<double>
GetKthTimes_ps(const std::map<int, EventData>& ev, int k) {
    std::vector<double> vals;
    vals.reserve(ev.size());
    for (auto& [id, e] : ev)
        if ((int)e.sortedTimes.size() >= k)
            vals.push_back(e.sortedTimes[k-1] * 1e3);  // ns → ps

    // Subtract median (robust center for skewed distributions)
    auto sv = vals;
    std::sort(sv.begin(), sv.end());
    double med = sv[sv.size()/2];
    for (double& t : vals) t -= med;
    return vals;
}

// ── Robust Gaussian fit using IQR-derived fit range ──────────────────────────
// Returns {sigma_gauss_ps, sigma_err_ps, sigma_iqr_ps}
struct FitResult { double sigGauss, sigErr, sigIQR; };

FitResult FitGaussIQR(const std::vector<double>& vals,
                       TH1D*& hOut,
                       const std::string& name, int color) {
    if (vals.empty()) return {0,0,0};

    // Empirical quantiles
    auto sv = vals;
    std::sort(sv.begin(), sv.end());
    int n = sv.size();
    auto at = [&](double q) { return sv[(int)(q*(n-1))]; };

    double q25 = at(0.25), q75 = at(0.75);
    double iqr  = q75 - q25;
    double sigIQR = iqr / 1.349;          // IQR → σ for Gaussian
    double q16  = at(0.16), q84 = at(0.84);
    double sigApprox = (q84 - q16) / 2.0; // ≈ σ from 68% quantile range

    // Histogram: ±5 IQR-sigmas around median (already subtracted → median~0)
    double hLo = -5.0 * sigIQR, hHi = 5.0 * sigIQR;
    int nbins = std::max(60, std::min(200, (int)(n/25)));
    hOut = new TH1D(name.c_str(),
                    ";#Deltat_{k} [ps];Events / bin",
                    nbins, hLo, hHi);
    hOut->SetLineColor(color);
    hOut->SetLineWidth(2);
    hOut->SetFillColorAlpha(color, 0.30);
    for (double v : vals) hOut->Fill(v);

    // Fit range: ±2.0 × sigApprox around peak (core Gaussian, ignore tails)
    double fitLo = -2.0 * sigApprox;
    double fitHi =  2.0 * sigApprox;
    TF1* f = new TF1((name + "_gaus").c_str(), "gaus", fitLo, fitHi);
    f->SetParameter(0, hOut->GetMaximum());
    f->SetParameter(1, 0.0);
    f->SetParameter(2, sigApprox);
    f->SetLineColor(color == kOrange+7 ? kOrange+4 : kBlue+3);
    f->SetLineWidth(2);
    hOut->Fit(f, "RQL0");  // R=range, Q=quiet, L=logL, 0=don't draw

    double sigGauss = std::abs(f->GetParameter(2));
    double sigErr   = f->GetParError(2);
    return {sigGauss, sigErr, sigIQR};
}

// ─────────────────────────────────────────────────────────────────────────────
void timing_resolution(const char* cap = "top") {

    gFaceFilter = std::string(cap) == "all" ? -1 :
                  std::string(cap) == "bot" ?  1 : 0;
    const char* capName = gFaceFilter<0 ? "top+bot caps" :
                          gFaceFilter==0 ? "top cap (IDs 0-3)" : "bot cap (IDs 4-7)";

    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetHistLineWidth(2);

    const int colV = kBlue+1;
    const int colT = kOrange+7;

    // ── Open files ────────────────────────────────────────────────────────────
    TFile* fV = TFile::Open("photon_hits_vikuiti.root",  "READ");
    TFile* fT = TFile::Open("photon_hits_tir_only.root", "READ");
    if (!fV || !fT) { printf("ERROR: cannot open ROOT files\n"); return; }

    TTree* tV = (TTree*)fV->Get("sipm_hits");
    TTree* tT = (TTree*)fT->Get("sipm_hits");

    printf("\n[timing_resolution]  face: %s\n", capName);
    printf("  Loading Vikuiti times ...\n");
    auto evV = LoadTimes(tV, gFaceFilter);
    printf("  Loading TIR-only times ...\n");
    auto evT = LoadTimes(tT, gFaceFilter);
    printf("  Events loaded: V=%zu  T=%zu\n", evV.size(), evT.size());

    // True mean Npe per event
    {
        long long nV = 0, nT = 0;
        for (auto& [id,e] : evV) nV += e.totalCount;
        for (auto& [id,e] : evT) nT += e.totalCount;
        printf("  Mean Npe/event (%s): Vikuiti=%.1f  TIR-only=%.1f  ratio=%.2f×\n",
               capName,
               (double)nV/evV.size(),
               (double)nT/evT.size(),
               (double)nV/(double)nT);
    }

    // ── Compute σ(k) for k=1..kMaxK ─────────────────────────────────────────
    std::vector<double> kVals(kMaxK),
        sigGV(kMaxK), sigGeV(kMaxK), sigIV(kMaxK),
        sigGT(kMaxK), sigGeT(kMaxK), sigIT(kMaxK);
    std::iota(kVals.begin(), kVals.end(), 1.0);

    std::map<int, TH1D*> hShowV, hShowT;

    printf("  Computing σ_gauss(k) and σ_IQR(k) for k=1..%d ...\n", kMaxK);
    for (int k = 1; k <= kMaxK; ++k) {
        auto vV = GetKthTimes_ps(evV, k);
        auto vT = GetKthTimes_ps(evT, k);

        TH1D* hV = nullptr, *hT = nullptr;
        auto rV = FitGaussIQR(vV, hV, "hV_k" + std::to_string(k), colV);
        auto rT = FitGaussIQR(vT, hT, "hT_k" + std::to_string(k), colT);

        sigGV[k-1]=rV.sigGauss; sigGeV[k-1]=rV.sigErr; sigIV[k-1]=rV.sigIQR;
        sigGT[k-1]=rT.sigGauss; sigGeT[k-1]=rT.sigErr; sigIT[k-1]=rT.sigIQR;

        bool show = std::find(kShowK.begin(),kShowK.end(),k) != kShowK.end();
        if (show) { hShowV[k]=hV; hShowT[k]=hT; }
        else       { delete hV; delete hT; }
    }

    printf("\n  k   σ_gauss_V   σ_IQR_V    σ_gauss_T   σ_IQR_T    ratio_gauss\n");
    for (int k : kShowK)
        printf("  %-3d  %7.2f ps  %7.2f ps   %7.2f ps  %7.2f ps   %.3f\n",
               k, sigGV[k-1], sigIV[k-1], sigGT[k-1], sigIT[k-1],
               sigGV[k-1]/sigGT[k-1]);
    printf("\n");

    // ── Find optimal k (minimum of σ_gauss) within k=2..30 ───────────────────
    // Search range avoids k=1 (dominated by Gaussian fit artifact for minimum)
    // and k>30 (complex multi-bounce regime).
    int kSearchLo = 2, kSearchHi = std::min(30, kMaxK);

    auto findMin = [&](const std::vector<double>& sig) {
        int kBest = kSearchLo;
        double sMin = sig[kSearchLo-1];
        for (int i = kSearchLo; i <= kSearchHi; ++i)
            if (sig[i-1] < sMin) { sMin = sig[i-1]; kBest = i; }
        return std::make_pair(kBest, sMin);
    };

    auto [kOptV, sigMinV] = findMin(sigGV);
    auto [kOptT, sigMinT] = findMin(sigGT);

    // Parabola refinement near minimum (±4 points, pol2 fit)
    auto refineMin = [&](int kOpt, const std::vector<double>& sig)
            -> std::pair<double,double> {
        int lo2 = std::max(kSearchLo, kOpt-4);
        int hi2 = std::min(kSearchHi, kOpt+4);
        TGraph* gp = new TGraph();
        for (int i = lo2; i <= hi2; ++i)
            gp->SetPoint(i-lo2, (double)i, sig[i-1]);
        TF1* fp = new TF1("fp_tmp", "pol2", lo2, hi2);
        gp->Fit(fp, "QR0");
        double a = fp->GetParameter(2);
        double b = fp->GetParameter(1);
        double kFine = (a > 0) ? -b/(2.*a) : (double)kOpt;
        double sFine = fp->Eval(kFine);
        delete gp;
        return {kFine, sFine};
    };

    auto [kOptV_f, sigOptV_f] = refineMin(kOptV, sigGV);
    auto [kOptT_f, sigOptT_f] = refineMin(kOptT, sigGT);

    printf("  Optimal k (minimum σ_gauss, parabola fit):\n");
    printf("    Vikuiti:   k_opt = %.1f  →  σ_opt = %.2f ps\n", kOptV_f, sigOptV_f);
    printf("    TIR-only:  k_opt = %.1f  →  σ_opt = %.2f ps\n", kOptT_f, sigOptT_f);
    printf("\n");

    // Write TSV for Beamer
    {
        std::ofstream tsv("timing_resolution.tsv");
        tsv << "k\tsigGauss_V\tsigIQR_V\tsigGauss_T\tsigIQR_T\n";
        for (int i = 0; i < kMaxK; ++i)
            tsv << (i+1) << "\t" << sigGV[i] << "\t" << sigIV[i]
                << "\t" << sigGT[i] << "\t" << sigIT[i] << "\n";
    }

    // ── Canvas 1: σ(k) vs k — two curves, both Gauss and IQR ────────────────
    {
        // Gauss graphs
        TGraphErrors* gV = new TGraphErrors(kMaxK);
        TGraphErrors* gT = new TGraphErrors(kMaxK);
        // IQR graphs (dashed, no error bars)
        TGraph* gVI = new TGraph(kMaxK);
        TGraph* gTI = new TGraph(kMaxK);
        for (int i = 0; i < kMaxK; ++i) {
            gV->SetPoint(i, kVals[i], sigGV[i]); gV->SetPointError(i,0,sigGeV[i]);
            gT->SetPoint(i, kVals[i], sigGT[i]); gT->SetPointError(i,0,sigGeT[i]);
            gVI->SetPoint(i, kVals[i], sigIV[i]);
            gTI->SetPoint(i, kVals[i], sigIT[i]);
        }
        gV->SetMarkerColor(colV);  gV->SetLineColor(colV);
        gV->SetMarkerStyle(20);    gV->SetMarkerSize(0.7);
        gT->SetMarkerColor(colT);  gT->SetLineColor(colT);
        gT->SetMarkerStyle(21);    gT->SetMarkerSize(0.7);
        gVI->SetLineColor(colV);   gVI->SetLineStyle(2);  gVI->SetLineWidth(1);
        gTI->SetLineColor(colT);   gTI->SetLineStyle(2);  gTI->SetLineWidth(1);

        TCanvas* c1 = new TCanvas("c_sigma","Timing resolution",900,600);
        c1->SetLeftMargin(0.12); c1->SetBottomMargin(0.12);

        double yAll[kMaxK*4];
        for (int i=0;i<kMaxK;i++){
            yAll[i]        = sigGV[i];
            yAll[kMaxK+i]  = sigGT[i];
            yAll[2*kMaxK+i]= sigIV[i];
            yAll[3*kMaxK+i]= sigIT[i];
        }
        double ymax = *std::max_element(yAll,yAll+4*kMaxK)*1.20;

        gV->SetTitle(Form(";k-th photon rank (1 = first);#sigma(t_{k}) [ps]  (%s)",capName));
        gV->GetYaxis()->SetRangeUser(0, ymax);
        gV->GetXaxis()->SetTitleSize(0.048);
        gV->GetYaxis()->SetTitleSize(0.048);
        gV->Draw("APE");
        gT->Draw("PE SAME");
        gVI->Draw("L SAME");
        gTI->Draw("L SAME");

        TLegend* leg = new TLegend(0.58,0.60,0.93,0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
        leg->AddEntry(gV,  "Vikuiti: #sigma_{Gauss}","pe");
        leg->AddEntry(gVI, "Vikuiti: #sigma_{IQR}","l");
        leg->AddEntry(gT,  "TIR-only: #sigma_{Gauss}","pe");
        leg->AddEntry(gTI, "TIR-only: #sigma_{IQR}","l");
        leg->Draw();

        // ── Mark k_opt on the plot ────────────────────────────────────────────
        // Vikuiti minimum: vertical dashed line + arrow + label
        {
            TLine* lV2 = new TLine(kOptV_f, 0, kOptV_f, sigOptV_f + 0.8);
            lV2->SetLineColor(colV); lV2->SetLineStyle(2); lV2->SetLineWidth(2);
            lV2->Draw();
            TArrow* arV = new TArrow(kOptV_f, sigOptV_f+1.5, kOptV_f, sigOptV_f+0.15,
                                     0.015, "|>");
            arV->SetFillColor(colV); arV->SetLineColor(colV); arV->SetLineWidth(2);
            arV->Draw();
            TLatex lxV; lxV.SetTextSize(0.032); lxV.SetTextColor(colV);
            lxV.SetTextAlign(22);
            lxV.DrawLatex(kOptV_f, sigOptV_f+2.1,
                Form("k_{opt}^{V}=%.0f, #sigma=%.1fps", kOptV_f, sigOptV_f));
        }
        // TIR-only minimum: same treatment shifted slightly to avoid overlap
        {
            double offT = (std::abs(kOptT_f - kOptV_f) < 3) ? 3.5 : 0.0;
            TLine* lT2 = new TLine(kOptT_f, 0, kOptT_f, sigOptT_f + 0.8);
            lT2->SetLineColor(colT); lT2->SetLineStyle(2); lT2->SetLineWidth(2);
            lT2->Draw();
            TArrow* arT = new TArrow(kOptT_f, sigOptT_f+1.5+offT,
                                     kOptT_f, sigOptT_f+0.15,
                                     0.015, "|>");
            arT->SetFillColor(colT); arT->SetLineColor(colT); arT->SetLineWidth(2);
            arT->Draw();
            TLatex lxT; lxT.SetTextSize(0.032); lxT.SetTextColor(colT+1);
            lxT.SetTextAlign(22);
            lxT.DrawLatex(kOptT_f, sigOptT_f+2.1+offT,
                Form("k_{opt}^{T}=%.0f, #sigma=%.1fps", kOptT_f, sigOptT_f));
        }

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.033);
        lat.DrawLatex(0.13,0.93,
            "EJ-228 cylinder  |  1 GeV #mu^{-}  |  5000 events  |  Gauss fit on core (#pm2#sigma_{68%})");

        c1->SaveAs("timing_sigma_vs_k.pdf");
        c1->SaveAs("timing_sigma_vs_k.png");
        printf("  Saved: timing_sigma_vs_k.pdf/png\n");
    }

    // ── Canvas 2: 6-panel distributions with Gaussian fits ───────────────────
    {
        TCanvas* c2 = new TCanvas("c_fits","Timing distributions", 1500, 900);
        c2->Divide(3, 2, 0.004, 0.005);

        for (int ip = 0; ip < (int)kShowK.size(); ++ip) {
            int k = kShowK[ip];
            c2->cd(ip+1);
            gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.16);
            gPad->SetTopMargin(0.10);

            TH1D* hV = hShowV[k];
            TH1D* hT = hShowT[k];
            if (!hV || !hT) continue;

            TF1* fV2 = hV->GetFunction(("hV_k"+std::to_string(k)+"_gaus").c_str());
            TF1* fT2 = hT->GetFunction(("hT_k"+std::to_string(k)+"_gaus").c_str());

            double ymax2 = std::max(hV->GetMaximum(), hT->GetMaximum()) * 1.25;
            hV->SetMaximum(ymax2); hT->SetMaximum(ymax2);
            hV->SetMinimum(0);

            // Align x-axes: use common range
            double xlo = std::min(hV->GetXaxis()->GetXmin(), hT->GetXaxis()->GetXmin());
            double xhi = std::max(hV->GetXaxis()->GetXmax(), hT->GetXaxis()->GetXmax());
            hV->GetXaxis()->SetRangeUser(xlo, xhi);

            hV->SetTitle(Form("k = %d;#Deltat_{k} [ps];Events / bin",k));
            hV->GetXaxis()->SetTitleSize(0.065);
            hV->GetYaxis()->SetTitleSize(0.060);
            hV->GetXaxis()->SetLabelSize(0.055);
            hV->GetYaxis()->SetLabelSize(0.053);
            hV->GetYaxis()->SetTitleOffset(1.3);
            hV->Draw("HIST");
            hT->Draw("HIST SAME");

            if (fV2) {
                fV2->SetLineWidth(2);
                fV2->SetLineColor(kBlue+3);
                fV2->SetLineStyle(1);
                fV2->Draw("SAME");
            }
            if (fT2) {
                fT2->SetLineWidth(2);
                fT2->SetLineColor(kOrange+4);
                fT2->SetLineStyle(2);
                fT2->Draw("SAME");
            }

            // σ labels (top-left box)
            TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.058);
            lat2.SetTextColor(kBlue+3);
            lat2.DrawLatex(0.20, 0.88, Form("#sigma_{V} = %.1f ps", sigGV[k-1]));
            lat2.SetTextColor(kOrange+4);
            lat2.DrawLatex(0.20, 0.79, Form("#sigma_{T} = %.1f ps", sigGT[k-1]));

            if (ip == 0) {
                TLegend* lg = new TLegend(0.52,0.72,0.94,0.88);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.053);
                lg->AddEntry(hV,"Vikuiti","f");
                lg->AddEntry(hT,"TIR-only","f");
                lg->Draw();
            }
        }

        // Overall title
        c2->cd(0);
        TLatex titLat; titLat.SetNDC(); titLat.SetTextSize(0.020);
        titLat.SetTextAlign(22);
        titLat.DrawLatex(0.5, 0.995,
            Form("t_{k} distributions | EJ-228 cylinder | %s | 5000 events | "
                 "Gauss fit on \\pm2#sigma_{68%%} core", capName));

        c2->SaveAs("timing_fits_panels.pdf");
        c2->SaveAs("timing_fits_panels.png");
        printf("  Saved: timing_fits_panels.pdf/png\n");
    }

    fV->Close(); fT->Close();
    printf("\nDone.\n\n");
}
