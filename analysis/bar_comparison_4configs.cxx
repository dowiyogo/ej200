// bar_comparison_4configs.cxx
// Unified comparison: TOP vs END readout, TIR vs Vikuiti
// Scintillators: EJ-204 (OPSC-101) and EJ-230 (OPSC-106)
//
// Data sources:
//  TOP+TIR,  TOP+Vik(skin): results/scan_resolution/  photon_hits_{ej204,ej230}_{tir,vik}_x*.root
//  END+Vik(air-gap):        results/scan_end_vikuiti/  photon_hits_{ej204,ej230}_x*.root
//  END+TIR(end-only):       results/scan_end_tir/      photon_hits_{ej204,ej230}_x*.root
//
// Usage: cd /home/rrios/ej200 && root -l -b -q 'analysis/bar_comparison_4configs.cxx'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <numeric>
#include <iostream>
#include <sstream>
#include <iomanip>

// ── Configuration ─────────────────────────────────────────────────────────────

static const std::vector<int> kPositions =
    {-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600};

static const int kKtop = 2;   // k-th photon for TOP analysis (same as original scan)
static const int kKend = 1;   // k-th photon for END analysis

struct Config {
    std::string id;       // unique key
    std::string label;    // plot legend
    std::string readout;  // "Top" or "End"
    std::string reflector;// "TIR" or "Vikuiti"
    std::string scint;    // "EJ-204" or "EJ-230"
    std::string rootdir;
    std::string tag;      // filename tag fragment
    int faceL, faceR;     // face_type for left/right (or -1 if single face)
    int faceTop;          // face_type for top (or -1)
    int color;
    int marker;
    int linestle;
};

// ── Gaussian fit from IQR ─────────────────────────────────────────────────────

bool gauss_fit(const std::vector<double>& vals, double& sigma, double& sigma_err,
               double& mean, double& mean_err) {
    if ((int)vals.size() < 15) return false;
    std::vector<double> v = vals;
    std::sort(v.begin(), v.end());
    double q25 = v[(int)(0.25*v.size())];
    double q75 = v[(int)(0.75*v.size())];
    double med  = v[v.size()/2];
    double sig  = (q75-q25)/1.349;
    if (sig <= 0) return false;
    TH1D h("hg","",100, med-4*sig, med+4*sig);
    for (double x : v) h.Fill(x);
    TF1 gaus("gf","gaus", med-2.5*sig, med+2.5*sig);
    gaus.SetParameters(h.GetMaximum(), med, sig);
    int st = h.Fit(&gaus,"RQN");
    if (st != 0 || gaus.GetParameter(2) <= 0) return false;
    sigma     = gaus.GetParameter(2) * 1e3;    // ns → ps
    sigma_err = gaus.GetParError(2)  * 1e3;
    mean      = gaus.GetParameter(1);           // ns
    mean_err  = gaus.GetParError(1);
    return true;
}

// ── Result per position ────────────────────────────────────────────────────────

struct PosResult {
    double x_mm      = 0;
    double sigma_ps  = 0;    // σ(T0) in ps  [or σ(k-th photon) for top]
    double sigma_err = 0;
    double npe       = 0;    // mean photons per event (L+R for end, total for top)
    bool   valid     = false;
};

// ── Analyse one file for END readout ─────────────────────────────────────────

PosResult analyse_end(const std::string& path, int x_mm, int kK) {
    PosResult res; res.x_mm = x_mm;
    TFile* f = TFile::Open(path.c_str(),"READ");
    if (!f || f->IsZombie()) { if(f){f->Close();delete f;} return res; }
    TTree* tree = nullptr;
    f->GetObject("sipm_hits", tree);
    if (!tree) { f->Close(); delete f; return res; }

    int ev_id=0, face=0; double tns=0;
    tree->SetBranchAddress("event_id",&ev_id);
    tree->SetBranchAddress("face_type",&face);
    tree->SetBranchAddress("time_ns",&tns);

    std::map<int,std::vector<double>> lhits, rhits;
    const Long64_t nent = tree->GetEntries();
    for (Long64_t i=0;i<nent;++i) {
        tree->GetEntry(i);
        if (face==0) lhits[ev_id].push_back(tns);
        else if (face==1) rhits[ev_id].push_back(tns);
    }
    f->Close(); delete f;

    long long nL=0, nR=0;
    for (auto& [e,v]:lhits) nL += v.size();
    for (auto& [e,v]:rhits) nR += v.size();

    // estimate total events from lhits + rhits map union
    std::set<int> all_evs;
    for (auto& [e,v]:lhits) all_evs.insert(e);
    for (auto& [e,v]:rhits) all_evs.insert(e);
    double nevents = all_evs.empty() ? 1 : all_evs.size();

    res.npe = (nL + nR) / nevents;

    std::vector<double> T0_vals;
    for (auto& [eid, lv] : lhits) {
        auto ir = rhits.find(eid);
        if (ir == rhits.end()) continue;
        auto& rv = ir->second;
        if ((int)lv.size()<kK || (int)rv.size()<kK) continue;
        std::nth_element(lv.begin(), lv.begin()+kK-1, lv.end());
        std::nth_element(rv.begin(), rv.begin()+kK-1, rv.end());
        T0_vals.push_back((lv[kK-1]+rv[kK-1])/2.0);
    }

    double mean=0, merr=0;
    if (gauss_fit(T0_vals, res.sigma_ps, res.sigma_err, mean, merr))
        res.valid = true;
    else if (!T0_vals.empty()) {
        // fallback: RMS
        double s=0; for(double v:T0_vals) s+=v;
        mean = s/T0_vals.size();
        double rms=0; for(double v:T0_vals) rms+=(v-mean)*(v-mean);
        res.sigma_ps  = std::sqrt(rms/T0_vals.size())*1e3;
        res.sigma_err = 0;
        res.valid = (T0_vals.size() >= 10);
    }
    return res;
}

// ── Analyse one file for TOP readout ─────────────────────────────────────────

PosResult analyse_top(const std::string& path, int x_mm, int kK) {
    PosResult res; res.x_mm = x_mm;
    TFile* f = TFile::Open(path.c_str(),"READ");
    if (!f || f->IsZombie()) { if(f){f->Close();delete f;} return res; }
    TTree* tree = nullptr;
    f->GetObject("sipm_hits", tree);
    if (!tree) { f->Close(); delete f; return res; }

    int ev_id=0, face=0; double tns=0;
    tree->SetBranchAddress("event_id",&ev_id);
    tree->SetBranchAddress("face_type",&face);
    tree->SetBranchAddress("time_ns",&tns);

    std::map<int,std::vector<double>> top_hits;
    long long ntop=0;
    const Long64_t nent = tree->GetEntries();
    for (Long64_t i=0;i<nent;++i) {
        tree->GetEntry(i);
        if (face==2) { top_hits[ev_id].push_back(tns); ntop++; }
    }
    f->Close(); delete f;

    double nevents = top_hits.empty() ? 1 : top_hits.size();
    res.npe = ntop / nevents;

    std::vector<double> tk_vals;
    for (auto& [eid, tv] : top_hits) {
        if ((int)tv.size() < kK) continue;
        std::nth_element(tv.begin(), tv.begin()+kK-1, tv.end());
        tk_vals.push_back(tv[kK-1]);
    }

    double mean=0, merr=0;
    if (gauss_fit(tk_vals, res.sigma_ps, res.sigma_err, mean, merr))
        res.valid = true;
    return res;
}

// ── Style ─────────────────────────────────────────────────────────────────────

void apply_style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.07); gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.13); gStyle->SetPadLeftMargin(0.14);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

// ── Main ─────────────────────────────────────────────────────────────────────

void bar_comparison_4configs() {
    apply_style();

    // Config: {id, label, readout, reflector, scint, rootdir, tag,
    //          faceL, faceR, faceTop, color, marker, linestyle}
    // For TOP:  faceL=-1, faceR=-1, faceTop=2
    // For END:  faceL=0,  faceR=1,  faceTop=-1

    struct Cfg {
        std::string id, label, readout, reflector, scint, rootdir, filetag;
        bool isEnd;
        int color, marker, lstyle;
    };

    std::vector<Cfg> configs = {
        // EJ-204
        {"ej204_top_tir", "EJ-204 TOP+TIR",     "Top","TIR","EJ-204",
         "results/scan_resolution","ej204_tir",   false, kBlue+1,   24, 2},
        {"ej204_top_vik", "EJ-204 TOP+Vikuiti",  "Top","Vik","EJ-204",
         "results/scan_resolution","ej204_vik",   false, kBlue+1,   20, 1},
        {"ej204_end_tir", "EJ-204 END+TIR",      "End","TIR","EJ-204",
         "results/scan_end_tir","ej204",           true,  kAzure+2,  26, 2},
        {"ej204_end_vik", "EJ-204 END+Vikuiti",  "End","Vik","EJ-204",
         "results/scan_end_vikuiti","ej204",       true,  kAzure+2,  22, 1},
        // EJ-230
        {"ej230_top_tir", "EJ-230 TOP+TIR",     "Top","TIR","EJ-230",
         "results/scan_resolution","ej230_tir",   false, kRed+1,    25, 2},
        {"ej230_top_vik", "EJ-230 TOP+Vikuiti",  "Top","Vik","EJ-230",
         "results/scan_resolution","ej230_vik",   false, kRed+1,    21, 1},
        {"ej230_end_tir", "EJ-230 END+TIR",      "End","TIR","EJ-230",
         "results/scan_end_tir","ej230",           true,  kOrange+1, 27, 2},
        {"ej230_end_vik", "EJ-230 END+Vikuiti",  "End","Vik","EJ-230",
         "results/scan_end_vikuiti","ej230",       true,  kOrange+1, 23, 1},
    };

    struct CfgResult { Cfg cfg; std::vector<PosResult> pts; };
    std::vector<CfgResult> all;

    for (const auto& cfg : configs) {
        CfgResult cr; cr.cfg = cfg;
        std::cout << "\n=== " << cfg.label << " ===\n";
        bool missing_any = false;
        for (int x : kPositions) {
            std::ostringstream path;
            path << cfg.rootdir << "/photon_hits_" << cfg.filetag << "_x" << x << ".root";
            PosResult res;
            if (cfg.isEnd)
                res = analyse_end(path.str(), x, kKend);
            else
                res = analyse_top(path.str(), x, kKtop);
            if (!res.valid && res.npe == 0) { missing_any = true; }
            std::cout << "  x=" << std::setw(5) << x
                      << "mm  σ(T0)=" << std::fixed << std::setprecision(1)
                      << res.sigma_ps << " ps"
                      << "  Npe=" << std::setprecision(1) << res.npe
                      << (res.valid ? "" : "  [low-stat]") << "\n";
            cr.pts.push_back(res);
        }
        if (missing_any) std::cout << "  WARNING: some files missing for " << cfg.label << "\n";
        all.push_back(cr);
    }

    const std::string kOut = "build_t0minidaq";

    // ── Plot 1: σ(T0) vs x ────────────────────────────────────────────────────
    {
        TCanvas c("cSigma","",1000,620); c.SetGrid();
        // two-panel: EJ-204 left, EJ-230 right
        c.Divide(2,1,0.002,0.002);

        for (int panel=0; panel<2; panel++) {
            c.cd(panel+1);
            gPad->SetGrid();
            const char* sc = (panel==0)?"EJ-204":"EJ-230";
            auto* fr = gPad->DrawFrame(-680,0,680,500);
            fr->GetXaxis()->SetTitle("Position  x  (mm)");
            fr->GetYaxis()->SetTitle("#sigma(T_{0})  (ps)");
            fr->GetXaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetTitleSize(0.052);
            fr->GetXaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetLabelSize(0.044);
            TLatex tit; tit.SetNDC(); tit.SetTextFont(42); tit.SetTextSize(0.055);
            tit.DrawLatex(0.18,0.91,sc);

            auto* leg = new TLegend(0.17,0.62,0.70,0.89);
            leg->SetTextFont(42); leg->SetTextSize(0.040);

            for (auto& cr : all) {
                if (cr.cfg.scint != sc) continue;
                std::vector<double> xs,ys,yes;
                for (auto& p : cr.pts) if (p.valid && p.sigma_ps>0) {
                    xs.push_back(p.x_mm); ys.push_back(p.sigma_ps); yes.push_back(p.sigma_err);
                }
                if (xs.empty()) continue;
                auto* g = new TGraphErrors((int)xs.size(),xs.data(),ys.data(),nullptr,yes.data());
                g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
                g->SetLineStyle(cr.cfg.lstyle); g->SetLineWidth(2);
                g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(0.9);
                g->Draw("LP same");
                leg->AddEntry(g, cr.cfg.label.c_str(), "lp");
            }
            leg->Draw();
        }
        c.SaveAs((kOut+"/comp_sigma_T0_vs_x.png").c_str());
        c.SaveAs((kOut+"/comp_sigma_T0_vs_x.pdf").c_str());
        std::cout << "\nSaved comp_sigma_T0_vs_x\n";
    }

    // ── Plot 2: Npe vs x ─────────────────────────────────────────────────────
    {
        TCanvas c("cNpe","",1000,620); c.SetGrid();
        c.Divide(2,1,0.002,0.002);

        for (int panel=0; panel<2; panel++) {
            c.cd(panel+1);
            gPad->SetGrid();
            gPad->SetLogy();
            const char* sc = (panel==0)?"EJ-204":"EJ-230";
            auto* fr = gPad->DrawFrame(-680,0.05,680,5000);
            fr->GetXaxis()->SetTitle("Position  x  (mm)");
            fr->GetYaxis()->SetTitle("Mean N_{pe} per event");
            fr->GetXaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetTitleSize(0.052);
            fr->GetXaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetLabelSize(0.044);
            TLatex tit; tit.SetNDC(); tit.SetTextFont(42); tit.SetTextSize(0.055);
            tit.DrawLatex(0.18,0.91,sc);

            auto* leg = new TLegend(0.17,0.12,0.72,0.42);
            leg->SetTextFont(42); leg->SetTextSize(0.040);

            for (auto& cr : all) {
                if (cr.cfg.scint != sc) continue;
                std::vector<double> xs,ys;
                for (auto& p : cr.pts) { xs.push_back(p.x_mm); ys.push_back(std::max(p.npe,0.01)); }
                auto* g = new TGraph((int)xs.size(),xs.data(),ys.data());
                g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
                g->SetLineStyle(cr.cfg.lstyle); g->SetLineWidth(2);
                g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(0.9);
                g->Draw("LP same");
                leg->AddEntry(g, cr.cfg.label.c_str(), "lp");
            }
            leg->Draw();
        }
        c.SaveAs((kOut+"/comp_npe_vs_x.png").c_str());
        c.SaveAs((kOut+"/comp_npe_vs_x.pdf").c_str());
        std::cout << "Saved comp_npe_vs_x\n";
    }

    // ── Summary table ─────────────────────────────────────────────────────────
    std::cout << "\n";
    std::cout << std::string(80,'=') << "\n";
    std::cout << "SUMMARY at x=0\n";
    std::cout << std::left
              << std::setw(28) << "Config"
              << std::setw(12) << "Npe@x=0"
              << std::setw(16) << "σ(T0)@x=0 (ps)"
              << "\n";
    std::cout << std::string(56,'-') << "\n";
    for (auto& cr : all) {
        auto it = std::find_if(cr.pts.begin(), cr.pts.end(),
                               [](const PosResult& p){ return p.x_mm == 0; });
        if (it == cr.pts.end()) continue;
        std::cout << std::left << std::setw(28) << cr.cfg.label
                  << std::setw(12) << std::fixed << std::setprecision(0) << it->npe
                  << (it->valid ? std::to_string((int)std::round(it->sigma_ps))+" ps"
                                : "< stat")
                  << "\n";
    }
    std::cout << std::string(80,'=') << "\n";
    std::cout << "All plots saved to " << kOut << "/\n";
}
