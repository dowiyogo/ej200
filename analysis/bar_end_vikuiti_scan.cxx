// bar_end_vikuiti_scan.cxx
// Analysis of End-only Vikuiti scan (EJ-200, EJ-204, EJ-230)
//
// Per position and scintillator:
//   (1) σ(T0) — timing resolution: T0 = (t_L + t_R)/2, position-independent
//   (2) σ(Δt) — time difference t_R - t_L → position resolution
//   (3) v_eff from mean(Δt) vs x linear fit
//   (4) Mean N_pe at each end
//
// Uses k-th first photon at left/right ends per event.
//
// Usage:  root -l -b -q 'bar_end_vikuiti_scan.cxx'

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

static const std::vector<int> kPositions = {-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600};
static const int kKend = 1;   // k-th photon at each end (1=first)

struct Config {
    std::string tag;
    std::string label;
    int color;
    int marker;
};

static const std::vector<Config> kConfigs = {
    {"ej200", "EJ-200", kGreen+2,   20},
    {"ej204", "EJ-204", kBlue+1,    21},
    {"ej230", "EJ-230", kRed+1,     22},
};

static const std::string kRootDir = "results/scan_end_vikuiti";
static const std::string kOutDir  = "build_end_vikuiti";

// ── Per-position result ────────────────────────────────────────────────────────

struct PosResult {
    double x_mm      = 0;
    double sigma_T0_ps = 0;    // σ((t_L+t_R)/2) in ps
    double sigma_T0_err = 0;
    double sigma_dt_ps = 0;    // σ(t_R - t_L) in ps → position resolution
    double sigma_dt_err = 0;
    double mean_dt_ns  = 0;    // mean(t_R - t_L) for v_eff
    double mean_dt_err = 0;
    double npe_left    = 0;    // mean photons at left end per event
    double npe_right   = 0;
};

// ── Helper: Gaussian fit from IQR ─────────────────────────────────────────────

bool gauss_fit(const std::vector<double>& vals, double& sigma, double& sigma_err) {
    if (vals.size() < 20) return false;
    std::vector<double> v = vals;
    std::sort(v.begin(), v.end());
    double q25 = v[(int)(0.25 * v.size())];
    double q75 = v[(int)(0.75 * v.size())];
    double med  = v[v.size() / 2];
    double sig  = (q75 - q25) / 1.349;
    if (sig <= 0) return false;
    TH1D h("hg","",120, med - 4*sig, med + 4*sig);
    for (double x : v) h.Fill(x);
    TF1 gaus("gf","gaus", med - 2.5*sig, med + 2.5*sig);
    gaus.SetParameters(h.GetMaximum(), med, sig);
    int st = h.Fit(&gaus, "RQN");
    if (st != 0 || gaus.GetParameter(2) <= 0) return false;
    sigma     = gaus.GetParameter(2) * 1e3;   // ns → ps
    sigma_err = gaus.GetParError(2)  * 1e3;
    return true;
}

// ── Analyse one ROOT file ──────────────────────────────────────────────────────

PosResult analyse_file(const std::string& path, int x_mm) {
    PosResult res;
    res.x_mm = x_mm;

    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "  MISSING: " << path << "\n";
        if (f) { f->Close(); delete f; }
        return res;
    }
    TTree* tree = nullptr;
    f->GetObject("sipm_hits", tree);
    if (!tree) { std::cerr << "  No tree in: " << path << "\n"; f->Close(); delete f; return res; }

    int    ev_id   = 0;
    int    face    = 0;
    double time_ns = 0;
    tree->SetBranchAddress("event_id", &ev_id);
    tree->SetBranchAddress("face_type", &face);
    tree->SetBranchAddress("time_ns",   &time_ns);

    // Collect k-th photon per event per side
    std::map<int, std::vector<double>> left_hits, right_hits;

    const Long64_t nent = tree->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);
        if (face == 0) left_hits[ev_id].push_back(time_ns);
        else if (face == 1) right_hits[ev_id].push_back(time_ns);
    }
    f->Close(); delete f;

    // Build T0 and Δt distributions
    std::vector<double> T0_vals, dt_vals;
    long long n_left_sum = 0, n_right_sum = 0;

    for (auto& [eid, lhits] : left_hits) {
        n_left_sum += (int)lhits.size();
        auto ir = right_hits.find(eid);
        if (ir == right_hits.end()) continue;
        auto& rhits = ir->second;
        n_right_sum += (int)rhits.size();

        if ((int)lhits.size() < kKend || (int)rhits.size() < kKend) continue;

        std::nth_element(lhits.begin(), lhits.begin() + kKend - 1, lhits.end());
        std::nth_element(rhits.begin(), rhits.begin() + kKend - 1, rhits.end());

        double tL = lhits[kKend - 1];
        double tR = rhits[kKend - 1];
        T0_vals.push_back((tL + tR) / 2.0);
        dt_vals.push_back(tR - tL);
    }

    res.npe_left  = left_hits.empty()  ? 0 : (double)n_left_sum  / left_hits.size();
    res.npe_right = right_hits.empty() ? 0 : (double)n_right_sum / right_hits.size();

    gauss_fit(T0_vals, res.sigma_T0_ps, res.sigma_T0_err);
    gauss_fit(dt_vals, res.sigma_dt_ps, res.sigma_dt_err);

    if (!dt_vals.empty()) {
        double sum = std::accumulate(dt_vals.begin(), dt_vals.end(), 0.0);
        double mean = sum / dt_vals.size();
        double var = 0;
        for (double d : dt_vals) var += (d-mean)*(d-mean);
        res.mean_dt_ns  = mean;
        res.mean_dt_err = std::sqrt(var/dt_vals.size()) / std::sqrt((double)dt_vals.size());
    }

    return res;
}

// ── Style ─────────────────────────────────────────────────────────────────────

void apply_style() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.06); gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.12); gStyle->SetPadLeftMargin(0.13);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

void bar_label(TCanvas* c) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(0.031); l.SetTextColor(kGray+1);
    l.DrawLatex(0.14, 0.965,
        "End-only + Vikuiti 3M  |  1400#times60#times10 mm^{3}  |  5000 ev/pt  |  no jitter");
}

// ── Main ─────────────────────────────────────────────────────────────────────

void bar_end_vikuiti_scan() {
    apply_style();

    struct ConfigResult { Config cfg; std::vector<PosResult> pts; };
    std::vector<ConfigResult> all;

    for (const auto& cfg : kConfigs) {
        ConfigResult cr; cr.cfg = cfg;
        std::cout << "\n=== " << cfg.label << " ===\n";
        for (int x : kPositions) {
            std::ostringstream path;
            path << kRootDir << "/photon_hits_" << cfg.tag << "_x" << x << ".root";
            auto res = analyse_file(path.str(), x);
            std::cout << "  x=" << std::setw(5) << x
                      << "mm  σ(T0)=" << std::fixed << std::setprecision(1) << res.sigma_T0_ps
                      << " ps  σ(Δt)=" << res.sigma_dt_ps
                      << " ps  Δt=" << res.mean_dt_ns << " ns"
                      << "  Npe=" << (int)(res.npe_left + res.npe_right) << "\n";
            cr.pts.push_back(res);
        }
        all.push_back(cr);
    }

    // ── Plot 1: σ(T0) vs x ────────────────────────────────────────────────────
    {
        TCanvas c("cT0","",800,560); c.SetGrid();
        auto* fr = c.DrawFrame(-680,0,680,120);
        fr->GetXaxis()->SetTitle("Position  x  (mm)");
        fr->GetYaxis()->SetTitle(Form("#sigma(T_{0})  [k=%d first photon]  (ps)", kKend));
        fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
        fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
        auto* leg = new TLegend(0.55,0.73,0.95,0.93);
        leg->SetTextFont(42); leg->SetTextSize(0.038);
        for (auto& cr : all) {
            std::vector<double> xs,ys,yes;
            for (auto& p : cr.pts) if (p.sigma_T0_ps>0) {
                xs.push_back(p.x_mm); ys.push_back(p.sigma_T0_ps); yes.push_back(p.sigma_T0_err);
            }
            if (xs.empty()) continue;
            auto* g = new TGraphErrors((int)xs.size(),xs.data(),ys.data(),nullptr,yes.data());
            g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2); g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(1.0);
            g->Draw("LP same");
            leg->AddEntry(g,cr.cfg.label.c_str(),"lp");
        }
        leg->Draw(); bar_label(&c);
        c.SaveAs((kOutDir+"/end_vik_sigma_T0_vs_x.png").c_str());
        c.SaveAs((kOutDir+"/end_vik_sigma_T0_vs_x.pdf").c_str());
        std::cout<<"\nSaved end_vik_sigma_T0_vs_x\n";
    }

    // ── Plot 2: σ(Δt) vs x → position resolution ──────────────────────────────
    {
        TCanvas c("cDt","",800,560); c.SetGrid();
        auto* fr = c.DrawFrame(-680,0,680,600);
        fr->GetXaxis()->SetTitle("Position  x  (mm)");
        fr->GetYaxis()->SetTitle(Form("#sigma(#Deltat)  [k=%d]  (ps)", kKend));
        fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
        fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
        auto* leg = new TLegend(0.55,0.73,0.95,0.93);
        leg->SetTextFont(42); leg->SetTextSize(0.038);
        for (auto& cr : all) {
            std::vector<double> xs,ys,yes;
            for (auto& p : cr.pts) if (p.sigma_dt_ps>0) {
                xs.push_back(p.x_mm); ys.push_back(p.sigma_dt_ps); yes.push_back(p.sigma_dt_err);
            }
            if (xs.empty()) continue;
            auto* g = new TGraphErrors((int)xs.size(),xs.data(),ys.data(),nullptr,yes.data());
            g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2); g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(1.0);
            g->Draw("LP same");
            leg->AddEntry(g,cr.cfg.label.c_str(),"lp");
        }
        leg->Draw(); bar_label(&c);
        c.SaveAs((kOutDir+"/end_vik_sigma_dt_vs_x.png").c_str());
        c.SaveAs((kOutDir+"/end_vik_sigma_dt_vs_x.pdf").c_str());
        std::cout<<"Saved end_vik_sigma_dt_vs_x\n";
    }

    // ── Plot 3: Mean Δt vs x → group velocity ────────────────────────────────
    {
        TCanvas c("cGV","",800,560); c.SetGrid();
        auto* fr = c.DrawFrame(-680,-15,680,15);
        fr->GetXaxis()->SetTitle("Position  x  (mm)");
        fr->GetYaxis()->SetTitle("#Deltat = t_{R}^{(k)} #minus t_{L}^{(k)}  (ns)");
        fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
        fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
        auto* leg = new TLegend(0.14,0.73,0.55,0.93);
        leg->SetTextFont(42); leg->SetTextSize(0.038);
        TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.030);
        double ytxt = 0.62;
        for (auto& cr : all) {
            std::vector<double> xs,ys,yes;
            for (auto& p : cr.pts) if (p.npe_left>0&&p.npe_right>0) {
                xs.push_back(p.x_mm); ys.push_back(p.mean_dt_ns); yes.push_back(p.mean_dt_err);
            }
            if (xs.size()<3) continue;
            auto* g = new TGraphErrors((int)xs.size(),xs.data(),ys.data(),nullptr,yes.data());
            g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2); g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(1.0);
            g->Draw("LP same");
            TF1 lin("lin","[0]*x+[1]",-660,660);
            g->Fit(&lin,"QN");
            double slope = lin.GetParameter(0);
            double veff  = -2.0/slope;
            double neff  = 299.792/veff;
            lin.SetLineColor(cr.cfg.color); lin.SetLineStyle(2); lin.SetLineWidth(1);
            lin.DrawClone("same");
            leg->AddEntry(g,cr.cfg.label.c_str(),"lp");
            lat.SetNDC(true); lat.SetTextColor(cr.cfg.color);
            lat.DrawLatex(0.14,ytxt,
                Form("%s: v_{eff}=%.0f mm/ns  n_{eff}=%.2f",
                     cr.cfg.label.c_str(),veff,neff));
            ytxt -= 0.055;
        }
        leg->Draw(); bar_label(&c);
        c.SaveAs((kOutDir+"/end_vik_groupvel.png").c_str());
        c.SaveAs((kOutDir+"/end_vik_groupvel.pdf").c_str());
        std::cout<<"Saved end_vik_groupvel\n";
    }

    // ── Plot 4: Mean Npe vs x ────────────────────────────────────────────────
    {
        TCanvas c("cNpe","",800,560); c.SetGrid();
        auto* fr = c.DrawFrame(-680,0,680,6000);
        fr->GetXaxis()->SetTitle("Position  x  (mm)");
        fr->GetYaxis()->SetTitle("Mean N_{pe} at end SiPMs (L+R) per event");
        fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
        fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
        auto* leg = new TLegend(0.55,0.73,0.95,0.93);
        leg->SetTextFont(42); leg->SetTextSize(0.038);
        for (auto& cr : all) {
            std::vector<double> xs,ys;
            for (auto& p : cr.pts) {
                xs.push_back(p.x_mm); ys.push_back(p.npe_left+p.npe_right);
            }
            auto* g = new TGraph((int)xs.size(),xs.data(),ys.data());
            g->SetLineColor(cr.cfg.color); g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2); g->SetMarkerStyle(cr.cfg.marker); g->SetMarkerSize(1.0);
            g->Draw("LP same");
            leg->AddEntry(g,cr.cfg.label.c_str(),"lp");
        }
        leg->Draw(); bar_label(&c);
        c.SaveAs((kOutDir+"/end_vik_npe_vs_x.png").c_str());
        c.SaveAs((kOutDir+"/end_vik_npe_vs_x.pdf").c_str());
        std::cout<<"Saved end_vik_npe_vs_x\n";
    }

    std::cout<<"\nAll plots saved to "<<kOutDir<<"/\n";
}
