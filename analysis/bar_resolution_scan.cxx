// bar_resolution_scan.cxx
// Analyse position-scan ROOT files to extract:
//   (1) Timing resolution σ(k) vs position — from top SiPMs (face_type=2)
//   (2) Group velocity  v_eff  — from end SiPMs (face_type=0,1)
//
// Usage:  root -l -b -q 'bar_resolution_scan.cxx'
//
// Input:  results/scan_resolution/photon_hits_{tag}_x{pos}.root
//         tag ∈ {ej204_tir, ej230_tir, ej204_vik, ej230_vik}
//         pos ∈ {-600,..,600} mm
//
// Output: build_t0minidaq/bar_scan_sigma_vs_x.{png,pdf}
//         build_t0minidaq/bar_scan_groupvel.{png,pdf}
//         build_t0minidaq/bar_scan_npe_vs_x.{png,pdf}

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

// ── Configuration ────────────────────────────────────────────────────────────

static const std::vector<int> kPositions = {-600,-500,-400,-300,-200,-100,0,100,200,300,400,500,600};
static const int kKopt = 2;   // k-th photon used for σ(k) — change if needed

struct Config {
    std::string tag;
    std::string label;
    int color;
    int style;     // line style (1=solid, 2=dashed)
};

static const std::vector<Config> kConfigs = {
    {"ej204_tir", "EJ-204  TIR",      kBlue+1,  1},
    {"ej230_tir", "EJ-230  TIR",      kRed+1,   1},
    {"ej204_vik", "EJ-204  Vikuiti",  kBlue+1,  2},
    {"ej230_vik", "EJ-230  Vikuiti",  kRed+1,   2},
};

static const std::string kRootDir  = "results/scan_resolution";
static const std::string kOutDir   = "build_t0minidaq";

// ── Per-position result ───────────────────────────────────────────────────────

struct PosResult {
    double x_mm   = 0;
    double sigma_ps = 0;   // σ(k) in picoseconds; 0 = failed/missing
    double sigma_err_ps = 0;
    double dt_ns  = 0;     // mean(t_R) – mean(t_L) in ns for end SiPMs
    double dt_err = 0;
    double npe_top = 0;    // mean photons reaching top SiPMs
    double npe_end = 0;    // mean photons reaching either end SiPM
};

// ── Helper: compute σ(k) from face_type==2 hits ──────────────────────────────

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
    f->GetObject("photon_hits", tree);
    if (!tree) {
        std::cerr << "  No tree in: " << path << "\n";
        f->Close(); delete f;
        return res;
    }

    int    ev_id    = 0;
    int    face     = 0;
    double time_ns  = 0;
    tree->SetBranchAddress("event_id", &ev_id);
    tree->SetBranchAddress("face_type", &face);
    tree->SetBranchAddress("time_ns",  &time_ns);

    // Collect hits per event
    std::map<int, std::vector<double>> top_hits;   // face_type==2
    std::map<int, double>              t_left;      // face_type==0  first photon
    std::map<int, double>              t_right;     // face_type==1  first photon

    const Long64_t nent = tree->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);
        if (face == 2) {
            top_hits[ev_id].push_back(time_ns);
        } else if (face == 0) {
            auto it = t_left.find(ev_id);
            if (it == t_left.end() || time_ns < it->second) t_left[ev_id] = time_ns;
        } else if (face == 1) {
            auto it = t_right.find(ev_id);
            if (it == t_right.end() || time_ns < it->second) t_right[ev_id] = time_ns;
        }
    }
    f->Close(); delete f;

    // ── σ(k) from top SiPMs ──
    std::vector<double> k_times;
    k_times.reserve(top_hits.size());
    for (auto& [eid, hits] : top_hits) {
        if ((int)hits.size() < kKopt) continue;
        std::nth_element(hits.begin(), hits.begin() + kKopt - 1, hits.end());
        k_times.push_back(hits[kKopt - 1]);
    }

    res.npe_top = (double)k_times.size();   // events with >=k photons

    if (k_times.size() > 20) {
        // IQR-based Gaussian fit
        std::sort(k_times.begin(), k_times.end());
        double q25 = k_times[(int)(0.25 * k_times.size())];
        double q75 = k_times[(int)(0.75 * k_times.size())];
        double med  = k_times[k_times.size() / 2];
        double iqr  = q75 - q25;
        double sig_est = iqr / 1.349;   // IQR ≈ 1.349σ for Gaussian

        TH1D h("hk", "", 120, med - 4*sig_est, med + 4*sig_est);
        for (double t : k_times) h.Fill(t);
        TF1 gaus("gfit", "gaus", med - 2.5*sig_est, med + 2.5*sig_est);
        gaus.SetParameters(h.GetMaximum(), med, sig_est);
        int fitstat = h.Fit(&gaus, "RQN");
        if (fitstat == 0 && gaus.GetParameter(2) > 0) {
            res.sigma_ps    = gaus.GetParameter(2) * 1e3;   // ns → ps
            res.sigma_err_ps = gaus.GetParError(2) * 1e3;
        }
    }

    // ── Group velocity: Δt = t_right - t_left ──
    std::vector<double> dts;
    for (auto& [eid, tl] : t_left) {
        auto ir = t_right.find(eid);
        if (ir != t_right.end()) dts.push_back(ir->second - tl);
    }
    if (!dts.empty()) {
        double sum = std::accumulate(dts.begin(), dts.end(), 0.0);
        double mean = sum / dts.size();
        double var  = 0;
        for (double d : dts) var += (d - mean) * (d - mean);
        res.dt_ns  = mean;
        res.dt_err = std::sqrt(var / dts.size()) / std::sqrt((double)dts.size());
        res.npe_end = (double)dts.size();
    }

    return res;
}

// ── Style helpers ─────────────────────────────────────────────────────────────

void apply_style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
}

TLegend* make_legend(double x1, double y1, double x2, double y2) {
    auto* leg = new TLegend(x1, y1, x2, y2);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    return leg;
}

void label_bar(TCanvas* c) {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.033);
    lat.SetTextColor(kGray+1);
    lat.DrawLatex(0.14, 0.965, "SHiP EJ-bar  1400 #times 60 #times 10 mm^{3}   |   5000 ev/point   |   no electronics jitter");
}

// ── Main ─────────────────────────────────────────────────────────────────────

void bar_resolution_scan() {
    apply_style();

    // Load all results
    struct ConfigResult {
        Config cfg;
        std::vector<PosResult> pts;
    };
    std::vector<ConfigResult> all;

    for (const auto& cfg : kConfigs) {
        ConfigResult cr;
        cr.cfg = cfg;
        std::cout << "\n=== " << cfg.label << " ===\n";
        for (int x : kPositions) {
            std::ostringstream path;
            path << kRootDir << "/photon_hits_" << cfg.tag
                 << "_x" << x << ".root";
            auto res = analyse_file(path.str(), x);
            std::cout << "  x=" << std::setw(5) << x
                      << "mm  σ(" << kKopt << ")=" << std::fixed << std::setprecision(1)
                      << res.sigma_ps << " ps  Δt=" << res.dt_ns << " ns\n";
            cr.pts.push_back(res);
        }
        all.push_back(cr);
    }

    // ── Plot 1: σ(k) vs x ────────────────────────────────────────────────────
    {
        TCanvas c("cSigma", "", 800, 560);
        c.SetGrid();

        auto* frame = c.DrawFrame(-680, 0, 680, 500);
        frame->GetXaxis()->SetTitle("Position along bar  x  (mm)");
        frame->GetYaxis()->SetTitle(Form("#sigma(k=%d)  (ps)", kKopt));
        frame->GetXaxis()->SetTitleSize(0.048);
        frame->GetYaxis()->SetTitleSize(0.048);
        frame->GetXaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetLabelSize(0.040);

        auto* leg = make_legend(0.14, 0.73, 0.50, 0.93);

        for (auto& cr : all) {
            std::vector<double> xs, ys, yes;
            for (auto& p : cr.pts) {
                if (p.sigma_ps > 0) {
                    xs.push_back(p.x_mm);
                    ys.push_back(p.sigma_ps);
                    yes.push_back(p.sigma_err_ps);
                }
            }
            if (xs.empty()) continue;

            auto* g = new TGraphErrors((int)xs.size(), xs.data(), ys.data(),
                                        nullptr, yes.data());
            g->SetLineColor(cr.cfg.color);
            g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2);
            g->SetLineStyle(cr.cfg.style);
            g->SetMarkerStyle(cr.cfg.style == 1 ? 20 : 24);
            g->SetMarkerSize(0.9);
            g->Draw("LP same");
            leg->AddEntry(g, cr.cfg.label.c_str(), "lp");
        }
        leg->Draw();
        label_bar(&c);

        c.SaveAs((kOutDir + "/bar_scan_sigma_vs_x.png").c_str());
        c.SaveAs((kOutDir + "/bar_scan_sigma_vs_x.pdf").c_str());
        std::cout << "\nSaved bar_scan_sigma_vs_x.{png,pdf}\n";
    }

    // ── Plot 2: Group velocity — Δt vs x ────────────────────────────────────
    // Fit Δt = -2x/(v_eff) + b  [x in mm, Δt in ns]
    // → slope = -2/v_eff  (v_eff in mm/ns)
    // → v_eff = c/n_eff, n_eff = c / v_eff
    {
        TCanvas c("cGV", "", 800, 560);
        c.SetGrid();

        auto* frame = c.DrawFrame(-680, -15, 680, 15);
        frame->GetXaxis()->SetTitle("Position along bar  x  (mm)");
        frame->GetYaxis()->SetTitle("#Deltat = t_{R} #minus t_{L}  (ns)");
        frame->GetXaxis()->SetTitleSize(0.048);
        frame->GetYaxis()->SetTitleSize(0.048);
        frame->GetXaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetLabelSize(0.040);

        auto* leg = make_legend(0.14, 0.73, 0.55, 0.93);

        TLatex lat;
        lat.SetTextFont(42);
        lat.SetTextSize(0.030);
        double ytext = 0.62;

        for (auto& cr : all) {
            std::vector<double> xs, ys, yes;
            for (auto& p : cr.pts) {
                if (p.npe_end > 0) {
                    xs.push_back(p.x_mm);
                    ys.push_back(p.dt_ns);
                    yes.push_back(p.dt_err);
                }
            }
            if (xs.size() < 3) continue;

            auto* g = new TGraphErrors((int)xs.size(), xs.data(), ys.data(),
                                        nullptr, yes.data());
            g->SetLineColor(cr.cfg.color);
            g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2);
            g->SetLineStyle(cr.cfg.style);
            g->SetMarkerStyle(cr.cfg.style == 1 ? 20 : 24);
            g->SetMarkerSize(0.9);
            g->Draw("LP same");

            // Linear fit
            TF1 lin("lin", "[0]*x + [1]", -660, 660);
            g->Fit(&lin, "QN");
            double slope = lin.GetParameter(0);   // ns/mm
            // v_eff = -2 / slope  [mm/ns] (slope is negative)
            double veff_mm_ns = -2.0 / slope;
            double c_mm_ns    = 299.792;           // speed of light in mm/ns
            double n_eff      = c_mm_ns / veff_mm_ns;

            lin.SetLineColor(cr.cfg.color);
            lin.SetLineStyle(cr.cfg.style);
            lin.SetLineWidth(1);
            lin.DrawClone("same");

            std::ostringstream label;
            label << cr.cfg.label << ":  v_{eff} = "
                  << std::fixed << std::setprecision(0) << veff_mm_ns
                  << " mm/ns  (n_{eff} = "
                  << std::fixed << std::setprecision(2) << n_eff << ")";
            leg->AddEntry(g, cr.cfg.label.c_str(), "lp");
            lat.SetNDC(true);
            lat.SetTextColor(cr.cfg.color);
            lat.DrawLatex(0.14, ytext,
                          Form("%s:  v_{eff} = %.0f mm/ns  (n_{eff} = %.2f)",
                               cr.cfg.label.c_str(), veff_mm_ns, n_eff));
            ytext -= 0.055;
        }
        leg->Draw();
        label_bar(&c);

        c.SaveAs((kOutDir + "/bar_scan_groupvel.png").c_str());
        c.SaveAs((kOutDir + "/bar_scan_groupvel.pdf").c_str());
        std::cout << "Saved bar_scan_groupvel.{png,pdf}\n";
    }

    // ── Plot 3: Mean Npe (top SiPMs) vs x ────────────────────────────────────
    {
        TCanvas c("cNpe", "", 800, 560);
        c.SetGrid();

        auto* frame = c.DrawFrame(-680, 0, 680, 800);
        frame->GetXaxis()->SetTitle("Position along bar  x  (mm)");
        frame->GetYaxis()->SetTitle("Events with #geq k top photons (/ 5000)");
        frame->GetXaxis()->SetTitleSize(0.048);
        frame->GetYaxis()->SetTitleSize(0.048);
        frame->GetXaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetLabelSize(0.040);

        auto* leg = make_legend(0.14, 0.73, 0.50, 0.93);

        for (auto& cr : all) {
            std::vector<double> xs, ys;
            for (auto& p : cr.pts) {
                xs.push_back(p.x_mm);
                ys.push_back(p.npe_top);
            }
            if (xs.empty()) continue;

            auto* g = new TGraph((int)xs.size(), xs.data(), ys.data());
            g->SetLineColor(cr.cfg.color);
            g->SetMarkerColor(cr.cfg.color);
            g->SetLineWidth(2);
            g->SetLineStyle(cr.cfg.style);
            g->SetMarkerStyle(cr.cfg.style == 1 ? 20 : 24);
            g->SetMarkerSize(0.9);
            g->Draw("LP same");
            leg->AddEntry(g, cr.cfg.label.c_str(), "lp");
        }
        leg->Draw();
        label_bar(&c);

        c.SaveAs((kOutDir + "/bar_scan_npe_vs_x.png").c_str());
        c.SaveAs((kOutDir + "/bar_scan_npe_vs_x.pdf").c_str());
        std::cout << "Saved bar_scan_npe_vs_x.{png,pdf}\n";
    }

    std::cout << "\nAll plots saved to " << kOutDir << "/\n";
}
