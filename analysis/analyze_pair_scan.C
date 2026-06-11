// analysis/analyze_pair_scan.C
// EXEC_08: per-position position-reconstruction from two adjacent top SiPMs.
//
// Reads per-photon sipm_hits TTrees from results/pairscan_<date>/*.root.
// For each position x_true, per event:
//   t_A = time of the N_PE_THRESH-th photon in channel ID_A  (leading-edge @ 4 PE)
//   t_B = time of the N_PE_THRESH-th photon in channel ID_B
//   DeltaT = t_A - t_B   (ps)
//   R      = ln(npe_A / npe_B)
// Applies threshold cut (both channels >= N_PE_THRESH).
// Fits Gaussian to DeltaT and R per position.
// Fits linear calibration: mu_DeltaT vs x_true -> slope m = 2/v_eff.
// Reconstructs x_rec and computes sigma_x per position.
//
// Usage:
//   root -l -b -q 'analysis/analyze_pair_scan.C("results/pairscan_2026-06-11")'
//
// Outputs in results_dir/analysis/:
//   pairscan_summary.csv
//   pairscan_deltaT_vs_x.pdf + .png
//   pairscan_sigma_vs_x.pdf + .png
//   pairscan_R_vs_x.pdf + .png
//   pairscan_efficiency.pdf + .png
//
// Pair and geometry (DetectorConstruction.cc:45-49, DetectorConstruction.hh:37-38):
//   ID_A = 28, x_A = -452 mm  (local_idx=12)
//   ID_B = 29, x_B = -432 mm  (local_idx=13)
//   modal pitch p = 20 mm
//   jitter label: INTRINSIC (no SPTR hook; /sipm/jitterSigma=0 ns)
//   leading-edge threshold: 4 PE equivalent (congruent_sum4_timing.C:218-221 convention)

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <map>
#include <string>
#include <vector>

// ── Configuration ─────────────────────────────────────────────────────────────
static const int    ID_A          = 28;
static const int    ID_B          = 29;
static const double X_A           = -452.0;  // mm
static const double X_B           = -432.0;  // mm
static const double X_MID         = -442.0;  // mm
static const double PITCH         =   20.0;  // mm
static const int    N_PE_THRESH   =    4;    // leading-edge: 4-PE order statistic
static const char*  JITTER_LABEL  = "INTRINSIC";

// ── Helper: list *.root files in directory ────────────────────────────────────
std::vector<std::string> ListRootFiles(const char* dir_path) {
    std::vector<std::string> files;
    DIR* dir = opendir(dir_path);
    if (!dir) {
        fprintf(stderr, "[ERROR] Cannot open directory: %s\n", dir_path);
        return files;
    }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string name(entry->d_name);
        if (name.size() > 5 && name.substr(name.size() - 5) == ".root"
            && name.find("pairscan_x") == 0) {
            files.push_back(std::string(dir_path) + "/" + name);
        }
    }
    closedir(dir);
    std::sort(files.begin(), files.end());
    return files;
}

// ── Helper: extract x_true from filename ─────────────────────────────────────
double ExtractXTrue(const std::string& path) {
    // pattern: pairscan_x<sign><digits>.<digit>mm.root
    auto slash = path.rfind('/');
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    // find 'x' after 'pairscan_'
    auto xpos = base.find("pairscan_x");
    if (xpos == std::string::npos) return 0.0;
    auto after_x = base.substr(xpos + 10);  // skip "pairscan_x"
    auto mm_pos  = after_x.find("mm.root");
    if (mm_pos == std::string::npos) return 0.0;
    return std::stod(after_x.substr(0, mm_pos));
}

// ── Per-position analysis struct ──────────────────────────────────────────────
struct PosResult {
    double x_true;
    double mu_dt;       // Gaussian mean of DeltaT [ps]
    double mu_dt_err;
    double sigma_dt;    // Gaussian sigma of DeltaT [ps]
    double sigma_dt_err;
    double sigma_x_dt;  // position resolution from DeltaT [mm]
    double mu_r;        // Gaussian mean of R = ln(npe_A/npe_B)
    double mu_r_err;
    double sigma_r;
    double sigma_r_err;
    double sigma_x_r;   // position resolution from R [mm]
    double efficiency;  // fraction of events passing both-channel threshold
    double npe_a;       // mean NPE channel A
    double npe_b;       // mean NPE channel B
    int    n_events;
    int    n_passed;
};

// ── Save canvas as PDF + PNG ──────────────────────────────────────────────────
void SaveCanvas(TCanvas* c, const std::string& base) {
    c->SaveAs((base + ".pdf").c_str());
    c->SaveAs((base + ".png").c_str());
}

// ── Main ──────────────────────────────────────────────────────────────────────
void analyze_pair_scan(const char* results_dir = "results/pairscan_2026-06-11") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Create analysis output directory
    std::string ana_dir = std::string(results_dir) + "/analysis";
    gSystem->mkdir(ana_dir.c_str(), /*recursive=*/kTRUE);

    auto files = ListRootFiles(results_dir);
    if (files.empty()) {
        fprintf(stderr, "[ERROR] No pairscan_x*.root files found in %s\n", results_dir);
        return;
    }
    printf("[INFO] Found %zu position files in %s\n", files.size(), results_dir);

    std::vector<PosResult> results;

    for (const auto& fpath : files) {
        const double x_true = ExtractXTrue(fpath);

        TFile* f = TFile::Open(fpath.c_str(), "READ");
        if (!f || f->IsZombie()) {
            fprintf(stderr, "[WARN] Cannot open %s — skipping\n", fpath.c_str());
            continue;
        }

        TTree* tree = dynamic_cast<TTree*>(f->Get("sipm_hits"));
        if (!tree) {
            fprintf(stderr, "[WARN] No sipm_hits TTree in %s\n", fpath.c_str());
            f->Close();
            continue;
        }

        // Read branches
        int    ev_id    = 0;
        int    gid      = 0;
        double time_ns  = 0.0;

        tree->SetBranchAddress("event_id",  &ev_id);
        tree->SetBranchAddress("global_id", &gid);
        tree->SetBranchAddress("time_ns",   &time_ns);

        // Accumulate per-event hits for channels A and B
        std::map<int, std::vector<double>> times_a, times_b;
        Long64_t nentries = tree->GetEntries();

        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            if (gid == ID_A) times_a[ev_id].push_back(time_ns);
            else if (gid == ID_B) times_b[ev_id].push_back(time_ns);
        }

        // Build histograms
        TH1D h_dt ("h_dt",  "h_dt",  120, -3.0, 3.0);   // DeltaT in ns -> convert to ps later
        TH1D h_r  ("h_r",   "h_r",   80,  -4.0, 4.0);

        int n_events = 0, n_passed = 0;
        double sum_npe_a = 0, sum_npe_b = 0;

        // Collect all event IDs that appear in either channel
        std::set<int> all_events;
        for (auto& [ev, _] : times_a) all_events.insert(ev);
        for (auto& [ev, _] : times_b) all_events.insert(ev);
        n_events = (int)all_events.size();

        for (int ev : all_events) {
            auto& ta = times_a[ev];
            auto& tb = times_b[ev];
            sum_npe_a += (double)ta.size();
            sum_npe_b += (double)tb.size();

            if ((int)ta.size() < N_PE_THRESH || (int)tb.size() < N_PE_THRESH)
                continue;

            // Sort and take N_PE_THRESH-th photon (0-based index N_PE_THRESH-1)
            std::sort(ta.begin(), ta.end());
            std::sort(tb.begin(), tb.end());
            double t_a = ta[N_PE_THRESH - 1];
            double t_b = tb[N_PE_THRESH - 1];

            h_dt.Fill(t_a - t_b);  // in ns

            double npe_a = (double)ta.size();
            double npe_b = (double)tb.size();
            if (npe_a > 0 && npe_b > 0)
                h_r.Fill(std::log(npe_a / npe_b));

            n_passed++;
        }

        PosResult res;
        res.x_true    = x_true;
        res.n_events  = n_events;
        res.n_passed  = n_passed;
        res.efficiency = (n_events > 0) ? (double)n_passed / n_events : 0.0;
        res.npe_a     = (n_events > 0) ? sum_npe_a / n_events : 0.0;
        res.npe_b     = (n_events > 0) ? sum_npe_b / n_events : 0.0;

        // Fit DeltaT
        if (h_dt.GetEntries() > 5) {
            TF1 fg("fg_dt", "gaus", h_dt.GetXaxis()->GetXmin(), h_dt.GetXaxis()->GetXmax());
            fg.SetParameters(h_dt.GetMaximum(), h_dt.GetMean(), h_dt.GetRMS());
            h_dt.Fit(&fg, "QNR");
            res.mu_dt       = fg.GetParameter(1) * 1000.0;   // ns -> ps
            res.mu_dt_err   = fg.GetParError(1)  * 1000.0;
            res.sigma_dt    = std::abs(fg.GetParameter(2)) * 1000.0;
            res.sigma_dt_err = fg.GetParError(2) * 1000.0;
        } else {
            res.mu_dt = res.mu_dt_err = res.sigma_dt = res.sigma_dt_err = 0.0;
        }

        // Fit R
        if (h_r.GetEntries() > 5) {
            TF1 fg_r("fg_r", "gaus", h_r.GetXaxis()->GetXmin(), h_r.GetXaxis()->GetXmax());
            fg_r.SetParameters(h_r.GetMaximum(), h_r.GetMean(), h_r.GetRMS());
            h_r.Fit(&fg_r, "QNR");
            res.mu_r       = fg_r.GetParameter(1);
            res.mu_r_err   = fg_r.GetParError(1);
            res.sigma_r    = std::abs(fg_r.GetParameter(2));
            res.sigma_r_err = fg_r.GetParError(2);
        } else {
            res.mu_r = res.mu_r_err = res.sigma_r = res.sigma_r_err = 0.0;
        }

        // sigma_x will be filled after global fit
        res.sigma_x_dt = res.sigma_x_r = 0.0;

        printf("[POS] x=%+.1f mm: %d evts, ε=%.2f, <NPE_A>=%.1f, <NPE_B>=%.1f, "
               "μ_Δt=%.1f ps, σ_Δt=%.1f ps\n",
               x_true, n_events, res.efficiency, res.npe_a, res.npe_b,
               res.mu_dt, res.sigma_dt);

        results.push_back(res);
        f->Close();
    }

    if (results.empty()) {
        fprintf(stderr, "[ERROR] No valid position results.\n");
        return;
    }

    // ── Global linear calibration: mu_DeltaT vs x_true ────────────────────────
    int n = (int)results.size();
    std::vector<double> x_arr(n), mu_dt_arr(n), mu_dt_err_arr(n);
    std::vector<double> mu_r_arr(n), mu_r_err_arr(n);
    for (int i = 0; i < n; i++) {
        x_arr[i]       = results[i].x_true;
        mu_dt_arr[i]   = results[i].mu_dt;
        mu_dt_err_arr[i] = (results[i].mu_dt_err > 0) ? results[i].mu_dt_err : 1.0;
        mu_r_arr[i]    = results[i].mu_r;
        mu_r_err_arr[i] = (results[i].mu_r_err > 0) ? results[i].mu_r_err : 0.001;
    }

    // Fit mu_DeltaT = m*x + b  (slope m in ps/mm)
    TGraphErrors g_dt_cal(n, x_arr.data(), mu_dt_arr.data(), nullptr, mu_dt_err_arr.data());
    TF1 f_lin("f_lin_dt", "pol1", x_arr.front() - 2, x_arr.back() + 2);
    g_dt_cal.Fit(&f_lin, "QNR");
    double slope_dt  = f_lin.GetParameter(1);   // ps/mm
    double slope_err = f_lin.GetParError(1);
    double v_eff     = 0.0;
    if (std::abs(slope_dt) > 1e-6) {
        v_eff = 2000.0 / slope_dt;  // cm/ns: 2 / (slope[ps/mm]) * 10 [mm/cm] / 1000 [ps/ns] * 10
        // v_eff [cm/ns] = 2 [no unit] / (slope [ps/mm]) * (1 ns/1000 ps) * (10 mm/1 cm)
        //              = 2 * 10 / (slope * 1000) * 10 = 2*10*10/(slope*1000)
        // Actually: slope [ps/mm] = dDeltaT/dx, DeltaT = (x_A-x)*(2/v_eff)*1000
        // so v_eff [mm/ns] = 2000/slope_dt [ps/mm], then /10 for cm/ns
        v_eff = (std::abs(slope_dt) > 1e-6) ? 2000.0 / std::abs(slope_dt) / 10.0 : 0.0;
    }

    // Fit mu_R = k*x + b_r (slope k related to 2/lambda_eff)
    TGraphErrors g_r_cal(n, x_arr.data(), mu_r_arr.data(), nullptr, mu_r_err_arr.data());
    TF1 f_lin_r("f_lin_r", "pol1", x_arr.front() - 2, x_arr.back() + 2);
    g_r_cal.Fit(&f_lin_r, "QNR");
    double slope_r   = f_lin_r.GetParameter(1);  // 1/mm
    double lambda_eff = 0.0;
    if (std::abs(slope_r) > 1e-9)
        lambda_eff = 2.0 / std::abs(slope_r) / 10.0;  // cm

    printf("\n[CALIB] Linear calibration:\n");
    printf("  DeltaT slope m = %.3f ± %.3f ps/mm  ->  v_eff = %.2f cm/ns\n",
           slope_dt, slope_err, v_eff);
    printf("  R slope k = %.4f ± %.4f 1/mm  ->  lambda_eff = %.1f cm\n",
           slope_r, f_lin_r.GetParError(1), lambda_eff);

    // ── Compute sigma_x for each position ─────────────────────────────────────
    for (auto& r : results) {
        if (std::abs(slope_dt) > 1e-6)
            r.sigma_x_dt = r.sigma_dt / std::abs(slope_dt);  // mm
        if (std::abs(slope_r) > 1e-9)
            r.sigma_x_r  = r.sigma_r  / std::abs(slope_r);   // mm
    }

    // ── Write CSV ─────────────────────────────────────────────────────────────
    std::string csv_path = ana_dir + "/pairscan_summary.csv";
    FILE* fcsv = fopen(csv_path.c_str(), "w");
    if (fcsv) {
        fprintf(fcsv, "x_true_mm,mu_dt_ps,mu_dt_err_ps,sigma_dt_ps,sigma_dt_err_ps,"
                      "sigma_x_dt_mm,mu_r,mu_r_err,sigma_r,sigma_r_err,"
                      "sigma_x_r_mm,efficiency,npe_a,npe_b,n_events,n_passed\n");
        for (const auto& r : results) {
            fprintf(fcsv,
                    "%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                    "%.4f,%.4f,%.4f,%.4f,%.3f,"
                    "%.4f,%.2f,%.2f,%d,%d\n",
                    r.x_true,
                    r.mu_dt, r.mu_dt_err, r.sigma_dt, r.sigma_dt_err, r.sigma_x_dt,
                    r.mu_r, r.mu_r_err, r.sigma_r, r.sigma_r_err, r.sigma_x_r,
                    r.efficiency, r.npe_a, r.npe_b, r.n_events, r.n_passed);
        }
        fclose(fcsv);
        printf("[OUT] Wrote %s\n", csv_path.c_str());
    }

    // ── Figure 1: mu_DeltaT vs x_true (calibration) ───────────────────────────
    {
        TCanvas c1("c1", "c1", 900, 600);
        c1.SetGrid();

        TGraphErrors* g = new TGraphErrors(n, x_arr.data(), mu_dt_arr.data(),
                                           nullptr, mu_dt_err_arr.data());
        g->SetTitle(TString::Format(
            "Pair scan IDs (%d,%d), p=%.0f mm, Top, %s;"
            "True muon position (mm);"
            "#mu_{#Deltat} (ps)",
            ID_A, ID_B, PITCH, JITTER_LABEL));
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlue + 1);
        g->GetYaxis()->SetTitle("#mu_{#Deltat} = #langle t_{A} - t_{B} #rangle (ps)");
        g->Draw("AP");

        TF1* fit = new TF1("fit_dt", "pol1", x_arr.front() - 2, x_arr.back() + 2);
        fit->SetLineColor(kRed);
        g_dt_cal.GetFunction("f_lin_dt")->Copy(*fit);
        fit->DrawClone("same");

        TLatex lx;
        lx.SetNDC();
        lx.SetTextSize(0.033);
        lx.DrawLatex(0.55, 0.25,
            TString::Format("slope = %.2f #pm %.2f ps/mm", slope_dt, slope_err));
        lx.DrawLatex(0.55, 0.20,
            TString::Format("v_{eff} = %.2f cm/ns", v_eff));

        SaveCanvas(&c1, ana_dir + "/pairscan_deltaT_vs_x");
    }

    // ── Figure 2: sigma_Dt and sigma_x vs x_true ──────────────────────────────
    {
        std::vector<double> sigma_dt_arr(n), sigma_dt_err_arr(n);
        std::vector<double> sigma_x_arr(n);
        for (int i = 0; i < n; i++) {
            sigma_dt_arr[i]  = results[i].sigma_dt;
            sigma_dt_err_arr[i] = results[i].sigma_dt_err;
            sigma_x_arr[i]   = results[i].sigma_x_dt;
        }

        TCanvas c2("c2", "c2", 900, 600);
        c2.Divide(1, 2);

        c2.cd(1);
        TGraphErrors* gs = new TGraphErrors(n, x_arr.data(), sigma_dt_arr.data(),
                                            nullptr, sigma_dt_err_arr.data());
        gs->SetTitle(TString::Format(
            "Pair scan (%d,%d) %s;True muon position (mm);#sigma_{#Deltat} (ps)",
            ID_A, ID_B, JITTER_LABEL));
        gs->SetMarkerStyle(20);
        gs->SetMarkerColor(kBlue + 1);
        gs->Draw("AP");

        c2.cd(2);
        TGraph* gx = new TGraph(n, x_arr.data(), sigma_x_arr.data());
        gx->SetTitle(TString::Format(
            ";True muon position (mm);#sigma_{x}^{#Deltat} (mm)"));
        gx->SetMarkerStyle(21);
        gx->SetMarkerColor(kGreen + 2);
        gx->Draw("AP");

        SaveCanvas(&c2, ana_dir + "/pairscan_sigma_vs_x");
    }

    // ── Figure 3: R = ln(N_A/N_B) vs x_true ──────────────────────────────────
    {
        std::vector<double> sigma_r_arr(n), sigma_r_err_arr(n);
        std::vector<double> sigma_x_r_arr(n);
        for (int i = 0; i < n; i++) {
            sigma_r_arr[i]    = results[i].sigma_r;
            sigma_r_err_arr[i] = results[i].sigma_r_err;
            sigma_x_r_arr[i] = results[i].sigma_x_r;
        }

        TCanvas c3("c3", "c3", 900, 600);
        c3.Divide(1, 2);

        c3.cd(1);
        TGraphErrors* gr = new TGraphErrors(n, x_arr.data(), mu_r_arr.data(),
                                            nullptr, mu_r_err_arr.data());
        gr->SetTitle(TString::Format(
            "Pair scan (%d,%d) %s;True muon position (mm);#mu_{R} = #langle ln(N_{A}/N_{B}) #rangle",
            ID_A, ID_B, JITTER_LABEL));
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kMagenta + 1);
        gr->Draw("AP");

        TF1 f_r_draw("f_r_draw", "pol1", x_arr.front() - 2, x_arr.back() + 2);
        g_r_cal.GetFunction("f_lin_r")->Copy(f_r_draw);
        f_r_draw.SetLineColor(kRed);
        f_r_draw.DrawClone("same");

        TLatex lx;
        lx.SetNDC(); lx.SetTextSize(0.04);
        lx.DrawLatex(0.55, 0.25,
            TString::Format("#lambda_{eff} = %.1f cm", lambda_eff));

        c3.cd(2);
        TGraph* gxr = new TGraph(n, x_arr.data(), sigma_x_r_arr.data());
        gxr->SetTitle(";True muon position (mm);#sigma_{x}^{Q} (mm)");
        gxr->SetMarkerStyle(21);
        gxr->SetMarkerColor(kMagenta + 1);
        gxr->Draw("AP");

        SaveCanvas(&c3, ana_dir + "/pairscan_R_vs_x");
    }

    // ── Figure 4: efficiency vs x_true ────────────────────────────────────────
    {
        std::vector<double> eff_arr(n);
        for (int i = 0; i < n; i++) eff_arr[i] = results[i].efficiency * 100.0;

        TCanvas c4("c4", "c4", 900, 400);
        c4.SetGrid();
        TGraph* ge = new TGraph(n, x_arr.data(), eff_arr.data());
        ge->SetTitle(TString::Format(
            "Pair scan (%d,%d) %s, threshold %d PE;"
            "True muon position (mm);"
            "Efficiency #varepsilon (%%)",
            ID_A, ID_B, JITTER_LABEL, N_PE_THRESH));
        ge->SetMarkerStyle(20);
        ge->SetMarkerColor(kOrange + 1);
        ge->GetYaxis()->SetRangeUser(0, 105);
        ge->Draw("APL");
        SaveCanvas(&c4, ana_dir + "/pairscan_efficiency");
    }

    // ── Summary printout ──────────────────────────────────────────────────────
    printf("\n=== Pair Scan EXEC_08 Summary ===\n");
    printf("  Pair: IDs (%d, %d), x_A=%.0f mm, x_B=%.0f mm, x_mid=%.0f mm\n",
           ID_A, ID_B, X_A, X_B, X_MID);
    printf("  Pitch p = %.0f mm\n", PITCH);
    printf("  Leading-edge threshold: %d PE\n", N_PE_THRESH);
    printf("  Jitter: %s\n", JITTER_LABEL);
    printf("  Calibration: slope_DeltaT = %.3f ps/mm, v_eff = %.2f cm/ns\n",
           slope_dt, v_eff);
    printf("  Calibration: slope_R = %.5f /mm, lambda_eff = %.1f cm\n",
           slope_r, lambda_eff);
    printf("  Positions processed: %d\n", n);

    // Mean sigma_x
    double sum_sx = 0, sum_sxr = 0;
    int cnt = 0;
    for (const auto& r : results) {
        if (r.sigma_x_dt > 0) { sum_sx += r.sigma_x_dt; cnt++; }
        if (r.sigma_x_r > 0)    sum_sxr += r.sigma_x_r;
    }
    if (cnt > 0) {
        printf("  Mean sigma_x(DeltaT) = %.2f mm\n", sum_sx / cnt);
        printf("  Mean sigma_x(R)      = %.2f mm\n", sum_sxr / cnt);
    }
    printf("  CSV: %s\n", csv_path.c_str());
    printf("=================================\n");
}
