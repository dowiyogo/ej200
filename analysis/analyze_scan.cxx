// analyze_scan.cxx
// Longitudinal scan analysis: Npe vs x  and  v_eff from Δt slope.
//
// Each position file is processed in a separate thread (one TFile per thread,
// which is safe with ROOT::EnableImplicitMT). Results are merged on the main
// thread for fitting with TGraph + TF1.
//
// Usage:
//   ./analyze_scan <run_dir> [OPTIONS]
//
// Options:
//   --threads  N   number of worker threads  (default: nCores)
//   --nevents  N   events simulated per position  (default: auto from tree)
//   --out      F   output ROOT file  (default: <run_dir>/scan_analysis.root)
//   --tsv      F   output TSV file   (default: <run_dir>/scan_analysis.tsv)

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <ROOT/RDataFrame.hxx>   // for ROOT::EnableImplicitMT

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <future>
#include <thread>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <set>

namespace fs = std::filesystem;

// ── Physics ────────────────────────────────────────────────────────────────
constexpr double kNbar      = 1.58;
constexpr double kCLight    = 29.9792458;       // cm/ns
constexpr double kVgroup    = kCLight / kNbar;  // 18.97 cm/ns
constexpr double kVeffLow   = 14.0;             // cm/ns physical lower bound
constexpr double kNpeEndMin = 5.0;              // PE/event minimum at END faces

// ── Per-position result ─────────────────────────────────────────────────────
struct PosResult {
    double x_mm      = 0.0;
    int    n_sim     = 0;    // events simulated (may be estimated from unique ids)
    int    n_hits    = 0;    // total hits in tree

    // Npe
    long   hits_L    = 0;
    long   hits_R    = 0;
    long   hits_T    = 0;
    double npe_L     = 0.0;
    double npe_R     = 0.0;
    double npe_T     = 0.0;

    // Timing
    int    n_paired  = 0;
    double mean_dt   = 0.0;  // mean(t_R_first - t_L_first) [ns]
    double rms_dt    = 0.0;  // sample std dev of dt [ns]
    double err_dt    = 0.0;  // standard error of mean [ns]

    std::string path;
    bool   ok = false;
};

// ── Single-file analysis (runs in worker thread) ────────────────────────────
PosResult analyzeFile(const std::string& path, int n_events_override) {
    PosResult res;
    res.path = path;

    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[FAIL] cannot open: " << path << "\n";
        return res;
    }

    TTree* tree = nullptr;
    f->GetObject("sipm_hits", tree);
    if (!tree || tree->GetEntries() == 0) {
        // empty / missing tree is not necessarily an error (0-hit run)
        res.ok = true;
        f->Close();
        return res;
    }

    // Read only the branches we need
    tree->SetBranchStatus("*",         0);
    tree->SetBranchStatus("event_id",  1);
    tree->SetBranchStatus("face_type", 1);
    tree->SetBranchStatus("time_ns",   1);
    tree->SetBranchStatus("gun_x_mm",  1);

    Int_t     eid;
    Int_t     face;
    Double_t  t, gx;
    tree->SetBranchAddress("event_id",  &eid);
    tree->SetBranchAddress("face_type", &face);
    tree->SetBranchAddress("time_ns",   &t);
    tree->SetBranchAddress("gun_x_mm",  &gx);

    // Accumulators (per-event first arrival)
    std::unordered_map<Int_t, double> tL, tR;   // event_id → t_first
    tL.reserve(8192);
    tR.reserve(8192);
    long hits_L = 0, hits_R = 0, hits_T = 0;
    double x_pos = -9999.0;
    std::set<Int_t> seen;

    const Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        seen.insert(eid);
        if (x_pos < -9000.0) x_pos = gx;

        switch (face) {
          case 0:  // END left
            ++hits_L;
            {   auto [it, ins] = tL.emplace(eid, t);
                if (!ins && t < it->second) it->second = t;   }
            break;
          case 1:  // END right
            ++hits_R;
            {   auto [it, ins] = tR.emplace(eid, t);
                if (!ins && t < it->second) it->second = t;   }
            break;
          default: // TOP or other
            ++hits_T;
            break;
        }
    }

    // Number of events: prefer caller override, else unique ids in tree
    int n_ev = (n_events_override > 0) ? n_events_override : (int)seen.size();

    // Δt: iterate over events with both END faces hit
    double sum1 = 0.0, sum2 = 0.0;
    int np = 0;
    for (auto& [id, t_left] : tL) {   // NOLINT(build/include_what_you_use)
        auto it = tR.find(id);
        if (it != tR.end()) {
            double dt = it->second - t_left;
            sum1 += dt;
            sum2 += dt * dt;
            ++np;
        }
    }

    res.x_mm     = x_pos;
    res.n_sim    = n_ev;
    res.n_hits   = (int)N;
    res.hits_L   = hits_L;
    res.hits_R   = hits_R;
    res.hits_T   = hits_T;
    res.npe_L    = n_ev > 0 ? (double)hits_L / n_ev : 0.0;
    res.npe_R    = n_ev > 0 ? (double)hits_R / n_ev : 0.0;
    res.npe_T    = n_ev > 0 ? (double)hits_T / n_ev : 0.0;
    res.n_paired = np;
    if (np > 0) {
        res.mean_dt = sum1 / np;
        res.rms_dt  = (np > 1)
            ? std::sqrt((sum2 - sum1*sum1/np) / (np - 1))
            : 0.0;
        res.err_dt  = (np > 0) ? res.rms_dt / std::sqrt(np) : 0.0;
    }
    res.ok = true;

    f->Close();
    return res;
}

// ── Thread-pool via std::async ──────────────────────────────────────────────
std::vector<PosResult> processInParallel(
    const std::vector<std::string>& files,
    int nThreads,
    int nEventsOverride)
{
    std::vector<PosResult> results;
    results.reserve(files.size());

    // Submit in batches of nThreads
    size_t i = 0;
    while (i < files.size()) {
        std::vector<std::future<PosResult>> batch;
        for (int t = 0; t < nThreads && i < files.size(); ++t, ++i)
            batch.push_back(std::async(std::launch::async,
                                       analyzeFile, files[i], nEventsOverride));
        for (auto& fut : batch)
            results.push_back(fut.get());
    }

    std::sort(results.begin(), results.end(),
              [](const PosResult& a, const PosResult& b){ return a.x_mm < b.x_mm; });
    return results;
}

// ── Helpers ─────────────────────────────────────────────────────────────────
static std::string check(bool ok) { return ok ? "  [OK]  " : "  [FAIL]"; }

static void printTable(const std::vector<PosResult>& res) {
    std::cout << "\n";
    std::cout << std::left  << std::setw(10) << "x [mm]"
              << std::right << std::setw(8)  << "N_ev"
              << std::setw(10) << "Npe_L"
              << std::setw(10) << "Npe_R"
              << std::setw(10) << "Npe_TOP"
              << std::setw(10) << "N_pair"
              << std::setw(12) << "<Δt> [ns]"
              << std::setw(10) << "σ(Δt)"
              << "\n";
    std::cout << std::string(80, '-') << "\n";
    for (const auto& r : res) {
        std::cout << std::fixed << std::setprecision(1)
                  << std::left  << std::setw(10) << r.x_mm
                  << std::right << std::setw(8)  << r.n_sim
                  << std::setw(10) << r.npe_L
                  << std::setw(10) << r.npe_R
                  << std::setw(10) << r.npe_T
                  << std::setw(10) << r.n_paired
                  << std::setw(12) << r.mean_dt
                  << std::setw(10) << r.rms_dt
                  << "\n";
    }
}

// ── Main ────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <run_dir> [--threads N] [--nevents N] [--out F] [--tsv F]\n";
        return 1;
    }

    std::string run_dir       = argv[1];
    int         nThreads      = (int)std::thread::hardware_concurrency();
    int         nEventsOvr    = 0;           // 0 = auto
    std::string out_root;
    std::string out_tsv;

    for (int i = 2; i < argc - 1; ++i) {
        std::string a = argv[i];
        if      (a == "--threads")  nThreads   = std::stoi(argv[++i]);
        else if (a == "--nevents")  nEventsOvr = std::stoi(argv[++i]);
        else if (a == "--out")      out_root   = argv[++i];
        else if (a == "--tsv")      out_tsv    = argv[++i];
    }
    if (out_root.empty()) out_root = run_dir + "/scan_analysis.root";
    if (out_tsv.empty())  out_tsv  = run_dir + "/scan_analysis.tsv";
    nThreads = std::max(1, nThreads);

    // ROOT multi-thread mode (makes TFile::Open thread-safe)
    ROOT::EnableImplicitMT(nThreads);

    // Discover ROOT files
    std::vector<std::string> files;
    for (auto& e : fs::recursive_directory_iterator(run_dir)) {
        if (e.path().extension() == ".root" &&
            e.path().filename().string().find("photon_hits") != std::string::npos)
            files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());

    if (files.empty()) {
        std::cerr << "No ROOT files found under " << run_dir << "\n";
        return 1;
    }

    std::cout << "══════════════════════════════════════════════════════════\n"
              << "  EJ-204 Scan Analysis (C++ / ROOT " << gROOT->GetVersion() << ")\n"
              << "  Run dir : " << run_dir  << "\n"
              << "  Files   : " << files.size() << " positions\n"
              << "  Threads : " << nThreads << "\n"
              << "══════════════════════════════════════════════════════════\n";

    auto t0 = std::chrono::steady_clock::now();
    auto results = processInParallel(files, nThreads, nEventsOvr);
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    long total_hits = 0;
    for (auto& r : results) total_hits += r.n_hits;
    std::cout << "  Processed " << total_hits << " hits in "
              << std::fixed << std::setprecision(1) << elapsed << " s\n";

    printTable(results);

    // ── v_eff fit ───────────────────────────────────────────────────────────
    std::cout << "\n── v_eff from Δt slope ────────────────────────────────────\n";

    std::vector<double> xs, dts, errs;
    for (const auto& r : results) {
        if (r.n_paired < 10) continue;   // skip positions with too few pairs
        xs.push_back(r.x_mm);
        dts.push_back(r.mean_dt);
        errs.push_back(r.err_dt > 0 ? r.err_dt : 0.001);
    }

    double v_eff = -1.0, slope = 0.0, chi2ndf = -1.0;
    TGraphErrors* gDt = nullptr;
    TF1*          fPol1 = nullptr;

    if ((int)xs.size() >= 4) {
        gDt = new TGraphErrors((int)xs.size(),
                               xs.data(), dts.data(), nullptr, errs.data());
        gDt->SetName("g_dt_vs_x");
        gDt->SetTitle(";Gun X [mm]; #Deltat = t_{R,first} - t_{L,first} [ns]");
        gDt->SetMarkerStyle(20);

        fPol1 = new TF1("f_pol1", "pol1", xs.front() - 20, xs.back() + 20);
        auto fitRes = gDt->Fit(fPol1, "QS");   // quiet, return status
        slope   = fPol1->GetParameter(1);       // ns/mm
        chi2ndf = (fitRes->Ndf() > 0)
                  ? fitRes->Chi2() / fitRes->Ndf() : -1.0;

        if (std::abs(slope) > 1e-10) {
            double v_mm_ns = -2.0 / slope;
            v_eff = v_mm_ns / 10.0;             // cm/ns
        }

        std::cout << "  Slope  a   = " << std::scientific << std::setprecision(5)
                  << slope << " ns/mm\n"
                  << "  Intercept  = " << fPol1->GetParameter(0) << " ns\n"
                  << std::fixed << std::setprecision(3)
                  << "  v_eff      = " << v_eff << " cm/ns"
                  << "   [expected " << kVeffLow << " – " << kVgroup << " cm/ns]\n"
                  << "  χ²/ndf     = " << std::setprecision(2) << chi2ndf << "\n"
                  << "  N positions used in fit: " << xs.size() << "\n";
    } else {
        std::cout << "  Too few paired positions for fit (" << xs.size() << ")\n";
    }

    // ── Mean Npe across positions ───────────────────────────────────────────
    double sum_L = 0, sum_R = 0, sum_T = 0;
    int n_valid = 0;
    for (auto& r : results) {
        if (!r.ok) continue;
        sum_L += r.npe_L;
        sum_R += r.npe_R;
        sum_T += r.npe_T;
        ++n_valid;
    }
    double mean_L = n_valid > 0 ? sum_L / n_valid : 0.0;
    double mean_R = n_valid > 0 ? sum_R / n_valid : 0.0;
    double mean_T = n_valid > 0 ? sum_T / n_valid : 0.0;

    // ── Verdict ─────────────────────────────────────────────────────────────
    bool npe_ok  = mean_L >= kNpeEndMin && mean_R >= kNpeEndMin;
    bool veff_ok = v_eff > 0 && v_eff >= kVeffLow && v_eff <= kVgroup * 1.15;
    bool ok_all  = npe_ok && veff_ok;

    std::cout << "\n── Verdict ─────────────────────────────────────────────────\n"
              << std::fixed << std::setprecision(1)
              << "  Npe END-left  : " << mean_L << " PE/ev" << check(mean_L >= kNpeEndMin) << "\n"
              << "  Npe END-right : " << mean_R << " PE/ev" << check(mean_R >= kNpeEndMin) << "\n"
              << "  Npe TOP       : " << mean_T << " PE/ev\n"
              << std::setprecision(2)
              << "  v_eff         : " << v_eff  << " cm/ns" << check(veff_ok) << "\n"
              << (ok_all ? "\n  ✓  SCAN VALID — physics confirmed\n"
                         : "\n  ✗  CHECK RESULTS\n")
              << "══════════════════════════════════════════════════════════\n";

    // ── Write outputs ───────────────────────────────────────────────────────
    // TSV
    {
        std::ofstream tsv(out_tsv);
        tsv << "x_mm\tn_events\tnpe_left\tnpe_right\tnpe_top"
               "\tn_paired\tmean_dt_ns\trms_dt_ns\n";
        for (const auto& r : results) {
            tsv << std::fixed << std::setprecision(3)
                << r.x_mm    << "\t" << r.n_sim  << "\t"
                << r.npe_L   << "\t" << r.npe_R  << "\t" << r.npe_T  << "\t"
                << r.n_paired<< "\t" << r.mean_dt << "\t" << r.rms_dt << "\n";
        }
        std::cout << "  TSV  → " << out_tsv << "\n";
    }

    // ROOT
    {
        TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

        // Npe graphs
        std::vector<double> vx, vnL, vnR, vnT;
        for (auto& r : results) {
            if (!r.ok) continue;
            vx.push_back(r.x_mm);
            vnL.push_back(r.npe_L);
            vnR.push_back(r.npe_R);
            vnT.push_back(r.npe_T);
        }
        auto* gL = new TGraph((int)vx.size(), vx.data(), vnL.data());
        gL->SetName("g_npe_left");
        gL->SetTitle(";Gun X [mm]; Npe END-left [PE/ev]");

        auto* gR = new TGraph((int)vx.size(), vx.data(), vnR.data());
        gR->SetName("g_npe_right");
        gR->SetTitle(";Gun X [mm]; Npe END-right [PE/ev]");

        auto* gT = new TGraph((int)vx.size(), vx.data(), vnT.data());
        gT->SetName("g_npe_top");
        gT->SetTitle(";Gun X [mm]; Npe TOP [PE/ev]");

        gL->Write(); gR->Write(); gT->Write();
        if (gDt)   { gDt->Write();  }
        if (fPol1) { fPol1->Write(); }
        fout->Close();
        std::cout << "  ROOT → " << out_root << "\n";
    }

    return ok_all ? 0 : 1;
}
