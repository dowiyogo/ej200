// tn_order_statistics.C
// ROOT C++ macro: per-channel <t_n> order-statistics for the EXEC_07 scan.
//
// Reads the per-photon sipm_hits TTree (branches: event_id, global_id, time_ns)
// for 7 key beam positions, computes mean and RMS of t_n (n=1..20) for each
// individual channel, and writes the statistics to a CSV file for auditing.
//
// Usage (batch):
//   root -l -b -q 'tn_order_statistics.C("data_dir","output_dir")'
//
// Output: <output_dir>/exec12b_tn_stats.csv
//   Columns: x_mm, ch_id, n, count, mean_ns, rms_ns, mean_npe
//
// High-quality deck figures are produced by exec12b_tn_dispersion.py (Python).
// This macro provides an independent ROOT-native audit of the statistics.
//
// Data source: sipm_hits TTree (analysis/exec07_photon_budget.py:337-342)
//   Branches used: event_id (int), global_id (int), time_ns (double)

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

void tn_order_statistics(
    const char* data_dir  = "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000",
    const char* output_dir = "analysis/exec07")
{
    const int KEY_POSITIONS[] = {-690, -650, -400, 0, 400, 650, 690};
    const int N_POSITIONS = 7;
    const int N_MAX       = 20;
    const int N_EVENTS    = 2000;
    const int MIN_ENTRIES = 10;

    TString csv_path = TString::Format("%s/exec12b_tn_stats.csv", output_dir);
    FILE* fout = fopen(csv_path.Data(), "w");
    if (!fout) {
        fprintf(stderr, "[ERROR] Cannot open %s for writing\n", csv_path.Data());
        return;
    }
    fprintf(fout, "x_mm,ch_id,n,count,mean_ns,rms_ns,mean_npe\n");

    for (int ip = 0; ip < N_POSITIONS; ip++) {
        int x_mm = KEY_POSITIONS[ip];
        TString fname = TString::Format("%s/photon_hits_x%dmm.root", data_dir, x_mm);

        TFile* f = TFile::Open(fname.Data(), "READ");
        if (!f || f->IsZombie()) {
            fprintf(stderr, "[WARN] Cannot open %s -- skipping\n", fname.Data());
            continue;
        }

        TTree* tree = dynamic_cast<TTree*>(f->Get("sipm_hits"));
        if (!tree) {
            fprintf(stderr, "[WARN] No sipm_hits TTree in %s\n", fname.Data());
            f->Close();
            continue;
        }

        // Read all entries into flat vectors
        int    ev_id = 0, ch_id = 0;
        double t_ns  = 0.0;
        tree->SetBranchAddress("event_id",  &ev_id);
        tree->SetBranchAddress("global_id", &ch_id);
        tree->SetBranchAddress("time_ns",   &t_ns);

        // ev_ch_map[event_id][channel_id] -> vector of times
        std::map<int, std::map<int, std::vector<double>>> ev_ch_map;
        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            if (ch_id < 0 || ch_id >= 86) continue;
            ev_ch_map[ev_id][ch_id].push_back(t_ns);
        }

        // Per-channel accumulators
        double sum_tn[86][N_MAX]  = {};
        double sum_tn2[86][N_MAX] = {};
        int    cnt_tn[86][N_MAX]  = {};
        double sum_npe[86]        = {};

        for (auto& [ev, ch_map] : ev_ch_map) {
            for (auto& [ch, times] : ch_map) {
                if (ch < 0 || ch >= 86) continue;
                std::sort(times.begin(), times.end());
                sum_npe[ch] += (double)times.size();
                for (int n = 0; n < N_MAX && n < (int)times.size(); n++) {
                    sum_tn[ch][n]  += times[n];
                    sum_tn2[ch][n] += times[n] * times[n];
                    cnt_tn[ch][n]++;
                }
            }
        }

        // Write CSV rows for channels with enough entries
        for (int ch = 0; ch < 86; ch++) {
            double mean_npe = sum_npe[ch] / N_EVENTS;
            for (int n = 0; n < N_MAX; n++) {
                int cnt = cnt_tn[ch][n];
                if (cnt < MIN_ENTRIES) continue;
                double mean_v = sum_tn[ch][n]  / cnt;
                double mean2  = sum_tn2[ch][n] / cnt;
                double var    = mean2 - mean_v * mean_v;
                double rms    = (var > 0) ? std::sqrt(var) : 0.0;
                fprintf(fout, "%d,%d,%d,%d,%.6f,%.6f,%.4f\n",
                        x_mm, ch, n + 1, cnt, mean_v, rms, mean_npe);
            }
        }

        printf("[OK] x=%+d mm: processed %lld entries\n", x_mm, nentries);
        f->Close();
    }

    fclose(fout);
    printf("Wrote %s\n", csv_path.Data());
}
