// grouped_resolution.C — ROOT/CINT Sum-of-N timing emulation.
// Usage:
//   root -l -q 'analysis/grouped_resolution.C("photon_hits_run*.root",4)'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

namespace {

using Group = std::vector<int>;
using Grouping = std::vector<Group>;
using EventKey = std::pair<int, int>; // rounded x, event_id

struct FitResult {
    double mu = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigmaErr = std::numeric_limits<double>::quiet_NaN();
    int n = 0;
};

std::map<std::string, Grouping> BuildGroupings() {
    return {
        {"single_end_left",  {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}}},
        {"sum4_end",         {{0,1,2,3}, {4,5,6,7}, {8,9,10,11}, {12,13,14,15}}},
        {"sum8_end",         {{0,1,2,3,4,5,6,7}, {8,9,10,11,12,13,14,15}}},
        {"full_end",         {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}}},
        {"single_top",       {{16}, {17}, {18}, {19}, {20}, {21}, {22}, {23}, {24}, {25},
                              {26}, {27}, {28}, {29}, {30}, {31}, {32}, {33}, {34}, {35}}},
        {"sum4_top",         {{16,17,18,19}, {20,21,22,23}, {24,25,26,27},
                              {28,29,30,31}, {32,33,34,35}}},
        {"full_top",         {{16,17,18,19,20,21,22,23,24,25,
                              26,27,28,29,30,31,32,33,34,35}}},
    };
}

double KthPhotonTime(std::vector<double> times, int k) {
    if (static_cast<int>(times.size()) < k) return std::numeric_limits<double>::quiet_NaN();
    std::nth_element(times.begin(), times.begin() + (k - 1), times.end());
    return times[k - 1];
}

double Mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double Std(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double m = Mean(v);
    double s2 = 0.0;
    for (double x : v) s2 += (x - m) * (x - m);
    return std::sqrt(s2 / static_cast<double>(v.size() - 1));
}

double Percentile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double pos = q * (v.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    const double frac = pos - lo;
    if (lo == hi) return v[lo];
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

FitResult FitFpt(const std::vector<double>& times, const std::string& hname) {
    FitResult out;
    out.n = static_cast<int>(times.size());
    if (times.size() < 10) return out;

    double lo = Percentile(times, 0.01);
    double hi = Percentile(times, 0.99);
    const double margin = std::max(0.5 * (hi - lo), 0.5);
    lo -= margin;
    hi += margin;

    TH1D h(hname.c_str(), "", 40, lo, hi);
    for (double t : times) h.Fill(t);

    TF1 f((hname + "_gaus").c_str(), "gaus", lo, hi);
    f.SetParameters(h.GetMaximum(), Mean(times), std::max(Std(times), 0.1));
    f.SetParLimits(2, 1e-4, hi - lo);
    const int status = h.Fit(&f, "QNR");
    if (status == 0) {
        out.mu = f.GetParameter(1);
        out.sigma = std::abs(f.GetParameter(2));
        out.sigmaErr = f.GetParError(2);
    } else {
        out.mu = Mean(times);
        out.sigma = Std(times);
        out.sigmaErr = out.sigma / std::sqrt(static_cast<double>(times.size()));
    }
    return out;
}

bool Contains(const Group& g, int id) {
    return std::find(g.begin(), g.end(), id) != g.end();
}

} // namespace

void grouped_resolution(const char* pattern = "photon_hits_merged.root",
                        int threshold = 4,
                        const char* outPdf = "grouped_resolution_root.pdf") {
    gStyle->SetOptStat(0);
    const auto groupings = BuildGroupings();

    TChain chain("sipm_hits");
    chain.Add(pattern);
    TTreeReader reader(&chain);
    TTreeReaderValue<int> event_id(reader, "event_id");
    TTreeReaderValue<int> global_id(reader, "global_id");
    TTreeReaderValue<double> time_ns(reader, "time_ns");
    TTreeReaderValue<double> gun_x_mm(reader, "gun_x_mm");

    std::map<std::string, std::map<EventKey, std::map<int, std::vector<double>>>> byGrouping;

    while (reader.Next()) {
        const int x = static_cast<int>(std::lround(*gun_x_mm));
        const EventKey key{x, *event_id};
        for (const auto& [name, groups] : groupings) {
            for (size_t ig = 0; ig < groups.size(); ++ig) {
                if (Contains(groups[ig], *global_id)) {
                    byGrouping[name][key][static_cast<int>(ig)].push_back(*time_ns);
                }
            }
        }
    }

    std::map<std::string, std::map<int, std::vector<double>>> triggersByX;
    for (const auto& [name, events] : byGrouping) {
        for (const auto& [key, perGroup] : events) {
            double earliest = std::numeric_limits<double>::infinity();
            for (const auto& [ig, times] : perGroup) {
                const double kth = KthPhotonTime(times, threshold);
                if (std::isfinite(kth)) earliest = std::min(earliest, kth);
            }
            if (std::isfinite(earliest)) triggersByX[name][key.first].push_back(earliest);
        }
    }

    TCanvas c("c_grouped_resolution", "Grouped resolution", 1100, 600);
    c.SetGrid();
    TMultiGraph mg;
    TLegend leg(0.62, 0.58, 0.88, 0.88);

    int color = 1;
    std::ofstream csv("grouped_resolution_root.csv");
    csv << "grouping,x_mm,mu_ns,sigma_ns,sigma_err,n_events\n";
    std::vector<TGraphErrors*> graphs;

    for (const auto& [name, byX] : triggersByX) {
        auto* g = new TGraphErrors();
        int ip = 0;
        for (const auto& [x, times] : byX) {
            const auto fit = FitFpt(times, "h_" + name + "_" + std::to_string(x));
            if (!std::isfinite(fit.sigma) || fit.sigma <= 0) continue;
            g->SetPoint(ip, x, fit.sigma * 1e3);
            g->SetPointError(ip, 0.0, fit.sigmaErr * 1e3);
            csv << name << "," << x << "," << fit.mu << "," << fit.sigma << ","
                << fit.sigmaErr << "," << fit.n << "\n";
            ++ip;
        }
        if (ip == 0) {
            delete g;
            continue;
        }
        g->SetName(("g_" + name).c_str());
        g->SetLineColor(color);
        g->SetMarkerColor(color);
        g->SetMarkerStyle(20 + (color % 10));
        g->SetLineWidth(2);
        mg.Add(g, "LP");
        leg.AddEntry(g, name.c_str(), "lp");
        graphs.push_back(g);
        ++color;
        if (color == 10) color = 30;
    }
    csv.close();

    mg.SetTitle(Form("Grouped SiPM timing resolution (threshold=%d p.e.);Muon x [mm];#sigma_{t} [ps]",
                     threshold));
    mg.Draw("A");
    leg.Draw();
    c.SaveAs(outPdf);
}
