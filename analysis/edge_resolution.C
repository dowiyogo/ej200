// edge_resolution.C — ROOT/CINT companion for edge_resolution.py.
// Usage from repository root:
//   root -l -q 'analysis/edge_resolution.C("photon_hits_run*.root")'

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"

namespace {

struct EventKey {
    int event;
    double x;
    bool operator<(const EventKey& other) const {
        if (event != other.event) return event < other.event;
        return x < other.x;
    }
};

struct EventFace {
    int n = 0;
    double fpt = 1.0e99;
};

double FitSigma(const std::vector<double>& times) {
    if (times.size() < 10) return std::nan("");
    auto sorted = times;
    std::sort(sorted.begin(), sorted.end());
    const double lo0 = sorted[std::max<size_t>(0, sorted.size() / 100)];
    const double hi0 = sorted[std::min(sorted.size() - 1, sorted.size() * 99 / 100)];
    const double margin = std::max(0.5 * (hi0 - lo0), 0.5);

    TH1D h("h_edge_fpt", "h_edge_fpt", 40, lo0 - margin, hi0 + margin);
    for (double t : times) h.Fill(t);
    TF1 f("f_edge_gaus", "gaus", lo0 - margin, hi0 + margin);
    h.Fit(&f, "QNR");
    return std::abs(f.GetParameter(2));
}

} // namespace

void edge_resolution(const char* pattern = "photon_hits_run*.root",
                     const char* outPdf = "edge_resolution.pdf") {
    gStyle->SetOptStat(0);

    TChain chain("sipm_hits");
    chain.Add(pattern);

    int event_id = 0;
    int face_type = 0;
    double time_ns = 0.0;
    double gun_x_mm = 0.0;

    chain.SetBranchAddress("event_id", &event_id);
    chain.SetBranchAddress("face_type", &face_type);
    chain.SetBranchAddress("time_ns", &time_ns);
    chain.SetBranchAddress("gun_x_mm", &gun_x_mm);

    std::map<double, std::map<EventKey, EventFace>> endByX;
    std::map<double, std::map<EventKey, EventFace>> topByX;

    const Long64_t n = chain.GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        EventKey key{static_cast<int>(i / 1000000000LL) + event_id, gun_x_mm};
        auto& store = (face_type == 2) ? topByX[gun_x_mm] : endByX[gun_x_mm];
        auto& ev = store[key];
        ev.n += 1;
        ev.fpt = std::min(ev.fpt, time_ns);
    }

    std::vector<double> xEnd, sEnd, eEnd, xTop, sTop, eTop;
    for (const auto& [x, events] : endByX) {
        std::vector<double> times;
        for (const auto& [key, ev] : events) times.push_back(ev.fpt);
        const double sigma = FitSigma(times);
        if (!std::isnan(sigma)) {
            xEnd.push_back(x);
            sEnd.push_back(1000.0 * sigma);
            eEnd.push_back(1000.0 * sigma / std::sqrt(std::max<size_t>(1, times.size())));
        }
    }
    for (const auto& [x, events] : topByX) {
        std::vector<double> times;
        for (const auto& [key, ev] : events) times.push_back(ev.fpt);
        const double sigma = FitSigma(times);
        if (!std::isnan(sigma)) {
            xTop.push_back(x);
            sTop.push_back(1000.0 * sigma);
            eTop.push_back(1000.0 * sigma / std::sqrt(std::max<size_t>(1, times.size())));
        }
    }

    TCanvas c("c_edge_resolution", "Edge timing resolution", 900, 550);
    TGraphErrors gEnd(xEnd.size(), xEnd.data(), sEnd.data(), nullptr, eEnd.data());
    TGraphErrors gTop(xTop.size(), xTop.data(), sTop.data(), nullptr, eTop.data());
    gEnd.SetTitle("Edge timing resolution;Muon x [mm];#sigma_{t} [ps]");
    gEnd.SetMarkerStyle(20);
    gEnd.SetLineColor(kViolet + 2);
    gEnd.SetMarkerColor(kViolet + 2);
    gTop.SetMarkerStyle(21);
    gTop.SetLineColor(kGreen + 2);
    gTop.SetMarkerColor(kGreen + 2);
    gEnd.Draw("AP");
    gTop.Draw("P SAME");

    TLegend leg(0.62, 0.72, 0.88, 0.88);
    leg.AddEntry(&gEnd, "End SiPMs", "lp");
    leg.AddEntry(&gTop, "Top SiPMs", "lp");
    leg.Draw();

    c.SaveAs(outPdf);
}
