#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSystem.h"

namespace {
double dCFDTime(const std::vector<double>& hitTimes, double fraction = 0.14) {
    if (hitTimes.empty()) return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> times = hitTimes;
    std::sort(times.begin(), times.end());

    const double amplitude = static_cast<double>(times.size());
    const double threshold = fraction * amplitude;
    if (threshold <= 1.0) return times.front();

    const int idx = static_cast<int>(std::ceil(threshold)) - 1;
    if (idx <= 0 || idx >= static_cast<int>(times.size())) {
        return times.back();
    }

    const double y0 = static_cast<double>(idx);
    const double y1 = static_cast<double>(idx + 1);
    const double t0 = times[idx - 1];
    const double t1 = times[idx];
    if (std::abs(t1 - t0) < 1e-12) return t1;
    return t0 + (threshold - y0) * (t1 - t0) / (y1 - y0);
}
} // namespace

void analyze_dCFD(const char* firstFile = "photon_hits_run000.root",
                  int faceToAnalyze = 0,
                  double fraction = 0.14,
                  double electronicsSigmaPs = 30.0) {
    TChain *chain = new TChain("sipm_hits");
    if (!gSystem->AccessPathName(firstFile)) {
        chain->Add(firstFile);
    } else {
        for (int i = 0; i <= 20; ++i) {
            TString name = Form("photon_hits_run%03d.root", i);
            if (!gSystem->AccessPathName(name)) {
                chain->Add(name);
            }
        }
    }

    if (chain->GetEntries() == 0) {
        std::cerr << "Error: No se encontraron entradas en sipm_hits." << std::endl;
        return;
    }

    Int_t event_id = 0;
    Int_t face_type = 0;
    Double_t time_ns = 0.0;
    chain->SetBranchAddress("event_id", &event_id);
    chain->SetBranchAddress("face_type", &face_type);
    chain->SetBranchAddress("time_ns", &time_ns);

    Long64_t nEntries = chain->GetEntries();
    std::map<Int_t, std::vector<Double_t>> eventHits;

    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);
        if (face_type == faceToAnalyze) {
            eventHits[event_id].push_back(time_ns);
        }
    }

    TH1D *hTriggerTime = new TH1D(
        "hTriggerTime",
        Form("Tiempo de trigger dCFD %.0f%%;t_{dCFD} [ns];Eventos / bin", fraction * 100.0),
        800, 0, 50);

    TRandom3 rng(12345);
    std::cout << "Procesando " << nEntries << " hits; eventos en face " << faceToAnalyze
              << ": " << eventHits.size() << std::endl;

    for (auto const& [id, times] : eventHits) {
        if (times.empty()) continue;

        const double tCFD = dCFDTime(times, fraction);
        if (std::isnan(tCFD)) continue;
        const double tMeasured = tCFD + rng.Gaus(0.0, electronicsSigmaPs * 1e-3);
        hTriggerTime->Fill(tMeasured);
    }

    TCanvas *c = new TCanvas("c", "Timing Resolution", 800, 600);
    hTriggerTime->Draw();
    hTriggerTime->Fit("gaus", "Q");

    TF1* fit = hTriggerTime->GetFunction("gaus");
    if (!fit) {
        std::cerr << "No se pudo ajustar una Gaussiana al histograma." << std::endl;
        return;
    }

    const Double_t res = fit->GetParameter(2);
    std::cout << "\nResolucion temporal dCFD estimada: " << res * 1000.0 << " ps" << std::endl;
}
