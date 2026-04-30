#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"

namespace {

enum TimeMarkMode {
    kFirstPhoton = 0,
    kNthPhoton = 1,
    kMeanFirstN = 2
};

struct SensorStats {
    int face_type = -1;
    int global_id = -1;
    int local_id = -1;
    int n_events = 0;
    double mean_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_ns = std::numeric_limits<double>::quiet_NaN();
};

const char* FaceName(int face) {
    if (face == 0) return "end-left";
    if (face == 1) return "end-right";
    if (face == 2) return "top";
    return "unknown";
}

const char* ModeName(int mode) {
    if (mode == kFirstPhoton) return "first-photon";
    if (mode == kNthPhoton) return "nth-photon";
    if (mode == kMeanFirstN) return "mean-first-n";
    return "unknown";
}

Double_t ComputeTimeMark(const std::vector<Double_t>& times, Int_t mode, size_t nParam) {
    if (times.empty()) return -1.0;

    std::vector<Double_t> sorted = times;
    std::sort(sorted.begin(), sorted.end());

    if (mode == kFirstPhoton) return sorted.front();

    if (mode == kNthPhoton) {
        if (nParam == 0 || sorted.size() < nParam) return -1.0;
        return sorted[nParam - 1];
    }

    if (mode == kMeanFirstN) {
        const size_t nUse = std::min(sorted.size(), nParam);
        if (nUse == 0) return -1.0;
        Double_t sum = 0.0;
        for (size_t i = 0; i < nUse; ++i) sum += sorted[i];
        return sum / static_cast<Double_t>(nUse);
    }

    return -1.0;
}

bool AddScanFiles(TChain& chain) {
    bool filesFound = false;
    const std::vector<TString> searchDirs = {"", "build/", "../build/"};

    for (int i = 0; i <= 20; ++i) {
        bool addedThisRun = false;
        for (const auto& dir : searchDirs) {
            TString name = Form("%sphoton_hits_run%03d.root", dir.Data(), i);
            if (!gSystem->AccessPathName(name)) {
                chain.Add(name);
                filesFound = true;
                addedThisRun = true;
                break;
            }
        }
        if (!addedThisRun) {
            std::cerr << "Aviso: no se encontro photon_hits_run"
                      << Form("%03d", i)
                      << ".root en ., build/ o ../build/" << std::endl;
        }
    }
    return filesFound;
}

double Mean(const std::vector<double>& v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum / static_cast<double>(v.size());
}

double StdDev(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double mu = Mean(v);
    double s2 = 0.0;
    for (double x : v) {
        const double d = x - mu;
        s2 += d * d;
    }
    return std::sqrt(s2 / static_cast<double>(v.size() - 1));
}

}  // namespace

void SiPMRankingScan_RMS(Int_t mode = kFirstPhoton, Int_t nParam = 1, Int_t minEvents = 50) {
    TChain chain("sipm_hits");
    if (!AddScanFiles(chain)) {
        std::cerr << "Error: no se encontraron archivos del scan." << std::endl;
        return;
    }

    std::cout << "SiPMRankingScan: modo = " << ModeName(mode)
              << ", nParam = " << nParam
              << ", minEvents = " << minEvents << std::endl;
    std::cout << "SiPMRankingScan: entradas totales en el chain = "
              << chain.GetEntries() << std::endl;

    Int_t event_id = -1;
    Int_t face_type = -1;
    Int_t global_id = -1;
    Int_t local_id = -1;
    Double_t time_ns = 0.0;
    Double_t gun_x_mm = 0.0;

    chain.SetBranchAddress("event_id", &event_id);
    chain.SetBranchAddress("face_type", &face_type);
    chain.SetBranchAddress("global_id", &global_id);
    chain.SetBranchAddress("local_id", &local_id);
    chain.SetBranchAddress("time_ns", &time_ns);
    chain.SetBranchAddress("gun_x_mm", &gun_x_mm);

    using EventHits = std::map<int, std::vector<Double_t>>;
    using XMap = std::map<int, EventHits>;
    XMap hitsBySensorAndX;
    std::map<int, std::pair<int, int>> sensorMeta;
    std::vector<int> xPositions;

    for (Long64_t i = 0; i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);
        const int x = static_cast<int>(std::lround(gun_x_mm));
        const int sensorKey = global_id;
        hitsBySensorAndX[sensorKey][x * 1000000 + event_id].push_back(time_ns);
        sensorMeta[sensorKey] = {face_type, local_id};
    }

    for (const auto& kv : hitsBySensorAndX) {
        for (const auto& xv : kv.second) {
            const int x = xv.first / 1000000;
            if (std::find(xPositions.begin(), xPositions.end(), x) == xPositions.end()) {
                xPositions.push_back(x);
            }
        }
    }
    std::sort(xPositions.begin(), xPositions.end());

    std::vector<int> sensors;
    for (const auto& kv : sensorMeta) sensors.push_back(kv.first);
    std::sort(sensors.begin(), sensors.end());

    if (xPositions.empty() || sensors.empty()) {
        std::cerr << "Error: no se pudieron reconstruir posiciones o sensores." << std::endl;
        return;
    }

    TH2D* hMean = new TH2D("hMeanTimeMap",
                           Form("Tiempo medio por SiPM (%s, n=%d); X del muon [mm]; global_id del SiPM",
                                ModeName(mode), nParam),
                           static_cast<int>(xPositions.size()),
                           xPositions.front() - 32.5,
                           xPositions.back() + 32.5,
                           static_cast<int>(sensors.size()),
                           sensors.front() - 0.5,
                           sensors.back() + 0.5);
    TH2D* hSigma = new TH2D("hSigmaTimeMap",
                            Form("Resolucion temporal por SiPM (%s, n=%d); X del muon [mm]; global_id del SiPM",
                                 ModeName(mode), nParam),
                            static_cast<int>(xPositions.size()),
                            xPositions.front() - 32.5,
                            xPositions.back() + 32.5,
                            static_cast<int>(sensors.size()),
                            sensors.front() - 0.5,
                            sensors.back() + 0.5);
    TH2D* hCounts = new TH2D("hCountMap",
                             "Eventos utiles por SiPM; X del muon [mm]; global_id del SiPM",
                             static_cast<int>(xPositions.size()),
                             xPositions.front() - 32.5,
                             xPositions.back() + 32.5,
                             static_cast<int>(sensors.size()),
                             sensors.front() - 0.5,
                             sensors.back() + 0.5);

    std::map<int, std::vector<SensorStats>> rankingByX;
    std::map<int, std::vector<double>> bestSigmaByFaceX;

    for (int sensor : sensors) {
        const int face = sensorMeta[sensor].first;
        const int local = sensorMeta[sensor].second;

        std::map<int, std::vector<double>> marksByX;
        for (const auto& ev : hitsBySensorAndX[sensor]) {
            const int x = ev.first / 1000000;
            const double tMark = ComputeTimeMark(ev.second, mode, static_cast<size_t>(nParam));
            if (tMark >= 0.0) marksByX[x].push_back(tMark);
        }

        for (const auto& xv : marksByX) {
            const int x = xv.first;
            const auto& marks = xv.second;
            const int n = static_cast<int>(marks.size());
            const double mu = Mean(marks);
            const double sigma = StdDev(marks);

            hCounts->Fill(x, sensor, n);
            if (n >= minEvents) {
                hMean->Fill(x, sensor, mu);
                hSigma->Fill(x, sensor, sigma * 1000.0);

                rankingByX[x].push_back(SensorStats{face, sensor, local, n, mu, sigma});
                bestSigmaByFaceX[face * 100000 + x].push_back(sigma * 1000.0);
            }
        }
    }

    std::cout << "\nRanking por posicion X\n";
    for (int x : xPositions) {
        auto& stats = rankingByX[x];
        if (stats.empty()) {
            std::cout << "  x = " << x << " mm: sin sensores con n >= "
                      << minEvents << std::endl;
            continue;
        }

        auto byMean = stats;
        std::sort(byMean.begin(), byMean.end(),
                  [](const SensorStats& a, const SensorStats& b) {
                      return a.mean_ns < b.mean_ns;
                  });

        auto bySigma = stats;
        std::sort(bySigma.begin(), bySigma.end(),
                  [](const SensorStats& a, const SensorStats& b) {
                      return a.sigma_ns < b.sigma_ns;
                  });

        std::cout << "  x = " << x << " mm" << std::endl;
        std::cout << "    Earliest 3 sensors:" << std::endl;
        for (size_t i = 0; i < std::min<size_t>(3, byMean.size()); ++i) {
            const auto& s = byMean[i];
            std::cout << "      " << i + 1
                      << ". gid=" << s.global_id
                      << " (" << FaceName(s.face_type) << ", local=" << s.local_id << ")"
                      << ", mean=" << s.mean_ns << " ns"
                      << ", sigma=" << s.sigma_ns * 1000.0 << " ps"
                      << ", n=" << s.n_events << std::endl;
        }
        std::cout << "    Best 3 sigma sensors:" << std::endl;
        for (size_t i = 0; i < std::min<size_t>(3, bySigma.size()); ++i) {
            const auto& s = bySigma[i];
            std::cout << "      " << i + 1
                      << ". gid=" << s.global_id
                      << " (" << FaceName(s.face_type) << ", local=" << s.local_id << ")"
                      << ", sigma=" << s.sigma_ns * 1000.0 << " ps"
                      << ", mean=" << s.mean_ns << " ns"
                      << ", n=" << s.n_events << std::endl;
        }
    }

    std::vector<double> xLeft, yLeft, xRight, yRight, xTop, yTop;
    for (int x : xPositions) {
        double bestLeft = std::numeric_limits<double>::infinity();
        double bestRight = std::numeric_limits<double>::infinity();
        double bestTop = std::numeric_limits<double>::infinity();

        for (const auto& s : rankingByX[x]) {
            const double sigmaPs = s.sigma_ns * 1000.0;
            if (s.face_type == 0) bestLeft = std::min(bestLeft, sigmaPs);
            if (s.face_type == 1) bestRight = std::min(bestRight, sigmaPs);
            if (s.face_type == 2) bestTop = std::min(bestTop, sigmaPs);
        }

        if (std::isfinite(bestLeft)) {
            xLeft.push_back(x);
            yLeft.push_back(bestLeft);
        }
        if (std::isfinite(bestRight)) {
            xRight.push_back(x);
            yRight.push_back(bestRight);
        }
        if (std::isfinite(bestTop)) {
            xTop.push_back(x);
            yTop.push_back(bestTop);
        }
    }

    TCanvas* cMaps = new TCanvas("cSiPMRankingMaps", "SiPM ranking maps", 1800, 600);
    cMaps->Divide(3, 1);
    cMaps->cd(1);
    hMean->Draw("COLZ TEXT");
    cMaps->cd(2);
    hSigma->Draw("COLZ TEXT");
    cMaps->cd(3);
    hCounts->Draw("COLZ TEXT");

    TCanvas* cBest = new TCanvas("cBestFaceSigma", "Best sigma by face", 1000, 550);
    TGraph* gLeft = new TGraph(static_cast<int>(xLeft.size()), xLeft.data(), yLeft.data());
    TGraph* gRight = new TGraph(static_cast<int>(xRight.size()), xRight.data(), yRight.data());
    TGraph* gTop = new TGraph(static_cast<int>(xTop.size()), xTop.data(), yTop.data());

    gLeft->SetTitle("Mejor resolucion temporal por cara; X del muon [mm]; mejor #sigma_{t} por cara [ps]");
    gLeft->SetLineColor(kBlue + 1);
    gLeft->SetMarkerColor(kBlue + 1);
    gLeft->SetMarkerStyle(20);
    gLeft->SetLineWidth(2);
    gLeft->SetMinimum(0.0);
    gLeft->Draw("ALP");

    gRight->SetLineColor(kRed + 1);
    gRight->SetMarkerColor(kRed + 1);
    gRight->SetMarkerStyle(21);
    gRight->SetLineWidth(2);
    gRight->Draw("LP SAME");

    gTop->SetLineColor(kGreen + 2);
    gTop->SetMarkerColor(kGreen + 2);
    gTop->SetMarkerStyle(22);
    gTop->SetLineWidth(2);
    gTop->Draw("LP SAME");

    TLegend* leg = new TLegend(0.62, 0.67, 0.90, 0.88);
    leg->AddEntry(gLeft, "Best end-left SiPM", "lp");
    leg->AddEntry(gRight, "Best end-right SiPM", "lp");
    leg->AddEntry(gTop, "Best top SiPM", "lp");
    leg->Draw();

    if (gLeft->GetHistogram()) {
        TLine* leftEnd = new TLine(-700.0, 0.0, -700.0, gLeft->GetHistogram()->GetMaximum());
        leftEnd->SetLineStyle(7);
        leftEnd->SetLineColor(kGray + 2);
        leftEnd->Draw();

        TLine* rightEnd = new TLine(700.0, 0.0, 700.0, gLeft->GetHistogram()->GetMaximum());
        rightEnd->SetLineStyle(7);
        rightEnd->SetLineColor(kGray + 2);
        rightEnd->Draw();
    }

    TLatex label;
    label.SetNDC(true);
    label.SetTextSize(0.035);
    label.DrawLatex(0.12, 0.92, Form("Ranking scan por SiPM individual (%s, n=%d)", ModeName(mode), nParam));
}
