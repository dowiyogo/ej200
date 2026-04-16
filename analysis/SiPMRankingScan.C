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
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
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
    double sigma_core_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_err_ns = std::numeric_limits<double>::quiet_NaN();
    double rms_ns = std::numeric_limits<double>::quiet_NaN();
};

struct FitSummary {
    double mean_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_core_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_err_ns = std::numeric_limits<double>::quiet_NaN();
    double rms_ns = std::numeric_limits<double>::quiet_NaN();
    int n = 0;
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

FitSummary FitCoreSigma(const std::vector<double>& marks, const TString& name) {
    FitSummary out;
    out.n = static_cast<int>(marks.size());
    if (marks.empty()) return out;

    out.mean_ns = Mean(marks);
    out.rms_ns = StdDev(marks);
    if (marks.size() < 20 || !(out.rms_ns > 0.0)) {
        out.sigma_core_ns = out.rms_ns;
        out.sigma_err_ns = (marks.size() > 1)
            ? out.rms_ns / std::sqrt(2.0 * (marks.size() - 1.0))
            : 0.0;
        return out;
    }

    const double lo = out.mean_ns - 4.0 * out.rms_ns;
    const double hi = out.mean_ns + 4.0 * out.rms_ns;
    TH1D h(name, name, 80, lo, hi);
    h.SetDirectory(nullptr);
    for (double t : marks) h.Fill(t);

    const double fitLo = out.mean_ns - 1.5 * out.rms_ns;
    const double fitHi = out.mean_ns + 1.5 * out.rms_ns;
    TF1 fit(name + "_gaus", "gaus", fitLo, fitHi);
    fit.SetParameters(std::max(1.0, h.GetMaximum()), out.mean_ns, std::max(0.02, out.rms_ns));
    fit.SetParLimits(2, 1e-4, 5.0 * out.rms_ns);

    const int fitStatus = h.Fit(&fit, "QNR");
    if (fitStatus == 0 && std::isfinite(fit.GetParameter(2)) && fit.GetParameter(2) > 0.0) {
        out.mean_ns = fit.GetParameter(1);
        out.sigma_core_ns = std::abs(fit.GetParameter(2));
        out.sigma_err_ns = fit.GetParError(2);
    } else {
        out.sigma_core_ns = out.rms_ns;
        out.sigma_err_ns = out.rms_ns / std::sqrt(2.0 * (marks.size() - 1.0));
    }
    return out;
}

template <typename T>
bool Contains(const std::vector<T>& v, const T& value) {
    return std::find(v.begin(), v.end(), value) != v.end();
}

}  // namespace

void SiPMRankingScan(Int_t mode = kFirstPhoton, Int_t nParam = 1, Int_t minEvents = 50) {
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

    using EventKey = std::pair<int, int>;
    using EventMap = std::map<EventKey, std::vector<Double_t>>;
    std::map<int, EventMap> hitsBySensor;
    std::map<int, std::pair<int, int>> sensorMeta;
    std::vector<int> xPositions;

    for (Long64_t i = 0; i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);
        const int x = static_cast<int>(std::lround(gun_x_mm));
        hitsBySensor[global_id][{x, event_id}].push_back(time_ns);
        sensorMeta[global_id] = {face_type, local_id};
        if (!Contains(xPositions, x)) xPositions.push_back(x);
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
                            Form("Resolucion del nucleo temporal por SiPM (%s, n=%d); X del muon [mm]; global_id del SiPM",
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

    for (int sensor : sensors) {
        const int face = sensorMeta[sensor].first;
        const int local = sensorMeta[sensor].second;

        std::map<int, std::vector<double>> marksByX;
        for (const auto& ev : hitsBySensor[sensor]) {
            const int x = ev.first.first;
            const double tMark = ComputeTimeMark(ev.second, mode, static_cast<size_t>(nParam));
            if (tMark >= 0.0) marksByX[x].push_back(tMark);
        }

        for (const auto& xv : marksByX) {
            const int x = xv.first;
            const FitSummary fit = FitCoreSigma(xv.second, Form("h_sensor_%d_x_%d", sensor, x));

            hCounts->SetBinContent(hCounts->FindBin(x, sensor), fit.n);
            if (fit.n >= minEvents) {
                hMean->SetBinContent(hMean->FindBin(x, sensor), fit.mean_ns);
                hSigma->SetBinContent(hSigma->FindBin(x, sensor), fit.sigma_core_ns * 1000.0);
                rankingByX[x].push_back(
                    SensorStats{face, sensor, local, fit.n, fit.mean_ns, fit.sigma_core_ns, fit.sigma_err_ns, fit.rms_ns});
            }
        }
    }

    std::cout << "\nRanking por posicion X\n";
    for (int x : xPositions) {
        auto stats = rankingByX[x];
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
                      return a.sigma_core_ns < b.sigma_core_ns;
                  });

        std::cout << "  x = " << x << " mm" << std::endl;
        std::cout << "    Earliest 3 sensors:" << std::endl;
        for (size_t i = 0; i < std::min<size_t>(3, byMean.size()); ++i) {
            const auto& s = byMean[i];
            std::cout << "      " << i + 1
                      << ". gid=" << s.global_id
                      << " (" << FaceName(s.face_type) << ", local=" << s.local_id << ")"
                      << ", mean=" << s.mean_ns << " ns"
                      << ", sigma_core=" << s.sigma_core_ns * 1000.0 << " ps"
                      << ", rms=" << s.rms_ns * 1000.0 << " ps"
                      << ", n=" << s.n_events << std::endl;
        }
        std::cout << "    Best 3 sigma-core sensors:" << std::endl;
        for (size_t i = 0; i < std::min<size_t>(3, bySigma.size()); ++i) {
            const auto& s = bySigma[i];
            std::cout << "      " << i + 1
                      << ". gid=" << s.global_id
                      << " (" << FaceName(s.face_type) << ", local=" << s.local_id << ")"
                      << ", sigma_core=" << s.sigma_core_ns * 1000.0 << " ps"
                      << ", rms=" << s.rms_ns * 1000.0 << " ps"
                      << ", mean=" << s.mean_ns << " ns"
                      << ", n=" << s.n_events << std::endl;
        }
    }

    std::vector<double> xLeft, yLeft, xRight, yRight, xTop, yTop;
    std::vector<double> xBest1, yBest1, xBest2, yBest2, xBest3, yBest3;

    for (int x : xPositions) {
        auto stats = rankingByX[x];
        if (stats.empty()) continue;

        double bestLeft = std::numeric_limits<double>::infinity();
        double bestRight = std::numeric_limits<double>::infinity();
        double bestTop = std::numeric_limits<double>::infinity();

        std::sort(stats.begin(), stats.end(),
                  [](const SensorStats& a, const SensorStats& b) {
                      return a.sigma_core_ns < b.sigma_core_ns;
                  });

        for (const auto& s : stats) {
            const double sigmaPs = s.sigma_core_ns * 1000.0;
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

        for (int nTake = 1; nTake <= 3; ++nTake) {
            if (static_cast<int>(stats.size()) < nTake) continue;

            std::vector<int> chosen;
            for (int i = 0; i < nTake; ++i) chosen.push_back(stats[i].global_id);

            std::map<int, std::vector<double>> combinedByEvent;
            for (int sensor : chosen) {
                for (const auto& ev : hitsBySensor[sensor]) {
                    if (ev.first.first != x) continue;
                    const double tMark = ComputeTimeMark(ev.second, mode, static_cast<size_t>(nParam));
                    if (tMark < 0.0) continue;
                    combinedByEvent[ev.first.second].push_back(tMark);
                }
            }

            std::vector<double> comboMarks;
            for (const auto& ev : combinedByEvent) {
                if (static_cast<int>(ev.second.size()) != nTake) continue;
                comboMarks.push_back(Mean(ev.second));
            }

            const FitSummary fit = FitCoreSigma(comboMarks, Form("h_combo_%d_x_%d", nTake, x));
            if (!(fit.n >= minEvents) || !(fit.sigma_core_ns > 0.0)) continue;

            if (nTake == 1) {
                xBest1.push_back(x);
                yBest1.push_back(fit.sigma_core_ns * 1000.0);
            } else if (nTake == 2) {
                xBest2.push_back(x);
                yBest2.push_back(fit.sigma_core_ns * 1000.0);
            } else {
                xBest3.push_back(x);
                yBest3.push_back(fit.sigma_core_ns * 1000.0);
            }
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

    gLeft->SetTitle("Mejor resolucion temporal por cara; X del muon [mm]; mejor #sigma_{core} por cara [ps]");
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

    TCanvas* cCombo = new TCanvas("cBestCombinations", "Best sensor combinations", 1000, 550);
    TGraph* gBest1 = new TGraph(static_cast<int>(xBest1.size()), xBest1.data(), yBest1.data());
    TGraph* gBest2 = new TGraph(static_cast<int>(xBest2.size()), xBest2.data(), yBest2.data());
    TGraph* gBest3 = new TGraph(static_cast<int>(xBest3.size()), xBest3.data(), yBest3.data());

    gBest1->SetTitle("Mejor combinacion de sensores por posicion; X del muon [mm]; #sigma_{core} [ps]");
    gBest1->SetLineColor(kBlack);
    gBest1->SetMarkerColor(kBlack);
    gBest1->SetMarkerStyle(20);
    gBest1->SetLineWidth(2);
    gBest1->SetMinimum(0.0);
    gBest1->Draw("ALP");

    gBest2->SetLineColor(kOrange + 7);
    gBest2->SetMarkerColor(kOrange + 7);
    gBest2->SetMarkerStyle(21);
    gBest2->SetLineWidth(2);
    gBest2->Draw("LP SAME");

    gBest3->SetLineColor(kViolet + 1);
    gBest3->SetMarkerColor(kViolet + 1);
    gBest3->SetMarkerStyle(22);
    gBest3->SetLineWidth(2);
    gBest3->Draw("LP SAME");

    TLegend* legCombo = new TLegend(0.60, 0.68, 0.90, 0.88);
    legCombo->AddEntry(gBest1, "Mejor 1 SiPM", "lp");
    legCombo->AddEntry(gBest2, "Promedio de mejores 2 SiPMs", "lp");
    legCombo->AddEntry(gBest3, "Promedio de mejores 3 SiPMs", "lp");
    legCombo->Draw();

    TLatex label;
    label.SetNDC(true);
    label.SetTextSize(0.035);
    label.DrawLatex(0.12, 0.92, Form("Ranking scan por SiPM individual (%s, n=%d)", ModeName(mode), nParam));
}
