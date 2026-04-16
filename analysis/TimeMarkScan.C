#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TF1.h"
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

const char* ModeName(int mode) {
    switch (mode) {
        case kFirstPhoton: return "first-photon";
        case kNthPhoton: return "nth-photon";
        case kMeanFirstN: return "mean-first-n";
        default: return "unknown";
    }
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

}  // namespace

void TimeMarkScan(Int_t mode = kMeanFirstN, Int_t nParam = 3, bool drawFits = true) {
    TChain chain("sipm_hits");
    if (!AddScanFiles(chain)) {
        std::cerr << "Error: no se encontraron archivos del scan." << std::endl;
        return;
    }

    std::cout << "TimeMarkScan: modo = " << ModeName(mode)
              << ", nParam = " << nParam << std::endl;
    std::cout << "TimeMarkScan: entradas totales en el chain = "
              << chain.GetEntries() << std::endl;

    Int_t event_id = -1;
    Int_t face_type = -1;
    Double_t time_ns = 0.0;
    Double_t gun_x_mm = 0.0;

    chain.SetBranchAddress("event_id", &event_id);
    chain.SetBranchAddress("face_type", &face_type);
    chain.SetBranchAddress("time_ns", &time_ns);
    chain.SetBranchAddress("gun_x_mm", &gun_x_mm);

    using HitPair = std::pair<std::vector<Double_t>, std::vector<Double_t>>;
    std::map<Int_t, std::map<Int_t, HitPair>> eventsByX;

    for (Long64_t i = 0; i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);
        const Int_t xRounded = static_cast<Int_t>(std::lround(gun_x_mm));
        if (face_type == 0) {
            eventsByX[xRounded][event_id].first.push_back(time_ns);
        } else if (face_type == 1) {
            eventsByX[xRounded][event_id].second.push_back(time_ns);
        }
    }

    if (eventsByX.empty()) {
        std::cerr << "Error: no hay eventos con face_type 0/1 en el chain." << std::endl;
        return;
    }

    TGraphErrors* grRes = new TGraphErrors();
    grRes->SetName("gTimeMarkResolution");
    std::vector<TH1D*> fitHists;
    int pointIdx = 0;

    for (const auto& [xPos, events] : eventsByX) {
        TH1D hDiff(Form("hDiff_x%d", xPos),
                   Form("x = %d mm; (t_{L}^{mark} - t_{R}^{mark}) / 2 [ns]; eventos", xPos),
                   160, -8.0, 8.0);
        hDiff.SetDirectory(nullptr);

        int nAccepted = 0;
        for (const auto& [evId, hits] : events) {
            (void)evId;
            const Double_t tL = ComputeTimeMark(hits.first, mode, static_cast<size_t>(nParam));
            const Double_t tR = ComputeTimeMark(hits.second, mode, static_cast<size_t>(nParam));
            if (tL < 0.0 || tR < 0.0) continue;
            hDiff.Fill((tL - tR) / 2.0);
            ++nAccepted;
        }

        if (hDiff.GetEntries() < 30) {
            std::cout << "  x = " << xPos
                      << " mm descartado: solo "
                      << hDiff.GetEntries()
                      << " entradas validas" << std::endl;
            continue;
        }

        const Double_t mean = hDiff.GetMean();
        const Double_t rms = hDiff.GetRMS();
        if (!(rms > 0.0)) continue;

        TF1 fit(Form("fit_x%d", xPos), "gaus", mean - 2.0 * rms, mean + 2.0 * rms);
        fit.SetParameters(hDiff.GetMaximum(), mean, rms);
        fit.SetParLimits(2, 1e-6, 5.0 * rms);
        const Int_t fitStatus = hDiff.Fit(&fit, "QNR");

        Double_t sigmaNs = rms;
        Double_t sigmaErrNs = (hDiff.GetEntries() > 1)
            ? rms / std::sqrt(2.0 * (hDiff.GetEntries() - 1.0))
            : 0.0;

        if (fitStatus == 0 && std::isfinite(fit.GetParameter(2)) && fit.GetParameter(2) > 0.0) {
            sigmaNs = std::abs(fit.GetParameter(2));
            sigmaErrNs = fit.GetParError(2);
        }

        if (drawFits) {
            TH1D* hClone = static_cast<TH1D*>(hDiff.Clone(Form("hDraw_x%d", xPos)));
            hClone->SetDirectory(nullptr);
            TF1* fitClone = static_cast<TF1*>(fit.Clone(Form("fitDraw_x%d", xPos)));
            hClone->GetListOfFunctions()->Clear();
            if (fitStatus == 0) hClone->GetListOfFunctions()->Add(fitClone);
            fitHists.push_back(hClone);
        }

        grRes->SetPoint(pointIdx, xPos, sigmaNs * 1000.0);
        grRes->SetPointError(pointIdx, 0.0, sigmaErrNs * 1000.0);

        std::cout << "  x = " << xPos
                  << " mm, accepted = " << nAccepted
                  << ", mean = " << mean
                  << " ns, sigma = " << sigmaNs * 1000.0
                  << " ps" << std::endl;
        ++pointIdx;
    }

    std::cout << "TimeMarkScan: puntos cargados en la grafica = "
              << pointIdx << std::endl;

    if (pointIdx == 0) {
        std::cerr << "Error: no se pudo construir ningun punto del scan." << std::endl;
        return;
    }

    TCanvas* cScan = new TCanvas("cTimeMarkScan", "Time-mark scan", 900, 550);
    grRes->SetTitle(Form("Resolucion temporal por time mark (%s, n=%d); Posicion X [mm]; #sigma_{t} [ps]",
                         ModeName(mode), nParam));
    grRes->SetMarkerStyle(20);
    grRes->SetMarkerColor(kAzure + 2);
    grRes->SetLineColor(kAzure + 2);
    grRes->SetLineWidth(2);
    grRes->SetMinimum(0.0);
    grRes->Draw("APL");

    TLine* target = new TLine(-700.0, 100.0, 700.0, 100.0);
    target->SetLineColor(kRed + 1);
    target->SetLineStyle(7);
    target->Draw();

    if (drawFits && !fitHists.empty()) {
        const int nPads = static_cast<int>(fitHists.size());
        const int nCols = 3;
        const int nRows = (nPads + nCols - 1) / nCols;
        TCanvas* cFits = new TCanvas("cTimeMarkFits", "Time-mark histograms and fits",
                                     1500, 430 * nRows);
        cFits->Divide(nCols, nRows);

        for (int i = 0; i < nPads; ++i) {
            cFits->cd(i + 1);
            fitHists[i]->SetLineColor(kBlack);
            fitHists[i]->SetMarkerStyle(20);
            fitHists[i]->SetMarkerSize(0.5);
            fitHists[i]->Draw("E");

            TF1* fit = fitHists[i]->GetListOfFunctions()->GetSize() > 0
                ? dynamic_cast<TF1*>(fitHists[i]->GetListOfFunctions()->At(0))
                : nullptr;

            if (fit) {
                fit->SetLineColor(kRed + 1);
                fit->SetLineWidth(2);
                fit->Draw("same");

                TLatex label;
                label.SetNDC(true);
                label.SetTextSize(0.045);
                label.DrawLatex(0.13, 0.85, Form("#sigma = %.0f ps", fit->GetParameter(2) * 1000.0));
            }
        }
    }
}
