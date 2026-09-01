#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"

// Función para obtener un tiempo de trigger estable
Double_t GetStableTime(const std::vector<Double_t>& times, size_t nPhotons = 1) {
    if (times.empty()) return -1.0;

    const size_t nUse = std::min(times.size(), nPhotons);
    std::vector<Double_t> sortedTimes = times;
    std::sort(sortedTimes.begin(), sortedTimes.end());

    Double_t sum = 0;
    for (size_t i = 0; i < nUse; ++i) sum += sortedTimes[i];
    return sum / static_cast<Double_t>(nUse);
}

void ResolutionScan_v2(size_t nPhotons = 1, bool drawFits = true) {
    TChain *chain = new TChain("sipm_hits");
    bool filesFound = false;

    const std::vector<TString> searchDirs = {"", "build/", "../build/"};
    for (int i = 0; i <= 20; ++i) {
        bool addedThisRun = false;
        for (const auto& dir : searchDirs) {
            TString name = Form("%sphoton_hits_run%03d.root", dir.Data(), i);
            if (!gSystem->AccessPathName(name)) {
                chain->Add(name);
                filesFound = true;
                addedThisRun = true;
                break;
            }
        }
        if (!addedThisRun) {
            std::cerr << "Aviso: no se encontro photon_hits_run"
                      << Form("%03d", i) << ".root en ., build/ o ../build/\n";
        }
    }

    if (!filesFound) {
        std::cerr << "Error: No se encontraron archivos run000...020" << std::endl;
        return;
    }

    std::cout << "ResolutionScan_v2: entradas totales en el chain = "
              << chain->GetEntries() << std::endl;

    Int_t event_id, face_type;
    Double_t time_ns, gun_x;
    chain->SetBranchAddress("event_id", &event_id);
    chain->SetBranchAddress("face_type", &face_type);
    chain->SetBranchAddress("time_ns", &time_ns);
    chain->SetBranchAddress("gun_x_mm", &gun_x);

    // Agrupamos: PosicionX redondeada [mm] -> Evento -> {Tiempos_L, Tiempos_R}
    std::map<Int_t, std::map<Int_t, std::pair<std::vector<Double_t>, std::vector<Double_t>>>> data;
    
    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        const Int_t gun_x_rounded = static_cast<Int_t>(std::lround(gun_x));
        if (face_type == 0) data[gun_x_rounded][event_id].first.push_back(time_ns);
        else if (face_type == 1) data[gun_x_rounded][event_id].second.push_back(time_ns);
    }

    if (data.size() == 1) {
        std::cerr << "Aviso: solo hay una posicion X en el ntuple (gun_x_mm = "
                  << data.begin()->first
                  << " mm). Si esperabas un scan, vuelve a generar los ROOT con "
                  << "/muon/gunX en lugar de /gun/position." << std::endl;
    }

    TGraphErrors *grRes = new TGraphErrors();
    std::vector<TH1D*> fitHists;
    int pIdx = 0;

    for (auto const& [x, events] : data) {
        TString hName = Form("hDiff_x%d", x);
        TH1D hDiff(hName, hName, 200, -10, 10);
        hDiff.SetDirectory(nullptr);

        for (auto const& [id, hits] : events) {
            Double_t tL = GetStableTime(hits.first, nPhotons);
            Double_t tR = GetStableTime(hits.second, nPhotons);
            
            if (tL >= 0.0 && tR >= 0.0) hDiff.Fill((tL - tR) / 2.0);
        }

        if (hDiff.GetEntries() > 30) {
            Double_t mean = hDiff.GetMean();
            Double_t rms = hDiff.GetRMS();
            if (rms <= 0.0) continue;

            Double_t sigma    = rms;
            Double_t sigmaErr = (hDiff.GetEntries() > 1)
                ? rms / std::sqrt(2.0 * (hDiff.GetEntries() - 1.0))
                : 0.0;

            hDiff.Fit("gaus", "Q", "", mean - 2*rms, mean + 2*rms);
            
            TF1 *f = hDiff.GetFunction("gaus");
            if (f && std::isfinite(f->GetParameter(2)) && f->GetParameter(2) > 0.0) {
                sigma    = f->GetParameter(2);
                sigmaErr = f->GetParError(2);
            }

            TH1D* hClone = nullptr;
            if (drawFits) {
                hClone = static_cast<TH1D*>(hDiff.Clone(Form("hDiff_draw_x%d", x)));
                hClone->SetDirectory(nullptr);
                hClone->SetTitle(Form("x = %d mm; (t_{L} - t_{R}) / 2 [ns]; eventos", x));
                fitHists.push_back(hClone);
            }

            grRes->SetPoint(pIdx, x, sigma * 1000.0);
            grRes->SetPointError(pIdx, 0, sigmaErr * 1000.0);
            std::cout << "  x = " << x
                      << " mm, entries = " << hDiff.GetEntries()
                      << ", mean = " << mean
                      << " ns, sigma = " << sigma * 1000.0
                      << " ps" << std::endl;
            pIdx++;
        } else {
            std::cout << "  x = " << x
                      << " mm descartado: solo "
                      << hDiff.GetEntries()
                      << " entradas" << std::endl;
        }
    }

    std::cout << "ResolutionScan_v2: puntos cargados en la grafica = "
              << pIdx << std::endl;
    if (pIdx == 0) {
        std::cerr << "Error: no se pudo construir ningun punto. "
                  << "Revisa los mensajes por posicion de arriba." << std::endl;
        return;
    }

    TCanvas *c = new TCanvas("c", "Scan de Resolucion", 800, 500);
    grRes->SetTitle("Resolucion Temporal Intrinsica; Posicion X [mm]; #sigma_{t} [ps]");
    grRes->SetMarkerStyle(20);
    grRes->SetMarkerColor(kAzure+2);
    grRes->SetLineColor(kAzure+2);
    grRes->SetMinimum(0.0);
    grRes->Draw("APL");
    
    TLine *target = new TLine(-700, 50, 700, 50);
    target->SetLineStyle(7); target->SetLineColor(kRed); target->Draw();

    if (drawFits && !fitHists.empty()) {
        const int nPads = static_cast<int>(fitHists.size());
        const int nCols = 3;
        const int nRows = (nPads + nCols - 1) / nCols;

        TCanvas *cFits = new TCanvas("cFits", "Histograms and fits by position", 1500, 450 * nRows);
        cFits->Divide(nCols, nRows);

        for (int i = 0; i < nPads; ++i) {
            cFits->cd(i + 1);
            fitHists[i]->SetLineColor(kBlack);
            fitHists[i]->SetMarkerStyle(20);
            fitHists[i]->SetMarkerSize(0.6);
            fitHists[i]->Draw("E");

            TF1* fit = fitHists[i]->GetFunction("gaus");
            if (fit) {
                fit->SetLineColor(kRed + 1);
                fit->SetLineWidth(2);
                fit->Draw("same");

                TLatex label;
                label.SetNDC(true);
                label.SetTextSize(0.045);
                label.DrawLatex(0.14, 0.85, Form("#sigma = %.0f ps", fit->GetParameter(2) * 1000.0));
            } else {
                TLatex label;
                label.SetNDC(true);
                label.SetTextSize(0.045);
                label.DrawLatex(0.14, 0.85, "sin ajuste gaussiano");
            }
        }
    }
}
