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

void ResolutionScan_v2(size_t nPhotons = 1) {
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

            hDiff.Fit("gaus", "Q", "", mean - 2*rms, mean + 2*rms);
            
            TF1 *f = hDiff.GetFunction("gaus");
            if (f) {
                grRes->SetPoint(pIdx, x, f->GetParameter(2) * 1000.0);
                grRes->SetPointError(pIdx, 0, f->GetParError(2) * 1000.0);
                pIdx++;
            }
        }
    }

    TCanvas *c = new TCanvas("c", "Scan de Resolucion", 800, 500);
    grRes->SetTitle("Resolucion Temporal Intrinsica; Posicion X [mm]; #sigma_{t} [ps]");
    grRes->SetMarkerStyle(20);
    grRes->SetMarkerColor(kAzure+2);
    grRes->Draw("APL");
    
    TLine *target = new TLine(-700, 50, 700, 50);
    target->SetLineStyle(7); target->SetLineColor(kRed); target->Draw();
}
