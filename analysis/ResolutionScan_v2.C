#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"

// Función para obtener un tiempo de trigger estable
Double_t GetStableTime(std::vector<Double_t>& times, size_t nPhotons = 5) {
    if (times.size() < nPhotons) return -1.0;
    std::sort(times.begin(), times.end());
    Double_t sum = 0;
    for(size_t i=0; i<nPhotons; ++i) sum += times[i];
    return sum / nPhotons;
}

void ResolutionScan_v2() {
    TChain *chain = new TChain("sipm_hits");
    bool filesFound = false;
    for (int i = 0; i <= 20; ++i) {
        TString name = Form("photon_hits_run%03d.root", i);
        if (!gSystem->AccessPathName(name)) {
            chain->Add(name);
            filesFound = true;
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

    // Agrupamos: Posicion -> Evento -> {Tiempos_L, Tiempos_R}
    std::map<Double_t, std::map<Int_t, std::pair<std::vector<Double_t>, std::vector<Double_t>>>> data;
    
    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        if (face_type == 0) data[gun_x][event_id].first.push_back(time_ns);
        else if (face_type == 1) data[gun_x][event_id].second.push_back(time_ns);
    }

    TGraphErrors *grRes = new TGraphErrors();
    int pIdx = 0;

    for (auto const& [x, events] : data) {
        TH1D hDiff("h", "", 100, -5, 5); // Ventana de 10ns para la diferencia

        for (auto const& [id, hits] : events) {
            Double_t tL = GetStableTime(const_cast<std::vector<Double_t>&>(hits.first));
            Double_t tR = GetStableTime(const_cast<std::vector<Double_t>&>(hits.second));
            
            if (tL > 0 && tR > 0) hDiff.Fill((tL - tR) / 2.0);
        }

        if (hDiff.GetEntries() > 30) {
            // Ajuste inteligente: centrar en el RMS para evitar fallos
            Double_t mean = hDiff.GetMean();
            Double_t rms = hDiff.GetRMS();
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