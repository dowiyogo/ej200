#include <iostream>
#include <vector>
#include <map>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

void ResolutionScan() {
    // 1. Configuración del Chain para los 21 puntos del scan
    TChain *chain = new TChain("sipm_hits");
    for (int i = 0; i <= 20; ++i) {
        chain->Add(Form("photon_hits_run%03d.root", i));
    }

    // Variables de lectura
    Int_t event_id, face_type;
    Double_t time_ns, gun_x;
    chain->SetBranchAddress("event_id", &event_id);
    chain->SetBranchAddress("face_type", &face_type);
    chain->SetBranchAddress("time_ns", &time_ns);
    chain->SetBranchAddress("gun_x_mm", &gun_x); // Usamos gun_x para separar los grupos

    // 2. Estructura para organizar: PosicionX -> EventID -> {tiempos_left, tiempos_right}
    std::map<Double_t, std::map<Int_t, std::pair<std::vector<Double_t>, std::vector<Double_t>>>> data;

    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);
        if (face_type == 0) data[gun_x][event_id].first.push_back(time_ns);
        if (face_type == 1) data[gun_x][event_id].second.push_back(time_ns);
    }

    // 3. Procesamiento por posición
    TGraphErrors *grRes = new TGraphErrors();
    int pointIdx = 0;

    for (auto const& [xPos, events] : data) {
        TString hName = Form("hDiff_x%.0f", xPos);
        // Histograma de la diferencia de tiempos para esta posición X
        TH1D *hDiff = new TH1D(hName, hName, 200, -10, 10);

        for (auto const& [evId, hits] : events) {
            if (hits.first.empty() || hits.second.empty()) continue;

            // Trigger por "Primer Fotón" en cada lado (lo más básico)
            Double_t tL = *std::min_element(hits.first.begin(), hits.first.end());
            Double_t tR = *std::min_element(hits.second.begin(), hits.second.end());
            
            hDiff->Fill((tL - tR) / 2.0); 
        }

        if (hDiff->GetEntries() > 20) {
            hDiff->Fit("gaus", "Q");
            TF1 *fit = hDiff->GetFunction("gaus");
            if (fit) {
                Double_t sigma = fit->GetParameter(2);
                Double_t sigmaErr = fit->GetParError(2);
                
                grRes->SetPoint(pointIdx, xPos, sigma * 1000.0); // Convertir a ps
                grRes->SetPointError(pointIdx, 0, sigmaErr * 1000.0);
                pointIdx++;
            }
        }
        delete hDiff;
    }

    // 4. Graficar Resolución vs X
    TCanvas *c2 = new TCanvas("c2", "Analisis de Resolucion por Posicion", 900, 600);
    grRes->SetTitle("Resolucion Temporal vs Posicion del Beam; Posicion X [mm]; #sigma_{t} [ps]");
    grRes->SetMarkerStyle(21);
    grRes->SetMarkerColor(kBlue+2);
    grRes->Draw("APL");
    
    // Línea de referencia del objetivo SHiP (50 ps)
    TLine *line = new TLine(-700, 50, 700, 50);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw();
}