#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"

void analyze_dCFD(const char* firstFile = "photon_hits_run000.root") {
    // 1. Configuración del TChain para procesar del run000 al run020
    TChain *chain = new TChain("sipm_hits"); // El nombre correcto según tu código
    for (int i = 0; i <= 20; ++i) {
        TString name = Form("photon_hits_run%03d.root", i);
        // Solo agregamos si el archivo existe
        if (!gSystem->AccessPathName(name)) {
            chain->Add(name);
        }
    }

    if (chain->GetEntries() == 0) {
        std::cerr << "Error: No se encontraron entradas. Asegurate de que los archivos run000...run020 esten aqui." << std::endl;
        return;
    }

    // 2. Variables del TTree (basado en tu SiPMSD.cc)
    Int_t event_id, face_type;
    Double_t time_ns;
    chain->SetBranchAddress("event_id", &event_id);
    chain->SetBranchAddress("face_type", &face_type);
    chain->SetBranchAddress("time_ns", &time_ns);

    // 3. Histograma para el tiempo de trigger (por ejemplo, cara Left)
    TH1D *hTriggerTime = new TH1D("hTriggerTime", "Tiempo de Trigger (Primeros fotones);t_{trigger} [ns];Cuentas", 1000, 0, 50);

    // 4. Lógica de análisis
    // Como no tienes waveforms digitalizadas, simularemos el trigger 
    // usando el tiempo del N-ésimo fotón o el promedio de los primeros.
    
    std::cout << "Procesando " << chain->GetEntries() << " hits de fotones..." << std::endl;

    // Para este análisis, agruparemos por evento para encontrar el "tiempo de llegada"
    // (Esto es lo más cercano a un dCFD sin tener la waveform completa)
    Long64_t nEntries = chain->GetEntries();
    std::map<Int_t, std::vector<Double_t>> eventHits;

    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);
        // Analizaremos solo la cara Left (face_type == 0) para este ejemplo
        if (face_type == 0) {
            eventHits[event_id].push_back(time_ns);
        }
    }

    for (auto const& [id, times] : eventHits) {
        if (times.empty()) continue;
        
        // Ordenamos los tiempos de arribo de los fotones del evento
        std::vector<Double_t> sortedTimes = times;
        std::sort(sortedTimes.begin(), sortedTimes.end());

        // Simulamos un trigger simple: el tiempo del 5to fotón detectado 
        // (Esto emula un umbral de baja intensidad)
        if (sortedTimes.size() >= 5) {
            hTriggerTime->Fill(sortedTimes[4]); 
        }
    }

    // 5. Visualización
    TCanvas *c = new TCanvas("c", "Timing Resolution", 800, 600);
    hTriggerTime->Draw();
    hTriggerTime->Fit("gaus");

    Double_t res = hTriggerTime->GetFunction("gaus")->GetParameter(2);
    std::cout << "\nResolucion temporal estimada: " << res * 1000.0 << " ps" << std::endl;
}