#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

/**
 * Macro para CINT/CLING: Analiza los hits de SiPM del Timing Detector.
 * Uso: root -l 'analyze_hits.C("photon_hits_run000.root")'
 */
void analyze_hits(const char* filename = "photon_hits_run000.root") {
    // 1. Configuración de estilo para publicaciones
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetLabelSize(0.035, "XYZ");

    // 2. Abrir el archivo y obtener el árbol
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        printf("Error: No se pudo abrir el archivo %s\n", filename);
        return;
    }

    TTree *tree = (TTree*)file->Get("sipm_hits");
    if (!tree) {
        printf("Error: No se encontró el TTree 'sipm_hits'\n");
        return;
    }

    // 3. Variables para lectura (basadas en la estructura de tu simulación)
    Int_t face_type, global_id;
    Double_t time_ns, wl_nm;
    
    tree->SetBranchAddress("face_type", &face_type);
    tree->SetBranchAddress("global_id", &global_id);
    tree->SetBranchAddress("time_ns",   &time_ns);
    tree->SetBranchAddress("wl_nm",     &wl_nm);

    // 4. Definición de histogramas
    // Separamos el tiempo de arribo por cara (End Left, End Right, Top)
    TH1D *hTimeLeft  = new TH1D("hTimeLeft",  "Arrival Time;Time [ns];Cuentas", 100, 0, 40);
    TH1D *hTimeRight = new TH1D("hTimeRight", "Arrival Time;Time [ns];Cuentas", 100, 0, 40);
    TH1D *hTimeTop   = new TH1D("hTimeTop",   "Arrival Time;Time [ns];Cuentas", 100, 0, 40);
    
    TH1D *hWavelength = new TH1D("hWavelength", "Detected Spectrum;Wavelength [nm];Cuentas", 80, 350, 550);
    TH1D *hOccupancy  = new TH1D("hOccupancy",  "SiPM Occupancy;Global ID;Hits", 60, 0, 60);

    // Colores consistentes con tu análisis previo
    hTimeLeft->SetLineColor(kBlue+1);
    hTimeRight->SetLineColor(kRed+1);
    hTimeTop->SetLineColor(kGreen+2);

    // 5. Loop de eventos
    Long64_t nEntries = tree->GetEntries();
    std::cout << "--> Analizando " << nEntries << " fotones detectados en " << filename << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        hWavelength->Fill(wl_nm);
        hOccupancy->Fill(global_id);

        // Clasificación por cara (FaceType definido en DetectorConstruction.cc)
        if (face_type == 0)      hTimeLeft->Fill(time_ns);
        else if (face_type == 1) hTimeRight->Fill(time_ns);
        else if (face_type == 2) hTimeTop->Fill(time_ns);
    }

    // 6. Visualización
    TCanvas *c1 = new TCanvas("c1", "SHiP Timing Detector Analysis", 1400, 900);
    c1->Divide(2, 2);

    // Panel 1: Espectro de tiempo comparado
    c1->cd(1);
    gPad->SetGrid();
    hTimeLeft->Draw("hist");
    hTimeRight->Draw("hist same");
    hTimeTop->Draw("hist same");
    
    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(hTimeLeft, "End Left (-X)", "l");
    leg->AddEntry(hTimeRight, "End Right (+X)", "l");
    leg->AddEntry(hTimeTop, "Top (+Y)", "l");
    leg->Draw();

    // Panel 2: Espectro óptico (PDE ya aplicado en SiPMSD.cc)
    c1->cd(2);
    hWavelength->SetFillColor(kAzure-9);
    hWavelength->Draw("hist");

    // Panel 3: Ocupación de los SiPMs
    c1->cd(3);
    hOccupancy->SetFillColor(kOrange-3);
    hOccupancy->Draw("hist");
    
    // Panel 4: Información resumida
    c1->cd(4);
    hOccupancy->Draw("hbar"); // Vista horizontal para identificar IDs

    c1->Update();
    std::cout << "--> Análisis finalizado exitosamente." << std::endl;
}