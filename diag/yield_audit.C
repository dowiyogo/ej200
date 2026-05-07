// diag/yield_audit.C — auditoria de yield de fotones detectados
// Uso: root -l -q 'diag/yield_audit.C("photon_hits_run000.root")'

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

void yield_audit(const char* fname = "photon_hits_run000.root") {
    gStyle->SetOptStat(1);

    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: no pude abrir " << fname << "\n";
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(f->Get("sipm_hits"));
    if (!tree) {
        std::cerr << "ERROR: no hay TTree 'sipm_hits'\n";
        return;
    }

    Int_t event_id, face_type, global_id;
    Double_t time_ns, wl_nm, energy_eV;
    tree->SetBranchAddress("event_id",  &event_id);
    tree->SetBranchAddress("face_type", &face_type);
    tree->SetBranchAddress("global_id", &global_id);
    tree->SetBranchAddress("time_ns",   &time_ns);
    tree->SetBranchAddress("wl_nm",     &wl_nm);
    tree->SetBranchAddress("energy_eV", &energy_eV);

    TH1D* hNpe = new TH1D("hNpe",
                          "Fotones detectados por evento;N_{pe};Eventos",
                          200, 0, 200);
    TH1D* hTime = new TH1D("hTime",
                           "Tiempo de arribo;t [ns];Fotones",
                           200, 0, 50);
    TH1D* hWL = new TH1D("hWL",
                         "Longitud de onda detectada;lambda [nm];Fotones",
                         100, 350, 550);
    TH1D* hSiPMId = new TH1D("hSiPMId",
                             "Ocupacion por SiPM;Global ID;Fotones",
                             40, 0, 40);

    std::map<Int_t, Int_t> npeByEvent;
    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        hTime->Fill(time_ns);
        hWL->Fill(wl_nm);
        hSiPMId->Fill(global_id);
        npeByEvent[event_id]++;
    }

    Int_t nEvents = static_cast<Int_t>(npeByEvent.size());
    double sumNpe = 0;
    Int_t nZeroEvents = 0;

    for (const auto& [eid, n] : npeByEvent) {
        hNpe->Fill(n);
        sumNpe += n;
        if (n == 0) ++nZeroEvents;
    }

    std::cout << "\n=== YIELD AUDIT REPORT ===\n";
    std::cout << "  Archivo            : " << fname << "\n";
    std::cout << "  Entradas en TTree  : " << nEntries << "\n";
    std::cout << "  Eventos con hits   : " << nEvents << "\n";
    std::cout << "  Media N_pe/evento  : "
              << (nEvents > 0 ? sumNpe / nEvents : 0.0) << "\n";
    std::cout << "  Max N_pe/evento    : "
              << hNpe->GetXaxis()->GetBinCenter(hNpe->GetMaximumBin()) << "\n";

    std::map<int, long long> byFace;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        byFace[face_type]++;
    }
    std::cout << "  Por cara:\n";
    const char* faceNames[] = {"End-left", "End-right", "Top"};
    for (auto& [face, cnt] : byFace) {
        if (face >= 0 && face <= 2) {
            std::cout << "    " << faceNames[face] << " : " << cnt
                      << " (" << 100.0 * cnt / nEntries << "%)\n";
        }
    }

    const double expectedPerEvent = 2.0 * 10000.0;
    const double detEff = (nEvents > 0 && expectedPerEvent > 0)
                              ? 100.0 * (sumNpe / nEvents) / expectedPerEvent
                              : 0.0;
    std::cout << "  Estimado generados/evento : ~" << expectedPerEvent << "\n";
    std::cout << "  Eficiencia total estimada : ~" << detEff << " %\n";
    std::cout << "  REFERENCIA Betancourt 2020: ~175 ph/evt (ambos ends)\n";
    std::cout << "  DEFICIT factor            : ~"
              << (sumNpe > 0 ? 175.0 / (sumNpe / nEvents) : -1.0) << "x\n";
    std::cout << "==========================\n\n";

    TCanvas* c = new TCanvas("c_audit", "Yield Audit", 1200, 800);
    c->Divide(2, 2);

    c->cd(1); hNpe->Draw(); hNpe->SetFillColor(kAzure - 9);
    c->cd(2); hTime->Draw(); hTime->SetLineColor(kBlue + 1);
    c->cd(3); hWL->Draw(); hWL->SetFillColor(kGreen - 9);
    c->cd(4); hSiPMId->Draw(); hSiPMId->SetFillColor(kOrange - 3);

    c->SaveAs("diag/yield_audit.pdf");
    std::cout << "  -> diag/yield_audit.pdf\n";
}
