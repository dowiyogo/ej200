#include "RunAction.hh"

#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4OpticalParameters.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

RunAction::RunAction() {
    // ── Register accumulables (solo en master thread) ────────────────────────
    // En modo MT, G4AccumulableManager registra en el master y los workers
    // lo heredan automáticamente vía Merge(). Registrar en cada worker
    // genera warnings de "already registered".
    if (IsMaster()) {
        auto* accMgr = G4AccumulableManager::Instance();
        accMgr->Register(fNEndLeft);
        accMgr->Register(fNEndRight);
        accMgr->Register(fNTop);
        accMgr->Register(fNEventsWithHits);
    }

    // ── Define ROOT ntuple (llamado en master y workers; G4Analysis lo maneja) ──
    auto* am = G4AnalysisManager::Instance();
    am->SetVerboseLevel(0);
    am->SetDefaultFileType("root");
    am->SetNtupleMerging(true);
    am->SetFileName("photon_hits");

    // Ntuple 0: one row per detected photon
    am->CreateNtuple("sipm_hits", "Detected optical photons in all SiPMs");
    am->CreateNtupleIColumn("event_id");   // 0
    am->CreateNtupleIColumn("face_type");  // 1  (0=end_left, 1=end_right, 2=top)
    am->CreateNtupleIColumn("global_id");  // 2  (0–35, unique across all SiPMs)
    am->CreateNtupleIColumn("local_id");   // 3  (index within face: 0-7 or 0-19)
    am->CreateNtupleDColumn("time_ns");    // 4
    am->CreateNtupleDColumn("energy_eV");  // 5
    am->CreateNtupleDColumn("wl_nm");      // 6
    am->CreateNtupleDColumn("pde");        // 7
    am->CreateNtupleDColumn("x_mm");       // 8
    am->CreateNtupleDColumn("y_mm");       // 9
    am->CreateNtupleDColumn("z_mm");       // 10
    am->CreateNtupleDColumn("gun_x_mm");   // 11  posicion x del muon primario [mm]
    am->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
    G4AccumulableManager::Instance()->Reset();
    auto* am = G4AnalysisManager::Instance();
    if (!am->IsOpenFile())       // <-- solo abrir si no está ya abierto
        am->OpenFile();
    G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(true);
}

void RunAction::EndOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Merge();

    auto* am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();

    const G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    G4String outFile = am->GetFileName();
    if (outFile.size() < 5 || outFile.substr(outFile.size() - 5) != ".root")
        outFile += ".root";

    G4cout
        << "\n=== EJ-200 Bar Run Summary ==="
        << "\n  Events run            : " << nEvents
        << "\n  Events with ≥1 hit    : " << fNEventsWithHits.GetValue()
        << "\n  End-left  photons     : " << fNEndLeft.GetValue()
        << "\n  End-right photons     : " << fNEndRight.GetValue()
        << "\n  Top SiPM  photons     : " << fNTop.GetValue()
        << "\n  ROOT output           : " << outFile
        << "\n==============================\n"
        << G4endl;
}
