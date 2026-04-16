#include "RunAction.hh"
#include <sstream>
#include <iomanip>
#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4OpticalParameters.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

RunAction::RunAction() {
    if (IsMaster()) {
        auto* accMgr = G4AccumulableManager::Instance();
        accMgr->Register(fNEndLeft);
        accMgr->Register(fNEndRight);
        accMgr->Register(fNTop);
        accMgr->Register(fNEventsWithHits);
    }

    auto* am = G4AnalysisManager::Instance();
    am->SetVerboseLevel(0);
    am->SetDefaultFileType("root");
    am->SetNtupleMerging(true);

    am->CreateNtuple("sipm_hits", "Detected optical photons in all SiPMs");
    am->CreateNtupleIColumn("event_id");
    am->CreateNtupleIColumn("face_type");
    am->CreateNtupleIColumn("global_id");
    am->CreateNtupleIColumn("local_id");
    am->CreateNtupleDColumn("time_ns");
    am->CreateNtupleDColumn("time_raw_ns");
    am->CreateNtupleDColumn("sptr_jitter_ns");
    am->CreateNtupleDColumn("energy_eV");
    am->CreateNtupleDColumn("wl_nm");
    am->CreateNtupleDColumn("pde");
    am->CreateNtupleDColumn("x_mm");
    am->CreateNtupleDColumn("y_mm");
    am->CreateNtupleDColumn("z_mm");
    am->CreateNtupleDColumn("gun_x_mm");
    am->CreateNtupleIColumn("sensor_model_id");
    am->CreateNtupleDColumn("sensor_sptr_sigma_ns");
    am->CreateNtupleDColumn("electronics_sigma_ns");
    am->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Reset();

    auto* am = G4AnalysisManager::Instance();

    std::ostringstream fname;
    fname << "photon_hits_run"
          << std::setw(3) << std::setfill('0')
          << run->GetRunID();

    am->SetFileName(fname.str());
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
        << "\n  Run ID                : " << run->GetRunID()
        << "\n  Events run            : " << nEvents
        << "\n  Events with ≥1 hit    : " << fNEventsWithHits.GetValue()
        << "\n  End-left  photons     : " << fNEndLeft.GetValue()
        << "\n  End-right photons     : " << fNEndRight.GetValue()
        << "\n  Top SiPM  photons     : " << fNTop.GetValue()
        << "\n  ROOT output           : " << outFile
        << "\n==============================\n"
        << G4endl;
}
