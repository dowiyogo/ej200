#include "RunAction.hh"
#include "SteppingAction.hh"

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
        accMgr->Register(fNScintPhotons);
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
    am->CreateNtupleDColumn("energy_eV");
    am->CreateNtupleDColumn("wl_nm");
    am->CreateNtupleDColumn("pde");
    am->CreateNtupleDColumn("x_mm");
    am->CreateNtupleDColumn("y_mm");
    am->CreateNtupleDColumn("z_mm");
    am->CreateNtupleDColumn("gun_x_mm");
    am->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Reset();
    BoundaryCensus::Reset();

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

    const G4int nSc  = fNScintPhotons.GetValue();
    const G4int nDet = fNEndLeft.GetValue() + fNEndRight.GetValue() + fNTop.GetValue();
    const G4double eff = (nSc > 0) ? 100.0 * nDet / nSc : 0.0;

    G4cout
        << "\n=== EJ-200 Bar Run Summary ==="
        << "\n  Run ID                : " << run->GetRunID()
        << "\n  Events run            : " << nEvents
        << "\n  Events with ≥1 hit    : " << fNEventsWithHits.GetValue()
        << "\n  End-left  photons     : " << fNEndLeft.GetValue()
        << "\n  End-right photons     : " << fNEndRight.GetValue()
        << "\n  Top SiPM  photons     : " << fNTop.GetValue()
        << "\n  Scint photons generated: " << nSc
        << "\n  Total photons detected : " << nDet
        << "\n  Detection efficiency   : " << std::fixed << std::setprecision(4)
        << eff << " %"
        << "\n  ROOT output           : " << outFile
        << "\n==============================\n"
        << G4endl;

    if (IsMaster()) {
        G4cout
            << "\n=== Boundary Census (diagnostic) ==="
            << "\n  Bar -> reflector panel   : " << BoundaryCensus::GetBarToMylar()
            << "\n  Bar -> World (escaped)   : " << BoundaryCensus::GetMylarToWorld()
            << "\n  Bar -> Bar (TIR/refl)    : " << BoundaryCensus::GetMylarReflected()
            << "\n  Bar -> SiPM (entering)   : " << BoundaryCensus::GetMylarToSiPM()
            << "\n  Killed in WorldLV        : " << BoundaryCensus::GetKilledWorld()
            << "\n====================================\n"
            << G4endl;
    }
}
