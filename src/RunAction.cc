#include "RunAction.hh"
#include "EventAction.hh"
#ifdef EJ200_ENABLE_DIAGNOSTICS
#include "TrackingAction.hh"
#endif
#ifdef EJ200_ENABLE_DIAGNOSTICS
#include "BoundaryCensus.hh"
#endif
#include "Randomize.hh"
#include "G4Threading.hh"
#include "DetectorConstruction.hh"
#include "SiPMModel.hh"
#include "SteppingAction.hh"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicsVector.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

namespace {
void LogActiveScintillator() {
    auto* detector = dynamic_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detector == nullptr) return;

    const auto* material = detector->GetActiveScintillatorMaterial();
    const auto* mpt = material ? material->GetMaterialPropertiesTable() : nullptr;
    if (mpt == nullptr) return;

    const auto* attenuation = mpt->GetProperty("ABSLENGTH");
    const auto* emission = mpt->GetProperty("SCINTILLATIONCOMPONENT1");
    G4double peakWavelength = 0.0;
    if (emission != nullptr && emission->GetVectorLength() > 0) {
        std::size_t peakIndex = 0;
        for (std::size_t i = 1; i < emission->GetVectorLength(); ++i) {
            if ((*emission)[i] > (*emission)[peakIndex]) peakIndex = i;
        }
        peakWavelength = 1239.84193 * eV * nm / emission->Energy(peakIndex);
    }

    G4cout
        << "\n=== Active Scintillator Baseline ==="
        << "\n  SSLG4 code            : " << detector->GetScintillatorCode()
        << "\n  Yield                 : "
        << mpt->GetConstProperty("SCINTILLATIONYIELD") * MeV << " ph/MeV"
        << "\n  Rise time tau_r       : "
        << mpt->GetConstProperty("SCINTILLATIONRISETIME1") / ns << " ns"
        << "\n  Decay time tau_d      : "
        << mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1") / ns << " ns"
        << "\n  Attenuation length    : "
        << (attenuation ? (*attenuation)[0] / cm : 0.0) << " cm"
        << "\n  Emission peak lambda  : " << peakWavelength / nm << " nm"
        << "\n  Finite rise time      : "
        << (G4OpticalParameters::Instance()->GetScintFiniteRiseTime() ? "enabled" : "DISABLED")
        << "\n====================================="
        << G4endl;
}

void LogReadoutConfiguration() {
    auto* detector = dynamic_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detector == nullptr) return;

    // --- Resolve BarSkin reflector (replaces old per-face border surfaces) ---
    G4double reflectivity = 0.0;
    G4String skinLabel = "OPEN/UNDEFINED";
    const auto& sipmSurfaces = detector->GetSiPMSurfaces();
    if (!sipmSurfaces.empty()) {
        auto* barPV = sipmSurfaces.begin()->second->GetVolume1();
        auto* barLV = barPV ? barPV->GetLogicalVolume() : nullptr;
        auto* skin  = barLV ? G4LogicalSkinSurface::GetSurface(barLV) : nullptr;
        if (skin != nullptr) {
            skinLabel = "reflective/BarSkin";
            auto* optical = dynamic_cast<G4OpticalSurface*>(skin->GetSurfaceProperty());
            auto* mpt  = optical ? optical->GetMaterialPropertiesTable() : nullptr;
            auto* prop = mpt ? mpt->GetProperty("REFLECTIVITY") : nullptr;
            if (prop != nullptr && prop->GetVectorLength() > 0)
                reflectivity = (*prop)[0];
        }
    }

    const auto faceState = [&skinLabel](G4bool instrumented) -> G4String {
        return instrumented ? G4String("instrumented") : skinLabel;
    };

    G4cout
        << "\n=== Active Readout / Wrapping Configuration ==="
        << "\n  Readout configuration : " << detector->GetReadoutConfiguration()
        << "\n  SiPM model            : " << detector->GetSiPMModel()
        << "\n  SiPM PDE file         : " << SiPMModel::DataFilePath(detector->GetSiPMModel())
        << "\n  -X face               : " << faceState(detector->IsEndInstrumented())
        << "\n  +X face               : " << faceState(detector->IsEndInstrumented())
        << "\n  -Y face               : " << faceState(false)
        << "\n  +Y face               : " << faceState(detector->IsTopInstrumented())
        << "\n  -Z face               : " << faceState(false)
        << "\n  +Z face               : " << faceState(false)
        << "\n  Reflector R           : " << reflectivity
        << "\n  Active End SiPMs      : " << detector->GetNActiveEndSiPMs()
        << " (L=" << detector->GetNActiveEndSiPMs() / 2
        << ", R=" << detector->GetNActiveEndSiPMs() / 2 << ")"
        << "\n  Active Top SiPMs      : " << detector->GetNActiveTopSiPMs()
        << "\n=============================================="
        << G4endl;
}
} // namespace

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
    am->CreateNtupleIColumn("track_id");
    am->FinishNtuple();
    EventAction::BookEnergyObservations();
    EventAction::BookFirstEncounters();
    EventAction::BookSiPMObservations();
}

void RunAction::BeginOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Reset();
    if (IsMaster()) BoundaryCensus::ResetSiPMEntries();
#ifdef EJ200_ENABLE_DIAGNOSTICS
    // EXEC_26: reinicio único y copia del motor, sin consumir números aleatorios.
    if (IsMaster()) {
        BoundaryCensus::Reset(); // EXEC_30: shared atomics reset once, before workers.
        BoundaryCensus::Instance().Reset();
        TerminalCensus::Reset();
    }
    const auto rngPrefix = "rng_run" + std::to_string(run->GetRunID()) +
                           "_thread" + std::to_string(G4Threading::G4GetThreadId());
    G4Random::saveEngineStatus((rngPrefix + "_begin.rndm").c_str());
    G4cout << "EXEC_26 RNG engine: " << G4Random::getTheEngine()->name()
           << "; state: " << rngPrefix << "_begin.rndm" << G4endl;

#endif
    auto* am = G4AnalysisManager::Instance();

    std::ostringstream fname;
    fname << "photon_hits_run"
          << std::setw(3) << std::setfill('0')
          << run->GetRunID();

    am->SetFileName(fname.str());
    am->OpenFile();

    G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(true);
    if (IsMaster()) {
        LogActiveScintillator();
        LogReadoutConfiguration();
    }
}

void RunAction::EndOfRunAction(const G4Run* run) {
    G4AccumulableManager::Instance()->Merge();

    auto* am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();

#ifdef EJ200_ENABLE_DIAGNOSTICS
    // EXEC_26: el maestro exporta tras finalizar los trabajadores; cada hilo guarda su motor.
    const auto rngPrefix = "rng_run" + std::to_string(run->GetRunID()) +
                           "_thread" + std::to_string(G4Threading::G4GetThreadId());
    G4Random::saveEngineStatus((rngPrefix + "_end.rndm").c_str());
    if (IsMaster()) {
        TerminalCensus::Write("terminal_fates_run" + std::to_string(run->GetRunID()));
        BoundaryCensus::Instance().Write(
            "boundary_census_run" + std::to_string(run->GetRunID()) + ".csv");
    }

#endif
    if (!IsMaster()) return; // EXEC_30: aggregate summary only after worker merges.
    const G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    G4String outFile = am->GetFileName();
    if (outFile.size() < 5 || outFile.substr(outFile.size() - 5) != ".root")
        outFile += ".root";

    const G4int nSc  = fNScintPhotons.GetValue();
    const G4int nDet = fNEndLeft.GetValue() + fNEndRight.GetValue() + fNTop.GetValue();
    const G4double eff = (nSc > 0) ? 100.0 * nDet / nSc : 0.0;

    G4cout
        << "\n=== EJ Scintillator Bar Run Summary ==="
        << "\n  Run ID                : " << run->GetRunID()
        << "\n  Events run            : " << nEvents
        << "\n  Events with ≥1 hit    : " << fNEventsWithHits.GetValue()
        << "\n  End-left  photons     : " << fNEndLeft.GetValue()
        << "\n  End-right photons     : " << fNEndRight.GetValue()
        << "\n  Top SiPM  photons     : " << fNTop.GetValue()
        << "\n  Scint photons generated: " << nSc
        << "\n  Bar -> SiPM (entering)   : " << BoundaryCensus::GetMylarToSiPM()
        << "\n  Total photons detected : " << nDet
        << "\n  Detection efficiency   : " << std::fixed << std::setprecision(4)
        << eff << " %"
        << "\n  ROOT output           : " << outFile
        << "\n==============================\n"
        << G4endl;

#ifdef EJ200_ENABLE_DIAGNOSTICS
    if (IsMaster()) {
        G4cout
            << "\n=== Boundary Census (diagnostic) ==="
            << "\n  Bar -> reflector panel   : " << BoundaryCensus::GetBarToMylar()
            << "\n  Bar -> World (escaped)   : " << BoundaryCensus::GetMylarToWorld()
            << "\n  Bar -> Bar (TIR/refl)    : " << BoundaryCensus::GetMylarReflected()
            << "\n  Killed in WorldLV        : " << BoundaryCensus::GetKilledWorld()
            << "\n  Spared World reflection : " << BoundaryCensus::GetSparedWorldReflection()
            << "\n====================================\n"
            << G4endl;
    }
#endif
}
