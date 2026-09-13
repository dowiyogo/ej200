#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4EventManager.hh"
#ifdef EJ200_ENABLE_DIAGNOSTICS
#include "TrackingAction.hh"
#endif
#ifdef EJ200_ENABLE_DIAGNOSTICS
#include "BoundaryCensus.hh"
#endif
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4TouchableHandle.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <atomic>

// hc constant for wavelength calculation [eV·nm]
static constexpr G4double kHC_eVnm = 1239.84193;  // eV·nm

// Required run observation, independent of optional diagnostic maps.
// Historical name retained: counts optical boundary encounters toward a SiPM,
// not unique photons or a proof that transmission/detection occurred.
namespace {
    std::atomic<long long> gMylarToSiPM{0};
}
namespace BoundaryCensus {
    long long GetMylarToSiPM() { return gMylarToSiPM.load(); }
    void ResetSiPMEntries() { gMylarToSiPM = 0; }
}

// Contadores de diagnóstico de frontera — acumulan durante toda la corrida.
// Se usan std::atomic para thread-safety en MT builds.
#ifdef EJ200_ENABLE_DIAGNOSTICS
namespace {
    std::atomic<long long> gBarToMylar{0};     // Bar -> reflector panel
    std::atomic<long long> gMylarToWorld{0};   // Bar -> World escape
    std::atomic<long long> gMylarReflected{0}; // Boundary reflection back to BarLV
    std::atomic<long long> gKilledWorld{0};    // Kills en WorldLV por SteppingAction
    std::atomic<long long> gSparedWorldReflection{0}; // Reflexiones que este guard no mata

    bool IsBarLV(const G4String& name) {
        return name == "BarLV";
    }

    bool IsReflectorLV(const G4String& name) {
        return name == "ReflectorYMinusLV" || name == "ReflectorXLV" ||
               name == "ReflectorZLV";
    }
}

namespace BoundaryCensus {
    long long GetBarToMylar()     { return gBarToMylar.load(); }
    long long GetMylarToWorld()   { return gMylarToWorld.load(); }
    long long GetMylarReflected() { return gMylarReflected.load(); }
    long long GetKilledWorld()    { return gKilledWorld.load(); }
    long long GetSparedWorldReflection() { return gSparedWorldReflection.load(); }
    void Reset() {
        gBarToMylar = 0;
        gMylarToWorld = 0;
        gMylarReflected = 0;
        gKilledWorld = 0;
        gSparedWorldReflection = 0;
    }
}

#endif
void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* event = dynamic_cast<EventAction*>(
        G4EventManager::GetEventManager()->GetUserEventAction());
    if (event) event->ObserveEnergy(step);
    static G4ThreadLocal G4OpBoundaryProcess* boundary_process = nullptr;
    // EXEC_27: no reutilizar el estado de una frontera anterior en otro paso.
    G4OpBoundaryProcessStatus boundary_status = Undefined;
    // Required by the escape guard, independent of optional diagnostics.
    if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::Definition() &&
        step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        if (!boundary_process) {
            auto* pv = step->GetTrack()->GetDefinition()
                           ->GetProcessManager()->GetProcessList();
            for (G4int i = 0; i < static_cast<G4int>(pv->size()); ++i) {
                if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                    boundary_process = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
                    break;
                }
            }
        }
        if (boundary_process) {
            boundary_status = boundary_process->GetStatus();
#ifdef EJ200_ENABLE_DIAGNOSTICS
            const auto* pre_pv = step->GetPreStepPoint()->GetPhysicalVolume();
            const auto* post_pv = step->GetPostStepPoint()->GetPhysicalVolume();
            if (pre_pv && post_pv) {
                BoundaryCensus::Instance().Record({
                    pre_pv->GetName(),
                    step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                    post_pv->GetName(),
                    step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                    static_cast<G4int>(boundary_status)
                });
            }
#endif
        }
    }

    auto* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;
    if (event) event->ObserveFirstEncounter(step, static_cast<G4int>(boundary_status));

    // ── Diagnóstico: contar fotones de centelleo en su primer step ───────────
    // Solo en step 1 para contar cada fotón exactamente una vez.
    // El proceso creador "Scintillation" indica origen correcto.
    if (track->GetCurrentStepNumber() == 1) {
        const G4VProcess* creator = track->GetCreatorProcess();
        if (creator != nullptr &&
            creator->GetProcessName() == "Scintillation") {
            auto* ra = dynamic_cast<RunAction*>(
                const_cast<G4UserRunAction*>(
                    G4RunManager::GetRunManager()->GetUserRunAction()));
            if (ra != nullptr) ra->AddScintPhoton();
        }
    }

    // Required SiPM observation and optional boundary diagnostics.
    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        const G4String preVolName =
            (step->GetPreStepPoint()->GetPhysicalVolume())
                ? step->GetPreStepPoint()->GetPhysicalVolume()
                      ->GetLogicalVolume()->GetName()
                : "NULL";
        const G4String postVolName =
            (step->GetPostStepPoint()->GetPhysicalVolume())
                ? step->GetPostStepPoint()->GetPhysicalVolume()
                      ->GetLogicalVolume()->GetName()
                : "NULL";

#ifdef EJ200_ENABLE_DIAGNOSTICS
        if (IsBarLV(preVolName) && IsReflectorLV(postVolName))
            ++gBarToMylar;

        if (IsBarLV(preVolName) && postVolName == "WorldLV")
            ++gMylarToWorld;

        if (IsBarLV(preVolName) && IsBarLV(postVolName))
            ++gMylarReflected;

#endif
        if (preVolName == "BarLV" &&
            (postVolName == "EndSiPMLV" || postVolName == "TopSiPMLV"))
            ++gMylarToSiPM;
    }

    // ── Wavelength filter ────────────────────────────────────────────────────
    // Kill photons outside the SiPM sensitivity window (300–900 nm).
    // EJ-200 emits in 380–500 nm so this filter has negligible effect on
    // scintillation photons, but efficiently removes any IR secondaries or
    // photons generated by Cherenkov outside the PDE range, avoiding
    // unnecessary tracking.
    const G4double energy = track->GetKineticEnergy();
    if (energy > 0.0) {
        // energy is in G4 internal units (MeV); eV = 1e-6, nm = 1 mm * 1e-6
        const G4double wl_nm = kHC_eVnm / (energy / eV);
        if (wl_nm < 300.0 || wl_nm > 900.0) {
#ifdef EJ200_ENABLE_DIAGNOSTICS
            TerminalCensus::MarkKill(track, "wavelength_filter");
#endif
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    // ── Geometry escape guard ────────────────────────────────────────────────
    auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();

    // Kill if outside the world entirely (safety net).
    if (postVol == nullptr) {
#ifdef EJ200_ENABLE_DIAGNOSTICS
        TerminalCensus::MarkKill(track, "null_post_volume");
#endif
        track->SetTrackStatus(fStopAndKill);
        return;
    }

    // Solo matar si el fotón realmente cruzó al mundo.
    // TIR y reflexión especular devuelven el fotón: no son escape.
    if (postVol->GetLogicalVolume()->GetName() == "WorldLV") {
        const auto st = boundary_status;
        const bool reflected = (st == TotalInternalReflection ||
                                st == FresnelReflection ||
                                st == LambertianReflection ||
                                st == SpikeReflection ||
                                st == BackScattering);
        if (reflected) {
#ifdef EJ200_ENABLE_DIAGNOSTICS
            ++gSparedWorldReflection;
#endif
        } else {
#ifdef EJ200_ENABLE_DIAGNOSTICS
            ++gKilledWorld;
#endif
#ifdef EJ200_ENABLE_DIAGNOSTICS
            TerminalCensus::MarkKill(track, "world_guard");
#endif
            track->SetTrackStatus(fStopAndKill);
        }
    }
}
