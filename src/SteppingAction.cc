#include "DisplayTrackingAction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <atomic>

// hc constant for wavelength calculation [eV·nm]
static constexpr G4double kHC_eVnm = 1239.84193;  // eV·nm

// Contadores de diagnóstico de frontera — acumulan durante toda la corrida.
// Se usan std::atomic para thread-safety en MT builds.
namespace {
    std::atomic<long long> gBarToMylar{0};     // Bar -> reflector panel
    std::atomic<long long> gMylarToWorld{0};   // Bar -> World escape
    std::atomic<long long> gMylarReflected{0}; // Boundary reflection back to BarLV
    std::atomic<long long> gMylarToSiPM{0};    // Bar -> SiPM detector
    std::atomic<long long> gKilledWorld{0};    // Kills en WorldLV por SteppingAction

    bool IsBarLV(const G4String& name) {
        return name == "BarLV";
    }

    bool IsReflectorLV(const G4String& name) {
        return name == "ReflectorYMinusLV" || name == "ReflectorXLV" ||
               name == "ReflectorZLV";
    }

    G4OpBoundaryProcess* GetOpBoundaryProcess() {
        static thread_local G4OpBoundaryProcess* boundary = nullptr;
        if (boundary != nullptr) return boundary;

        auto* pm = G4OpticalPhoton::Definition()->GetProcessManager();
        if (pm == nullptr) return nullptr;

        G4ProcessVector* post = pm->GetPostStepProcessVector(typeDoIt);
        if (post == nullptr) return nullptr;

        const size_t nProc = post->size();
        for (size_t i = 0; i < nProc; ++i) {
            auto* proc = (*post)[i];
            auto* asBoundary = dynamic_cast<G4OpBoundaryProcess*>(proc);
            if (asBoundary != nullptr) {
                boundary = asBoundary;
                return boundary;
            }
        }
        return nullptr;
    }
}

namespace BoundaryCensus {
    long long GetBarToMylar()     { return gBarToMylar.load(); }
    long long GetMylarToWorld()   { return gMylarToWorld.load(); }
    long long GetMylarReflected() { return gMylarReflected.load(); }
    long long GetMylarToSiPM()    { return gMylarToSiPM.load(); }
    long long GetKilledWorld()    { return gKilledWorld.load(); }
    void Reset() {
        gBarToMylar = 0;
        gMylarToWorld = 0;
        gMylarReflected = 0;
        gMylarToSiPM = 0;
        gKilledWorld = 0;
    }
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;

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

    // ── Diagnóstico de frontera ──────────────────────────────────────────────
    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        auto* boundary = GetOpBoundaryProcess();
        if (boundary != nullptr) {
            DisplayDiagnostics::RecordBoundaryStatus(
                track,
                boundary->GetStatus(),
                step->GetPostStepPoint()->GetPosition());
        }

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

        if (IsBarLV(preVolName) && IsReflectorLV(postVolName))
            ++gBarToMylar;

        if (IsBarLV(preVolName) && postVolName == "WorldLV")
            ++gMylarToWorld;

        if (IsBarLV(preVolName) && IsBarLV(postVolName))
            ++gMylarReflected;

        if (IsBarLV(preVolName) &&
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
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    // ── Geometry escape guard ────────────────────────────────────────────────
    auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();

    // Kill if outside the world entirely (safety net).
    if (postVol == nullptr) {
        track->SetTrackStatus(fStopAndKill);
        return;
    }

    // Kill photons that reach the air (world) volume. Photons entering a SiPM
    // land in EndSiPMLV or TopSiPMLV, so this correctly spares detected photons.
    if (postVol->GetLogicalVolume()->GetName() == "WorldLV") {
        ++gKilledWorld;
        track->SetTrackStatus(fStopAndKill);
    }
}
