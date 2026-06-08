#include "SteppingAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <atomic>
#include <cmath>
#include <unordered_set>

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

    std::atomic<long long> gGenerated{0};
    std::atomic<long long> gGeneratedScintillation{0};
    std::atomic<long long> gGeneratedCherenkov{0};
    std::atomic<long long> gReachedEndFace{0};
    std::atomic<long long> gEnteredEndSiPM{0};
    std::atomic<long long> gEnteredTopSiPM{0};
    std::atomic<long long> gDetectedEnd{0};
    std::atomic<long long> gDetectedTop{0};
    std::atomic<long long> gRejectedPDEEnd{0};
    std::atomic<long long> gRejectedPDETop{0};
    std::atomic<long long> gBulkAbsorption{0};
    std::atomic<long long> gSurfaceAbsorption{0};
    std::atomic<long long> gEscapedWorld{0};
    std::atomic<long long> gWavelengthKilled{0};
    std::atomic<long long> gTotalInternalReflection{0};
    std::atomic<long long> gBoundaryReflection{0};

    thread_local std::unordered_set<unsigned long long> gEndFaceSeen;

    bool IsBarLV(const G4String& name) {
        return name == "BarLV";
    }

    bool IsReflectorLV(const G4String& name) {
        return name == "ReflectorYMinusLV" || name == "ReflectorXLV" ||
               name == "ReflectorZLV";
    }

    G4OpBoundaryProcess* BoundaryProcess() {
        thread_local G4OpBoundaryProcess* process = nullptr;
        if (process != nullptr) return process;
        auto* manager = G4OpticalPhoton::Definition()->GetProcessManager();
        if (manager == nullptr) return nullptr;
        auto* processes = manager->GetPostStepProcessVector(typeDoIt);
        for (G4int i = 0; i < processes->entries(); ++i) {
            process = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]);
            if (process != nullptr) return process;
        }
        return nullptr;
    }

    unsigned long long TrackKey(const G4Track* track) {
        const auto* run = G4RunManager::GetRunManager()->GetCurrentRun();
        const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
        const unsigned long long runId = run ? static_cast<unsigned>(run->GetRunID()) : 0;
        const unsigned long long eventId = event ? static_cast<unsigned>(event->GetEventID()) : 0;
        return (runId << 48) ^ (eventId << 24) ^ static_cast<unsigned>(track->GetTrackID());
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

namespace PhotonBudget {
    long long GetGenerated()               { return gGenerated.load(); }
    long long GetGeneratedScintillation()  { return gGeneratedScintillation.load(); }
    long long GetGeneratedCherenkov()      { return gGeneratedCherenkov.load(); }
    long long GetReachedEndFace()          { return gReachedEndFace.load(); }
    long long GetEnteredEndSiPM()          { return gEnteredEndSiPM.load(); }
    long long GetEnteredTopSiPM()          { return gEnteredTopSiPM.load(); }
    long long GetDetectedEnd()             { return gDetectedEnd.load(); }
    long long GetDetectedTop()             { return gDetectedTop.load(); }
    long long GetRejectedPDEEnd()          { return gRejectedPDEEnd.load(); }
    long long GetRejectedPDETop()          { return gRejectedPDETop.load(); }
    long long GetBulkAbsorption()          { return gBulkAbsorption.load(); }
    long long GetSurfaceAbsorption()       { return gSurfaceAbsorption.load(); }
    long long GetEscapedWorld()            { return gEscapedWorld.load(); }
    long long GetWavelengthKilled()        { return gWavelengthKilled.load(); }
    long long GetTotalInternalReflection() { return gTotalInternalReflection.load(); }
    long long GetBoundaryReflection()      { return gBoundaryReflection.load(); }
    void AddDetectedEnd()                  { ++gDetectedEnd; }
    void AddDetectedTop()                  { ++gDetectedTop; }
    void AddRejectedPDEEnd()               { ++gRejectedPDEEnd; }
    void AddRejectedPDETop()               { ++gRejectedPDETop; }
    void Reset() {
        gGenerated = 0;
        gGeneratedScintillation = 0;
        gGeneratedCherenkov = 0;
        gReachedEndFace = 0;
        gEnteredEndSiPM = 0;
        gEnteredTopSiPM = 0;
        gDetectedEnd = 0;
        gDetectedTop = 0;
        gRejectedPDEEnd = 0;
        gRejectedPDETop = 0;
        gBulkAbsorption = 0;
        gSurfaceAbsorption = 0;
        gEscapedWorld = 0;
        gWavelengthKilled = 0;
        gTotalInternalReflection = 0;
        gBoundaryReflection = 0;
    }
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;

    // Count every optical photon exactly once and retain its production mode.
    if (track->GetCurrentStepNumber() == 1) {
        ++gGenerated;
        const G4VProcess* creator = track->GetCreatorProcess();
        const G4String creatorName = creator ? creator->GetProcessName() : "";
        if (creatorName == "Scintillation") {
            auto* ra = dynamic_cast<RunAction*>(
                const_cast<G4UserRunAction*>(
                    G4RunManager::GetRunManager()->GetUserRunAction()));
            if (ra != nullptr) ra->AddScintPhoton();
            ++gGeneratedScintillation;
        } else if (creatorName == "Cerenkov") {
            ++gGeneratedCherenkov;
        }
    }

    const auto* postProcess = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (postProcess != nullptr && postProcess->GetProcessName() == "OpAbsorption") {
        ++gBulkAbsorption;
    }

    // ── Diagnóstico de frontera ──────────────────────────────────────────────
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

        if (IsBarLV(preVolName) && IsReflectorLV(postVolName))
            ++gBarToMylar;

        if (IsBarLV(preVolName) && postVolName == "WorldLV")
            ++gMylarToWorld;

        if (IsBarLV(preVolName) && IsBarLV(postVolName))
            ++gMylarReflected;

        if (IsBarLV(preVolName) &&
            (postVolName == "EndSiPMLV" || postVolName == "TopSiPMLV"))
            ++gMylarToSiPM;

        if (IsBarLV(preVolName) &&
            (postVolName == "EndSiPMLV" ||
             std::abs(step->GetPostStepPoint()->GetPosition().x()) >= 699.999 * mm) &&
            gEndFaceSeen.insert(TrackKey(track)).second) {
            ++gReachedEndFace;
        }
        if (IsBarLV(preVolName) && postVolName == "EndSiPMLV") ++gEnteredEndSiPM;
        if (IsBarLV(preVolName) && postVolName == "TopSiPMLV") ++gEnteredTopSiPM;

        auto* boundary = BoundaryProcess();
        if (boundary != nullptr) {
            const auto status = boundary->GetStatus();
            if (status == Absorption) ++gSurfaceAbsorption;
            if (status == TotalInternalReflection) ++gTotalInternalReflection;
            if (status == FresnelReflection || status == LambertianReflection ||
                status == LobeReflection || status == SpikeReflection ||
                status == BackScattering) {
                ++gBoundaryReflection;
            }
        }
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
            ++gWavelengthKilled;
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
        ++gEscapedWorld;
        track->SetTrackStatus(fStopAndKill);
    }
}
