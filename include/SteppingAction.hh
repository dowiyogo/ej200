#pragma once
#include "G4UserSteppingAction.hh"

// --------------------------------------------------------------------------
// SteppingAction — optional per-step diagnostics.
//
// Currently minimal: kills optical photons that escape the world volume
// (should not happen with wrapping, but acts as a safety net).
// Can be extended to add per-step scoring or debugging output.
// --------------------------------------------------------------------------
class SteppingAction : public G4UserSteppingAction {
  public:
    SteppingAction()  = default;
    ~SteppingAction() = default;

    void UserSteppingAction(const G4Step*) override;
};

// Acceso externo a los contadores de diagnóstico de frontera.
// Definidos en SteppingAction.cc; thread-safe (std::atomic).
// Historical accessor names are kept for compatibility with diag scripts.
namespace BoundaryCensus {
    long long GetBarToMylar();
    long long GetMylarToWorld();
    long long GetMylarReflected();
    long long GetMylarToSiPM();
    long long GetKilledWorld();
    void      Reset();
}

// Photon-budget tallies for diagnostic runs. Crossings and terminal fates are
// reported separately because a photon can cross/reflect many times.
namespace PhotonBudget {
    long long GetGenerated();
    long long GetGeneratedScintillation();
    long long GetGeneratedCherenkov();
    long long GetReachedEndFace();
    long long GetEnteredEndSiPM();
    long long GetEnteredTopSiPM();
    long long GetDetectedEnd();
    long long GetDetectedTop();
    long long GetRejectedPDEEnd();
    long long GetRejectedPDETop();
    long long GetBulkAbsorption();
    long long GetSurfaceAbsorption();
    long long GetEscapedWorld();
    long long GetWavelengthKilled();
    long long GetTotalInternalReflection();
    long long GetBoundaryReflection();
    void AddDetectedEnd();
    void AddDetectedTop();
    void AddRejectedPDEEnd();
    void AddRejectedPDETop();
    void Reset();
}
