#pragma once
#include "G4UserSteppingAction.hh"

// Optical wavelength and escape guards, scintillation scoring, and optional census.
// Boundary status lookup remains active with diagnostics OFF to preserve TIR.
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
    // Required physics smoke observation, available with diagnostics OFF.
    long long GetMylarToSiPM();
    void ResetSiPMEntries();
}

#ifdef EJ200_ENABLE_DIAGNOSTICS
namespace BoundaryCensus {
    long long GetBarToMylar();
    long long GetMylarToWorld();
    long long GetMylarReflected();
    long long GetKilledWorld();
    long long GetSparedWorldReflection();
    void      Reset();
}

#endif
