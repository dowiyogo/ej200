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
