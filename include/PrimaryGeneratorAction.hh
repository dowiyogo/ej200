#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class G4GenericMessenger;

// --------------------------------------------------------------------------
// PrimaryGeneratorAction
//
// Default: vertical mu- through the bar centre (X=0, angle=0°).
//   Position : (0, 0, +60 mm) — 55 mm above the +Z bar face
//   Direction: (0, 0, -1)     — straight down
//   Energy   : 1 GeV
//
// ── New UI commands (available after /run/initialize) ──────────────────────
//
//  /muon/angle <deg>
//      Incidence angle θ from the -Z vertical axis, in the XZ plane.
//      Direction becomes (sinθ, 0, -cosθ).
//      Positive θ tilts toward +X.  Example: /muon/angle 30 deg
//
//  /muon/gunX <x> mm
//      Set gun X position directly in mm.
//
// The standard /gun/* commands still work for particle type, energy, and
// position. If no /muon/gunX override is active, the
// current /gun/position X coordinate is preserved.
// --------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    // Called by messenger commands
    void SetAngleDeg      (G4double deg);
    void SetGunXmm        (G4double xMm);

  private:
    G4ParticleGun       fGun;
    G4GenericMessenger* fMessenger   = nullptr;

    G4double fAngleDeg      = 0.0;    // incidence angle [degrees]
    G4double fGunX          = 0.0;    // direct X override [G4 internal = mm]
    G4bool   fUseDirectGunX = false;  // true after /muon/gunX
};
