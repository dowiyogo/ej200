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
//  /muon/midpointSiPMs <i> <j>
//      Place the gun X at the geometric midpoint between top SiPMs i and j
//      (0-based indices, e.g. "9 10" for the central gap).
//      Reads the current pitch from DetectorConstruction at event time, so
//      it remains valid after /det/topSiPMPitch changes.
//
//  /muon/gunX <x> mm
//      Set gun X position directly in mm. This clears midpoint mode.
//
// The standard /gun/* commands still work for particle type, energy, and
// position. If no /muon/gunX or /muon/midpointSiPMs override is active, the
// current /gun/position X coordinate is preserved.
// --------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    // Called by messenger commands
    void SetAngleDeg      (G4double deg);
    void SetMidpointSiPMs (G4String indices);  // format: "i j"
    void SetGunXmm        (G4double xMm);

  private:
    G4ParticleGun       fGun;
    G4GenericMessenger* fMessenger   = nullptr;

    G4double fAngleDeg      = 0.0;    // incidence angle [degrees]
    G4double fGunX          = 0.0;    // direct X override [G4 internal = mm]
    G4bool   fUseDirectGunX = false;  // true after /muon/gunX
    G4int    fMidSiPM1      = -1;     // >= 0 → use midpoint mode
    G4int    fMidSiPM2   = -1;
};
