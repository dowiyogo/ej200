#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class G4GenericMessenger;

// --------------------------------------------------------------------------
// PrimaryGeneratorAction — muon gun for EJ-228 cylinder
//
// Default: horizontal mu- through the cylinder centre from +X.
//   Position : (+80 mm, 0, 0)
//   Direction: (−1, 0, 0)        — along −X through the mantle
//   Energy   : 1 GeV             — minimum-ionizing muon
//
// ── UI commands ────────────────────────────────────────────────────────────
//  /muon/gunY <y> mm    Transverse Y position scan (default 0 mm).
//  /muon/gunZ <z> mm    Axial Z position scan along cylinder axis (default 0 mm).
//                       Range: within ±kCylHalfH = ±12.5 mm.
//  /muon/angle <deg>    Tilt from −X direction in the XZ plane.
//                       0° = horizontal; positive = tilts toward +Z.
//
// The standard /gun/particle and /gun/energy commands also work.
// --------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    void SetGunYmm   (G4double y);
    void SetGunZmm   (G4double z);
    void SetAngleDeg (G4double deg);

  private:
    G4ParticleGun       fGun;
    G4GenericMessenger* fMessenger = nullptr;

    G4double fGunYmm   = 0.0;
    G4double fGunZmm   = 0.0;
    G4double fAngleDeg = 0.0;
};
