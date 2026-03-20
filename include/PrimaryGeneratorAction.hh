#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4Event;

// --------------------------------------------------------------------------
// PrimaryGeneratorAction
//
// Default: vertical mu- (PDG ID –13) a y=+60 mm (sobre la cara +Y del bar,
// que termina en y=+30 mm), dirigido hacia abajo (0, -1, 0), energia 4 GeV.
// Con esta posicion el muon entra por la cara superior y atraviesa los 60 mm
// completos del bar antes de salir por la cara inferior.
//
// Configurable via macros:
//   /gun/particle        mu-
//   /gun/energy          4 GeV
//   /gun/position        0 60 0 mm
//   /gun/direction       0 -1 0
//
// Para el scan longitudinal (scan.mac):
//   /gun/position {xpos} 60 0 mm   ← y=60 mm siempre, x varía
//   /control/loop scan.mac xpos -650 650 70
//   (scan.mac sets /gun/position {xpos} 0 0 mm)
// --------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event*) override;

  private:
    G4ParticleGun fGun;
};
