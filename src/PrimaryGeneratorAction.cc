#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4MuonMinus.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fGun(1)  // 1 primary per event
{
    fGun.SetParticleDefinition(G4MuonMinus::Definition());
    fGun.SetParticleEnergy(1.0 * GeV);
    // Y = +60 mm: por encima de la cara superior del bar (kBarHalfY = 30 mm).
    // El muon entra por la cara +Y, atraviesa los 60 mm completos del bar
    // y sale por la cara -Y.  X = Z = 0 → centro longitudinal y transversal.
    fGun.SetParticlePosition({0.0, 60.0 * mm, 0.0});
    fGun.SetParticleMomentumDirection({0.0, -1.0, 0.0}); // vertically downward
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    fGun.GeneratePrimaryVertex(event);
}
