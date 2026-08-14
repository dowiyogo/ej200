#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4MuonMinus.hh"
#include "G4SystemOfUnits.hh"

// ---------------------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fGun(1)
{
    fGun.SetParticleDefinition(G4MuonMinus::Definition());
    fGun.SetParticleEnergy(1.0 * GeV);
    // Horizontal muon entering the cylinder from +X, through the centre.
    fGun.SetParticlePosition({80.0 * mm, 0.0, 0.0});
    fGun.SetParticleMomentumDirection({-1.0, 0.0, 0.0});

    fMessenger = new G4GenericMessenger(this, "/muon/", "Muon gun control");

    // /muon/gunY <y> mm — transverse scan
    {
        auto& cmd = fMessenger->DeclareMethodWithUnit(
            "gunY", "mm", &PrimaryGeneratorAction::SetGunYmm,
            "Set gun Y position [mm] for transverse scan (default 0).");
        cmd.SetParameterName("y", false);
    }

    // /muon/gunZ <z> mm — axial scan along cylinder axis
    {
        auto& cmd = fMessenger->DeclareMethodWithUnit(
            "gunZ", "mm", &PrimaryGeneratorAction::SetGunZmm,
            "Set gun Z position [mm] for axial scan. Range: ±12.5 mm (cylinder half-height).");
        cmd.SetParameterName("z", false);
    }

    // /muon/angle <deg> — tilt from -X in the XZ plane
    {
        auto& cmd = fMessenger->DeclareMethod(
            "angle", &PrimaryGeneratorAction::SetAngleDeg,
            "Tilt angle from −X direction in XZ plane [deg]. 0=horizontal (default).");
        cmd.SetParameterName("deg", false);
        cmd.SetRange("deg >= -90 && deg <= 90");
        cmd.SetDefaultValue("0");
    }
}

// ---------------------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fMessenger; }

void PrimaryGeneratorAction::SetGunYmm   (G4double y) { fGunYmm   = y; }
void PrimaryGeneratorAction::SetGunZmm   (G4double z) { fGunZmm   = z; }
void PrimaryGeneratorAction::SetAngleDeg (G4double d) { fAngleDeg = d; }

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    // Direction: angle θ from -X axis in the XZ plane.
    //   θ=0 → (-1, 0, 0)   (horizontal, default)
    //   θ>0 → tilts toward +Z
    const G4double theta = fAngleDeg * CLHEP::deg;
    fGun.SetParticlePosition({80.0 * mm, fGunYmm, fGunZmm});
    fGun.SetParticleMomentumDirection({-std::cos(theta), 0.0, std::sin(theta)});
    fGun.GeneratePrimaryVertex(event);
}
