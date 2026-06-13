#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4MuonMinus.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"


// ---------------------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fGun(1)   // 1 primary per event
{
    fGun.SetParticleDefinition(G4MuonMinus::Definition());
    fGun.SetParticleEnergy(1.0 * GeV);
    // Default: gun 55 mm above the +Z bar face (kBarHalfZ = 5 mm),
    // straight down in Z.
    fGun.SetParticlePosition({0.0, 0.0, 60.0 * mm});
    fGun.SetParticleMomentumDirection({0.0, 0.0, -1.0});

    // ── UI messenger ─────────────────────────────────────────────────────────
    fMessenger = new G4GenericMessenger(this, "/muon/", "Muon gun control");

    // /muon/angle <deg>  — incidence angle in the XZ plane
    {
        auto& cmd = fMessenger->DeclareMethod(
            "angle",
            &PrimaryGeneratorAction::SetAngleDeg,
            "Set muon incidence angle from vertical (-Z), in degrees [XZ plane].\n"
            "  0 deg = vertical (default).  Positive = tilts toward +X.\n"
            "  Example: /muon/angle 30");
        cmd.SetParameterName("deg", false);
        cmd.SetRange("deg >= -90 && deg <= 90");
        cmd.SetDefaultValue("0");
    }

    // /muon/gunX <x> mm  — direct X override
    {
        auto& cmd = fMessenger->DeclareMethodWithUnit(
            "gunX", "mm",
            &PrimaryGeneratorAction::SetGunXmm,
            "Set gun X position directly [mm].\n"
            "  Example: /muon/gunX 0 mm");
        (void)cmd;
        cmd.SetParameterName("x", false);
    }
}

// ---------------------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fMessenger;
}

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::SetAngleDeg(G4double deg) {
    fAngleDeg = deg;
}

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::SetGunXmm(G4double xMm) {
    fGunX          = xMm;   // xMm already in G4 internal units (mm=1)
    fUseDirectGunX = true;
}

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    // ── Resolve X position ───────────────────────────────────────────────────
    const G4ThreeVector basePos = fGun.GetParticlePosition();
    G4double gunX = basePos.x();

    if (fUseDirectGunX) {
        gunX = fGunX;
    }

    // ── Resolve momentum direction from angle ────────────────────────────────
    // θ is measured from -Z axis in the XZ plane.
    // Direction: (sin θ, 0, -cos θ)
    const G4double theta = fAngleDeg * CLHEP::deg;
    const G4double sinT  = std::sin(theta);
    const G4double cosT  = std::cos(theta);

    // Preserve the current /gun/position Y and Z coordinates.
    // Only X is overridden by /muon/gunX.
    fGun.SetParticlePosition({gunX, basePos.y(), basePos.z()});
    fGun.SetParticleMomentumDirection({sinT, 0.0, -cosT});

    fGun.GeneratePrimaryVertex(event);
}
