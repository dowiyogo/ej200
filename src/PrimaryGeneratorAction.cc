#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4MuonMinus.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <sstream>
#include <stdexcept>

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

    // /muon/midpointSiPMs <i> <j>  — gun X at midpoint of two top SiPMs
    {
        fMessenger->DeclareMethod(
            "midpointSiPMs",
            &PrimaryGeneratorAction::SetMidpointSiPMs,
            "Place gun X at midpoint of two top SiPMs (0-based indices).\n"
            "  The pitch is read from DetectorConstruction at run time.\n"
            "  Example: /muon/midpointSiPMs 9 10");
    }

    // /muon/gunX <x> mm  — direct X override
    {
        auto& cmd = fMessenger->DeclarePropertyWithUnit(
            "gunX", "mm", fGunX,
            "Set gun X position directly [mm].\n"
            "  Clears any midpointSiPMs setting.\n"
            "  Example: /muon/gunX 0 mm");
        (void)cmd;
        // Changing fGunX clears the midpoint mode
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
void PrimaryGeneratorAction::SetMidpointSiPMs(G4String indices) {
    // Parse "i j" from the macro argument string
    std::istringstream ss(indices);
    G4int a = -1, b = -1;
    if (!(ss >> a >> b) || a < 0 || b < 0) {
        G4cerr << "[PrimaryGeneratorAction] /muon/midpointSiPMs: "
               << "invalid argument \"" << indices
               << "\". Expected two non-negative integers, e.g. \"9 10\".\n";
        return;
    }
    fMidSiPM1 = a;
    fMidSiPM2 = b;
}

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::SetGunXmm(G4double xMm) {
    fGunX    = xMm;   // xMm already in G4 internal units (mm=1) from DeclarePropertyWithUnit
    fMidSiPM1 = -1;   // clear midpoint mode
    fMidSiPM2 = -1;
}

// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    // ── Resolve X position ───────────────────────────────────────────────────
    G4double gunX = fGunX;

    if (fMidSiPM1 >= 0 && fMidSiPM2 >= 0) {
        // Look up the current pitch and count from DetectorConstruction
        const auto* dc = dynamic_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());

        if (dc != nullptr) {
            const G4double pitch  = dc->GetTopSiPMPitch();
            const G4int    nTotal = dc->GetNTopSiPMs();
            const G4double x1 = DetectorConstruction::TopSiPMCenterX(
                                     fMidSiPM1, pitch, nTotal);
            const G4double x2 = DetectorConstruction::TopSiPMCenterX(
                                     fMidSiPM2, pitch, nTotal);
            gunX = 0.5 * (x1 + x2);
        } else {
            G4cerr << "[PrimaryGeneratorAction] Cannot resolve DetectorConstruction "
                      "for midpoint calculation; using fGunX.\n";
        }
    }

    // ── Resolve momentum direction from angle ────────────────────────────────
    // θ is measured from -Z axis in the XZ plane.
    // Direction: (sin θ, 0, -cos θ)
    const G4double theta = fAngleDeg * CLHEP::deg;
    const G4double sinT  = std::sin(theta);
    const G4double cosT  = std::cos(theta);

    // Gun position: above the bar (+Z face is at kBarHalfZ = 5 mm).
    // Keep Y=0 (centre of bar width) and place the gun 55 mm above the top face.
    // X is set to gunX; Z stays at +60 mm regardless of angle so the muon
    // always enters through the top face even at large angles.
    fGun.SetParticlePosition({gunX, 0.0, 60.0 * mm});
    fGun.SetParticleMomentumDirection({sinT, 0.0, -cosT});

    fGun.GeneratePrimaryVertex(event);
}
