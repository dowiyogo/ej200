#include "SiPMSD.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

// ---------------------------------------------------------------------------
// PDE table — Hamamatsu S13360-6025 (25 µm cell, 6×6 mm²)
// 33 points, 300–940 nm.  Source: Hamamatsu S13360 series datasheet.
// ---------------------------------------------------------------------------
namespace {
    constexpr G4int    kNPDE   = 33;
    constexpr G4double kWL_nm[kNPDE] = {
        300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500,
        520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940
    };
    constexpr G4double kPDE[kNPDE] = {
        0.000, 0.050, 0.120, 0.180, 0.260, 0.330, 0.380, 0.400,
        0.405, 0.403, 0.390, 0.370, 0.340, 0.310, 0.280, 0.240,
        0.210, 0.180, 0.150, 0.120, 0.100, 0.080, 0.060, 0.050,
        0.040, 0.030, 0.025, 0.020, 0.015, 0.010, 0.008, 0.004, 0.001
    };
} // namespace

// ---------------------------------------------------------------------------
SiPMSD::SiPMSD(const G4String& name)
    : G4VSensitiveDetector(name)
{
    // Build PDE curve in ascending energy order
    const G4double hc = 1239.84193 * eV * nm;
    for (G4int i = 0; i < kNPDE; ++i)
        fPDEVec.InsertValues(hc / (kWL_nm[i] * nm), kPDE[i]);

    // ── UI messenger ─────────────────────────────────────────────────────────
    fMessenger = new G4GenericMessenger(this, "/sipm/", "SiPM detector control");

    auto& cmd = fMessenger->DeclareMethodWithUnit(
        "jitterSigma", "ns",
        &SiPMSD::SetJitterSigma,
        "Set electronic time-jitter sigma [ns].\n"
        "  Default: 0.020 ns (= 20 ps).\n"
        "  Example: /sipm/jitterSigma 0.050 ns");
    cmd.SetParameterName("sigma", false);
    cmd.SetRange("sigma >= 0");
}

SiPMSD::~SiPMSD() {
    delete fMessenger;
}

// ---------------------------------------------------------------------------
G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    auto* track = step->GetTrack();

    // Only optical photons entering from a geometry boundary
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

    auto* pre = step->GetPreStepPoint();
    if (pre->GetStepStatus() != fGeomBoundary) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    // ── Identify SiPM by copy number ────────────────────────────────────────
    // G4Track::GetVolume() always returns the current volume of the track,
    // avoiding the pre-step touchable ambiguity at geometry boundaries.
    auto* pv = track->GetVolume();
    if (pv == nullptr) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }
    const G4int globalId = pv->GetCopyNo();

    // ── Photon kinematics ────────────────────────────────────────────────────
    const G4double energy    = pre->GetKineticEnergy();
    const G4double energy_eV = energy / eV;
    const G4double wl_nm     = (energy_eV > 0.0) ? (1239.84193 / energy_eV) : 0.0;

    // ── Electronic time jitter ───────────────────────────────────────────────
    // Simulate the timing resolution of the readout electronics by smearing
    // the photon arrival time with a Gaussian of zero mean and sigma = fJitterSigma.
    // G4RandGauss::shoot(mean, sigma) draws from the CLHEP Gaussian RNG.
    const G4double jitter  = G4RandGauss::shoot(0.0, fJitterSigma);
    const G4double time_ns = (track->GetGlobalTime() + jitter) / ns;

    const G4ThreeVector pos = pre->GetPosition();

    // ── PDE Bernoulli trial ──────────────────────────────────────────────────
    const G4double pde      = GetPDE(energy);
    const G4bool   detected = (G4UniformRand() < pde);

    if (detected) {
        const G4int eventId =
            G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

        auto* ea = dynamic_cast<EventAction*>(
            G4EventManager::GetEventManager()->GetUserEventAction());
        if (ea != nullptr) {
            const G4int face = DetectorConstruction::FaceType(globalId);
            if      (face == 0) ea->AddEndLeftHit();
            else if (face == 1) ea->AddEndRightHit();
            else                ea->AddTopHit();
        }

        auto* am = G4AnalysisManager::Instance();
        am->FillNtupleIColumn(0, 0, eventId);
        am->FillNtupleIColumn(0, 1, DetectorConstruction::FaceType(globalId));
        am->FillNtupleIColumn(0, 2, globalId);
        am->FillNtupleIColumn(0, 3, DetectorConstruction::LocalId(globalId));
        am->FillNtupleDColumn(0, 4, time_ns);          // ← jitter-smeared time
        am->FillNtupleDColumn(0, 5, energy_eV);
        am->FillNtupleDColumn(0, 6, wl_nm);
        am->FillNtupleDColumn(0, 7, pde);
        am->FillNtupleDColumn(0, 8, pos.x() / mm);
        am->FillNtupleDColumn(0, 9, pos.y() / mm);
        am->FillNtupleDColumn(0, 10, pos.z() / mm);
        const G4double gunX = ea ? ea->GetGunXmm() : 0.0;
        am->FillNtupleDColumn(0, 11, gunX);
        am->AddNtupleRow(0);
    }

    track->SetTrackStatus(fStopAndKill);
    return detected;
}

// ---------------------------------------------------------------------------
G4double SiPMSD::GetPDE(G4double energy) const {
    return fPDEVec.Value(energy);
}
