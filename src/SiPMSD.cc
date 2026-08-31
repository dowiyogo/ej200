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
SiPMSD::SiPMSD(const G4String& name)
    : G4VSensitiveDetector(name)
{
    SetModel(fModel);

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
    auto* post = step->GetPostStepPoint();
    if (pre->GetStepStatus() != fGeomBoundary &&
        post->GetStepStatus() != fGeomBoundary) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    // ── Identify SiPM by copy number ────────────────────────────────────────
    // BoundaryInvokeSD calls the SD while the track is still in the bar, so
    // the post-step physical volume is the authoritative SiPM placement.
    auto* pv = post->GetPhysicalVolume();
    if (pv == nullptr ||
        (pv->GetLogicalVolume()->GetName() != "EndSiPMLV" &&
         pv->GetLogicalVolume()->GetName() != "TopSiPMLV" &&
         pv->GetLogicalVolume()->GetName() != "SparseTopSiPMLV")) {
        pv = track->GetVolume();
    }
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

    const G4ThreeVector pos = post->GetPosition();
    const G4double pde = GetPDE(energy);
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
    am->FillNtupleDColumn(0, 4, time_ns);
    am->FillNtupleDColumn(0, 5, energy_eV);
    am->FillNtupleDColumn(0, 6, wl_nm);
    am->FillNtupleDColumn(0, 7, pde);
    am->FillNtupleDColumn(0, 8, pos.x() / mm);
    am->FillNtupleDColumn(0, 9, pos.y() / mm);
    am->FillNtupleDColumn(0, 10, pos.z() / mm);
    const G4double gunX = ea ? ea->GetGunXmm() : 0.0;
    am->FillNtupleDColumn(0, 11, gunX);
    am->AddNtupleRow(0);

    track->SetTrackStatus(fStopAndKill);
    return true;
}

// ---------------------------------------------------------------------------
G4double SiPMSD::GetPDE(G4double energy) const {
    if (energy <= 0.0) return 0.0;
    return SiPMModel::InterpolatePDE(fPDECurve, 1239.84193 / (energy / eV));
}

void SiPMSD::SetModel(const G4String& model) {
    const G4String canonical = SiPMModel::CanonicalName(model);
    if (canonical.empty()) return;
    fModel = canonical;
    fPDECurve = SiPMModel::LoadPDECurve(fModel);
}
