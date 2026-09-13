#include "EventAction.hh"
#include "PhysicalObservation.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

namespace {
G4int VolumeId(const G4VPhysicalVolume* volume) {
    if (!volume) return 0;
    const auto& name = volume->GetName();
    if (name == "BarPV") return 1;
    if (name == "AirGapYMinusPV") return 2;
    if (name == "AirGapZPlusPV") return 3;
    if (name == "AirGapZMinusPV") return 4;
    if (name == "EndSiPMLeft_PV") return 5;
    if (name == "EndSiPMRight_PV") return 6;
    if (name == "TopSiPMPV") return 7;
    if (name == "WorldPV") return 8;
    return 0;
}
}

void EventAction::BookFirstEncounters() {
    auto* am = G4AnalysisManager::Instance();
    am->CreateNtuple("first_bar_encounters", "One first physical bar encounter per optical track");
    for (const auto* name : {"event_id", "track_id", "source", "pre_copy", "post_copy",
                            "pre_volume_id", "post_volume_id", "outcome", "exiting_bar"})
        am->CreateNtupleIColumn(name);
    am->CreateNtupleFColumn("cos_incidence");
#ifdef EJ200_ENABLE_DIAGNOSTICS
    for (const auto* name : {"boundary_status", "normal_valid",
                            "normal_orientation_valid"}) am->CreateNtupleIColumn(name);
    for (const auto* name : {"normal_norm", "energy_eV"}) am->CreateNtupleDColumn(name);
#endif
    am->FinishNtuple();
}

void EventAction::ObserveFirstEncounter(const G4Step* step, G4int status) {
    if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary ||
        status == Undefined || status == NotAtBoundary || status == StepTooSmall) return;
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    const auto* prePV = pre->GetPhysicalVolume();
    const auto* postPV = post->GetPhysicalVolume();
    if (!prePV || !postPV || prePV == postPV) return;
    const bool exiting = prePV->GetLogicalVolume()->GetName() == "BarLV";
    const bool entering = postPV->GetLogicalVolume()->GetName() == "BarLV";
    if (!exiting && !entering) return;
    if (!fFirstEncounterTracks.insert(step->GetTrack()->GetTrackID()).second) return;

    // Equivalent solid normal, oriented from pre to post without looking at the
    // momentum sign. For entry into a daughter, negate that daughter's outward
    // normal. For exit from a daughter or bar, use the pre solid's outward normal.
    const bool intoDaughter = postPV->GetMotherLogical() == prePV->GetLogicalVolume();
    const auto& touch = intoDaughter ? post->GetTouchableHandle() : pre->GetTouchableHandle();
    const auto* surfacePV = intoDaughter ? postPV : prePV;
    const auto transform = touch->GetHistory()->GetTopTransform();
    const auto localPoint = transform.TransformPoint(post->GetPosition());
    const auto* solid = surfacePV->GetLogicalVolume()->GetSolid();
    auto localNormal = solid->SurfaceNormal(localPoint);
    if (intoDaughter) localNormal = -localNormal;
    const auto normal = transform.Inverse().TransformAxis(localNormal);
    const auto cosine = static_cast<G4float>(pre->GetMomentumDirection().dot(normal));
    auto* am = G4AnalysisManager::Instance();
    const int vals[] = {G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID(),
        step->GetTrack()->GetTrackID(), PhysicalObservation::Source(step->GetTrack()),
        prePV->GetCopyNo(), postPV->GetCopyNo(), VolumeId(prePV), VolumeId(postPV),
        PhysicalObservation::Outcome(status), exiting};
    for (int i = 0; i < 9; ++i) am->FillNtupleIColumn(2, i, vals[i]);
    am->FillNtupleFColumn(2, 9, cosine);
#ifdef EJ200_ENABLE_DIAGNOSTICS
    // These raw toolkit and normal-validation fields are intentionally absent
    // from production output. GetGlobalExitNormal is avoided because it writes
    // navigator caches in Geant4 11.4.
    const auto outward = intoDaughter ? -localNormal : localNormal;
    constexpr G4double probe = 1.e-5*mm;
    const bool orientationValid = solid->Inside(localPoint+probe*outward) != kInside &&
                                  solid->Inside(localPoint-probe*outward) != kOutside;
    const bool valid = solid->Inside(localPoint) == kSurface;
    am->FillNtupleIColumn(2, 10, status);
    am->FillNtupleIColumn(2, 11, valid);
    am->FillNtupleIColumn(2, 12, orientationValid);
    am->FillNtupleDColumn(2, 13, normal.mag());
    am->FillNtupleDColumn(2, 14, pre->GetKineticEnergy()/eV);
#endif
    am->AddNtupleRow(2);
}
