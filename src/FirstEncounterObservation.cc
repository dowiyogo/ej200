#include "EventAction.hh"
#include "PhysicalObservation.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

void EventAction::BookFirstEncounters() {
    auto* am = G4AnalysisManager::Instance();
    am->CreateNtuple("first_bar_encounters", "One first physical bar encounter per optical track");
    for (const auto* name : {"event_id", "track_id", "source", "pre_copy", "post_copy",
                            "boundary_status", "outcome", "exiting_bar", "normal_valid",
                            "navigator_valid"}) am->CreateNtupleIColumn(name);
    am->CreateNtupleSColumn("pre_volume");
    am->CreateNtupleSColumn("post_volume");
    for (const auto* name : {"cos_incidence", "normal_norm", "normal_dot_navigator", "energy_eV"})
        am->CreateNtupleDColumn(name);
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
    // No navigator relocation and no normal sign correction: read-only cross-check.
    G4bool navigatorValid = false;
    const auto navigatorNormal = G4TransportationManager::GetTransportationManager()
        ->GetNavigatorForTracking()->GetGlobalExitNormal(post->GetPosition(), &navigatorValid);
    const bool valid = solid->Inside(localPoint) == kSurface;
    const auto cosine = pre->GetMomentumDirection().dot(normal);
    auto* am = G4AnalysisManager::Instance();
    const int vals[] = {G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID(),
        step->GetTrack()->GetTrackID(), PhysicalObservation::Source(step->GetTrack()),
        prePV->GetCopyNo(), postPV->GetCopyNo(), status, PhysicalObservation::Outcome(status),
        exiting, valid, navigatorValid};
    for (int i = 0; i < 10; ++i) am->FillNtupleIColumn(2, i, vals[i]);
    am->FillNtupleSColumn(2, 10, prePV->GetName());
    am->FillNtupleSColumn(2, 11, postPV->GetName());
    am->FillNtupleDColumn(2, 12, cosine);
    am->FillNtupleDColumn(2, 13, normal.mag());
    am->FillNtupleDColumn(2, 14, navigatorValid ? normal.dot(navigatorNormal) : -2.);
    am->FillNtupleDColumn(2, 15, pre->GetKineticEnergy()/eV);
    am->AddNtupleRow(2);
}
