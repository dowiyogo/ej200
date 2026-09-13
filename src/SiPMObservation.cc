#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicalObservation.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

namespace {
std::uint64_t Key(G4int track, G4int sensor) {
    return (std::uint64_t(static_cast<std::uint32_t>(track)) << 32) |
           static_cast<std::uint32_t>(sensor);
}
}

void EventAction::BookSiPMObservations() {
    auto* am = G4AnalysisManager::Instance();
    am->CreateNtuple("sipm_event_counts", "Independent accepted incidents and SD detections per active SiPM/event");
    for (const auto* name : {"event_id", "global_id", "face_type", "incident", "detected",
        "detected_unique", "matched_detected", "unmatched_detected", "incident_scint",
        "reflected_encounters_excluded", "duplicate_incident", "surface_detection",
        "surface_absorption", "transmitted", "unknown_surface_pde"}) am->CreateNtupleIColumn(name);
    for (const auto* name : {"expected_surface_pde_sum", "surface_pde_variance_sum", "incident_wavelength_nm_sum"})
        am->CreateNtupleDColumn(name);
    am->FinishNtuple();
}

void EventAction::BeginSiPMObservations() {
    fIncidentKeys.clear();
    fDetectionKeys.clear();
    fSiPMObservations.clear();
    const auto* detector = dynamic_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detector) for (const auto& item : detector->GetSiPMSurfaces())
        fSiPMObservations.emplace(item.first, SiPMObservation{});
}

void EventAction::ObserveSiPMDetection(G4int trackId, G4int globalId) {
    auto& record = fSiPMObservations[globalId];
    ++record.detected;
    if (fDetectionKeys.insert(Key(trackId, globalId)).second) ++record.detectedUnique;
}

void EventAction::ObserveSiPMIncident(const G4Step* step, G4int status) {
    if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;
    const auto* pre = step->GetPreStepPoint()->GetPhysicalVolume();
    const auto* post = step->GetPostStepPoint()->GetPhysicalVolume();
    if (!pre || !post || pre == post) return;
    const auto& name = post->GetLogicalVolume()->GetName();
    if (name != "EndSiPMLV" && name != "TopSiPMLV") return;
    auto& record = fSiPMObservations[post->GetCopyNo()];
    if (PhysicalObservation::Outcome(status) == 1) ++record.reflected;
    if (!PhysicalObservation::AcceptedIncident(status)) return;
    if (!fIncidentKeys.insert(Key(step->GetTrack()->GetTrackID(), post->GetCopyNo())).second) {
        ++record.duplicateIncident;
        return;
    }
    ++record.incident;
    if (PhysicalObservation::Source(step->GetTrack()) == 1) ++record.incidentScint;
    if (status == Detection) ++record.surfaceDetection;
    else if (status == Absorption) ++record.absorption;
    else ++record.transmitted;
    const auto energy = step->GetPreStepPoint()->GetKineticEnergy();
    if (energy > 0.) record.wavelengthSum += 1239.84193/(energy/eV);
    // Actual directed boundary property, not an efficiency reconstructed from hits.
    const auto* border = G4LogicalBorderSurface::GetSurface(pre, post);
    const auto* surface = border ? dynamic_cast<const G4OpticalSurface*>(border->GetSurfaceProperty()) : nullptr;
    const auto* mpt = surface ? surface->GetMaterialPropertiesTable() : nullptr;
    const auto* efficiency = mpt ? mpt->GetProperty("EFFICIENCY") : nullptr;
    if (efficiency) {
        const auto p = efficiency->Value(energy);
        record.expectedPDE += p;
        record.variancePDE += p*(1.-p);
    } else ++record.unknownPDE;
}

void EventAction::EndSiPMObservations(G4int eventId) {
    // Boundary SD callbacks may precede stepping observations. Join only now.
    for (auto key : fDetectionKeys) if (fIncidentKeys.count(key))
        ++fSiPMObservations[static_cast<std::uint32_t>(key)].matched;
    auto* am = G4AnalysisManager::Instance();
    for (const auto& item : fSiPMObservations) {
        const auto& r = item.second;
        const int vals[] = {eventId, item.first, DetectorConstruction::FaceType(item.first),
            r.incident, r.detected, r.detectedUnique, r.matched, r.detectedUnique-r.matched,
            r.incidentScint, r.reflected, r.duplicateIncident, r.surfaceDetection,
            r.absorption, r.transmitted, r.unknownPDE};
        for (int i = 0; i < 15; ++i) am->FillNtupleIColumn(3, i, vals[i]);
        am->FillNtupleDColumn(3, 15, r.expectedPDE);
        am->FillNtupleDColumn(3, 16, r.variancePDE);
        am->FillNtupleDColumn(3, 17, r.wavelengthSum);
        am->AddNtupleRow(3);
    }
}
