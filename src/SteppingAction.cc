#include "SteppingAction.hh"

#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* track = step->GetTrack();

    // Kill optical photons that have left the world volume (safety net).
    // In normal operation the wrapping surface reflects them back; this
    // handles edge cases (e.g. photon emitted very close to a boundary).
    if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
        if (step->GetPostStepPoint()->GetPhysicalVolume() == nullptr) {
            track->SetTrackStatus(fStopAndKill);
        }
    }
}
