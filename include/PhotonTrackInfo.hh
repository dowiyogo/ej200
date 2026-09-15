#pragma once

#include "G4VUserTrackInformation.hh"
#include "globals.hh"

// Minimal track-owned state for EXEC_46. Geant4 deletes it with the track.
class PhotonTrackInfo final : public G4VUserTrackInformation {
  public:
    PhotonTrackInfo() = default;
    ~PhotonTrackInfo() override = default;

    void AddBoundaryEncounter();
    G4int GetBoundaryEncounters() const { return fBoundaryEncounters; }

  private:
    G4int fBoundaryEncounters = 0;
};
