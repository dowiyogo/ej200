#pragma once

#include "G4OpBoundaryProcess.hh"
#include "G4ThreeVector.hh"
#include "G4UserTrackingAction.hh"

class G4GenericMessenger;
class G4Track;

class DisplayTrackingAction : public G4UserTrackingAction {
  public:
    DisplayTrackingAction();
    ~DisplayTrackingAction() override;

    void PreUserTrackingAction(const G4Track*) override;

  private:
    G4GenericMessenger* fMessenger = nullptr;
    int                 fMaxOpticalTrajectories = 100;
    bool                fStoreMuon = true;
    bool                fStoreOtherParticles = false;
};

namespace DisplayDiagnostics {
void BeginEvent();
void EndEvent(int eventId);

void RecordBoundaryStatus(const G4Track* track,
                          G4OpBoundaryProcessStatus status,
                          const G4ThreeVector& position);

int GetConfiguredMaxOpticalTrajectories();
int GetOpticalTracksSeen();
int GetOpticalTrajectoriesStored();
}
