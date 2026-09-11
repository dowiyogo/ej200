#pragma once
#include "G4UserTrackingAction.hh"
#include "globals.hh"
class G4Track;
namespace TerminalCensus {
void MarkKill(const G4Track*, const G4String& reason);
void Reset();
void Write(const G4String& prefix);
}
class TrackingAction : public G4UserTrackingAction {
public:
    void PreUserTrackingAction(const G4Track*) override;
    void PostUserTrackingAction(const G4Track*) override;
};
