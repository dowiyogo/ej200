#pragma once
#include "G4Accumulable.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction {
  public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction  (const G4Run*) override;

    // Called by EventAction at end-of-event with per-event photon counts
    void AddEndLeft (G4int n) { fNEndLeft  += n; }
    void AddEndRight(G4int n) { fNEndRight += n; }
    void AddTop     (G4int n) { fNTop      += n; }
    void FlushEvent (G4bool anyHit) { if (anyHit) fNEventsWithHits += 1; }

  private:
    G4Accumulable<G4int> fNEndLeft       {"NEndLeft",        0};
    G4Accumulable<G4int> fNEndRight      {"NEndRight",       0};
    G4Accumulable<G4int> fNTop           {"NTop",            0};
    G4Accumulable<G4int> fNEventsWithHits{"NEventsWithHits", 0};
};
