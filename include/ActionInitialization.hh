#pragma once
#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization {
  public:
    ActionInitialization()  = default;
    ~ActionInitialization() = default;

    // Called for master thread (MT builds): only RunAction
    void BuildForMaster() const override;
    // Called for worker threads (or single-thread)
    void Build()          const override;
};
