#include "ActionInitialization.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
    auto* run   = new RunAction();
    auto* event = new EventAction(run);

    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(run);
    SetUserAction(event);
    SetUserAction(new SteppingAction());
}
