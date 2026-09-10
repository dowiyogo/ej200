#include "ActionInitialization.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#ifdef EJ200_ENABLE_DIAGNOSTICS
#include "TrackingAction.hh"
#endif

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
#ifdef EJ200_ENABLE_DIAGNOSTICS
    SetUserAction(new TrackingAction());
#endif
}
