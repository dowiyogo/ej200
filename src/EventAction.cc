#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* ra) : fRunAction(ra) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fNEndLeft  = 0;
    fNEndRight = 0;
    fNTop      = 0;

    // Extraer posicion x del vertice primario.
    // G4ParticleGun siempre crea exactamente un G4PrimaryVertex, por lo
    // que GetPrimaryVertex(0) es seguro.  En caso de que no exista
    // (no deberia ocurrir), se deja fGunXmm = 0.
    const G4PrimaryVertex* vtx = event ? event->GetPrimaryVertex(0) : nullptr;
    fGunXmm = vtx ? vtx->GetX0() / mm : 0.0;
}

void EventAction::EndOfEventAction(const G4Event*) {
    if (!fRunAction) return;

    // Flush per-event counters (incremented by SiPMSD) into run accumulables.
    fRunAction->AddEndLeft (fNEndLeft);
    fRunAction->AddEndRight(fNEndRight);
    fRunAction->AddTop     (fNTop);
    fRunAction->FlushEvent (fNEndLeft > 0 || fNEndRight > 0 || fNTop > 0);
}
