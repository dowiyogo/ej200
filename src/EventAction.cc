#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4OpticalPhoton.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

EventAction::EventAction(RunAction* ra) : fRunAction(ra) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fNEndLeft  = 0;
    fNEndRight = 0;
    fNTop      = 0;
    fEdep = fNonIonizing = fOpticalEdep = 0.;
    fProducedScint = fProducedOptical = 0;
    fFirstEncounterTracks.clear();

    // Extraer posicion x del vertice primario.
    // G4ParticleGun siempre crea exactamente un G4PrimaryVertex, por lo
    // que GetPrimaryVertex(0) es seguro.  En caso de que no exista
    // (no deberia ocurrir), se deja fGunXmm = 0.
    const G4PrimaryVertex* vtx = event ? event->GetPrimaryVertex(0) : nullptr;
    fGunXmm = vtx ? vtx->GetX0() / mm : 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    auto* am = G4AnalysisManager::Instance();
    const G4int ints[] = {event->GetEventID(), fProducedScint, fProducedOptical,
                         fNEndLeft, fNEndRight, fNTop};
    for (G4int i = 0; i < 6; ++i) am->FillNtupleIColumn(1, i, ints[i]);
    const G4double values[] = {fEdep/MeV, fNonIonizing/MeV,
                              (fEdep-fNonIonizing)/MeV, fOpticalEdep/MeV, fGunXmm};
    for (G4int i = 0; i < 5; ++i) am->FillNtupleDColumn(1, 6+i, values[i]);
    am->AddNtupleRow(1);
    if (!fRunAction) return;

    // Flush per-event counters (incremented by SiPMSD) into run accumulables.
    fRunAction->AddEndLeft (fNEndLeft);
    fRunAction->AddEndRight(fNEndRight);
    fRunAction->AddTop     (fNTop);
    fRunAction->FlushEvent (fNEndLeft > 0 || fNEndRight > 0 || fNTop > 0);
}

void EventAction::BookEnergyObservations() {
    auto* am = G4AnalysisManager::Instance();
    am->CreateNtuple("event_observables", "Energy and generated photons in BarLV per generated event");
    for (const auto* name : {"event_id", "produced_scint", "produced_optical",
                            "detected_left", "detected_right", "detected_top"})
        am->CreateNtupleIColumn(name);
    for (const auto* name : {"edep_total_MeV", "edep_nonionizing_MeV",
                            "edep_ionizing_MeV", "edep_optical_MeV", "gun_x_mm"})
        am->CreateNtupleDColumn(name);
    am->FinishNtuple();
}

void EventAction::ObserveEnergy(const G4Step* step) {
    const auto* pv = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv || pv->GetLogicalVolume()->GetName() != "BarLV") return;
    fEdep += step->GetTotalEnergyDeposit();
    fNonIonizing += step->GetNonIonizingEnergyDeposit();
    if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::Definition())
        fOpticalEdep += step->GetTotalEnergyDeposit();
    // Observe secondaries at creation, including ones not subsequently tracked.
    const auto* secondaries = step->GetSecondaryInCurrentStep();
    if (!secondaries) return;
    for (const auto* child : *secondaries) {
        if (child->GetDefinition() != G4OpticalPhoton::Definition()) continue;
        ++fProducedOptical;
        const auto* creator = child->GetCreatorProcess();
        if (creator && creator->GetProcessName() == "Scintillation") ++fProducedScint;
    }
}
