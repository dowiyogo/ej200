#pragma once
#include "G4UserEventAction.hh"
#include "globals.hh"
#include <unordered_set>
#include <map>
#include <cstdint>

class RunAction;
class G4Step;

// Contadores por evento — llenados por SiPMSD, volcados a RunAction al
// final del evento.
//
// Tambien almacena la posicion x del vertice primario (gun_x_mm) para que
// SiPMSD pueda escribirla en cada fila del ntuple.  Esto es esencial para
// el analisis de resolucion temporal vs posicion longitudinal (scan.mac).
class EventAction : public G4UserEventAction {
  public:
    explicit EventAction(RunAction* ra);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction  (const G4Event*) override;

    void AddEndLeftHit()  { ++fNEndLeft;  }
    void AddEndRightHit() { ++fNEndRight; }
    void AddTopHit()      { ++fNTop;      }

    // Posicion x del muon primario en este evento [mm].
    // Extraida del G4PrimaryVertex en BeginOfEventAction y consumida
    // por SiPMSD::ProcessHits para llenar la columna gun_x_mm del ntuple.
    G4double GetGunXmm() const { return fGunXmm; }
    // Passive production observations; no RNG or tracking mutations.
    static void BookEnergyObservations();
    void ObserveEnergy(const G4Step*);
    static void BookFirstEncounters();
    void ObserveFirstEncounter(const G4Step*, G4int boundaryStatus);
    static void BookSiPMObservations();
    void BeginSiPMObservations();
    void EndSiPMObservations(G4int eventId);
    void ObserveSiPMIncident(const G4Step*, G4int boundaryStatus);
    void ObserveSiPMDetection(G4int trackId, G4int globalId);
    void RegisterDetectedTrackId(G4int trackId);

  private:
    RunAction* fRunAction = nullptr;
    G4int      fNEndLeft  = 0;
    G4int      fNEndRight = 0;
    G4int      fNTop      = 0;
    G4double   fGunXmm    = 0.0;
    G4double fEdep = 0., fNonIonizing = 0., fOpticalEdep = 0.;
    G4int fProducedScint = 0, fProducedOptical = 0;
    std::unordered_set<G4int> fFirstEncounterTracks;
    struct SiPMObservation {
        G4int incident = 0, detected = 0, detectedUnique = 0, matched = 0;
        G4int incidentScint = 0, reflected = 0, duplicateIncident = 0;
        G4int surfaceDetection = 0, absorption = 0, transmitted = 0, unknownPDE = 0;
        G4double expectedPDE = 0., variancePDE = 0., wavelengthSum = 0.;
    };
    std::map<G4int, SiPMObservation> fSiPMObservations;
    std::unordered_set<std::uint64_t> fIncidentKeys, fDetectionKeys;
    std::unordered_set<G4int> fDetectedTrackIds;
};
