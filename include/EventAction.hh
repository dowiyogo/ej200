#pragma once
#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

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

  private:
    RunAction* fRunAction = nullptr;
    G4int      fNEndLeft  = 0;
    G4int      fNEndRight = 0;
    G4int      fNTop      = 0;
    G4double   fGunXmm    = 0.0;
};
