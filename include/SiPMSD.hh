#pragma once
#include "G4MaterialPropertyVector.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

// --------------------------------------------------------------------------
// SiPMSD — sensitive detector para todos los volúmenes SiPM
//           (end-left, end-right, top).
//
// La PDE se aplica aquí via trial de Bernoulli.  La interpolación usa
// G4MaterialPropertyVector::Value(), que es el mismo motor que G4 usa
// internamente para todas las propiedades ópticas (interpolación lineal
// en orden ascendente de energía).
//
// ¿Por qué no DETECTIONEFFICIENCY en el MPT de la superficie?
//   G4OpBoundaryProcess en G4 11.x mata los fotones etiquetados con
//   DETECTIONEFFICIENCY *antes* de que entren al volumen SiPM, por lo
//   que el SD nunca sería llamado.  Mantener la PDE aquí da control total.
//
// Convención de copy-number (fijada en DetectorConstruction):
//   0 –  7  →  end-left  (cara -X)
//   8 – 15  →  end-right (cara +X)
//  16 – 35  →  top       (cara +Y)
// --------------------------------------------------------------------------
class SiPMSD : public G4VSensitiveDetector
{
  public:
    explicit SiPMSD(const G4String& name);
    ~SiPMSD() override = default;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

  private:
    // Interpola la PDE a la energía del fotón [unidades internas de G4]
    G4double GetPDE(G4double energy) const;

    // PDE como función de energía del fotón, llenada en orden ascendente.
    // G4MaterialPropertyVector::Value() provee interpolación lineal con
    // clamping en los extremos de la tabla.
    G4MaterialPropertyVector fPDEVec;
};
