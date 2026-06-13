#pragma once
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "SiPMModel.hh"

#include <vector>

class G4GenericMessenger;

// --------------------------------------------------------------------------
// SiPMSD — sensitive detector for the two END SiPM arrays.
//
// PDE is applied by G4OpBoundaryProcess from the selected surface EFFICIENCY
// curve. This SD records only photons accepted by that boundary process.
//
// Electronic jitter:
//   The registered hit time is  t_hit = t_global + Gauss(0, σ_jitter).
//   Default σ_jitter = 20 ps.  Configurable via:
//     /sipm/jitterSigma <value> ns   (e.g. /sipm/jitterSigma 0.020 ns)
//
// Copy-number convention (set by DetectorConstruction):
//   0 –  7  →  end-left  (face -X)
//   8 – 15  →  end-right (face +X)
// --------------------------------------------------------------------------
class SiPMSD : public G4VSensitiveDetector
{
  public:
    explicit SiPMSD(const G4String& name);
    ~SiPMSD() override;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

    // Configurable jitter sigma (G4 internal time units)
    void     SetJitterSigma(G4double sigma) { fJitterSigma = sigma; }
    G4double GetJitterSigma() const         { return fJitterSigma; }
    void     SetModel(const G4String& model);
    G4String GetModel() const { return fModel; }

  private:
    G4double GetPDE(G4double energy) const;

    G4String fModel = "AFBR-S4N66P024M";
    std::vector<SiPMModel::PDEPoint> fPDECurve;

    // Jitter sigma in G4 internal units (nanoseconds scale).
    // Default: 20 ps = 0.020 ns.
    G4double fJitterSigma = 20.0e-3;  // 1 G4 time unit = 1 ns → 20e-3 = 20 ps

    G4GenericMessenger* fMessenger = nullptr;
};
