#pragma once
#include "G4MaterialPropertyVector.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4GenericMessenger;

// --------------------------------------------------------------------------
// SiPMSD — sensitive detector for all SiPM volumes (end-left, end-right, top).
//
// PDE is applied via a Bernoulli trial with linear interpolation over the
// Hamamatsu S13360-6025 table (300–940 nm, 33 points).
//
// Electronic jitter:
//   The registered hit time is  t_hit = t_global + Gauss(0, σ_jitter).
//   Default σ_jitter = 20 ps.  Configurable via:
//     /sipm/jitterSigma <value> ns   (e.g. /sipm/jitterSigma 0.020 ns)
//
// Copy-number convention (set by DetectorConstruction):
//   0 –  7  →  end-left  (face -X)
//   8 – 15  →  end-right (face +X)
//  16 …    →  top        (face +Y)
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

  private:
    G4double GetPDE(G4double energy) const;

    G4MaterialPropertyVector fPDEVec;

    // Jitter sigma in G4 internal units (nanoseconds scale).
    // Default: 20 ps = 0.020 ns.
    G4double fJitterSigma = 20.0e-3;  // 1 G4 time unit = 1 ns → 20e-3 = 20 ps

    G4GenericMessenger* fMessenger = nullptr;
};
