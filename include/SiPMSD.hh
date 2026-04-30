#pragma once
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include <vector>

class G4GenericMessenger;

// --------------------------------------------------------------------------
// SiPMSD — sensitive detector for all SiPM volumes (end-left, end-right, top).
//
// PDE is applied via a Bernoulli trial with linear interpolation over the
// selected sensor model curve.
//
// Electronic jitter:
//   The registered hit time is  t_hit = t_global + Gauss(0, σ_jitter).
//   Configurable via:
//     /sipm/sptrSigma <value> ns   (e.g. /sipm/sptrSigma 0.030 ns)
//
// Copy-number convention (set by DetectorConstruction):
//   0 –  7  →  end-left  (face -X)
//   8 – 15  →  end-right (face +X)
//  16 …    →  top        (face +Y)
// --------------------------------------------------------------------------
class SiPMSD : public G4VSensitiveDetector
{
  public:
    enum class SensorModel : G4int {
        kBroadcomAFBRS4N66P024M = 0,
        kHamamatsuS13360        = 1,
        kHamamatsuS14160_10um   = 2,
        kHamamatsuS14160_15um   = 3
    };

    explicit SiPMSD(const G4String& name);
    ~SiPMSD() override;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

    void     SetSPTRSigma(G4double sigma)        { fSPTRSigma = sigma; }
    G4double GetSPTRSigma() const                { return fSPTRSigma; }

    void     SetElectronicsSigma(G4double sigma) { fElectronicsSigma = sigma; }
    G4double GetElectronicsSigma() const         { return fElectronicsSigma; }

    void        SetSensorModelByName(const G4String& modelName);
    SensorModel GetSensorModel() const { return fSensorModel; }
    G4int       GetSensorModelId() const;
    G4String    GetSensorModelName() const;

  private:
    G4double GetPDE(G4double energy) const;
    void     ApplySensorDefaults(SensorModel model);
    void     LoadPDECurve(SensorModel model);

    std::vector<G4double> fPDEWavelengthNm;
    std::vector<G4double> fPDEValue;

    SensorModel fSensorModel      = SensorModel::kBroadcomAFBRS4N66P024M;
    G4double    fSPTRSigma        = 30.0e-3;
    G4double    fElectronicsSigma = 30.0e-3;

    G4GenericMessenger* fMessenger = nullptr;
};
