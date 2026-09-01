#pragma once

#include "globals.hh"

#include <vector>

namespace SiPMModel {

struct PDEPoint {
    G4double wavelengthNm;
    G4double efficiency;
};

G4String CanonicalName(const G4String& requested);
G4String DataFilePath(const G4String& canonicalModel);
std::vector<PDEPoint> LoadPDECurve(const G4String& canonicalModel);
G4double InterpolatePDE(const std::vector<PDEPoint>& curve, G4double wavelengthNm);

} // namespace SiPMModel
