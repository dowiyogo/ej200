#pragma once
#include "G4OpBoundaryProcess.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

namespace PhysicalObservation {
// Engine-neutral result: 0 unknown, 1 reflection, 2 transmission,
// 3 absorption, 4 detection. Raw toolkit status is stored alongside it.
inline int Outcome(int status) {
    switch (status) {
    case FresnelRefraction: case Transmission: case SameMaterial:
    case CoatedDielectricRefraction: case CoatedDielectricFrustratedTransmission:
        return 2;
    case TotalInternalReflection: case FresnelReflection: case LambertianReflection:
    case LobeReflection: case SpikeReflection: case BackScattering:
    case CoatedDielectricReflection: return 1;
    case Absorption: case NoRINDEX: return 3;
    case Detection: return 4;
    default: return 0;
    }
}
inline int Source(const G4Track* track) {
    const auto* creator = track->GetCreatorProcess();
    if (!creator) return 0; // primary
    if (creator->GetProcessName() == "Scintillation") return 1;
    if (creator->GetProcessName() == "Cerenkov") return 2;
    return 3;
}
inline bool AcceptedIncident(int status) {
    return status == Detection || status == Absorption || Outcome(status) == 2;
}
} // namespace PhysicalObservation
