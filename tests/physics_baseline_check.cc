#include "DetectorConstruction.hh"
#include "Materials.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4PhysicsVector.hh"
#include "G4SystemOfUnits.hh"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
struct Expected {
    const char* name;
    G4Material* (*create)();
    G4double yield;
    G4double rise;
    G4double decay;
    G4double attenuation;
    G4double peakWavelength;
};

[[noreturn]] void Fail(const std::string& message) {
    std::cerr << "physics_baseline_check FAILED: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

void RequireEqual(const std::string& label, G4double actual, G4double expected) {
    if (actual != expected) {
        std::cerr << "physics_baseline_check FAILED: " << label
                  << " actual=" << actual << " expected=" << expected << '\n';
        std::exit(EXIT_FAILURE);
    }
}

G4double PeakWavelength(const G4MaterialPropertyVector* emission) {
    if (emission == nullptr || emission->GetVectorLength() == 0) {
        Fail("missing SCINTILLATIONCOMPONENT1");
    }
    std::size_t peakIndex = 0;
    for (std::size_t i = 1; i < emission->GetVectorLength(); ++i) {
        if ((*emission)[i] > (*emission)[peakIndex]) peakIndex = i;
    }
    return 1239.84193 * eV * nm / emission->Energy(peakIndex);
}

void CheckMaterial(const Expected& expected) {
    auto* material = expected.create();
    auto* mpt = material->GetMaterialPropertiesTable();
    if (mpt == nullptr) Fail(std::string(expected.name) + " has no material properties table");

    RequireEqual(std::string(expected.name) + " yield",
                 mpt->GetConstProperty("SCINTILLATIONYIELD"), expected.yield);
    RequireEqual(std::string(expected.name) + " rise",
                 mpt->GetConstProperty("SCINTILLATIONRISETIME1"), expected.rise);
    RequireEqual(std::string(expected.name) + " decay",
                 mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1"), expected.decay);
    RequireEqual(std::string(expected.name) + " resolution scale",
                 mpt->GetConstProperty("RESOLUTIONSCALE"), 1.0);
    RequireEqual(std::string(expected.name) + " prompt fraction",
                 mpt->GetConstProperty("SCINTILLATIONYIELD1"), 1.0);

    if (mpt->GetConstProperty("SCINTILLATIONRISETIME1") <= 0.0) {
        Fail(std::string(expected.name) + " rise time is not positive");
    }

    auto* attenuation = mpt->GetProperty("ABSLENGTH");
    if (attenuation == nullptr || attenuation->GetVectorLength() == 0) {
        Fail(std::string(expected.name) + " has no ABSLENGTH");
    }
    for (std::size_t i = 0; i < attenuation->GetVectorLength(); ++i) {
        RequireEqual(std::string(expected.name) + " attenuation",
                     (*attenuation)[i], expected.attenuation);
    }

    RequireEqual(std::string(expected.name) + " emission peak",
                 PeakWavelength(mpt->GetProperty("SCINTILLATIONCOMPONENT1")),
                 expected.peakWavelength);
}
} // namespace

int main() {
    Materials::EnableFiniteScintillationRiseTime();
    if (!G4OpticalParameters::Instance()->GetScintFiniteRiseTime()) {
        Fail("finite scintillation rise time is disabled");
    }

    const Expected materials[] = {
        {"EJ-204", Materials::CreateEJ204, 10400.0 / MeV, 0.7 * ns, 1.8 * ns,
         160.0 * cm, 408.0 * nm},
        {"EJ-200", Materials::CreateEJ200, 10000.0 / MeV, 0.9 * ns, 2.1 * ns,
         380.0 * cm, 425.0 * nm},
        {"EJ-230", Materials::CreateEJ230, 9700.0 / MeV, 0.5 * ns, 1.5 * ns,
         120.0 * cm, 391.0 * nm},
    };
    for (const auto& material : materials) CheckMaterial(material);

    DetectorConstruction detector;
    if (detector.GetScintillatorMaterial() != "EJ204") {
        Fail("default active scintillator is not EJ204");
    }

    std::cout << "physics_baseline_check PASSED\n";
    return EXIT_SUCCESS;
}
