#include "DetectorConstruction.hh"
#include "Materials.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4PhysicsVector.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
[[noreturn]] void Fail(const std::string& message) {
    std::cerr << "sslg4_properties_check FAILED: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

void RequireNear(const std::string& label, G4double actual, G4double expected,
                 G4double tolerance) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "sslg4_properties_check FAILED: " << label
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
} // namespace

int main() {
    Materials::EnableFiniteScintillationRiseTime();
    if (!G4OpticalParameters::Instance()->GetScintFiniteRiseTime()) {
        Fail("finite scintillation rise time is disabled");
    }

    DetectorConstruction detector;
    detector.SetScintillatorCode("OPSC-101");
    detector.Construct();
    if (detector.GetScintillatorCode() != "OPSC-101") {
        Fail("default active scintillator is not OPSC-101");
    }

    auto* material = detector.GetActiveScintillatorMaterial();
    if (material == nullptr || material->GetName() != "opsc-101") {
        Fail("bar material was not created by OrganicScintillatorFactory");
    }
    auto* mpt = material->GetMaterialPropertiesTable();
    if (mpt == nullptr) Fail("OPSC-101 has no material properties table");

    RequireNear("yield", mpt->GetConstProperty("SCINTILLATIONYIELD"),
                10400.0 / MeV, 1.0e-9 / MeV);
    RequireNear("rise", mpt->GetConstProperty("SCINTILLATIONRISETIME1"),
                0.7 * ns, 1.0e-12 * ns);
    RequireNear("decay", mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1"),
                1.8 * ns, 1.0e-12 * ns);

    auto* attenuation = mpt->GetProperty("ABSLENGTH");
    if (attenuation == nullptr || attenuation->GetVectorLength() == 0) {
        Fail("missing ABSLENGTH");
    }
    for (std::size_t i = 0; i < attenuation->GetVectorLength(); ++i) {
        RequireNear("ABSLENGTH", (*attenuation)[i], 160.0 * cm, 1.0e-9 * cm);
    }

    auto* rindex = mpt->GetProperty("RINDEX");
    if (rindex == nullptr || rindex->GetVectorLength() == 0) Fail("missing RINDEX");
    for (std::size_t i = 0; i < rindex->GetVectorLength(); ++i) {
        RequireNear("RINDEX", (*rindex)[i], 1.58, 1.0e-12);
    }

    RequireNear("emission peak",
                PeakWavelength(mpt->GetProperty("SCINTILLATIONCOMPONENT1")),
                408.8 * nm, 2.0 * nm);

    std::cout << "\n=== Effective SSLG4 OPSC-101 MPT dump ===\n";
    mpt->DumpTable();
    std::cout << "sslg4_properties_check PASSED\n";
    return EXIT_SUCCESS;
}
