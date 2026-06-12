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

    // EXEC_13 / EJ-230: verify OPSC-106 properties against datasheet
    // Eljen EJ-228/EJ-230, Rev. Aug 2023
    DetectorConstruction detector;
    detector.SetScintillatorCode("OPSC-106");
    detector.Construct();
    if (detector.GetScintillatorCode() != "OPSC-106") {
        Fail("active scintillator code is not OPSC-106");
    }

    auto* material = detector.GetActiveScintillatorMaterial();
    if (material == nullptr || material->GetName() != "opsc-106") {
        Fail("bar material was not created by OrganicScintillatorFactory for opsc-106");
    }
    auto* mpt = material->GetMaterialPropertiesTable();
    if (mpt == nullptr) Fail("OPSC-106 has no material properties table");

    // yield = 9700 ph/MeV (datasheet)
    RequireNear("yield", mpt->GetConstProperty("SCINTILLATIONYIELD"),
                9700.0 / MeV, 1.0e-9 / MeV);
    // rise time = 0.5 ns (datasheet)
    RequireNear("rise", mpt->GetConstProperty("SCINTILLATIONRISETIME1"),
                0.5 * ns, 1.0e-12 * ns);
    // decay time = 1.5 ns (datasheet)
    RequireNear("decay", mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1"),
                1.5 * ns, 1.0e-12 * ns);

    // attenuation length = 120 cm (datasheet)
    auto* attenuation = mpt->GetProperty("ABSLENGTH");
    if (attenuation == nullptr || attenuation->GetVectorLength() == 0) {
        Fail("missing ABSLENGTH");
    }
    for (std::size_t i = 0; i < attenuation->GetVectorLength(); ++i) {
        RequireNear("ABSLENGTH", (*attenuation)[i], 120.0 * cm, 1.0e-9 * cm);
    }

    // refractive index = 1.58 (datasheet)
    auto* rindex = mpt->GetProperty("RINDEX");
    if (rindex == nullptr || rindex->GetVectorLength() == 0) Fail("missing RINDEX");
    for (std::size_t i = 0; i < rindex->GetVectorLength(); ++i) {
        RequireNear("RINDEX", (*rindex)[i], 1.58, 1.0e-12);
    }

    // emission peak ≈ 391 nm (datasheet); SSLG4 data peaks at ~390.5 nm
    RequireNear("emission peak",
                PeakWavelength(mpt->GetProperty("SCINTILLATIONCOMPONENT1")),
                391.0 * nm, 2.0 * nm);

    std::cout << "\n=== Effective SSLG4 OPSC-106 (EJ-230) MPT dump ===\n";
    mpt->DumpTable();
    std::cout << "sslg4_properties_check PASSED (OPSC-106/EJ-230)\n";
    return EXIT_SUCCESS;
}
