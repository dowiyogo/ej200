#include "DetectorConstruction.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>

namespace {
[[noreturn]] void Fail(const std::string& message) {
    std::cerr << "readout_config_check FAILED: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string& message) {
    if (!condition) Fail(message);
}

void RequireBarSkin(const DetectorConstruction& detector) {
    const auto& sipmSurfaces = detector.GetSiPMSurfaces();
    Require(!sipmSurfaces.empty(), "no SiPM surfaces");
    auto* barPV = sipmSurfaces.begin()->second->GetVolume1();
    Require(barPV != nullptr, "SiPM border has no volume1");
    auto* barLV = barPV->GetLogicalVolume();
    Require(barLV != nullptr, "BarPV has no logical volume");

    auto* skin = G4LogicalSkinSurface::GetSurface(barLV);
    Require(skin != nullptr, "BarLV has no reflector skin surface");
    auto* optical = dynamic_cast<G4OpticalSurface*>(skin->GetSurfaceProperty());
    Require(optical != nullptr, "Bar skin is not an optical surface");
    Require(optical->GetName() == "BarSkinReflector",
        "Bar skin does not use BarSkinReflector");
    Require(detector.GetReflectorSurfaces().empty(),
        "reflector panel border surfaces should be removed");
}

void RequireSiPMSurfaceProperties(const DetectorConstruction& detector) {
    const auto& surfaces = detector.GetSiPMSurfaces();
    Require(!surfaces.empty(), "no SiPM surfaces");
    auto* optical = dynamic_cast<G4OpticalSurface*>(
        surfaces.begin()->second->GetSurfaceProperty());
    Require(optical != nullptr, "SiPM border is not an optical surface");
    auto* mpt = optical->GetMaterialPropertiesTable();
    Require(mpt != nullptr, "SiPM surface has no MPT");
    auto* efficiency = mpt->GetProperty("EFFICIENCY");
    auto* reflectivity = mpt->GetProperty("REFLECTIVITY");
    Require(efficiency != nullptr, "SiPM surface has no EFFICIENCY");
    Require(reflectivity != nullptr, "SiPM surface has no REFLECTIVITY");
    const G4double energy420 = 1239.84193 * eV * nm / (420.0 * nm);
    Require(std::abs(efficiency->Value(energy420) - 0.63) < 1.0e-6,
            "Broadcom PDE is not 63% at 420 nm");
    Require(std::abs(reflectivity->Value(energy420)) < 1.0e-12,
            "SiPM reflectivity is not zero");
}

void RequireEndTopPlacements(const DetectorConstruction& detector) {
    const auto& surfaces = detector.GetSiPMSurfaces();
    Require(surfaces.size() == DetectorConstruction::kNTotalSiPMs,
            "EndTop does not have 86 SiPM border surfaces");

    std::set<G4int> copyNumbers;
    for (G4int globalId = 0; globalId < DetectorConstruction::kNTotalSiPMs; ++globalId) {
        const auto found = surfaces.find(globalId);
        Require(found != surfaces.end(), "missing global ID " + std::to_string(globalId));
        auto* physical = found->second->GetVolume2();
        Require(physical != nullptr, "surface has no SiPM physical volume");
        Require(physical->GetCopyNo() == globalId,
                "copy number differs from global ID " + std::to_string(globalId));
        Require(copyNumbers.insert(physical->GetCopyNo()).second,
                "duplicate copy number " + std::to_string(globalId));
        if (globalId >= 16) {
            const G4int localId = globalId - 16;
            const G4double expected = DetectorConstruction::TopSiPMCenterX(localId);
            Require(std::abs(physical->GetTranslation().x() - expected) < 0.01 * mm,
                    "wrong Top x for global ID " + std::to_string(globalId));
            Require(std::abs(physical->GetTranslation().z()) < 0.01 * mm,
                    "wrong Top z for global ID " + std::to_string(globalId));
        }
    }
    Require(std::abs(DetectorConstruction::TopSiPMCenterX(35) -
                     DetectorConstruction::TopSiPMCenterX(34) - 24.0 * mm) < 0.01 * mm,
            "central Top pair is not separated by 24 mm");
}

void CheckEnd() {
    DetectorConstruction detector;
    Require(detector.GetReadoutConfiguration() == "End", "default configuration is not End");
    detector.Construct();

    Require(detector.IsEndInstrumented(), "End config does not instrument the ends");
    Require(!detector.IsTopInstrumented(), "End config instruments Top");
    Require(detector.GetNActiveEndSiPMs() == 16, "End config does not activate 16 End SiPMs");
    Require(detector.GetNActiveTopSiPMs() == 0, "End config activates Top SiPMs");
    Require(detector.GetSiPMSurfaces().size() == 16, "End config has wrong SiPM surface count");
    RequireBarSkin(detector);
    RequireSiPMSurfaceProperties(detector);
}

void CheckTop() {
    DetectorConstruction detector;
    detector.SetReadoutConfiguration("Top");
    detector.Construct();

    Require(!detector.IsEndInstrumented(), "Top config instruments End");
    Require(detector.IsTopInstrumented(), "Top config does not instrument Top");
    Require(detector.GetNActiveEndSiPMs() == 0, "Top config activates End SiPMs");
    Require(detector.GetNActiveTopSiPMs() == detector.GetNTopSiPMs(),
            "Top config has wrong active Top SiPM count");
    Require(detector.GetNActiveTopSiPMs() == 70, "Top config does not activate 70 Top SiPMs");
    Require(detector.GetSiPMSurfaces().size() ==
                static_cast<std::size_t>(detector.GetNTopSiPMs()),
            "Top config has wrong SiPM surface count");
    RequireBarSkin(detector);
    RequireSiPMSurfaceProperties(detector);
}

void CheckEndTop() {
    DetectorConstruction detector;
    detector.SetReadoutConfiguration("EndTop");
    detector.Construct();

    Require(detector.IsEndInstrumented() && detector.IsTopInstrumented(),
            "EndTop does not instrument both readouts");
    Require(detector.GetNActiveEndSiPMs() == 16, "EndTop does not activate 16 End SiPMs");
    Require(detector.GetNActiveTopSiPMs() == 70, "EndTop does not activate 70 Top SiPMs");
    RequireBarSkin(detector);
    RequireEndTopPlacements(detector);
    RequireSiPMSurfaceProperties(detector);
}
} // namespace

int main() {
    CheckEnd();
    CheckTop();
    CheckEndTop();
    std::cout << "readout_config_check PASSED\n";
    return EXIT_SUCCESS;
}
