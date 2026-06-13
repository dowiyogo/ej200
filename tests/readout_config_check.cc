#include "DetectorConstruction.hh"

#include "G4LogicalBorderSurface.hh"
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

void RequireReflector(const DetectorConstruction& detector, const G4String& face,
                      bool expected) {
    const auto& surfaces = detector.GetReflectorSurfaces();
    const auto found = surfaces.find(face);
    Require((found != surfaces.end()) == expected,
            std::string(face) + (expected ? " is not wrapped" : " is unexpectedly wrapped"));
    if (!expected) return;

    auto* optical = dynamic_cast<G4OpticalSurface*>(found->second->GetSurfaceProperty());
    Require(optical != nullptr, std::string(face) + " reflector is not an optical surface");
    Require(optical->GetName() == "BarSkinReflector",
            std::string(face) + " does not use BarSkinReflector");
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

void RequireSolidYPlusReflector(const DetectorConstruction& detector) {
    const auto& surfaces = detector.GetReflectorSurfaces();
    const auto found = surfaces.find("+Y");
    Require(found != surfaces.end(), "+Y reflector is missing");
    const G4String entityType =
        found->second->GetVolume2()->GetLogicalVolume()->GetSolid()->GetEntityType();
    Require(entityType == "G4Box", "+Y reflector is unexpectedly perforated");
}

void RequireEndPlacements(const DetectorConstruction& detector) {
    const auto& surfaces = detector.GetSiPMSurfaces();
    Require(surfaces.size() == DetectorConstruction::kNTotalSiPMs,
            "end-only geometry does not have 16 SiPM border surfaces");

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
    }
}

void CheckEndOnly() {
    DetectorConstruction detector;
    detector.Construct();

    Require(detector.GetNActiveEndSiPMs() == 16, "geometry does not activate 16 End SiPMs");
    Require(detector.GetSiPMSurfaces().size() == 16, "geometry has wrong SiPM surface count");
    RequireReflector(detector, "+Y", true);
    RequireSolidYPlusReflector(detector);
    RequireReflector(detector, "-X", false);
    RequireReflector(detector, "+X", false);
    RequireEndPlacements(detector);
    RequireSiPMSurfaceProperties(detector);
}
} // namespace

int main() {
    CheckEndOnly();
    std::cout << "readout_config_check PASSED\n";
    return EXIT_SUCCESS;
}
