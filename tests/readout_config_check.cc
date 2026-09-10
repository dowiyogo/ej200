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

void RequirePanelReflectors(const DetectorConstruction& detector) {
    const auto& sipmSurfaces = detector.GetSiPMSurfaces();
    Require(!sipmSurfaces.empty(), "no SiPM surfaces");
    auto* barPV = sipmSurfaces.begin()->second->GetVolume1();
    Require(barPV != nullptr, "SiPM border has no volume1");
    Require(G4LogicalSkinSurface::GetSurface(barPV->GetLogicalVolume()) == nullptr,
            "explicit air-gap geometry must not have a bar skin reflector");

    const auto& panels = detector.GetReflectorSurfaces();
    Require(panels.size() == (detector.IsTopInstrumented() ? 3u : 4u),
            "wrong lateral reflector count");
    Require(panels.count("+Y") == (detector.IsTopInstrumented() ? 0u : 1u),
            "YPlus reflector must be present only without TOP sensors");
    for (const auto& entry : panels) {
        auto* border = entry.second;
        auto* air = border->GetVolume1();
        auto* wrap = border->GetVolume2();
        Require(air && wrap, "reflector border has missing volume");
        Require(air->GetName().find("AirGap") == 0 &&
                wrap->GetName().find("Vikuiti") == 0,
                "reflector border must point from air to wrap");
        auto* optical = dynamic_cast<G4OpticalSurface*>(border->GetSurfaceProperty());
        Require(optical && optical->GetType() == dielectric_metal &&
                optical->GetFinish() == polished, "wrong air-wrap reflector model");
        auto* mpt = optical->GetMaterialPropertiesTable();
        auto* reflectivity = mpt ? mpt->GetProperty("REFLECTIVITY") : nullptr;
        Require(reflectivity && reflectivity->GetVectorLength() > 0,
                "reflector has no reflectivity spectrum");
        for (std::size_t i = 0; i < reflectivity->GetVectorLength(); ++i)
            Require(std::abs((*reflectivity)[i] - 0.98) < 1e-12,
                    "reflector reflectivity differs from 0.98");
        auto* reverse = G4LogicalBorderSurface::GetSurface(wrap, air);
        Require(reverse && reverse->GetSurfaceProperty() == optical,
                "reverse wrap-air reflector border missing or inconsistent");
        auto* barAir = G4LogicalBorderSurface::GetSurface(barPV, air);
        auto* airBar = G4LogicalBorderSurface::GetSurface(air, barPV);
        Require(barAir && airBar, "bidirectional bar-air border missing");
        auto* barOptical = dynamic_cast<G4OpticalSurface*>(barAir->GetSurfaceProperty());
        Require(barOptical && barOptical->GetType() == dielectric_dielectric &&
                barOptical->GetFinish() == polished &&
                airBar->GetSurfaceProperty() == barOptical,
                "bar-air interface must retain polished dielectric model for TIR");
    }
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
    RequirePanelReflectors(detector);
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
    RequirePanelReflectors(detector);
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
    RequirePanelReflectors(detector);
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
