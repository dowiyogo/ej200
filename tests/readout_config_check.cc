#include "DetectorConstruction.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include <cstdlib>
#include <iostream>
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
    const auto& surfaces = detector.GetSiPMSurfaces();
    Require(!surfaces.empty(), "no SiPM border surfaces");
    auto* barPV = surfaces.begin()->second->GetVolume1();
    Require(barPV != nullptr, "SiPM border has no bar physical volume");
    auto* barLV = barPV->GetLogicalVolume();
    Require(barLV != nullptr, "BarPV has no logical volume");

    auto* skin = G4LogicalSkinSurface::GetSurface(barLV);
    Require(skin != nullptr, "BarLV has no reflector skin");
    auto* optical = dynamic_cast<G4OpticalSurface*>(skin->GetSurfaceProperty());
    Require(optical != nullptr, "Bar skin is not an optical surface");
    Require(optical->GetName() == "BarSkinReflector",
        "Bar skin does not use BarSkinReflector");
    Require(detector.GetReflectorSurfaces().empty(),
        "reflector panel border surfaces should be removed");
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
    Require(detector.GetSiPMSurfaces().size() ==
                static_cast<std::size_t>(detector.GetNTopSiPMs()),
            "Top config has wrong SiPM surface count");
    RequireBarSkin(detector);
}

void CheckEndTop() {
    DetectorConstruction detector;
    detector.SetReadoutConfiguration("EndTop");
    detector.Construct();

    Require(detector.IsEndInstrumented() && detector.IsTopInstrumented(),
            "EndTop does not instrument both readouts");
    RequireBarSkin(detector);
}
} // namespace

int main() {
    CheckEnd();
    CheckTop();
    CheckEndTop();
    std::cout << "readout_config_check PASSED\n";
    return EXIT_SUCCESS;
}
