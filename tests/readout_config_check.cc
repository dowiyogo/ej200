#include "DetectorConstruction.hh"

#include "G4LogicalBorderSurface.hh"
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

void CheckEnd() {
    DetectorConstruction detector;
    Require(detector.GetReadoutConfiguration() == "End", "default configuration is not End");
    detector.Construct();

    Require(detector.IsEndInstrumented(), "End config does not instrument the ends");
    Require(!detector.IsTopInstrumented(), "End config instruments Top");
    Require(detector.GetNActiveEndSiPMs() == 16, "End config does not activate 16 End SiPMs");
    Require(detector.GetNActiveTopSiPMs() == 0, "End config activates Top SiPMs");
    Require(detector.GetSiPMSurfaces().size() == 16, "End config has wrong SiPM surface count");
    RequireReflector(detector, "+Y", true);
    RequireReflector(detector, "-X", false);
    RequireReflector(detector, "+X", false);
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
    RequireReflector(detector, "+Y", false);
    RequireReflector(detector, "-X", true);
    RequireReflector(detector, "+X", true);
}

void CheckEndTop() {
    DetectorConstruction detector;
    detector.SetReadoutConfiguration("EndTop");
    detector.Construct();

    Require(detector.IsEndInstrumented() && detector.IsTopInstrumented(),
            "EndTop does not instrument both readouts");
    RequireReflector(detector, "+Y", false);
    RequireReflector(detector, "-X", false);
    RequireReflector(detector, "+X", false);
}
} // namespace

int main() {
    CheckEnd();
    CheckTop();
    CheckEndTop();
    std::cout << "readout_config_check PASSED\n";
    return EXIT_SUCCESS;
}
