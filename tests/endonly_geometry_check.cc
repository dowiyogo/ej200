#include "DetectorConstruction.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
[[noreturn]] void Fail(const std::string& message) {
    std::cerr << "endonly_geometry_check FAILED: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

void Require(bool condition, const std::string& message) {
    if (!condition) Fail(message);
}
} // namespace

int main() {
    DetectorConstruction detector;
    detector.SetReadoutConfiguration("EndTop");
    detector.Construct();

    Require(detector.GetTopSurface() == "mylar", "default top surface is not mylar");
    Require(detector.GetNActiveEndSiPMs() == 16, "mylar mode does not keep exactly 16 END channels");
    Require(detector.GetNActiveTopSiPMs() == 0, "mylar mode registers active TOP channels");
    Require(detector.GetSiPMSurfaces().size() == 16, "mylar mode does not expose exactly 16 SiPM surfaces");
    for (G4int id = 0; id < 16; ++id) {
        Require(detector.GetSiPMSurfaces().count(id) == 1,
                "missing END channel " + std::to_string(id));
    }
    for (const auto& item : detector.GetSiPMSurfaces()) {
        Require(item.first < 16, "TOP channel registered in mylar mode");
    }

    auto* skin = detector.GetBarSkinSurface();
    Require(skin != nullptr, "bar logical skin surface is missing");
    auto* optical = dynamic_cast<G4OpticalSurface*>(skin->GetSurfaceProperty());
    Require(optical != nullptr, "bar skin surface is not an optical surface");
    Require(optical->GetName() == "MylarReflector",
            "bar skin surface does not use MylarReflector in default mylar mode");

    std::cout << "endonly_geometry_check PASSED: 16 END, 0 TOP, BarLV skin = MylarReflector\n";
    return EXIT_SUCCESS;
}
