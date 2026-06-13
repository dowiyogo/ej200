#include "DetectorConstruction.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

#include <cstdlib>
#include <iostream>
#include <set>
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

    const auto& reflectors = detector.GetReflectorSurfaces();
    const std::set<G4String> nonEndFaces = {"-Y", "+Y", "-Z", "+Z"};
    G4SurfaceProperty* shared = nullptr;
    for (const auto& face : nonEndFaces) {
        const auto found = reflectors.find(face);
        Require(found != reflectors.end(), std::string(face) + " is not closed by Mylar");
        auto* optical = dynamic_cast<G4OpticalSurface*>(found->second->GetSurfaceProperty());
        Require(optical != nullptr && optical->GetName() == "MylarReflector",
                std::string(face) + " does not use MylarReflector");
        if (shared == nullptr) shared = optical;
        Require(shared == optical, "non-END faces do not share one uniform Mylar surface");
    }
    const auto plusY = reflectors.at("+Y")->GetVolume2()->GetLogicalVolume()->GetSolid();
    Require(plusY->GetEntityType() == "G4Box", "+Y Mylar panel contains TOP windows");

    std::cout << "endonly_geometry_check PASSED: 16 END, 0 TOP, +Y/-Y/+Z/-Z closed\n";
    return EXIT_SUCCESS;
}
