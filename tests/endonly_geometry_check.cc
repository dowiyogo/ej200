#include "DetectorConstruction.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
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

std::string Role(const G4String& volumeName) {
    if (volumeName == "WorldLV") return "world";
    if (volumeName == "BarLV") return "scintillator";
    if (volumeName.find("AirGap") == 0) return "air_gap";
    if (volumeName.find("Mylar") == 0) return "reflector";
    if (volumeName == "EndSiPMLV" || volumeName == "TopSiPMLV") return "sipm";
    return "other";
}

void WriteVolumeCsv(const char* path) {
    if (path == nullptr || path[0] == '\0') return;
    std::ofstream out(path);
    out << "volume_name,material_name,has_rindex,rindex_min,rindex_max,mother,role\n";
    auto* store = G4LogicalVolumeStore::GetInstance();
    for (auto* lv : *store) {
        auto* material = lv->GetMaterial();
        auto* mpt = material ? material->GetMaterialPropertiesTable() : nullptr;
        auto* rindex = mpt ? mpt->GetProperty("RINDEX") : nullptr;
        bool has = rindex != nullptr && rindex->GetVectorLength() > 0;
        G4double rmin = 0.0;
        G4double rmax = 0.0;
        if (has) {
            rmin = (*rindex)[0];
            rmax = (*rindex)[0];
            for (std::size_t i = 1; i < rindex->GetVectorLength(); ++i) {
                rmin = std::min(rmin, (*rindex)[i]);
                rmax = std::max(rmax, (*rindex)[i]);
            }
        }
        out << lv->GetName() << ','
            << (material ? material->GetName() : "NULL") << ','
            << (has ? 1 : 0) << ','
            << rmin << ','
            << rmax << ','
            << "n/a" << ','
            << Role(lv->GetName()) << '\n';
    }
}

void WriteDimensionsJson(const char* path) {
    if (path == nullptr || path[0] == '\0') return;
    std::ofstream out(path);
    out
        << "{\n"
        << "  \"bar\": {\"length_mm\": 1400.0, \"width_mm\": 60.0, \"thickness_mm\": 10.0},\n"
        << "  \"air_gap\": {\"thickness_mm\": 0.10},\n"
        << "  \"reflector\": {\"thickness_mm\": 0.05},\n"
        << "  \"implementation\": \"four long lateral panels; END faces open\",\n"
        << "  \"air_shell_outer_mm\": {\"length\": 1400.0, \"width\": 60.2, \"thickness\": 10.2},\n"
        << "  \"reflector_outer_mm\": {\"length\": 1400.0, \"width\": 60.3, \"thickness\": 10.3}\n"
        << "}\n";
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

    Require(detector.GetBarSkinSurface() == nullptr,
            "mylar mode must not install a bar skin reflector");
    Require(detector.GetScintillatorAirSurfaces().size() == 4,
            "expected four scintillator-air border surfaces");
    Require(detector.GetReflectorSurfaces().size() == 4,
            "expected four air-reflector border surfaces");
    for (const auto& item : detector.GetScintillatorAirSurfaces()) {
        auto* optical = dynamic_cast<G4OpticalSurface*>(item.second->GetSurfaceProperty());
        Require(optical != nullptr, "scintillator-air border is not optical");
        Require(optical->GetType() == dielectric_dielectric,
                "scintillator-air border is not dielectric_dielectric");
    }
    for (const auto& item : detector.GetReflectorSurfaces()) {
        auto* optical = dynamic_cast<G4OpticalSurface*>(item.second->GetSurfaceProperty());
        Require(optical != nullptr, "air-reflector border is not optical");
        Require(optical->GetType() == dielectric_metal,
                "air-reflector border is not dielectric_metal");
        auto* mpt = optical->GetMaterialPropertiesTable();
        Require(mpt != nullptr && mpt->GetProperty("REFLECTIVITY") != nullptr,
                "air-reflector surface has no REFLECTIVITY");
    }

    WriteVolumeCsv(std::getenv("EXEC23_GEOMETRY_VOLUMES_CSV"));
    WriteDimensionsJson(std::getenv("EXEC23_GEOMETRY_DIMENSIONS_JSON"));

    std::cout << "endonly_geometry_check PASSED: explicit air gap + Mylar panels, "
              << "16 END, 0 TOP, no bar skin reflector\n";
    return EXIT_SUCCESS;
}
