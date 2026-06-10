#include "SiPMModel.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifndef EJ200_DATA_DIR
#define EJ200_DATA_DIR "data"
#endif

namespace SiPMModel {
namespace {
G4String Normalized(G4String value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char ch) {
                                   return ch == '-' || ch == '_' || ch == ' ';
                               }),
                value.end());
    return value;
}
} // namespace

G4String CanonicalName(const G4String& requested) {
    const G4String normalized = Normalized(requested);
    if (normalized == "AFBRS4N66P024M" || normalized == "BROADCOM") {
        return "AFBR-S4N66P024M";
    }
    return "";
}

G4String DataFilePath(const G4String& canonicalModel) {
    const G4String model = CanonicalName(canonicalModel);
    if (model.empty()) {
        throw std::runtime_error("unsupported SiPM model: " + std::string(canonicalModel));
    }

    const char* configuredRoot = std::getenv("EJ200_DATA_DIR");
    const G4String root = configuredRoot ? configuredRoot : EJ200_DATA_DIR;
    return root + "/sipm/" + model + "_pde.txt";
}

std::vector<PDEPoint> LoadPDECurve(const G4String& canonicalModel) {
    const G4String path = DataFilePath(canonicalModel);
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("unable to open SiPM PDE file: " + std::string(path));
    }

    std::vector<PDEPoint> curve;
    std::string line;
    while (std::getline(input, line)) {
        const auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;

        std::istringstream values(line);
        PDEPoint point{};
        if (!(values >> point.wavelengthNm >> point.efficiency)) {
            throw std::runtime_error("invalid SiPM PDE row in " + std::string(path));
        }
        if (point.efficiency < 0.0 || point.efficiency > 1.0) {
            throw std::runtime_error("SiPM PDE outside [0,1] in " + std::string(path));
        }
        if (!curve.empty() && point.wavelengthNm <= curve.back().wavelengthNm) {
            throw std::runtime_error("SiPM wavelengths are not strictly ascending in " +
                                     std::string(path));
        }
        curve.push_back(point);
    }

    if (curve.size() < 2) {
        throw std::runtime_error("SiPM PDE file has fewer than two points: " +
                                 std::string(path));
    }
    return curve;
}

G4double InterpolatePDE(const std::vector<PDEPoint>& curve, G4double wavelengthNm) {
    if (curve.empty() || wavelengthNm < curve.front().wavelengthNm ||
        wavelengthNm > curve.back().wavelengthNm) {
        return 0.0;
    }

    const auto upper = std::lower_bound(
        curve.begin(), curve.end(), wavelengthNm,
        [](const PDEPoint& point, G4double wavelength) {
            return point.wavelengthNm < wavelength;
        });
    if (upper == curve.begin()) return upper->efficiency;
    if (upper == curve.end()) return curve.back().efficiency;

    const auto lower = upper - 1;
    const G4double fraction =
        (wavelengthNm - lower->wavelengthNm) /
        (upper->wavelengthNm - lower->wavelengthNm);
    return lower->efficiency +
           fraction * (upper->efficiency - lower->efficiency);
}

} // namespace SiPMModel
