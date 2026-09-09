#include "BoundaryCensus.hh"
#include "G4OpBoundaryProcess.hh"
#include <fstream>
#include <set>
#include <stdexcept>

namespace BoundaryCensus {
namespace {
// Nombres de todos los estados públicos de Geant4 11.4, sin agrupar reflexiones.
const char* StatusName(G4int status) {
    switch (status) {
        case Undefined: return "Undefined";
        case Transmission: return "Transmission";
        case FresnelRefraction: return "FresnelRefraction";
        case FresnelReflection: return "FresnelReflection";
        case TotalInternalReflection: return "TotalInternalReflection";
        case LambertianReflection: return "LambertianReflection";
        case LobeReflection: return "LobeReflection";
        case SpikeReflection: return "SpikeReflection";
        case BackScattering: return "BackScattering";
        case Absorption: return "Absorption";
        case Detection: return "Detection";
        case NotAtBoundary: return "NotAtBoundary";
        case SameMaterial: return "SameMaterial";
        case StepTooSmall: return "StepTooSmall";
        case NoRINDEX: return "NoRINDEX";
        case PolishedLumirrorAirReflection: return "PolishedLumirrorAirReflection";
        case PolishedLumirrorGlueReflection: return "PolishedLumirrorGlueReflection";
        case PolishedAirReflection: return "PolishedAirReflection";
        case PolishedTeflonAirReflection: return "PolishedTeflonAirReflection";
        case PolishedTiOAirReflection: return "PolishedTiOAirReflection";
        case PolishedTyvekAirReflection: return "PolishedTyvekAirReflection";
        case PolishedVM2000AirReflection: return "PolishedVM2000AirReflection";
        case PolishedVM2000GlueReflection: return "PolishedVM2000GlueReflection";
        case EtchedLumirrorAirReflection: return "EtchedLumirrorAirReflection";
        case EtchedLumirrorGlueReflection: return "EtchedLumirrorGlueReflection";
        case EtchedAirReflection: return "EtchedAirReflection";
        case EtchedTeflonAirReflection: return "EtchedTeflonAirReflection";
        case EtchedTiOAirReflection: return "EtchedTiOAirReflection";
        case EtchedTyvekAirReflection: return "EtchedTyvekAirReflection";
        case EtchedVM2000AirReflection: return "EtchedVM2000AirReflection";
        case EtchedVM2000GlueReflection: return "EtchedVM2000GlueReflection";
        case GroundLumirrorAirReflection: return "GroundLumirrorAirReflection";
        case GroundLumirrorGlueReflection: return "GroundLumirrorGlueReflection";
        case GroundAirReflection: return "GroundAirReflection";
        case GroundTeflonAirReflection: return "GroundTeflonAirReflection";
        case GroundTiOAirReflection: return "GroundTiOAirReflection";
        case GroundTyvekAirReflection: return "GroundTyvekAirReflection";
        case GroundVM2000AirReflection: return "GroundVM2000AirReflection";
        case GroundVM2000GlueReflection: return "GroundVM2000GlueReflection";
        case Dichroic: return "Dichroic";
        case CoatedDielectricReflection: return "CoatedDielectricReflection";
        case CoatedDielectricRefraction: return "CoatedDielectricRefraction";
        case CoatedDielectricFrustratedTransmission: return "CoatedDielectricFrustratedTransmission";
        default: return "UnknownStatus";
    }
}
// Escapado CSV: los nombres se preservan aunque contengan comillas o comas.
G4String CsvQuote(const G4String& value) {
    G4String result = "\"";
    for (char ch : value) {
        if (ch == '"') result += '"';
        result += ch;
    }
    return result + "\"";
}
}

Census& Instance() {
    static Census instance;
    return instance;
}
void Census::Record(const BoundaryKey& key) {
    const std::lock_guard<std::mutex> lock(mutex_);
    ++counts_[key];
}
void Census::Reset() {
    const std::lock_guard<std::mutex> lock(mutex_);
    counts_.clear();
}
void Census::Write(const G4String& filename) const {
    const std::lock_guard<std::mutex> lock(mutex_);
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open boundary census CSV: " + filename);
    out << "pre,pre_copy,post,post_copy,status_int,status_name,count\n";
    // Para cada par observado se exportan también los estados con cero entradas.
    std::set<std::tuple<G4String, G4int, G4String, G4int>> pairs;
    for (const auto& entry : counts_)
        pairs.emplace(entry.first.pre, entry.first.pre_copy,
                      entry.first.post, entry.first.post_copy);
    const auto writeRow = [&out](const BoundaryKey& key, G4long count) {
        out << CsvQuote(key.pre) << ',' << key.pre_copy << ','
            << CsvQuote(key.post) << ',' << key.post_copy << ','
            << key.status << ',' << StatusName(key.status) << ',' << count << '\n';
    };
    for (const auto& pair : pairs) {
        for (G4int status = Undefined; status <= CoatedDielectricFrustratedTransmission; ++status) {
            BoundaryKey key{std::get<0>(pair), std::get<1>(pair),
                            std::get<2>(pair), std::get<3>(pair), status};
            auto it = counts_.find(key);
            writeRow(key, it == counts_.end() ? 0 : it->second);
        }
    }
    // Un valor desconocido nunca se descarta silenciosamente.
    for (const auto& entry : counts_)
        if (entry.first.status < Undefined || entry.first.status > CoatedDielectricFrustratedTransmission)
            writeRow(entry.first, entry.second);
    out.flush();
    if (!out) throw std::runtime_error("Cannot write boundary census CSV: " + filename);
}
} // namespace BoundaryCensus
