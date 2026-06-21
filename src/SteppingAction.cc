#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

static constexpr G4double kHC_eVnm = 1239.84193;

namespace {
std::atomic<long long> gBarToMylar{0};
std::atomic<long long> gMylarToWorld{0};
std::atomic<long long> gMylarReflected{0};
std::atomic<long long> gMylarToSiPM{0};
std::atomic<long long> gKilledWorld{0};
std::atomic<long long> gExec26ScintAirSurfaceLossKilled{0};

struct SpeedStats {
    long long nSteps = 0;
    long long nSuper = 0;
    long long nMissingRindex = 0;
    double maxRatio = 0.0;
};

struct TrackHistory {
    int nScintAirReflections = 0;
    int nAirReflectorReflections = 0;
    int nTotalInternalReflections = 0;
    bool visitedAirGap = false;
    bool visitedReflectorBoundary = false;
};

struct DetectionRow {
    int eventId = -1;
    int trackId = -1;
    int globalId = -1;
    int detectedEndSide = -1;
    int nScintAirReflections = 0;
    int nAirReflectorReflections = 0;
    int nTotalInternalReflections = 0;
    bool visitedAirGap = false;
    bool visitedReflectorBoundary = false;
    double firstDetectionTimeNs = 0.0;
    double gunXmm = 0.0;
};

std::mutex gAuditMutex;
std::map<std::string, long long> gVolumeSteps;
std::map<std::string, long long> gBoundaryStatusCounts;
std::map<std::string, SpeedStats> gSpeedStats;
std::unordered_map<long long, TrackHistory> gTrackHistory;
std::vector<DetectionRow> gDetections;

bool IsBarLV(const G4String& name) { return name == "BarLV"; }
bool IsAirGapLV(const G4String& name) { return name.find("AirGap") == 0; }
bool IsReflectorLV(const G4String& name) { return name.find("Mylar") == 0; }
bool IsSiPMLV(const G4String& name) {
    return name == "EndSiPMLV" || name == "TopSiPMLV";
}
bool IsScintillatorMaterial(const G4String& name) {
    const std::string s = name;
    return s == "OPSC-101" || s == "opsc-101" || s == "EJ-204" ||
           s == "ej-204" || s == "EJ-200" || s == "ej-200";
}
bool IsReflectionStatus(G4OpBoundaryProcessStatus status) {
    return status == FresnelReflection ||
           status == TotalInternalReflection ||
           status == LambertianReflection ||
           status == LobeReflection ||
           status == SpikeReflection ||
           status == BackScattering;
}
long long TrackKey(int eventId, int trackId) {
    return (static_cast<long long>(eventId) << 32) ^
           static_cast<unsigned int>(trackId);
}
std::string StatusName(G4OpBoundaryProcessStatus status) {
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
        case CoatedDielectricFrustratedTransmission:
            return "CoatedDielectricFrustratedTransmission";
    }
    return "Unknown";
}
G4OpBoundaryProcess* BoundaryProcess() {
    static G4OpBoundaryProcess* process = nullptr;
    if (process != nullptr) return process;
    auto* manager = G4OpticalPhoton::Definition()->GetProcessManager();
    if (manager == nullptr) return nullptr;
    auto* processes = manager->GetProcessList();
    for (G4int i = 0; i < processes->size(); ++i) {
        process = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]);
        if (process != nullptr) return process;
    }
    return nullptr;
}
int CurrentEventId() {
    auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    return event ? event->GetEventID() : -1;
}
double CurrentGunXmm() {
    auto* ea = dynamic_cast<EventAction*>(
        G4EventManager::GetEventManager()->GetUserEventAction());
    return ea ? ea->GetGunXmm() : 0.0;
}
std::string OutputDir() {
    const char* dir = std::getenv("SHIP_OPTICAL_AUDIT_DIR");
    return (dir != nullptr && dir[0] != '\0') ? std::string(dir) : std::string();
}
std::string RunTag() {
    const char* tag = std::getenv("SHIP_OPTICAL_AUDIT_TAG");
    return (tag != nullptr && tag[0] != '\0') ? std::string(tag) : std::string("run");
}
bool FileExists(const std::string& path) {
    std::ifstream in(path);
    return in.good();
}
double Exec26ScintAirSurfaceLoss() {
    static const double value = [] {
        const char* raw = std::getenv("EXEC26_SCINT_AIR_SURFACE_LOSS");
        if (raw == nullptr || raw[0] == '\0') return 0.0;
        char* end = nullptr;
        const double parsed = std::strtod(raw, &end);
        if (end == raw || !std::isfinite(parsed)) return 0.0;
        return std::max(0.0, std::min(1.0, parsed));
    }();
    return value;
}
bool ApplyExec26ScintAirSurfaceLoss(const G4Step* step) {
    const double loss = Exec26ScintAirSurfaceLoss();
    if (loss <= 0.0) return false;

    auto* pre = step->GetPreStepPoint();
    auto* post = step->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary) return false;

    auto* preVol = pre->GetPhysicalVolume();
    auto* postVol = post->GetPhysicalVolume();
    if (preVol == nullptr || postVol == nullptr) return false;

    const G4String preName = preVol->GetLogicalVolume()->GetName();
    const G4String postName = postVol->GetLogicalVolume()->GetName();
    if (!IsBarLV(preName)) return false;

    auto* boundary = BoundaryProcess();
    const auto status = boundary ? boundary->GetStatus() : Undefined;
    const bool scintAirEncounter =
        IsAirGapLV(postName) ||
        (IsBarLV(postName) && IsReflectionStatus(status));
    if (!scintAirEncounter) return false;
    if (G4UniformRand() >= loss) return false;

    {
        std::lock_guard<std::mutex> lock(gAuditMutex);
        gBoundaryStatusCounts["BarLV->ScintAirSurfaceLoss|Killed"] += 1;
    }
    ++gExec26ScintAirSurfaceLossKilled;
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return true;
}
void AuditStep(const G4Step* step) {
    auto* pre = step->GetPreStepPoint();
    auto* post = step->GetPostStepPoint();
    auto* track = step->GetTrack();
    auto* preVol = pre->GetPhysicalVolume();
    auto* postVol = post->GetPhysicalVolume();
    auto* material = pre->GetMaterial();
    const G4String preName =
        preVol ? preVol->GetLogicalVolume()->GetName() : "NULL";
    const G4String postName =
        postVol ? postVol->GetLogicalVolume()->GetName() : "NULL";
    const G4String matName = material ? material->GetName() : "NULL";
    const int eventId = CurrentEventId();
    const long long key = TrackKey(eventId, track->GetTrackID());

    G4double nPhase = 0.0;
    bool hasRindex = false;
    if (material != nullptr) {
        auto* mpt = material->GetMaterialPropertiesTable();
        auto* rindex = mpt ? mpt->GetProperty("RINDEX") : nullptr;
        if (rindex != nullptr && rindex->GetVectorLength() > 0) {
            nPhase = rindex->Value(track->GetKineticEnergy());
            hasRindex = nPhase > 0.0;
        }
    }
    const double ratio = hasRindex
        ? track->GetVelocity() / (c_light / nPhase)
        : 0.0;

    std::lock_guard<std::mutex> lock(gAuditMutex);
    gVolumeSteps[std::string(preName)] += 1;
    auto& speed = gSpeedStats[std::string(preName) + "|" + std::string(matName)];
    speed.nSteps += 1;
    if (!hasRindex) {
        speed.nMissingRindex += 1;
    } else {
        speed.maxRatio = std::max(speed.maxRatio, ratio);
        if (ratio > 1.001) speed.nSuper += 1;
    }

    auto& history = gTrackHistory[key];
    if (IsAirGapLV(preName) || IsAirGapLV(postName)) history.visitedAirGap = true;

    if (post->GetStepStatus() == fGeomBoundary) {
        auto* boundary = BoundaryProcess();
        const auto status = boundary ? boundary->GetStatus() : Undefined;
        const std::string statusName = StatusName(status);
        const std::string boundaryName = std::string(preName) + "->" + std::string(postName);
        gBoundaryStatusCounts[boundaryName + "|" + statusName] += 1;

        if (IsBarLV(preName) && IsAirGapLV(postName)) ++gBarToMylar;
        if (IsBarLV(preName) && postName == "WorldLV") ++gMylarToWorld;
        if (IsBarLV(preName) && IsBarLV(postName)) ++gMylarReflected;
        if (IsBarLV(preName) && IsSiPMLV(postName)) ++gMylarToSiPM;

        if ((IsBarLV(preName) && IsAirGapLV(postName)) ||
            (IsAirGapLV(preName) && IsBarLV(postName))) {
            if (IsReflectionStatus(status)) history.nScintAirReflections += 1;
        }
        if ((IsAirGapLV(preName) && IsReflectorLV(postName)) ||
            (IsReflectorLV(preName) && IsAirGapLV(postName))) {
            history.visitedReflectorBoundary = true;
            if (IsReflectionStatus(status)) history.nAirReflectorReflections += 1;
        }
        if (status == TotalInternalReflection) {
            history.nTotalInternalReflections += 1;
        }
    }
}
} // namespace

namespace BoundaryCensus {
long long GetBarToMylar()     { return gBarToMylar.load(); }
long long GetMylarToWorld()   { return gMylarToWorld.load(); }
long long GetMylarReflected() { return gMylarReflected.load(); }
long long GetMylarToSiPM()    { return gMylarToSiPM.load(); }
long long GetKilledWorld()    { return gKilledWorld.load(); }
void Reset() {
    gBarToMylar = 0;
    gMylarToWorld = 0;
    gMylarReflected = 0;
    gMylarToSiPM = 0;
    gKilledWorld = 0;
    gExec26ScintAirSurfaceLossKilled = 0;
}
} // namespace BoundaryCensus

namespace OpticalAudit {
void Reset() {
    std::lock_guard<std::mutex> lock(gAuditMutex);
    gVolumeSteps.clear();
    gBoundaryStatusCounts.clear();
    gSpeedStats.clear();
    gTrackHistory.clear();
    gDetections.clear();
}

void RecordDetection(int eventId, int trackId, int globalId,
                     double timeNs, double gunXmm) {
    const long long key = TrackKey(eventId, trackId);
    std::lock_guard<std::mutex> lock(gAuditMutex);
    const auto found = gTrackHistory.find(key);
    const TrackHistory history = found != gTrackHistory.end()
        ? found->second
        : TrackHistory{};
    DetectionRow row;
    row.eventId = eventId;
    row.trackId = trackId;
    row.globalId = globalId;
    row.detectedEndSide = (globalId < 8) ? 0 : ((globalId < 16) ? 1 : 2);
    row.nScintAirReflections = history.nScintAirReflections;
    row.nAirReflectorReflections = history.nAirReflectorReflections;
    row.nTotalInternalReflections = history.nTotalInternalReflections;
    row.visitedAirGap = history.visitedAirGap;
    row.visitedReflectorBoundary = history.visitedReflectorBoundary;
    row.firstDetectionTimeNs = timeNs;
    row.gunXmm = gunXmm;
    gDetections.push_back(row);
}

void WriteOutputs() {
    const std::string dir = OutputDir();
    if (dir.empty()) return;
    const std::string tag = RunTag();
    std::lock_guard<std::mutex> lock(gAuditMutex);

    {
        const std::string path = dir + "/optical_boundary_counts.csv";
        const bool needHeader = !FileExists(path);
        std::ofstream out(path, std::ios::app);
        if (needHeader) out << "run_tag,boundary,status,count\n";
        for (const auto& item : gBoundaryStatusCounts) {
            const auto split = item.first.find('|');
            out << tag << ','
                << item.first.substr(0, split) << ','
                << item.first.substr(split + 1) << ','
                << item.second << '\n';
        }
    }
    {
        const std::string path = dir + "/optical_volume_counts.csv";
        const bool needHeader = !FileExists(path);
        std::ofstream out(path, std::ios::app);
        if (needHeader) out << "run_tag,volume,n_steps\n";
        for (const auto& item : gVolumeSteps) {
            out << tag << ',' << item.first << ',' << item.second << '\n';
        }
    }
    {
        const std::string path = dir + "/speed_audit_steps.csv";
        const bool needHeader = !FileExists(path);
        std::ofstream out(path, std::ios::app);
        if (needHeader) {
            out << "run_tag,volume,material,n_steps,max_ratio,"
                << "superluminal_fraction,missing_rindex\n";
        }
        for (const auto& item : gSpeedStats) {
            const auto split = item.first.find('|');
            const auto& stats = item.second;
            const double frac = stats.nSteps > 0
                ? static_cast<double>(stats.nSuper) / stats.nSteps
                : 0.0;
            out << tag << ','
                << item.first.substr(0, split) << ','
                << item.first.substr(split + 1) << ','
                << stats.nSteps << ','
                << std::setprecision(12) << stats.maxRatio << ','
                << frac << ','
                << stats.nMissingRindex << '\n';
        }
    }
    {
        const std::string path = dir + "/detected_photon_history.csv";
        const bool needHeader = !FileExists(path);
        std::ofstream out(path, std::ios::app);
        if (needHeader) {
            out << "run_tag,event_id,track_id,global_id,detected_end_side,"
                << "n_reflections_scint_air,n_reflections_air_reflector,"
                << "n_total_internal_reflections,visited_air_gap,"
                << "visited_reflector_boundary,first_detection_time_ns,gun_x_mm\n";
        }
        for (const auto& row : gDetections) {
            out << tag << ','
                << row.eventId << ','
                << row.trackId << ','
                << row.globalId << ','
                << row.detectedEndSide << ','
                << row.nScintAirReflections << ','
                << row.nAirReflectorReflections << ','
                << row.nTotalInternalReflections << ','
                << (row.visitedAirGap ? 1 : 0) << ','
                << (row.visitedReflectorBoundary ? 1 : 0) << ','
                << std::setprecision(12) << row.firstDetectionTimeNs << ','
                << row.gunXmm << '\n';
        }
    }
    {
        const std::string path = dir + "/" + tag + "_optical_history_summary.json";
        std::ofstream out(path);
        out << "{\n"
            << "  \"run_tag\": \"" << tag << "\",\n"
            << "  \"n_detected_rows\": " << gDetections.size() << ",\n"
            << "  \"n_tracked_histories\": " << gTrackHistory.size() << ",\n"
            << "  \"n_boundary_categories\": " << gBoundaryStatusCounts.size() << "\n"
            << "}\n";
    }
}
} // namespace OpticalAudit

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* track = step->GetTrack();
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;

    AuditStep(step);
    if (ApplyExec26ScintAirSurfaceLoss(step)) return;

    if (track->GetCurrentStepNumber() == 1) {
        const G4VProcess* creator = track->GetCreatorProcess();
        if (creator != nullptr && creator->GetProcessName() == "Scintillation") {
            auto* ra = dynamic_cast<RunAction*>(
                const_cast<G4UserRunAction*>(
                    G4RunManager::GetRunManager()->GetUserRunAction()));
            if (ra != nullptr) ra->AddScintPhoton();
        }
    }

    const G4double energy = track->GetKineticEnergy();
    if (energy > 0.0) {
        const G4double wl_nm = kHC_eVnm / (energy / eV);
        if (wl_nm < 300.0 || wl_nm > 900.0) {
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    if (postVol == nullptr) {
        track->SetTrackStatus(fStopAndKill);
        return;
    }
    if (postVol->GetLogicalVolume()->GetName() == "WorldLV") {
        ++gKilledWorld;
        track->SetTrackStatus(fStopAndKill);
    }
}
