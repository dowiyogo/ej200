#include "DisplayTrackingAction.hh"

#include "G4GenericMessenger.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
thread_local int gConfiguredMaxOpticalTrajectories = 100;
thread_local int gOpticalOrdinal = 0;
thread_local int gOpticalSeen = 0;
thread_local int gOpticalStored = 0;
thread_local std::unordered_set<int> gDisplayedOpticalTrackIds;

struct BoundaryMark {
    int eventId = -1;
    int trackId = -1;
    int order = 0;
    G4OpBoundaryProcessStatus status = Undefined;
    double x_mm = 0.0;
    double y_mm = 0.0;
    double z_mm = 0.0;
};

thread_local std::unordered_map<int, std::vector<BoundaryMark>>
    gBoundaryMarksByTrack;

const char* BoundaryStatusToString(G4OpBoundaryProcessStatus status) {
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
        default: return "Unknown";
    }
}

const char* BoundaryClassLabel(G4OpBoundaryProcessStatus status) {
    switch (status) {
        case FresnelReflection:
        case TotalInternalReflection:
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
        case BackScattering:
        case PolishedLumirrorAirReflection:
        case PolishedLumirrorGlueReflection:
        case PolishedAirReflection:
        case PolishedTeflonAirReflection:
        case PolishedTiOAirReflection:
        case PolishedTyvekAirReflection:
        case PolishedVM2000AirReflection:
        case PolishedVM2000GlueReflection:
        case EtchedLumirrorAirReflection:
        case EtchedLumirrorGlueReflection:
        case EtchedAirReflection:
        case EtchedTeflonAirReflection:
        case EtchedTiOAirReflection:
        case EtchedTyvekAirReflection:
        case EtchedVM2000AirReflection:
        case EtchedVM2000GlueReflection:
        case GroundLumirrorAirReflection:
        case GroundLumirrorGlueReflection:
        case GroundAirReflection:
        case GroundTeflonAirReflection:
        case GroundTiOAirReflection:
        case GroundTyvekAirReflection:
        case GroundVM2000AirReflection:
        case GroundVM2000GlueReflection:
        case CoatedDielectricReflection:
            return "reflection";
        case Transmission:
        case FresnelRefraction:
        case CoatedDielectricRefraction:
        case CoatedDielectricFrustratedTransmission:
            return "refraction_or_transmission";
        case Absorption:
        case Detection:
            return "termination";
        default:
            return "other";
    }
}
}  // namespace

DisplayTrackingAction::DisplayTrackingAction() {
    fMessenger = new G4GenericMessenger(this, "/display/",
                                        "Display trajectory controls");

    auto& maxOptCmd =
        fMessenger->DeclareProperty("maxOpticalTrajectories",
                                    fMaxOpticalTrajectories,
                                    "Maximum optical photon trajectories to store per event");
    maxOptCmd.SetParameterName("max", false);
    maxOptCmd.SetRange("max>=0");

    fMessenger->DeclareProperty("storeMuon", fStoreMuon,
                                "Store primary muon trajectory");
    fMessenger->DeclareProperty("storeOtherParticles", fStoreOtherParticles,
                                "Store non-muon/non-optical trajectories");
}

DisplayTrackingAction::~DisplayTrackingAction() {
    delete fMessenger;
}

void DisplayTrackingAction::PreUserTrackingAction(const G4Track* track) {
    if (track == nullptr || fpTrackingManager == nullptr) return;

    gConfiguredMaxOpticalTrajectories = std::max(0, fMaxOpticalTrajectories);

    const auto* particle = track->GetParticleDefinition();
    const int trackId = track->GetTrackID();
    const bool isOptical = (particle == G4OpticalPhoton::Definition());
    const bool isMuon = (std::abs(particle->GetPDGEncoding()) == 13);

    bool storeTrajectory = false;
    if (isOptical) {
        ++gOpticalSeen;
        ++gOpticalOrdinal;
        if (gOpticalOrdinal <= gConfiguredMaxOpticalTrajectories) {
            storeTrajectory = true;
            ++gOpticalStored;
            gDisplayedOpticalTrackIds.insert(trackId);
        }
    }
    else if (isMuon) {
        storeTrajectory = fStoreMuon;
    }
    else {
        storeTrajectory = fStoreOtherParticles;
    }

    fpTrackingManager->SetStoreTrajectory(storeTrajectory ? 1 : 0);
}

namespace DisplayDiagnostics {
void BeginEvent() {
    gOpticalOrdinal = 0;
    gOpticalSeen = 0;
    gOpticalStored = 0;
    gDisplayedOpticalTrackIds.clear();
    gBoundaryMarksByTrack.clear();
}

void EndEvent(int eventId) {
    std::vector<BoundaryMark> marks;
    marks.reserve(gBoundaryMarksByTrack.size() * 2U);

    for (auto& kv : gBoundaryMarksByTrack) {
        for (auto& mark : kv.second) {
            mark.eventId = eventId;
            marks.push_back(mark);
        }
    }

    std::sort(marks.begin(), marks.end(), [](const BoundaryMark& a,
                                             const BoundaryMark& b) {
        if (a.trackId != b.trackId) return a.trackId < b.trackId;
        return a.order < b.order;
    });

    char path[128];
    std::snprintf(path, sizeof(path), "scatter_points_event%03d.csv", eventId);

    std::ofstream out(path);
    if (!out) {
        G4cerr << "[DisplayDiagnostics] Could not write " << path << G4endl;
        return;
    }

    out << "event_id,track_id,interaction_order,status,status_class,marker_label,x_mm,y_mm,z_mm\n";
    for (const auto& mark : marks) {
        out << mark.eventId << ','
            << mark.trackId << ','
            << mark.order << ','
            << BoundaryStatusToString(mark.status) << ','
            << BoundaryClassLabel(mark.status) << ','
            << (mark.order == 1 ? "scatter1" : "scatter2") << ','
            << std::fixed << std::setprecision(6)
            << mark.x_mm << ','
            << mark.y_mm << ','
            << mark.z_mm << '\n';
    }

    G4cout << "[DisplayDiagnostics] Wrote " << path
           << " with " << marks.size() << " marker rows" << G4endl;
}

void RecordBoundaryStatus(const G4Track* track,
                          G4OpBoundaryProcessStatus status,
                          const G4ThreeVector& position) {
    if (track == nullptr) return;
    const int trackId = track->GetTrackID();
    if (gDisplayedOpticalTrackIds.find(trackId) == gDisplayedOpticalTrackIds.end()) {
        return;
    }

    auto& marks = gBoundaryMarksByTrack[trackId];
    if (marks.size() >= 2U) return;

    BoundaryMark mark;
    mark.trackId = trackId;
    mark.order = static_cast<int>(marks.size()) + 1;
    mark.status = status;
    mark.x_mm = position.x() / mm;
    mark.y_mm = position.y() / mm;
    mark.z_mm = position.z() / mm;
    marks.push_back(mark);
}

int GetConfiguredMaxOpticalTrajectories() {
    return gConfiguredMaxOpticalTrajectories;
}

int GetOpticalTracksSeen() {
    return gOpticalSeen;
}

int GetOpticalTrajectoriesStored() {
    return gOpticalStored;
}
}  // namespace DisplayDiagnostics
