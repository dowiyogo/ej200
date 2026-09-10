#include "TrackingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <tuple>

namespace {
using Key = std::tuple<G4String, G4String, G4String>;
std::mutex censusMutex;
std::map<Key, long long> allFates, scintFates;
std::map<G4String, long long> audit;
// Event/track identities deduplicate terminal callbacks without changing track data.
thread_local G4int eventId = -1, runId = -1;
thread_local std::set<G4int> started, finished;
thread_local std::map<G4int, G4String> killFlags;
bool IsScint(const G4Track* t) {
    return t->GetCreatorProcess() && t->GetCreatorProcess()->GetProcessName() == "Scintillation";
}
G4String Quote(const G4String& s) {
    G4String r = "\"";
    for (char c : s) { if (c == '"') r += '"'; r += c; }
    return r + "\"";
}
void WriteMap(const G4String& path, const std::map<Key, long long>& counts) {
    std::ofstream out(path);
    out << "process,volume,kill_reason,count\n";
    for (const auto& entry : counts)
        out << Quote(std::get<0>(entry.first)) << ',' << Quote(std::get<1>(entry.first))
            << ',' << Quote(std::get<2>(entry.first)) << ',' << entry.second << '\n';
    out.flush();
    if (!out) throw std::runtime_error("Cannot write terminal census: " + path);
}
}
namespace TerminalCensus {
void MarkKill(const G4Track* track, const G4String& reason) {
    killFlags[track->GetTrackID()] = reason;
}
void Reset() {
    const std::lock_guard<std::mutex> lock(censusMutex);
    allFates.clear(); scintFates.clear(); audit.clear();
    for (const auto* key : {"started_optical", "started_scintillation", "terminal_optical",
             "terminal_scintillation", "duplicate_terminal", "nonterminal_post", "unknown_process",
             "unknown_volume", "terminal_without_start"}) audit[key] = 0;
}
void Write(const G4String& prefix) {
    const std::lock_guard<std::mutex> lock(censusMutex);
    WriteMap(prefix + "_all.csv", allFates);
    WriteMap(prefix + "_scintillation.csv", scintFates);
    std::ofstream out(prefix + "_audit.csv");
    out << "metric,count\n";
    for (const auto& entry : audit) out << entry.first << ',' << entry.second << '\n';
    out.flush();
    if (!out) throw std::runtime_error("Cannot write terminal audit: " + prefix);
}
}
void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;
    const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    const G4int id = event ? event->GetEventID() : -1;
    const auto* run = G4RunManager::GetRunManager()->GetCurrentRun();
    const G4int rid = run ? run->GetRunID() : -1;
    if (id != eventId || rid != runId) {
        eventId = id; runId = rid; started.clear(); finished.clear(); killFlags.clear();
    }
    if (started.insert(track->GetTrackID()).second) {
        killFlags.erase(track->GetTrackID());
        const std::lock_guard<std::mutex> lock(censusMutex);
        ++audit["started_optical"];
        if (IsScint(track)) ++audit["started_scintillation"];
    }
}
void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;
    const std::lock_guard<std::mutex> lock(censusMutex);
    if (track->GetTrackStatus() != fStopAndKill && track->GetTrackStatus() != fKillTrackAndSecondaries) {
        ++audit["nonterminal_post"]; return;
    }
    if (!finished.insert(track->GetTrackID()).second) { ++audit["duplicate_terminal"]; return; }
    if (!started.count(track->GetTrackID())) ++audit["terminal_without_start"];
    const auto* step = track->GetStep();
    const auto* post = step ? step->GetPostStepPoint() : nullptr;
    const auto* process = post ? post->GetProcessDefinedStep() : nullptr;
    const auto* pv = track->GetVolume();
    const G4String processName = process ? process->GetProcessName() : "unknown";
    const G4String volumeName = pv && pv->GetLogicalVolume() ? pv->GetLogicalVolume()->GetName() : "unknown";
    const auto flag = killFlags.find(track->GetTrackID());
    const G4String reason = flag == killFlags.end() ? "none" : flag->second;
    const Key key{processName, volumeName, reason};
    ++allFates[key]; ++audit["terminal_optical"];
    if (IsScint(track)) { ++scintFates[key]; ++audit["terminal_scintillation"]; }
    if (!process) ++audit["unknown_process"];
    if (volumeName == "unknown") ++audit["unknown_volume"];
    killFlags.erase(track->GetTrackID());
}
