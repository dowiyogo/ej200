#include "TrackingAction.hh"
#include "G4Event.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
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
std::map<G4int, long long> barBoundaryAll, barBoundaryScint;
using StatusKey = std::tuple<G4String, G4String, G4String, G4int>;
std::map<StatusKey, long long> terminalStatusAll, terminalStatusScint;
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
    barBoundaryAll.clear(); barBoundaryScint.clear();
    terminalStatusAll.clear(); terminalStatusScint.clear();
    for (const auto* key : {"started_optical", "started_scintillation", "terminal_optical",
             "terminal_scintillation", "duplicate_terminal", "nonterminal_post", "unknown_process",
             "unknown_volume", "terminal_without_start"}) audit[key] = 0;
}
void Write(const G4String& prefix) {
    const std::lock_guard<std::mutex> lock(censusMutex);
    WriteMap(prefix + "_all.csv", allFates);
    WriteMap(prefix + "_scintillation.csv", scintFates);
    const auto writeBoundary = [&prefix](const G4String& suffix, const std::map<G4int, long long>& counts) {
        std::ofstream out(prefix + suffix);
        out << "process,volume,kill_reason,status_int,count\n";
        for (const auto& item : counts)
            out << "Transportation,BarLV,none," << item.first << ',' << item.second << '\n';
        out.flush();
        if (!out) throw std::runtime_error("Cannot write terminal boundary states: " + prefix);
    };
    writeBoundary("_bar_boundary_all.csv", barBoundaryAll);
    writeBoundary("_bar_boundary_scintillation.csv", barBoundaryScint);
    const auto writeStatuses = [&prefix](const G4String& suffix, const std::map<StatusKey, long long>& counts) {
        std::ofstream out(prefix + suffix);
        out << "process,volume,kill_reason,status_int,count\n";
        for (const auto& item : counts)
            out << Quote(std::get<0>(item.first)) << ',' << Quote(std::get<1>(item.first)) << ','
                << Quote(std::get<2>(item.first)) << ',' << std::get<3>(item.first) << ',' << item.second << '\n';
        out.flush();
        if (!out) throw std::runtime_error("Cannot write terminal status ledger: " + prefix);
    };
    writeStatuses("_states_all.csv", terminalStatusAll);
    writeStatuses("_states_scintillation.csv", terminalStatusScint);
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
    // EXEC_30: current final-step status, never a stale previous boundary.
    static G4ThreadLocal G4OpBoundaryProcess* boundary = nullptr;
    if (!boundary) {
        auto* processes = track->GetDefinition()->GetProcessManager()->GetProcessList();
        for (G4int i = 0; i < static_cast<G4int>(processes->size()); ++i)
            if ((*processes)[i]->GetProcessName() == "OpBoundary") {
                boundary = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]); break;
            }
    }
    const G4int status = post && post->GetStepStatus() == fGeomBoundary
        ? (boundary ? static_cast<G4int>(boundary->GetStatus()) : static_cast<G4int>(Undefined))
        : static_cast<G4int>(NotAtBoundary);
    const StatusKey statusKey{processName, volumeName, reason, status};
    ++terminalStatusAll[statusKey];
    if (IsScint(track)) ++terminalStatusScint[statusKey];
    if (processName == "Transportation" && volumeName == "BarLV" && reason == "none") {
        ++barBoundaryAll[status];
        if (IsScint(track)) ++barBoundaryScint[status];
    }
    ++allFates[key]; ++audit["terminal_optical"];
    if (IsScint(track)) { ++scintFates[key]; ++audit["terminal_scintillation"]; }
    if (!process) ++audit["unknown_process"];
    if (volumeName == "unknown") ++audit["unknown_volume"];
    killFlags.erase(track->GetTrackID());
}
