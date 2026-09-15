#include "PhotonTrackInfo.hh"

#include "G4Exception.hh"

#include <limits>

namespace {
constexpr G4int kMaximumBoundaryEncounters =
    std::numeric_limits<G4int>::max();
}

void PhotonTrackInfo::AddBoundaryEncounter() {
    if (fBoundaryEncounters < 0 ||
        fBoundaryEncounters == kMaximumBoundaryEncounters) {
        G4ExceptionDescription message;
        message << "Invalid boundary-encounter counter value "
                << fBoundaryEncounters << ".";
        G4Exception("PhotonTrackInfo::AddBoundaryEncounter",
                    "EXEC46_INVALID_BOUNDARY_COUNTER", FatalException, message);
    }
    ++fBoundaryEncounters;
}
