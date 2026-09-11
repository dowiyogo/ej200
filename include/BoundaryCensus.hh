#pragma once
#include "globals.hh"
#include <map>
#include <mutex>
#include <tuple>

namespace BoundaryCensus {
// Clave direccional: volumen físico y copia a cada lado, más estado observado.
struct BoundaryKey {
    G4String pre;
    G4int pre_copy;
    G4String post;
    G4int post_copy;
    G4int status;
    bool operator<(const BoundaryKey& other) const {
        return std::tie(pre, pre_copy, post, post_copy, status) <
               std::tie(other.pre, other.pre_copy, other.post, other.post_copy, other.status);
    }
};

// Un único contador compartido; el mutex protege la acumulación entre hilos.
// La API de los contadores históricos permanece en el mismo espacio de nombres.
class Census {
  public:
    void Record(const BoundaryKey& key);
    void Reset();
    void Write(const G4String& filename) const;
    Census(const Census&) = delete;
    Census& operator=(const Census&) = delete;
  private:
    Census() = default;
    friend Census& Instance();
    mutable std::mutex mutex_;
    std::map<BoundaryKey, G4long> counts_;
};
Census& Instance();
} // namespace BoundaryCensus
