#pragma once
#include "G4LogicalBorderSurface.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// --------------------------------------------------------------------------
// DetectorConstruction — EJ-228 cylinder (feat/ej228-cylinder)
//
// Volume hierarchy:
//   WorldPV  (G4_AIR, 200×200×200 mm³)
//     ├─ CylPV        EJ-228, r=12.5 mm, h=25 mm  (axis along Z)
//     ├─ AirGapPV     annular,  r=12.5–12.6 mm, h=25 mm (mantle only)
//     ├─ VikuitiPV    annular,  r=12.6–12.7 mm, h=25 mm (Vikuiti 3M ESR)
//     ├─ TopSiPM_PV × 4   z=+12.75 mm, 2×2 grid (global IDs 0–3)
//     └─ BotSiPM_PV × 4   z=−12.75 mm, 2×2 grid (global IDs 4–7)
//
// Optical surfaces:
//   CylPV → AirGapPV         polished dielectric-dielectric (TIR at r=12.5 mm)
//   AirGapPV → VikuitiPV     Vikuiti specular  R=0.98
//   CylPV → each SiPM PV     dielectric-metal, EFFICIENCY (AFBR-S4N66P024M PDE)
//
// face_type encoding (matches SiPMSD convention):
//   0 = top cap  (+Z, global IDs 0–3)
//   1 = bottom cap (−Z, global IDs 4–7)
//
// GROUPVEL fix: WorldLV air gets GROUPVEL = c/n_bar to eliminate the
// G4Transportation velocity-aliasing bug after TIR bounces.
// --------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    static constexpr G4int kNTopSiPMs    = 4;
    static constexpr G4int kNBottomSiPMs = 4;
    static constexpr G4int kNTotalSiPMs  = 8;

    DetectorConstruction()  = default;
    ~DetectorConstruction() = default;

    G4VPhysicalVolume* Construct()           override;
    void               ConstructSDandField() override;

    // face_type 0=top (+Z), 1=bottom (−Z)
    static G4int FaceType(G4int globalId) { return globalId < kNTopSiPMs ? 0 : 1; }
    static G4int LocalId (G4int globalId) {
        return globalId < kNTopSiPMs ? globalId : globalId - kNTopSiPMs;
    }

    // ── Interface stubs kept for RunAction compatibility ──────────────────
    G4Material* GetActiveScintillatorMaterial() const { return fScintMat; }
    G4String    GetScintillatorCode()    const { return "EJ-228"; }
    G4String    GetSiPMModel()           const { return fSiPMModel; }
    G4String    GetReadoutConfiguration()const { return "CylTopBot"; }
    G4bool      IsEndInstrumented()      const { return true;  }  // top+bot caps
    G4bool      IsTopInstrumented()      const { return false; }
    G4int       GetNActiveEndSiPMs()     const { return kNTotalSiPMs; }
    G4int       GetNActiveTopSiPMs()     const { return 0; }
    G4int       GetNTopSiPMs()           const { return 0; }

    // Kept so RunAction can iterate surfaces for diagnostics
    const std::map<G4int, G4LogicalBorderSurface*>& GetSiPMSurfaces() const
    { return fSiPMSurfaces; }
    const std::map<G4String, G4LogicalBorderSurface*>& GetReflectorSurfaces() const
    { return fReflectorSurfaces; }

    // Unused bar helpers — stubs so old call sites still compile
    static G4double TopSiPMCenterX(G4int) { return 0.0; }
    void SetScintillatorCode(G4String)     {}
    void SetSiPMModel(G4String)            {}
    void SetReadoutConfiguration(G4String) {}
    void SetEdgeWrapMode(G4String)         {}

  private:
    G4LogicalVolume* fSiPMLV  = nullptr;
    G4Material*      fScintMat = nullptr;
    G4String         fSiPMModel = "AFBR-S4N66P024M";

    std::map<G4int,    G4LogicalBorderSurface*> fSiPMSurfaces;
    std::map<G4String, G4LogicalBorderSurface*> fReflectorSurfaces;  // empty
};
