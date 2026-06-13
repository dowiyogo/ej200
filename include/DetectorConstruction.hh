#pragma once
#include "G4LogicalBorderSurface.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;
class G4GenericMessenger;

// --------------------------------------------------------------------------
// DetectorConstruction
//
// Volume hierarchy:
//   WorldPV (air)
//     ├─ BarPV (SSLG4 OPSC-106/EJ-230 by default — EXEC_13 campaign)
//     │   ├─ EndSiPMLeft_PV  × 8   (global IDs  0– 7)
//     │   ├─ EndSiPMRight_PV × 8   (global IDs  8–15)
//     └─ Reflector*PV panels on every non-END face
//
// SiPM coupling volumes are BarLV daughters with BarPV->SiPM border surfaces.
//
// UI commands:
//   /det/scintillator OPSC-106|OPSC-101|OPSC-100
//   /sipm/model AFBR-S4N66P024M
//   /det/edgeWrap mylar|air|black  — accepted for legacy macros; no-op
//
// Border surfaces: BarPV → each SiPM physical volume; BarPV → reflector panels.
// --------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    static constexpr G4int kNEndSiPMs = 8;   // per side (8×1 array)
    static constexpr G4int kNTotalSiPMs = 2 * kNEndSiPMs;

    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct()           override;
    void               ConstructSDandField() override;

    // Border-surface map (globalId → surface) — kept for external consumers
    const std::map<G4int, G4LogicalBorderSurface*>& GetSiPMSurfaces() const
    { return fSiPMSurfaces; }
    const std::map<G4String, G4LogicalBorderSurface*>& GetReflectorSurfaces() const
    { return fReflectorSurfaces; }

    // Active scintillator selector — triggers geometry rebuild
    void        SetScintillatorCode(G4String code);
    G4String    GetScintillatorCode() const { return fScintCode; }
    G4Material* GetActiveScintillatorMaterial() const { return fActiveScintillator; }

    void     SetSiPMModel(G4String model);
    G4String GetSiPMModel() const { return fSiPMModel; }

    G4int GetNActiveEndSiPMs() const { return 2 * kNEndSiPMs; }

    // Legacy edge-wrap command — retained as a no-op for old macros
    void SetEdgeWrapMode(G4String mode);

    // Helpers for analysis (independent of pitch/count)
    static G4int FaceType(G4int globalId);  // 0=end_left, 1=end_right
    static G4int LocalId (G4int globalId);  // index within face

  private:
    G4LogicalVolume*   fEndSiPMLV  = nullptr;
    G4VPhysicalVolume* fBarPhys          = nullptr;

    G4String           fScintCode = "OPSC-106";
    G4String           fSiPMModel = "AFBR-S4N66P024M";
    G4Material*         fActiveScintillator = nullptr;

    G4GenericMessenger* fMessenger = nullptr;
    G4GenericMessenger* fSiPMMessenger = nullptr;

    std::map<G4int, G4LogicalBorderSurface*> fSiPMSurfaces;
    std::map<G4String, G4LogicalBorderSurface*> fReflectorSurfaces;
};
