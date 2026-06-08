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
//     ├─ BarPV (selectable EJ-204/EJ-200/EJ-230; default EJ-204)
//     │   ├─ EndSiPMLeft_PV  × 8   (global IDs  0– 7)
//     │   ├─ EndSiPMRight_PV × 8   (global IDs  8–15)
//     │   └─ TopSiPMPV       × N   (global IDs 16…15+N)
//     └─ Reflector*PV panels (optical border surfaces, R=0.98)
//
// The Mylar wrap volume was removed in fix/geometry-bar-in-world.
// Optical reflection is now handled by explicit reflector panels on the
// non-SiPM faces; SiPMs are BarLV daughters with BarPV→SiPM border surfaces.
//
// N (number of top SiPMs) is computed from the configurable pitch so that all
// SiPMs remain inside the bar footprint.  Default pitch = 70 mm → N = 20.
//
// UI command (available after /run/initialize):
//   /det/topSiPMPitch <value> mm   — triggers ReinitializeGeometry()
//   /det/scintillatorMaterial EJ204|EJ200|EJ230
//   /det/readout End|Top|EndTop
//   /det/edgeWrap mylar|air|black  — accepted for legacy macros; no-op
//
// Border surfaces: BarPV → each SiPM physical volume; BarPV → reflector panels.
// --------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    static constexpr G4int kNEndSiPMs = 8;   // per side (8×1 array)

    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct()           override;
    void               ConstructSDandField() override;

    // Border-surface map (globalId → surface) — kept for external consumers
    const std::map<G4int, G4LogicalBorderSurface*>& GetSiPMSurfaces() const
    { return fSiPMSurfaces; }
    const std::map<G4String, G4LogicalBorderSurface*>& GetReflectorSurfaces() const
    { return fReflectorSurfaces; }

    // X-centre of top SiPM with given index (0-based) for a given pitch [G4 units]
    static G4double TopSiPMCenterX(G4int idx, G4double pitch, G4int nTotal);

    // Configurable top-SiPM pitch — triggers geometry rebuild
    void     SetTopSiPMPitch(G4double pitchMm);  // receives value in mm from macro
    G4double GetTopSiPMPitch() const { return fTopSiPMPitch; }
    G4int    GetNTopSiPMs()    const { return fNTopSiPMs; }

    // Active scintillator selector — triggers geometry rebuild
    void        SetScintillatorMaterial(G4String materialName);
    G4String    GetScintillatorMaterial() const { return fScintillatorName; }
    G4Material* GetActiveScintillatorMaterial() const { return fActiveScintillator; }

    // Readout/wrapping selector — triggers geometry rebuild
    void     SetReadoutConfiguration(G4String configuration);
    G4String GetReadoutConfiguration() const { return fReadoutConfiguration; }
    G4bool   IsEndInstrumented() const;
    G4bool   IsTopInstrumented() const;
    G4int    GetNActiveEndSiPMs() const
    { return IsEndInstrumented() ? 2 * kNEndSiPMs : 0; }
    G4int    GetNActiveTopSiPMs() const
    { return IsTopInstrumented() ? fNTopSiPMs : 0; }

    // Legacy edge-wrap command — retained as a no-op for old macros
    void SetEdgeWrapMode(G4String mode);

    // Helpers for analysis (independent of pitch/count)
    static G4int FaceType(G4int globalId);  // 0=end_left, 1=end_right, 2=top
    static G4int LocalId (G4int globalId);  // index within face

  private:
    // Compute how many top SiPMs fit inside the bar at the requested pitch
    static G4int ComputeNTopSiPMs(G4double pitch);

    G4LogicalVolume*   fEndSiPMLV  = nullptr;
    G4LogicalVolume*   fTopSiPMLV  = nullptr;
    G4VPhysicalVolume* fBarPhys          = nullptr;

    G4double           fTopSiPMPitch = 70.0;  // G4 internal units (1 = 1 mm)
    G4int              fNTopSiPMs    = 20;
    G4String           fScintillatorName = "EJ204";
    G4String           fReadoutConfiguration = "End";
    G4Material*         fActiveScintillator = nullptr;

    G4GenericMessenger* fMessenger = nullptr;

    std::map<G4int, G4LogicalBorderSurface*> fSiPMSurfaces;
    std::map<G4String, G4LogicalBorderSurface*> fReflectorSurfaces;
};
