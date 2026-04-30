#pragma once
#include "G4LogicalBorderSurface.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GenericMessenger;

// --------------------------------------------------------------------------
// DetectorConstruction
//
// Geometry: 1.4 m × 60 mm × 10 mm EJ-200/EJ-230 bar wrapped in a thin Mylar film
// (25 µm, n=1.65) implemented as a nested physical volume.  The Mylar layer
// provides reflection via TIR at the Mylar–air interface without requiring a
// G4OpticalSurface on the outer boundary, improving simulation speed.
//
// Volume hierarchy:
//   WorldPV (air)
//     ├─ WrapCapLeftPV   (configurable edge-cap material)
//     │    └─ BarCapLeftPV
//     ├─ WrapCenterPV    (Mylar)
//     │    └─ BarCenterPV
//     └─ WrapCapRightPV  (configurable edge-cap material)
//          └─ BarCapRightPV
//          (bar material is EJ-230 by default; EJ-200 still selectable)
//     ├─ EndSiPMLeft_PV  × 8   (global IDs  0– 7)
//     ├─ EndSiPMRight_PV × 8   (global IDs  8–15)
//     └─ TopSiPMPV       × N   (global IDs 16…15+N)
//
// N (number of top SiPMs) is computed from the configurable pitch so that all
// SiPMs remain inside the bar footprint.  Default pitch = 70 mm → N = 20.
//
// UI command (available after /run/initialize):
//   /det/topSiPMPitch <value> mm   — triggers ReinitializeGeometry()
//   /det/scintillatorMaterial EJ230|EJ200
//   /det/edgeWrap mylar|air|black  — edge cap material, default mylar
//
// Border surfaces: the adjacent wrap segment → each SiPM physical volume
// (not BarPV → SiPM, because the wrap layer is geometrically between bar and SiPM).
// --------------------------------------------------------------------------
enum class EdgeWrapMode { Mylar, Air, BlackTape };

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

    // X-centre of top SiPM with given index (0-based) for a given pitch [G4 units]
    static G4double TopSiPMCenterX(G4int idx, G4double pitch, G4int nTotal);

    // Configurable top-SiPM pitch — triggers geometry rebuild
    void     SetTopSiPMPitch(G4double pitchMm);  // receives value in mm from macro
    G4double GetTopSiPMPitch() const { return fTopSiPMPitch; }
    G4int    GetNTopSiPMs()    const { return fNTopSiPMs; }

    // Scintillator material selector. Default is EJ230 for the current study.
    void     SetScintillatorMaterial(G4String materialName);
    G4String GetScintillatorMaterial() const { return fScintillatorName; }

    // Configurable edge wrapping material — triggers geometry rebuild
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
    G4VPhysicalVolume* fWrapPhys         = nullptr;  // central Mylar envelope alias
    G4VPhysicalVolume* fWrapCenterPhys   = nullptr;
    G4VPhysicalVolume* fWrapCapLeftPhys  = nullptr;
    G4VPhysicalVolume* fWrapCapRightPhys = nullptr;

    G4double           fTopSiPMPitch = 70.0;  // G4 internal units (1 = 1 mm)
    G4int              fNTopSiPMs    = 20;
    G4String           fScintillatorName = "EJ230";
    EdgeWrapMode       fEdgeWrapMode = EdgeWrapMode::Mylar;

    G4GenericMessenger* fMessenger = nullptr;

    std::map<G4int, G4LogicalBorderSurface*> fSiPMSurfaces;
};
