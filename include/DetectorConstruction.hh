#pragma once
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SystemOfUnits.hh"
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
//     ├─ BarPV (SSLG4 OPSC-101/EJ-204 by default)
//     │   ├─ EndSiPMLeft_PV  × 8   (global IDs  0– 7)
//     │   ├─ EndSiPMRight_PV × 8   (global IDs  8–15)
//     │   └─ TopSiPMPV       × 70  (global IDs 16–85, sipm mode only)
//     ├─ AirGap*PV       × 4   (long lateral faces only, END faces open)
//     └─ Mylar*PV        × 4   (external reflector panels behind air gap)
//
// In mylar mode (default), explicit air-gap and Mylar panels wrap only the four
// long lateral faces.
// In sipm mode, the skin uses Materials::CreateBarSkinReflector.
// SiPM coupling volumes are BarLV daughters with explicit BarPV→SiPM border surfaces.
//
// UI commands:
//   /det/scintillator OPSC-101|OPSC-100
//   /det/readout End|Top|EndTop
//   /sipm/model AFBR-S4N66P024M
//   /ship/geom/topSurface mylar|sipm
//   /ship/geom/mylar/reflectivity <0..1>
//   /ship/geom/mylar/specularLobe <0..1>
//   /ship/geom/mylar/sigmaAlpha <angle>
//   /ship/geom/mylar/finish polished|ground
//   /ship/geom/mylar/airGap <length>
//   /ship/geom/sipmEfficiencyMode nominal|zero|one
//   /det/edgeWrap mylar|air|black  — accepted for legacy macros; no-op
//
// Border surfaces: BarPV → each SiPM physical volume.
// Reflector: G4LogicalSkinSurface on BarLV.
// --------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    static constexpr G4int kNEndSiPMs = 8;   // per side (8 single-element placements)
    static constexpr G4int kNTopSiPMs = 70;
    static constexpr G4int kNTotalSiPMs = 2 * kNEndSiPMs + kNTopSiPMs;

    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct()           override;
    void               ConstructSDandField() override;

    // Surface accessors
    const std::map<G4int, G4LogicalBorderSurface*>& GetSiPMSurfaces() const
    { return fSiPMSurfaces; }
    G4LogicalSkinSurface* GetBarSkinSurface() const { return fBarSkinSurface; }
    const std::map<G4String, G4LogicalBorderSurface*>& GetReflectorSurfaces() const
    { return fReflectorSurfaces; }
    const std::map<G4String, G4LogicalBorderSurface*>& GetScintillatorAirSurfaces() const
    { return fScintillatorAirSurfaces; }

    // X-centre of top SiPM local index 0..69 [G4 units].
    static G4double TopSiPMCenterX(G4int idx);
    G4int GetNTopSiPMs() const { return kNTopSiPMs; }

    // Active scintillator selector — triggers geometry rebuild
    void        SetScintillatorCode(G4String code);
    G4String    GetScintillatorCode() const { return fScintCode; }
    G4Material* GetActiveScintillatorMaterial() const { return fActiveScintillator; }

    void     SetSiPMModel(G4String model);
    G4String GetSiPMModel() const { return fSiPMModel; }

    // Readout/wrapping selector — triggers geometry rebuild
    void     SetReadoutConfiguration(G4String configuration);
    G4String GetReadoutConfiguration() const { return fReadoutConfiguration; }
    G4bool   IsEndInstrumented() const;
    G4bool   IsTopInstrumented() const;
    G4int    GetNActiveEndSiPMs() const
    { return IsEndInstrumented() ? 2 * kNEndSiPMs : 0; }
    G4int    GetNActiveTopSiPMs() const
    { return IsTopInstrumented() ? kNTopSiPMs : 0; }

    void     SetTopSurface(G4String mode);
    G4String GetTopSurface() const { return fTopSurface; }
    void     SetMylarReflectivity(G4double value);
    void     SetMylarSpecularLobe(G4double value);
    void     SetMylarSigmaAlpha(G4double value);
    void     SetMylarFinish(G4String value);
    void     SetMylarAirGap(G4double value);
    void     SetSiPMEfficiencyMode(G4String value);
    G4double GetMylarReflectivity() const { return fMylarReflectivity; }
    G4double GetMylarSpecularLobe() const { return fMylarSpecularLobe; }
    G4double GetMylarSigmaAlpha() const { return fMylarSigmaAlpha; }
    G4String GetMylarFinish() const { return fMylarFinish; }
    G4double GetMylarAirGap() const { return fMylarAirGap; }
    G4String GetSiPMEfficiencyMode() const { return fSiPMEfficiencyMode; }

    // Legacy edge-wrap command — retained as a no-op for old macros
    void SetEdgeWrapMode(G4String mode);

    // Helpers for analysis (independent of pitch/count)
    static G4int FaceType(G4int globalId);  // 0=end_left, 1=end_right, 2=top
    static G4int LocalId (G4int globalId);  // index within face

  private:
    G4LogicalVolume*    fEndSiPMLV      = nullptr;
    G4LogicalVolume*    fTopSiPMLV      = nullptr;
    G4VPhysicalVolume*  fBarPhys        = nullptr;
    G4LogicalSkinSurface* fBarSkinSurface = nullptr;

    G4String           fScintCode = "OPSC-101";
    G4String           fSiPMModel = "AFBR-S4N66P024M";
    G4String           fReadoutConfiguration = "End";
    G4String           fTopSurface = "mylar";
    G4String           fSiPMEfficiencyMode = "nominal";
    G4double           fMylarReflectivity = 0.95;
    G4double           fMylarSpecularLobe = 1.0;
    G4double           fMylarSigmaAlpha = 0.1 * CLHEP::deg;
    G4String           fMylarFinish = "polished";
    G4double           fMylarAirGap = 0.10 * CLHEP::mm;
    G4Material*         fActiveScintillator = nullptr;

    G4GenericMessenger* fMessenger = nullptr;
    G4GenericMessenger* fSiPMMessenger = nullptr;
    G4GenericMessenger* fGeometryMessenger = nullptr;
    G4GenericMessenger* fMylarMessenger = nullptr;

    std::map<G4int, G4LogicalBorderSurface*> fSiPMSurfaces;
    std::map<G4String, G4LogicalBorderSurface*> fScintillatorAirSurfaces;
    std::map<G4String, G4LogicalBorderSurface*> fReflectorSurfaces;
};
