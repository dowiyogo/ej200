#pragma once
#include "G4LogicalBorderSurface.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>

class G4LogicalVolume;
class G4VPhysicalVolume;

// --------------------------------------------------------------------------
// DetectorConstruction
//
// Geometry: 1.4 m × 60 mm × 10 mm EJ-200 bar wrapped in reflective material.
//
// SiPMs:
//   End SiPMs — two 8×1 arrays (8 SiPMs, 6×6 mm² each) on the ±X end faces.
//               global IDs  0– 7  (left, –X side)
//                           8–15  (right, +X side)
//
//   Top SiPMs — kNTopSiPMs individual 6×6 mm² SiPMs on the +Y face,
//               distributed uniformly along the FULL bar length.
//               global IDs 16 … (15 + kNTopSiPMs)
//
// Surfaces:
//   bar → world:    reflective wrapping (dielectric_metal, R=0.95)
//   bar → sipmPhys: optical coupling + DETECTIONEFFICIENCY table
//
// The SiPM material (CreateSiPMCoupling, n=1.58) matches the bar refractive
// index so there is no TIR at the bar–SiPM interface.
// --------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    static constexpr G4int kNEndSiPMs   = 8;    // per side (8×1 array)
    static constexpr G4int kNTopSiPMs   = 20;   // along full bar top face
    static constexpr G4int kNSiPMsTotal = 2 * kNEndSiPMs + kNTopSiPMs; // 36

    DetectorConstruction()  = default;
    ~DetectorConstruction() = default;

    G4VPhysicalVolume* Construct()        override;
    void               ConstructSDandField() override;

    // Border-surface map (globalId → surface) — consumed by SiPMSD
    const std::map<G4int, G4LogicalBorderSurface*>& GetSiPMSurfaces() const
    { return fSiPMSurfaces; }

    // X-centre of top SiPM with given index (0 … kNTopSiPMs-1)
    static G4double TopSiPMCenterX(G4int idx);

    // Helpers for analysis
    static G4int FaceType(G4int globalId);  // 0=end_left, 1=end_right, 2=top
    static G4int LocalId(G4int globalId);   // index within face

  private:
    G4LogicalVolume*  fEndSiPMLV = nullptr; // shared LV for all 16 end SiPMs
    G4LogicalVolume*  fTopSiPMLV = nullptr; // shared LV for all 20 top SiPMs
    G4VPhysicalVolume* fBarPhys  = nullptr;

    std::map<G4int, G4LogicalBorderSurface*> fSiPMSurfaces; // built in Construct()
};
