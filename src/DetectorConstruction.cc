#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalSkinSurface.hh"

// ── Geometry constants ───────────────────────────────────────────────────────
static constexpr G4double kBarHalfX = 700.0 * mm;  // 1.4 m total length
static constexpr G4double kBarHalfY =  30.0 * mm;  // 60 mm width
static constexpr G4double kBarHalfZ =   5.0 * mm;  // 10 mm height

// End SiPM (on ±X face): thin slab, 6×6 mm² active area
static constexpr G4double kEndHalfX = 0.25 * mm;  // thickness
static constexpr G4double kEndHalfY = 3.0  * mm;  // 6 mm in Y
static constexpr G4double kEndHalfZ = 3.0  * mm;  // 6 mm in Z
static constexpr G4double kEndPitch = 2*kEndHalfY + 1.5*mm; // 7.5 mm pitch

// Top SiPM (on +Y face): thin slab, 6×6 mm² active area
static constexpr G4double kTopHalfX = 3.0  * mm;  // 6 mm in X
static constexpr G4double kTopHalfY = 0.25 * mm;  // thickness
static constexpr G4double kTopHalfZ = 3.0  * mm;  // 6 mm in Z

// ── Static helpers ────────────────────────────────────────────────────────────
G4double DetectorConstruction::TopSiPMCenterX(G4int idx) {
    // 20 SiPMs uniformly from x = -665 mm to x = +665 mm, step = 70 mm.
    // Span: ±665 mm  →  well inside the ±700 mm bar.
    return (-665.0 + idx * 70.0) * mm;
}

G4int DetectorConstruction::FaceType(G4int globalId) {
    if (globalId < kNEndSiPMs)           return 0; // end left
    if (globalId < 2 * kNEndSiPMs)       return 1; // end right
    return 2;                                       // top
}

G4int DetectorConstruction::LocalId(G4int globalId) {
    if (globalId < kNEndSiPMs)     return globalId;
    if (globalId < 2*kNEndSiPMs)   return globalId - kNEndSiPMs;
    return globalId - 2*kNEndSiPMs;
}

// ── Construct ─────────────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();

    // ── Materials ────────────────────────────────────────────────────────────
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* barMat   = Materials::CreateEJ200();
    G4Material* sipmMat  = Materials::CreateSiPMCoupling();

    // World air must have RINDEX for optical photon tracking
    {
        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", {2.0*eV, 4.0*eV}, {1.0, 1.0});
        worldMat->SetMaterialPropertiesTable(mpt);
    }

    // ── Optical surfaces ─────────────────────────────────────────────────────
    auto* barSurface = Materials::CreateBarSurface();
    auto* sipmSurface = Materials::CreateSiPMSurface();

    // ── World volume ─────────────────────────────────────────────────────────
    auto* worldSolid = new G4Box("WorldSolid", 1.6*m, 0.25*m, 0.25*m);
    auto* worldLV    = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPhys  = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                         nullptr, false, 0, true);

    // ── Scintillating bar ────────────────────────────────────────────────────
    auto* barSolid = new G4Box("BarSolid", kBarHalfX, kBarHalfY, kBarHalfZ);
    auto* barLV    = new G4LogicalVolume(barSolid, barMat, "BarLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.1, 0.3, 1.0, 0.25));
        va->SetForceSolid(true);
        barLV->SetVisAttributes(va);
    }
    fBarPhys = new G4PVPlacement(nullptr, {}, barLV, "BarPV", worldLV,
                                 false, 0, true);

    // Reflective wrapping: bar → world (applies to all bar–world interfaces
    // except where SiPM volumes are placed, since those have their own border
    // surfaces with higher priority)
    new G4LogicalSkinSurface("BarSkin", barLV, barSurface);

    // ── End SiPMs — shared logical volume, 16 placements ─────────────────────
    // 8×1 array on each ±X face; 8 SiPMs stacked along Y.
    // Centres: y = −26.25, −18.75, −11.25, −3.75, +3.75, +11.25, +18.75, +26.25 mm
    // (span 58.5 mm centred on y=0; fits within ±30 mm bar width)

    auto* endSolid = new G4Box("EndSiPMSolid", kEndHalfX, kEndHalfY, kEndHalfZ);
    fEndSiPMLV     = new G4LogicalVolume(endSolid, sipmMat, "EndSiPMLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.0, 0.85, 0.2, 0.6));
        va->SetForceSolid(true);
        fEndSiPMLV->SetVisAttributes(va);
    }

    for (G4int i = 0; i < kNEndSiPMs; ++i) {
        // Centres: (i - 3.5) * 7.5 mm → −26.25, −18.75 … +18.75, +26.25 mm
        // Max edge = 26.25 + 3.0 = 29.25 mm  <  bar half-Y = 30 mm  ✓
        const G4double cy = (i - 3.5) * kEndPitch;

        // Left side (–X): global IDs 0..7, copy numbers 0..7
        const G4int leftId = i;
        auto* leftPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(-(kBarHalfX + kEndHalfX), cy, 0.0),
            fEndSiPMLV, "EndSiPMLeft_PV", worldLV, false, leftId, true);
        fSiPMSurfaces[leftId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(leftId),
            fBarPhys, leftPhys, sipmSurface);

        // Right side (+X): global IDs 8..15, copy numbers 8..15
        const G4int rightId = i + kNEndSiPMs;
        auto* rightPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(+(kBarHalfX + kEndHalfX), cy, 0.0),
            fEndSiPMLV, "EndSiPMRight_PV", worldLV, false, rightId, true);
        fSiPMSurfaces[rightId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(rightId),
            fBarPhys, rightPhys, sipmSurface);
    }

    // ── Top SiPMs — shared logical volume, kNTopSiPMs placements ─────────────
    // Uniformly distributed along the FULL bar length (both halves).
    // Each is a 6×6 mm² slab flush against the +Y bar face.

    auto* topSolid = new G4Box("TopSiPMSolid", kTopHalfX, kTopHalfY, kTopHalfZ);
    fTopSiPMLV     = new G4LogicalVolume(topSolid, sipmMat, "TopSiPMLV");
    {
        auto* va = new G4VisAttributes(G4Colour(1.0, 0.2, 0.0, 0.6));
        va->SetForceSolid(true);
        fTopSiPMLV->SetVisAttributes(va);
    }

    for (G4int i = 0; i < kNTopSiPMs; ++i) {
        const G4int globalId = 2 * kNEndSiPMs + i; // 16..35
        const G4double cx    = TopSiPMCenterX(i);

        auto* topPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(cx, kBarHalfY + kTopHalfY, 0.0),
            fTopSiPMLV, "TopSiPMPV", worldLV, false, globalId, true);

        fSiPMSurfaces[globalId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(globalId),
            fBarPhys, topPhys, sipmSurface);
    }

    return worldPhys;
}

// ── ConstructSDandField ───────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField() {
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* sipmSD = new SiPMSD("SiPMSD");
    sdManager->AddNewDetector(sipmSD);

    SetSensitiveDetector(fEndSiPMLV, sipmSD);
    SetSensitiveDetector(fTopSiPMLV, sipmSD);
}
