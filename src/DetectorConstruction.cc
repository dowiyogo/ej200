#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "CLHEP/Units/PhysicalConstants.h"

// ── Cylinder geometry constants ───────────────────────────────────────────────
static constexpr G4double kCylRadius    = 12.5 * mm;  // 25 mm diameter
static constexpr G4double kCylHalfH    = 12.5 * mm;  // 25 mm height
static constexpr G4double kAirGapThick  = 0.10 * mm;  // air gap between scint and Vikuiti
static constexpr G4double kWrapThick    = 0.10 * mm;  // Vikuiti 3M ESR wrap

// SiPMs: AFBR-S4N66P024M, 6×6 mm² active area, 2×2 grid per cap
// Grid centres at (±kSiPMOff, ±kSiPMOff) — touching at x=0 and y=0, fits in r=12.5
static constexpr G4double kSiPMHalf    = 3.0  * mm;  // half-side of 6 mm active area
static constexpr G4double kSiPMThick   = 0.25 * mm;  // half-thickness of SiPM slab
static constexpr G4double kSiPMOff     = 3.0  * mm;  // grid offset: centres at (±3, ±3) mm

// ── FaceType / LocalId ───────────────────────────────────────────────────────
// (static inline in header — no definition needed here)

// ── Construct ────────────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();

    // ── World (air) ──────────────────────────────────────────────────────────
    auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    {
        // GROUPVEL fix: after TIR bounces on the mantle, G4Transportation would
        // pick up the air velocity (c). Setting GROUPVEL = c/n_bar in the world
        // material prevents superluminal steps (same fix applied in feat/endtop-sslg4).
        // WorldLV photons are killed immediately by TIR so there is no physical
        // side-effect; the fix only corrects the velocity used during boundary crossing.
        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX",  {2.0*eV, 4.0*eV}, {1.0, 1.0});
        const G4double vg = CLHEP::c_light / 1.58;
        mpt->AddProperty("GROUPVEL", {2.0*eV, 4.0*eV}, {vg, vg});
        worldMat->SetMaterialPropertiesTable(mpt);
    }

    auto* worldBox = new G4Box("WorldBox", 100.*mm, 100.*mm, 100.*mm);
    auto* worldLV  = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPV  = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                       nullptr, false, 0, true);

    // ── Scintillator cylinder (EJ-228) ───────────────────────────────────────
    fScintMat      = Materials::CreateEJ228();
    auto* cylSolid = new G4Tubs("CylSolid", 0., kCylRadius, kCylHalfH,
                                 0., CLHEP::twopi);
    auto* cylLV    = new G4LogicalVolume(cylSolid, fScintMat, "CylLV");
    auto* visScint = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.35));
    visScint->SetForceSolid(false);
    cylLV->SetVisAttributes(visScint);
    auto* cylPV    = new G4PVPlacement(nullptr, {}, cylLV, "CylPV",
                                       worldLV, false, 0, true);

    // ── Air gap (annular tube, mantle only) ──────────────────────────────────
    const G4double rAirIn  = kCylRadius;
    const G4double rAirOut = kCylRadius + kAirGapThick;
    auto* airSolid = new G4Tubs("AirGapSolid", rAirIn, rAirOut, kCylHalfH,
                                 0., CLHEP::twopi);
    auto* airLV    = new G4LogicalVolume(airSolid, worldMat, "AirGapLV");
    airLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* airPV    = new G4PVPlacement(nullptr, {}, airLV, "AirGapPV",
                                       worldLV, false, 0, true);

    // ── Vikuiti 3M ESR wrap (annular shell, mantle only) ────────────────────
    const G4double rWrapIn  = rAirOut;
    const G4double rWrapOut = rWrapIn + kWrapThick;
    auto* wrapMat  = Materials::CreateMylar();   // structural host; reflectance via surface
    auto* wrapSolid = new G4Tubs("WrapSolid", rWrapIn, rWrapOut, kCylHalfH,
                                  0., CLHEP::twopi);
    auto* wrapLV   = new G4LogicalVolume(wrapSolid, wrapMat, "VikuitiLV");
    auto* visWrap  = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85, 0.6));
    visWrap->SetForceSolid(true);
    wrapLV->SetVisAttributes(visWrap);
    auto* wrapPV   = new G4PVPlacement(nullptr, {}, wrapLV, "VikuitiPV",
                                       worldLV, false, 0, true);

    // ── Optical surface: scintillator ↔ air gap (polished → TIR) ────────────
    // Geant4 applies Fresnel equations automatically with RINDEX(EJ-228)=1.58
    // and RINDEX(air)=1.0, giving θ_c = arcsin(1/1.58) = 39.3°.
    new G4LogicalBorderSurface("CylAirSurf", cylPV, airPV,
                               Materials::CreateBarSurface());

    // ── Optical surface: air gap ↔ Vikuiti (specular metal R=0.98) ──────────
    new G4LogicalBorderSurface("AirVikuitiSurf", airPV, wrapPV,
                               Materials::CreateVikuitiSurface());

    // ── SiPMs (shared logical volume, distinct physical placements) ──────────
    auto* sipmMat   = Materials::CreateSiPMCoupling();  // n=1.58, eliminates TIR at cap face
    auto* sipmSolid = new G4Box("SiPMSolid", kSiPMHalf, kSiPMHalf, kSiPMThick);
    fSiPMLV         = new G4LogicalVolume(sipmSolid, sipmMat, "SiPMLV");
    auto* visSiPM   = new G4VisAttributes(G4Colour(0.1, 0.1, 0.9, 0.9));
    visSiPM->SetForceSolid(true);
    fSiPMLV->SetVisAttributes(visSiPM);

    auto* sipmSurf = Materials::CreateSiPMSurface(fSiPMModel);

    // 2×2 grid: centres at (±kSiPMOff, ±kSiPMOff)
    // Adjacent SiPMs touch (not overlap) at x=0 and y=0 — Geant4 allows this.
    const G4double ox = kSiPMOff, oy = kSiPMOff;
    const G4double gridXY[4][2] = {{-ox,-oy},{-ox,+oy},{+ox,-oy},{+ox,+oy}};

    G4cout << "\n[DetectorConstruction] EJ-228 cylinder — SiPM layout:\n";

    // Top cap: face_type 0, global IDs 0–3, at z = +(kCylHalfH + kSiPMThick)
    const G4double zTop = kCylHalfH + kSiPMThick;
    for (G4int i = 0; i < kNTopSiPMs; ++i) {
        G4ThreeVector pos(gridXY[i][0], gridXY[i][1], zTop);
        auto* pv = new G4PVPlacement(nullptr, pos, fSiPMLV,
                                     "TopSiPM_PV_" + std::to_string(i),
                                     worldLV, false, i, true);
        fSiPMSurfaces[i] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(i), cylPV, pv, sipmSurf);
        G4cout << "  Top  ID " << i << " at ("
               << pos.x()/mm << ", " << pos.y()/mm << ", " << pos.z()/mm << ") mm\n";
    }

    // Bottom cap: face_type 1, global IDs 4–7, at z = -(kCylHalfH + kSiPMThick)
    const G4double zBot = -(kCylHalfH + kSiPMThick);
    for (G4int i = 0; i < kNBottomSiPMs; ++i) {
        const G4int gid = kNTopSiPMs + i;
        G4ThreeVector pos(gridXY[i][0], gridXY[i][1], zBot);
        auto* pv = new G4PVPlacement(nullptr, pos, fSiPMLV,
                                     "BotSiPM_PV_" + std::to_string(gid),
                                     worldLV, false, gid, true);
        fSiPMSurfaces[gid] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(gid), cylPV, pv, sipmSurf);
        G4cout << "  Bot  ID " << gid << " at ("
               << pos.x()/mm << ", " << pos.y()/mm << ", " << pos.z()/mm << ") mm\n";
    }

    G4cout << "  Vikuiti R=0.98, air gap " << kAirGapThick/mm << " mm, wrap "
           << kWrapThick/mm << " mm\n"
           << "  TIR angle θ_c = arcsin(1/1.58) = 39.3° on mantle\n\n";

    return worldPV;
}

// ── ConstructSDandField ───────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField() {
    auto* sipmSD = new SiPMSD("SiPMSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
    SetSensitiveDetector(fSiPMLV, sipmSD);
}
