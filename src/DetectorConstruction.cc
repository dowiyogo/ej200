#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"

// ── Geometry constants ────────────────────────────────────────────────────────
static constexpr G4double kBarHalfX   = 700.0 * mm;  // 1.4 m total length
static constexpr G4double kBarHalfY   =  30.0 * mm;  // 60 mm width
static constexpr G4double kBarHalfZ   =   5.0 * mm;  // 10 mm height

// Mylar wrap thickness (25 µm film)
static constexpr G4double kMylarThick = 0.025 * mm;

// End SiPM (on ±X face): thin slab, 6×6 mm² active area
static constexpr G4double kEndHalfX   = 0.25 * mm;
static constexpr G4double kEndHalfY   = 3.0  * mm;  // 6 mm in Y
static constexpr G4double kEndHalfZ   = 3.0  * mm;  // 6 mm in Z
static constexpr G4double kEndPitch   = 2*kEndHalfY + 1.5*mm; // 7.5 mm

// Top SiPM (on +Y face): thin slab, 6×6 mm² active area
static constexpr G4double kTopHalfX   = 3.0  * mm;  // 6 mm in X
static constexpr G4double kTopHalfY   = 0.25 * mm;  // thickness
static constexpr G4double kTopHalfZ   = 3.0  * mm;  // 6 mm in Z

// ── Static helpers ────────────────────────────────────────────────────────────

// Maximum X available for SiPM centres (must stay inside bar footprint)
static constexpr G4double kTopSiPMMaxX = kBarHalfX - kTopHalfX;  // 697 mm

G4int DetectorConstruction::ComputeNTopSiPMs(G4double pitch) {
    // Number of SiPMs that can be placed symmetrically within ±kTopSiPMMaxX
    // at the requested pitch.  Always returns an even number ≥ 2.
    if (pitch <= 0.0) return 2;
    G4int n = static_cast<G4int>(std::floor(2.0 * kTopSiPMMaxX / pitch)) + 1;
    if (n < 2)        n = 2;
    if (n % 2 != 0)   n -= 1;   // enforce symmetry (even count)
    return n;
}

// X-centre of top SiPM index idx (0-based) given nTotal SiPMs at pitch
G4double DetectorConstruction::TopSiPMCenterX(G4int idx, G4double pitch,
                                               G4int nTotal) {
    // Symmetric distribution: centres at (idx - (nTotal-1)/2) * pitch
    return (idx - 0.5 * (nTotal - 1)) * pitch;
}

G4int DetectorConstruction::FaceType(G4int globalId) {
    if (globalId < kNEndSiPMs)       return 0;  // end left
    if (globalId < 2 * kNEndSiPMs)   return 1;  // end right
    return 2;                                    // top
}

G4int DetectorConstruction::LocalId(G4int globalId) {
    if (globalId < kNEndSiPMs)     return globalId;
    if (globalId < 2*kNEndSiPMs)   return globalId - kNEndSiPMs;
    return globalId - 2*kNEndSiPMs;
}

// ── Constructor / Destructor ──────────────────────────────────────────────────
DetectorConstruction::DetectorConstruction() {
    fNTopSiPMs = ComputeNTopSiPMs(fTopSiPMPitch);

    // UI messenger — available after construction, before /run/initialize
    fMessenger = new G4GenericMessenger(this, "/det/", "Detector geometry control");

    auto& cmd = fMessenger->DeclareMethodWithUnit(
        "topSiPMPitch", "mm",
        &DetectorConstruction::SetTopSiPMPitch,
        "Set distance between adjacent top SiPMs [mm]. Triggers geometry rebuild.\n"
        "  Valid range: ~10 mm … 73 mm (limited by bar length and SiPM count).\n"
        "  Common values: 40 (4 cm), 50 (5 cm), 70 (7 cm, default).");
    cmd.SetParameterName("pitch", false);
    cmd.SetRange("pitch > 0");
}

DetectorConstruction::~DetectorConstruction() {
    delete fMessenger;
}

// ── SetTopSiPMPitch ───────────────────────────────────────────────────────────
void DetectorConstruction::SetTopSiPMPitch(G4double pitchMm) {
    // pitchMm is already in G4 internal units (mm=1) because DeclareMethodWithUnit
    // performs the conversion before calling this function.
    fTopSiPMPitch = pitchMm;
    fNTopSiPMs    = ComputeNTopSiPMs(fTopSiPMPitch);
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// ── Construct ─────────────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();

    // ── Materials ────────────────────────────────────────────────────────────
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* mylarMat = Materials::CreateMylar();
    G4Material* barMat   = Materials::CreateEJ200();
    G4Material* sipmMat  = Materials::CreateSiPMCoupling();

    // Air RINDEX required for optical-photon tracking
    {
        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", {2.0*eV, 4.0*eV}, {1.0, 1.0});
        worldMat->SetMaterialPropertiesTable(mpt);
    }

    // SiPM coupling surface (bar–SiPM boundary, polished dielectric–dielectric)
    auto* sipmSurface = Materials::CreateSiPMSurface();

    // ── World volume ─────────────────────────────────────────────────────────
    auto* worldSolid = new G4Box("WorldSolid", 1.6*m, 0.25*m, 0.25*m);
    auto* worldLV    = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPhys  = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                         nullptr, false, 0, true);

    // ── Mylar wrap — thin envelope around the bar ────────────────────────────
    // The 25 µm Mylar layer (n=1.65) acts as a passive reflector:
    //   • Photons from bar (n=1.58) → Mylar: small Fresnel reflection, mostly transmitted
    //   • Photons from Mylar → air (n=1.0): TIR for angles > arcsin(1/1.65) ≈ 37.3°
    // This replaces the G4LogicalSkinSurface on the bar, avoiding the overhead
    // of per-step surface lookups on every bar–world boundary crossing.
    const G4double wHX = kBarHalfX + kMylarThick;
    const G4double wHY = kBarHalfY + kMylarThick;
    const G4double wHZ = kBarHalfZ + kMylarThick;

    auto* wrapSolid = new G4Box("WrapSolid", wHX, wHY, wHZ);
    auto* wrapLV    = new G4LogicalVolume(wrapSolid, mylarMat, "WrapLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.15));
        va->SetForceSolid(true);
        wrapLV->SetVisAttributes(va);
    }
    fWrapPhys = new G4PVPlacement(nullptr, {}, wrapLV, "WrapPV", worldLV,
                                  false, 0, true);

    // ── Scintillating bar — placed inside the Mylar wrap ─────────────────────
    auto* barSolid = new G4Box("BarSolid", kBarHalfX, kBarHalfY, kBarHalfZ);
    auto* barLV    = new G4LogicalVolume(barSolid, barMat, "BarLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.1, 0.3, 1.0, 0.25));
        va->SetForceSolid(true);
        barLV->SetVisAttributes(va);
    }
    // Bar is the daughter of WrapLV — centred at origin inside the Mylar shell
    fBarPhys = new G4PVPlacement(nullptr, {}, barLV, "BarPV", wrapLV,
                                 false, 0, true);

    // Note: no G4LogicalSkinSurface on barLV — the Mylar volume handles reflection.

    // ── End SiPMs — 8×1 array on each ±X face ────────────────────────────────
    // Placed in WorldLV, flush against the outer face of the Mylar wrap.
    // x = ±(wHX + kEndHalfX)  ← just outside the wrap
    auto* endSolid = new G4Box("EndSiPMSolid", kEndHalfX, kEndHalfY, kEndHalfZ);
    fEndSiPMLV     = new G4LogicalVolume(endSolid, sipmMat, "EndSiPMLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.0, 0.85, 0.2, 0.6));
        va->SetForceSolid(true);
        fEndSiPMLV->SetVisAttributes(va);
    }

    fSiPMSurfaces.clear();

    for (G4int i = 0; i < kNEndSiPMs; ++i) {
        const G4double cy = (i - 3.5) * kEndPitch;

        const G4int leftId = i;
        auto* leftPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(-(wHX + kEndHalfX), cy, 0.0),
            fEndSiPMLV, "EndSiPMLeft_PV", worldLV, false, leftId, true);
        // Border surface: WrapPV → SiPM (photon exits Mylar into SiPM coupling)
        fSiPMSurfaces[leftId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(leftId),
            fWrapPhys, leftPhys, sipmSurface);

        const G4int rightId = i + kNEndSiPMs;
        auto* rightPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(+(wHX + kEndHalfX), cy, 0.0),
            fEndSiPMLV, "EndSiPMRight_PV", worldLV, false, rightId, true);
        fSiPMSurfaces[rightId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(rightId),
            fWrapPhys, rightPhys, sipmSurface);
    }

    // ── Top SiPMs — N SiPMs along the +Y face, configurable pitch ────────────
    // N = fNTopSiPMs is computed from fTopSiPMPitch (default 70 mm → 20 SiPMs).
    // Placed in WorldLV, flush against the +Y outer face of the Mylar wrap.
    // y = wHY + kTopHalfY  ← just above the wrap top face
    auto* topSolid = new G4Box("TopSiPMSolid", kTopHalfX, kTopHalfY, kTopHalfZ);
    fTopSiPMLV     = new G4LogicalVolume(topSolid, sipmMat, "TopSiPMLV");
    {
        auto* va = new G4VisAttributes(G4Colour(1.0, 0.2, 0.0, 0.6));
        va->SetForceSolid(true);
        fTopSiPMLV->SetVisAttributes(va);
    }

    for (G4int i = 0; i < fNTopSiPMs; ++i) {
        const G4int    globalId = 2 * kNEndSiPMs + i;  // 16, 17, …
        const G4double cx       = TopSiPMCenterX(i, fTopSiPMPitch, fNTopSiPMs);

        auto* topPhys = new G4PVPlacement(
            nullptr,
            G4ThreeVector(cx, wHY + kTopHalfY, 0.0),
            fTopSiPMLV, "TopSiPMPV", worldLV, false, globalId, true);

        fSiPMSurfaces[globalId] = new G4LogicalBorderSurface(
            "SiPMSurf_" + std::to_string(globalId),
            fWrapPhys, topPhys, sipmSurface);
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
