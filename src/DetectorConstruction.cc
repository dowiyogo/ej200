#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "OrganicScintillatorFactory.hh"
#include "SiPMModel.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

// ── Geometry constants ────────────────────────────────────────────────────────
static constexpr G4double kBarHalfX   = 700.0 * mm;  // 1.4 m total length
static constexpr G4double kBarHalfY   =  30.0 * mm;  // 60 mm width
static constexpr G4double kBarHalfZ   =   5.0 * mm;  // 10 mm height

// End SiPM (on ±X face): thin slab, 6×6 mm² active area
static constexpr G4double kEndHalfX   = 0.25 * mm;
static constexpr G4double kEndHalfY   = 3.0  * mm;  // 6 mm in Y
static constexpr G4double kEndHalfZ   = 3.0  * mm;  // 6 mm in Z
static constexpr G4double kEndPitch   = 2*kEndHalfY + 1.5*mm; // 7.5 mm

// ── Static helpers ────────────────────────────────────────────────────────────
G4int DetectorConstruction::FaceType(G4int globalId) {
    if (globalId < kNEndSiPMs)       return 0;  // end left
    return 1;                                    // end right
}

G4int DetectorConstruction::LocalId(G4int globalId) {
    if (globalId < kNEndSiPMs)     return globalId;
    return globalId - kNEndSiPMs;
}

// ── Constructor / Destructor ──────────────────────────────────────────────────
DetectorConstruction::DetectorConstruction() {
    // UI messenger — available after construction, before /run/initialize
    fMessenger = new G4GenericMessenger(this, "/det/", "Detector geometry control");

    auto& matCmd = fMessenger->DeclareMethod(
        "scintillator",
        &DetectorConstruction::SetScintillatorCode,
        "Set SSLG4 organic scintillator: OPSC-106 (EJ-230, default), OPSC-101 (EJ-204), or OPSC-100 (EJ-200).");
    matCmd.SetParameterName("code", false);

    auto& legacyMatCmd = fMessenger->DeclareMethod(
        "scintillatorMaterial",
        &DetectorConstruction::SetScintillatorCode,
        "Legacy alias for /det/scintillator.");
    legacyMatCmd.SetParameterName("code", false);

    fMessenger->DeclareMethod(
        "edgeWrap",
        &DetectorConstruction::SetEdgeWrapMode,
        "Legacy no-op. The Mylar wrap volume was removed; reflection is handled\n"
        "by explicit reflector panels around BarLV.");

    fSiPMMessenger = new G4GenericMessenger(this, "/sipm/", "SiPM model control");
    auto& sipmModelCmd = fSiPMMessenger->DeclareMethod(
        "model", &DetectorConstruction::SetSiPMModel,
        "Set SiPM model. Default and currently supported: AFBR-S4N66P024M.");
    sipmModelCmd.SetParameterName("model", false);
}

DetectorConstruction::~DetectorConstruction() {
    delete fSiPMMessenger;
    delete fMessenger;
}

void DetectorConstruction::SetScintillatorCode(G4String code) {
    std::transform(code.begin(), code.end(), code.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    code.erase(std::remove_if(code.begin(), code.end(),
                              [](unsigned char ch) { return ch == '_' || ch == ' '; }),
               code.end());
    if (code == "EJ204") code = "OPSC-101";
    if (code == "EJ200") code = "OPSC-100";
    if (code == "EJ230") code = "OPSC-106";

    if (code != "OPSC-101" && code != "OPSC-100" && code != "OPSC-106") {
        G4cerr << "[DetectorConstruction] Unknown /det/scintillator \""
               << code << "\". Use OPSC-106 (EJ-230), OPSC-101 (EJ-204), or OPSC-100 (EJ-200).\n";
        return;
    }
    if (code == fScintCode) return;

    fScintCode = code;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

void DetectorConstruction::SetSiPMModel(G4String model) {
    const G4String canonical = SiPMModel::CanonicalName(model);
    if (canonical.empty()) {
        G4cerr << "[DetectorConstruction] Unknown /sipm/model \"" << model
               << "\". Use AFBR-S4N66P024M.\n";
        return;
    }
    if (canonical == fSiPMModel) return;
    fSiPMModel = canonical;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

// ── SetEdgeWrapMode ──────────────────────────────────────────────────────────
void DetectorConstruction::SetEdgeWrapMode(G4String mode) {
    G4cout << "[DetectorConstruction] /det/edgeWrap " << mode
           << " ignored: WrapLV was removed in fix/geometry-bar-in-world.\n";
}

// ── Construct ─────────────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();

    // ── Materials ────────────────────────────────────────────────────────────
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    fActiveScintillator =
        OrganicScintillatorFactory::GetInstance()->Get(fScintCode, true);
    if (fScintCode == "OPSC-101") {
        auto* mpt = fActiveScintillator->GetMaterialPropertiesTable();
        if (mpt != nullptr) {
            auto* attenuation = mpt->GetProperty("ABSLENGTH");
            G4bool overrideAttenuation =
                attenuation == nullptr || attenuation->GetVectorLength() == 0;
            if (!overrideAttenuation) {
                for (std::size_t i = 0; i < attenuation->GetVectorLength(); ++i) {
                    if (std::abs((*attenuation)[i] - 160.0 * cm) > 1.0e-9 * cm) {
                        overrideAttenuation = true;
                    }
                }
            }
            if (overrideAttenuation) {
                mpt->RemoveProperty("ABSLENGTH");
                mpt->AddProperty("ABSLENGTH",
                                 {1.5 * eV, 6.5 * eV},
                                 {160.0 * cm, 160.0 * cm});
                G4cout << "[DetectorConstruction] OPSC-101 ABSLENGTH overridden to 160 cm."
                       << G4endl;
            }
            if (!mpt->ConstPropertyExists("SCINTILLATIONRISETIME1") ||
                std::abs(mpt->GetConstProperty("SCINTILLATIONRISETIME1") -
                         0.7 * ns) > 1.0e-12 * ns) {
                mpt->RemoveConstProperty("SCINTILLATIONRISETIME1");
                mpt->AddConstProperty("SCINTILLATIONRISETIME1", 0.7 * ns);
                G4cout << "[DetectorConstruction] OPSC-101 rise time overridden to 0.7 ns."
                       << G4endl;
            }
        }
    }
    // Enforce EJ-230 datasheet values (Rev. Aug 2023): att=120 cm, rise=0.5 ns.
    // SSLG4 opsc-106 data already matches; this block guards against upstream drift.
    if (fScintCode == "OPSC-106") {
        auto* mpt = fActiveScintillator->GetMaterialPropertiesTable();
        if (mpt != nullptr) {
            auto* attenuation = mpt->GetProperty("ABSLENGTH");
            G4bool overrideAttenuation =
                attenuation == nullptr || attenuation->GetVectorLength() == 0;
            if (!overrideAttenuation) {
                for (std::size_t i = 0; i < attenuation->GetVectorLength(); ++i) {
                    if (std::abs((*attenuation)[i] - 120.0 * cm) > 1.0e-9 * cm) {
                        overrideAttenuation = true;
                    }
                }
            }
            if (overrideAttenuation) {
                mpt->RemoveProperty("ABSLENGTH");
                mpt->AddProperty("ABSLENGTH",
                                 {1.5 * eV, 6.5 * eV},
                                 {120.0 * cm, 120.0 * cm});
                G4cout << "[DetectorConstruction] OPSC-106 ABSLENGTH overridden to 120 cm."
                       << G4endl;
            }
            if (!mpt->ConstPropertyExists("SCINTILLATIONRISETIME1") ||
                std::abs(mpt->GetConstProperty("SCINTILLATIONRISETIME1") -
                         0.5 * ns) > 1.0e-12 * ns) {
                mpt->RemoveConstProperty("SCINTILLATIONRISETIME1");
                mpt->AddConstProperty("SCINTILLATIONRISETIME1", 0.5 * ns);
                G4cout << "[DetectorConstruction] OPSC-106 rise time overridden to 0.5 ns."
                       << G4endl;
            }
        }
    }
    G4Material* barMat   = fActiveScintillator;
    G4Material* sipmMat  = Materials::CreateSiPMCoupling();

    // Air RINDEX required for optical-photon tracking
    {
        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", {2.0*eV, 4.0*eV}, {1.0, 1.0});
        worldMat->SetMaterialPropertiesTable(mpt);
    }

    // SiPM detection surface: selected PDE as EFFICIENCY and zero reflectivity.
    auto* sipmSurface = Materials::CreateSiPMSurface(fSiPMModel);

    // ── World volume ─────────────────────────────────────────────────────────
    auto* worldSolid = new G4Box("WorldSolid", 1.6*m, 0.25*m, 0.25*m);
    auto* worldLV    = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPhys  = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                         nullptr, false, 0, true);

    // ── Scintillating bar — direct WorldLV daughter ──────────────────────────
    auto* barSolid = new G4Box("BarSolid", kBarHalfX, kBarHalfY, kBarHalfZ);
    auto* barLV    = new G4LogicalVolume(barSolid, barMat, "BarLV");
    {
        auto* va = new G4VisAttributes(G4Colour(0.1, 0.3, 1.0, 0.25));
        va->SetForceSolid(true);
        barLV->SetVisAttributes(va);
    }
    fBarPhys = new G4PVPlacement(nullptr, {}, barLV, "BarPV",
                                 worldLV, false, 0, true);

    fEndSiPMLV = nullptr;
    fSiPMSurfaces.clear();
    fReflectorSurfaces.clear();

    // Reflective panels cover every non-instrumented face. Explicit sibling
    // volumes avoid global skin surfaces shadowing active SiPM faces.
    auto* reflector = Materials::CreateMylarReflector();
    const G4double foilHalfT = 0.5 * um;
    auto* reflYMinusSolid = new G4Box("ReflectorYMinusSolid",
                                      kBarHalfX, foilHalfT, kBarHalfZ);
    auto* reflZSolid = new G4Box("ReflectorZSolid",
                                 kBarHalfX, kBarHalfY, foilHalfT);
    auto* reflYMinusLV = new G4LogicalVolume(reflYMinusSolid, worldMat,
                                             "ReflectorYMinusLV");
    auto* reflZLV = new G4LogicalVolume(reflZSolid, worldMat, "ReflectorZLV");
    reflYMinusLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    reflZLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto* reflYMinusPhys = new G4PVPlacement(
        nullptr, {0.0, -(kBarHalfY + foilHalfT), 0.0},
        reflYMinusLV, "ReflectorYMinusPV", worldLV, false, 0, true);
    auto* reflZMinusPhys = new G4PVPlacement(
        nullptr, {0.0, 0.0, -(kBarHalfZ + foilHalfT)},
        reflZLV, "ReflectorZMinusPV", worldLV, false, 3, true);
    auto* reflZPlusPhys = new G4PVPlacement(
        nullptr, {0.0, 0.0, +(kBarHalfZ + foilHalfT)},
        reflZLV, "ReflectorZPlusPV", worldLV, false, 4, true);

    fReflectorSurfaces["-Y"] = new G4LogicalBorderSurface(
        "BarReflector_YMinus", fBarPhys, reflYMinusPhys, reflector);
    fReflectorSurfaces["-Z"] = new G4LogicalBorderSurface(
        "BarReflector_ZMinus", fBarPhys, reflZMinusPhys, reflector);
    fReflectorSurfaces["+Z"] = new G4LogicalBorderSurface(
        "BarReflector_ZPlus", fBarPhys, reflZPlusPhys, reflector);

    auto* reflYPlusSolid = new G4Box(
        "ReflectorYPlusBaseSolid", kBarHalfX, foilHalfT, kBarHalfZ);
    auto* reflYPlusLV = new G4LogicalVolume(
        reflYPlusSolid, worldMat, "ReflectorYPlusLV");
    reflYPlusLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* reflYPlusPhys = new G4PVPlacement(
        nullptr, {0.0, +(kBarHalfY + foilHalfT), 0.0},
        reflYPlusLV, "ReflectorYPlusPV", worldLV, false, 5, true);
    fReflectorSurfaces["+Y"] = new G4LogicalBorderSurface(
        "BarReflector_YPlus", fBarPhys, reflYPlusPhys, reflector);
    {
        auto* endSolid = new G4Box("EndSiPMSolid", kEndHalfX, kEndHalfY, kEndHalfZ);
        fEndSiPMLV = new G4LogicalVolume(endSolid, sipmMat, "EndSiPMLV");
        auto* va = new G4VisAttributes(G4Colour(0.0, 0.85, 0.2, 0.6));
        va->SetForceSolid(true);
        fEndSiPMLV->SetVisAttributes(va);

        for (G4int i = 0; i < kNEndSiPMs; ++i) {
            const G4double cy = (i - 3.5) * kEndPitch;
            const G4int leftId = i;
            auto* leftPhys = new G4PVPlacement(
                nullptr, G4ThreeVector(-(kBarHalfX - kEndHalfX), cy, 0.0),
                fEndSiPMLV, "EndSiPMLeft_PV", barLV, false, leftId, true);
            fSiPMSurfaces[leftId] = new G4LogicalBorderSurface(
                "SiPMSurf_" + std::to_string(leftId), fBarPhys, leftPhys, sipmSurface);

            const G4int rightId = i + kNEndSiPMs;
            auto* rightPhys = new G4PVPlacement(
                nullptr, G4ThreeVector(+(kBarHalfX - kEndHalfX), cy, 0.0),
                fEndSiPMLV, "EndSiPMRight_PV", barLV, false, rightId, true);
            fSiPMSurfaces[rightId] = new G4LogicalBorderSurface(
                "SiPMSurf_" + std::to_string(rightId), fBarPhys, rightPhys, sipmSurface);
        }
    }

    return worldPhys;
}

// ── ConstructSDandField ───────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField() {
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* sipmSD = dynamic_cast<SiPMSD*>(
        sdManager->FindSensitiveDetector("SiPMSD", false));
    if (sipmSD == nullptr) {
        sipmSD = new SiPMSD("SiPMSD");
        sdManager->AddNewDetector(sipmSD);
    }
    sipmSD->SetModel(fSiPMModel);

    if (fEndSiPMLV != nullptr) SetSensitiveDetector(fEndSiPMLV, sipmSD);
}
