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
#include "G4SubtractionSolid.hh"
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

// Top SiPM (on +Y face): thin slab, 6×6 mm² active area
static constexpr G4double kTopHalfX   = 3.0  * mm;  // 6 mm in X
static constexpr G4double kTopHalfY   = 0.25 * mm;  // thickness
static constexpr G4double kTopHalfZ   = 3.0  * mm;  // 6 mm in Z

// ── Static helpers ────────────────────────────────────────────────────────────
G4double DetectorConstruction::TopSiPMCenterX(G4int idx) {
    if (idx < 0 || idx >= kNTopSiPMs) return 0.0;
    if (idx < 35) return (-692.0 + 20.0 * idx) * mm;
    return (12.0 + 20.0 * (idx - 35)) * mm;
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

    auto& readoutCmd = fMessenger->DeclareMethod(
        "readout",
        &DetectorConstruction::SetReadoutConfiguration,
        "Set readout/wrapping configuration: End (default), Top, or EndTop.");
    readoutCmd.SetParameterName("configuration", false);

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

void DetectorConstruction::SetReadoutConfiguration(G4String configuration) {
    std::transform(configuration.begin(), configuration.end(), configuration.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    configuration.erase(
        std::remove_if(configuration.begin(), configuration.end(),
                       [](unsigned char ch) { return ch == '-' || ch == '_' || ch == ' '; }),
        configuration.end());

    G4String canonical;
    if (configuration == "END") canonical = "End";
    else if (configuration == "TOP") canonical = "Top";
    else if (configuration == "ENDTOP") canonical = "EndTop";
    else {
        G4cerr << "[DetectorConstruction] Unknown /det/readout \"" << configuration
               << "\". Use End, Top, or EndTop.\n";
        return;
    }
    if (canonical == fReadoutConfiguration) return;

    fReadoutConfiguration = canonical;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

G4bool DetectorConstruction::IsEndInstrumented() const {
    return fReadoutConfiguration == "End" || fReadoutConfiguration == "EndTop";
}

G4bool DetectorConstruction::IsTopInstrumented() const {
    return fReadoutConfiguration == "Top" || fReadoutConfiguration == "EndTop";
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
    fTopSiPMLV = nullptr;
    fSiPMSurfaces.clear();
    fReflectorSurfaces.clear();

    // Reflective panels cover every non-instrumented face. Explicit sibling
    // volumes avoid global skin surfaces shadowing active SiPM faces.
    auto* reflector = Materials::CreateBarSkinReflector();
    const G4double foilHalfT = 0.5 * um;
    auto* reflYMinusSolid = new G4Box("ReflectorYMinusSolid",
                                      kBarHalfX, foilHalfT, kBarHalfZ);
    auto* reflZSolid = new G4Box("ReflectorZSolid",
                                 kBarHalfX, kBarHalfY, foilHalfT);
    auto* reflXSolid = new G4Box("ReflectorXSolid",
                                 foilHalfT, kBarHalfY, kBarHalfZ);
    auto* reflYMinusLV = new G4LogicalVolume(reflYMinusSolid, worldMat,
                                             "ReflectorYMinusLV");
    auto* reflZLV = new G4LogicalVolume(reflZSolid, worldMat, "ReflectorZLV");
    auto* reflXLV = new G4LogicalVolume(reflXSolid, worldMat, "ReflectorXLV");
    reflYMinusLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    reflZLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    reflXLV->SetVisAttributes(G4VisAttributes::GetInvisible());

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

    G4VSolid* reflYPlusSolid = new G4Box(
        "ReflectorYPlusBaseSolid", kBarHalfX, foilHalfT, kBarHalfZ);
    if (IsTopInstrumented()) {
        auto* window = new G4Box(
            "ReflectorYPlusWindowSolid", kTopHalfX, 2.0 * foilHalfT, kTopHalfZ);
        for (G4int i = 0; i < kNTopSiPMs; ++i) {
            reflYPlusSolid = new G4SubtractionSolid(
                "ReflectorYPlusCut_" + std::to_string(i),
                reflYPlusSolid, window, nullptr,
                G4ThreeVector(TopSiPMCenterX(i), 0.0, 0.0));
        }
    }
    auto* reflYPlusLV = new G4LogicalVolume(
        reflYPlusSolid, worldMat, "ReflectorYPlusLV");
    reflYPlusLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* reflYPlusPhys = new G4PVPlacement(
        nullptr, {0.0, +(kBarHalfY + foilHalfT), 0.0},
        reflYPlusLV, "ReflectorYPlusPV", worldLV, false, 5, true);
    fReflectorSurfaces["+Y"] = new G4LogicalBorderSurface(
        "BarReflector_YPlus", fBarPhys, reflYPlusPhys, reflector);
    if (!IsEndInstrumented()) {
        auto* reflXMinusPhys = new G4PVPlacement(
            nullptr, {-(kBarHalfX + foilHalfT), 0.0, 0.0},
            reflXLV, "ReflectorXMinusPV", worldLV, false, 1, true);
        auto* reflXPlusPhys = new G4PVPlacement(
            nullptr, {+(kBarHalfX + foilHalfT), 0.0, 0.0},
            reflXLV, "ReflectorXPlusPV", worldLV, false, 2, true);
        fReflectorSurfaces["-X"] = new G4LogicalBorderSurface(
            "BarReflector_XMinus", fBarPhys, reflXMinusPhys, reflector);
        fReflectorSurfaces["+X"] = new G4LogicalBorderSurface(
            "BarReflector_XPlus", fBarPhys, reflXPlusPhys, reflector);
    }

    if (IsEndInstrumented()) {
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

    if (IsTopInstrumented()) {
        auto* topSolid = new G4Box("TopSiPMSolid", kTopHalfX, kTopHalfY, kTopHalfZ);
        fTopSiPMLV = new G4LogicalVolume(topSolid, sipmMat, "TopSiPMLV");
        auto* va = new G4VisAttributes(G4Colour(1.0, 0.2, 0.0, 0.6));
        va->SetForceSolid(true);
        fTopSiPMLV->SetVisAttributes(va);

        for (G4int i = 0; i < kNTopSiPMs; ++i) {
            const G4int globalId = 2 * kNEndSiPMs + i;
            const G4double cx = TopSiPMCenterX(i);

            auto* topPhys = new G4PVPlacement(
                nullptr, G4ThreeVector(cx, kBarHalfY - kTopHalfY, 0.0),
                fTopSiPMLV, "TopSiPMPV", barLV, false, globalId, true);

            fSiPMSurfaces[globalId] = new G4LogicalBorderSurface(
                "SiPMSurf_" + std::to_string(globalId), fBarPhys, topPhys, sipmSurface);
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
    if (fTopSiPMLV != nullptr) SetSensitiveDetector(fTopSiPMLV, sipmSD);
}
