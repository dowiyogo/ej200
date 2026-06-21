#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "OrganicScintillatorFactory.hh"
#include "SiPMModel.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <string>

// ── Geometry constants ────────────────────────────────────────────────────────
static constexpr G4double kBarHalfX   = 700.0 * mm;  // 1.4 m total length
static constexpr G4double kBarHalfY   =  30.0 * mm;  // 60 mm width
static constexpr G4double kBarHalfZ   =   5.0 * mm;  // 10 mm height
static constexpr G4double kAirGapThickness = 0.10 * mm;
static constexpr G4double kReflectorThickness = 0.05 * mm;

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
        "Set SSLG4 organic scintillator: OPSC-101 (EJ-204, default) or OPSC-100 (EJ-200).");
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
        "Legacy no-op. Reflection is handled by explicit air-gap and Mylar panels.");

    fSiPMMessenger = new G4GenericMessenger(this, "/sipm/", "SiPM model control");
    auto& sipmModelCmd = fSiPMMessenger->DeclareMethod(
        "model", &DetectorConstruction::SetSiPMModel,
        "Set SiPM model. Default and currently supported: AFBR-S4N66P024M.");
    sipmModelCmd.SetParameterName("model", false);

    fGeometryMessenger =
        new G4GenericMessenger(this, "/ship/geom/", "End-only geometry control");
    auto& topSurfaceCmd = fGeometryMessenger->DeclareMethod(
        "topSurface", &DetectorConstruction::SetTopSurface,
        "Set top-face treatment: mylar (default, no Top SiPMs) or sipm (base EndTop behavior).");
    topSurfaceCmd.SetParameterName("mode", false);

    fMylarMessenger =
        new G4GenericMessenger(this, "/ship/geom/mylar/", "Tunable Mylar reflector");
    auto& reflectivityCmd = fMylarMessenger->DeclareMethod(
        "reflectivity", &DetectorConstruction::SetMylarReflectivity,
        "Set uniform aluminized-Mylar reflectivity.");
    reflectivityCmd.SetParameterName("value", false);
    reflectivityCmd.SetRange("value >= 0.0 && value <= 1.0");

    auto& lobeCmd = fMylarMessenger->DeclareMethod(
        "specularLobe", &DetectorConstruction::SetMylarSpecularLobe,
        "Set unified-model specular-lobe fraction; the remainder is Lambertian.");
    lobeCmd.SetParameterName("value", false);
    lobeCmd.SetRange("value >= 0.0 && value <= 1.0");

    auto& sigmaCmd = fMylarMessenger->DeclareMethodWithUnit(
        "sigmaAlpha", "deg", &DetectorConstruction::SetMylarSigmaAlpha,
        "Set unified-model microfacet roughness sigma alpha.");
    sigmaCmd.SetParameterName("angle", false);
    sigmaCmd.SetRange("angle >= 0.0");
}

DetectorConstruction::~DetectorConstruction() {
    delete fMylarMessenger;
    delete fGeometryMessenger;
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

    if (code != "OPSC-101" && code != "OPSC-100") {
        G4cerr << "[DetectorConstruction] Unknown /det/scintillator \""
               << code << "\". Use OPSC-101 or OPSC-100.\n";
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
    return fTopSurface == "sipm" &&
           (fReadoutConfiguration == "Top" || fReadoutConfiguration == "EndTop");
}

void DetectorConstruction::SetTopSurface(G4String mode) {
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (mode != "mylar" && mode != "sipm") {
        G4cerr << "[DetectorConstruction] Unknown /ship/geom/topSurface \""
               << mode << "\". Use mylar or sipm.\n";
        return;
    }
    if (mode == fTopSurface) return;
    fTopSurface = mode;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

void DetectorConstruction::SetMylarReflectivity(G4double value) {
    if (value == fMylarReflectivity) return;
    fMylarReflectivity = value;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

void DetectorConstruction::SetMylarSpecularLobe(G4double value) {
    if (value == fMylarSpecularLobe) return;
    fMylarSpecularLobe = value;
    if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
}

void DetectorConstruction::SetMylarSigmaAlpha(G4double value) {
    if (value == fMylarSigmaAlpha) return;
    fMylarSigmaAlpha = value;
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
    if (fActiveScintillator->GetMaterialPropertiesTable() == nullptr) {
        fActiveScintillator->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
    }
    {
        auto* mpt = fActiveScintillator->GetMaterialPropertiesTable();
        auto* rindex = mpt ? mpt->GetProperty("RINDEX") : nullptr;
        if (rindex == nullptr || rindex->GetVectorLength() == 0) {
            mpt->AddProperty("RINDEX",
                             {1.5 * eV, 6.5 * eV},
                             {1.58, 1.58});
            G4cout << "[DetectorConstruction] " << fScintCode
                   << " RINDEX set to 1.58 for optical tracking." << G4endl;
        }
    }
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
    G4Material* barMat   = fActiveScintillator;
    G4Material* airGapMat = worldMat;
    G4Material* reflectorMat = Materials::CreateMylar();
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
    fScintillatorAirSurfaces.clear();
    fReflectorSurfaces.clear();
    fBarSkinSurface = nullptr;

    if (fTopSurface == "mylar") {
        auto* scintAirSurface = Materials::CreateBarSurface();
        auto* airReflectorSurface = Materials::CreateMylarReflector(
            fMylarReflectivity, fMylarSpecularLobe, fMylarSigmaAlpha);

        struct PanelSpec {
            G4String key;
            G4String stem;
            G4ThreeVector airHalf;
            G4ThreeVector airCenter;
            G4ThreeVector reflectorHalf;
            G4ThreeVector reflectorCenter;
        };

        const G4double g = kAirGapThickness;
        const G4double m = kReflectorThickness;
        const PanelSpec panels[] = {
            {"+Y", "YPlus",
             {kBarHalfX, 0.5*g, kBarHalfZ + g},
             {0.0, kBarHalfY + 0.5*g, 0.0},
             {kBarHalfX, 0.5*m, kBarHalfZ + g + m},
             {0.0, kBarHalfY + g + 0.5*m, 0.0}},
            {"-Y", "YMinus",
             {kBarHalfX, 0.5*g, kBarHalfZ + g},
             {0.0, -(kBarHalfY + 0.5*g), 0.0},
             {kBarHalfX, 0.5*m, kBarHalfZ + g + m},
             {0.0, -(kBarHalfY + g + 0.5*m), 0.0}},
            {"+Z", "ZPlus",
             {kBarHalfX, kBarHalfY, 0.5*g},
             {0.0, 0.0, kBarHalfZ + 0.5*g},
             {kBarHalfX, kBarHalfY + g, 0.5*m},
             {0.0, 0.0, kBarHalfZ + g + 0.5*m}},
            {"-Z", "ZMinus",
             {kBarHalfX, kBarHalfY, 0.5*g},
             {0.0, 0.0, -(kBarHalfZ + 0.5*g)},
             {kBarHalfX, kBarHalfY + g, 0.5*m},
             {0.0, 0.0, -(kBarHalfZ + g + 0.5*m)}},
        };

        for (std::size_t i = 0; i < std::size(panels); ++i) {
            const auto& panel = panels[i];
            auto* airSolid = new G4Box(
                G4String("AirGap") + panel.stem + "Solid",
                panel.airHalf.x(), panel.airHalf.y(), panel.airHalf.z());
            auto* airLV = new G4LogicalVolume(
                airSolid, airGapMat, G4String("AirGap") + panel.stem + "LV");
            auto* airVA = new G4VisAttributes(G4Colour(0.6, 0.8, 1.0, 0.12));
            airVA->SetForceSolid(true);
            airLV->SetVisAttributes(airVA);
            auto* airPV = new G4PVPlacement(
                nullptr, panel.airCenter, airLV, G4String("AirGap") + panel.stem + "PV",
                worldLV, false, static_cast<G4int>(i), true);

            auto* reflSolid = new G4Box(
                G4String("Mylar") + panel.stem + "Solid",
                panel.reflectorHalf.x(), panel.reflectorHalf.y(), panel.reflectorHalf.z());
            auto* reflLV = new G4LogicalVolume(
                reflSolid, reflectorMat, G4String("Mylar") + panel.stem + "LV");
            auto* reflVA = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85, 0.35));
            reflVA->SetForceSolid(true);
            reflLV->SetVisAttributes(reflVA);
            auto* reflPV = new G4PVPlacement(
                nullptr, panel.reflectorCenter, reflLV, G4String("Mylar") + panel.stem + "PV",
                worldLV, false, static_cast<G4int>(i), true);

            fScintillatorAirSurfaces[panel.key] = new G4LogicalBorderSurface(
                G4String("ScintillatorToAirGap_") + panel.stem,
                fBarPhys, airPV, scintAirSurface);
            new G4LogicalBorderSurface(
                G4String("AirGapToScintillator_") + panel.stem,
                airPV, fBarPhys, scintAirSurface);

            fReflectorSurfaces[panel.key] = new G4LogicalBorderSurface(
                G4String("AirGapToReflector_") + panel.stem,
                airPV, reflPV, airReflectorSurface);
            new G4LogicalBorderSurface(
                G4String("ReflectorToAirGap_") + panel.stem,
                reflPV, airPV, airReflectorSurface);
        }
    } else {
        auto* reflector = Materials::CreateBarSkinReflector();
        fBarSkinSurface = new G4LogicalSkinSurface("BarSkin", barLV, reflector);
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
