#include "SiPMSD.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4GenericMessenger.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cctype>

// ---------------------------------------------------------------------------
// PDE curves.
//
// Broadcom AFBR-S4N66P024M:
//   Digitized from AFBR-S4N66P024M-DS105 Figure 6.  The datasheet table gives
//   spectral range 250-900 nm and peak PDE 63% at 420 nm (12 V overvoltage).
//
// Hamamatsu S14160:
//   Digitized from s14160-1310ps_etc_kapd1070e Figure "Photon detection
//   efficiency vs. wavelength".  Datasheet table gives peak PDE 18% for
//   10 um pixels and 32% for 15 um pixels at 460 nm.
//
// Hamamatsu S13360:
//   Legacy curve kept for backward compatibility with previous analyses.
//
// SPTR defaults below are simulation defaults.  The supplied datasheets provide
// PDE/electrical curves but do not quote a numeric SPTR sigma for every model.
// ---------------------------------------------------------------------------
namespace {
    constexpr G4double kAFBRS4N66P024M_WL_nm[] = {
        250, 260, 270, 280, 290, 300, 320, 340, 360, 380, 400,
        420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620,
        640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840,
        860, 880, 900
    };
    constexpr G4double kAFBRS4N66P024M_PDE[] = {
        0.055, 0.065, 0.100, 0.200, 0.320, 0.420, 0.430, 0.500,
        0.520, 0.600, 0.630, 0.630, 0.610, 0.580, 0.550, 0.520,
        0.480, 0.440, 0.410, 0.380, 0.340, 0.310, 0.290, 0.270,
        0.250, 0.230, 0.205, 0.180, 0.155, 0.135, 0.115, 0.095,
        0.080, 0.065, 0.050
    };

    constexpr G4double kS14160_WL_nm[] = {
        290, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480,
        500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900
    };
    constexpr G4double kS14160_10um_PDE[] = {
        0.000, 0.050, 0.080, 0.110, 0.120, 0.130, 0.150, 0.165,
        0.180, 0.185, 0.183, 0.180, 0.175, 0.165, 0.150, 0.135,
        0.125, 0.115, 0.105, 0.095, 0.088, 0.080, 0.070, 0.060,
        0.050, 0.043, 0.036, 0.030, 0.025, 0.020, 0.015, 0.010
    };
    constexpr G4double kS14160_15um_PDE[] = {
        0.000, 0.080, 0.130, 0.170, 0.210, 0.220, 0.260, 0.290,
        0.320, 0.320, 0.320, 0.310, 0.290, 0.270, 0.250, 0.235,
        0.220, 0.200, 0.180, 0.160, 0.145, 0.130, 0.110, 0.095,
        0.080, 0.068, 0.058, 0.050, 0.043, 0.037, 0.032, 0.028
    };

    constexpr G4double kS13360_WL_nm[] = {
        300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500,
        520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940
    };
    constexpr G4double kS13360_PDE[] = {
        0.000, 0.050, 0.120, 0.180, 0.260, 0.330, 0.380, 0.400,
        0.405, 0.403, 0.390, 0.370, 0.340, 0.310, 0.280, 0.240,
        0.210, 0.180, 0.150, 0.120, 0.100, 0.080, 0.060, 0.050,
        0.040, 0.030, 0.025, 0.020, 0.015, 0.010, 0.008, 0.004, 0.001
    };

    template <size_t N>
    constexpr G4int SizeOf(const G4double (&)[N]) {
        return static_cast<G4int>(N);
    }
} // namespace

// ---------------------------------------------------------------------------
SiPMSD::SiPMSD(const G4String& name)
    : G4VSensitiveDetector(name)
{
    ApplySensorDefaults(fSensorModel);

    // ── UI messenger ─────────────────────────────────────────────────────────
    fMessenger = new G4GenericMessenger(this, "/sipm/", "SiPM detector control");

    auto& sensorCmd = fMessenger->DeclareMethod(
        "sensorModel",
        &SiPMSD::SetSensorModelByName,
        "Set SiPM model. Supported values:\n"
        "  Broadcom, AFBR-S4N66P024M\n"
        "  S14160-6010PS, S14160-3010PS, S14160-1310PS\n"
        "  S14160-6015PS, S14160-3015PS, S14160-1315PS\n"
        "  Hamamatsu, S13360, S13360-6025 (legacy)");
    sensorCmd.SetParameterName("model", false);

    auto& sptrCmd = fMessenger->DeclareMethodWithUnit(
        "sptrSigma", "ns",
        &SiPMSD::SetSPTRSigma,
        "Override the SiPM SPTR sigma [ns].\n"
        "  Defaults: Broadcom = 0.030 ns, Hamamatsu models = 0.150 ns.");
    sptrCmd.SetParameterName("sigma", false);
    sptrCmd.SetRange("sigma >= 0");

    auto& elecCmd = fMessenger->DeclareMethodWithUnit(
        "electronicsSigma", "ns",
        &SiPMSD::SetElectronicsSigma,
        "Set FastIC+ electronics jitter RMS [ns].\n"
        "  Default: 0.030 ns (= 30 ps).");
    elecCmd.SetParameterName("sigma", false);
    elecCmd.SetRange("sigma >= 0");
}

SiPMSD::~SiPMSD() {
    delete fMessenger;
}

// ---------------------------------------------------------------------------
G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    auto* track = step->GetTrack();

    // Only optical photons entering from a geometry boundary
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

    auto* pre = step->GetPreStepPoint();
    if (pre->GetStepStatus() != fGeomBoundary) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    // ── Identify SiPM by copy number ────────────────────────────────────────
    // G4Track::GetVolume() always returns the current volume of the track,
    // avoiding the pre-step touchable ambiguity at geometry boundaries.
    auto* pv = track->GetVolume();
    if (pv == nullptr) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }
    const G4int globalId = pv->GetCopyNo();

    // ── Photon kinematics ────────────────────────────────────────────────────
    const G4double energy    = pre->GetKineticEnergy();
    const G4double energy_eV = energy / eV;
    const G4double wl_nm     = (energy_eV > 0.0) ? (1239.84193 / energy_eV) : 0.0;

    const G4double rawTimeNs = track->GetGlobalTime() / ns;

    const G4ThreeVector pos = pre->GetPosition();

    // ── PDE Bernoulli trial ──────────────────────────────────────────────────
    const G4double pde      = GetPDE(energy);
    const G4bool   detected = (G4UniformRand() < pde);

    if (detected) {
        const G4double sptrJitterNs = G4RandGauss::shoot(0.0, fSPTRSigma) / ns;
        const G4double timeNs       = rawTimeNs + sptrJitterNs;
        const G4int eventId =
            G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

        auto* ea = dynamic_cast<EventAction*>(
            G4EventManager::GetEventManager()->GetUserEventAction());
        if (ea != nullptr) {
            const G4int face = DetectorConstruction::FaceType(globalId);
            if      (face == 0) ea->AddEndLeftHit();
            else if (face == 1) ea->AddEndRightHit();
            else                ea->AddTopHit();
        }

        auto* am = G4AnalysisManager::Instance();
        am->FillNtupleIColumn(0, 0, eventId);
        am->FillNtupleIColumn(0, 1, DetectorConstruction::FaceType(globalId));
        am->FillNtupleIColumn(0, 2, globalId);
        am->FillNtupleIColumn(0, 3, DetectorConstruction::LocalId(globalId));
        am->FillNtupleDColumn(0, 4, timeNs);
        am->FillNtupleDColumn(0, 5, rawTimeNs);
        am->FillNtupleDColumn(0, 6, sptrJitterNs);
        am->FillNtupleDColumn(0, 7, energy_eV);
        am->FillNtupleDColumn(0, 8, wl_nm);
        am->FillNtupleDColumn(0, 9, pde);
        am->FillNtupleDColumn(0, 10, pos.x() / mm);
        am->FillNtupleDColumn(0, 11, pos.y() / mm);
        am->FillNtupleDColumn(0, 12, pos.z() / mm);
        const G4double gunX = ea ? ea->GetGunXmm() : 0.0;
        am->FillNtupleDColumn(0, 13, gunX);
        am->FillNtupleIColumn(0, 14, GetSensorModelId());
        am->FillNtupleDColumn(0, 15, fSPTRSigma / ns);
        am->FillNtupleDColumn(0, 16, fElectronicsSigma / ns);
        am->AddNtupleRow(0);
    }

    track->SetTrackStatus(fStopAndKill);
    return detected;
}

// ---------------------------------------------------------------------------
G4double SiPMSD::GetPDE(G4double energy) const {
    if (energy <= 0.0 || fPDEWavelengthNm.empty()) return 0.0;

    const G4double wlNm = 1239.84193 / (energy / eV);
    if (wlNm < fPDEWavelengthNm.front() || wlNm > fPDEWavelengthNm.back()) {
        return 0.0;
    }

    auto upper = std::lower_bound(
        fPDEWavelengthNm.begin(), fPDEWavelengthNm.end(), wlNm);
    if (upper == fPDEWavelengthNm.begin()) return std::clamp(fPDEValue.front(), 0.0, 1.0);
    if (upper == fPDEWavelengthNm.end())   return std::clamp(fPDEValue.back(),  0.0, 1.0);

    const size_t hi = static_cast<size_t>(upper - fPDEWavelengthNm.begin());
    const size_t lo = hi - 1;
    const G4double x0 = fPDEWavelengthNm[lo];
    const G4double x1 = fPDEWavelengthNm[hi];
    const G4double y0 = fPDEValue[lo];
    const G4double y1 = fPDEValue[hi];
    const G4double t = (wlNm - x0) / (x1 - x0);
    return std::clamp(y0 + t * (y1 - y0), 0.0, 1.0);
}

void SiPMSD::SetSensorModelByName(const G4String& modelName) {
    G4String lowered = modelName;
    std::transform(
        lowered.begin(), lowered.end(), lowered.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    lowered.erase(
        std::remove_if(
            lowered.begin(), lowered.end(),
            [](unsigned char ch) {
                return ch == ' ' || ch == '_' || ch == '-';
            }),
        lowered.end());

    if (lowered == "broadcom" || lowered == "afbr" ||
        lowered == "afbrs4n66p024m" || lowered == "nuvmt" ||
        lowered == "nuvmt14m" || lowered == "broadcomnuvmt14m") {
        ApplySensorDefaults(SensorModel::kBroadcomAFBRS4N66P024M);
        return;
    }

    if (lowered == "hamamatsu" || lowered == "s13360" || lowered == "s133606025") {
        ApplySensorDefaults(SensorModel::kHamamatsuS13360);
        return;
    }

    if (lowered == "s141601310ps" || lowered == "s141603010ps" ||
        lowered == "s141606010ps" || lowered == "s1416010um" ||
        lowered == "hamamatsus1416010um") {
        ApplySensorDefaults(SensorModel::kHamamatsuS14160_10um);
        return;
    }

    if (lowered == "s141601315ps" || lowered == "s141603015ps" ||
        lowered == "s141606015ps" || lowered == "s1416015um" ||
        lowered == "hamamatsus1416015um" || lowered == "s14160") {
        ApplySensorDefaults(SensorModel::kHamamatsuS14160_15um);
        return;
    }

    G4cerr << "[SiPMSD] Unknown /sipm/sensorModel '" << modelName
           << "'. Keeping " << GetSensorModelName() << G4endl;
}

void SiPMSD::LoadPDECurve(SensorModel model) {
    const G4double* wl  = nullptr;
    const G4double* pde = nullptr;
    G4int n = 0;

    switch (model) {
        case SensorModel::kBroadcomAFBRS4N66P024M:
            wl = kAFBRS4N66P024M_WL_nm;
            pde = kAFBRS4N66P024M_PDE;
            n = SizeOf(kAFBRS4N66P024M_WL_nm);
            break;
        case SensorModel::kHamamatsuS14160_10um:
            wl = kS14160_WL_nm;
            pde = kS14160_10um_PDE;
            n = SizeOf(kS14160_WL_nm);
            break;
        case SensorModel::kHamamatsuS14160_15um:
            wl = kS14160_WL_nm;
            pde = kS14160_15um_PDE;
            n = SizeOf(kS14160_WL_nm);
            break;
        case SensorModel::kHamamatsuS13360:
            wl = kS13360_WL_nm;
            pde = kS13360_PDE;
            n = SizeOf(kS13360_WL_nm);
            break;
    }

    fPDEWavelengthNm.assign(wl, wl + n);
    fPDEValue.assign(pde, pde + n);
}

void SiPMSD::ApplySensorDefaults(SensorModel model) {
    fSensorModel = model;
    fSPTRSigma   = (model == SensorModel::kBroadcomAFBRS4N66P024M) ? 30.0e-3 : 150.0e-3;
    LoadPDECurve(model);
    G4cout << "[SiPMSD] Sensor model: " << GetSensorModelName()
           << ", SPTR sigma = " << (fSPTRSigma / ns) << " ns" << G4endl;
}

G4int SiPMSD::GetSensorModelId() const {
    return static_cast<G4int>(fSensorModel);
}

G4String SiPMSD::GetSensorModelName() const {
    if (fSensorModel == SensorModel::kHamamatsuS13360) {
        return "Hamamatsu S13360-6025 (legacy)";
    }
    if (fSensorModel == SensorModel::kHamamatsuS14160_10um) {
        return "Hamamatsu S14160-6010PS/3010PS/1310PS";
    }
    if (fSensorModel == SensorModel::kHamamatsuS14160_15um) {
        return "Hamamatsu S14160-6015PS/3015PS/1315PS";
    }
    return "Broadcom AFBR-S4N66P024M";
}