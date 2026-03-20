#include "SiPMSD.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

// ---------------------------------------------------------------------------
// Tabla PDE — Hamamatsu S13360-6025 (celda 25 µm, 6×6 mm²)
// 33 puntos, 300–940 nm.  Fuente: datasheet Hamamatsu S13360 series.
// ---------------------------------------------------------------------------
namespace {
    constexpr G4int    kNPDE   = 33;
    constexpr G4double kWL_nm[kNPDE] = {
        300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500,
        520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940
    };
    constexpr G4double kPDE[kNPDE] = {
        0.000, 0.050, 0.120, 0.180, 0.260, 0.330, 0.380, 0.400,
        0.405, 0.403, 0.390, 0.370, 0.340, 0.310, 0.280, 0.240,
        0.210, 0.180, 0.150, 0.120, 0.100, 0.080, 0.060, 0.050,
        0.040, 0.030, 0.025, 0.020, 0.015, 0.010, 0.008, 0.004, 0.001
    };
} // namespace

// ---------------------------------------------------------------------------
SiPMSD::SiPMSD(const G4String& name)
    : G4VSensitiveDetector(name)
{
    // Construir la curva PDE usando G4MaterialPropertyVector::InsertValues().
    // InsertValues() mantiene orden ascendente de energía automáticamente.
    // hc = 1239.84193 eV·nm
    const G4double hc = 1239.84193 * eV * nm;
    for (G4int i = 0; i < kNPDE; ++i) {
        fPDEVec.InsertValues(hc / (kWL_nm[i] * nm), kPDE[i]);
    }
}

// ---------------------------------------------------------------------------
G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    auto* track = step->GetTrack();

    // Solo fotones opticos entrando desde una frontera geometrica
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

    auto* pre = step->GetPreStepPoint();
    if (pre->GetStepStatus() != fGeomBoundary) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }

    // ── Identificar el SiPM por copy number ─────────────────────────────────
    //
    // FIX: se usa G4Track::GetVolume() en lugar de
    //      pre->GetTouchableHandle()->GetCopyNumber().
    //
    // En pasos que empiezan en una frontera (fGeomBoundary), el touchable
    // del punto pre puede devolver el volumen ANTERIOR (la barra, copyNo=0)
    // en ciertas implementaciones de Geant4 11.x, haciendo que todos los
    // hits aparezcan en SiPM 0 (extremo izquierdo) y el otro extremo
    // quede aparentemente vacio.
    //
    // G4Track::GetVolume() siempre devuelve el volumen actual del track
    // (el SiPM recien entrado) de forma inequivoca en cualquier version.
    auto* pv = track->GetVolume();
    if (pv == nullptr) {
        track->SetTrackStatus(fStopAndKill);
        return false;
    }
    const G4int globalId = pv->GetCopyNo();

    // Cinematica del foton
    const G4double energy     = pre->GetKineticEnergy();
    const G4double energy_eV  = energy / eV;
    const G4double wl_nm      = (energy_eV > 0.0) ? (1239.84193 / energy_eV) : 0.0;
    const G4double time_ns    = track->GetGlobalTime() / ns;
    const G4ThreeVector pos   = pre->GetPosition();

    // Trial de Bernoulli con PDE
    // G4MaterialPropertyVector::Value() interpola linealmente con clamping
    const G4double pde      = GetPDE(energy);
    const G4bool   detected = (G4UniformRand() < pde);

    if (detected) {
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
        am->FillNtupleDColumn(0, 4, time_ns);
        am->FillNtupleDColumn(0, 5, energy_eV);
        am->FillNtupleDColumn(0, 6, wl_nm);
        am->FillNtupleDColumn(0, 7, pde);
        am->FillNtupleDColumn(0, 8, pos.x() / mm);
        am->FillNtupleDColumn(0, 9, pos.y() / mm);
        am->FillNtupleDColumn(0, 10, pos.z() / mm);
        // Posicion x del muon primario: extraida del EventAction en
        // BeginOfEventAction desde el G4PrimaryVertex.  Es constante para
        // todos los fotones del mismo evento y permite segmentar el analisis
        // por posicion longitudinal sin recurrir a metadatos externos.
        const G4double gunX = ea ? ea->GetGunXmm() : 0.0;
        am->FillNtupleDColumn(0, 11, gunX);
        am->AddNtupleRow(0);
    }

    // Matar el foton en cualquier caso (plano absorbente)
    track->SetTrackStatus(fStopAndKill);
    return detected;
}

// ---------------------------------------------------------------------------
G4double SiPMSD::GetPDE(G4double energy) const
{
    // G4MaterialPropertyVector::Value() hace interpolacion lineal con clamping:
    //   energy < E_min  -> PDE[0]   (= 0.000 en 940 nm, la mas baja)
    //   energy > E_max  -> PDE[N-1] (= 0.000 en 300 nm, la mas alta)
    // Nota: la tabla esta en orden ascendente de energia (longitudes de onda
    // descendentes), como exige G4PhysicsVector.
    return fPDEVec.Value(energy);
}
