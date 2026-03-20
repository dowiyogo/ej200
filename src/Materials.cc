#include "Materials.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

namespace Materials {

// ---------------------------------------------------------------------------
G4Material* CreateEJ200() {
    auto* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // Espectro de emision del EJ-200 — digitalizado del datasheet de Eljen
    // Technology (figura "EJ-200 Emission Spectrum").
    //
    // El espectro cubre 380–500 nm con pico en 425 nm.  Se usaron 22 puntos
    // para capturar el pico estrecho y las colas con fidelidad suficiente.
    //
    // ORDEN: longitud de onda DESCENDENTE (equivale a energia ASCENDENTE),
    // como exige G4MaterialPropertiesTable para propiedades vectoriales.
    //
    // Valores de amplitud relativa leidos del grafico del datasheet:
    //   λ(nm) | amplitud
    //   500   | 0.010
    //   490   | 0.020
    //   480   | 0.050
    //   470   | 0.100
    //   465   | 0.145
    //   460   | 0.200
    //   455   | 0.280
    //   450   | 0.370
    //   445   | 0.460
    //   440   | 0.560
    //   435   | 0.680
    //   432   | 0.780
    //   430   | 0.860
    //   428   | 0.930
    //   426   | 0.980
    //   425   | 1.000   <- pico
    //   423   | 0.970
    //   420   | 0.880
    //   415   | 0.680
    //   410   | 0.460
    //   405   | 0.280
    //   400   | 0.150
    //   395   | 0.070
    //   390   | 0.030
    //   385   | 0.010
    //   380   | 0.002
    const G4int nEm = 26;
    G4double wl_nm[nEm]  = {500, 490, 480, 470, 465, 460, 455, 450, 445, 440,
                             435, 432, 430, 428, 426, 425, 423, 420, 415, 410,
                             405, 400, 395, 390, 385, 380};
    G4double relOut[nEm] = {0.010, 0.020, 0.050, 0.100, 0.145, 0.200, 0.280, 0.370,
                             0.460, 0.560, 0.680, 0.780, 0.860, 0.930, 0.980, 1.000,
                             0.970, 0.880, 0.680, 0.460, 0.280, 0.150, 0.070, 0.030,
                             0.010, 0.002};

    const G4double hc = 1239.84193 * eV * nm; // h·c in eV·nm

    // Convertir a energias foton (wl_nm esta en orden DESCENDENTE, por lo
    // que photonE queda en orden ASCENDENTE, como requiere G4).
    G4double photonE[nEm], spectrum[nEm];
    for (G4int i = 0; i < nEm; ++i) {
        photonE[i] = hc / (wl_nm[i] * nm);
        spectrum[i] = relOut[i];
    }

    // Optical properties of EJ-200
    const G4int nOpt = 4;
    G4double eOpt[nOpt]  = {2.0*eV, 2.6*eV, 3.1*eV, 4.0*eV};
    G4double rIdx[nOpt]  = {1.58,   1.58,   1.58,   1.58};
    G4double absL[nOpt]  = {3.8*m,  3.8*m,  3.8*m,  3.8*m};
    //G4double rayl[nOpt]  = {1.5*m,  1.5*m,  1.5*m,  1.5*m};

    auto* mpt = new G4MaterialPropertiesTable();

    // Refractive index, bulk absorption, Rayleigh scattering
    mpt->AddProperty("RINDEX",    eOpt, rIdx, nOpt);
    mpt->AddProperty("ABSLENGTH", eOpt, absL, nOpt);
    //mpt->AddProperty("RAYLEIGH",  eOpt, rayl, nOpt);

    // Scintillation (G4 v11+ property names)
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonE, spectrum, nEm);
    mpt->AddConstProperty("SCINTILLATIONYIELD",        10000.0 / MeV);
    mpt->AddConstProperty("RESOLUTIONSCALE",           1.0);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1",       1.0);

    mat->SetMaterialPropertiesTable(mpt);
    return mat;
}

// ---------------------------------------------------------------------------
G4Material* CreateSiPMCoupling() {
    // Use SiO2 as the base material and override RINDEX to match the bar.
    // This models a perfect optical-coupling compound between bar and SiPM,
    // eliminating TIR at the bar–SiPM interface (arcsin(1.0/1.58) ≈ 39° → 90°).
    auto* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    const G4int n = 4;
    G4double e[n] = {2.0*eV, 2.6*eV, 3.1*eV, 4.0*eV};
    G4double r[n] = {1.58,   1.58,   1.58,   1.58};  // matched to EJ-200

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, r, n);
    mat->SetMaterialPropertiesTable(mpt);
    return mat;
}

// ---------------------------------------------------------------------------
G4OpticalSurface* CreateBarSurface() {
    // Interfaz barra-aire: dielectric_dielectric polished.
    // Geant4 aplica las ecuaciones de Fresnel automaticamente.
    // Para angulos > arcsin(1/1.58) = 39.3 deg ocurre TIR.
    // No se necesita MPT: el RINDEX de los materiales lo define todo.
    auto* surf = new G4OpticalSurface("BarSurface");
    surf->SetType(dielectric_dielectric);
    surf->SetModel(unified);
    surf->SetFinish(polished);
    surf->SetSigmaAlpha(0.0);
    return surf;
}

// ---------------------------------------------------------------------------
G4OpticalSurface* CreateSiPMSurface() {
    // SiPM detection surface.
    // dielectric_dielectric | polished: photon enters the SiPM volume via
    // Snell's law (no TIR since n_SiPM = n_bar = 1.58 via coupling material).
    // DETECTIONEFFICIENCY is read by SiPMSD::ProcessHits and applied manually.
    // PDE data from Hamamatsu S13360-6025 (or equivalent 6 mm SiPM, 25 μm cell).

    auto* surf = new G4OpticalSurface("SiPMSurface");
    surf->SetType(dielectric_dielectric);
    surf->SetModel(unified);
    surf->SetFinish(polished);
    surf->SetSigmaAlpha(0.0);

    const G4double hc = 1239.84193 * eV * nm;

    const G4int n = 33;
    G4double wl_nm[n] = {
        300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500,
        520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900,
        920, 940
    };
    G4double pde[n] = {
        0.000, 0.050, 0.120, 0.180, 0.260, 0.330, 0.380, 0.400,
        0.405, 0.403, 0.390, 0.370, 0.340, 0.310, 0.280, 0.240,
        0.210, 0.180, 0.150, 0.120, 0.100, 0.080, 0.060, 0.050,
        0.040, 0.030, 0.025, 0.020, 0.015, 0.010, 0.008, 0.004, 0.001
    };

    // Convert wavelengths to photon energies; G4 requires ascending energy order.
    G4double energy[n], detEff[n];
    for (G4int i = 0; i < n; ++i) {
        energy[n - 1 - i]  = hc / wl_nm[i];
        detEff[n - 1 - i]  = pde[i];
    }

    auto* mpt = new G4MaterialPropertiesTable();
    // No DETECTIONEFFICIENCY here — PDE is applied manually in SiPMSD.
    // The surface is purely optical-contact: polished dielectric_dielectric
    // with n_SiPM = n_bar = 1.58 ensures full transmission (no TIR, no loss).
    surf->SetMaterialPropertiesTable(mpt);
    return surf;
}

} // namespace Materials
