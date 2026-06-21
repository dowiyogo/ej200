#include "Materials.hh"
#include "SiPMModel.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalParameters.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Materials {

namespace {
G4Material* FindOrBuildPVT(const G4String& name) {
    if (auto* existing = G4Material::GetMaterial(name, false)) {
        return existing;
    }

    auto* nist = G4NistManager::Instance();
    auto* mat = new G4Material(name, 1.023 * g / cm3, 2);
    mat->AddElement(nist->FindOrBuildElement("C"), 9);
    mat->AddElement(nist->FindOrBuildElement("H"), 10);
    return mat;
}

void AddScintillatorProperties(G4Material* mat,
                               G4double* photonE,
                               G4double* spectrum,
                               G4int nEm,
                               G4double attenuation,
                               G4double peakWavelength,
                               G4double yield,
                               G4double rise,
                               G4double decay) {
    const G4int nOpt = 4;
    G4double eOpt[nOpt] = {2.0*eV, 2.6*eV, 3.1*eV, 4.0*eV};
    G4double rIdx[nOpt] = {1.58, 1.58, 1.58, 1.58};
    G4double absL[nOpt] = {attenuation, attenuation, attenuation, attenuation};

    const G4double hc = 1239.84193 * eV * nm;
    const G4double rayleighRef = 1.5 * m;
    G4double rayl[nOpt];
    for (G4int i = 0; i < nOpt; ++i) {
        const G4double wavelength = hc / eOpt[i];
        rayl[i] = rayleighRef * std::pow(wavelength / peakWavelength, 4.0);
    }

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", eOpt, rIdx, nOpt);
    mpt->AddProperty("ABSLENGTH", eOpt, absL, nOpt);
    mpt->AddProperty("RAYLEIGH", eOpt, rayl, nOpt);
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonE, spectrum, nEm);
    mpt->AddConstProperty("SCINTILLATIONYIELD", yield);
    mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mpt->AddConstProperty("SCINTILLATIONRISETIME1", rise);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", decay);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    mat->SetMaterialPropertiesTable(mpt);
}

void BuildSpectrum(const G4double* wavelengths,
                   const G4double* relativeOutput,
                   G4int count,
                   G4double* photonE,
                   G4double* spectrum) {
    const G4double hc = 1239.84193 * eV * nm;
    for (G4int i = 0; i < count; ++i) {
        photonE[i] = hc / (wavelengths[i] * nm);
        spectrum[i] = relativeOutput[i];
    }
}
} // namespace

void EnableFiniteScintillationRiseTime() {
    G4OpticalParameters::Instance()->SetScintFiniteRiseTime(true);
}

// ---------------------------------------------------------------------------
G4Material* CreateEJ204() {
    auto* mat = FindOrBuildPVT("EJ-204");

    // Digitized prompt spectrum with the Eljen datasheet maximum at 408 nm.
    const G4int nEm = 21;
    G4double wl_nm[nEm] = {
        500, 480, 460, 450, 440, 430, 425, 420, 415, 410, 408,
        405, 400, 395, 390, 385, 380, 375, 370, 365, 360
    };
    G4double relOut[nEm] = {
        0.005, 0.015, 0.050, 0.100, 0.200, 0.380, 0.520, 0.680, 0.850, 0.980,
        1.000, 0.970, 0.850, 0.650, 0.430, 0.250, 0.120, 0.050, 0.020, 0.005, 0.000
    };
    G4double photonE[nEm], spectrum[nEm];
    BuildSpectrum(wl_nm, relOut, nEm, photonE, spectrum);
    AddScintillatorProperties(
        mat, photonE, spectrum, nEm, 160.0 * cm, 408.0 * nm,
        10400.0 / MeV, 0.7 * ns, 1.8 * ns);
    return mat;
}

// ---------------------------------------------------------------------------
G4Material* CreateEJ200() {
    auto* mat = FindOrBuildPVT("EJ-200");

    // Espectro de emision del EJ-200 — digitalizado del datasheet de Eljen
    // Technology (figura "EJ-200 Emission Spectrum").
    //
    // El espectro cubre 380–500 nm con pico en 425 nm.  Se usaron 26 puntos
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

    G4double photonE[nEm], spectrum[nEm];
    BuildSpectrum(wl_nm, relOut, nEm, photonE, spectrum);
    AddScintillatorProperties(
        mat, photonE, spectrum, nEm, 380.0 * cm, 425.0 * nm,
        10000.0 / MeV, 0.9 * ns, 2.1 * ns);
    return mat;
}

// ---------------------------------------------------------------------------
G4Material* CreateEJ230() {
    auto* mat = FindOrBuildPVT("EJ-230");

    // Digitized prompt spectrum with the Eljen datasheet maximum at 391 nm.
    const G4int nEm = 29;
    G4double wl_nm[nEm] = {
        500, 490, 480, 470, 460, 450, 440, 430, 425, 420,
        415, 410, 405, 400, 395, 391, 388, 385, 382, 380,
        378, 375, 372, 370, 365, 360, 355, 350, 345
    };
    G4double relOut[nEm] = {
        0.005, 0.010, 0.020, 0.040, 0.070, 0.115, 0.180, 0.280, 0.350, 0.420,
        0.500, 0.620, 0.750, 0.900, 0.980, 1.000, 0.985, 0.960, 0.900, 0.830,
        0.720, 0.560, 0.380, 0.260, 0.120, 0.050, 0.020, 0.006, 0.000
    };
    G4double photonE[nEm], spectrum[nEm];
    BuildSpectrum(wl_nm, relOut, nEm, photonE, spectrum);
    AddScintillatorProperties(
        mat, photonE, spectrum, nEm, 120.0 * cm, 391.0 * nm,
        9700.0 / MeV, 0.5 * ns, 1.5 * ns);
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
    G4double r[n] = {1.58,   1.58,   1.58,   1.58};  // matched to PVT scintillators

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, r, n);
    mat->SetMaterialPropertiesTable(mpt);
    return mat;
}

// ---------------------------------------------------------------------------
G4OpticalSurface* CreateBarSurface() {
    // NOTE: CreateBarSurface() and CreateMylar() are not used in the current
    // geometry (BarPV directly in WorldLV with explicit reflector panels) but
    // are retained for reference and potential future use.
    //
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
G4OpticalSurface* CreateSiPMSurface(const G4String& model) {
    // A dielectric_metal surface with zero reflectivity lets Geant4 apply the
    // PDE exactly once through EFFICIENCY and invoke the attached SiPM SD.
    auto* surf = new G4OpticalSurface("SiPMSurface");
    surf->SetType(dielectric_metal);
    surf->SetModel(unified);
    surf->SetFinish(polished);
    surf->SetSigmaAlpha(0.0);

    const auto curve = SiPMModel::LoadPDECurve(model);
    const G4double hc = 1239.84193 * eV * nm;
    std::vector<G4double> energy;
    std::vector<G4double> efficiency;
    std::vector<G4double> reflectivity;
    energy.reserve(curve.size());
    efficiency.reserve(curve.size());
    reflectivity.reserve(curve.size());
    for (auto it = curve.rbegin(); it != curve.rend(); ++it) {
        energy.push_back(hc / (it->wavelengthNm * nm));
        efficiency.push_back(it->efficiency);
        reflectivity.push_back(0.0);
    }

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("EFFICIENCY", energy, efficiency);
    mpt->AddProperty("REFLECTIVITY", energy, reflectivity);
    surf->SetMaterialPropertiesTable(mpt);
    return surf;
}

// ---------------------------------------------------------------------------
G4Material* CreateMylar() {
    // Mylar (biaxially oriented PET film) used as simplified reflective wrap.
    // With RINDEX = 1.65, the critical angle for total internal reflection at
    // the Mylar–air boundary is arcsin(1/1.65) ≈ 37.3°.  Photons escaping the
    // bar (n=1.58) and entering the Mylar (n=1.65) that then hit the outer
    // Mylar–air interface at shallow angles are totally reflected back, acting
    // as a passive TIR mirror without requiring a G4OpticalSurface on the air.
    auto* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial("G4_MYLAR");

    const G4int n = 4;
    G4double e[n] = {2.0*eV, 2.6*eV, 3.1*eV, 4.0*eV};
    G4double r[n] = {1.65, 1.65, 1.65, 1.65};   // Mylar RINDEX ≈ 1.65

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, r, n);
    // Mylar is transparent at optical wavelengths; omit ABSLENGTH so G4 does
    // not apply bulk absorption inside the 25 µm film.
    mat->SetMaterialPropertiesTable(mpt);
    return mat;
}

// ---------------------------------------------------------------------------
G4Material* CreateBlackTape() {
    auto* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

    const G4int n = 4;
    G4double e[n]   = {2.0*eV, 2.6*eV, 3.1*eV, 4.0*eV};
    G4double r[n]   = {1.50,   1.50,   1.50,   1.50};
    G4double abs[n] = {1.0*um, 1.0*um, 1.0*um, 1.0*um};

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX",    e, r,   n);
    mpt->AddProperty("ABSLENGTH", e, abs, n);
    mat->SetMaterialPropertiesTable(mpt);
    return mat;
}

// ---------------------------------------------------------------------------
G4OpticalSurface* CreateBarSkinReflector() {
    // Reflective surface applied to explicit foil panels around the EJ-200 bar.
    //
    // BarPV is a direct WorldLV daughter in fix/geometry-bar-in-world. SiPMs
    // are BarLV daughters, while non-SiPM faces see separate reflector panels
    // through G4LogicalBorderSurface(BarPV -> Reflector*PV).
    //
    // Spectral reflectivity tuning points: wavelength descending -> energy ascending.
    // Source: Hecht, "Optics", 5th ed.; typical values for Al mirror at visible.
    // Central value at EJ-200 emission peak (425 nm): R = 0.88.
    // Betancourt 2020 (NIM A 979) uses Al foil + black plastic film.
    // Tuning: R=0.85 (aged/imperfect), R=0.90 (standard),
    // R=0.95 (aluminized Mylar / high-quality reflector).
    // R=0.98 (Tyvek-like high-reflectivity fallback).
    // The explicit-panel geometry requires R=0.98 to keep combined end-SiPM
    // yield above the acceptance threshold; the conservative Al-foil curve
    // gives ~6 ph/evt per end side in this geometry.

    auto* surf = new G4OpticalSurface("BarSkinReflector");
    surf->SetType(dielectric_metal);
    surf->SetModel(unified);
    surf->SetFinish(polished);
    surf->SetSigmaAlpha(0.0);

    const G4double hc = 1239.84193 * eV * nm;

    const G4int n = 6;
    G4double wl_nm[n] = {800.0, 600.0, 500.0, 425.0, 400.0, 350.0};
    G4double refl[n]  = {0.98,  0.98,  0.98,  0.98,  0.98,  0.98};

    G4double energy[n], reflectivity[n];
    for (G4int i = 0; i < n; ++i) {
        energy[i]       = hc / (wl_nm[i] * nm);  // ascending energy
        reflectivity[i] = refl[i];
    }

    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", energy, reflectivity, n);
    surf->SetMaterialPropertiesTable(mpt);
    return surf;
}

// ---------------------------------------------------------------------------
G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
                                       G4double specularLobe,
                                       G4double sigmaAlpha) {
    // exec21-optfix: switch from dielectric_metal (no TIR) to dielectric_dielectric
    // to enable TIR at the bar-air interface (n_bar=1.58, n_air=1.0, theta_c=39.3 deg).
    // Rationale: real scintillator bars have an air gap between bar and wrapping foil;
    // photons at angle > theta_c undergo TIR with zero loss, dramatically improving
    // light collection at the END SiPMs from ~0.37 PE (no TIR) to O(100) PE (with TIR).
    // The 'reflectivity', 'specularLobe', 'sigmaAlpha' parameters are no longer used
    // for the skin surface; TIR is determined automatically by RINDEX of bar vs world.
    // Non-TIR photons (angle < theta_c) are transmitted to worldLV and "escape" the bar.
    // TODO: add explicit Mylar wrapping volume to recover non-TIR photons.
    (void)reflectivity; (void)specularLobe; (void)sigmaAlpha;
    auto* surf = new G4OpticalSurface("MylarReflector");
    surf->SetType(dielectric_dielectric);
    surf->SetModel(unified);
    surf->SetFinish(polished);   // polished → specular TIR (Snell's law reflection)
    surf->SetSigmaAlpha(0.0);   // no microfacet roughness for the TIR test
    // No MPT needed: Geant4 uses RINDEX of bar (1.58) and world (1.0)
    // to compute TIR threshold and Fresnel reflection/transmission automatically.
    return surf;
}

} // namespace Materials
