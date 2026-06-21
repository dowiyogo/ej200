#pragma once
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

// --------------------------------------------------------------------------
// Materials — centralised optical material and surface factory (from Repo2)
// --------------------------------------------------------------------------
namespace Materials {

// Enable Geant4's biexponential scintillation timing sampler so each
// material's SCINTILLATIONRISETIME1 property is honored.
void EnableFiniteScintillationRiseTime();

// EJ-204 fast-timing plastic scintillator (PVT base).
// Emission peak 408 nm, n=1.58, yield=10 400 ph/MeV,
// rise=0.7 ns, decay=1.8 ns, bulk attenuation length=160 cm.
G4Material* CreateEJ204();

// EJ-200 plastic scintillator (PVT base).
// Full emission spectrum (peak 425 nm), n=1.58,
// yield=10 000 ph/MeV, rise=0.9 ns, decay=2.1 ns, bulk att. length=3.8 m.
G4Material* CreateEJ200();

// EJ-230 fast-timing plastic scintillator (PVT base).
// Emission peak 391 nm, n=1.58, yield=9700 ph/MeV,
// rise=0.5 ns, decay=1.5 ns, bulk attenuation length=120 cm.
G4Material* CreateEJ230();

// Mylar wrapping film (G4_MYLAR base).
// RINDEX = 1.65 across the optical range.  Used as a thin physical-volume
// wrapper around the scintillator bar to provide reflection via TIR and
// Fresnel equations without requiring a G4OpticalSurface at the air boundary.
G4Material* CreateMylar();

// Reflective wrapping surface for explicit reflector panels.
// Uses dielectric_metal with idealized high reflectivity R(λ) = 0.99.
G4OpticalSurface* CreateBarSkinReflector();

// Tunable aluminized-Mylar reflector for the end-only geometry.
// The unified-model remainder after specularLobe is Lambertian.
G4OpticalSurface* CreateMylarReflector(G4double reflectivity,
                                       G4double specularLobe,
                                       G4double sigmaAlpha,
                                       const G4String& finish);

// Black absorbing tape — high-density polymer with extremely short absorption
// length and low refractive index, used to simulate "black tape" / "vinyl"
// edge wrapping that absorbs all incident photons.
G4Material* CreateBlackTape();

// Optical coupling material for SiPM volumes (G4_SILICON_DIOXIDE base).
// RINDEX set to 1.58 to match the bar and eliminate TIR at the bar–SiPM
// interface — equivalent to perfect optical grease coupling.
G4Material* CreateSiPMCoupling();

// Reflective wrapping surface (Tyvek-like).
// dielectric_metal | groundfrontpainted | R = 0.95
// Apply as G4LogicalBorderSurface(bar → world) to model all non-SiPM faces.
G4OpticalSurface* CreateBarSurface();

// SiPM optical surface. The selected model's PDE file is installed as
// EFFICIENCY(lambda), with REFLECTIVITY(lambda)=0.
G4OpticalSurface* CreateSiPMSurface(const G4String& model,
                                    const G4String& efficiencyMode = "nominal");

} // namespace Materials
