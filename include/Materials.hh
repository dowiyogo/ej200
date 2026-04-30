#pragma once
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

// --------------------------------------------------------------------------
// Materials — centralised optical material and surface factory (from Repo2)
// --------------------------------------------------------------------------
namespace Materials {

// EJ-200 plastic scintillator (G4_PLASTIC_SC_VINYLTOLUENE base).
// Full 20-point emission spectrum (peak 425 nm), n=1.58,
// yield=10 000 ph/MeV, tau1=2.1 ns, bulk att. length=3.8 m.
G4Material* CreateEJ200();

// EJ-230 fast-timing plastic scintillator (polyvinyltoluene base).
// Emission spectrum from Eljen EJ-228/EJ-230 datasheet (peak 391 nm),
// n=1.58, yield=9700 ph/MeV, rise=0.5 ns, decay=1.5 ns,
// bulk attenuation length=120 cm.
G4Material* CreateEJ230();

// Mylar wrapping film (G4_MYLAR base).
// RINDEX = 1.65 across the optical range.  Used as a thin physical-volume
// wrapper around the scintillator bar to provide reflection via TIR and
// Fresnel equations without requiring a G4OpticalSurface at the air boundary.
G4Material* CreateMylar();

// Optical coupling material for SiPM volumes (G4_SILICON_DIOXIDE base).
// RINDEX set to 1.58 to match the bar and eliminate TIR at the bar–SiPM
// interface — equivalent to perfect optical grease coupling.
G4Material* CreateSiPMCoupling();

// Reflective wrapping surface (Tyvek-like).
// dielectric_metal | groundfrontpainted | R = 0.95
// Apply as G4LogicalBorderSurface(bar → world) to model all non-SiPM faces.
G4OpticalSurface* CreateBarSurface();

// SiPM optical surface.
// dielectric_dielectric | polished; PDE is applied manually in SiPMSD from
// the selected sensor datasheet curve.
// Apply as G4LogicalBorderSurface(bar → sipmPhys) for every SiPM volume.
G4OpticalSurface* CreateSiPMSurface();

} // namespace Materials
