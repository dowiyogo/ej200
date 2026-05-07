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

// Mylar wrapping film (G4_MYLAR base).
// RINDEX = 1.65 across the optical range.  Used as a thin physical-volume
// wrapper around the scintillator bar to provide reflection via TIR and
// Fresnel equations without requiring a G4OpticalSurface at the air boundary.
G4Material* CreateMylar();

// Reflective wrapping surface (aluminum foil + black film, Betancourt 2020).
// dielectric_metal | polished | REFLECTIVITY = 0.85–0.92 (Al foil, 350–800 nm).
// Apply as G4LogicalBorderSurface(WrapPV → WorldPV) so photons hitting the
// outer face of any wrap segment are reflected back into the bar.
// At each wavelength, the photon is reflected with probability R(λ) or
// absorbed (lost) with probability 1-R(λ). No refraction into WorldPV.
G4OpticalSurface* CreateWrapReflectorSurface();

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

// SiPM optical surface.
// dielectric_dielectric | polished | DETECTIONEFFICIENCY = PDE(λ)
// PDE table from Hamamatsu S13360 series (33 points, 300–940 nm).
// Apply as G4LogicalBorderSurface(bar → sipmPhys) for every SiPM volume.
G4OpticalSurface* CreateSiPMSurface();

} // namespace Materials
