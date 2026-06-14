#!/usr/bin/env python3
"""
validate_group_velocity.py

Validates that Geant4's effective photon group velocity in EJ-204
equals c/n = 299.7925/1.58 ≈ 189.742 mm/ns, as expected from a flat
RINDEX with no explicit GROUPVEL property set.

HOW TO USE
----------
Two modes:

  Mode A — theoretical check (no simulation needed, always runnable):
    python validate_group_velocity.py --theoretical

  Mode B — simulation output check (requires a ROOT file from a
    dedicated single-photon axial propagation run):
    python validate_group_velocity.py --root path/to/singlephoton_run.root

EXPECTED SIMULATION SETUP FOR MODE B
--------------------------------------
Run a Geant4 macro that:
  1. Generates a single 408 nm photon at x = -700 mm (left end face centre)
     with momentum direction +X (along bar axis).
  2. Records the hit time t_hit (ns) at the right SiPM (x = +700 mm).
  3. Repeats for N events; stores t_hit in a TTree branch "t_hit_ns".
  4. Expected: t_hit ≈ 1400 mm / 189.742 mm/ns ≈ 7.38 ns
     (if no scattering or absorption — only valid for photons that
      travel straight without reflection).

ROOT tree expected: tree name "photon_hits", branch "t_hit_ns".
"""

from __future__ import annotations

import argparse
import math
import sys


# ---------------------------------------------------------------------------
# Physical constants and simulation parameters
# ---------------------------------------------------------------------------
C_MM_NS = 299.792458      # speed of light, mm/ns
N_EJ204 = 1.58            # RINDEX from Materials.cc line 41
LAMBDA_PEAK_NM = 408.0    # peak emission, Materials.cc line 99
BAR_HALF_LENGTH_MM = 700.0 # half-length of EJ-204 bar in SSLG4 geometry
RAYLEIGH_REF_M = 1.5      # Materials.cc line 45
ABSLENGTH_CM = 160.0      # Materials.cc line 99


def theoretical_check() -> int:
    """Mode A: derive and print expected group velocity, no simulation needed."""
    v_g = C_MM_NS / N_EJ204
    t_traverse = (2 * BAR_HALF_LENGTH_MM) / v_g  # full bar length

    # Rayleigh at peak emission
    rayleigh_peak_mm = RAYLEIGH_REF_M * 1000.0  # 1.5 m → 1500 mm = 150 cm

    # Mean free path dominated by min(ABSLENGTH, Rayleigh) at peak
    abslength_mm = ABSLENGTH_CM * 10.0
    mean_free_path_mm = 1.0 / (1.0 / abslength_mm + 1.0 / rayleigh_peak_mm)

    print("=" * 60)
    print("EJ-204 photon group velocity — theoretical check")
    print("=" * 60)
    print(f"\nSimulation parameters (from src/Materials.cc):")
    print(f"  RINDEX           = {N_EJ204} (flat, no dispersion)")
    print(f"  ABSLENGTH        = {ABSLENGTH_CM} cm  (bulk, applied to optical path s)")
    print(f"  RAYLEIGH at peak = {rayleigh_peak_mm:.0f} mm = {rayleigh_peak_mm/10:.0f} cm")
    print(f"  No GROUPVEL property set")
    print()
    print(f"Derived quantities:")
    print(f"  v_g = c / n = {C_MM_NS:.6f} / {N_EJ204} = {v_g:.4f} mm/ns")
    print(f"  t_traverse (full bar, straight) = {2*BAR_HALF_LENGTH_MM:.0f} mm / {v_g:.4f} mm/ns")
    print(f"            = {t_traverse:.4f} ns")
    print(f"  Effective mean free path at peak = {mean_free_path_mm:.1f} mm")
    print(f"    (combined ABSLENGTH + Rayleigh; actual path is zig-zag)")
    print()
    print("Expected Geant4 behaviour:")
    print(f"  A 408 nm photon emitted at x = −700 mm along +X arrives at")
    print(f"  x = +700 mm at t ≈ {t_traverse:.4f} ns IF it travels straight")
    print(f"  without reflection (i.e., hitting neither +Y/−Y nor +Z/−Z walls).")
    print(f"  With reflections the path length s > 1400 mm, so t_hit > {t_traverse:.4f} ns.")
    print()
    print("Consistency check:")
    v_g_check = C_MM_NS / N_EJ204
    tol = 0.001  # mm/ns
    ok = abs(v_g_check - 189.742) < tol
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] c/n = {v_g_check:.4f} mm/ns  (want 189.742 ± {tol} mm/ns)")
    print("=" * 60)
    return 0 if ok else 1


def simulation_check(root_path: str) -> int:
    """Mode B: read a ROOT file from a single-photon axial run, verify t_hit."""
    try:
        import ROOT  # type: ignore
    except ImportError:
        print("ERROR: PyROOT not available. Run Mode A (--theoretical) instead.")
        return 1

    f = ROOT.TFile.Open(root_path, "READ")
    if not f or f.IsZombie():
        print(f"ERROR: Cannot open {root_path}")
        return 1

    tree = f.Get("photon_hits")
    if not tree:
        print("ERROR: TTree 'photon_hits' not found in file")
        f.Close()
        return 1

    t_hits = [entry.t_hit_ns for entry in tree if hasattr(entry, "t_hit_ns")]
    if not t_hits:
        print("ERROR: Branch 't_hit_ns' not found or empty")
        f.Close()
        return 1

    f.Close()

    # Filter direct photons: t < 1.5 × t_straight (remove multiply-scattered)
    v_g = C_MM_NS / N_EJ204
    t_straight = (2 * BAR_HALF_LENGTH_MM) / v_g
    direct = [t for t in t_hits if t < 1.5 * t_straight]

    if not direct:
        print(f"ERROR: No photon hits found with t < {1.5*t_straight:.2f} ns")
        return 1

    t_mean = sum(direct) / len(direct)
    t_std  = math.sqrt(sum((t - t_mean)**2 for t in direct) / len(direct))
    v_meas = (2 * BAR_HALF_LENGTH_MM) / t_mean

    print("=" * 60)
    print("EJ-204 group velocity — simulation check")
    print("=" * 60)
    print(f"  Total hits: {len(t_hits)}   Direct (t < {1.5*t_straight:.2f} ns): {len(direct)}")
    print(f"  Mean t_hit = {t_mean:.4f} ± {t_std:.4f} ns")
    print(f"  v_measured = 1400 mm / {t_mean:.4f} ns = {v_meas:.4f} mm/ns")
    print(f"  v_expected = {v_g:.4f} mm/ns")

    tol_frac = 0.01   # 1% tolerance
    ok = abs(v_meas - v_g) / v_g < tol_frac
    status = "PASS" if ok else "FAIL"
    dev = (v_meas - v_g) / v_g * 100
    print(f"  [{status}] deviation = {dev:+.2f}%  (want < ±{tol_frac*100:.0f}%)")
    print("=" * 60)
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--theoretical", action="store_true",
                       help="Theoretical check only (no simulation needed)")
    group.add_argument("--root", metavar="FILE",
                       help="Path to ROOT file from single-photon axial run")
    args = ap.parse_args()

    if args.theoretical:
        return theoretical_check()
    else:
        return simulation_check(args.root)


if __name__ == "__main__":
    sys.exit(main())
