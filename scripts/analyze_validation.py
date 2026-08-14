#!/usr/bin/env python3
"""
Validation analysis for feat/endtop-sslg4 physics check.

Checks:
  1. Npe at END SiPMs per event (expect >10; old physics gave ~0.3)
  2. Effective optical speed v_eff from timing slope (expect 15-19 cm/ns)
  3. Top SiPM Npe (sanity: expect non-zero for EndTop readout)

Usage:
  python3 scripts/analyze_validation.py runs/validation_YYYYMMDD_HHMMSS
"""

from __future__ import annotations
import sys, pathlib, textwrap
import numpy as np
import uproot

# ── Physical constants / expectations ────────────────────────────────────────
BAR_HALF_X_MM   = 700.0          # bar half-length along x
N_BAR           = 1.58           # refractive index EJ-204
C_LIGHT_CM_NS   = 29.9792458     # c in cm/ns
V_GROUP_EXP     = C_LIGHT_CM_NS / N_BAR   # ~18.97 cm/ns  (upper bound: bare GROUPVEL)
V_EFF_EXP_LOW   = 14.0           # cm/ns  lower: scattering, geometry losses
V_EFF_EXP_HIGH  = V_GROUP_EXP    # cm/ns  upper: bare GROUPVEL (no scattering)
NPE_END_MIN     = 5.0            # PE/event minimum at each END face
NPE_END_OLD     = 0.37           # PE/event with old dielectric_metal physics (for reference)

def load_run(run_dir: pathlib.Path) -> dict[str, np.ndarray]:
    """Load all ROOT files in run_dir into a single hit table."""
    files = sorted(run_dir.rglob("*.root"))
    if not files:
        sys.exit(f"ERROR: no ROOT files found under {run_dir}")

    cols: dict[str, list] = {k: [] for k in
        ["event_id", "face_type", "global_id", "time_ns", "gun_x_mm"]}

    for f in files:
        with uproot.open(f) as rf:
            tree = rf["sipm_hits"]
            arr = tree.arrays(list(cols.keys()), library="np")
            for k in cols:
                cols[k].append(arr[k])

    return {k: np.concatenate(v) for k, v in cols.items()}

def first_photon_times(hits: dict, face: int
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Return (gun_x_mm, first_photon_time_ns) per event for a given face."""
    mask = hits["face_type"] == face
    ev   = hits["event_id"][mask]
    t    = hits["time_ns"][mask]
    x    = hits["gun_x_mm"][mask]

    # group by event_id, take min time
    events = np.unique(ev)
    t_first = np.empty(len(events))
    x_first = np.empty(len(events))
    for i, eid in enumerate(events):
        sel = ev == eid
        t_first[i] = t[sel].min()
        x_first[i] = x[sel][0]
    return x_first, t_first

def npe_per_event(hits: dict, face: int) -> dict[float, float]:
    """Mean Npe per event grouped by gun_x_mm position."""
    mask = hits["face_type"] == face
    ev   = hits["event_id"][mask]
    x    = hits["gun_x_mm"][mask]

    result = {}
    for pos in np.unique(x):
        pmask   = x == pos
        n_hits  = len(ev[pmask])
        n_evts  = len(np.unique(hits["event_id"][hits["gun_x_mm"] == pos]))
        result[float(pos)] = n_hits / max(n_evts, 1)
    return result

def veff_from_slope(x_mm: np.ndarray, dt_ns: np.ndarray) -> tuple[float, float]:
    """
    Fit Δt = a·x + b  where  a = -2 / v_eff  (x in mm, Δt in ns).
    Returns (v_eff_cm_ns, fit_residual_std_ns).
    """
    A = np.column_stack([x_mm, np.ones_like(x_mm)])
    result = np.linalg.lstsq(A, dt_ns, rcond=None)
    a, _ = result[0]
    dt_fit = A @ result[0]
    residuals = dt_ns - dt_fit
    # a = -2/v_eff  (x in mm → need /10 for cm)
    v_eff_mm_ns = -2.0 / a  if abs(a) > 1e-9 else float("nan")
    v_eff_cm_ns = v_eff_mm_ns / 10.0
    return v_eff_cm_ns, float(np.std(residuals))

def check(label: str, value: float, low: float, high: float,
          unit: str = "", fmt: str = ".2f") -> str:
    ok = low <= value <= high
    status = "✓ OK  " if ok else "✗ FAIL"
    return (f"  {status}  {label}: {value:{fmt}} {unit}"
            f"  [expected {low:{fmt}}–{high:{fmt}} {unit}]")

# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:
    if len(sys.argv) < 2:
        sys.exit("Usage: analyze_validation.py <run_dir>")
    run_dir = pathlib.Path(sys.argv[1])

    print("=" * 62)
    print(f"  Validation Analysis — feat/endtop-sslg4")
    print(f"  Run dir: {run_dir.name}")
    print("=" * 62)

    hits = load_run(run_dir)
    positions = np.unique(hits["gun_x_mm"])
    n_events_total = len(np.unique(hits["event_id"]))
    print(f"\n  Hits loaded  : {len(hits['event_id']):,}")
    print(f"  Events total : {n_events_total:,}")
    print(f"  Positions    : {len(positions)} "
          f"({positions[0]:.0f} .. {positions[-1]:.0f} mm)")

    # ── 1. Npe per face ──────────────────────────────────────────────────────
    print("\n── 1. Npe per event ─────────────────────────────────────────────")
    print(f"  (reference: old dielectric_metal gave ~{NPE_END_OLD:.2f} PE/event at END)")

    npe_left  = npe_per_event(hits, face=0)
    npe_right = npe_per_event(hits, face=1)
    npe_top   = npe_per_event(hits, face=2)

    mean_left  = np.mean(list(npe_left.values()))  if npe_left  else 0.0
    mean_right = np.mean(list(npe_right.values())) if npe_right else 0.0
    mean_top   = np.mean(list(npe_top.values()))   if npe_top   else 0.0

    print(check("END-left  mean Npe", mean_left,  NPE_END_MIN, 200, "PE/ev"))
    print(check("END-right mean Npe", mean_right, NPE_END_MIN, 200, "PE/ev"))
    print(f"  ─        TOP       mean Npe: {mean_top:.1f} PE/ev  (sanity check)")

    # ── 2. v_eff from timing ─────────────────────────────────────────────────
    print("\n── 2. Effective optical speed (group velocity check) ────────────")
    print("  Method: slope of (t_right_first − t_left_first) vs x")
    print(f"  Expected v_eff: {V_EFF_EXP_LOW:.1f}–{V_EFF_EXP_HIGH:.1f} cm/ns")
    print(f"  (c_light would give {C_LIGHT_CM_NS:.1f} cm/ns → bug if v_eff >> {V_EFF_EXP_HIGH:.1f})")

    x_L, t_L = first_photon_times(hits, face=0)
    x_R, t_R = first_photon_times(hits, face=1)

    # Match events that have BOTH faces hit
    ev_all = hits["event_id"]
    x_all  = hits["gun_x_mm"]
    dt_rows: list[tuple[float, float]] = []

    for pos in positions:
        pmask = x_all == pos
        evs_at_pos = np.unique(ev_all[pmask])
        for eid in evs_at_pos:
            eL = (hits["face_type"] == 0) & (hits["event_id"] == eid)
            eR = (hits["face_type"] == 1) & (hits["event_id"] == eid)
            if eL.any() and eR.any():
                dt_rows.append((pos, hits["time_ns"][eR].min() - hits["time_ns"][eL].min()))

    if len(dt_rows) < 4:
        print("  WARNING: too few events with both END faces hit — skip v_eff fit")
    else:
        dt_arr  = np.array(dt_rows)
        x_fit   = dt_arr[:, 0]
        dt_fit  = dt_arr[:, 1]

        # Show per-position mean Δt
        print(f"\n  {'Position':>10}  {'N paired':>8}  {'mean Δt (ns)':>12}")
        for pos in positions:
            sel = x_fit == pos
            if sel.any():
                print(f"  {pos:>+10.0f}mm  {sel.sum():>8d}  {dt_fit[sel].mean():>+12.3f}")

        v_eff, residual_std = veff_from_slope(x_fit, dt_fit)
        print()
        print(check("v_eff (from Δt slope)",
                    v_eff, V_EFF_EXP_LOW, V_EFF_EXP_HIGH, "cm/ns"))
        print(f"  ─        fit residual σ : {residual_std:.3f} ns")

        if v_eff > C_LIGHT_CM_NS * 0.85:
            print("\n  ⚠ WARNING: v_eff ≈ c_light → GROUP VELOCITY BUG present in EXEC_23 geometry")
        elif v_eff > V_EFF_EXP_HIGH * 1.15:
            print("\n  ⚠ WARNING: v_eff slightly above GROUPVEL — check geometry")
        else:
            print("\n  ✓ Group velocity looks physical")

    # ── 3. Summary verdict ───────────────────────────────────────────────────
    print("\n── 3. Verdict ───────────────────────────────────────────────────")
    npe_ok    = (mean_left >= NPE_END_MIN) and (mean_right >= NPE_END_MIN)
    veff_ok   = V_EFF_EXP_LOW <= v_eff <= V_EFF_EXP_HIGH * 1.15 if len(dt_rows) >= 4 else None

    if npe_ok and (veff_ok is None or veff_ok):
        print("  ✓ PHYSICS VALID — ready for production scan")
        print("    Next: ./scripts/run_t0minidaq_endtop_scan_5000.sh")
    elif not npe_ok:
        print("  ✗ FAIL — Npe too low at END SiPMs")
        print(f"    END-left: {mean_left:.2f}, END-right: {mean_right:.2f} PE/ev")
        print("    Check: TIR fix in CreateBarSkinReflector (dielectric_dielectric?)")
    else:
        print("  ✗ FAIL — v_eff out of physical range → group velocity bug")
        print("    Check: air-gap panel geometry in DetectorConstruction.cc")

    print("=" * 62)

if __name__ == "__main__":
    main()
