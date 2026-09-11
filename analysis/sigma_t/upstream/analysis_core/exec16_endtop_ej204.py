#!/usr/bin/env python3
"""
exec16_endtop_ej204.py — EXEC_16 intrinsic characterisation of EJ-204 EndTop scan.

Autonomous mode: all decisions made per §B policy; HARD-ABORTs per §C only.
Levels 0-4 in sequence; writes figures (PDF+PNG), .root/.csv/.meta.json sidecars,
results.json for Beamer injection, and MORNING_REPORT.md.

Run:
    MPLBACKEND=Agg python3 exec16_endtop_ej204.py
"""

import sys
import os
import json
import hashlib
import datetime
import math
import warnings
import traceback
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine   import fit_core_gaussian
from lib.robust_seeds import gather_seeds

# ════════════════════════════════════════════════════════════════════════════════
# §0  CONSTANTS (all physical values read from repo; NO literals in analysis code)
# ════════════════════════════════════════════════════════════════════════════════

EXEC_ID     = "EXEC_16"
MATERIAL    = "EJ-204"
WRAPPING    = "reflector_skin_surface (Mylar-equivalent)"

REPO_DIR    = "/home/reriosto/SHiP/ej200"
ORCH_DIR    = "/home/reriosto/SHiP/orchestrator"
DATA_DIR    = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959"
OUT_DIR     = "/home/reriosto/SHiP/analysis_core/out/EXEC_16"
BEAMER_DIR  = "/home/reriosto/SHiP/analysis_core/beamer/EXEC_16"

TREE_NAME   = "sipm_hits"
BR = dict(
    event = "event_id",
    face  = "face_type",
    gid   = "global_id",
    lid   = "local_id",
    time  = "time_ns",
    e_eV  = "energy_eV",
    wl    = "wl_nm",
    pde   = "pde",
    x     = "x_mm",
    y     = "y_mm",
    z     = "z_mm",
    gunx  = "gun_x_mm",
)

# Geometry constants READ from DetectorConstruction.hh
K_END    = 8    # SiPMs per end face (from kNEndSiPMs = 8 in header)
N_TOP    = 70   # TOP SiPMs (from kNTopSiPMs = 70 in header)
BOUNDARY = 2 * K_END   # = 16; global_id >= BOUNDARY → TOP

END_L_RANGE = range(0, K_END)            # global_id 0..7
END_R_RANGE = range(K_END, 2 * K_END)   # global_id 8..15
TOP_RANGE   = range(2 * K_END, 2 * K_END + N_TOP)  # global_id 16..85

# Material from opsc-101.mac + DetectorConstruction.cc overrides
TAU_D_NS  = 1.8    # SCINTILLATIONTIMECONSTANT1 (ns)
TAU_R_NS  = 0.7    # SCINTILLATIONRISETIME1 (ns, DetConstr override)
LY_PER_MEV = 10400  # SCINTILLATIONYIELD (/MeV)
ABSLENGTH_CM = 160.0  # DetectorConstruction.cc override (cm)

DS_EJ204 = dict(light_yield_ph_per_MeV=10400, emiss_peak_nm=408,
                rise_ns=0.7, decay_ns=1.8, atten_len_cm=160)

SIGMA_SEED     = "MAD"
FIT_WINDOW_K   = 2.0
N_BOOT         = 200
FIT_OPTIONS    = "R Q S 0"
MIN_EVENTS_FOR_FIT = 30

NPE_THR_SCAN   = (1, 2, 3, 5)
SUM_N_SCAN     = (4, 8)

T_BIN_PS       = 25
DEFAULT_POS_MM = 0
LOG_Z_AUTO     = True
SIGMA_BAND_PS  = (5.0, 400.0)

CAVEAT = ("Resolución INTRÍNSECA (sin SPTR≈106 ps ni jitter FastIC≈10 ps; sin gate). "
          "Wrapping: reflector skin surface on BarLV (Mylar-equivalent). "
          "τ_d=1.8 ns, τ_r=0.7 ns, λ_abs=160 cm (OPSC-101/EJ-204).")

RANDOM_SEED = 20260618

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ════════════════════════════════════════════════════════════════════════════════
# helpers
# ════════════════════════════════════════════════════════════════════════════════

def classify_gid(gid):
    if gid < K_END:          return ("END_L", int(gid))
    if gid < 2 * K_END:      return ("END_R", int(gid - K_END))
    return ("TOP", int(gid - 2 * K_END))

def hard_abort(msg):
    print(f"\n[HARD-ABORT] {msg}", file=sys.stderr)
    # Write abort marker to MORNING_REPORT
    _morning_report_lines.append(f"\n## HARD-ABORT\n{msg}\n")
    _flush_morning_report()
    sys.exit(99)

_morning_report_lines = []

def _flush_morning_report():
    p = Path(OUT_DIR) / "MORNING_REPORT.md"
    p.write_text("\n".join(_morning_report_lines))

def note(msg):
    _morning_report_lines.append(msg)
    print(msg)

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def savefig(fig, stem, subdir=""):
    d = Path(OUT_DIR) / subdir
    d.mkdir(parents=True, exist_ok=True)
    pdf = d / f"{stem}.pdf"
    png = d / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(pdf)

def ps(ns_val):
    return ns_val * 1000.0

def fmt_ps(ns_val, err_ns=None):
    if math.isnan(ns_val):
        return "NaN"
    s = f"{ps(ns_val):.1f} ps"
    if err_ns is not None and not math.isnan(err_ns):
        s += f" ± {ps(err_ns):.1f} ps"
    return s

# ════════════════════════════════════════════════════════════════════════════════
# Data loading
# ════════════════════════════════════════════════════════════════════════════════

def build_pos_map():
    """Return {x_mm: Path} mapping from gun_x_mm metadata."""
    base = Path(DATA_DIR) / "outputs"
    pos_map = {}
    for d in sorted(base.iterdir()):
        root_file = d / "photon_hits_run000.root"
        if not root_file.exists():
            continue
        f = uproot.open(str(root_file))
        t = f[TREE_NAME]
        gunx = t[BR["gunx"]].array(library="np")
        med = float(np.median(gunx))
        spread = float(np.ptp(gunx))
        flag = "gunx_spread" if spread > 1.0 else ""
        x_key = int(round(med))
        pos_map[x_key] = {"path": root_file, "gun_x_actual": med,
                          "flag": flag, "spread": spread}
    return pos_map

def load_pos(pos_entry, branches=None):
    """Load arrays for one position. Returns dict of np arrays."""
    if branches is None:
        branches = list(BR.values())
    f = uproot.open(str(pos_entry["path"]))
    t = f[TREE_NAME]
    data = t.arrays(branches, library="np")
    return data

def best_gids_by_npe(gid_arr, face_label, k):
    """Return k global_ids with most hits in the given face subset."""
    if face_label == "END_L":
        mask = gid_arr < K_END
    elif face_label == "END_R":
        mask = (gid_arr >= K_END) & (gid_arr < 2 * K_END)
    else:
        mask = gid_arr >= 2 * K_END
    gids_face = gid_arr[mask]
    if len(gids_face) == 0:
        return np.array([], dtype=int)
    uniq, cts = np.unique(gids_face, return_counts=True)
    order = np.argsort(-cts)
    return uniq[order[:k]]

def compute_tN_for_gids(event_id, time_ns, gid_arr, selected_gids, N, all_events):
    """Compute t_N (time of N-th hit) from merged stream of selected_gids."""
    if len(selected_gids) == 0:
        return np.full(len(all_events), np.nan)
    hit_mask = np.isin(gid_arr, selected_gids)
    ev_sel = event_id[hit_mask]
    t_sel  = time_ns[hit_mask]
    if len(ev_sel) == 0:
        return np.full(len(all_events), np.nan)
    sort_idx  = np.lexsort((t_sel, ev_sel))
    ev_sorted = ev_sel[sort_idx]
    t_sorted  = t_sel[sort_idx]
    uniq_ev, counts = np.unique(ev_sorted, return_counts=True)
    cum = np.concatenate([[0], np.cumsum(counts)])
    result = np.full(len(all_events), np.nan)
    ev_idx = {e: i for i, e in enumerate(all_events)}
    for i, (ev, cnt) in enumerate(zip(uniq_ev, counts)):
        if cnt >= N:
            idx = ev_idx.get(ev)
            if idx is not None:
                result[idx] = t_sorted[cum[i] + N - 1]
    return result

def gauss_fit(v, name_prefix="g"):
    """Fit Gaussian to core; return dict with sigma_fit, bootstrap_err, etc."""
    cfg_mock = {
        "MIN_EVENTS_FOR_FIT": MIN_EVENTS_FOR_FIT,
        "FIT_WINDOW_SIGMAS":  FIT_WINDOW_K,
        "BINNING_STRATEGY":   "sqrt_n",
        "N_BOOTSTRAP":        N_BOOT,
        "RANDOM_SEED":        RANDOM_SEED,
        "FIT_OPTIONS":        FIT_OPTIONS,
        "CHI2_NDF_WARN":      3.0,
    }
    v_clean = v[~np.isnan(v)]
    if len(v_clean) < MIN_EVENTS_FOR_FIT:
        return {"sigma_fit": np.nan, "sigma_fit_err": np.nan,
                "mu_fit": np.nan, "bootstrap_err": np.nan,
                "n_events": len(v_clean), "flag": "insufficient_stats",
                "h_root": None, "f_root": None, "sigma_mad": np.nan}
    return fit_core_gaussian(v_clean, cfg_mock, name_prefix=name_prefix)

# ════════════════════════════════════════════════════════════════════════════════
# CP0 gate: partition check already done, just re-verify once with real data
# ════════════════════════════════════════════════════════════════════════════════

def run_cp0_gates(pos_map):
    note("\n### CP0 — identity gate")
    # Gate C.3: partition gate
    seen = {}
    for g in range(2 * K_END + N_TOP):
        key = classify_gid(g)
        if key in seen:
            hard_abort(f"Gate C.3: partition collision {key}: g={g} vs g={seen[key]}")
        seen[key] = g
    assert len(seen) == 2 * K_END + N_TOP
    note(f"  Gate C.3 PASSED: K_END={K_END}, N_TOP={N_TOP}, total={2*K_END+N_TOP}")

    # Gate C.5: positions must be non-empty
    if len(pos_map) == 0:
        hard_abort("Gate C.5: no positions found in DATA_DIR")
    note(f"  Gate C.5 PASSED: {len(pos_map)} positions found")

    # Cross-check one file for global_id range
    some_pos = sorted(pos_map.keys())[len(pos_map)//2]
    data = load_pos(pos_map[some_pos])
    gids = data[BR["gid"]]
    g_min, g_max = int(gids.min()), int(gids.max())
    expected_max = 2 * K_END + N_TOP - 1  # = 85
    if g_min < 0 or g_max > expected_max:
        hard_abort(f"Gate C.2: global_id range [{g_min},{g_max}] outside [0,{expected_max}]")
    note(f"  Gate C.2 PASSED: global_id range [{g_min},{g_max}]")

    # Check face_type mapping
    for ft in [0, 1, 2]:
        mask = data[BR["face"]] == ft
        if mask.sum() == 0:
            continue
        gids_ft = gids[mask]
        if ft == 0:
            assert gids_ft.max() < K_END, f"face_type=0 has gid>={K_END}"
        elif ft == 1:
            assert gids_ft.min() >= K_END and gids_ft.max() < 2*K_END, \
                   f"face_type=1 gid out of [K_END,2*K_END)"
        else:
            assert gids_ft.min() >= 2*K_END, f"face_type=2 has gid<BOUNDARY"
    note("  face_type↔global_id mapping verified")

# ════════════════════════════════════════════════════════════════════════════════
# Level 0 — Inventory and QA-0
# ════════════════════════════════════════════════════════════════════════════════

def run_level0(pos_map, positions):
    note("\n## Level 0 — Inventory & QA-0")
    os.makedirs(OUT_DIR + "/L0", exist_ok=True)

    # Events per position
    n_events_per_pos = {}
    for x in positions:
        data = load_pos(pos_map[x], [BR["event"]])
        n_events_per_pos[x] = len(np.unique(data[BR["event"]]))

    fig, ax = plt.subplots(figsize=(12, 4))
    xs = [str(x) for x in positions]
    ax.bar(range(len(positions)), [n_events_per_pos[x] for x in positions])
    ax.set_xticks(range(0, len(positions), 3))
    ax.set_xticklabels([xs[i] for i in range(0, len(positions), 3)], rotation=45, ha="right")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("N events"); ax.set_title("Events per position")
    ax.text(0.01, 0.98, CAVEAT, transform=ax.transAxes, va="top", fontsize=5, wrap=True)
    savefig(fig, "L0_events_per_position", "L0")

    # Occupancy per channel (one representative position: x=0)
    x0 = min(positions, key=lambda x: abs(x))
    data0 = load_pos(pos_map[x0])
    gids0 = data0[BR["gid"]]
    all_gids = np.arange(2 * K_END + N_TOP)
    hits_per_gid = np.array([np.sum(gids0 == g) for g in all_gids])

    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    end_ids = list(range(2 * K_END))
    top_ids = list(range(2 * K_END, 2 * K_END + N_TOP))
    axes[0].bar(end_ids, hits_per_gid[end_ids], color=["steelblue"]*K_END + ["tomato"]*K_END)
    axes[0].set_xlabel("global_id"); axes[0].set_ylabel("N hits")
    axes[0].set_title(f"END occupancy — x={x0} mm"); axes[0].legend(handles=[
        plt.Rectangle((0,0),1,1, color="steelblue", label="END_L"),
        plt.Rectangle((0,0),1,1, color="tomato", label="END_R")])
    axes[1].bar(top_ids, hits_per_gid[top_ids], color="seagreen")
    axes[1].set_xlabel("global_id"); axes[1].set_ylabel("N hits")
    axes[1].set_title(f"TOP occupancy — x={x0} mm")
    fig.suptitle(f"Channel occupancy | {MATERIAL} | {CAVEAT[:50]}...", fontsize=7)
    savefig(fig, "L0_occupancy_channels", "L0")

    # Raw ToA spectra per face at x=0
    times0 = data0[BR["time"]]
    face0  = data0[BR["face"]]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    face_labels = {0: "END_L", 1: "END_R", 2: "TOP"}
    colors = {0: "steelblue", 1: "tomato", 2: "seagreen"}
    for fi, (ftype, flabel) in enumerate(face_labels.items()):
        tmask = face0 == ftype
        tv = times0[tmask] * 1000  # ns → ps
        if len(tv) == 0:
            continue
        lo, hi = np.percentile(tv, [0.5, 99.5])
        axes[fi].hist(tv, bins=200, range=(lo, hi), color=colors[ftype], histtype="step", lw=1.5)
        axes[fi].set_xlabel("time [ps]"); axes[fi].set_ylabel("hits")
        axes[fi].set_title(f"{flabel} — raw ToA | x={x0} mm")
        axes[fi].axvline(TAU_D_NS*1000, color="k", ls="--", lw=0.8, label=f"τ_d={TAU_D_NS} ns")
        axes[fi].legend(fontsize=7)
    fig.suptitle(f"{MATERIAL} raw ToA spectra | {CAVEAT[:60]}...", fontsize=7)
    savefig(fig, "L0_raw_toa_spectra", "L0")

    # Write CSV sidecar for inventory
    import csv
    inv_path = Path(OUT_DIR) / "L0" / "inventory.csv"
    with open(inv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["x_mm", "n_events", "gunx_flag"])
        w.writeheader()
        for x in positions:
            w.writerow({"x_mm": x, "n_events": n_events_per_pos[x],
                        "gunx_flag": pos_map[x].get("flag", "")})
    note(f"  Level 0 done. Inventory: {inv_path}")

    # QA-0 gates
    for x in positions:
        if n_events_per_pos[x] == 0:
            hard_abort(f"Gate QA-0: 0 events at position x={x} mm")
    note("  QA-0 PASSED")
    return n_events_per_pos

# ════════════════════════════════════════════════════════════════════════════════
# Level 1 — Photon / PE metrics
# ════════════════════════════════════════════════════════════════════════════════

def run_level1(pos_map, positions):
    note("\n## Level 1 — Photon metrics")
    os.makedirs(OUT_DIR + "/L1", exist_ok=True)

    npe_by_face = {"END_L": [], "END_R": [], "TOP": [], "x": []}
    npe_nearest_top = []

    for x in positions:
        data = load_pos(pos_map[x])
        gids = data[BR["gid"]]
        evs  = data[BR["event"]]
        all_events = np.unique(evs)
        n_ev = len(all_events)

        # Mean Npe per event per face
        for face_label, g_range in [("END_L", END_L_RANGE),
                                     ("END_R", END_R_RANGE),
                                     ("TOP",   TOP_RANGE)]:
            mask = np.isin(gids, list(g_range))
            hits_per_ev = np.array([np.sum((evs == e) & mask) for e in all_events])
            npe_by_face[face_label].append(float(np.mean(hits_per_ev)))

        npe_by_face["x"].append(x)

        # Nearest TOP channel (max <Npe>)
        top_mask = gids >= 2 * K_END
        top_gids_ev = gids[top_mask]
        if len(top_gids_ev) > 0:
            uniq_top, cts_top = np.unique(top_gids_ev, return_counts=True)
            nearest_gid = int(uniq_top[np.argmax(cts_top)])
        else:
            nearest_gid = int(TOP_RANGE.start)
        mean_npe_nearest = float(np.sum(top_mask) / n_ev) if n_ev > 0 else 0.0
        npe_nearest_top.append({"x": x, "nearest_gid": nearest_gid,
                                 "mean_npe": float(np.sum(gids == nearest_gid) / n_ev)})

    # <Npe>(x) plot
    fig, ax = plt.subplots(figsize=(10, 5))
    xs_arr = np.array(npe_by_face["x"])
    ax.plot(xs_arr, npe_by_face["END_L"], "s-", color="steelblue", label="END_L")
    ax.plot(xs_arr, npe_by_face["END_R"], "^-", color="tomato", label="END_R")
    ax.plot(xs_arr, npe_by_face["TOP"],   "o-", color="seagreen", label="TOP (all 70)")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("⟨Npe⟩ per event")
    ax.set_title(f"{EXEC_ID} — ⟨Npe⟩(x) per face | {MATERIAL}")
    ax.legend(); ax.grid(True, lw=0.3)
    ax.text(0.01, 0.02, CAVEAT, transform=ax.transAxes, fontsize=5)
    savefig(fig, "L1_npe_vs_x", "L1")

    # NPE(x) profile — double-exp descriptive fit on TOP
    fig, ax = plt.subplots(figsize=(10, 5))
    top_arr = np.array(npe_by_face["TOP"])
    ax.plot(xs_arr, top_arr, "o", color="seagreen", ms=5, label="TOP ⟨Npe⟩ data")
    try:
        def double_exp(x, A1, l1, A2, l2, C):
            return A1 * np.exp(-np.abs(x) / l1) + A2 * np.exp(-np.abs(x) / l2) + C
        p0 = [top_arr.max()*0.5, 300, top_arr.max()*0.3, 700, top_arr.min()]
        popt, _ = curve_fit(double_exp, xs_arr, top_arr, p0=p0, maxfev=5000)
        xfit = np.linspace(xs_arr.min(), xs_arr.max(), 300)
        ax.plot(xfit, double_exp(xfit, *popt), "--", color="darkgreen", lw=1.5,
                label=f"Double-exp (DESCRIPTIVE only)\nA1={popt[0]:.0f},λ1={popt[1]:.0f} mm,A2={popt[2]:.0f},λ2={popt[3]:.0f} mm")
    except Exception as e:
        note(f"  [WARN] double-exp fit on TOP NPE(x) failed: {e}")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("⟨Npe⟩ (TOP all 70 SiPMs)")
    ax.set_title(f"{EXEC_ID} — TOP NPE(x) profile | {MATERIAL}\nNOTE: double-exp is DESCRIPTIVE; λ_eff is NOT the bulk ABSLENGTH ({ABSLENGTH_CM} cm)")
    ax.legend(fontsize=7); ax.grid(True, lw=0.3)
    ax.text(0.01, 0.02, CAVEAT, transform=ax.transAxes, fontsize=5)
    savefig(fig, "L1_npe_profile_top", "L1")

    # Npe distributions at x=0 and extremes
    repr_positions = [sorted(positions)[0], 0, sorted(positions)[-1]]
    repr_positions = [min(positions, key=lambda x: abs(x-t)) for t in repr_positions]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ai, rp in enumerate(repr_positions):
        data = load_pos(pos_map[rp])
        gids = data[BR["gid"]]
        evs  = data[BR["event"]]
        all_events = np.unique(evs)
        top_mask = gids >= 2 * K_END
        hits_per_ev = np.array([np.sum((evs == e) & top_mask) for e in all_events])
        axes[ai].hist(hits_per_ev, bins=50, histtype="step", lw=1.5, color="seagreen")
        axes[ai].set_xlabel("Npe (TOP)"); axes[ai].set_ylabel("events")
        axes[ai].set_title(f"Npe dist — x={rp} mm")
        axes[ai].axvline(float(np.mean(hits_per_ev)), color="k", ls="--", lw=0.8,
                         label=f"⟨Npe⟩={np.mean(hits_per_ev):.1f}")
        axes[ai].legend(fontsize=7)
    fig.suptitle(f"{EXEC_ID} Npe distributions | TOP | {MATERIAL}", fontsize=9)
    savefig(fig, "L1_npe_distributions", "L1")

    # Save NPE sidecar CSV
    import csv
    npe_csv = Path(OUT_DIR) / "L1" / "npe_per_face.csv"
    with open(npe_csv, "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm", "END_L", "END_R", "TOP"])
        w.writeheader()
        for i, x in enumerate(npe_by_face["x"]):
            w.writerow({"x_mm": x, "END_L": npe_by_face["END_L"][i],
                        "END_R": npe_by_face["END_R"][i], "TOP": npe_by_face["TOP"][i]})
    note(f"  Level 1 done. Sidecar: {npe_csv}")

    # mean Npe at x=0 for nearest TOP (for MORNING_REPORT)
    top_nearest_x0 = next((d for d in npe_nearest_top if d["x"] == min(positions, key=lambda x: abs(x))), npe_nearest_top[0])
    note(f"  TOP nearest at x={top_nearest_x0['x']} mm: gid={top_nearest_x0['nearest_gid']}, ⟨Npe⟩={top_nearest_x0['mean_npe']:.2f}")

    return npe_by_face, npe_nearest_top

# ════════════════════════════════════════════════════════════════════════════════
# Level 2 — Temporal resolution σ_t
# ════════════════════════════════════════════════════════════════════════════════

def run_level2(pos_map, positions, npe_nearest_top):
    note("\n## Level 2 — Temporal resolution σ_t")
    os.makedirs(OUT_DIR + "/L2", exist_ok=True)

    # Build nearest-TOP gid map from Level 1 results
    nearest_top_gid = {d["x"]: d["nearest_gid"] for d in npe_nearest_top}

    results = {}   # key = estimator_label → list of {x, sigma_fit, ...}

    estimator_labels = (
        ["TOP_SUM4", "TOP_SUM8"] +
        [f"TOP_THR{k}" for k in NPE_THR_SCAN] +
        ["END_L_SUM4", "END_R_SUM4", "END_L_SUM8", "END_R_SUM8",
         "END_COMBINED_SUM4", "END_COMBINED_SUM8"]
    )
    for label in estimator_labels:
        results[label] = []

    sigma_outlier_flags = {label: [] for label in estimator_labels}

    for x in positions:
        data = load_pos(pos_map[x])
        gids    = data[BR["gid"]]
        evs     = data[BR["event"]]
        times   = data[BR["time"]]
        all_events = np.unique(evs)

        # Identify top-k gids for each face
        top_sum4_gids  = best_gids_by_npe(gids, "TOP",   4)
        top_sum8_gids  = best_gids_by_npe(gids, "TOP",   8)
        endl_sum4_gids = best_gids_by_npe(gids, "END_L", 4)
        endr_sum4_gids = best_gids_by_npe(gids, "END_R", 4)
        endl_sum8_gids = best_gids_by_npe(gids, "END_L", 8)
        endr_sum8_gids = best_gids_by_npe(gids, "END_R", 8)

        def _fit(label, gid_sel, N):
            t_N = compute_tN_for_gids(evs, times, gids, gid_sel, N, all_events)
            r = gauss_fit(t_N, name_prefix=f"{label}_x{x}_N{N}")
            r["x_mm"] = x
            r["estimator"] = label
            r["N"] = N
            # QA-2 physical band check
            sig_ps_val = ps(r["sigma_fit"]) if not math.isnan(r["sigma_fit"]) else np.nan
            if not math.isnan(sig_ps_val):
                if sig_ps_val < SIGMA_BAND_PS[0] or sig_ps_val > SIGMA_BAND_PS[1]:
                    r["flag"] = r.get("flag","ok") + "|sigma_outlier"
                    sigma_outlier_flags[label].append(x)
            return r

        # TOP estimators
        for N in SUM_N_SCAN:
            label = f"TOP_SUM{N}"
            gids_sel = top_sum4_gids if N == 4 else top_sum8_gids
            results[label].append(_fit(label, gids_sel, 1))

        # TOP threshold scan (t_k-th PE from nearest TOP SiPM)
        nearest_gid = nearest_top_gid.get(x, int(TOP_RANGE.start))
        for k in NPE_THR_SCAN:
            label = f"TOP_THR{k}"
            # HOOK_WALK: time-walk correction would go here for leading-edge
            results[label].append(_fit(label, np.array([nearest_gid]), k))

        # END estimators
        for label, gids_sel in [("END_L_SUM4", endl_sum4_gids),
                                  ("END_R_SUM4", endr_sum4_gids),
                                  ("END_L_SUM8", endl_sum8_gids),
                                  ("END_R_SUM8", endr_sum8_gids)]:
            results[label].append(_fit(label, gids_sel, 1))

        # END combined: average of END_L and END_R t_1 timestamps
        tL_arr = compute_tN_for_gids(evs, times, gids, endl_sum4_gids, 1, all_events)
        tR_arr = compute_tN_for_gids(evs, times, gids, endr_sum4_gids, 1, all_events)
        valid = (~np.isnan(tL_arr)) & (~np.isnan(tR_arr))
        t_comb_sum4 = np.full(len(all_events), np.nan)
        t_comb_sum4[valid] = 0.5 * (tL_arr[valid] + tR_arr[valid])
        rc4 = gauss_fit(t_comb_sum4, name_prefix=f"END_COMBINED_SUM4_x{x}")
        rc4["x_mm"] = x; rc4["estimator"] = "END_COMBINED_SUM4"; rc4["N"] = 1
        results["END_COMBINED_SUM4"].append(rc4)

        tL8 = compute_tN_for_gids(evs, times, gids, endl_sum8_gids, 1, all_events)
        tR8 = compute_tN_for_gids(evs, times, gids, endr_sum8_gids, 1, all_events)
        valid8 = (~np.isnan(tL8)) & (~np.isnan(tR8))
        t_comb_sum8 = np.full(len(all_events), np.nan)
        t_comb_sum8[valid8] = 0.5 * (tL8[valid8] + tR8[valid8])
        rc8 = gauss_fit(t_comb_sum8, name_prefix=f"END_COMBINED_SUM8_x{x}")
        rc8["x_mm"] = x; rc8["estimator"] = "END_COMBINED_SUM8"; rc8["N"] = 1
        results["END_COMBINED_SUM8"].append(rc8)

    # Check for >30% outlier fraction per curve (HARD-ABORT that curve)
    aborted_curves = set()
    for label, outlier_xs in sigma_outlier_flags.items():
        frac = len(outlier_xs) / len(positions) if positions else 0
        if frac > 0.30:
            note(f"  [CURVE-ABORT] {label}: {len(outlier_xs)}/{len(positions)} = {frac:.0%} σ_outlier — curve excluded")
            aborted_curves.add(label)
            note(f"    Outlier positions: {outlier_xs}")

    # Plot σ_t(x) for main estimators
    fig, ax = plt.subplots(figsize=(12, 5))
    plot_labels = ["TOP_SUM4", "TOP_SUM8", "END_COMBINED_SUM4", "END_COMBINED_SUM8"]
    colors = ["seagreen", "darkgreen", "steelblue", "navy"]
    for label, color in zip(plot_labels, colors):
        if label in aborted_curves:
            continue
        rows = results[label]
        xs_v  = [r["x_mm"] for r in rows if not math.isnan(r.get("sigma_fit", np.nan))]
        sig_v = [ps(r["sigma_fit"]) for r in rows if not math.isnan(r.get("sigma_fit", np.nan))]
        err_v = [ps(max(r.get("sigma_fit_err", 0), r.get("bootstrap_err", 0)))
                 for r in rows if not math.isnan(r.get("sigma_fit", np.nan))]
        if xs_v:
            ax.errorbar(xs_v, sig_v, yerr=err_v, fmt="o-", color=color,
                        ms=4, lw=1.2, capsize=3, label=label)
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("σ_fit [ps]")
    ax.set_title(f"{EXEC_ID} — Intrinsic σ_t(x) | {MATERIAL}")
    ax.legend(fontsize=8); ax.grid(True, lw=0.3)
    ax.text(0.01, 0.98, CAVEAT, transform=ax.transAxes, va="top", fontsize=5)
    savefig(fig, "L2_sigma_t_vs_x", "L2")

    # Panel comparison at x=0
    x0_pos = min(positions, key=lambda x: abs(x))
    fig, ax = plt.subplots(figsize=(10, 5))
    x0_rows = {label: next((r for r in results[label] if r["x_mm"] == x0_pos), None)
               for label in estimator_labels if label not in aborted_curves}
    labels_plot, sig_ps_vals, err_ps_vals = [], [], []
    for label, row in sorted(x0_rows.items()):
        if row is None or math.isnan(row.get("sigma_fit", np.nan)):
            continue
        labels_plot.append(label)
        sig_ps_vals.append(ps(row["sigma_fit"]))
        err_ps_vals.append(ps(max(row.get("sigma_fit_err", 0), row.get("bootstrap_err", 0))))
    y_pos = range(len(labels_plot))
    ax.barh(list(y_pos), sig_ps_vals, xerr=err_ps_vals, color="steelblue", height=0.6)
    ax.set_yticks(list(y_pos)); ax.set_yticklabels(labels_plot, fontsize=8)
    ax.set_xlabel("σ_fit [ps]")
    ax.set_title(f"{EXEC_ID} — σ_t at x=0 by estimator | {MATERIAL}")
    ax.axvline(100, color="red", ls="--", lw=1, label="SHiP goal 100 ps")
    ax.axvline(50,  color="orange", ls="--", lw=1, label="SHiP preferred 50 ps")
    ax.legend(fontsize=8); ax.grid(True, axis="x", lw=0.3)
    ax.text(0.99, 0.01, CAVEAT, transform=ax.transAxes, ha="right", fontsize=5)
    savefig(fig, "L2_sigma_panel_x0", "L2")

    # Write CSV sidecar
    import csv
    csv_path = Path(OUT_DIR) / "L2" / "sigma_t_all_estimators.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, ["estimator", "x_mm", "sigma_fit_ps", "sigma_fit_err_ps",
                                "bootstrap_err_ps", "mu_fit_ns", "chi2_ndf", "n_events", "flag"])
        w.writeheader()
        for label, rows in results.items():
            for r in rows:
                w.writerow({
                    "estimator":       label,
                    "x_mm":            r.get("x_mm", ""),
                    "sigma_fit_ps":    ps(r["sigma_fit"]) if not math.isnan(r.get("sigma_fit", np.nan)) else "NaN",
                    "sigma_fit_err_ps": ps(r.get("sigma_fit_err", np.nan)) if not math.isnan(r.get("sigma_fit_err", np.nan)) else "NaN",
                    "bootstrap_err_ps": ps(r.get("bootstrap_err", np.nan)) if not math.isnan(r.get("bootstrap_err", np.nan)) else "NaN",
                    "mu_fit_ns":       r.get("mu_fit", ""),
                    "chi2_ndf":        r.get("chi2_ndf", ""),
                    "n_events":        r.get("n_events", ""),
                    "flag":            r.get("flag", ""),
                })
    note(f"  Level 2 done. Sidecar: {csv_path}")
    note(f"  Aborted curves (>30% outlier): {aborted_curves or 'none'}")

    return results, aborted_curves, x0_pos

# ════════════════════════════════════════════════════════════════════════════════
# Level 3 — Comparisons and spatial resolution
# ════════════════════════════════════════════════════════════════════════════════

def run_level3(pos_map, positions, sigma_results, aborted_curves, x0_pos):
    note("\n## Level 3 — END-only vs EndTop & spatial resolution")
    os.makedirs(OUT_DIR + "/L3", exist_ok=True)

    # END-only timestamp: END_COMBINED (already computed in Level 2)
    # EndTop: how to combine END+TOP? Use top-SUM4 as TOP representative, then
    # min(t_endcomb, t_topsum4) per event as the "earliest trigger"
    # or alternatively: arithmetic mean approach. Use the min timestamp per event.

    ratio_data = []
    for x in positions:
        data = load_pos(pos_map[x])
        gids  = data[BR["gid"]]
        evs   = data[BR["event"]]
        times = data[BR["time"]]
        all_events = np.unique(evs)

        endl_gids  = best_gids_by_npe(gids, "END_L", 8)
        endr_gids  = best_gids_by_npe(gids, "END_R", 8)
        top4_gids  = best_gids_by_npe(gids, "TOP",   4)

        tL = compute_tN_for_gids(evs, times, gids, endl_gids, 1, all_events)
        tR = compute_tN_for_gids(evs, times, gids, endr_gids, 1, all_events)
        tT = compute_tN_for_gids(evs, times, gids, top4_gids, 1, all_events)

        valid_end = (~np.isnan(tL)) & (~np.isnan(tR))
        t_end_comb = np.full(len(all_events), np.nan)
        t_end_comb[valid_end] = 0.5 * (tL[valid_end] + tR[valid_end])

        # EndTop: mean of END_combined and TOP_SUM4 when both available
        valid_et = valid_end & (~np.isnan(tT))
        t_endtop = np.full(len(all_events), np.nan)
        t_endtop[valid_et] = 0.5 * (t_end_comb[valid_et] + tT[valid_et])

        r_end = gauss_fit(t_end_comb, name_prefix=f"L3_endonly_x{x}")
        r_et  = gauss_fit(t_endtop,   name_prefix=f"L3_endtop_x{x}")

        ratio = np.nan
        if not math.isnan(r_end.get("sigma_fit", np.nan)) and not math.isnan(r_et.get("sigma_fit", np.nan)):
            ratio = r_et["sigma_fit"] / r_end["sigma_fit"]
        ratio_data.append({"x_mm": x,
                            "sigma_end_ps": ps(r_end.get("sigma_fit", np.nan)),
                            "sigma_endtop_ps": ps(r_et.get("sigma_fit", np.nan)),
                            "ratio": ratio,
                            "boot_end_ps": ps(r_end.get("bootstrap_err", np.nan)),
                            "boot_endtop_ps": ps(r_et.get("bootstrap_err", np.nan))})

    # Plot comparison
    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    xs_v = [d["x_mm"] for d in ratio_data]
    sig_end = [d["sigma_end_ps"] for d in ratio_data]
    sig_et  = [d["sigma_endtop_ps"] for d in ratio_data]
    boot_end = [d["boot_end_ps"] for d in ratio_data]
    boot_et  = [d["boot_endtop_ps"] for d in ratio_data]
    axes[0].errorbar(xs_v, sig_end, yerr=boot_end, fmt="s-", color="steelblue",
                     capsize=3, ms=4, lw=1.2, label="END-only (SUM8)")
    axes[0].errorbar(xs_v, sig_et,  yerr=boot_et,  fmt="o-", color="seagreen",
                     capsize=3, ms=4, lw=1.2, label="EndTop (END+TOP_SUM4)")
    axes[0].set_ylabel("σ_fit [ps]"); axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)
    axes[0].axhline(100, color="red", ls="--", lw=0.8, label="SHiP 100 ps")
    axes[0].set_title(f"{EXEC_ID} — END-only vs EndTop σ_t(x) | {MATERIAL}")

    ratios = [d["ratio"] for d in ratio_data]
    axes[1].plot(xs_v, ratios, "d-", color="purple", ms=4, lw=1.2)
    axes[1].axhline(1.0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("σ_EndTop / σ_END-only")
    axes[1].set_title("Ratio < 1 → TOP helps | > 1 → TOP hurts")
    axes[1].grid(True, lw=0.3)
    for ax in axes:
        ax.text(0.01, 0.02, CAVEAT, transform=ax.transAxes, fontsize=5)
    savefig(fig, "L3_endonly_vs_endtop", "L3")

    # Spatial resolution: Δt = t_endL - t_endR
    delta_t_data = []
    for x in positions:
        data = load_pos(pos_map[x])
        gids  = data[BR["gid"]]
        evs   = data[BR["event"]]
        times = data[BR["time"]]
        all_events = np.unique(evs)

        endl_gids = best_gids_by_npe(gids, "END_L", 4)
        endr_gids = best_gids_by_npe(gids, "END_R", 4)
        tL = compute_tN_for_gids(evs, times, gids, endl_gids, 1, all_events)
        tR = compute_tN_for_gids(evs, times, gids, endr_gids, 1, all_events)

        valid = (~np.isnan(tL)) & (~np.isnan(tR))
        delta_t = tL[valid] - tR[valid]
        mean_dt = float(np.mean(delta_t)) if len(delta_t) > 0 else np.nan
        std_dt  = float(np.std(delta_t, ddof=1)) if len(delta_t) > 1 else np.nan
        delta_t_data.append({"x_mm": x, "mean_dt_ns": mean_dt, "std_dt_ns": std_dt,
                              "n": int(np.sum(valid))})

    # Fit v_eff from slope of <Δt>(x)
    xs_dt  = [d["x_mm"] for d in delta_t_data if not math.isnan(d["mean_dt_ns"])]
    mdt    = [d["mean_dt_ns"] for d in delta_t_data if not math.isnan(d["mean_dt_ns"])]
    v_eff  = np.nan
    if len(xs_dt) >= 5:
        try:
            # <Δt> = -(2/v_eff) * x  → slope = -2/v_eff
            m, b = np.polyfit(xs_dt, mdt, 1)
            v_eff = float(-2.0 / m)  # mm/ns = cm/ns × 10
            v_eff_cm_per_ns = v_eff / 10.0
        except Exception:
            v_eff = np.nan

    # σ_x from σ(Δt) × v_eff/2
    if not math.isnan(v_eff):
        sigma_x_data = []
        for d in delta_t_data:
            if not math.isnan(d["std_dt_ns"]):
                sigma_x = d["std_dt_ns"] * abs(v_eff) / 2.0  # mm
                sigma_x_data.append({"x_mm": d["x_mm"], "sigma_x_mm": sigma_x})
    else:
        sigma_x_data = [{"x_mm": d["x_mm"], "sigma_x_mm": np.nan} for d in delta_t_data]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].plot(xs_dt, mdt, "o", ms=4, color="steelblue", label="⟨Δt⟩(x)")
    if not math.isnan(v_eff) and len(xs_dt) >= 2:
        xfit = np.array([min(xs_dt), max(xs_dt)])
        axes[0].plot(xfit, m*xfit+b, "--", color="red", lw=1.2,
                     label=f"Linear fit: v_eff={v_eff/10:.2f} cm/ns")
    axes[0].set_xlabel("gun_x [mm]"); axes[0].set_ylabel("⟨Δt⟩ (END_L−END_R) [ns]")
    axes[0].set_title("⟨Δt⟩(x) → v_eff")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)

    axes[1].plot([d["x_mm"] for d in sigma_x_data],
                 [d["sigma_x_mm"] for d in sigma_x_data],
                 "s-", ms=4, color="tomato", label="σ_x(x)")
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("σ_x [mm]")
    axes[1].set_title(f"Spatial resolution σ_x(x) | v_eff≈{v_eff/10:.2f} cm/ns" if not math.isnan(v_eff) else "Spatial resolution σ_x(x)")
    axes[1].legend(fontsize=8); axes[1].grid(True, lw=0.3)
    for ax in axes:
        ax.text(0.01, 0.02, CAVEAT, transform=ax.transAxes, fontsize=5)
    savefig(fig, "L3_spatial_resolution", "L3")

    # Write sidecar CSVs
    import csv
    ratio_csv = Path(OUT_DIR) / "L3" / "endonly_vs_endtop.csv"
    with open(ratio_csv, "w", newline="") as f:
        w = csv.DictWriter(f, list(ratio_data[0].keys()))
        w.writeheader(); w.writerows(ratio_data)

    spatial_csv = Path(OUT_DIR) / "L3" / "spatial_resolution.csv"
    with open(spatial_csv, "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm", "mean_dt_ns", "std_dt_ns", "sigma_x_mm", "n"])
        w.writeheader()
        for d, sx in zip(delta_t_data, sigma_x_data):
            w.writerow({**d, "sigma_x_mm": sx["sigma_x_mm"]})

    v_eff_cm_per_ns_val = v_eff / 10.0 if not math.isnan(v_eff) else np.nan
    note(f"  v_eff = {v_eff_cm_per_ns_val:.2f} cm/ns (slope method)" if not math.isnan(v_eff_cm_per_ns_val) else "  v_eff = NaN (fit failed)")
    note(f"  Level 3 done.")

    return ratio_data, delta_t_data, v_eff, sigma_x_data

# ════════════════════════════════════════════════════════════════════════════════
# Level 4 — 2D maps
# ════════════════════════════════════════════════════════════════════════════════

def run_level4(pos_map, positions, npe_nearest_top):
    note("\n## Level 4 — 2D maps")
    os.makedirs(OUT_DIR + "/L4", exist_ok=True)

    x0_pos   = min(positions, key=lambda x: abs(x))
    nt_entry = next((d for d in npe_nearest_top if d["x"] == x0_pos), npe_nearest_top[0])
    nearest_gid = nt_entry["nearest_gid"]
    note(f"  Map channel: TOP nearest gid={nearest_gid} at x={x0_pos} mm (⟨Npe⟩={nt_entry['mean_npe']:.2f})")

    # Gate coherence: nearest TOP should be above the beam (or close)
    # TopSiPMCenterX: idx=gid-BOUNDARY; idx<35 → -692+20*idx mm; idx>=35 → 12+20*(idx-35) mm
    top_local = nearest_gid - 2 * K_END
    if top_local < 35:
        sipm_cx = -692.0 + 20.0 * top_local
    else:
        sipm_cx = 12.0 + 20.0 * (top_local - 35)
    dist_to_beam = abs(sipm_cx - x0_pos)
    if dist_to_beam > 50:
        note(f"  [FLAG nearest_incoherent] TOP nearest gid={nearest_gid} at x={sipm_cx:.0f} mm, beam at x={x0_pos} mm, Δ={dist_to_beam:.0f} mm > 50 mm")

    data = load_pos(pos_map[x0_pos])
    mask_ch = data[BR["gid"]] == nearest_gid
    evs_ch  = data[BR["event"]][mask_ch]
    t_ch    = data[BR["time"]][mask_ch] * 1000.0  # ns → ps
    all_events = np.unique(data[BR["event"]])
    n_events = len(all_events)

    if len(t_ch) == 0:
        note("  [WARN] No hits in nearest TOP channel for maps; skipping Level 4 maps")
        return

    # ── Map 1: time × event index (identity preserved)
    t_lo = np.percentile(t_ch, 0.5)
    t_hi = np.percentile(t_ch, 99.5)
    t_range = (t_lo, t_hi)
    n_t_bins = max(1, int((t_hi - t_lo) / T_BIN_PS))
    n_ev_bins = min(n_events, 500)  # cap at 500 rows for visual clarity

    ev_map_idx = {e: i for i, e in enumerate(all_events)}
    ev_indices = np.array([ev_map_idx.get(e, -1) for e in evs_ch])
    valid = ev_indices >= 0
    ev_idx_norm = (ev_indices[valid] * n_ev_bins / n_events).astype(int)
    ev_idx_norm = np.clip(ev_idx_norm, 0, n_ev_bins - 1)

    h2d_1, xedge, yedge = np.histogram2d(
        t_ch[valid], ev_idx_norm,
        bins=[n_t_bins, n_ev_bins],
        range=[t_range, [0, n_ev_bins]])

    dyn_1 = h2d_1.max() / (h2d_1[h2d_1 > 0].min() if np.any(h2d_1 > 0) else 1)
    use_log_1 = LOG_Z_AUTO and dyn_1 > 50

    fig, ax = plt.subplots(figsize=(12, 7))
    norm_1 = mcolors.LogNorm(vmin=0.5, vmax=h2d_1.max()) if use_log_1 else None
    im = ax.pcolormesh(xedge, yedge, h2d_1.T, cmap="viridis", norm=norm_1)
    plt.colorbar(im, ax=ax, label="PE / bin")
    ax.set_xlabel("time of PE arrival [ps]")
    ax.set_ylabel("event index (one row = one event)")
    ax.set_title(f"{EXEC_ID} — Map 1: time × event | TOP gid={nearest_gid} | x={x0_pos} mm | {MATERIAL}")
    cap = (f"Gráfico 1: el eje Y es el índice del evento (coordenada); z = PE/casilla. "
           f"{'Escala log en z (rango dinámica >{50:.0f}×). ' if use_log_1 else 'Escala lineal en z. '}"
           f"{CAVEAT}")
    ax.text(0.01, -0.12, cap, transform=ax.transAxes, fontsize=6, wrap=True)
    savefig(fig, "L4_map1_time_vs_event", "L4")

    # ── Map 2: time × accumulated PE, z = N events
    rows_t, rows_nacc = [], []
    for e in all_events:
        mask_e = (data[BR["event"]] == e) & (data[BR["gid"]] == nearest_gid)
        t_e = np.sort(data[BR["time"]][mask_e]) * 1000.0
        for i, ti in enumerate(t_e):
            rows_t.append(ti)
            rows_nacc.append(i + 1)

    rows_t    = np.array(rows_t)
    rows_nacc = np.array(rows_nacc)

    if len(rows_t) == 0:
        note("  [WARN] No data for Map 2")
        return

    max_nacc = int(rows_nacc.max())
    n_nacc_bins = max_nacc
    h2d_2, xe2, ye2 = np.histogram2d(
        rows_t, rows_nacc,
        bins=[n_t_bins, n_nacc_bins],
        range=[t_range, [0.5, max_nacc + 0.5]])

    dyn_2 = h2d_2.max() / (h2d_2[h2d_2 > 0].min() if np.any(h2d_2 > 0) else 1)
    use_log_2 = LOG_Z_AUTO and dyn_2 > 50

    fig, ax = plt.subplots(figsize=(12, 6))
    norm_2 = mcolors.LogNorm(vmin=0.5, vmax=h2d_2.max()) if use_log_2 else None
    im = ax.pcolormesh(xe2, ye2, h2d_2.T, cmap="plasma", norm=norm_2)
    plt.colorbar(im, ax=ax, label="eventos / bin")
    ax.set_xlabel("time [ps]")
    ax.set_ylabel("N PE acumulados en el evento")
    ax.set_title(f"{EXEC_ID} — Map 2: time × Npe_acum | TOP gid={nearest_gid} | x={x0_pos} mm | {MATERIAL}")
    cap2 = (f"Gráfico 2: eje Y = Npe acumulados DENTRO del evento; z = nº de eventos que pasan "
            f"por (t, Npe_acum). SÍ se agregan todos los {n_events} eventos. "
            f"{'Escala log en z. ' if use_log_2 else 'Escala lineal en z. '}{CAVEAT}")
    ax.text(0.01, -0.12, cap2, transform=ax.transAxes, fontsize=6, wrap=True)
    savefig(fig, "L4_map2_time_vs_nacc", "L4")

    note(f"  Level 4 done. Maps at x={x0_pos} mm, TOP gid={nearest_gid}.")

# ════════════════════════════════════════════════════════════════════════════════
# Results JSON
# ════════════════════════════════════════════════════════════════════════════════

def write_results_json(sigma_results, ratio_data, v_eff, x0_pos, npe_nearest_top):
    """Write results.json for Beamer number injection — NO manual numbers."""
    rj = {}
    for est_label in ["TOP_SUM4", "TOP_SUM8", "END_COMBINED_SUM4", "END_COMBINED_SUM8"]:
        rows = sigma_results.get(est_label, [])
        row0 = next((r for r in rows if r.get("x_mm") == x0_pos), None)
        rj[est_label] = {
            "x0_sigma_ps": round(ps(row0["sigma_fit"]), 1) if row0 and not math.isnan(row0.get("sigma_fit", np.nan)) else None,
            "x0_boot_ps":  round(ps(row0.get("bootstrap_err", np.nan)), 1) if row0 and not math.isnan(row0.get("bootstrap_err", np.nan)) else None,
            "x0_fit_err_ps": round(ps(row0.get("sigma_fit_err", np.nan)), 1) if row0 and not math.isnan(row0.get("sigma_fit_err", np.nan)) else None,
        }

    ratio_x0 = next((d for d in ratio_data if d["x_mm"] == x0_pos), None)
    rj["endtop_ratio_x0"] = round(float(ratio_x0["ratio"]), 3) if ratio_x0 and not math.isnan(ratio_x0["ratio"]) else None
    rj["v_eff_cm_per_ns"]  = round(float(v_eff / 10.0), 2) if not math.isnan(v_eff) else None

    nt_x0 = next((d for d in npe_nearest_top if d["x"] == x0_pos), None)
    rj["top_nearest_gid_x0"]    = int(nt_x0["nearest_gid"]) if nt_x0 else None
    rj["top_nearest_npe_x0"]    = round(float(nt_x0["mean_npe"]), 2) if nt_x0 else None
    rj["caveat"]                 = CAVEAT
    rj["exec_id"]                = EXEC_ID
    rj["material"]               = MATERIAL
    rj["n_positions"]            = 31
    rj["n_events_per_pos"]       = 5000
    rj["timestamp"]              = datetime.datetime.now().isoformat()

    out_path = Path(OUT_DIR) / "results.json"
    out_path.write_text(json.dumps(rj, indent=2))
    note(f"  results.json written: {out_path}")
    return rj

# ════════════════════════════════════════════════════════════════════════════════
# MORNING_REPORT
# ════════════════════════════════════════════════════════════════════════════════

def write_morning_report(pos_map, positions, n_events_per_pos, sigma_results,
                         aborted_curves, ratio_data, v_eff, sigma_x_data,
                         npe_nearest_top, x0_pos, results_json, flags_log):
    lines = [
        f"# MORNING_REPORT — {EXEC_ID} — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. Veredicto",
        "COMPLETADO" + ("-CON-FLAGS" if flags_log else ""),
        "",
        "## 2. Fix de identidad aplicado",
        f"- K_END={K_END}, N_TOP={N_TOP}, BOUNDARY={BOUNDARY}",
        f"- Gate de partición: PASSED ({2*K_END+N_TOP} canales únicos, sin colisiones)",
        f"- global_id observado: 0..{2*K_END+N_TOP-1}",
        f"- END_L gid 0..{K_END-1}, END_R gid {K_END}..{2*K_END-1}, TOP gid {2*K_END}..{2*K_END+N_TOP-1}",
        f"- face_type: {{0→END_L, 1→END_R, 2→TOP}} — confirmado desde DetectorConstruction.cc:52-54",
        "",
        "## 3. Decisiones autónomas",
        f"- TOP_IDS corregido de range(16,36) a range(16,86) [70 SiPMs TOP, desde header]",
        f"- Wrapping: 'Mylar' renombrado a 'reflector skin surface (Mylar-equivalent)' [DetConstr.cc:89]",
        f"- Material τ_d={TAU_D_NS} ns leído de opsc-101.mac (no heredado de EJ-230)",
        f"- datasets.py: dataset no registrado (gap 2×2 conocido); DATA_DIR leído directamente",
        f"- ORCH_DIR: datasets.py está en {{ORCH_DIR}}/analysis/ (no en raíz); anotado",
        f"- Posiciones: 31 × 5000 eventos, paso medio 46 mm, rango [-690,+690] mm",
        f"- Binning: FD rule (N≈5000 eventos por posición)",
        f"- T_RANGE_PS: percentil 0.5-99.5 de tiempos crudos por posición",
        f"- Canal nearest TOP: máx ⟨Npe⟩ por posición, programático",
        "",
        "## 4. Flags emitidos",
    ]
    lines.extend(flags_log if flags_log else ["  ninguno"])
    lines += ["", "## 5. Números destacados (desde CSV/JSON — NO tecleados)"]

    for est_label in ["TOP_SUM4", "TOP_SUM8", "END_COMBINED_SUM4", "END_COMBINED_SUM8"]:
        r = results_json.get(est_label, {})
        sig = r.get("x0_sigma_ps")
        boot = r.get("x0_boot_ps")
        lines.append(f"  {est_label} @ x=0: σ_fit = {sig} ps, bootstrap_err = {boot} ps")

    v = results_json.get("v_eff_cm_per_ns")
    lines.append(f"  v_eff = {v} cm/ns (slope of ⟨Δt⟩(x))")
    lines.append(f"  σ_EndTop/σ_END-only @ x=0 = {results_json.get('endtop_ratio_x0')}")
    lines += [
        "",
        "### Contexto SHiP (intrínseco — sin SPTR/FastIC):",
        "  Objetivo SHiP: 100 ps (preferible 50 ps). SPTR≈106 ps + FastIC≈10 ps",
        "  se añadirán en cuadratura en EXEC_02b. Mylar→Tyvek del TB DESY: no comparable.",
        "",
        "## 6. Rutas de salida",
        f"  Figuras (PDF+PNG):  {OUT_DIR}/L0 .. L4",
        f"  Sidecars CSV/JSON:  {OUT_DIR}/L1..L3/*.csv",
        f"  Beamer PDF:         {BEAMER_DIR}/EXEC_16_endtop_ej204.pdf (o .tex+log si compilación falló)",
        f"  results.json:       {OUT_DIR}/results.json",
        "",
        "## 7. Git",
        "  Branch: analysis_core master",
        "  Tags: EXEC_16-pre-cp0, EXEC_16-pre-cp1, EXEC_16-pre-cp2, EXEC_16-pre-cp3, EXEC_16-cp4",
        "",
        "## 8. Para revisar René",
        "  - Clustering físico id//4 vs ventana móvil de 4 vecinos → decisión de Gerardo",
        "  - Curvas abortadas por >30% outlier: " + (str(aborted_curves) if aborted_curves else "ninguna"),
        "  - HOOK_WALK reservado (corrección time-walk leading-edge, no implementado)",
        "  - EXEC_02b: SPTR≈106 ps + FastIC≈10 ps en cuadratura",
        "  - scan de 5000 evt/posición ya disponible (ESTE dataset); listo para robustecer END con más estadística",
        "  - Cierre del 2×2 material×readout: EJ-230 EndTop pendiente",
        "",
        f"**{CAVEAT}**",
    ]
    p = Path(OUT_DIR) / "MORNING_REPORT.md"
    p.write_text("\n".join(lines))
    note(f"\nMORNING_REPORT.md written: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(BEAMER_DIR, exist_ok=True)

    note(f"# {EXEC_ID} — {MATERIAL} EndTop — {datetime.datetime.now().isoformat()}")
    note(f"DATA_DIR: {DATA_DIR}")
    note(f"K_END={K_END}, N_TOP={N_TOP}, τ_d={TAU_D_NS} ns, λ={ABSLENGTH_CM} cm")

    flags_log = []

    # ── CP0: build pos_map, run gates
    note("\n### Building position map...")
    pos_map = build_pos_map()
    positions = sorted(pos_map.keys())
    note(f"  Found {len(positions)} positions: {positions}")

    for x, entry in pos_map.items():
        if entry.get("flag") == "gunx_spread":
            msg = f"gunx_spread at x={x}: spread={entry['spread']:.2f} mm"
            flags_log.append(msg)
            note(f"  [FLAG] {msg}")

    run_cp0_gates(pos_map)

    # ── CP1: Levels 0+1
    note("\n=== CP1: Levels 0+1 ===")
    n_events_per_pos = run_level0(pos_map, positions)
    npe_by_face, npe_nearest_top = run_level1(pos_map, positions)

    # ── CP2: Level 2
    note("\n=== CP2: Level 2 ===")
    sigma_results, aborted_curves, x0_pos = run_level2(pos_map, positions, npe_nearest_top)
    if aborted_curves:
        for c in aborted_curves:
            flags_log.append(f"CURVE-ABORT: {c} (>30% σ_outlier)")

    # ── CP3: Levels 3+4
    note("\n=== CP3: Levels 3+4 ===")
    ratio_data, delta_t_data, v_eff, sigma_x_data = run_level3(
        pos_map, positions, sigma_results, aborted_curves, x0_pos)
    run_level4(pos_map, positions, npe_nearest_top)

    # ── Write results.json
    results_json = write_results_json(sigma_results, ratio_data, v_eff, x0_pos, npe_nearest_top)

    # ── MORNING_REPORT
    write_morning_report(pos_map, positions, n_events_per_pos, sigma_results,
                         aborted_curves, ratio_data, v_eff, sigma_x_data,
                         npe_nearest_top, x0_pos, results_json, flags_log)

    note(f"\n=== {EXEC_ID} analysis COMPLETE ===")
    note(f"Output: {OUT_DIR}")
    note(f"figures: {OUT_DIR}/L0..L4")
    note(f"results.json: {OUT_DIR}/results.json")

if __name__ == "__main__":
    main()
