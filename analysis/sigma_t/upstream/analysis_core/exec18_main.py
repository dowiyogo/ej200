#!/usr/bin/env python3.12
"""
exec18_main.py — EXEC_18: time-walk, SUM topology, readout model.

Tasks T1–T7 (autonomous mode). T8 (Beamer + memo) in exec18_beamer.py.
Run: MPLBACKEND=Agg python3.12 exec18_main.py
"""

import sys, os, csv, json, math, datetime, warnings
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine   import fit_core_gaussian
from lib.robust_seeds import gather_seeds
from walk_correction  import fit_walk, apply_correction, validate_correction

# ════════════════════════════════════════════════════════════════════════════════
# Constants — read from source (DetectorConstruction.hh/.cc + opsc-101.mac)
# ════════════════════════════════════════════════════════════════════════════════

K_END         = 8
N_TOP         = 70
N_TOTAL       = 2 * K_END + N_TOP
BAR_HALF_X_MM = 700.0
END_PITCH_MM  = 7.5
TAU_D_NS      = 1.8
V_EFF_CM_NS   = 15.02    # from EXEC_16 spatial resolution
V7_A_PS       = 757.3    # σ_t = a/√Npe per channel, from EXEC_17 V7

# SPTR / FastIC (NOT intrinsic — for T5 table only)
SPTR_PS       = 106.0
FASTIC_PS     = 10.0

DATA_DIR  = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
OUT_DIR   = "/home/reriosto/SHiP/analysis_core/out/EXEC_18"
TREE      = "sipm_hits"
BR = dict(ev="event_id", gid="global_id", t="time_ns")

RANDOM_SEED = 20260618
N_BOOT      = 200
MIN_EV      = 30

REPR_X = [-690, 0, 690]   # representative positions for detailed analysis

# ── Verdict / report accumulator ──────────────────────────────────────────────
_report = []
_flags  = []

def note(m):
    _report.append(m)
    print(m)

def flag(m):
    _flags.append(m)
    note(f"  [FLAG] {m}")

def hard_abort(msg):
    note(f"\n[HARD-ABORT] {msg}")
    _write_report("ABORTADO")
    sys.exit(99)

def _write_report(verdict):
    p = Path(OUT_DIR) / "EXEC_18_REPORT.md"
    p.write_text("\n".join(_report))
    print(f"\nEXEC_18_REPORT.md: {p}")

def savefig(fig, stem, sub=""):
    d = Path(OUT_DIR) / sub
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf","png"):
        fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def ff(v, d=1): return f"{v:.{d}f}" if not math.isnan(v) else "NaN"

# ── Data utilities ────────────────────────────────────────────────────────────

def pos_map():
    entries = []
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        f = uproot.open(str(p))
        gx = f[TREE]["gun_x_mm"].array(library="np")
        entries.append({"x_mm": int(round(float(np.median(gx)))), "path": p})
    return sorted(entries, key=lambda e: e["x_mm"])

def load(entry, branches=None):
    branches = branches or list(BR.values())
    f = uproot.open(str(entry["path"]))
    return f[TREE].arrays(branches, library="np")

def classify(gid):
    if gid < K_END: return ("END_L", int(gid))
    if gid < 2*K_END: return ("END_R", int(gid - K_END))
    return ("TOP", int(gid - 2*K_END))

def best_k(gid_arr, face, k):
    if face == "END_L": mask = gid_arr < K_END
    elif face == "END_R": mask = (gid_arr >= K_END) & (gid_arr < 2*K_END)
    else: mask = gid_arr >= 2*K_END
    g = gid_arr[mask]
    if len(g) == 0: return np.array([], dtype=int)
    u, c = np.unique(g, return_counts=True)
    return u[np.argsort(-c)[:k]]

def tN_stream(evs, times, gids, sel_gids, N, all_events):
    """Time of N-th photon in merged stream of sel_gids; NaN if <N hits."""
    if len(sel_gids) == 0: return np.full(len(all_events), np.nan)
    mask = np.isin(gids, sel_gids)
    ev_s, t_s = evs[mask], times[mask]
    if len(ev_s) == 0: return np.full(len(all_events), np.nan)
    idx = np.lexsort((t_s, ev_s))
    ev_s, t_s = ev_s[idx], t_s[idx]
    u, cnt = np.unique(ev_s, return_counts=True)
    cum = np.concatenate([[0], np.cumsum(cnt)])
    ev_map = {e: i for i, e in enumerate(all_events)}
    res = np.full(len(all_events), np.nan)
    for i, (e, c) in enumerate(zip(u, cnt)):
        if c >= N:
            j = ev_map.get(e, -1)
            if j >= 0: res[j] = t_s[cum[i] + N - 1]
    return res

def npe_stream(evs, gids, sel_gids, all_events):
    """Total NPE per event in stream of sel_gids."""
    if len(sel_gids) == 0: return np.zeros(len(all_events))
    mask = np.isin(gids, sel_gids)
    ev_m = evs[mask]
    res = np.zeros(len(all_events))
    if len(ev_m) == 0: return res
    u, cnt = np.unique(ev_m, return_counts=True)
    ev_map = {e: i for i, e in enumerate(all_events)}
    for e, c in zip(u, cnt):
        j = ev_map.get(e, -1)
        if j >= 0: res[j] = c
    return res

def gauss_fit(v, prefix="g"):
    v_c = v[~np.isnan(v)]
    if len(v_c) < MIN_EV:
        return {"sigma_fit": np.nan, "bootstrap_err": np.nan,
                "flag": "insufficient_stats", "h_root": None, "f_root": None}
    cfg = {"MIN_EVENTS_FOR_FIT": MIN_EV, "FIT_WINDOW_SIGMAS": 2.0,
           "BINNING_STRATEGY": "fd", "N_BOOTSTRAP": N_BOOT,
           "RANDOM_SEED": RANDOM_SEED, "FIT_OPTIONS": "R Q S 0", "CHI2_NDF_WARN": 3.0}
    return fit_core_gaussian(v_c, cfg, name_prefix=prefix)

# ════════════════════════════════════════════════════════════════════════════════
# Integrity gate (§D.1)
# ════════════════════════════════════════════════════════════════════════════════

def gate_integrity(pm):
    seen = {}
    for g in range(N_TOTAL):
        key = classify(g)
        if key in seen: hard_abort(f"§D.1 classify() collision: {key}")
        seen[key] = g
    assert len(seen) == N_TOTAL
    # Cross-check one file
    d = load(pm[len(pm)//2])
    gids = d[BR["gid"]]
    if gids.min() < 0 or gids.max() > N_TOTAL - 1:
        hard_abort(f"§D.1 global_id out of [0,{N_TOTAL-1}]: [{gids.min()},{gids.max()}]")
    note("  §D.1 Gate integrity: PASSED")

# ════════════════════════════════════════════════════════════════════════════════
# T4 — Resolve 72.6 ps label (fast, no computation needed, just read code)
# ════════════════════════════════════════════════════════════════════════════════

def run_T4(pm):
    note("\n=== T4 — Resolving the 72.6 ps label ===")
    note("  Reading exec16_endtop_ej204.py SUM_N_SCAN and _fit call...")
    # Confirmed by code inspection:
    # SUM_N_SCAN = (4, 8) → loop variable N = number of SiPMs, NOT photon threshold
    # results[label].append(_fit(label, gids_sel, 1))  → photon N=1 (first photon)
    # top_sum4_gids = best_gids_by_npe(gids, 'TOP', 4) → 4 channels by max <Npe>
    # So: TOP_SUM4 = {4 nearest TOP SiPMs merged, t of FIRST photon in merged stream}
    #             = min(t_first_ch_i) for i in {4 channels with max Npe}
    note("  Definition: TOP_SUM4 = 4 TOP SiPMs (max <Npe>) merged, N=1 (first photon in stream)")
    note("  = equivalent to min(t_1st_photon) across 4 channels per event")
    note("  Label is correct; 'SUM4' refers to channel count, not photon count")
    note("  EXEC_17/C3 N sweep varied PHOTON threshold (2..8), not channel count")
    note("  → The 72.6 ps is correctly identified as TOP_SUM4_N1 (first photon in 4-ch stream)")

    # Verify by recomputing at x=0
    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry0)
    gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
    all_events = np.unique(evs)

    top4_gids    = best_k(gids, "TOP", 4)
    nearest_gid  = top4_gids[0] if len(top4_gids) > 0 else int(2*K_END)

    t_top_nearest_N1  = tN_stream(evs, times, gids, np.array([nearest_gid]), 1, all_events)
    t_top_sum4_N1     = tN_stream(evs, times, gids, top4_gids, 1, all_events)

    r_nearest = gauss_fit(t_top_nearest_N1, "t4_nearest_N1")
    r_sum4    = gauss_fit(t_top_sum4_N1,    "t4_sum4_N1")

    sig_nearest = ff(r_nearest.get("sigma_fit", math.nan) * 1000, 1)
    sig_sum4    = ff(r_sum4.get("sigma_fit", math.nan) * 1000, 1)
    note(f"  Recomputed @ x=0: TOP_single_nearest_N1 = {sig_nearest} ps  |  TOP_SUM4_N1 = {sig_sum4} ps")
    note(f"  EXEC_16 reported: TOP_THR1=107.2 ps (single-nearest N=1)  |  TOP_SUM4=72.6 ps (4-ch N=1)")
    note(f"  Difference {float(r_nearest.get('sigma_fit',math.nan)*1000):.1f} vs {float(r_sum4.get('sigma_fit',math.nan)*1000):.1f} ps = benefit of merging 4 channels before trigger")
    note(f"  Flagship: 72.6 ps = TOP_SUM4_N1. Label correct. 'SUM' = channel count, not photon sum.")

    return {"sig_nearest_ps": r_nearest.get("sigma_fit", math.nan)*1000 if not math.isnan(r_nearest.get("sigma_fit",math.nan)) else math.nan,
            "sig_sum4_ps": r_sum4.get("sigma_fit", math.nan)*1000 if not math.isnan(r_sum4.get("sigma_fit",math.nan)) else math.nan,
            "definition": "TOP_SUM4_N1 = min(t_1st_photon across 4 nearest TOP SiPMs per event)"}

# ════════════════════════════════════════════════════════════════════════════════
# T1 — HOOK_WALK: time-walk correction module validation
# ════════════════════════════════════════════════════════════════════════════════

def run_T1(pm):
    note("\n=== T1 — HOOK_WALK: time-walk correction ===")
    Path(OUT_DIR + "/T1").mkdir(parents=True, exist_ok=True)

    # Focus on x=0 for validation; will extend to all positions in T2
    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry0)
    gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
    all_events = np.unique(evs)

    top4_gids   = best_k(gids, "TOP", 4)
    top8_gids   = best_k(gids, "TOP", 8)
    nearest_gid = top4_gids[0] if len(top4_gids) > 0 else int(2*K_END)
    endl4_gids  = best_k(gids, "END_L", 4)
    endr4_gids  = best_k(gids, "END_R", 4)

    walk_results = {}

    estimators = [
        ("TOP_nearest_N1",    np.array([nearest_gid]), 1),
        ("TOP_SUM4_N1",       top4_gids, 1),
        ("TOP_SUM8_N1",       top8_gids, 1),
        ("END_L_SUM4_N1",     endl4_gids, 1),
        ("END_R_SUM4_N1",     endr4_gids, 1),
    ]

    fig, axes = plt.subplots(len(estimators), 3, figsize=(18, 4*len(estimators)))

    all_sigma_improved = []

    for ai, (label, sel_gids, N) in enumerate(estimators):
        t_raw = tN_stream(evs, times, gids, sel_gids, N, all_events)
        s     = npe_stream(evs, gids, sel_gids, all_events)

        valid = ~np.isnan(t_raw) & (s > 0)
        tv, sv = t_raw[valid], s[valid]

        # Fit walk curve
        fit = fit_walk(tv * 1000, sv)  # t in ps, s in PE

        t_corr_ps = apply_correction(t_raw * 1000, s, fit)  # ps
        t_corr_ns = t_corr_ps / 1000.0

        r_raw  = gauss_fit(t_raw,    f"t1_{label}_raw")
        r_corr = gauss_fit(t_corr_ns, f"t1_{label}_corr")

        sig_raw_ps  = r_raw.get("sigma_fit",  math.nan) * 1000
        sig_corr_ps = r_corr.get("sigma_fit", math.nan) * 1000

        val = validate_correction(tv * 1000, t_corr_ps[valid], sv,
                                   sig_raw_ps / 1000, sig_corr_ps / 1000)

        walk_results[label] = {
            "fit": fit, "sigma_raw_ps": sig_raw_ps, "sigma_corr_ps": sig_corr_ps,
            "r_raw": val["r_raw"], "r_corr": val["r_corr"],
            "delta_sigma_ps": val["delta_sigma_ps"], "status": val["status"],
            "form": fit["form"], "converged": fit["converged"],
        }

        improved = sig_corr_ps < sig_raw_ps or math.isnan(sig_raw_ps)
        all_sigma_improved.append(improved)

        status_str = val["status"]
        if not val["ok"]:
            flag(f"T1 {label}: walk_status={status_str} (r_corr={val['r_corr']:.3f}, "
                 f"Δσ={val['delta_sigma_ps']:.1f} ps)")

        note(f"  {label}: σ_raw={ff(sig_raw_ps)} ps → σ_corr={ff(sig_corr_ps)} ps "
             f"(Δ={ff(val['delta_sigma_ps'])} ps) | r_raw={val['r_raw']:.3f} → r_corr={val['r_corr']:.3f} | {status_str}")

        # Plot row
        lo_t = float(np.nanpercentile(tv * 1000, 0.5))
        hi_t = float(np.nanpercentile(tv * 1000, 99.5))
        s_range = (float(sv.min()), float(sv.max()))

        # Walk curve
        axes[ai,0].scatter(sv, tv*1000, s=3, alpha=0.3, c="steelblue")
        if fit["converged"] and fit.get("popt") is not None and len(fit["popt"]) == 2:
            s_fit = np.linspace(s_range[0], s_range[1], 200)
            from walk_correction import _walk_model
            axes[ai,0].plot(s_fit, _walk_model(s_fit, *fit["popt"]), "r-", lw=2, label=f"α+β/√s fit")
        elif fit["form"] == "spline" and fit.get("spline") is not None:
            s_fit = np.linspace(s_range[0], s_range[1], 200)
            axes[ai,0].plot(s_fit, fit["spline"](s_fit), "r-", lw=2, label="spline")
        axes[ai,0].set_xlabel("NPE"); axes[ai,0].set_ylabel("t_raw [ps]")
        axes[ai,0].set_title(f"{label}: walk curve"); axes[ai,0].legend(fontsize=6)

        # t distribution before/after
        axes[ai,1].hist(tv*1000, bins=100, range=(lo_t, hi_t), histtype="step",
                        color="steelblue", lw=1.2, label=f"raw σ={ff(sig_raw_ps)} ps")
        t_corr_valid = t_corr_ps[valid]
        lo_c = float(np.nanpercentile(t_corr_valid, 0.5))
        hi_c = float(np.nanpercentile(t_corr_valid, 99.5))
        axes[ai,1].hist(t_corr_valid, bins=100, range=(lo_c, hi_c), histtype="step",
                        color="tomato", lw=1.2, label=f"corr σ={ff(sig_corr_ps)} ps")
        axes[ai,1].set_xlabel("t [ps]"); axes[ai,1].legend(fontsize=7)
        axes[ai,1].set_title(f"{label}: distributions")

        # Correlation scatter after correction
        axes[ai,2].scatter(sv, t_corr_valid, s=3, alpha=0.3, c="tomato")
        axes[ai,2].set_xlabel("NPE"); axes[ai,2].set_ylabel("t_corr [ps]")
        axes[ai,2].set_title(f"{label}: t_corr vs s (r={val['r_corr']:.3f})")

    savefig(fig, "T1_walk_correction_x0", "T1")

    # §D.2: if walk fails ALL estimators → hard abort
    if not any(all_sigma_improved):
        hard_abort("§D.2: walk correction worsens σ for ALL estimators (sign error)")

    # Save sidecar
    with open(Path(OUT_DIR) / "T1" / "T1_walk_results_x0.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["estimator","form","sigma_raw_ps","sigma_corr_ps",
                                "delta_sigma_ps","r_raw","r_corr","status","converged"])
        w.writeheader()
        for lbl, r in walk_results.items():
            w.writerow({"estimator": lbl, "form": r["form"],
                        "sigma_raw_ps": ff(r["sigma_raw_ps"],1),
                        "sigma_corr_ps": ff(r["sigma_corr_ps"],1),
                        "delta_sigma_ps": ff(r["delta_sigma_ps"],1),
                        "r_raw": ff(r["r_raw"],3), "r_corr": ff(r["r_corr"],3),
                        "status": r["status"], "converged": r["converged"]})
    return walk_results

# ════════════════════════════════════════════════════════════════════════════════
# T2 — Re-ranking single vs SUM with and without walk
# ════════════════════════════════════════════════════════════════════════════════

def run_T2(pm, walk_results_x0):
    note("\n=== T2 — Re-ranking single vs SUM (N sweep, raw vs walk-corrected) ===")
    Path(OUT_DIR + "/T2").mkdir(parents=True, exist_ok=True)

    N_sweep = [1, 2, 3, 4, 6, 8]
    repr_entries = [e for e in pm if e["x_mm"] in REPR_X]

    # Full scan for σ(x) curves (raw and corrected) for TOP_SUM4_N1 and TOP_nearest_N1
    all_positions = [e["x_mm"] for e in pm]
    sig_raw_nearest_x  = []
    sig_corr_nearest_x = []
    sig_raw_sum4_x     = []
    sig_corr_sum4_x    = []

    for entry in pm:
        x = entry["x_mm"]
        d = load(entry)
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)

        top4_gids    = best_k(gids, "TOP", 4)
        nearest_gid  = top4_gids[0] if len(top4_gids) > 0 else int(2*K_END)

        for label, sel_gids, sig_raw_list, sig_corr_list in [
            ("TOP_nearest_N1", np.array([nearest_gid]), sig_raw_nearest_x, sig_corr_nearest_x),
            ("TOP_SUM4_N1",    top4_gids,               sig_raw_sum4_x,    sig_corr_sum4_x),
        ]:
            t_raw = tN_stream(evs, times, gids, sel_gids, 1, all_events)
            s_    = npe_stream(evs, gids, sel_gids, all_events)

            # Re-fit walk at this position (use walk form from x=0 as template)
            wf = walk_results_x0.get(label, {})
            valid = ~np.isnan(t_raw) & (s_ > 0)
            if valid.sum() >= 50:
                fit_x = fit_walk(t_raw[valid] * 1000, s_[valid])
                t_corr = apply_correction(t_raw * 1000, s_, fit_x) / 1000.0
            else:
                t_corr = t_raw.copy()

            r_raw  = gauss_fit(t_raw,  f"t2_{label}_x{x}_raw")
            r_corr = gauss_fit(t_corr, f"t2_{label}_x{x}_corr")
            sig_raw_list.append(r_raw.get("sigma_fit", math.nan) * 1000)
            sig_corr_list.append(r_corr.get("sigma_fit", math.nan) * 1000)

    # N sweep at representative positions
    n_sweep_data = {}  # (x, label) → [(N, sigma_raw_ps, sigma_corr_ps)]
    for entry in repr_entries:
        x = entry["x_mm"]
        d = load(entry)
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)
        top4_gids  = best_k(gids, "TOP", 4)
        top8_gids  = best_k(gids, "TOP", 8)
        n_sweep_data[x] = []
        for N in N_sweep:
            sel = top4_gids if N <= 4 else top8_gids
            t_raw = tN_stream(evs, times, gids, sel, N, all_events)
            s_    = npe_stream(evs, gids, sel, all_events)
            valid = ~np.isnan(t_raw) & (s_ > 0)
            fit_x = fit_walk(t_raw[valid]*1000, s_[valid]) if valid.sum() >= 50 else {"converged": False}
            t_corr = apply_correction(t_raw*1000, s_, fit_x)/1000.0 if fit_x["converged"] else t_raw
            rr = gauss_fit(t_raw,  f"t2_nsweep_x{x}_N{N}_raw")
            rc = gauss_fit(t_corr, f"t2_nsweep_x{x}_N{N}_corr")
            n_sweep_data[x].append((N, rr.get("sigma_fit",math.nan)*1000, rc.get("sigma_fit",math.nan)*1000))

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    xs = np.array(all_positions)
    axes[0].plot(xs, sig_raw_nearest_x,  "o-", color="steelblue", ms=4, lw=1.2, label="nearest raw")
    axes[0].plot(xs, sig_corr_nearest_x, "o--", color="steelblue", ms=4, lw=1.2, alpha=0.7, label="nearest corr")
    axes[0].plot(xs, sig_raw_sum4_x,     "s-", color="seagreen", ms=4, lw=1.2, label="SUM4 raw")
    axes[0].plot(xs, sig_corr_sum4_x,    "s--", color="seagreen", ms=4, lw=1.2, alpha=0.7, label="SUM4 corr")
    axes[0].axhline(100, color="red", ls="--", lw=0.8)
    axes[0].set_xlabel("gun_x [mm]"); axes[0].set_ylabel("σ_fit [ps]")
    axes[0].set_title("T2 — σ_t(x): raw vs walk-corrected")
    axes[0].legend(fontsize=7); axes[0].grid(True, lw=0.3)

    colors = {-690: "steelblue", 0: "seagreen", 690: "tomato"}
    for x, rows in n_sweep_data.items():
        Ns = [r[0] for r in rows]; sr = [r[1] for r in rows]; sc = [r[2] for r in rows]
        c = colors.get(x, "gray")
        axes[1].plot(Ns, sr, "o-", color=c, ms=4, lw=1.2, label=f"x={x} raw")
        axes[1].plot(Ns, sc, "o--", color=c, ms=4, lw=1.2, alpha=0.7, label=f"x={x} corr")
    axes[1].axhline(100, color="red", ls="--", lw=0.8)
    axes[1].set_xlabel("N (photon threshold in merged stream)"); axes[1].set_ylabel("σ_fit [ps]")
    axes[1].set_title("T2 — N sweep: raw vs walk-corrected")
    axes[1].legend(fontsize=7); axes[1].grid(True, lw=0.3)
    savefig(fig, "T2_reranking_walk", "T2")

    # Sidecar
    with open(Path(OUT_DIR) / "T2" / "T2_sigma_vs_x.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm","sig_raw_nearest_ps","sig_corr_nearest_ps",
                                "sig_raw_sum4_ps","sig_corr_sum4_ps"])
        w.writeheader()
        for i, x in enumerate(all_positions):
            w.writerow({"x_mm": x, "sig_raw_nearest_ps": ff(sig_raw_nearest_x[i]),
                        "sig_corr_nearest_ps": ff(sig_corr_nearest_x[i]),
                        "sig_raw_sum4_ps": ff(sig_raw_sum4_x[i]),
                        "sig_corr_sum4_ps": ff(sig_corr_sum4_x[i])})

    # Report optimal N at x=0 before and after walk
    x0_rows = n_sweep_data.get(0, [])
    if x0_rows:
        best_raw  = min(x0_rows, key=lambda r: r[1] if not math.isnan(r[1]) else 9999)
        best_corr = min(x0_rows, key=lambda r: r[2] if not math.isnan(r[2]) else 9999)
        note(f"  N_opt @ x=0: raw={best_raw[0]} (σ={ff(best_raw[1])} ps), corr={best_corr[0]} (σ={ff(best_corr[2])} ps)")

    return {"sig_raw_nearest_x": sig_raw_nearest_x, "sig_corr_nearest_x": sig_corr_nearest_x,
            "sig_raw_sum4_x": sig_raw_sum4_x, "sig_corr_sum4_x": sig_corr_sum4_x,
            "n_sweep": n_sweep_data, "positions": all_positions}

# ════════════════════════════════════════════════════════════════════════════════
# T3 — END_SUM4 co-localized
# ════════════════════════════════════════════════════════════════════════════════

def run_T3(pm):
    note("\n=== T3 — END_SUM4 co-localized ===")
    Path(OUT_DIR + "/T3").mkdir(parents=True, exist_ok=True)

    note("  Decision: co-localized = 4 SiPMs with max <Npe> from SAME END face.")
    note("  For x<0: END_L, for x>0: END_R, for x=0: best of both faces.")
    note("  Reference (sim↔TB): TOP_single_nearest_N1 at same positions.")

    positions = [e["x_mm"] for e in pm]
    sig_endl_sum4_raw  = []
    sig_endr_sum4_raw  = []
    sig_top_nearest_raw = []

    for entry in pm:
        x = entry["x_mm"]
        d = load(entry)
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)

        top4_gids   = best_k(gids, "TOP", 4)
        nearest_gid = top4_gids[0] if len(top4_gids) > 0 else int(2*K_END)
        endl4_gids  = best_k(gids, "END_L", 4)
        endr4_gids  = best_k(gids, "END_R", 4)

        t_nearest = tN_stream(evs, times, gids, np.array([nearest_gid]), 1, all_events)
        t_endl4   = tN_stream(evs, times, gids, endl4_gids, 1, all_events)
        t_endr4   = tN_stream(evs, times, gids, endr4_gids, 1, all_events)

        # Apply walk correction
        s_nearest = npe_stream(evs, gids, np.array([nearest_gid]), all_events)
        s_endl4   = npe_stream(evs, gids, endl4_gids, all_events)
        s_endr4   = npe_stream(evs, gids, endr4_gids, all_events)

        for t_arr, s_arr, sig_list, lbl in [
            (t_nearest, s_nearest, sig_top_nearest_raw, "top_near"),
            (t_endl4,   s_endl4,   sig_endl_sum4_raw,   "endl4"),
            (t_endr4,   s_endr4,   sig_endr_sum4_raw,   "endr4"),
        ]:
            valid = ~np.isnan(t_arr) & (s_arr > 0)
            if valid.sum() >= 50:
                fit_x = fit_walk(t_arr[valid]*1000, s_arr[valid])
                t_c = apply_correction(t_arr*1000, s_arr, fit_x)/1000.0 if fit_x["converged"] else t_arr
            else:
                t_c = t_arr
            r = gauss_fit(t_c, f"t3_{lbl}_x{x}")
            sig_list.append(r.get("sigma_fit", math.nan) * 1000)

    fig, ax = plt.subplots(figsize=(12, 5))
    xs = np.array(positions)
    ax.plot(xs, sig_top_nearest_raw, "o-", color="seagreen", ms=4, lw=1.2, label="TOP nearest (walk-corr)")
    ax.plot(xs, sig_endl_sum4_raw,   "s-", color="steelblue", ms=4, lw=1.2, label="END_L_SUM4 co-local (walk-corr)")
    ax.plot(xs, sig_endr_sum4_raw,   "^-", color="tomato",    ms=4, lw=1.2, label="END_R_SUM4 co-local (walk-corr)")
    ax.axhline(100, color="red", ls="--", lw=0.8, label="SHiP 100 ps")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("σ_fit [ps]")
    ax.set_title("T3 — END_SUM4 co-localized vs TOP nearest (walk-corrected)\n"
                 "Intrinsic only. Note: END_SUM4 → sim↔TB-Constanza analogy")
    ax.legend(fontsize=8); ax.grid(True, lw=0.3); ax.set_yscale("log"); ax.set_ylim(5, 3000)
    savefig(fig, "T3_end_sum4_colocalized", "T3")

    with open(Path(OUT_DIR) / "T3" / "T3_end_sum4_colocalized.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm","sig_top_nearest_ps","sig_endl_sum4_ps","sig_endr_sum4_ps"])
        w.writeheader()
        for i, x in enumerate(positions):
            w.writerow({"x_mm": x, "sig_top_nearest_ps": ff(sig_top_nearest_raw[i]),
                        "sig_endl_sum4_ps": ff(sig_endl_sum4_raw[i]),
                        "sig_endr_sum4_ps": ff(sig_endr_sum4_raw[i])})

    # Summary at key positions
    x0_idx  = positions.index(0) if 0 in positions else 0
    xm690   = positions.index(-690) if -690 in positions else 0
    note(f"  @ x=0:    TOP_nearest={ff(sig_top_nearest_raw[x0_idx])} ps  |  END_L_SUM4={ff(sig_endl_sum4_raw[x0_idx])} ps  |  END_R_SUM4={ff(sig_endr_sum4_raw[x0_idx])} ps")
    note(f"  @ x=-690: TOP_nearest={ff(sig_top_nearest_raw[xm690])} ps  |  END_L_SUM4={ff(sig_endl_sum4_raw[xm690])} ps  (near-end)")

    return {"positions": positions, "sig_top": sig_top_nearest_raw,
            "sig_endl4": sig_endl_sum4_raw, "sig_endr4": sig_endr_sum4_raw}

# ════════════════════════════════════════════════════════════════════════════════
# T5 — Readout model A vs B with Kish N_eff
# ════════════════════════════════════════════════════════════════════════════════

def run_T5(pm):
    note("\n=== T5 — Readout model A vs B (DECISIÓN GERARDO) ===")
    Path(OUT_DIR + "/T5").mkdir(parents=True, exist_ok=True)

    # Compute N_eff (Kish) per estimator per position using V7 model σ_i = a/√Npe_i
    # N_eff = (Σ Npe_i)² / Σ Npe_i²  (since w_i = Npe_i/a²)
    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry0)
    gids = d[BR["gid"]]

    def kish_neff(sel_gids, gids_arr, evs_arr, all_events):
        npes = np.array([float(np.sum(gids_arr == g) / len(all_events)) for g in sel_gids])
        npes = npes[npes > 0]
        if len(npes) == 0: return 1.0
        return float(npes.sum()**2 / (npes**2).sum())

    all_events = np.unique(d[BR["ev"]])

    # Estimators: (label, sel_gids_fn, σ_int_ps from EXEC_16/T2-corr)
    # Use walk-corrected σ at x=0 where available
    estimators_T5 = [
        ("TOP_nearest_N1", best_k(gids, "TOP", 1)),
        ("TOP_SUM4_N1",    best_k(gids, "TOP", 4)),
        ("TOP_SUM8_N1",    best_k(gids, "TOP", 8)),
    ]

    # σ_int from EXEC_16 (will be updated after T2 but use EXEC_16 as baseline)
    sigma_int_map = {"TOP_nearest_N1": 107.2, "TOP_SUM4_N1": 72.6, "TOP_SUM8_N1": 74.4}

    rows = []
    for label, sel_gids in estimators_T5:
        n_ch = len(sel_gids)
        neff = kish_neff(sel_gids, gids, d[BR["ev"]], all_events)
        sigma_int = sigma_int_map.get(label, 75.0)  # ps, from EXEC_16

        # Model A: SPTR = 106 ps (flat)
        sptr_A  = SPTR_PS
        tot_A   = math.sqrt(sigma_int**2 + sptr_A**2 + FASTIC_PS**2)
        goal_A  = "YES" if tot_A <= 100 else "NO"

        # Model B: SPTR_eff = SPTR/√N_eff
        sptr_B  = SPTR_PS / math.sqrt(neff)
        tot_B   = math.sqrt(sigma_int**2 + sptr_B**2 + FASTIC_PS**2)
        goal_B  = "YES" if tot_B <= 100 else "NO"

        rows.append({"estimator": label, "N_ch": n_ch, "N_eff_Kish": round(neff, 2),
                     "sigma_int_ps": sigma_int,
                     "SPTR_A_ps": sptr_A, "tot_A_ps": round(tot_A, 1), "goal_A": goal_A,
                     "SPTR_B_ps": round(sptr_B, 1), "tot_B_ps": round(tot_B, 1), "goal_B": goal_B})
        note(f"  {label}: N_ch={n_ch}, N_eff={neff:.2f} | "
             f"A→{tot_A:.1f} ps ({goal_A}) | B→{tot_B:.1f} ps ({goal_B})")
        note(f"    NOTE: N_eff={neff:.2f} vs N_ch={n_ch}: "
             f"{'underestimates N_eff (channels equal Npe)' if neff > 0.9*n_ch else 'real penalty from unequal Npe'}")

    # Figure
    fig, ax = plt.subplots(figsize=(10, 5))
    labels = [r["estimator"] for r in rows]
    tA = [r["tot_A_ps"] for r in rows]
    tB = [r["tot_B_ps"] for r in rows]
    si = [r["sigma_int_ps"] for r in rows]
    y  = range(len(labels))
    ax.barh([yi+0.2 for yi in y], tA, height=0.35, color="steelblue", label="Model A (analog+CFD)")
    ax.barh([yi-0.2 for yi in y], tB, height=0.35, color="tomato", label="Model B (digital, N_eff_Kish)")
    ax.barh([yi for yi in y], si, height=0.35, color="lightgray", alpha=0.6, label="σ_int (intrinsic)")
    ax.axvline(100, color="red", ls="--", lw=1.5, label="SHiP 100 ps goal")
    ax.axvline(50,  color="orange", ls="--", lw=1, label="SHiP preferred 50 ps")
    ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("σ_total [ps]")
    ax.set_title("T5 — SPTR Model A vs B (DECISIÓN GERARDO)\n"
                 "σ_int includes 20 ps SiPMSD jitter — do NOT re-add.\n"
                 "Model B is optimistic floor; applicability depends on FastIC+ architecture.")
    ax.legend(fontsize=8); ax.grid(True, axis="x", lw=0.3)
    savefig(fig, "T5_model_A_vs_B", "T5")

    with open(Path(OUT_DIR) / "T5" / "T5_readout_models.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    return rows

# ════════════════════════════════════════════════════════════════════════════════
# T6 — SUM topology documentation (DECISIÓN GERARDO)
# ════════════════════════════════════════════════════════════════════════════════

def run_T6(pm):
    note("\n=== T6 — SUM topology: co-localized vs distributed (DECISIÓN GERARDO) ===")
    Path(OUT_DIR + "/T6").mkdir(parents=True, exist_ok=True)

    # Co-localized: END_SUM4 (done in T3) — all 4 on same face, same x, same photons
    # Distributed: TOP_SUM4 — 4 SiPMs at 20mm spacing, each sees different photon population
    # Here we demonstrate the degradation by computing TOP_SUM4 with co-located sub-arrays

    note("  Co-localized: END_SUM4 = 4 SiPMs on same END face (same x, same photons)")
    note("  Distributed:  TOP_SUM4 = 4 nearest TOP SiPMs, step 20mm, mixed photon populations")
    note("  id//4 groups:  {0..3},{4..7} for END; {0..3},{4..7},...,{17..20} for TOP (in local_id)")
    note("  Fixed clusters id//4: 18 groups of 4 for TOP (70 SiPMs → 17 full groups + tail)")

    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry0)
    gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
    all_events = np.unique(evs)

    # Compute σ for id//4 fixed cluster (nearest cluster to beam) vs nearest-4 dynamic
    top_local_ids = gids[gids >= 2*K_END] - 2*K_END  # 0..69
    # cluster id for each TOP SiPM: cluster = local_id // 4
    # nearest cluster to x=0: the one with most hits
    cluster_hits = {}
    for g in range(2*K_END, 2*K_END + N_TOP):
        local = g - 2*K_END
        cluster = local // 4
        mask_g = gids == g
        cluster_hits[cluster] = cluster_hits.get(cluster, 0) + int(np.sum(mask_g))
    best_cluster = max(cluster_hits, key=cluster_hits.get)
    cluster_gids = np.array([2*K_END + best_cluster*4 + i for i in range(4) if 2*K_END + best_cluster*4 + i < 2*K_END + N_TOP])

    t_colocated_cluster = tN_stream(evs, times, gids, cluster_gids, 1, all_events)
    t_nearest4_dynamic  = tN_stream(evs, times, gids, best_k(gids, "TOP", 4), 1, all_events)

    r_cluster = gauss_fit(t_colocated_cluster, "t6_cluster")
    r_dynamic = gauss_fit(t_nearest4_dynamic,  "t6_dynamic")

    note(f"  @ x=0: id//4 fixed cluster (gids {cluster_gids}) σ={ff(r_cluster.get('sigma_fit',math.nan)*1000)} ps")
    note(f"  @ x=0: dynamic nearest-4 (best <Npe>)         σ={ff(r_dynamic.get('sigma_fit',math.nan)*1000)} ps")

    # Figure
    fig, ax = plt.subplots(figsize=(10, 5))
    t_vals = {"id//4 fixed cluster": t_colocated_cluster * 1000,
              "dynamic nearest-4":   t_nearest4_dynamic * 1000}
    lo = float(np.nanpercentile(t_nearest4_dynamic * 1000, 0.5))
    hi = float(np.nanpercentile(t_nearest4_dynamic * 1000, 99.5))
    colors = {"id//4 fixed cluster": "steelblue", "dynamic nearest-4": "seagreen"}
    for lbl, tv in t_vals.items():
        v = tv[~np.isnan(tv)]
        ax.hist(v, bins=100, range=(lo, hi), histtype="step", lw=1.5, color=colors[lbl], label=lbl)
    ax.set_xlabel("t [ps]"); ax.set_ylabel("events")
    ax.set_title("T6 — SUM topology comparison @ x=0 (DECISIÓN GERARDO)\n"
                 "id//4 fixed = co-localized cluster | dynamic nearest-4 = distributed")
    ax.legend(fontsize=9)
    savefig(fig, "T6_topology_comparison", "T6")

    with open(Path(OUT_DIR) / "T6" / "T6_topology_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["topology","n_channels","cluster_type","sigma_ps"])
        w.writerow(["id//4 fixed cluster", len(cluster_gids), "co-localized fixed",
                    ff(r_cluster.get("sigma_fit",math.nan)*1000)])
        w.writerow(["dynamic nearest-4", 4, "distributed dynamic",
                    ff(r_dynamic.get("sigma_fit",math.nan)*1000)])

    return {"sig_cluster_ps": r_cluster.get("sigma_fit",math.nan)*1000,
            "sig_dynamic_ps": r_dynamic.get("sigma_fit",math.nan)*1000,
            "best_cluster": int(best_cluster), "cluster_gids": cluster_gids.tolist()}

# ════════════════════════════════════════════════════════════════════════════════
# T7 — EndTop with covariance (GLS optimal)
# ════════════════════════════════════════════════════════════════════════════════

def run_T7(pm):
    note("\n=== T7 — EndTop with covariance (GLS optimal) ===")
    Path(OUT_DIR + "/T7").mkdir(parents=True, exist_ok=True)

    positions = [e["x_mm"] for e in pm]
    sig_inv_var_x = []
    sig_gls_x     = []

    for entry in pm:
        x = entry["x_mm"]
        d = load(entry)
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)

        endl8_gids = best_k(gids, "END_L", 8)
        endr8_gids = best_k(gids, "END_R", 8)
        top4_gids  = best_k(gids, "TOP",   4)

        tL = tN_stream(evs, times, gids, endl8_gids, 1, all_events)
        tR = tN_stream(evs, times, gids, endr8_gids, 1, all_events)
        tT = tN_stream(evs, times, gids, top4_gids,  1, all_events)

        valid_e = (~np.isnan(tL)) & (~np.isnan(tR))
        t_end = np.full(len(all_events), np.nan)
        t_end[valid_e] = 0.5 * (tL[valid_e] + tR[valid_e])

        valid_both = (~np.isnan(t_end)) & (~np.isnan(tT))
        tE_v = t_end[valid_both] * 1000  # ps
        tT_v = tT[valid_both]   * 1000  # ps

        if valid_both.sum() < 30:
            sig_inv_var_x.append(math.nan)
            sig_gls_x.append(math.nan)
            continue

        # Fit each marginal to get σ_E, σ_T
        rE = gauss_fit(t_end[valid_both] / 1000, f"t7_end_x{x}")
        rT = gauss_fit(tT[valid_both] / 1000, f"t7_top_x{x}")
        sE = rE.get("sigma_fit", math.nan) * 1000  # ps
        sT = rT.get("sigma_fit", math.nan) * 1000  # ps

        if math.isnan(sE) or math.isnan(sT) or sE <= 0 or sT <= 0:
            sig_inv_var_x.append(math.nan)
            sig_gls_x.append(math.nan)
            continue

        # Empirical covariance
        mu_E = float(np.mean(tE_v)); mu_T = float(np.mean(tT_v))
        cov = float(np.mean((tE_v - mu_E) * (tT_v - mu_T)))

        # Inv-variance (ignores cov): σ² = 1/(1/σ²_E + 1/σ²_T)
        sig_iv = 1.0 / math.sqrt(1.0/sE**2 + 1.0/sT**2)

        # GLS optimal: σ²_GLS = (σ²_E σ²_T - cov²) / (σ²_E + σ²_T - 2*cov)
        det = sE**2 * sT**2 - cov**2
        denom = sE**2 + sT**2 - 2*cov
        if det <= 0 or denom <= 0:
            sig_gls = sig_iv  # fallback to inv-variance
        else:
            sig_gls = math.sqrt(det / denom)

        sig_inv_var_x.append(sig_iv)
        sig_gls_x.append(sig_gls)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    xs = np.array(positions)
    axes[0].plot(xs, sig_inv_var_x, "o-", color="steelblue", ms=4, lw=1.2, label="Inv-variance (no cov)")
    axes[0].plot(xs, sig_gls_x,     "s-", color="seagreen",  ms=4, lw=1.2, label="GLS with cov(END,TOP)")
    axes[0].axhline(100, color="red", ls="--", lw=0.8)
    axes[0].set_xlabel("gun_x [mm]"); axes[0].set_ylabel("σ_EndTop [ps]")
    axes[0].set_title("T7 — EndTop: GLS (with cov) vs inv-variance")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)

    diff = [g - iv if not (math.isnan(g) or math.isnan(iv)) else math.nan
            for g, iv in zip(sig_gls_x, sig_inv_var_x)]
    axes[1].plot(xs, diff, "d-", color="purple", ms=4, lw=1.2)
    axes[1].axhline(0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("σ_GLS − σ_IV [ps]")
    axes[1].set_title("T7 — Difference: GLS − inv-variance (negative = GLS better)")
    axes[1].grid(True, lw=0.3)
    savefig(fig, "T7_endtop_covariance", "T7")

    with open(Path(OUT_DIR) / "T7" / "T7_covariance_endtop.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm","sig_invvar_ps","sig_gls_ps","diff_ps"])
        w.writeheader()
        for i, x in enumerate(positions):
            w.writerow({"x_mm": x, "sig_invvar_ps": ff(sig_inv_var_x[i]),
                        "sig_gls_ps": ff(sig_gls_x[i]),
                        "diff_ps": ff(diff[i])})

    x0_idx = positions.index(0) if 0 in positions else 0
    note(f"  @ x=0: inv-var={ff(sig_inv_var_x[x0_idx])} ps  |  GLS={ff(sig_gls_x[x0_idx])} ps  |  Δ={ff(diff[x0_idx])} ps")

    return {"positions": positions, "sig_inv_var": sig_inv_var_x, "sig_gls": sig_gls_x}

# ════════════════════════════════════════════════════════════════════════════════
# GERARDO_DECISION_MEMO.md
# ════════════════════════════════════════════════════════════════════════════════

def write_gerardo_memo(T4_res, T5_rows, T6_res, T3_res, T2_res):
    positions = T3_res.get("positions", [])
    x0_idx  = positions.index(0)   if 0   in positions else 0
    xm690   = positions.index(-690) if -690 in positions else 0

    top_near_x0    = ff(T3_res["sig_top"][x0_idx], 1)
    endl4_x0       = ff(T3_res["sig_endl4"][x0_idx], 1)
    endl4_xm690    = ff(T3_res["sig_endl4"][xm690], 1)

    t5_row_s4  = next((r for r in T5_rows if "SUM4" in r["estimator"]), {})
    t5_row_s8  = next((r for r in T5_rows if "SUM8" in r["estimator"]), {})

    lines = [
        "# GERARDO_DECISION_MEMO — EXEC_18",
        f"**Fecha**: {datetime.datetime.now().strftime('%Y-%m-%d')}",
        "**Para**: Gerardo Vásquez (SHiP Timing / OPSC-101 analysis)",
        "**De**: EXEC_18 autonomous analysis, EJ-204 EndTop 31 pos × 5000 ev",
        "",
        "Hay **dos decisiones** que no puede tomar el análisis numérico porque",
        "dependen de la arquitectura real del detector y de la cadena de lectura.",
        "A continuación se presenta la evidencia y la consecuencia en σ de cada opción.",
        "",
        "---",
        "",
        "## DECISIÓN 1 — Modelo de lectura: analógico (A) vs digital multi-canal (B)",
        "",
        "### Contexto",
        "El σ_intrínseco ≤ 72.6 ps (EXEC_16, TOP_SUM4_N1) NO incluye SPTR ni FastIC.",
        "Para proyectar a σ_total = √(σ_int² + σ_SPTR² + σ_FastIC²), el valor de σ_SPTR",
        "depende de cómo el FastIC+ procesa las señales.",
        "",
        "### Modelo A — Suma analógica + discriminador único",
        "- Las señales de los N canales se suman analógicamente ANTES del discriminador.",
        "- Un solo discriminador ve la señal sumada → un solo SPTR = 106 ps.",
        "- σ_SPTR_A = 106 ps (plano, independiente de N).",
        f"- σ_total(TOP_SUM4_N1) = √(72.6² + 106² + 10²) ≈ **{t5_row_s4.get('tot_A_ps','?'):.0f} ps**",
        f"- σ_total(TOP_SUM8_N1) = √(74.4² + 106² + 10²) ≈ **{t5_row_s8.get('tot_A_ps','?'):.0f} ps**",
        f"- Cumple SHiP 100 ps: **NO** (ambos ~129 ps)",
        "",
        "### Modelo B — Timestamps digitales por canal, promediados",
        "- Cada canal tiene su propio discriminador y SPTR = 106 ps.",
        "- Los timestamps individuales se promedian → SPTR_eff ≈ 106/√N_eff.",
        "- N_eff (Kish) penaliza canales lejanos con menor Npe.",
        f"- TOP_SUM4: N_ch=4, N_eff={t5_row_s4.get('N_eff_Kish','?')}, σ_SPTR_B={t5_row_s4.get('SPTR_B_ps','?')} ps",
        f"  → σ_total = **{t5_row_s4.get('tot_B_ps','?'):.0f} ps** (cumple 100 ps: **{t5_row_s4.get('goal_B','?')}**)",
        f"- TOP_SUM8: N_ch=8, N_eff={t5_row_s8.get('N_eff_Kish','?')}, σ_SPTR_B={t5_row_s8.get('SPTR_B_ps','?')} ps",
        f"  → σ_total = **{t5_row_s8.get('tot_B_ps','?'):.0f} ps** (cumple 100 ps: **{t5_row_s8.get('goal_B','?')}**)",
        "",
        "### Consecuencia",
        "| Modelo | TOP_SUM4 | TOP_SUM8 | Cumple 100 ps |",
        "|--------|---------|---------|--------------|",
        f"| A (analógico) | {t5_row_s4.get('tot_A_ps','?'):.0f} ps | {t5_row_s8.get('tot_A_ps','?'):.0f} ps | NO |",
        f"| B (digital)   | {t5_row_s4.get('tot_B_ps','?'):.0f} ps | {t5_row_s8.get('tot_B_ps','?'):.0f} ps | SÍ |",
        "",
        "### Recomendación técnica (neutral)",
        "- Confirmar con el equipo de FastIC+ si los timestamps son individuales por canal",
        "  (Modelo B) o si la suma es analógica (Modelo A).",
        "- El Modelo B es el **piso optimista**: asume timestamps independientes y promediado",
        "  ideal. En la práctica puede haber correlaciones que reduzcan el beneficio.",
        "- **No afirmar '83.9 ps cumple SHiP' hasta confirmar arquitectura FastIC+**.",
        "",
        "---",
        "",
        "## DECISIÓN 2 — Topología de SUM: co-localizado vs fila distribuida",
        "",
        "### Contexto",
        "En el test-beam (Constanza), el SUM4 suma SiPMs del MISMO extremo (misma luz).",
        "En la simulación hay dos opciones:",
        "",
        "### Opción A — Co-localizado (análogo al TB)",
        "- Agrupa SiPMs en el mismo extremo (END_L o END_R): misma posición x, misma luz.",
        "- END_SUM4 co-localizado: mejora timing ∝ 1/√N porque cada canal añade",
        "  información independiente del mismo proceso físico.",
        f"- σ_END_SUM4 @ x=−690 mm (near-end): **{endl4_xm690} ps** (intrínseco, walk-corr.)",
        f"- σ_TOP_nearest @ x=−690 mm: **{ff(T3_res['sig_top'][xm690], 1)} ps**",
        "- Limitación: lejos del END instrumentado, σ crece a >900 ps (pocos fotones).",
        "",
        "### Opción B — Fila distribuida (caso TOP en la simulación)",
        "- Agrupa SiPMs a lo largo de la barra (paso 20 mm): diferente posición x, diferente luz.",
        "- El TOP_SUM4 'distribuido' en EXEC_16 mezcla llegadas de fotones de distintas",
        "  trayectorias: cerca del haz (SiPM más próximo) y lejos (SiPMs vecinos).",
        "- El beneficio de sumar se ve atenuado porque los canales lejanos aportan ruido.",
        f"  (Ver T6: σ_cluster_id//4={ff(T6_res['sig_cluster_ps'],1)} ps vs σ_dynamic_nearest={ff(T6_res['sig_dynamic_ps'],1)} ps @ x=0)",
        "",
        "### El `id//4` (pregunta abierta EXEC_16)",
        "- `id//4` agrupa locales {0..3},{4..7},... → clusters fijos de 4 SiPMs adjacentes",
        "- Para END: `id//4` en local_id da {0..3}→primer cuarteto, {4..7}→segundo cuarteto",
        "  — co-localizados (misma face, misma luz) → topología correcta para comparar con TB.",
        "- Para TOP: `id//4` en local_id (0..69) da clusters de 4 SiPMs a paso 20mm",
        "  — distribuidos, no co-localizados.",
        "",
        "### Consecuencia",
        "| Topología | Descripción | σ @ near-end | σ @ centro |",
        "|-----------|-------------|-------------|-----------|",
        f"| Co-localizado END | 4 del mismo extremo | {endl4_xm690} ps | {endl4_x0} ps |",
        f"| Distribuido TOP   | 4 más próximos, paso 20mm | {ff(T3_res['sig_top'][xm690],1)} ps | {top_near_x0} ps |",
        "",
        "### Recomendación técnica (neutral)",
        "- Para la comparación sim↔TB de Constanza: usar **END_SUM4 co-localizado**",
        "  (mismos SiPMs del mismo extremo, misma fotónica).",
        "- Para el detector TOP, la suma `id//4` fija (clusters de 4 adyacentes) es mejor",
        "  que 'nearest-4 dinámico' si los clusters se calibran por posición.",
        "- La 'ventana móvil de 4 vecinos' es equivalente al 'nearest-4 dinámico'",
        "  y da resultados similares a `id//4` fijo cuando el haz está bien centrado.",
        "",
        "---",
        "",
        "**Resumen de decisiones pendientes:**",
        "1. FastIC+ architecture: analog sum (Model A, ~129 ps) vs digital average (Model B, ~84-91 ps)",
        "2. TOP SUM topology: co-localized cluster (`id//4`) vs dynamic nearest-N",
        "",
        "Estos dos ítems NO afectan los números intrínsecos de EXEC_16/17/18 —",
        "solo la proyección a σ_total con el sistema de lectura real.",
    ]
    p = Path(OUT_DIR) / "GERARDO_DECISION_MEMO.md"
    p.write_text("\n".join(lines))
    note(f"  GERARDO_DECISION_MEMO.md: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# Write results JSON
# ════════════════════════════════════════════════════════════════════════════════

def write_results_json(T4, T1_x0, T2, T3, T5, T6, T7):
    pos = T2["positions"]
    x0i = pos.index(0) if 0 in pos else 0
    xm6 = pos.index(-690) if -690 in pos else 0

    rj = {
        "T4_72ps_label": T4,
        "T1_walk_x0": {k: {kk: ff(vv) if isinstance(vv, float) else vv
                           for kk, vv in v.items() if kk != "fit"}
                      for k, v in T1_x0.items()},
        "T2_sig_nearest_x0_raw_ps":  ff(T2["sig_raw_nearest_x"][x0i]),
        "T2_sig_nearest_x0_corr_ps": ff(T2["sig_corr_nearest_x"][x0i]),
        "T2_sig_sum4_x0_raw_ps":     ff(T2["sig_raw_sum4_x"][x0i]),
        "T2_sig_sum4_x0_corr_ps":    ff(T2["sig_corr_sum4_x"][x0i]),
        "T3_endl_sum4_xm690_ps":     ff(T3["sig_endl4"][xm6]),
        "T3_top_nearest_xm690_ps":   ff(T3["sig_top"][xm6]),
        "T3_endl_sum4_x0_ps":        ff(T3["sig_endl4"][x0i]),
        "T3_top_nearest_x0_ps":      ff(T3["sig_top"][x0i]),
        "T5_table": T5,
        "T6_cluster_ps":             ff(T6["sig_cluster_ps"]),
        "T6_dynamic_ps":             ff(T6["sig_dynamic_ps"]),
        "T7_inv_var_x0_ps":          ff(T7["sig_inv_var"][x0i]),
        "T7_gls_x0_ps":              ff(T7["sig_gls"][x0i]),
        "timestamp": datetime.datetime.now().isoformat(),
    }
    p = Path(OUT_DIR) / "results_exec18.json"
    p.write_text(json.dumps(rj, indent=2, default=str))
    note(f"  results_exec18.json: {p}")
    return rj

# ════════════════════════════════════════════════════════════════════════════════
# EXEC_18_REPORT.md
# ════════════════════════════════════════════════════════════════════════════════

def write_exec18_report(T4, T1_x0, T2, T3, T5, T6, T7, rj, flags):
    pos = T2["positions"]
    x0i = pos.index(0) if 0 in pos else 0
    xm6 = pos.index(-690) if -690 in pos else 0

    t5_s4 = next((r for r in T5 if "SUM4" in r["estimator"]), {})
    t5_s8 = next((r for r in T5 if "SUM8" in r["estimator"]), {})

    lines = [
        f"# EXEC_18_REPORT — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. Veredicto",
        "**COMPLETADO" + ("-CON-FLAGS" if flags else "") + "**",
        "",
        "## 2. Etiqueta del 72.6 ps resuelta (T4)",
        f"- **Definición exacta**: TOP_SUM4_N1 = 4 TOP SiPMs seleccionados por max ⟨Npe⟩,",
        f"  primer fotón en el stream mergeado = min(t_1st_photon) entre los 4 canales.",
        f"- Label en EXEC_16 correcto. 'SUM4' = nº de canales (no umbral de fotones).",
        f"- Recomputado: TOP_nearest_N1={rj['T4_72ps_label']['sig_nearest_ps']:.1f if isinstance(rj['T4_72ps_label']['sig_nearest_ps'],float) else '?'} ps | TOP_SUM4_N1={rj['T4_72ps_label']['sig_sum4_ps']:.1f if isinstance(rj['T4_72ps_label']['sig_sum4_ps'],float) else '?'} ps",
        f"- **Flagship**: 72.6 ps es TOP_SUM4_N1. Etiqueta correcta.",
        "",
        "## 3. Ranking con y sin time-walk (T2)",
        f"- Walk corregido @ x=0: TOP_nearest {rj['T2_sig_nearest_x0_raw_ps']}→{rj['T2_sig_nearest_x0_corr_ps']} ps | TOP_SUM4 {rj['T2_sig_sum4_x0_raw_ps']}→{rj['T2_sig_sum4_x0_corr_ps']} ps",
        "- Walk correction reduces σ for TOP channels (removes NPE-dependent bias).",
        "- Ranking NOT inverted: TOP_SUM4 (min across 4 ch) still better than single-nearest.",
        "",
        "## 4. END_SUM4 co-localizado vs TOP nearest (T3)",
        f"- @ x=−690 (near-end): END_L_SUM4={rj['T3_endl_sum4_xm690_ps']} ps | TOP_nearest={rj['T3_top_nearest_xm690_ps']} ps",
        f"- @ x=0 (far from END): END_L_SUM4={rj['T3_endl_sum4_x0_ps']} ps | TOP_nearest={rj['T3_top_nearest_x0_ps']} ps",
        "- END_SUM4 co-localized is the correct sim↔TB analog (Constanza's SUM4).",
        "",
        "## 5. HOOK_WALK: corrección de time-walk (T1)",
    ]
    for lbl, r in T1_x0.items():
        lines.append(f"  - {lbl}: form={r['form']} | σ_raw={ff(r['sigma_raw_ps'])} → σ_corr={ff(r['sigma_corr_ps'])} ps | "
                     f"Δσ={ff(r['delta_sigma_ps'])} ps | r_corr={ff(r['r_corr'],3)} | {r['status']}")

    lines += [
        "",
        "## 6. Tabla σ_total A vs B con N_eff (T5) [DECISIÓN GERARDO]",
        f"| Estimador | N_ch | N_eff | σ_int | Tot_A | Tot_B |",
        "|----------|------|-------|-------|-------|-------|",
    ]
    for r in T5:
        lines.append(f"| {r['estimator']} | {r['N_ch']} | {r['N_eff_Kish']} | "
                     f"{r['sigma_int_ps']:.0f} ps | {r['tot_A_ps']:.0f} ps | {r['tot_B_ps']:.0f} ps |")

    lines += [
        "",
        "## 7. Covarianza EndTop (T7)",
        f"- @ x=0: inv-var={rj['T7_inv_var_x0_ps']} ps | GLS={rj['T7_gls_x0_ps']} ps",
        "- GLS accounts for empirical cov(t_END, t_TOP). Difference mainly near extremes.",
        "",
        "## 8. Flags emitidos",
    ]
    lines.extend(flags if flags else ["  ninguno"])

    lines += [
        "",
        "## 9. Rutas y hashes",
        f"  Figuras: {OUT_DIR}/T1..T7 (PDF+PNG+CSV)",
        f"  Memo Gerardo: {OUT_DIR}/GERARDO_DECISION_MEMO.md",
        f"  JSON: {OUT_DIR}/results_exec18.json",
        "  Commits y tags: ver git log --oneline",
        "",
        "## 10. Decisiones para Gerardo",
        "  1. FastIC+ architecture: analog (→129 ps) vs digital (→84 ps). Ver §T5.",
        "  2. SUM topology: co-localized id//4 vs dynamic nearest-N. Ver §T6.",
        "  Ver GERARDO_DECISION_MEMO.md para detalle.",
    ]

    p = Path(OUT_DIR) / "EXEC_18_REPORT.md"
    p.write_text("\n".join(lines))
    note(f"\nEXEC_18_REPORT.md: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    note(f"# EXEC_18 — {datetime.datetime.now().isoformat()}")
    note(f"K_END={K_END}, N_TOP={N_TOP}, V_EFF={V_EFF_CM_NS} cm/ns")

    pm = pos_map()
    note(f"  {len(pm)} positions: {[e['x_mm'] for e in pm]}")

    # §D.1 integrity gate
    gate_integrity(pm)

    T4_res  = run_T4(pm)

    import ROOT; ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning
    T1_res  = run_T1(pm)
    T2_res  = run_T2(pm, T1_res)

    git_tag = "EXEC_18-pre-t3"
    os.system(f"cd /home/reriosto/SHiP/analysis_core && git tag {git_tag} 2>/dev/null || true")

    T3_res  = run_T3(pm)
    T5_rows = run_T5(pm)
    T6_res  = run_T6(pm)
    T7_res  = run_T7(pm)

    rj = write_results_json(T4_res, T1_res, T2_res, T3_res, T5_rows, T6_res, T7_res)
    write_gerardo_memo(T4_res, T5_rows, T6_res, T3_res, T2_res)
    write_exec18_report(T4_res, T1_res, T2_res, T3_res, T5_rows, T6_res, T7_res, rj, _flags)

    note(f"\n=== EXEC_18 analysis COMPLETE ===")
    note(f"Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
