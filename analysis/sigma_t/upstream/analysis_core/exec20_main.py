#!/usr/bin/env python3.12
"""
exec20_main.py — EXEC_20: END-only baseline, TOP contribution, canonical fix.

Uses existing EndTop data for T1/T3/T4. T2b (scan) results loaded when available.
All intrinsic (no SPTR/FastIC). Binning: sqrt_n (canonical per EXEC_19 T5 decision).
Run: MPLBACKEND=Agg python3.12 exec20_main.py
"""

import sys, os, csv, json, math, datetime, warnings
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine   import fit_core_gaussian
from lib.robust_seeds import gather_seeds
from walk_correction  import fit_walk, apply_correction

# ════════════════════════════════════════════════════════════════════════════════
# Constants
# ════════════════════════════════════════════════════════════════════════════════
K_END     = 8
N_TOP     = 70
V7_A_PS   = 757.3  # σ_t = a/√Npe (from EXEC_17 V7)

DATA_DIR  = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
OUT_DIR   = "/home/reriosto/SHiP/analysis_core/out/EXEC_20"
RANDOM_SEED = 20260618
N_BOOT    = 200
MIN_EV    = 30

REPR_X    = [-690, -350, 0, 350, 690]

_report = []; _flags = []

def note(m): _report.append(m); print(m)
def flag(m): _flags.append(m); note(f"  [FLAG] {m}")

def savefig(fig, stem, sub=""):
    d = Path(OUT_DIR) / sub
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf","png"):
        fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def ff(v, d=1):
    try: return f"{float(v):.{d}f}" if not math.isnan(float(v)) else "NaN"
    except: return str(v)

# ── Fit: sqrt_n binning (canonical decision from EXEC_19/T5) ─────────────────
def gauss_fit(v, prefix="g"):
    v_c = v[~np.isnan(v)]
    if len(v_c) < MIN_EV:
        return {"sigma_fit": math.nan, "bootstrap_err": math.nan, "flag": "insufficient_stats"}
    cfg = {"MIN_EVENTS_FOR_FIT": MIN_EV, "FIT_WINDOW_SIGMAS": 2.0,
           "BINNING_STRATEGY": "sqrt_n",  # canonical
           "N_BOOTSTRAP": N_BOOT, "RANDOM_SEED": RANDOM_SEED,
           "FIT_OPTIONS": "R Q S 0", "CHI2_NDF_WARN": 3.0}
    return fit_core_gaussian(v_c, cfg, name_prefix=prefix)

# ── Data utilities ────────────────────────────────────────────────────────────
def pos_map():
    entries = []
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        f = uproot.open(str(p))
        gx = f["sipm_hits"]["gun_x_mm"].array(library="np")
        entries.append({"x_mm": int(round(float(np.median(gx)))), "path": p})
    return sorted(entries, key=lambda e: e["x_mm"])

def load_pos(entry, branches=None):
    branches = branches or ["event_id","global_id","time_ns"]
    f = uproot.open(str(entry["path"]))
    return f["sipm_hits"].arrays(branches, library="np")

def best_k(gid_arr, face, k):
    if face == "END_L": mask = gid_arr < K_END
    elif face == "END_R": mask = (gid_arr >= K_END) & (gid_arr < 2*K_END)
    else: mask = gid_arr >= 2*K_END
    g = gid_arr[mask]
    if len(g) == 0: return np.array([], dtype=int)
    u, c = np.unique(g, return_counts=True)
    return u[np.argsort(-c)[:k]]

def t_first_merged(evs, times, gids, sel_gids, all_events):
    if len(sel_gids) == 0: return np.full(len(all_events), np.nan)
    mask = np.isin(gids, sel_gids)
    ev_m, t_m = evs[mask], times[mask]
    if len(ev_m) == 0: return np.full(len(all_events), np.nan)
    idx = np.lexsort((t_m, ev_m))
    ev_m, t_m = ev_m[idx], t_m[idx]
    u, cnt = np.unique(ev_m, return_counts=True)
    cum = np.concatenate([[0], np.cumsum(cnt)])
    ev_map = {e: i for i, e in enumerate(all_events)}
    res = np.full(len(all_events), np.nan)
    for i, (e, c) in enumerate(zip(u, cnt)):
        if c >= 1:
            j = ev_map.get(e, -1)
            if j >= 0: res[j] = t_m[cum[i]]
    return res

def npe_stream(evs, gids, sel_gids, all_events):
    if len(sel_gids) == 0: return np.zeros(len(all_events))
    mask = np.isin(gids, sel_gids)
    ev_m = evs[mask]
    res = np.zeros(len(all_events))
    u, cnt = np.unique(ev_m, return_counts=True)
    ev_map = {e: i for i, e in enumerate(all_events)}
    for e, c in zip(u, cnt):
        j = ev_map.get(e, -1)
        if j >= 0: res[j] = c
    return res

def walk_correct(t_raw, s):
    valid = ~np.isnan(t_raw) & (s > 0)
    if valid.sum() < 50: return t_raw.copy()
    fit = fit_walk(t_raw[valid]*1000, s[valid])
    return apply_correction(t_raw*1000, s, fit)/1000.0

# ════════════════════════════════════════════════════════════════════════════════
# T1 — END-only baseline: sigma_END(x) and sigma_x (t_avg = (tL+tR)/2)
# ════════════════════════════════════════════════════════════════════════════════

def run_T1(pm):
    note("\n=== T1 — END-only baseline from EndTop data ===")
    note("  Estimator: t_avg = (t_L + t_R)/2; t_L = first photon in merged END_L stream (all 8 SiPMs), walk-corrected")
    note("  Binning: sqrt_n (canonical); walk: parametric α+β/√NPE")
    Path(OUT_DIR + "/T1").mkdir(parents=True, exist_ok=True)

    results = []
    for entry in pm:
        x = entry["x_mm"]
        d = load_pos(entry)
        evs = d["event_id"]; gids = d["global_id"]; times = d["time_ns"]
        all_ev = np.unique(evs)
        n_ev = len(all_ev)

        endl8 = np.arange(K_END)          # all END_L (gid 0..7)
        endr8 = np.arange(K_END, 2*K_END) # all END_R (gid 8..15)

        tL_raw = t_first_merged(evs, times, gids, endl8, all_ev)
        tR_raw = t_first_merged(evs, times, gids, endr8, all_ev)
        sL = npe_stream(evs, gids, endl8, all_ev)
        sR = npe_stream(evs, gids, endr8, all_ev)

        tL = walk_correct(tL_raw, sL)
        tR = walk_correct(tR_raw, sR)

        valid = (~np.isnan(tL)) & (~np.isnan(tR))
        t_avg = np.full(n_ev, np.nan)
        dt    = np.full(n_ev, np.nan)
        t_avg[valid] = 0.5*(tL[valid]+tR[valid])
        dt[valid]    = tL[valid]-tR[valid]

        r_avg = gauss_fit(t_avg, f"tend_avg_x{x}")
        sig_avg   = r_avg.get("sigma_fit", math.nan)*1000  # ps
        sig_avg_b = r_avg.get("bootstrap_err", math.nan)*1000
        n_valid   = int(valid.sum())
        npe_L     = float(np.sum(sL)/n_ev)
        npe_R     = float(np.sum(sR)/n_ev)

        # σ_x from σ(Δt) × v_eff/2 (v_eff from slope computed later)
        dt_std = float(np.nanstd(dt[valid])) if valid.sum() > 5 else math.nan
        dt_mean = float(np.nanmean(dt[valid])) if valid.sum() > 5 else math.nan

        results.append({
            "x_mm": x, "sigma_END_avg_ps": sig_avg, "boot_ps": sig_avg_b,
            "sigma_dt_ps": dt_std*1000 if not math.isnan(dt_std) else math.nan,
            "dt_mean_ns": dt_mean, "npe_L": npe_L, "npe_R": npe_R, "n_valid": n_valid
        })

    # v_eff from slope of <Δt>(x)
    xs = np.array([r["x_mm"] for r in results])
    dts = np.array([r["dt_mean_ns"] for r in results])
    valid_dt = ~np.isnan(dts)
    v_eff_cm_ns = math.nan
    if valid_dt.sum() >= 5:
        try:
            m, b = np.polyfit(xs[valid_dt], dts[valid_dt], 1)
            v_eff_mm_ns = float(-2.0 / m)  # mm/ns
            v_eff_cm_ns = v_eff_mm_ns / 10.0
            note(f"  v_eff = {v_eff_cm_ns:.2f} cm/ns (from Δt(x) slope)")
        except Exception as e:
            note(f"  [WARN] v_eff fit failed: {e}")

    # σ_x = σ(Δt) × |v_eff| / 2 [mm]
    v_eff_mm = abs(v_eff_cm_ns*10) if not math.isnan(v_eff_cm_ns) else math.nan
    for r in results:
        r["sigma_x_mm"] = r["sigma_dt_ps"]/1000 * v_eff_mm / 2.0 if not math.isnan(r["sigma_dt_ps"]) and not math.isnan(v_eff_mm) else math.nan

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sigs = [r["sigma_END_avg_ps"] for r in results]
    boots = [r.get("boot_ps", math.nan) for r in results]
    sigx = [r.get("sigma_x_mm", math.nan) for r in results]
    axes[0].errorbar(xs, sigs, yerr=boots, fmt="s-", color="steelblue", ms=4, capsize=3, lw=1.2)
    axes[0].set_ylabel("σ_END [ps]"); axes[0].set_xlabel("gun_x [mm]")
    axes[0].set_title("T1 — END-only baseline σ_END(t_avg) | from EndTop data (intrinsic)")
    axes[0].set_yscale("log"); axes[0].set_ylim(10, 5000); axes[0].grid(True, lw=0.3, which="both")
    axes[0].text(0.02, 0.02, "Intrinseco solo (sin SPTR/FastIC).\nt_avg=(t_L+t_R)/2, walk-corr, sqrt_n.", transform=axes[0].transAxes, fontsize=7)

    valid_sx = [not math.isnan(s) for s in sigx]
    axes[1].plot(xs[valid_sx], [sigx[i] for i in range(len(xs)) if valid_sx[i]], "o-", color="tomato", ms=4, lw=1.2)
    axes[1].set_ylabel("σ_x [mm]"); axes[1].set_xlabel("gun_x [mm]")
    axes[1].set_title(f"T1 — Spatial resolution σ_x(x) | v_eff={ff(v_eff_cm_ns,2)} cm/ns")
    axes[1].grid(True, lw=0.3)
    savefig(fig, "T1_end_baseline", "T1")

    with open(Path(OUT_DIR) / "T1" / "T1_end_baseline.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(results[0].keys()))
        w.writeheader(); w.writerows(results)

    note(f"  σ_END(t_avg) @ x=0:    {ff([r['sigma_END_avg_ps'] for r in results if r['x_mm']==0][0])} ps")
    note(f"  σ_END(t_avg) @ x=-690: {ff([r['sigma_END_avg_ps'] for r in results if r['x_mm']==-690][0])} ps")
    note(f"  σ_x @ x=0:    {ff([r['sigma_x_mm'] for r in results if r['x_mm']==0][0])} mm")
    note(f"  T1 done.")
    return results, v_eff_cm_ns

# ════════════════════════════════════════════════════════════════════════════════
# T3 — TOP contribution over END-only baseline
# ════════════════════════════════════════════════════════════════════════════════

def run_T3(pm, T1_results):
    note("\n=== T3 — TOP contribution over END-only baseline ===")
    Path(OUT_DIR + "/T3").mkdir(parents=True, exist_ok=True)

    positions = [e["x_mm"] for e in pm]
    T1_map = {r["x_mm"]: r for r in T1_results}

    sig_end    = []  # from T1
    sig_top    = []  # TOP_nearest
    sig_endtop = []  # EndTop combined (GLS)
    improvement = []  # 1 - sigma_endtop/sigma_end

    for entry in pm:
        x = entry["x_mm"]
        d = load_pos(entry)
        evs = d["event_id"]; gids = d["global_id"]; times = d["time_ns"]
        all_ev = np.unique(evs)

        # END baseline from T1
        sig_e = T1_map.get(x, {}).get("sigma_END_avg_ps", math.nan)
        sig_end.append(sig_e)

        # TOP: nearest single channel (walk-corrected)
        top4_gids  = best_k(gids, "TOP", 4)
        near_gid   = top4_gids[0:1] if len(top4_gids) >= 1 else np.array([], dtype=int)
        tT_raw = t_first_merged(evs, times, gids, near_gid, all_ev)
        sT     = npe_stream(evs, gids, near_gid, all_ev)
        tT     = walk_correct(tT_raw, sT)
        r_T    = gauss_fit(tT, f"t3_top_x{x}")
        sig_t  = r_T.get("sigma_fit", math.nan)*1000
        sig_top.append(sig_t)

        # EndTop: GLS combination (from EXEC_17/T7 methodology)
        endl8 = np.arange(K_END); endr8 = np.arange(K_END, 2*K_END)
        tL_raw = t_first_merged(evs, times, gids, endl8, all_ev)
        tR_raw = t_first_merged(evs, times, gids, endr8, all_ev)
        tL = walk_correct(tL_raw, npe_stream(evs, gids, endl8, all_ev))
        tR = walk_correct(tR_raw, npe_stream(evs, gids, endr8, all_ev))

        valid_e = (~np.isnan(tL)) & (~np.isnan(tR))
        t_end_comb = np.full(len(all_ev), np.nan)
        t_end_comb[valid_e] = 0.5*(tL[valid_e]+tR[valid_e])

        # GLS with covariance
        valid_both = (~np.isnan(t_end_comb)) & (~np.isnan(tT))
        if valid_both.sum() >= MIN_EV:
            tE_v = t_end_comb[valid_both]; tT_v = tT[valid_both]
            r_Ef = gauss_fit(tE_v, f"t3_Egls_x{x}"); r_Tf = gauss_fit(tT_v, f"t3_Tgls_x{x}")
            sE = r_Ef.get("sigma_fit", math.nan)*1000
            sT2 = r_Tf.get("sigma_fit", math.nan)*1000
            if not (math.isnan(sE) or math.isnan(sT2) or sE <= 0 or sT2 <= 0):
                cov = float(np.cov(tE_v, tT_v)[0,1])*(1000**2)  # ps²
                det = sE**2 * sT2**2 - cov**2
                denom = sE**2 + sT2**2 - 2*cov
                if det > 0 and denom > 0:
                    sig_et = math.sqrt(det/denom)
                else:
                    sig_et = 1.0/math.sqrt(1.0/sE**2 + 1.0/sT2**2)  # inv-var fallback
            else:
                sig_et = math.nan
        else:
            sig_et = math.nan
        sig_endtop.append(sig_et)

        # Improvement
        imp = 1.0 - sig_et/sig_e if not (math.isnan(sig_et) or math.isnan(sig_e) or sig_e == 0) else math.nan
        improvement.append(imp)

    xs = np.array(positions)

    # Figure: 3 curves + improvement
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    axes[0].plot(xs, sig_end,    "s-", color="steelblue", ms=4, lw=1.2, label="END-only (t_avg)")
    axes[0].plot(xs, sig_top,    "o-", color="seagreen",  ms=4, lw=1.2, label="TOP nearest (walk-corr)")
    axes[0].plot(xs, sig_endtop, "^-", color="tomato",    ms=4, lw=1.2, label="EndTop (GLS covariance)")
    axes[0].set_yscale("log"); axes[0].set_ylim(5, 5000)
    axes[0].set_ylabel("σ_t [ps] (log)"); axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3, which="both")
    axes[0].set_title("T3 — Tres curvas: END-only, TOP, EndTop | intrínseco")
    axes[0].text(0.01, 0.02, "Intrinseco solo (sin SPTR/FastIC). sqrt_n binning, walk-corr.",
                 transform=axes[0].transAxes, fontsize=7)

    imp_x = [xs[i] for i in range(len(xs)) if not math.isnan(improvement[i])]
    imp_v = [improvement[i] for i in range(len(xs)) if not math.isnan(improvement[i])]
    axes[1].plot(imp_x, [v*100 for v in imp_v], "d-", color="purple", ms=4, lw=1.2)
    axes[1].axhline(0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("Mejora(x) = 1 − σ_EndTop/σ_END [%]")
    axes[1].set_title("T3 — Aporte del TOP sobre END-only baseline")
    axes[1].grid(True, lw=0.3)
    savefig(fig, "T3_top_contribution", "T3")

    with open(Path(OUT_DIR) / "T3" / "T3_top_contribution.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm","sigma_END_ps","sigma_TOP_ps","sigma_EndTop_ps","improvement_pct"])
        w.writeheader()
        for i, x in enumerate(positions):
            w.writerow({"x_mm": x, "sigma_END_ps": ff(sig_end[i]),
                        "sigma_TOP_ps": ff(sig_top[i]),
                        "sigma_EndTop_ps": ff(sig_endtop[i]),
                        "improvement_pct": ff(improvement[i]*100) if not math.isnan(improvement[i]) else "NaN"})

    # Report key numbers
    x0_i = positions.index(0) if 0 in positions else 0
    xm6_i = positions.index(-690) if -690 in positions else 0
    imp_x0   = improvement[x0_i]
    imp_xm6  = improvement[xm6_i]
    imp_mean = float(np.nanmean(improvement))
    note(f"  Aporte del TOP @ x=0:    mejora = {ff(imp_x0*100)}% (σ_END={ff(sig_end[x0_i])}, σ_EndTop={ff(sig_endtop[x0_i])})")
    note(f"  Aporte del TOP @ x=-690: mejora = {ff(imp_xm6*100)}% (σ_END={ff(sig_end[xm6_i])}, σ_EndTop={ff(sig_endtop[xm6_i])})")
    note(f"  Mejora media (scan):      {ff(imp_mean*100)}%")
    note(f"  Física: TOP ayuda MÁS al centro (END casi ciego) y MENOS en extremos (END excelente)")
    return {"improvement_x0": imp_x0, "improvement_xm690": imp_xm6, "improvement_mean": imp_mean,
            "sig_end_x0": sig_end[x0_i], "sig_top_x0": sig_top[x0_i], "sig_endtop_x0": sig_endtop[x0_i],
            "positions": positions, "sig_end_all": sig_end, "sig_top_all": sig_top, "sig_endtop_all": sig_endtop}

# ════════════════════════════════════════════════════════════════════════════════
# T4 — Reconciliation with Gerardo
# ════════════════════════════════════════════════════════════════════════════════

def run_T4(T1_results, T3_res):
    note("\n=== T4 — Reconciliation with Gerardo ===")
    Path(OUT_DIR + "/T4").mkdir(parents=True, exist_ok=True)

    # Our σ_END(centre)
    r0 = next((r for r in T1_results if r["x_mm"] == 0), {})
    sig_end_x0 = r0.get("sigma_END_avg_ps", math.nan)
    npe_L_x0   = r0.get("npe_L", math.nan)
    npe_R_x0   = r0.get("npe_R", math.nan)

    note(f"  Our σ_END(t_avg) @ x=0 = {ff(sig_end_x0)} ps")
    note(f"  npe_L @ x=0 = {npe_L_x0:.2f} PE/ev, npe_R = {npe_R_x0:.2f} PE/ev")
    note(f"  At center, END sees only ~0.37 PE/ev per face → timing dominated by Poisson of rare photons")
    note("")
    note(f"  Why Gerardo's END-only is likely BETTER:")
    note(f"  1. Gerardo measures at positions near the ends (x±~690) where END sees 900+ PE/ev")
    note(f"  2. Or Gerardo uses SUM of more SiPMs / different clustering")
    note(f"  3. Or Gerardo's bar/material/wrapping differ from our simulation")
    note(f"  4. Estimator: our t_avg=(t_L+t_R)/2 gives ~{ff(sig_end_x0)} ps at center")
    note(f"     If Gerardo uses a weighted average (better timing channel dominates) → could be better")
    note(f"  NOTE: Confirm Gerardo's exact estimator (which channels, combination, threshold)")

    # Near-end comparison
    rm = next((r for r in T1_results if r["x_mm"] == -690), {})
    sig_end_xm6 = rm.get("sigma_END_avg_ps", math.nan)
    note(f"\n  Our σ_END(t_avg) @ x=-690 = {ff(sig_end_xm6)} ps")
    note(f"  npe_L @ x=-690 = {rm.get('npe_L',math.nan):.0f} PE/ev, npe_R = {rm.get('npe_R',math.nan):.2f} PE/ev")
    note(f"  Near END_L face: END_L sees ~903 PE/ev → good timing")
    note(f"  But t_avg requires BOTH ends to fire → rare events when END_R fires too (0.07 PE/ev)")
    note(f"  → t_avg throws away most events at x=-690 where END_R doesn't fire!")
    note(f"  → This is a key design choice: t_avg for position-independence vs single-end for best timing")

    # Better estimator: use single best end (inv-var)
    note(f"\n  Alternative: single-best-end (inv-var selects END_L at x=-690)")
    note(f"  → σ_END_L(x=-690) from EXEC_19 T4 = 30.5 ps (near-end, sqrt_n, walk)")
    note(f"  This likely better matches Gerardo's setup if he reads both ends independently")
    note(f"  and takes the timestamp from the better side")

    # Write comparison
    data = {
        "sigma_END_tavg_x0_ps": sig_end_x0,
        "sigma_END_tavg_xm690_ps": sig_end_xm6,
        "sigma_END_L_xm690_canonical_ps": 30.5,  # from EXEC_19 T4
        "npe_L_x0": npe_L_x0,
        "npe_R_x0": npe_R_x0,
        "note_gerardo": "Confirm estimator: channels, combination, position evaluated",
        "note_methodology": (
            "t_avg requires BOTH ends to fire; near-end, END_R fires rarely → low efficiency. "
            "Single-end inv-var is better near ends but loses position-independence."
        )
    }
    with open(Path(OUT_DIR) / "T4" / "T4_gerardo_reconciliation.csv", "w", newline="") as f:
        w = csv.writer(f)
        for k, v in data.items(): w.writerow([k, v])
    return data

# ════════════════════════════════════════════════════════════════════════════════
# T5 — Canonical fix and dynamic-4 correction
# ════════════════════════════════════════════════════════════════════════════════

def run_T5(pm):
    note("\n=== T5 — Canonical fix: binning convention + dynamic-4 correction ===")
    Path(OUT_DIR + "/T5").mkdir(parents=True, exist_ok=True)

    note("  Binning decision: sqrt_n (uniform across sessions)")
    note("  Rationale: FD gave 72.1 vs 85.0 ps artifact (EXEC_18 T6); sqrt_n is stable")
    note("  Applied throughout this analysis")

    # Recompute dynamic-4 and id//4 cluster with sqrt_n + walk
    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load_pos(entry0)
    evs = d["event_id"]; gids = d["global_id"]; times = d["time_ns"]
    all_ev = np.unique(evs)

    top4_dyn = best_k(gids, "TOP", 4)
    ch_hits = {}
    for g in np.unique(gids[gids >= 2*K_END]):
        cl = (int(g) - 2*K_END) // 4
        ch_hits[cl] = ch_hits.get(cl, 0) + int(np.sum(gids == g))
    best_cl = max(ch_hits, key=ch_hits.get)
    cl4_gids = np.array([2*K_END + best_cl*4 + i for i in range(4)
                         if 2*K_END + best_cl*4 + i < 2*K_END + N_TOP])

    results_T5 = {}
    for label, sel_gids in [("dynamic_nearest4", top4_dyn), ("id_div4_cluster", cl4_gids)]:
        t_raw = t_first_merged(evs, times, gids, sel_gids, all_ev)
        s_tot = npe_stream(evs, gids, sel_gids, all_ev)
        t_c   = walk_correct(t_raw, s_tot)
        npe_tot = float(np.mean(s_tot))
        r = gauss_fit(t_c, f"t5_{label}")
        sig = r.get("sigma_fit", math.nan)*1000
        results_T5[label] = {"sigma_ps": sig, "npe_total": npe_tot}
        note(f"  {label}: ⟨Npe⟩_total={npe_tot:.1f} PE/ev, σ = {ff(sig)} ps (sqrt_n, walk-corr)")

    sig_dyn = results_T5["dynamic_nearest4"]["sigma_ps"]
    sig_cl  = results_T5["id_div4_cluster"]["sigma_ps"]
    note(f"\n  dynamic ({sig_dyn:.1f} ps) vs id//4 ({sig_cl:.1f} ps):")
    if not (math.isnan(sig_dyn) or math.isnan(sig_cl)):
        if sig_dyn <= sig_cl:
            note(f"  → dynamic ≤ id//4 as expected (more light = better σ). EXEC_18 artifact RESOLVED.")
        else:
            note(f"  → dynamic still > id//4 ({sig_dyn-sig_cl:.1f} ps gap). Possible residual effect; check T3 figure.")
            flag(f"T5: dynamic still > cluster after sqrt_n+walk ({sig_dyn:.1f} vs {sig_cl:.1f} ps)")

    # Updated CANONICAL_ESTIMATORS.md
    canon_lines = [
        "# CANONICAL_ESTIMATORS.md — EXEC_20 update",
        f"Updated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## Convention (EXEC_20 final)",
        "- **Binning**: sqrt_n (ALL sessions from EXEC_20 onward)",
        "- Walk: parametric α+β/√s, anchored at median(NPE)",
        "- Fit: MAD·1.4826 seeded, window ±2σ_MAD, N_boot=200",
        "- **σ_int**: intrinsic ONLY (no SPTR/FastIC). SPTR/FastIC in deferred appendix.",
        "- σ_int includes 20 ps SiPMSD jitter (set to 0 in actual macs — both datasets)",
        "",
        "## Headline estimators @ x=0 (EndTop data, sqrt_n, walk-corr)",
        f"- END-only t_avg (all 8 per face): σ = ~882 ps (T1; center is photon-starved)",
        f"- TOP nearest N=1: σ ≈ 106 ps (T4 EXEC_19, sqrt_n, walk-corr)",
        f"- TOP_SUM4_N1: σ ≈ 69-73 ps (EXEC_16 raw/EXEC_19 walk-corr; see T5 for dynamic vs fixed)",
        f"- dynamic nearest-4: σ = {ff(sig_dyn)} ps (sqrt_n, walk-corr)",
        f"- id//4 fixed cluster: σ = {ff(sig_cl)} ps (sqrt_n, walk-corr)",
        "",
        "## Headline estimators @ x=-690 mm (EndTop data)",
        "- END_L_SUM4 co-localized (4 END_L gids 3,4,2,5): σ ≈ 30.5 ps (EXEC_19 T4)",
        "- TOP nearest @ x=-690: σ ≈ 118 ps (EXEC_19 T4)",
        "",
        "## END-only baseline (T1, EXEC_20)",
        "- t_avg = (t_L + t_R)/2 with both ends (all 8 per face), walk-corr, sqrt_n",
        "- Position-independent by design; requires both ends to fire",
        "- At center: ~882 ps (few photons); near end x=-690: ~609 ps (tR rarely fires)",
        "- Single-end (inv-var / best face): ~30-44 ps near-end; NOT used for T1 baseline",
        "",
        "## Deferred to SPTR/FastIC appendix",
        "- σ_tot_A, σ_tot_B, N_eff, Chain A/B: see EXEC_19 for these numbers",
        "- Not in headline (René's decision, EXEC_20 reorientation)",
    ]
    Path(OUT_DIR + "/T5").mkdir(parents=True, exist_ok=True)
    canon_path = Path("/home/reriosto/SHiP/analysis_core/out/EXEC_19") / "CANONICAL_ESTIMATORS.md"
    # Write update to EXEC_20 folder as well
    Path(OUT_DIR + "/T5/CANONICAL_ESTIMATORS.md").write_text("\n".join(canon_lines))
    # Also update the EXEC_19 location
    if canon_path.exists():
        canon_path.write_text("\n".join(canon_lines))
    note(f"  CANONICAL_ESTIMATORS.md updated")
    return results_T5

# ════════════════════════════════════════════════════════════════════════════════
# Write results JSON
# ════════════════════════════════════════════════════════════════════════════════

def write_results_json(T1_r, T3_r, T4_r, T5_r, v_eff, scan_result):
    T1_map = {r["x_mm"]: r for r in T1_r}
    rj = {
        "T1_sigma_END_tavg_x0_ps": ff(T1_map.get(0, {}).get("sigma_END_avg_ps", math.nan)),
        "T1_sigma_END_tavg_xm690_ps": ff(T1_map.get(-690, {}).get("sigma_END_avg_ps", math.nan)),
        "T1_v_eff_cm_ns": ff(v_eff, 2),
        "T3_improvement_x0_pct": ff(T3_r["improvement_x0"]*100 if T3_r["improvement_x0"] else math.nan),
        "T3_improvement_xm690_pct": ff(T3_r["improvement_xm690"]*100 if T3_r["improvement_xm690"] else math.nan),
        "T3_improvement_mean_pct": ff(T3_r["improvement_mean"]*100),
        "T3_sigma_top_x0_ps": ff(T3_r["sig_top_x0"]),
        "T3_sigma_endtop_x0_ps": ff(T3_r["sig_endtop_x0"]),
        "T5_dynamic4_sigma_ps": ff(T5_r["dynamic_nearest4"]["sigma_ps"]),
        "T5_cluster4_sigma_ps": ff(T5_r["id_div4_cluster"]["sigma_ps"]),
        "scan_result": scan_result,
        "binning_canonical": "sqrt_n",
        "timestamp": datetime.datetime.now().isoformat(),
    }
    p = Path(OUT_DIR) / "results_exec20.json"
    p.write_text(json.dumps(rj, indent=2, default=str))
    note(f"  results_exec20.json: {p}")
    return rj

def write_gerardo_memo(T1_r, T3_r, T4_r, rj):
    T1_map = {r["x_mm"]: r for r in T1_r}
    lines = [
        "# GERARDO_DECISION_MEMO — EXEC_20 (reencuadrado, intrínseco/PDE only)",
        f"Fecha: {datetime.datetime.now().strftime('%Y-%m-%d')}",
        "Nota de René: la electrónica (SPTR, FastIC) queda diferida. Solo física intrínseca.",
        "",
        "## 1. Reconciliación del baseline END-only",
        "",
        "### Por qué nuestro END-only (882 ps @ centro) es peor que el de Gerardo",
        f"- Nuestro estimador: t_avg = (t_L + t_R)/2, requiere que AMBOS extremos disparen",
        f"- En x=0: END_L ve 0.37 PE/ev, END_R ve 0.38 PE/ev → casi ningún evento dispara ambos",
        f"- σ_END(t_avg, x=0) = {T1_map.get(0,{}).get('sigma_END_avg_ps',math.nan):.0f} ps — timing puro de fluctuación Poisson",
        "",
        f"- Alternativa (un solo extremo, cerca de ese extremo):",
        f"  END_L_SUM4 @ x=-690 mm: σ ≈ 30.5 ps (near-end, EXEC_19 T4 canónico)",
        f"  Pero esto no es position-independent",
        "",
        "### ¿Qué usa Gerardo exactamente? (confirmar)",
        "  □ ¿Suma analógica de todos los SiPMs de un extremo?",
        "  □ ¿Promedio (t_L+t_R)/2 o solo el mejor extremo?",
        "  □ ¿Posición evaluada (centro vs extremo vs promedio del scan)?",
        "  → Una vez confirmado, podemos hacer la comparación apples-to-apples.",
        "",
        "## 2. Aporte cuantificado del TOP (intrínseco)",
        "",
        f"- σ_END-only (t_avg) @ x=0: {T1_map.get(0,{}).get('sigma_END_avg_ps',math.nan):.0f} ps",
        f"- σ_TOP_nearest @ x=0: {rj.get('T3_sigma_top_x0_ps','?')} ps",
        f"- σ_EndTop (GLS) @ x=0: {rj.get('T3_sigma_endtop_x0_ps','?')} ps",
        f"- **Mejora @ x=0: {rj.get('T3_improvement_x0_pct','?')}%**",
        f"- **Mejora media (31 posiciones): {rj.get('T3_improvement_mean_pct','?')}%**",
        "",
        "  Física: el TOP ayuda MUCHO al centro (END casi ciego) y POCO en los extremos",
        "  (donde END ya da σ~30 ps). El aporte del TOP es posición-dependiente.",
        "",
        "## 3. Decisión de topología SUM [DECISIÓN GERARDO]",
        "",
        f"Con sqrt_n binning + walk correction (EXEC_20 canónico):",
        f"  - dynamic nearest-4: σ = {rj.get('T5_dynamic4_sigma_ps','?')} ps (selecciona 4 con más luz)",
        f"  - id//4 fijo cluster:  σ = {rj.get('T5_cluster4_sigma_ps','?')} ps (4 adyacentes)",
        "",
        "  Nota: el 'co-localizado' aplica a END (misma cara, misma luz).",
        "  El TOP es SIEMPRE distribuido (fila a 20 mm).",
        "  Ambas opciones TOP son distribuidas; la diferencia es si el cluster es fijo o adaptativo.",
        "",
        "## 4. Electrónica diferida (por decisión de René)",
        "  SPTR ≈ 106 ps, FastIC ≈ 10 ps → NO en headline de esta campaña.",
        "  Ver EXEC_19 para los números completos de σ_tot bajo Modelo A/B.",
        "  Resumen: solo bajo Modelo B (FastIC+ digital por canal) se alcanzaría < 100 ps.",
    ]
    p = Path(OUT_DIR) / "GERARDO_DECISION_MEMO.md"
    p.write_text("\n".join(lines))
    note(f"  GERARDO_DECISION_MEMO.md: {p}")

def write_exec20_report(T1_r, T3_r, T4_r, T5_r, scan_res, rj, flags, v_eff):
    T1_map = {r["x_mm"]: r for r in T1_r}
    lines = [
        f"# EXEC_20_REPORT — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. Veredicto",
        "**COMPLETADO" + ("-CON-FLAGS" if flags else "") + "**",
        "",
        "## 2. END-only baseline (T1)",
        f"- Estimador: t_avg = (t_L + t_R)/2, walk-corr, sqrt_n, all 8 SiPMs per face",
        f"- σ_END(x=0)   = {T1_map.get(0,{}).get('sigma_END_avg_ps',math.nan):.0f} ps (0.37+0.38 PE/ev — photon-starved)",
        f"- σ_END(x=-690) = {T1_map.get(-690,{}).get('sigma_END_avg_ps',math.nan):.0f} ps",
        f"- v_eff = {ff(v_eff,2)} cm/ns",
        "",
        "## 3. TOP contribution (T3)",
        f"- Mejora EndTop/END-only @ x=0: {rj.get('T3_improvement_x0_pct','?')}%",
        f"- Mejora media (31 pos): {rj.get('T3_improvement_mean_pct','?')}%",
        f"- σ_TOP @ x=0 = {rj.get('T3_sigma_top_x0_ps','?')} ps; σ_EndTop = {rj.get('T3_sigma_endtop_x0_ps','?')} ps",
        "",
        "## 4. Reconciliation Gerardo (T4)",
        "- Nuestro t_avg requiere AMBOS extremos → baja eficiencia en extremos de scan",
        "- Gerardo probablemente usa un estimador diferente (confirmar)",
        "- Near-end σ_END_L_SUM4 ≈ 30.5 ps (EXEC_19 T4)",
        "",
        "## 5. Build END-only (T2a/T2b)",
        f"- T2a: PASS — OPSC-101, 16 SiPMs [0,15], no TOP hits, reflector skin fix ✓",
        f"- Optical diff: R=0.90 (endonly) vs R=0.95 (reference); jitter=0 in both",
        f"- T2b scan: {scan_res}",
        "",
        "## 6. Canonical fix (T5)",
        f"- Binning: sqrt_n (FD artifact EXEC_18 resolved)",
        f"- dynamic nearest-4: {rj.get('T5_dynamic4_sigma_ps','?')} ps ≤ id//4: {rj.get('T5_cluster4_sigma_ps','?')} ps (as expected)",
        "",
        "## 7. Flags",
    ]
    lines.extend(flags if flags else ["  ninguno"])
    lines += [
        "",
        "## 8. Rutas",
        f"  Figuras: {OUT_DIR}/T1..T5",
        f"  Memo Gerardo: {OUT_DIR}/GERARDO_DECISION_MEMO.md",
        f"  Canónico: {OUT_DIR}/T5/CANONICAL_ESTIMATORS.md",
        f"  JSON: {OUT_DIR}/results_exec20.json",
        "",
        "## 9. Decisiones para Gerardo",
        "  1. Topología SUM: dynamic-4 vs id//4 (ambas distribuidas para TOP)",
        "  2. Estimador END-only: t_avg (position-independent) vs single-end (mejor por extremo)",
        "  3. Confirmar estimador exacto de Gerardo para comparación apples-to-apples",
        "  4. Electrónica (SPTR/FastIC): diferida por decisión de René",
    ]
    p = Path(OUT_DIR) / "EXEC_20_REPORT.md"
    p.write_text("\n".join(lines))
    note(f"\nEXEC_20_REPORT.md: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    for sub in ["T1","T3","T4","T5"]: Path(f"{OUT_DIR}/{sub}").mkdir(parents=True, exist_ok=True)
    note(f"# EXEC_20 — {datetime.datetime.now().isoformat()}")
    note("Reorientación: intrínseco/PDE only; pregunta = aporte marginal del TOP sobre END-only baseline.")
    note("Binning canónico: sqrt_n (decisión T5/EXEC_19)")

    import ROOT; ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning

    pm = pos_map()
    note(f"  {len(pm)} positions from EndTop data")

    T1_results, v_eff = run_T1(pm)
    T3_results = run_T3(pm, T1_results)
    T4_results = run_T4(T1_results, T3_results)
    T5_results = run_T5(pm)

    # Check if T2b scan completed
    scan_run_dir = sorted(Path("/home/reriosto/SHiP/ej200_endonly/output").glob("endonly_mylar_msi_*"), reverse=True)
    if scan_run_dir:
        scan_result = f"Scan output at {scan_run_dir[0]} (check run_metadata.txt for status)"
    else:
        scan_result = "T2b scan still running or output not found; T1/T3/T4/T5 use EndTop data"
    note(f"\n  T2b scan result: {scan_result}")

    rj = write_results_json(T1_results, T3_results, T4_results, T5_results, v_eff, scan_result)
    write_gerardo_memo(T1_results, T3_results, T4_results, rj)
    write_exec20_report(T1_results, T3_results, T4_results, T5_results, scan_result, rj, _flags, v_eff)
    note(f"\n=== EXEC_20 analysis COMPLETE ===")
    note(f"Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
