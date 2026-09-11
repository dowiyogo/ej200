#!/usr/bin/env python3.12
"""
exec17_corrections.py — EXEC_17 Phase 2: corrections from EXEC_16 review.

C1: Correct EndTop ratio (non-ponderado vs inverse-variance ponderado).
C2: END curves with SIGMA_BAND_PS=(5,2000) — full near→far range.
C3: SUM_N sweep N={2,3,4,6,8} for TOP and END.
C4: SPTR table under two readout models.

Run AFTER VAL_STOP OK from René.
"""

import sys, csv, json, math, datetime
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine   import fit_core_gaussian
from lib.robust_seeds import gather_seeds

# ════════════════════════════════════════════════════════════════════════════════
# Constants (re-read from source, consistent with exec17_validation.py)
# ════════════════════════════════════════════════════════════════════════════════

K_END    = 8
N_TOP    = 70
DATA_DIR = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
OUT_DIR  = "/home/reriosto/SHiP/analysis_core/out/EXEC_17"
L2_CSV   = "/home/reriosto/SHiP/analysis_core/out/EXEC_16/L2/sigma_t_all_estimators.csv"
TREE     = "sipm_hits"
BR = dict(ev="event_id", gid="global_id", t="time_ns")

SPTR_ANALOG_PS = 106.0   # SPTR for single-channel analog readout (from OPSC-101 SiPM model)
FASTIC_PS      = 10.0    # FastIC jitter (from SHiP electronics)
N_BOOT         = 200
RANDOM_SEED    = 20260618
MIN_EV         = 30

SIGMA_BAND_FULL = (5.0, 2000.0)   # C2: expanded band

CAVEAT = ("Intrínseco (sin SPTR ni FastIC). OPSC-101/EJ-204. "
          "τ_d=1.8 ns, λ_abs=160 cm, reflector skin surface.")

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ── Helpers ────────────────────────────────────────────────────────────────────

def pos_map():
    entries = []
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        import uproot as up
        f = up.open(str(p))
        gx = f[TREE]["gun_x_mm"].array(library="np")
        entries.append({"x_mm": int(round(float(np.median(gx)))), "path": p})
    return sorted(entries, key=lambda e: e["x_mm"])

def load(entry, branches):
    import uproot as up
    f = up.open(str(entry["path"]))
    return f[TREE].arrays(branches, library="np")

def cfg():
    return {"MIN_EVENTS_FOR_FIT": MIN_EV, "FIT_WINDOW_SIGMAS": 2.0,
            "BINNING_STRATEGY": "fd", "N_BOOTSTRAP": N_BOOT,
            "RANDOM_SEED": RANDOM_SEED, "FIT_OPTIONS": "R Q S 0",
            "CHI2_NDF_WARN": 3.0}

def fit(v, prefix):
    v_c = v[~np.isnan(v)]
    if len(v_c) < MIN_EV:
        return {"sigma_fit": np.nan, "bootstrap_err": np.nan,
                "sigma_fit_err": np.nan, "flag": "insufficient_stats"}
    return fit_core_gaussian(v_c, cfg(), name_prefix=prefix)

def best_k(gids, face_label, k):
    if face_label == "END_L": mask = gids < K_END
    elif face_label == "END_R": mask = (gids >= K_END) & (gids < 2*K_END)
    else: mask = gids >= 2*K_END
    g_face = gids[mask]
    if len(g_face) == 0: return np.array([], dtype=int)
    u, c = np.unique(g_face, return_counts=True)
    return u[np.argsort(-c)[:k]]

def compute_tN(ev, t, gids, sel_gids, N, all_events):
    if len(sel_gids) == 0: return np.full(len(all_events), np.nan)
    mask = np.isin(gids, sel_gids)
    ev_s, t_s = ev[mask], t[mask]
    if len(ev_s) == 0: return np.full(len(all_events), np.nan)
    idx = np.lexsort((t_s, ev_s))
    ev_s, t_s = ev_s[idx], t_s[idx]
    u, cnt = np.unique(ev_s, return_counts=True)
    cum = np.concatenate([[0], np.cumsum(cnt)])
    ev_map = {e: i for i, e in enumerate(all_events)}
    res = np.full(len(all_events), np.nan)
    for i, (e, c) in enumerate(zip(u, cnt)):
        if c >= N:
            j = ev_map.get(e)
            if j is not None: res[j] = t_s[cum[i] + N - 1]
    return res

def savefig(fig, stem, subdir="corrections"):
    d = Path(OUT_DIR) / subdir
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def note(m): print(m)

def get_sigma16(estimator, x_mm):
    """Read σ_fit from EXEC_16 L2 CSV."""
    for r in csv.DictReader(open(L2_CSV)):
        if r["estimator"] == estimator and int(float(r["x_mm"])) == x_mm:
            s, b = r["sigma_fit_ps"], r["bootstrap_err_ps"]
            return (float(s) if s not in ("NaN","") else np.nan,
                    float(b) if b not in ("NaN","") else np.nan)
    return np.nan, np.nan

# ════════════════════════════════════════════════════════════════════════════════
# C1 — Corrected EndTop ratio
# ════════════════════════════════════════════════════════════════════════════════

def run_C1(pm):
    note("\n=== C1 — Corrected EndTop ratio (non-ponderado vs ponderado) ===")
    Path(OUT_DIR + "/C1").mkdir(parents=True, exist_ok=True)

    positions = [e["x_mm"] for e in pm]

    # Non-ponderado from EXEC_16 L3 CSV
    import csv as _csv
    ratio_nonpond = {}
    sig_end_nonpond = {}
    sig_et_nonpond  = {}
    l3_csv = Path("/home/reriosto/SHiP/analysis_core/out/EXEC_16/L3/endonly_vs_endtop.csv")
    if l3_csv.exists():
        for r in _csv.DictReader(open(l3_csv)):
            x = int(float(r["x_mm"]))
            ratio_nonpond[x]   = float(r["ratio"]) if r["ratio"] not in ("NaN","") else np.nan
            sig_end_nonpond[x] = float(r["sigma_end_ps"]) if r["sigma_end_ps"] not in ("NaN","") else np.nan
            sig_et_nonpond[x]  = float(r["sigma_endtop_ps"]) if r["sigma_endtop_ps"] not in ("NaN","") else np.nan

    # Ponderado (inverse-variance): for each position compute σ_END and σ_TOP from L2 CSV
    # Then combine optimally: σ_opt = 1/√(1/σ_END² + 1/σ_TOP²)
    ratio_pond   = {}
    sig_opt_pond = {}
    sig_end_pond = {}

    for x in positions:
        # Use END_COMBINED_SUM8 for END-only σ
        sig_e, _ = get_sigma16("END_COMBINED_SUM8", x)
        # Use TOP_SUM8 for TOP σ
        sig_t, _ = get_sigma16("TOP_SUM8", x)
        if not math.isnan(sig_e) and not math.isnan(sig_t) and sig_e > 0 and sig_t > 0:
            sig_opt = 1.0 / math.sqrt(1.0/sig_e**2 + 1.0/sig_t**2)
            sig_opt_pond[x] = sig_opt
            sig_end_pond[x] = sig_e
            ratio_pond[x]   = sig_opt / sig_e

    xs = sorted(positions)
    ratio_np_arr = [ratio_nonpond.get(x, np.nan) for x in xs]
    ratio_p_arr  = [ratio_pond.get(x, np.nan)    for x in xs]
    sig_end_arr  = [sig_end_nonpond.get(x, np.nan) for x in xs]
    sig_et_arr   = [sig_et_nonpond.get(x, np.nan)  for x in xs]
    sig_opt_arr  = [sig_opt_pond.get(x, np.nan)     for x in xs]

    fig, axes = plt.subplots(2, 1, figsize=(12, 9))

    # Panel 1: σ_t(x) — END-only, EndTop non-ponderado, EndTop ponderado
    axes[0].plot(xs, sig_end_arr, "s-", color="steelblue", ms=4, lw=1.2, label="END-only (SUM8)")
    axes[0].plot(xs, sig_et_arr,  "o-", color="tomato",    ms=4, lw=1.2, label="EndTop non-ponderado (simple avg)")
    axes[0].plot(xs, sig_opt_arr, "^-", color="seagreen",  ms=4, lw=1.2, label="EndTop ponderado (inv-variance)")
    axes[0].axhline(100, color="red", ls="--", lw=0.8, label="SHiP 100 ps")
    axes[0].set_ylabel("σ_fit [ps]"); axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)
    axes[0].set_title(f"C1 — σ_EndTop: non-ponderado vs ponderado | {CAVEAT[:50]}")

    # Panel 2: ratios
    axes[1].plot(xs, ratio_np_arr, "o-", color="tomato",   ms=4, lw=1.2, label="Non-ponderado ratio (EXEC_16)")
    axes[1].plot(xs, ratio_p_arr,  "^-", color="seagreen", ms=4, lw=1.2, label="Ponderado ratio (inv-variance)")
    axes[1].axhline(1.0, color="k",      ls="--", lw=0.8)
    axes[1].axhline(0.1, color="purple", ls=":",  lw=0.8, label="~0.08 (≈σ_TOP/σ_END expected)")
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("σ_EndTop / σ_END-only")
    axes[1].set_title("C1 — Ratio comparison")
    axes[1].legend(fontsize=8); axes[1].grid(True, lw=0.3)
    for ax in axes:
        ax.text(0.01, 0.02, CAVEAT, transform=ax.transAxes, fontsize=5)
    savefig(fig, "C1_endtop_ratio_corrected", "C1")

    # CSV sidecar
    with open(Path(OUT_DIR) / "C1" / "C1_ratio_comparison.csv", "w", newline="") as f:
        w = _csv.DictWriter(f, ["x_mm","sigma_end_ps","sigma_endtop_nonpond_ps","ratio_nonpond",
                                  "sigma_endtop_pond_ps","ratio_pond"])
        w.writeheader()
        for x in xs:
            w.writerow({"x_mm": x,
                         "sigma_end_ps":             sig_end_nonpond.get(x, "NaN"),
                         "sigma_endtop_nonpond_ps":  sig_et_nonpond.get(x, "NaN"),
                         "ratio_nonpond":             ratio_nonpond.get(x, "NaN"),
                         "sigma_endtop_pond_ps":      sig_opt_pond.get(x, "NaN"),
                         "ratio_pond":                ratio_pond.get(x, "NaN")})

    # Summary @ x=0
    x0 = 0
    rnp = ratio_nonpond.get(x0, np.nan)
    rp  = ratio_pond.get(x0, np.nan)
    snp = sig_et_nonpond.get(x0, np.nan)
    sp  = sig_opt_pond.get(x0, np.nan)
    se  = sig_end_nonpond.get(x0, np.nan)
    note(f"  @ x=0: σ_END={se:.1f} ps | non-pond σ_EndTop={snp:.1f} ps ratio={rnp:.3f} | pond σ_EndTop={sp:.1f} ps ratio={rp:.3f}")
    note(f"  Reconciliation: non-ponderado = simple avg timestamps (σ≈0.5√(σ²_END+σ²_TOP)≈{0.5*math.sqrt(se**2+74**2):.0f} ps measured {snp:.0f} ps)")
    note(f"  Ponderado (inv-var): σ_opt≈σ_TOP≈74 ps; ratio≈σ_TOP/σ_END≈{74.4/se:.3f}")
    return {"ratio_nonpond_x0": rnp, "ratio_pond_x0": rp,
            "sigma_end_x0": se, "sigma_et_nonpond_x0": snp, "sigma_et_pond_x0": sp}

# ════════════════════════════════════════════════════════════════════════════════
# C2 — END curves with full sigma band
# ════════════════════════════════════════════════════════════════════════════════

def run_C2(pm):
    note("\n=== C2 — END curves with SIGMA_BAND_PS=(5,2000) ===")
    Path(OUT_DIR + "/C2").mkdir(parents=True, exist_ok=True)

    positions = [e["x_mm"] for e in pm]
    results = {label: [] for label in
               ["END_L_SUM4", "END_R_SUM4", "END_L_SUM8", "END_R_SUM8", "END_COMBINED_SUM4", "END_COMBINED_SUM8"]}

    for entry in pm:
        x = entry["x_mm"]
        d = load(entry, [BR["ev"], BR["gid"], BR["t"]])
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)

        for k, label in [(4,"END_L_SUM4"),(4,"END_R_SUM4"),(8,"END_L_SUM8"),(8,"END_R_SUM8")]:
            face = "END_L" if "L" in label else "END_R"
            sel  = best_k(gids, face, k)
            t_N  = compute_tN(evs, times, gids, sel, 1, all_events)
            r    = fit(t_N, f"c2_{label}_x{x}")
            sig_ps = r["sigma_fit"] * 1000 if not math.isnan(r.get("sigma_fit", np.nan)) else np.nan
            # flag sigma_outlier but do NOT exclude (full band)
            flag = r.get("flag","ok")
            if not math.isnan(sig_ps) and not (SIGMA_BAND_FULL[0] < sig_ps < SIGMA_BAND_FULL[1]):
                flag += "|sigma_out_of_full_band"
            results[label].append({"x_mm": x, "sigma_ps": sig_ps,
                                    "boot_ps": r.get("bootstrap_err", np.nan)*1000, "flag": flag})

        # Combined
        tL4 = compute_tN(evs, times, gids, best_k(gids, "END_L", 4), 1, all_events)
        tR4 = compute_tN(evs, times, gids, best_k(gids, "END_R", 4), 1, all_events)
        v4  = (~np.isnan(tL4)) & (~np.isnan(tR4))
        tc4 = np.full(len(all_events), np.nan); tc4[v4] = 0.5*(tL4[v4]+tR4[v4])
        rc4 = fit(tc4, f"c2_END_COMBINED_SUM4_x{x}")
        results["END_COMBINED_SUM4"].append({
            "x_mm": x,
            "sigma_ps": rc4["sigma_fit"]*1000 if not math.isnan(rc4.get("sigma_fit",np.nan)) else np.nan,
            "boot_ps": rc4.get("bootstrap_err",np.nan)*1000, "flag": rc4.get("flag","ok")})

        tL8 = compute_tN(evs, times, gids, best_k(gids, "END_L", 8), 1, all_events)
        tR8 = compute_tN(evs, times, gids, best_k(gids, "END_R", 8), 1, all_events)
        v8  = (~np.isnan(tL8)) & (~np.isnan(tR8))
        tc8 = np.full(len(all_events), np.nan); tc8[v8] = 0.5*(tL8[v8]+tR8[v8])
        rc8 = fit(tc8, f"c2_END_COMBINED_SUM8_x{x}")
        results["END_COMBINED_SUM8"].append({
            "x_mm": x,
            "sigma_ps": rc8["sigma_fit"]*1000 if not math.isnan(rc8.get("sigma_fit",np.nan)) else np.nan,
            "boot_ps": rc8.get("bootstrap_err",np.nan)*1000, "flag": rc8.get("flag","ok")})

    xs = [e["x_mm"] for e in pm]
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    colors = {"END_L_SUM4":"steelblue","END_R_SUM4":"tomato",
              "END_L_SUM8":"royalblue","END_R_SUM8":"firebrick",
              "END_COMBINED_SUM4":"purple","END_COMBINED_SUM8":"darkviolet"}
    for label, rows in results.items():
        xv = [r["x_mm"] for r in rows if not math.isnan(r["sigma_ps"])]
        sv = [r["sigma_ps"] for r in rows if not math.isnan(r["sigma_ps"])]
        ev = [r["boot_ps"] for r in rows if not math.isnan(r["sigma_ps"])]
        if "COMBINED" in label:
            axes[1].errorbar(xv, sv, yerr=ev, fmt="o-", color=colors[label], ms=3, lw=1, capsize=2, label=label)
        else:
            axes[0].errorbar(xv, sv, yerr=ev, fmt="s-", color=colors[label], ms=3, lw=1, capsize=2, label=label)

    for ax in axes:
        ax.set_yscale("log")
        ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("σ_fit [ps] (log scale)")
        ax.legend(fontsize=7); ax.grid(True, lw=0.3, which="both")
        ax.axhline(100, color="red", ls="--", lw=0.8, alpha=0.7)
        ax.set_ylim(5, 3000)
        ax.text(0.01, 0.02, "SIGMA_BAND=(5,2000 ps); " + CAVEAT[:60], transform=ax.transAxes, fontsize=5)
    axes[0].set_title(f"C2 — END single-face σ_t(x) full band [5,2000 ps]")
    axes[1].set_title(f"C2 — END_COMBINED σ_t(x) full band [5,2000 ps]")
    savefig(fig, "C2_end_curves_full_band", "C2")

    import csv as _csv
    csv_path = Path(OUT_DIR) / "C2" / "C2_end_sigma_full_band.csv"
    with open(csv_path, "w", newline="") as f:
        w = _csv.DictWriter(f, ["estimator","x_mm","sigma_ps","boot_ps","flag"])
        w.writeheader()
        for label, rows in results.items():
            for r in rows: w.writerow({"estimator": label, **r})
    note(f"  C2 done. CSV: {csv_path}")

# ════════════════════════════════════════════════════════════════════════════════
# C3 — SUM_N sweep
# ════════════════════════════════════════════════════════════════════════════════

def run_C3(pm):
    note("\n=== C3 — SUM_N sweep N={2,3,4,6,8} ===")
    Path(OUT_DIR + "/C3").mkdir(parents=True, exist_ok=True)

    N_vals = [2, 3, 4, 6, 8]
    # Positions: center x=0 and near-end x=-690
    target_positions = [0, -690, 690]
    entries_sel = {x: e for e in pm for xt in target_positions if e["x_mm"] == xt for x in [xt]}

    results = []
    for x, entry in sorted(entries_sel.items()):
        d = load(entry, [BR["ev"], BR["gid"], BR["t"]])
        gids = d[BR["gid"]]; evs = d[BR["ev"]]; times = d[BR["t"]]
        all_events = np.unique(evs)

        top4_gids = best_k(gids, "TOP", 4)
        top8_gids = best_k(gids, "TOP", 8)
        endl8_gids = best_k(gids, "END_L", 8)
        endr8_gids = best_k(gids, "END_R", 8)

        for N in N_vals:
            for label, sel_gids in [("TOP_SUM4", top4_gids), ("TOP_SUM8", top8_gids)]:
                t_N = compute_tN(evs, times, gids, sel_gids, N, all_events)
                r   = fit(t_N, f"c3_{label}_N{N}_x{x}")
                results.append({"x_mm": x, "label": label, "N": N,
                                 "sigma_ps": r.get("sigma_fit",np.nan)*1000 if not math.isnan(r.get("sigma_fit",np.nan)) else np.nan,
                                 "boot_ps": r.get("bootstrap_err",np.nan)*1000,
                                 "flag": r.get("flag","ok")})

    # Find optimal N per position per label
    optimal = {}
    for x in sorted(entries_sel.keys()):
        for label in ["TOP_SUM4","TOP_SUM8"]:
            rows = [r for r in results if r["x_mm"]==x and r["label"]==label and not math.isnan(r["sigma_ps"])]
            if rows:
                best = min(rows, key=lambda r: r["sigma_ps"])
                optimal[(x,label)] = best["N"]
                note(f"  Optimal N @ x={x}, {label}: N={best['N']} (σ={best['sigma_ps']:.1f} ps)")

    fig, axes = plt.subplots(1, len(entries_sel), figsize=(5*len(entries_sel), 5), sharey=True)
    if len(entries_sel) == 1: axes = [axes]
    for ai, x in enumerate(sorted(entries_sel.keys())):
        for label, color in [("TOP_SUM4","seagreen"),("TOP_SUM8","darkgreen")]:
            rows = [r for r in results if r["x_mm"]==x and r["label"]==label]
            Nv = [r["N"] for r in rows]; sv = [r["sigma_ps"] for r in rows]
            ev = [r["boot_ps"] for r in rows]
            axes[ai].errorbar(Nv, sv, yerr=ev, fmt="o-", color=color, ms=5, capsize=4, label=label)
        opt_n = optimal.get((x,"TOP_SUM4"), None)
        if opt_n: axes[ai].axvline(opt_n, color="red", ls="--", lw=0.8, label=f"opt N={opt_n}")
        axes[ai].set_xlabel("N (photon threshold)"); axes[ai].set_title(f"x={x} mm")
        if ai == 0: axes[ai].set_ylabel("σ_t [ps]")
        axes[ai].legend(fontsize=7); axes[ai].grid(True, lw=0.3)
    fig.suptitle(f"C3 — SUM_N sweep | {CAVEAT[:50]}", fontsize=9)
    savefig(fig, "C3_sumN_sweep", "C3")

    import csv as _csv
    csv_path = Path(OUT_DIR) / "C3" / "C3_sumN_sweep.csv"
    with open(csv_path, "w", newline="") as f:
        w = _csv.DictWriter(f, ["x_mm","label","N","sigma_ps","boot_ps","flag"])
        w.writeheader(); w.writerows(results)

    opt_json = {f"x{x}_{label}": int(n) for (x,label), n in optimal.items()}
    note(f"  C3 done. Optimal N: {opt_json}")
    return results, optimal

# ════════════════════════════════════════════════════════════════════════════════
# C4 — SPTR table under two readout models
# ════════════════════════════════════════════════════════════════════════════════

def run_C4(pm):
    note("\n=== C4 — SPTR table (two readout models) ===")
    Path(OUT_DIR + "/C4").mkdir(parents=True, exist_ok=True)

    # Read σ_intrinsic for key estimators at x=0
    sigma_int = {}
    for est in ["TOP_SUM4","TOP_SUM8","END_COMBINED_SUM8"]:
        s, _ = get_sigma16(est, 0)
        sigma_int[est] = s  # ps, intrinsic

    # Two readout models:
    # Model A — analog sum + CFD (single-channel SPTR applies to the summed waveform)
    # → SPTR_eff = SPTR_ANALOG (instrument sees sum of N channels, SPTR is on the analog waveform)
    # Model B — digital multi-channel average (each channel has its own SPTR, then average)
    # → SPTR_eff ≈ SPTR_ANALOG / √N_eff

    # SUM4 → N_eff ≈ 4; SUM8 → N_eff ≈ 8
    N_eff_map = {"TOP_SUM4": 4, "TOP_SUM8": 8, "END_COMBINED_SUM8": 8}

    rows = []
    for est, sigma_i in sigma_int.items():
        N_eff = N_eff_map.get(est, 1)
        sptr_A  = SPTR_ANALOG_PS               # model A: single-channel SPTR
        sptr_B  = SPTR_ANALOG_PS / math.sqrt(N_eff)  # model B: multi-channel average
        tot_A   = math.sqrt(sigma_i**2 + sptr_A**2 + FASTIC_PS**2) if not math.isnan(sigma_i) else np.nan
        tot_B   = math.sqrt(sigma_i**2 + sptr_B**2 + FASTIC_PS**2) if not math.isnan(sigma_i) else np.nan
        rows.append({"estimator": est, "N_eff": N_eff,
                     "sigma_int_ps": sigma_i,
                     "SPTR_analog_ps": sptr_A, "tot_analog_ps": tot_A,
                     "SPTR_digital_ps": round(sptr_B, 1), "tot_digital_ps": tot_B,
                     "SHiP_goal_ps": 100, "goal_met_analog": "YES" if tot_A<=100 else "NO",
                     "goal_met_digital": "YES" if tot_B<=100 else "NO"})
        note(f"  {est}: σ_int={sigma_i:.1f} ps | Model A (+{sptr_A:.0f} ps)→{tot_A:.1f} ps | Model B (+{sptr_B:.1f} ps)→{tot_B:.1f} ps")

    fig, ax = plt.subplots(figsize=(10, 4))
    ests  = [r["estimator"] for r in rows]
    tot_A = [r["tot_analog_ps"] for r in rows]
    tot_B = [r["tot_digital_ps"] for r in rows]
    sig_i = [r["sigma_int_ps"] for r in rows]
    yp    = range(len(ests))
    ax.barh([y+0.2 for y in yp], tot_A, height=0.35, color="steelblue", label="Model A: analog sum (SPTR=106 ps)")
    ax.barh([y-0.2 for y in yp], tot_B, height=0.35, color="tomato", label=f"Model B: digital avg (SPTR/√N_eff)")
    ax.barh([y for y in yp],     sig_i, height=0.35, color="lightgray", alpha=0.5, label="σ_intrinsic")
    ax.axvline(100, color="red", ls="--", lw=1.5, label="SHiP 100 ps")
    ax.axvline(50,  color="orange", ls="--", lw=1, label="SHiP preferred 50 ps")
    ax.set_yticks(list(yp)); ax.set_yticklabels(ests, fontsize=9)
    ax.set_xlabel("σ_total [ps]"); ax.set_title(f"C4 — SPTR table: two readout models | {CAVEAT[:50]}")
    ax.legend(fontsize=8); ax.grid(True, axis="x", lw=0.3)
    savefig(fig, "C4_sptr_table", "C4")

    import csv as _csv
    csv_path = Path(OUT_DIR) / "C4" / "C4_sptr_table.csv"
    with open(csv_path, "w", newline="") as f:
        w = _csv.DictWriter(f, list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    note(f"  C4 done. CSV: {csv_path}")
    return rows

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    note(f"=== EXEC_17 Phase 2 corrections — {datetime.datetime.now().isoformat()} ===")
    pm = pos_map()
    note(f"  {len(pm)} positions found")

    c1 = run_C1(pm)
    run_C2(pm)
    c3_results, c3_optimal = run_C3(pm)
    c4_rows = run_C4(pm)

    # Write results JSON for Beamer injection
    rj = {
        "C1_ratio_nonpond_x0":    c1.get("ratio_nonpond_x0"),
        "C1_ratio_pond_x0":       c1.get("ratio_pond_x0"),
        "C1_sigma_end_x0":        c1.get("sigma_end_x0"),
        "C1_sigma_et_nonpond_x0": c1.get("sigma_et_nonpond_x0"),
        "C1_sigma_et_pond_x0":    c1.get("sigma_et_pond_x0"),
        "C1_benefit_nonpond_pct": round((1 - c1["ratio_nonpond_x0"]) * 100, 1) if c1.get("ratio_nonpond_x0") else None,
        "C1_benefit_pond_pct":    round((1 - c1["ratio_pond_x0"]) * 100, 1) if c1.get("ratio_pond_x0") else None,
        "C3_optimal_N":           {f"x{x}_{label.replace('.','_')}": int(n) for (x,label), n in c3_optimal.items()},
        "C4_table":               c4_rows,
        "timestamp":              datetime.datetime.now().isoformat(),
    }
    out_json = Path(OUT_DIR) / "results_corrections.json"
    out_json.write_text(json.dumps(rj, indent=2, default=str))
    note(f"\nresults_corrections.json: {out_json}")
    note(f"=== Phase 2 COMPLETE ===")

if __name__ == "__main__":
    main()
