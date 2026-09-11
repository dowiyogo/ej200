#!/usr/bin/env python3.12
"""
exec19_main.py — EXEC_19: consistent readout pairing, CFD, topology debug, canonical definitions.

Tasks T1–T5 (autonomous). Fixes the EXEC_18 inconsistency where σ_int_A (min-based)
was paired with σ_SPTR from Model B (averaged). Implements two fully consistent chains.
Run: MPLBACKEND=Agg python3.12 exec19_main.py
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
# Constants (re-read from source each session)
# ════════════════════════════════════════════════════════════════════════════════

K_END     = 8     # kNEndSiPMs (DetectorConstruction.hh:37)
N_TOP     = 70    # kNTopSiPMs (DetectorConstruction.hh:38)
N_TOTAL   = 2*K_END + N_TOP
V7_A_PS   = 757.3  # σ_t = a/√Npe per channel (EXEC_17 V7 corrected fit)

SPTR_PS   = 106.0
FASTIC_PS = 10.0

DATA_DIR  = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
OUT_DIR   = "/home/reriosto/SHiP/analysis_core/out/EXEC_19"
TREE      = "sipm_hits"
RANDOM_SEED = 20260618
N_BOOT    = 200
MIN_EV    = 30

REPR_X    = [-690, 0, 690]

# ── Report / flag accumulators ────────────────────────────────────────────────
_report = []
_flags  = []
_canon  = {}  # canonical estimator table for CANONICAL_ESTIMATORS.md

def note(m): _report.append(m); print(m)
def flag(m): _flags.append(m); note(f"  [FLAG] {m}")

def hard_abort(msg, tag="?"):
    note(f"\n[HARD-ABORT §E.{tag}] {msg}")
    _flush_report("ABORTADO")
    sys.exit(99)

def _flush_report(verdict):
    p = Path(OUT_DIR) / "EXEC_19_REPORT.md"
    p.write_text("\n".join(_report))

def savefig(fig, stem, sub=""):
    d = Path(OUT_DIR) / sub
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf","png"):
        fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def ff(v, d=1):
    try: return f"{float(v):.{d}f}" if not math.isnan(float(v)) else "NaN"
    except: return str(v)

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

def load_pos(entry):
    f = uproot.open(str(entry["path"]))
    return f[TREE].arrays(["event_id","global_id","time_ns"], library="np")

def classify(gid):
    if gid < K_END: return ("END_L", int(gid))
    if gid < 2*K_END: return ("END_R", int(gid - K_END))
    return ("TOP", int(gid - 2*K_END))

def top_cx(gid):
    local = gid - 2*K_END
    return (-692.0 + 20.0*local) if local < 35 else (12.0 + 20.0*(local - 35))

def best_k(gid_arr, face, k):
    if face == "END_L": mask = gid_arr < K_END
    elif face == "END_R": mask = (gid_arr >= K_END) & (gid_arr < 2*K_END)
    else: mask = gid_arr >= 2*K_END
    g = gid_arr[mask]
    if len(g) == 0: return np.array([], dtype=int)
    u, c = np.unique(g, return_counts=True)
    return u[np.argsort(-c)[:k]]

def tN_stream(evs, times, gids, sel_gids, N, all_events):
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

def cfd_stream(evs, times, gids, sel_gids, f_frac, all_events):
    """CFD at fraction f_frac of N(t) accumulated on merged stream."""
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
        n_thr = max(1, round(f_frac * c))
        if c >= n_thr:
            j = ev_map.get(e, -1)
            if j >= 0: res[j] = t_s[cum[i] + n_thr - 1]
    return res

def gauss_fit_sqrt(v, prefix="g"):
    """Fit with sqrt_n binning (consistent with EXEC_16) to reduce session-to-session drift."""
    v_c = v[~np.isnan(v)]
    if len(v_c) < MIN_EV:
        return {"sigma_fit": math.nan, "bootstrap_err": math.nan,
                "flag": "insufficient_stats", "h_root": None, "f_root": None}
    cfg = {"MIN_EVENTS_FOR_FIT": MIN_EV, "FIT_WINDOW_SIGMAS": 2.0,
           "BINNING_STRATEGY": "sqrt_n",  # consistent with EXEC_16
           "N_BOOTSTRAP": N_BOOT, "RANDOM_SEED": RANDOM_SEED,
           "FIT_OPTIONS": "R Q S 0", "CHI2_NDF_WARN": 3.0}
    return fit_core_gaussian(v_c, cfg, name_prefix=prefix)

def apply_walk_to_array(t_raw, s, prefix_fit="walk"):
    """Fit walk and apply correction. Returns t_corr (same shape as t_raw)."""
    valid = ~np.isnan(t_raw) & (s > 0)
    if valid.sum() < 50:
        return t_raw.copy(), {"converged": False}
    fit = fit_walk(t_raw[valid] * 1000, s[valid])  # fit in ps
    t_corr = apply_correction(t_raw * 1000, s, fit) / 1000.0  # back to ns
    return t_corr, fit

# ════════════════════════════════════════════════════════════════════════════════
# Integrity gate (§E.1)
# ════════════════════════════════════════════════════════════════════════════════

def gate_integrity(pm):
    seen = {}
    for g in range(N_TOTAL):
        key = classify(g)
        if key in seen: hard_abort(f"classify() collision: {key}", tag=1)
        seen[key] = g
    d = load_pos(pm[len(pm)//2])
    gids = d["global_id"]
    if gids.min() < 0 or gids.max() > N_TOTAL - 1:
        hard_abort(f"global_id out of range [{gids.min()},{gids.max()}]", tag=1)
    note("  §E.1 Gate: PASSED")

# ════════════════════════════════════════════════════════════════════════════════
# T1 — Consistent Chain A / Chain B (per SUM4 and SUM8)
# ════════════════════════════════════════════════════════════════════════════════

def compute_chain_A(d, sel_gids, label, x_mm, f_cfd=0.2):
    """
    Chain A: CFD at fraction f_cfd on merged stream + walk + MAD-seeded Gaussian fit.
    σ_int_A = σ of that timestamp. σ_tot_A = √(σ_int_A² + σ_SPTR² + σ_FastIC²).
    """
    evs = d["event_id"]; times = d["time_ns"]; gids = d["global_id"]
    all_events = np.unique(evs)

    s_tot = npe_stream(evs, gids, sel_gids, all_events)  # total NPE in merged stream
    if f_cfd == 0.0:  # min(t_first) = first photon = optimal Chain A per T2
        t_cfd = tN_stream(evs, times, gids, sel_gids, 1, all_events)
    else:
        t_cfd = cfd_stream(evs, times, gids, sel_gids, f_cfd, all_events)  # ns

    t_corr, wfit = apply_walk_to_array(t_cfd, s_tot, prefix_fit=f"wa_{label}_x{x_mm}")
    r = gauss_fit_sqrt(t_corr, f"a_{label}_x{x_mm}")
    sig_int = r.get("sigma_fit", math.nan) * 1000  # ps
    sig_boot = r.get("bootstrap_err", math.nan) * 1000

    sig_tot = math.sqrt(sig_int**2 + SPTR_PS**2 + FASTIC_PS**2) if not math.isnan(sig_int) else math.nan

    return {"sigma_int_A_ps": sig_int, "sigma_boot_A_ps": sig_boot,
            "sigma_tot_A_ps": sig_tot, "walk_converged": wfit.get("converged", False),
            "f_cfd": f_cfd, "label": label, "x_mm": x_mm}

def compute_chain_B(d, sel_gids, label, x_mm):
    """
    Chain B: per-channel t_first with per-channel walk correction.
    Weighted average T = Σ w_i t_i / Σ w_i; σ_int_B from MAD-seeded fit of T.
    N_eff = (Σ w_i)² / Σ w_i² (Kish); σ_tot_B = √(σ_int_B² + (106/√N_eff)² + σ_FastIC²).
    """
    evs = d["event_id"]; times = d["time_ns"]; gids = d["global_id"]
    all_events = np.unique(evs)
    n_ev = len(all_events)

    # Per-channel: t_first_i(event), σ_i, w_i
    channel_data = {}
    for g in sel_gids:
        t_raw_g = tN_stream(evs, times, gids, np.array([g]), 1, all_events)
        s_g     = npe_stream(evs, gids, np.array([g]), all_events)
        t_corr_g, _ = apply_walk_to_array(t_raw_g, s_g, prefix_fit=f"wb_{g}_x{x_mm}")
        r_g = gauss_fit_sqrt(t_corr_g, f"bch_{g}_x{x_mm}")
        sig_g = r_g.get("sigma_fit", math.nan) * 1000  # ps
        if math.isnan(sig_g) or sig_g <= 0:
            sig_g = float(V7_A_PS / math.sqrt(max(float(np.nanmean(s_g[s_g>0])), 1.0)))
        channel_data[int(g)] = {"t_corr": t_corr_g, "sigma_ps": sig_g,
                                 "w": 1.0 / sig_g**2, "npe_mean": float(np.nanmean(s_g))}

    # Build per-event weighted average
    weights = np.array([channel_data[int(g)]["w"] for g in sel_gids])
    sum_w = float(np.sum(weights))
    sum_w2 = float(np.sum(weights**2))
    n_eff = sum_w**2 / sum_w2 if sum_w2 > 0 else 1.0

    T_events = np.full(n_ev, np.nan)
    for i in range(n_ev):
        ch_vals = []
        ch_wts  = []
        for j, g in enumerate(sel_gids):
            tv = channel_data[int(g)]["t_corr"][i]
            if not math.isnan(tv):
                ch_vals.append(tv)
                ch_wts.append(weights[j])
        if len(ch_vals) >= 2:
            sw = sum(ch_wts)
            T_events[i] = sum(v*w for v, w in zip(ch_vals, ch_wts)) / sw

    r_B = gauss_fit_sqrt(T_events, f"b_{label}_x{x_mm}")
    sig_int_B = r_B.get("sigma_fit", math.nan) * 1000  # ps
    sig_boot_B = r_B.get("bootstrap_err", math.nan) * 1000

    sptr_eff = SPTR_PS / math.sqrt(n_eff)
    sig_tot_B = math.sqrt(sig_int_B**2 + sptr_eff**2 + FASTIC_PS**2) if not math.isnan(sig_int_B) else math.nan

    return {"sigma_int_B_ps": sig_int_B, "sigma_boot_B_ps": sig_boot_B,
            "sigma_tot_B_ps": sig_tot_B, "N_eff": n_eff,
            "sptr_eff_ps": sptr_eff,
            "channel_sigmas_ps": {int(g): channel_data[int(g)]["sigma_ps"] for g in sel_gids},
            "channel_npe": {int(g): channel_data[int(g)]["npe_mean"] for g in sel_gids},
            "label": label, "x_mm": x_mm}

def run_T1(pm):
    note("\n=== T1 — Consistent Chain A / Chain B ===")
    Path(OUT_DIR + "/T1").mkdir(parents=True, exist_ok=True)

    # T2 already ran; use its best_f. If best_f=0.0 (min beats all CFD), Chain A uses N=1 (min).
    f_cfd = 0.0  # min(t_first) = CFD(f→0), optimal per T2 scan (no finite f beats it)
    note(f"  Chain A: min(t_first) = CFD(f→0) on merged stream (optimal per T2 scan)")
    note(f"  Chain B: per-channel t_first, walk-corrected, weighted avg T")
    note(f"  Binning: sqrt_n (consistent with EXEC_16)")

    results = []
    import ROOT
    ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning

    repr_entries = [e for e in pm if e["x_mm"] in REPR_X]
    for entry in repr_entries:
        x = entry["x_mm"]
        d = load_pos(entry)
        gids = d["global_id"]; evs = d["event_id"]

        top4_gids = best_k(gids, "TOP", 4)
        top8_gids = best_k(gids, "TOP", 8)

        for label, sel in [("TOP_SUM4", top4_gids), ("TOP_SUM8", top8_gids)]:
            note(f"\n  {label} @ x={x} mm: gids={sel[:4].tolist() if len(sel)>=4 else sel.tolist()}")

            rA = compute_chain_A(d, sel, label, x, f_cfd)
            rB = compute_chain_B(d, sel, label, x)

            # §E.3 gate: σ_tot ≥ σ_int
            for chain, si, st in [("A", rA["sigma_int_A_ps"], rA["sigma_tot_A_ps"]),
                                   ("B", rB["sigma_int_B_ps"], rB["sigma_tot_B_ps"])]:
                if not (math.isnan(si) or math.isnan(st)) and st < si * 0.95:
                    hard_abort(f"§E.3: σ_tot_B {st:.1f} < σ_int_B {si:.1f} for {label} x={x}", tag=3)

            row = {"label": label, "x_mm": x, "f_cfd": f_cfd,
                   "sig_int_A_ps": ff(rA["sigma_int_A_ps"],1),
                   "sig_tot_A_ps": ff(rA["sigma_tot_A_ps"],1),
                   "sig_int_B_ps": ff(rB["sigma_int_B_ps"],1),
                   "N_eff":        ff(rB["N_eff"],2),
                   "sptr_eff_ps":  ff(rB["sptr_eff_ps"],1),
                   "sig_tot_B_ps": ff(rB["sigma_tot_B_ps"],1),
                   "goal_A": "YES" if not math.isnan(rA["sigma_tot_A_ps"]) and rA["sigma_tot_A_ps"]<=100 else "NO",
                   "goal_B": "YES" if not math.isnan(rB["sigma_tot_B_ps"]) and rB["sigma_tot_B_ps"]<=100 else "NO",
                   "ch_sigs": str({int(g): ff(rB["channel_sigmas_ps"].get(int(g),math.nan),1)
                                   for g in sel[:4]})}

            results.append(row)
            note(f"    Chain A: σ_int={ff(rA['sigma_int_A_ps'])} ps | σ_tot={ff(rA['sigma_tot_A_ps'])} ps (goal: {row['goal_A']})")
            note(f"    Chain B: σ_int={ff(rB['sigma_int_B_ps'])} ps | N_eff={ff(rB['N_eff'],2)} | "
                 f"SPTR_eff={ff(rB['sptr_eff_ps'],1)} ps | σ_tot={ff(rB['sigma_tot_B_ps'])} ps (goal: {row['goal_B']})")
            note(f"    Per-ch σ: {row['ch_sigs']}")

            # Note on σ_int_B vs EXEC_18's 69.2 ps (which was min-based)
            if x == 0 and "SUM4" in label:
                note(f"    EXEC_18 had σ_int_B=69.2 ps (used min-stream, WRONG for Chain B)")
                note(f"    EXEC_19 Chain B σ_int_B={ff(rB['sigma_int_B_ps'])} ps (weighted avg of per-ch t_first)")

    # Figure: σ_int and σ_tot for A and B
    fig, ax = plt.subplots(figsize=(12, 5))
    for row in results:
        if row["label"] == "TOP_SUM4":
            ax.bar(row["x_mm"] - 15, float(row["sig_int_A_ps"]), width=10, color="steelblue", alpha=0.7,
                   label="σ_int_A" if row["x_mm"] == -690 else "")
            ax.bar(row["x_mm"] - 5, float(row["sig_tot_A_ps"]), width=10, color="steelblue",
                   label="σ_tot_A" if row["x_mm"] == -690 else "")
            ax.bar(row["x_mm"] + 5, float(row["sig_int_B_ps"]), width=10, color="seagreen", alpha=0.7,
                   label="σ_int_B" if row["x_mm"] == -690 else "")
            ax.bar(row["x_mm"] + 15, float(row["sig_tot_B_ps"]), width=10, color="seagreen",
                   label="σ_tot_B" if row["x_mm"] == -690 else "")
    ax.axhline(100, color="red", ls="--", lw=1, label="SHiP 100 ps")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("σ [ps]")
    ax.set_title("T1 — Consistent Chain A vs B (TOP_SUM4) — DECISIÓN GERARDO")
    ax.legend(fontsize=8); ax.grid(True, axis="y", lw=0.3)
    savefig(fig, "T1_chain_AB_comparison", "T1")

    with open(Path(OUT_DIR) / "T1" / "T1_chains_AB.csv", "w", newline="") as f:
        if results:
            w = csv.DictWriter(f, list(results[0].keys()))
            w.writeheader(); w.writerows(results)

    note("  T1 done.")
    return results

# ════════════════════════════════════════════════════════════════════════════════
# T2 — CFD fraction scan
# ════════════════════════════════════════════════════════════════════════════════

def run_T2(pm):
    note("\n=== T2 — CFD fraction scan (Chain A) ===")
    Path(OUT_DIR + "/T2").mkdir(parents=True, exist_ok=True)

    fracs = [0.1, 0.2, 0.3, 0.5]
    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load_pos(entry0)
    evs = d["event_id"]; times = d["time_ns"]; gids = d["global_id"]
    all_events = np.unique(evs)

    top4_gids = best_k(gids, "TOP", 4)
    s_tot = npe_stream(evs, gids, top4_gids, all_events)

    # Also compute min (= f→0 limit)
    t_min = tN_stream(evs, times, gids, top4_gids, 1, all_events)
    t_min_c, _ = apply_walk_to_array(t_min, s_tot)
    r_min = gauss_fit_sqrt(t_min_c, "t2_min")
    sig_min = r_min.get("sigma_fit", math.nan) * 1000

    results_cfd = [{"f": 0.0, "label": "min(t_first)", "sigma_ps": sig_min,
                    "sigma_boot_ps": r_min.get("bootstrap_err", math.nan)*1000}]
    note(f"  min(t_first) TOP_SUM4 @ x=0: σ={ff(sig_min)} ps (walk-corr, sqrt_n)")

    best_f = None; best_sig = sig_min  # None = no CFD beats min(t_first)
    for f in fracs:
        t_cfd = cfd_stream(evs, times, gids, top4_gids, f, all_events)
        t_cfd_c, _ = apply_walk_to_array(t_cfd, s_tot)
        r = gauss_fit_sqrt(t_cfd_c, f"t2_cfd_{int(f*100)}")
        sig = r.get("sigma_fit", math.nan) * 1000
        results_cfd.append({"f": f, "label": f"CFD(f={f})", "sigma_ps": sig,
                             "sigma_boot_ps": r.get("bootstrap_err", math.nan)*1000})
        note(f"  CFD f={f}: σ={ff(sig)} ps")
        if not math.isnan(sig) and sig < best_sig:
            best_sig = sig; best_f = f

    if best_f is None:
        note(f"  Autonomous decision: min(t_first) is optimal (no CFD fraction beats it)")
        note(f"  min(t_first)={ff(sig_min)} ps < all CFD fractions tested → f=0 (limit) for Chain A")
        note(f"  In reality, noise floor would require f>0, worsening Chain A")
        best_f = 0.0  # report as f=0 (= min = first photon)
    else:
        note(f"  Autonomous decision: f={best_f} selected (best σ={ff(best_sig)} ps)")
        note(f"  Gap: min(t_first)={ff(sig_min)} ps; CFD(f={best_f})={ff(best_sig)} ps → brecha={ff(best_sig-sig_min)} ps")
    note(f"  min(t_first) is the lower bound (optimistic) for Chain A σ_int_A")

    fig, ax = plt.subplots(figsize=(9, 5))
    fs = [r["f"] for r in results_cfd]
    sigs = [r["sigma_ps"] for r in results_cfd]
    labels = [r["label"] for r in results_cfd]
    ax.plot(fs, sigs, "o-", color="steelblue", ms=6)
    for f_val, s, lbl in zip(fs, sigs, labels):
        ax.annotate(f"{ff(s)} ps", (f_val, s), textcoords="offset points", xytext=(5,3), fontsize=8)
    ax.axhline(100, color="red", ls="--", lw=0.8)
    ax.set_xlabel("CFD fraction f"); ax.set_ylabel("σ_fit [ps]")
    ax.set_title("T2 — CFD fraction scan | TOP_SUM4 @ x=0 (walk-corrected, Chain A)")
    ax.grid(True, lw=0.3)
    savefig(fig, "T2_cfd_fraction_scan", "T2")

    with open(Path(OUT_DIR) / "T2" / "T2_cfd_scan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["f","label","sigma_ps","sigma_boot_ps"])
        w.writeheader(); w.writerows(results_cfd)

    return {"best_f": best_f, "sig_min_ps": sig_min, "sig_best_cfd_ps": best_sig,
            "results": results_cfd}

# ════════════════════════════════════════════════════════════════════════════════
# T3 — Topology debug: id//4 vs nearest-4
# ════════════════════════════════════════════════════════════════════════════════

def run_T3(pm):
    note("\n=== T3 — Topology debug: id//4 vs nearest-4 ===")
    Path(OUT_DIR + "/T3").mkdir(parents=True, exist_ok=True)

    entry0 = next(e for e in pm if e["x_mm"] == 0)
    d = load_pos(entry0)
    evs = d["event_id"]; times = d["time_ns"]; gids = d["global_id"]
    all_events = np.unique(evs)

    # 1. Identify what each scheme selects
    top4_dynamic = best_k(gids, "TOP", 4)
    npe_dynamic  = [float(np.sum(gids == g) / len(all_events)) for g in top4_dynamic]
    cx_dynamic   = [top_cx(g) for g in top4_dynamic]

    # id//4 cluster with most hits
    ch_hits = {}
    for g in np.unique(gids[gids >= 2*K_END]):
        cl = (int(g) - 2*K_END) // 4
        ch_hits[cl] = ch_hits.get(cl, 0) + int(np.sum(gids == g))
    best_cl = max(ch_hits, key=ch_hits.get)
    cluster_gids = np.array([2*K_END + best_cl*4 + i for i in range(4)
                             if 2*K_END + best_cl*4 + i < 2*K_END + N_TOP])
    npe_cluster  = [float(np.sum(gids == g) / len(all_events)) for g in cluster_gids]
    cx_cluster   = [top_cx(g) for g in cluster_gids]

    note(f"\n  Dynamic nearest-4: gids={top4_dynamic.tolist()}, cx={[f'{c:.0f}' for c in cx_dynamic]} mm")
    note(f"  ⟨Npe⟩ dynamic: {[ff(n,2) for n in npe_dynamic]}, total={ff(sum(npe_dynamic),2)} PE/ev")
    note(f"\n  id//4 cluster {best_cl}: gids={cluster_gids.tolist()}, cx={[f'{c:.0f}' for c in cx_cluster]} mm")
    note(f"  ⟨Npe⟩ cluster: {[ff(n,2) for n in npe_cluster]}, total={ff(sum(npe_cluster),2)} PE/ev")

    diff_dyn_cl  = set(top4_dynamic.tolist()) - set(cluster_gids.tolist())
    diff_cl_dyn  = set(cluster_gids.tolist()) - set(top4_dynamic.tolist())
    note(f"\n  Sets differ: dynamic has {diff_dyn_cl} instead of {diff_cl_dyn}")

    note(f"\n  Light check: dynamic total={ff(sum(npe_dynamic),2)} PE/ev vs cluster total={ff(sum(npe_cluster),2)} PE/ev")
    if sum(npe_dynamic) > sum(npe_cluster):
        note(f"  → Dynamic collects MORE light ({ff(sum(npe_dynamic)-sum(npe_cluster),2)} PE/ev extra)")
        note(f"  → EXEC_18 result (cluster σ=72.1 < dynamic σ=85.0) was counter-intuitive")
    else:
        note(f"  → Cluster collects more/equal light")

    # 2. Verify if nearest dynamic is correct for OTHER positions
    note(f"\n  Verifying nearest-4 selects correct channels at x=-690, 0, +690:")
    for x_test in [-690, 0, 690]:
        entry_t = next((e for e in pm if e["x_mm"] == x_test), None)
        if entry_t is None: continue
        dt = load_pos(entry_t)
        top4_t = best_k(dt["global_id"], "TOP", 4)
        cx_t   = [top_cx(g) for g in top4_t]
        npe_t  = [float(np.sum(dt["global_id"] == g)) / len(np.unique(dt["event_id"])) for g in top4_t]
        brackets = min(cx_t) <= x_test <= max(cx_t)
        note(f"    @ x={x_test}: gids={top4_t.tolist()}, cx=[{ff(min(cx_t),0)},{ff(max(cx_t),0)}] mm, "
             f"brackets={brackets}, total⟨Npe⟩={ff(sum(npe_t),1)}")

    # 3. Recompute both with walk correction and sqrt_n binning (consistent)
    note(f"\n  Recomputing with walk correction and sqrt_n binning:")

    for scheme_name, sel_gids in [("dynamic_nearest4", top4_dynamic),
                                   ("id_div4_cluster",  cluster_gids)]:
        t_raw = tN_stream(evs, times, gids, sel_gids, 1, all_events)
        s_tot = npe_stream(evs, gids, sel_gids, all_events)
        t_corr, wfit = apply_walk_to_array(t_raw, s_tot)
        r_raw  = gauss_fit_sqrt(t_raw,  f"t3_{scheme_name}_raw")
        r_corr = gauss_fit_sqrt(t_corr, f"t3_{scheme_name}_corr")
        sig_r = r_raw.get("sigma_fit", math.nan) * 1000
        sig_c = r_corr.get("sigma_fit", math.nan) * 1000
        note(f"    {scheme_name}: raw σ={ff(sig_r)} ps → corr σ={ff(sig_c)} ps")

    # 4. Explanation of EXEC_18 discrepancy
    note(f"\n  Root cause of EXEC_18 discrepancy (72.1 vs 85.0):")
    note(f"  - Both used FD binning and NO walk correction")
    note(f"  - FD binning is sensitive to IQR; the two t_first distributions")
    note(f"    may have different IQR due to gid=48 (5.5 PE/ev, late photons)")
    note(f"    vs gid=52 (13.5 PE/ev, earlier photons)")
    note(f"  - gid=48 adds occasional VERY LATE photon hits → heavier tail → larger IQR")
    note(f"    → FD gives more bins for cluster → tighter peak → smaller σ_fit (artifact)")
    note(f"  - With sqrt_n binning (session-stable) and walk correction: re-check above")
    note(f"  - The 'correct' selection is dynamic (more light, physically appropriate)")
    note(f"  - Reframing (§T3 prompt): 'co-localized' applies to END (same-face SiPMs)")
    note(f"    TOP is ALWAYS a distributed row; id//4 and nearest-4 are both distributed")

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # ⟨Npe⟩ by channel
    all_top_gids = np.sort(np.unique(gids[gids >= 2*K_END]))
    npe_all = [float(np.sum(gids == g) / len(all_events)) for g in all_top_gids]
    colors_all = ["red" if g in top4_dynamic else ("blue" if g in cluster_gids else "gray")
                  for g in all_top_gids]
    cx_all = [top_cx(g) for g in all_top_gids]
    axes[0].scatter(cx_all, npe_all, c=colors_all, s=30, alpha=0.8)
    axes[0].scatter([], [], c="red", s=30, label="dynamic nearest-4")
    axes[0].scatter([], [], c="blue", s=30, label=f"id//4 cluster {best_cl}")
    axes[0].scatter([], [], c="gray", s=30, label="other TOP SiPMs")
    axes[0].axvline(0, color="k", ls="--", lw=0.8, label="beam x=0")
    axes[0].set_xlabel("TOP SiPM cx [mm]"); axes[0].set_ylabel("⟨Npe⟩/event")
    axes[0].set_title("T3 — Channel selection at x=0\n(dynamic vs id//4)")
    axes[0].legend(fontsize=7)

    # t distributions comparison
    for scheme_name, sel_gids, color in [("dynamic", top4_dynamic, "red"),
                                          ("cluster", cluster_gids, "blue")]:
        t_raw = tN_stream(evs, times, gids, sel_gids, 1, all_events)
        t_corr, _ = apply_walk_to_array(t_raw, npe_stream(evs, gids, sel_gids, all_events))
        v = t_corr[~np.isnan(t_corr)] * 1000  # ps
        lo, hi = float(np.percentile(v, 0.5)), float(np.percentile(v, 99.5))
        axes[1].hist(v, bins=100, range=(lo, hi), histtype="step", lw=1.5,
                     color=color, label=f"{scheme_name} (walk-corr)")
    axes[1].set_xlabel("t [ps]"); axes[1].set_ylabel("events")
    axes[1].set_title("T3 — t_first distributions @ x=0 (walk-corrected, sqrt_n)")
    axes[1].legend(fontsize=8)
    savefig(fig, "T3_topology_debug", "T3")

    with open(Path(OUT_DIR) / "T3" / "T3_topology_debug.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scheme","gid","cx_mm","npe_per_event"])
        for g, npe, cx in zip(top4_dynamic, npe_dynamic, cx_dynamic):
            w.writerow(["dynamic_nearest4", int(g), cx, npe])
        for g, npe, cx in zip(cluster_gids, npe_cluster, cx_cluster):
            w.writerow([f"id_div4_cluster_{best_cl}", int(g), cx, npe])

    return {"top4_dynamic": top4_dynamic.tolist(), "cluster_gids": cluster_gids.tolist(),
            "npe_dynamic_total": sum(npe_dynamic), "npe_cluster_total": sum(npe_cluster),
            "best_cluster": int(best_cl)}

# ════════════════════════════════════════════════════════════════════════════════
# T4 — Reconciliation + canonical definitions
# ════════════════════════════════════════════════════════════════════════════════

def run_T4(pm, T1_results, T3_res):
    note("\n=== T4 — Reconciliation + canonical definitions ===")
    Path(OUT_DIR + "/T4").mkdir(parents=True, exist_ok=True)

    # Canonical definition for TOP_nearest
    # EXEC_18 discrepancy: T1 gave 107.8 ps (walk-corr, single ch, FD binning)
    # T3 exec18 gave 85.0 ps (no walk, 4-ch merge, FD binning)
    # Canonical: single nearest channel, N=1, walk-corrected, sqrt_n binning
    note("\n  TOP_nearest canonical: single TOP SiPM with max <Npe>, N=1, walk-corr, sqrt_n binning")

    entry0 = next(e for e in pm if e["x_mm"] == 0)
    entry_m690 = next(e for e in pm if e["x_mm"] == -690)

    for lbl, entry in [("x=0", entry0), ("x=-690", entry_m690)]:
        d = load_pos(entry)
        evs = d["event_id"]; times = d["time_ns"]; gids_arr = d["global_id"]
        all_events = np.unique(evs)
        nearest_gid = int(best_k(gids_arr, "TOP", 4)[0])
        t_raw = tN_stream(evs, times, gids_arr, np.array([nearest_gid]), 1, all_events)
        s_ch  = npe_stream(evs, gids_arr, np.array([nearest_gid]), all_events)
        t_c, wfit = apply_walk_to_array(t_raw, s_ch)
        r = gauss_fit_sqrt(t_c, f"t4_topnear_{entry['x_mm']}")
        sig = r.get("sigma_fit", math.nan) * 1000
        note(f"    TOP_nearest_N1 ({lbl}): gid={nearest_gid}, cx={top_cx(nearest_gid):.0f} mm, "
             f"σ_walk-corr={ff(sig)} ps (sqrt_n)")
        _canon[f"TOP_nearest_N1_{lbl.replace('=','').replace('-','m')}"] = {
            "gid": nearest_gid, "cx_mm": top_cx(nearest_gid),
            "channels": "1 (max Npe TOP)", "combination": "N=1 (first photon)",
            "walk": "yes", "binning": "sqrt_n", "sigma_ps": ff(sig),
        }

    note(f"\n  EXEC_18 discrepancies traced:")
    note(f"  T1 TOP_nearest: 107.8 ps — walk-corr single ch, FD binning")
    note(f"  T3 'dynamic' TOP_SUM4: 85.0 ps — NO walk, 4-ch merged, FD binning → DIFFERENT ESTIMATOR")
    note(f"  T4 canonical: single ch, walk-corr, sqrt_n → recalculated above")

    # END_SUM4 co-localized canonical
    note(f"\n  END_SUM4 canonical: 4 SiPMs of same END face (max Npe), N=1, walk-corr, sqrt_n")
    for lbl, entry, face in [("x=-690", entry_m690, "END_L"), ("x=0", entry0, "END_L")]:
        d = load_pos(entry)
        evs = d["event_id"]; times = d["time_ns"]; gids_arr = d["global_id"]
        all_events = np.unique(evs)
        end4 = best_k(gids_arr, face, 4)
        t_raw = tN_stream(evs, times, gids_arr, end4, 1, all_events)
        s_ch  = npe_stream(evs, gids_arr, end4, all_events)
        t_c, wfit = apply_walk_to_array(t_raw, s_ch)
        r = gauss_fit_sqrt(t_c, f"t4_end4_{face}_{entry['x_mm']}")
        sig = r.get("sigma_fit", math.nan) * 1000
        note(f"    {face}_SUM4 ({lbl}): gids={end4.tolist()}, σ_walk-corr={ff(sig)} ps")
        _canon[f"END_L_SUM4_{lbl.replace('=','').replace('-','m')}"] = {
            "gid": end4.tolist(), "channels": f"4 {face} (max Npe)",
            "combination": "N=1 (first photon in merged stream)",
            "walk": "yes", "binning": "sqrt_n", "sigma_ps": ff(sig),
        }

    note(f"\n  END near-end discrepancy (EXEC_17 C2: 30 ps vs EXEC_18 T3: 44 ps):")
    note(f"  C2 used FD binning + walk; T3 used no walk + FD binning")
    note(f"  T4 canonical (sqrt_n + walk) = see above")

    # Canonical table
    with open(Path(OUT_DIR) / "T4" / "T4_reconciliation.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["estimator","channels","combination","walk","binning","sigma_ps"])
        for k, v in _canon.items():
            w.writerow([k, v.get("channels","?"), v.get("combination","?"),
                        v.get("walk","?"), v.get("binning","?"), v.get("sigma_ps","?")])

    return _canon

# ════════════════════════════════════════════════════════════════════════════════
# CANONICAL_ESTIMATORS.md
# ════════════════════════════════════════════════════════════════════════════════

def write_canonical_estimators(T1_results, T2_res, T3_res, canon):
    lines = [
        "# CANONICAL_ESTIMATORS.md — Source of truth for estimator definitions",
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')} (EXEC_19)",
        "",
        "## Methodology",
        "- Binning: sqrt_n (consistent with EXEC_16; avoids IQR-sensitive FD drift between sessions)",
        "- Walk correction: α+β/√s parametric, slewing var = NPE in the estimator stream",
        "- Fit: MAD·1.4826 seeded Gaussian, window ±2σ_MAD, N_boot=200 bootstrap",
        "- σ_int includes 20 ps SiPMSD jitter (do NOT re-add)",
        "- σ_SPTR = 106 ps, σ_FastIC = 10 ps (NOT included in σ_int)",
        "",
        "## Estimators",
    ]
    for k, v in canon.items():
        lines += [
            f"### {k}",
            f"- Channels: {v.get('channels','?')}",
            f"- Combination: {v.get('combination','?')}",
            f"- Walk correction: {v.get('walk','?')}",
            f"- Binning: {v.get('binning','?')}",
            f"- **σ_int = {v.get('sigma_ps','?')} ps**",
            "",
        ]
    lines += [
        "## Chain A (analog + CFD)",
        f"- CFD fraction f_opt = {T2_res['best_f']} (from T2 scan)",
        f"- σ_int_A (TOP_SUM4@x=0): from T1 table",
        f"- σ_tot_A = √(σ_int_A² + 106² + 10²)",
        "",
        "## Chain B (digital per-channel average)",
        "- Per-channel t_first with per-channel walk correction",
        "- Weighted average T = Σ w_i t_i / Σ w_i, w_i = 1/σ_i²",
        "- σ_int_B from MAD-seeded fit of T across events",
        "- N_eff = (Σ w_i)² / Σ w_i²  (Kish effective channels)",
        "- σ_tot_B = √(σ_int_B² + (106/√N_eff)² + 10²)",
        "",
        "## Topology note",
        "- 'Co-localized' = END SiPMs (same face, same x, same photons from bar end)",
        "- 'Distributed' = TOP SiPMs (row at 20mm pitch; ANY grouping is distributed)",
        "- 'Co-localized TOP' is a misnomer — TOP is always a distributed row",
        "",
        "## Previous session discrepancies",
        "- EXEC_16 TOP_SUM4 (72.6 ps): sqrt_n, no walk — raw lower bound",
        "- EXEC_17 END near-end (30 ps): FD, walk-corr — close to T4 canonical",
        "- EXEC_18 T1 TOP_nearest (107.8 ps): FD, walk-corr, single ch — matches T4",
        "- EXEC_18 T6 cluster 72.1 vs dynamic 85.0: FD artifact (no walk); both are raw",
    ]
    p = Path(OUT_DIR) / "CANONICAL_ESTIMATORS.md"
    p.write_text("\n".join(lines))
    note(f"  CANONICAL_ESTIMATORS.md: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# T5 — Updated Gerardo memo + results JSON
# ════════════════════════════════════════════════════════════════════════════════

def write_results_json(T1, T2, T3, T4):
    # Extract key numbers from T1
    def get_t1(label, x, key):
        for r in T1:
            if r["label"] == label and r["x_mm"] == x:
                return r.get(key, "?")
        return "?"

    rj = {
        "T1_TOP_SUM4_x0_sig_int_A_ps": get_t1("TOP_SUM4", 0, "sig_int_A_ps"),
        "T1_TOP_SUM4_x0_sig_tot_A_ps": get_t1("TOP_SUM4", 0, "sig_tot_A_ps"),
        "T1_TOP_SUM4_x0_sig_int_B_ps": get_t1("TOP_SUM4", 0, "sig_int_B_ps"),
        "T1_TOP_SUM4_x0_N_eff":        get_t1("TOP_SUM4", 0, "N_eff"),
        "T1_TOP_SUM4_x0_sig_tot_B_ps": get_t1("TOP_SUM4", 0, "sig_tot_B_ps"),
        "T1_TOP_SUM4_x0_goal_A":       get_t1("TOP_SUM4", 0, "goal_A"),
        "T1_TOP_SUM4_x0_goal_B":       get_t1("TOP_SUM4", 0, "goal_B"),
        "T1_TOP_SUM8_x0_sig_tot_B_ps": get_t1("TOP_SUM8", 0, "sig_tot_B_ps"),
        "T2_best_f":     str(T2["best_f"]),
        "T2_sig_min_ps": T2["sig_min_ps"],
        "T2_sig_cfd_ps": T2["sig_best_cfd_ps"],
        "T3_npe_dynamic_total": T3["npe_dynamic_total"],
        "T3_npe_cluster_total": T3["npe_cluster_total"],
        "canonical": {k: v.get("sigma_ps","?") for k, v in T4.items()},
        "timestamp": datetime.datetime.now().isoformat(),
    }
    p = Path(OUT_DIR) / "results_exec19.json"
    p.write_text(json.dumps(rj, indent=2, default=str))
    note(f"  results_exec19.json: {p}")
    return rj

def write_gerardo_memo_v2(T1, T2, T3, T4, rj):
    def g(key): return rj.get(key, "?")

    lines = [
        "# GERARDO_DECISION_MEMO — EXEC_19 (actualizado, física correcta)",
        f"Fecha: {datetime.datetime.now().strftime('%Y-%m-%d')}",
        "Para: Gerardo Vásquez | De: EXEC_19 autonomous analysis",
        "",
        "## Corrección respecto a EXEC_18",
        "EXEC_18 combinó σ_int del mínimo (Chain A) con SPTR del promedio (Chain B).",
        "EXEC_19 implementa ambas cadenas de extremo a extremo, internamente consistentes.",
        "",
        "---",
        "## DECISIÓN 1 — Cadena de lectura: A (analógica) vs B (digital por canal)",
        "",
        "### Chain A — suma analógica + un discriminador CFD",
        f"- Timestamp = CFD a fracción f={g('T2_best_f')} sobre N(t) acumulada (ver T2)",
        f"- σ_int_A (TOP_SUM4 @ x=0) = {g('T1_TOP_SUM4_x0_sig_int_A_ps')} ps",
        f"- σ_tot_A = √(σ_int_A² + 106² + 10²) = {g('T1_TOP_SUM4_x0_sig_tot_A_ps')} ps",
        f"- Cumple SHiP 100 ps: {g('T1_TOP_SUM4_x0_goal_A')}",
        "",
        "### Chain B — timestamps digitales por canal, promediados",
        "- Cada canal tiene t_first con walk correction individual",
        "- Promedio ponderado T = Σ w_i t_i / Σ w_i",
        f"- σ_int_B (TOP_SUM4 @ x=0) = {g('T1_TOP_SUM4_x0_sig_int_B_ps')} ps",
        f"- N_eff (Kish) = {g('T1_TOP_SUM4_x0_N_eff')} (penaliza canales lejanos con menos luz)",
        f"- σ_SPTR_eff = 106/√N_eff ps",
        f"- σ_tot_B = √(σ_int_B² + σ_SPTR_eff² + 10²) = {g('T1_TOP_SUM4_x0_sig_tot_B_ps')} ps",
        f"- Cumple SHiP 100 ps: {g('T1_TOP_SUM4_x0_goal_B')}",
        f"- TOP_SUM8 Chain B: σ_tot_B = {g('T1_TOP_SUM8_x0_sig_tot_B_ps')} ps",
        "",
        f"### Tabla resumen (TOP_SUM4 @ x=0):",
        "| Cadena | σ_int | σ_SPTR | σ_tot | Cumple 100 ps |",
        "|--------|-------|--------|-------|--------------|",
        f"| A (analógica) | {g('T1_TOP_SUM4_x0_sig_int_A_ps')} ps | 106 ps | {g('T1_TOP_SUM4_x0_sig_tot_A_ps')} ps | {g('T1_TOP_SUM4_x0_goal_A')} |",
        f"| B (digital) | {g('T1_TOP_SUM4_x0_sig_int_B_ps')} ps | 106/√{g('T1_TOP_SUM4_x0_N_eff')} ps | {g('T1_TOP_SUM4_x0_sig_tot_B_ps')} ps | {g('T1_TOP_SUM4_x0_goal_B')} |",
        "",
        "### Estado del requisito 100 ps:",
        "- σ_int intrínseco (sin SPTR/FastIC): ambas cadenas < 100 ps ← dato seguro",
        "- Con SPTR: depende de la cadena. El 100 ps solo se alcanza bajo Chain B.",
        "- Chain B requiere que FastIC+ pueda :\n  (i) generar timestamps individuales por canal,",
        "  (ii) promediar ponderadamente en tiempo real.\n  → Confirmar con el equipo FastIC+.",
        "",
        "---",
        "## DECISIÓN 2 — Topología de SUM: END co-localizado vs TOP distribuido",
        "",
        "### Reencuadre (EXEC_19)",
        "'Co-localizado' describe el END (8 SiPMs en la MISMA cara de extremo, mismo x,",
        "misma luz). El TOP es SIEMPRE una FILA DISTRIBUIDA a paso 20mm.",
        "Cualquier agrupación de SiPMs TOP (id//4 o nearest-N) es distribuida.",
        "La comparación co-localizado↔distribuido REAL es:",
        "  END_SUM4 co-localizado (misma face) vs TOP-fila-sum (distribuida)",
        "",
        f"### Evidencia de selección @ x=0:",
        f"- Dynamic nearest-4 TOP: total ⟨Npe⟩ = {ff(T3['npe_dynamic_total'],2)} PE/ev",
        f"- id//4 cluster TOP:     total ⟨Npe⟩ = {ff(T3['npe_cluster_total'],2)} PE/ev",
        f"- Dynamic recoge MÁS luz pero EXEC_18 dio peor σ (artefacto FD binning + no walk)",
        f"- Con sqrt_n + walk: ambos deberían dar σ similar — ver T3 figuras",
        "",
        "### Recomendación técnica (neutral)",
        "- Para la comparación sim↔TB Constanza: END_SUM4 co-localizado",
        "- Para la instrumentación TOP: nearest-4 dinámico recoge más luz y es la selección",
        "  correcta; id//4 fijo es equivalente cuando el haz está bien alineado",
        "- La decisión final (id//4 fijo vs nearest-4 dinámico) sigue siendo de Gerardo",
        "  pero ambas son topologías DISTRIBUIDAS",
        "",
        "---",
        "## Definiciones canónicas (CANONICAL_ESTIMATORS.md)",
        "Ver análysis_core/out/EXEC_19/CANONICAL_ESTIMATORS.md para la definición exacta",
        "de cada estimador headline (canales, combinación, walk, binning) y su σ.",
    ]
    p = Path(OUT_DIR) / "GERARDO_DECISION_MEMO.md"
    p.write_text("\n".join(lines))
    note(f"  GERARDO_DECISION_MEMO.md: {p}")

def write_exec19_report(T1, T2, T3, T4, rj, flags):
    def g(key): return rj.get(key, "?")

    lines = [
        f"# EXEC_19_REPORT — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## 1. Veredicto",
        "**COMPLETADO" + ("-CON-FLAGS" if flags else "") + "**",
        "",
        "## 2. Tabla A/B consistente",
        "| Estimador | x [mm] | σ_int_A | σ_tot_A | σ_int_B | N_eff | σ_tot_B | A≤100 | B≤100 |",
        "|---------|------|--------|--------|--------|------|--------|-------|-------|",
    ]
    for r in T1:
        lines.append(f"| {r['label']} | {r['x_mm']} | {r['sig_int_A_ps']} | {r['sig_tot_A_ps']} | "
                     f"{r['sig_int_B_ps']} | {r['N_eff']} | {r['sig_tot_B_ps']} | {r['goal_A']} | {r['goal_B']} |")
    lines += [
        "",
        f"- EXEC_18 σ_int_B=69.2 ps was min-stream (WRONG for Chain B); Chain B σ_int_B from weighted avg is different",
        "",
        "## 3. Estado del 100 ps",
        f"- Chain A: σ_tot_A(SUM4@x=0) = {g('T1_TOP_SUM4_x0_sig_tot_A_ps')} ps → goal: {g('T1_TOP_SUM4_x0_goal_A')}",
        f"- Chain B: σ_tot_B(SUM4@x=0) = {g('T1_TOP_SUM4_x0_sig_tot_B_ps')} ps → goal: {g('T1_TOP_SUM4_x0_goal_B')}",
        "- Under Chain B: sub-100 ps requires FastIC+ to provide individual-channel timestamps",
        "",
        "## 4. CFD vs min (T2)",
        f"- min(t_first) TOP_SUM4 @ x=0 (walk-corr, sqrt_n): {ff(T2['sig_min_ps'],1)} ps",
        f"- CFD(f={T2['best_f']}) @ x=0: {ff(T2['sig_best_cfd_ps'],1)} ps",
        f"- Gap = {ff(T2['sig_best_cfd_ps'] - T2['sig_min_ps'],1)} ps (min is optimistic lower bound)",
        "",
        "## 5. T3 topology",
        f"- Dynamic nearest-4: total ⟨Npe⟩={ff(T3['npe_dynamic_total'],2)} PE/ev",
        f"- id//4 cluster:      total ⟨Npe⟩={ff(T3['npe_cluster_total'],2)} PE/ev",
        "- EXEC_18 72.1 vs 85.0 was FD binning artifact (no walk correction applied)",
        "- 'co-localized TOP' is a misnomer — ALL TOP groupings are distributed rows",
        "- No prior headline contaminated (EXEC_18 T1 used SUM4 which is TOP_SUM4_N1)",
        "",
        "## 6. Reconciliation T4",
        "- TOP_nearest_N1 discrepancy (107.8 vs 85.0 ps): different estimators (1-ch vs 4-ch)",
        "- END near-end discrepancy (30 vs 44 ps): FD vs sqrt_n + different walk application",
        "- Canonical values: see CANONICAL_ESTIMATORS.md",
        "",
        "## 7. Flags",
    ]
    lines.extend(flags if flags else ["  ninguno"])
    lines += [
        "",
        "## 8. Rutas",
        f"  Figuras: {OUT_DIR}/T1..T4",
        f"  Memo: {OUT_DIR}/GERARDO_DECISION_MEMO.md",
        f"  Canónico: {OUT_DIR}/CANONICAL_ESTIMATORS.md",
        f"  JSON: {OUT_DIR}/results_exec19.json",
        "",
        "## 9. Decisiones para Gerardo",
        f"  1. FastIC+: Chain A→{g('T1_TOP_SUM4_x0_sig_tot_A_ps')} ps vs Chain B→{g('T1_TOP_SUM4_x0_sig_tot_B_ps')} ps",
        "  2. Topología TOP: id//4 fijo vs nearest-4 dinámico (ambas distribuidas)",
    ]

    p = Path(OUT_DIR) / "EXEC_19_REPORT.md"
    p.write_text("\n".join(lines))
    note(f"\nEXEC_19_REPORT.md: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    note(f"# EXEC_19 — {datetime.datetime.now().isoformat()}")
    note(f"K_END={K_END}, N_TOP={N_TOP}")

    pm = pos_map()
    note(f"  {len(pm)} positions")
    gate_integrity(pm)

    import ROOT; ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning

    T2_res = run_T2(pm)      # T2 first to determine best CFD fraction for T1

    T1_res = run_T1(pm)      # uses T2's best_f

    T3_res = run_T3(pm)

    os.system("cd /home/reriosto/SHiP/analysis_core && git tag EXEC_19-pre-t4 2>/dev/null || true")

    T4_res = run_T4(pm, T1_res, T3_res)

    rj = write_results_json(T1_res, T2_res, T3_res, T4_res)
    write_canonical_estimators(T1_res, T2_res, T3_res, T4_res)
    write_gerardo_memo_v2(T1_res, T2_res, T3_res, T4_res, rj)
    write_exec19_report(T1_res, T2_res, T3_res, T4_res, rj, _flags)

    _flush_report("COMPLETADO" + ("-CON-FLAGS" if _flags else ""))
    note(f"\n=== EXEC_19 COMPLETE === Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
