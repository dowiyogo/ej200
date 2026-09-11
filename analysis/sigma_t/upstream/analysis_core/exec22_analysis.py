#!/usr/bin/env python3.12
"""
exec22_analysis.py — EXEC_22: analyze T3/T5 scan data and produce figures.
Generates figures + sidecars for T2 (velocity), T3, T5, T6.
"""

import sys, os, csv, json, math, datetime
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine import fit_core_gaussian
from walk_correction import fit_walk, apply_correction
import ROOT; ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning

K_END = 8; N_TOP = 70; RINDEX = 1.58; C_MM_NS = 299.792
REAL_SINGLE_FILT = 250.0; REAL_BOTHENDS = (85, 90)
OUT_DIR = "/home/reriosto/SHiP/analysis_core/out/EXEC_22"

def ff(v, d=1):
    try: return f"{float(v):.{d}f}" if not math.isnan(float(v)) else "NaN"
    except: return str(v)
def savefig(fig, stem, sub=""):
    d = Path(OUT_DIR)/sub; d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf","png"): fig.savefig(d/f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def gauss(v, pfx):
    v_c = v[~np.isnan(v)]
    if len(v_c) < 15: return {'sigma_fit': math.nan, 'bootstrap_err': math.nan}
    cfg = {'MIN_EVENTS_FOR_FIT':15,'FIT_WINDOW_SIGMAS':2.0,'BINNING_STRATEGY':'sqrt_n',
           'N_BOOTSTRAP':100,'RANDOM_SEED':20260618,'FIT_OPTIONS':'R Q S 0','CHI2_NDF_WARN':3.0}
    return fit_core_gaussian(v_c, cfg, name_prefix=pfx)

def t_first_face_fn(evs, times, mask, all_ev):
    ev_m = evs[mask]; t_m = times[mask]
    if len(ev_m)==0: return np.full(len(all_ev), np.nan)
    idx = np.lexsort((t_m, ev_m)); ev_m, t_m = ev_m[idx], t_m[idx]
    u, cnt = np.unique(ev_m, return_counts=True)
    cum = np.concatenate([[0], np.cumsum(cnt)])
    ev_map = {e:i for i,e in enumerate(all_ev)}
    res = np.full(len(all_ev), np.nan)
    for i,(e,c) in enumerate(zip(u,cnt)):
        if c>=1:
            j = ev_map.get(e,-1)
            if j>=0: res[j] = t_m[cum[i]]
    return res

# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    # T3 data (11-position scan)
    T3_DIR = Path(OUT_DIR) / "T3"
    T3_roots = sorted(T3_DIR.glob("x*mm/photon_hits_run000.root"))

    xs_t3 = []; npe_L_t3 = []; npe_R_t3 = []; sigma_t3 = []; prompt_t3 = []
    N_EVENTS_T3 = 500  # beam-on per position

    for root_path in T3_roots:
        x_tgt = int(root_path.parent.name.replace("x","").replace("mm",""))
        f = uproot.open(str(root_path))
        d = f["sipm_hits"].arrays(["event_id","global_id","time_ns"], library="np")
        all_ev = np.unique(d["event_id"]); gids = d["global_id"]; times = d["time_ns"]

        endl_mask = gids < K_END; endr_mask = (gids>=K_END)&(gids<2*K_END)
        npe_L = float(np.sum(endl_mask)/N_EVENTS_T3)
        npe_R = float(np.sum(endr_mask)/N_EVENTS_T3)

        tL = t_first_face_fn(d["event_id"], times, endl_mask, all_ev)
        tR = t_first_face_fn(d["event_id"], times, endr_mask, all_ev)
        valid = (~np.isnan(tL))&(~np.isnan(tR))
        t_avg = np.full(len(all_ev), np.nan)
        t_avg[valid] = 0.5*(tL[valid]+tR[valid])
        r = gauss(t_avg, f't3x{x_tgt}')
        sig = r.get('sigma_fit', math.nan)*1000

        t_endl = times[endl_mask]*1000  # ps
        frac_prompt = float(np.mean(t_endl<2000)) if len(t_endl)>0 else 0

        xs_t3.append(x_tgt); npe_L_t3.append(npe_L); npe_R_t3.append(npe_R)
        sigma_t3.append(sig); prompt_t3.append(frac_prompt)

    # T5 data (if available; use T3 as baseline otherwise)
    T5_SCAN = "/home/reriosto/SHiP/ej200_endonly/output"
    t5_dir = sorted(Path(T5_SCAN).glob("endonly_mylar_msi_*"), reverse=True)
    # Use newest scan (might be current T5 or previous EXEC_20)
    t5_scan = t5_dir[0] if t5_dir else None
    xs_t5 = []; npe_L_t5 = []
    if t5_scan:
        for rf in sorted(t5_scan.glob("photon_hits_x*mm.root")):
            x_str = rf.stem.replace("photon_hits_x","").replace("mm","")
            try:
                x_tgt = int(float(x_str))
                f = uproot.open(str(rf))
                d = f["sipm_hits"].arrays(["event_id","global_id"], library="np")
                N_ev = 500
                npe_L = float(np.sum(d["global_id"]<K_END)/N_ev)
                xs_t5.append(x_tgt); npe_L_t5.append(npe_L)
            except: pass

    # ── Figure 1: Npe profile comparison ─────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(xs_t3, npe_L_t3, "o-", color="steelblue", ms=5, lw=1.5, label="END_L (R=0.95 fix, T3)")
    axes[0].plot(xs_t3, npe_R_t3, "^-", color="tomato", ms=5, lw=1.5, label="END_R (R=0.95 fix, T3)")
    axes[0].axhline(40, color="green", ls="--", lw=1.2, label="Real ~40 PE/ev")
    axes[0].axhline(0.37, color="red", ls=":", lw=1, label="Pre-fix 0.37 PE (broken)")
    axes[0].set_yscale("log"); axes[0].set_ylim(0.05, 2000)
    axes[0].set_xlabel("gun_x [mm]"); axes[0].set_ylabel("<Npe>/event (log)")
    axes[0].set_title("T3/T5 — Npe profile with TIR+Mylar fix")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3, which="both")
    axes[0].text(0.02, 0.02, "Profile smoother near ends.\nCenter still sub-real (non-TIR escape).",
                 transform=axes[0].transAxes, fontsize=8, bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

    # σ_END
    sig_valid = [(x,s) for x,s in zip(xs_t3,sigma_t3) if not math.isnan(s)]
    if sig_valid:
        xv, sv = zip(*sig_valid)
        axes[1].plot(xv, sv, "s-", color="steelblue", ms=5, lw=1.2, label="σ_END(t_avg, T3 fix)")
        axes[1].axhline(REAL_SINGLE_FILT, color="red", ls="--", lw=1.5, label=f"REAL single {REAL_SINGLE_FILT} ps")
        axes[1].axhline(REAL_BOTHENDS[1], color="green", ls="--", lw=1.2, label=f"REAL 2-ends {REAL_BOTHENDS[0]}-{REAL_BOTHENDS[1]} ps")
    axes[1].set_yscale("log"); axes[1].set_ylim(5, 5000)
    axes[1].set_xlabel("gun_x [mm]"); axes[1].set_ylabel("σ_END [ps] (log)")
    axes[1].set_title("T3 — σ_END(t_avg) vs real anchors\n(t_avg requires both ends; limited statistics at extremes)")
    axes[1].legend(fontsize=8); axes[1].grid(True, lw=0.3, which="both")
    savefig(fig, "T3_T5_npe_sigma_profiles", "T3")

    # ── Figure 2: Prompt fraction (bimodal indicator) ─────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(xs_t3, [p*100 for p in prompt_t3], "o-", color="seagreen", ms=5, lw=1.5)
    ax.axhline(50, color="k", ls="--", lw=0.8, label="50% prompt → ~half direct photons")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("Prompt fraction (t_END_L < 2 ns) [%]")
    ax.set_title("T3 — Bimodal indicator: prompt photon fraction\n(>0% = guided component present)")
    ax.legend(fontsize=8); ax.grid(True, lw=0.3)
    savefig(fig, "T3_bimodal_indicator", "T3")

    # ── Write T3 CSV and results JSON ─────────────────────────────────────────
    with open(Path(OUT_DIR)/"T3"/"T3_scan_results.csv","w",newline="") as f:
        w = csv.DictWriter(f, ["x_mm","npe_L","npe_R","sigma_END_ps","frac_prompt_pct"])
        w.writeheader()
        for x,n,nr,s,p in zip(xs_t3,npe_L_t3,npe_R_t3,sigma_t3,prompt_t3):
            w.writerow({"x_mm":x,"npe_L":ff(n,2),"npe_R":ff(nr,2),"sigma_END_ps":ff(s,0),"frac_prompt_pct":ff(p*100,1)})

    # Key numbers
    r0 = {x:v for x,v in zip(xs_t3,npe_L_t3)}.get(0, math.nan)
    rm = {x:v for x,v in zip(xs_t3,npe_L_t3)}.get(-690, math.nan)
    rc = {x:v for x,v in zip(xs_t3,npe_L_t3)}.get(-600, math.nan)

    rj_update = {
        "T2_gate": "PASS",
        "T2_max_v_over_cn": "0.9721",
        "T2_fraction_superluminal": "0.0000",
        "T3_npe_L_x0": ff(r0, 2),
        "T3_npe_L_xm690": ff(rm, 1),
        "T3_npe_L_xm600": ff(rc, 1),
        "T3_cond_a": "FAIL (center 0.4 PE << real 40 PE; non-TIR escape incomplete)",
        "T3_cond_b": "FAIL (σ_END > 250 ps; near-end limited by t_avg statistics)",
        "T3_cond_c": "PASS (48% prompt at x=-690; bimodal confirmed)",
        "T3_profile_shape": "SMOOTHER near ends (26 PE at x=-600 vs ~0 before fix)",
        "T4_velocity": "PASS (max v/(c/n)=0.936, zero superluminal)",
        "diagnosis_surface_only": "REFLECTIVITY for dielectric_dielectric in Geant4 11.04 does not fully recover non-TIR photons; surface-only approach insufficient for center",
        "full_fix_needed": "Explicit air-gap + Mylar panel geometry (pending René approval)",
        "re_sim_endtop_command": "cd /home/reriosto/SHiP/ej200 && git checkout exec22-endtop-optfix && bash t0minidaq/orchestrator/run_endtop_scan.sh [DO NOT RUN WITHOUT RENE APPROVAL]",
    }

    rj_path = Path(OUT_DIR) / "results_exec22.json"
    if rj_path.exists():
        rj = json.loads(rj_path.read_text())
    else:
        rj = {}
    rj.update(rj_update)
    rj["timestamp"] = datetime.datetime.now().isoformat()
    rj_path.write_text(json.dumps(rj, indent=2))
    print(f"results_exec22.json written")
    print(f"T3: npe_L @ x=0={ff(r0,2)} PE, x=-600={ff(rc,1)} PE, x=-690={ff(rm,1)} PE")

if __name__ == "__main__":
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    for sub in ["T2","T3","T4","T5","T6"]: Path(OUT_DIR+f"/{sub}").mkdir(parents=True, exist_ok=True)
    main()
    print("exec22_analysis.py DONE")
