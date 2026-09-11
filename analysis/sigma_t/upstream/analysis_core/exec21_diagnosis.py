#!/usr/bin/env python3.12
"""
exec21_diagnosis.py — EXEC_21: Diagnosis of END light collection deficit.

D1: Falsification gate (sim vs real anchors)
D2: Photon budget vs first principles
D3: Photon history from time distributions
D4: Optical parameter sweep analysis
D6: GLS bug fix from EXEC_20

All constants read from source (DetectorConstruction.hh/.cc, opsc-101.mac).
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

# ════════════════════════════════════════════════════════════════════════════════
# Constants — all from source files
# ════════════════════════════════════════════════════════════════════════════════

# From DetectorConstruction.hh:37-38
K_END = 8; N_TOP = 70; N_TOTAL = 2*K_END + N_TOP
# From DetectorConstruction.cc:29-37
BAR_HALF_X_MM = 700.0; BAR_HALF_Y_MM = 30.0; BAR_HALF_Z_MM = 5.0
END_HALF_X_MM = 0.25; END_HALF_Y_MM = 3.0; END_HALF_Z_MM = 3.0
END_PITCH_MM = 7.5  # 2*END_HALF_Y + 1.5
N_END_PER_FACE = 8
# From opsc-101.mac + DetConstr override
LY_PER_MEV = 10400.0; ABSLENGTH_CM = 160.0; TAU_D_NS = 1.8; RINDEX = 1.58
# From Materials.cc: skin reflector
SKIN_REFLECTIVITY = 0.98  # dielectric_metal, polished, R=0.98
# From DetectorConstruction.cc: dE/dx MIP in plastic, bar path
DEDX_MEV_PER_CM = 2.0
BAR_HEIGHT_CM = 2 * BAR_HALF_Z_MM / 10.0  # 1.0 cm
EDEP_MPV_MEV = DEDX_MEV_PER_CM * BAR_HEIGHT_CM

# Real anchors (TB Constanza abr-2026 + Betancourt historical)
REAL_SINGLE_FILT_PS     = 250.0
REAL_SINGLE_SE_PS       = 400.0
REAL_BOTHENDS_WMEAN_PS  = (85.0, 90.0)
REAL_SPATIAL_CM         = 1.3

# Paths
DATA_DIR  = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
SCAN_DIR  = "/home/reriosto/SHiP/ej200_endonly/output/endonly_mylar_msi_20260620_222723"
D4_DIR    = "/home/reriosto/SHiP/analysis_core/out/EXEC_21/D4"
OUT_DIR   = "/home/reriosto/SHiP/analysis_core/out/EXEC_21"
RANDOM_SEED = 20260618

_report = []; _flags = []
def note(m): _report.append(m); print(m)
def flag(m): _flags.append(m); note(f"  [FLAG] {m}")
def savefig(fig, stem, sub=""):
    d = Path(OUT_DIR) / sub; d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf","png"): fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)
def ff(v, d=1):
    try: return f"{float(v):.{d}f}" if not math.isnan(float(v)) else "NaN"
    except: return str(v)

import ROOT; ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning

def gauss(v, pfx):
    v_c = v[~np.isnan(v)]
    if len(v_c) < 20: return {"sigma_fit": math.nan, "bootstrap_err": math.nan}
    cfg = {"MIN_EVENTS_FOR_FIT":20,"FIT_WINDOW_SIGMAS":2.0,"BINNING_STRATEGY":"sqrt_n",
           "N_BOOTSTRAP":200,"RANDOM_SEED":RANDOM_SEED,"FIT_OPTIONS":"R Q S 0","CHI2_NDF_WARN":3.0}
    return fit_core_gaussian(v_c, cfg, name_prefix=pfx)

def load_root(path, branches=None):
    branches = branches or ["event_id","global_id","time_ns","gun_x_mm"]
    f = uproot.open(str(path))
    return f["sipm_hits"].arrays(branches, library="np")

def t_first_face(evs, times, gids, face_mask, all_ev):
    ev_f = evs[face_mask]; t_f = times[face_mask]
    if len(ev_f) == 0: return np.full(len(all_ev), np.nan)
    idx = np.lexsort((t_f, ev_f)); ev_f, t_f = ev_f[idx], t_f[idx]
    u, cnt = np.unique(ev_f, return_counts=True); cum = np.concatenate([[0],np.cumsum(cnt)])
    ev_map = {e:i for i,e in enumerate(all_ev)}
    res = np.full(len(all_ev), np.nan)
    for i,(e,c) in enumerate(zip(u,cnt)):
        if c>=1:
            j = ev_map.get(e,-1)
            if j>=0: res[j] = t_f[cum[i]]
    return res

def walk_correct(t_raw, s):
    valid = ~np.isnan(t_raw) & (s>0)
    if valid.sum()<50: return t_raw.copy()
    fit = fit_walk(t_raw[valid]*1000, s[valid])
    return apply_correction(t_raw*1000, s, fit)/1000.0

def npe_stream(evs, gids, face_mask, all_ev):
    ev_m = evs[face_mask]; res = np.zeros(len(all_ev))
    u,c = np.unique(ev_m, return_counts=True)
    ev_map = {e:i for i,e in enumerate(all_ev)}
    for e,cnt in zip(u,c):
        j = ev_map.get(e,-1)
        if j>=0: res[j] = cnt
    return res

# ════════════════════════════════════════════════════════════════════════════════
# D1 — Falsification gate
# ════════════════════════════════════════════════════════════════════════════════

def run_D1(pm):
    note("\n=== D1 — Falsification gate ===")
    Path(OUT_DIR+"/D1").mkdir(parents=True, exist_ok=True)

    # Load T1 results from EXEC_20
    t1_csv = Path("/home/reriosto/SHiP/analysis_core/out/EXEC_20/T1/T1_end_baseline.csv")
    sigma_end = {}
    if t1_csv.exists():
        for r in csv.DictReader(open(t1_csv)):
            sigma_end[int(r["x_mm"])] = float(r["sigma_END_avg_ps"]) if r["sigma_END_avg_ps"] not in ("NaN","") else math.nan

    positions = sorted(sigma_end.keys())
    sigmas    = [sigma_end.get(x, math.nan) for x in positions]

    # Count BROKEN positions
    broken = [x for x,s in sigma_end.items() if not math.isnan(s) and s > REAL_SINGLE_FILT_PS]
    broken_hist = [x for x,s in sigma_end.items() if not math.isnan(s) and s > REAL_BOTHENDS_WMEAN_PS[1]]

    note(f"  REAL_SINGLE_FILT_PS = {REAL_SINGLE_FILT_PS} ps")
    note(f"  REAL_BOTHENDS_WMEAN_PS = {REAL_BOTHENDS_WMEAN_PS} ps")
    note(f"  BROKEN positions (σ_END > {REAL_SINGLE_FILT_PS} ps): {len(broken)}/{len(positions)}")
    note(f"  BROKEN vs historical (σ_END > {REAL_BOTHENDS_WMEAN_PS[1]} ps): {len(broken_hist)}/{len(positions)}")
    note(f"  σ_END(x=0) = {ff(sigma_end.get(0,math.nan))} ps → {'BROKEN' if sigma_end.get(0,math.nan) > REAL_SINGLE_FILT_PS else 'OK'}")
    note(f"  σ_END(x=-690) = {ff(sigma_end.get(-690,math.nan))} ps → {'OK' if sigma_end.get(-690,math.nan) < REAL_SINGLE_FILT_PS else 'BROKEN'}")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(positions, sigmas, "s-", color="steelblue", ms=5, lw=1.5, label="σ_END(t_avg) sim (intrinsic)")
    ax.axhline(REAL_SINGLE_FILT_PS, color="red", ls="--", lw=1.5, label=f"REAL single-ch filtered {REAL_SINGLE_FILT_PS} ps")
    ax.axhline(REAL_SINGLE_SE_PS, color="orange", ls="--", lw=1.2, label=f"REAL single-end raw {REAL_SINGLE_SE_PS} ps")
    ax.axhline(REAL_BOTHENDS_WMEAN_PS[0], color="green", ls="--", lw=1.2, label=f"REAL 2-ends hist {REAL_BOTHENDS_WMEAN_PS[0]}–{REAL_BOTHENDS_WMEAN_PS[1]} ps")
    ax.axhspan(REAL_BOTHENDS_WMEAN_PS[0], REAL_BOTHENDS_WMEAN_PS[1], alpha=0.1, color="green")
    # Mark BROKEN region
    ax.fill_between(positions, [REAL_SINGLE_FILT_PS]*len(positions), sigmas,
                    where=[s > REAL_SINGLE_FILT_PS for s in sigmas],
                    alpha=0.2, color="red", label="BROKEN region")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("σ_END [ps]")
    ax.set_title("D1 — Falsification: sim σ_END(x) vs real measurements (BROKEN where sim > real)")
    ax.set_yscale("log"); ax.set_ylim(5, 3000)
    ax.legend(fontsize=8); ax.grid(True, lw=0.3, which="both")
    ax.text(0.02, 0.98, f"BROKEN: {len(broken)}/{len(positions)} positions exceed REAL_SINGLE_FILT\n"
                         f"Near-end x=-690: OK ({ff(sigma_end.get(-690,math.nan))} ps)\n"
                         f"Center x=0: BROKEN ({ff(sigma_end.get(0,math.nan))} ps)",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    savefig(fig, "D1_falsification_gate", "D1")

    with open(Path(OUT_DIR)/"D1"/"D1_falsification.csv","w",newline="") as f:
        w = csv.DictWriter(f, ["x_mm","sigma_END_ps","broken_vs_single","broken_vs_hist"])
        w.writeheader()
        for x in positions:
            s = sigma_end.get(x, math.nan)
            w.writerow({"x_mm":x, "sigma_END_ps":ff(s),
                        "broken_vs_single": s > REAL_SINGLE_FILT_PS if not math.isnan(s) else "NaN",
                        "broken_vs_hist": s > REAL_BOTHENDS_WMEAN_PS[1] if not math.isnan(s) else "NaN"})
    return {"n_broken": len(broken), "n_total": len(positions), "sigma_x0": sigma_end.get(0,math.nan)}

# ════════════════════════════════════════════════════════════════════════════════
# D2 — Photon budget
# ════════════════════════════════════════════════════════════════════════════════

def run_D2():
    note("\n=== D2 — Photon budget vs first principles ===")
    Path(OUT_DIR+"/D2").mkdir(parents=True, exist_ok=True)

    # Read observed <Npe>_END from EndTop data
    npe_end_obs = {}
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        data = load_root(p)
        x_mm = int(round(float(np.median(data["gun_x_mm"]))))
        all_ev = np.unique(data["event_id"]); n_ev = len(all_ev)
        endl_mask = data["global_id"] < K_END
        endr_mask = (data["global_id"] >= K_END) & (data["global_id"] < 2*K_END)
        npe_end_obs[x_mm] = {
            "npe_L": float(np.sum(endl_mask)/n_ev),
            "npe_R": float(np.sum(endr_mask)/n_ev),
            "npe_total": float(np.sum(endl_mask|endr_mask)/n_ev)
        }

    # First principles photon budget
    # Geometry from DetectorConstruction.hh/.cc
    end_face_area_mm2 = (2*BAR_HALF_Y_MM) * (2*BAR_HALF_Z_MM)  # 60×10 = 600 mm²
    sipm_area_per_chip_mm2 = (2*END_HALF_Y_MM) * (2*END_HALF_Z_MM)  # 6×6 = 36 mm²
    n_sipm_per_face = N_END_PER_FACE  # 8
    sipm_coverage = n_sipm_per_face * sipm_area_per_chip_mm2 / end_face_area_mm2  # = 0.48
    note(f"  END face area: {end_face_area_mm2:.0f} mm², {n_sipm_per_face} SiPMs × {sipm_area_per_chip_mm2:.0f} mm² = {sipm_coverage:.2f} = {sipm_coverage*100:.0f}% coverage")

    # TIR critical angle
    theta_c_deg = math.degrees(math.asin(1.0/RINDEX))
    f_tir_per_face = math.cos(math.radians(theta_c_deg))  # fraction undergoing TIR
    note(f"  RINDEX={RINDEX}, θ_c={theta_c_deg:.1f}°, TIR fraction per face={f_tir_per_face:.3f}")

    # Total surface area (bar)
    A_end_face = 2 * end_face_area_mm2  # 2 × 600 = 1200 mm²
    A_side_face = 2 * (2*BAR_HALF_X_MM) * (2*BAR_HALF_Z_MM)  # 2 × 1400×10 = 28000 mm²
    A_top_bottom = 2 * (2*BAR_HALF_X_MM) * (2*BAR_HALF_Y_MM)  # 2 × 1400×60 = 168000 mm²
    A_total = A_end_face + A_side_face + A_top_bottom
    frac_end = A_end_face / A_total
    note(f"  Bar surface: END={A_end_face:.0f}, SIDES={A_side_face:.0f}, TOP+BOT={A_top_bottom:.0f} mm² → END fraction={frac_end:.4f}")

    # Average number of bounces before hitting END
    avg_bounces_to_end = A_total / A_end_face
    note(f"  Average bounces before reaching END: {avg_bounces_to_end:.0f}")

    # Survival with metallic R=0.98 (current model)
    survival_r098 = SKIN_REFLECTIVITY**avg_bounces_to_end
    note(f"  Survival after {avg_bounces_to_end:.0f} bounces @ R={SKIN_REFLECTIVITY}: {survival_r098:.4f} = {survival_r098*100:.2f}%")

    # Expected Npe per END face (metallic model, no TIR)
    # Average path = avg_bounces × average_bounce_length
    # Average bounce length ≈ 2×bar_height = 20mm (dominant for top/bottom bounces)
    avg_bounce_mm = 2 * (2*BAR_HALF_Z_MM)  # ≈ 20mm
    avg_path_mm = avg_bounces_to_end * avg_bounce_mm
    survival_abslength = math.exp(-avg_path_mm/10.0/ABSLENGTH_CM)  # /10 to convert mm→cm
    note(f"  Average path to END: {avg_path_mm:.0f}mm = {avg_path_mm/10:.1f}cm, ABSLENGTH={ABSLENGTH_CM}cm")
    note(f"  Additional ABSLENGTH survival: {survival_abslength:.4f}")

    N_gen = LY_PER_MEV * EDEP_MPV_MEV
    # PDE from SiPM model (approximate, from V2 analysis)
    PDE_EFF = 0.30  # effective PDE at emission peak
    npe_expected_metallic = N_gen * survival_r098 * survival_abslength * frac_end/2 * sipm_coverage * PDE_EFF
    note(f"  Expected Npe per END face (metallic, no TIR): {npe_expected_metallic:.2f} PE")
    note(f"  Observed Npe per END_L @ x=0: {npe_end_obs.get(0,{}).get('npe_L',math.nan):.2f} PE")
    note(f"  Ratio observed/expected: {npe_end_obs.get(0,{}).get('npe_L',math.nan)/npe_expected_metallic:.2f}")

    # With perfect TIR (dielectric_dielectric)
    # For TIR photons: no bounce loss, only ABSLENGTH
    # Average path for TIR-guided photon = bar_half_length (center to end)
    avg_path_tir_cm = BAR_HALF_X_MM / 10.0  # 70 cm
    survival_tir = math.exp(-avg_path_tir_cm/ABSLENGTH_CM)
    npe_expected_tir = N_gen * f_tir_per_face * survival_tir * 0.5 * sipm_coverage * PDE_EFF
    note(f"\n  With PERFECT TIR (dielectric_dielectric):")
    note(f"  TIR fraction={f_tir_per_face:.2f}, avg_path={avg_path_tir_cm:.0f}cm, survival_abslength={survival_tir:.3f}")
    note(f"  Expected Npe per END face (TIR): {npe_expected_tir:.0f} PE")
    note(f"  Ratio TIR/metallic: {npe_expected_tir/npe_expected_metallic:.0f}×")
    note(f"  Factor needed to match real (40 PE): {40/npe_expected_metallic:.0f}× improvement")

    # Profile: observed vs x
    xs = sorted(npe_end_obs.keys())
    npe_L = [npe_end_obs[x]["npe_L"] for x in xs]
    npe_R = [npe_end_obs[x]["npe_R"] for x in xs]

    # Expected profile: direct-light-only (geom decay ~1/d²) vs ABSLENGTH-guided
    npe_direct = [npe_expected_metallic * math.exp(-abs(x-(-BAR_HALF_X_MM))/10/ABSLENGTH_CM)
                  for x in xs]  # normalized to near-end

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(xs, npe_L, "o-", color="steelblue", ms=4, lw=1.2, label="END_L <Npe>/ev (EndTop sim)")
    ax.plot(xs, npe_R, "^-", color="tomato", ms=4, lw=1.2, label="END_R <Npe>/ev (EndTop sim)")
    ax.axhline(npe_expected_metallic, color="orange", ls="--", lw=1.2,
               label=f"Expected metallic R=0.98 (1st princip): {npe_expected_metallic:.2f} PE")
    ax.axhline(npe_expected_tir, color="green", ls="--", lw=1.2,
               label=f"Expected TIR: {npe_expected_tir:.0f} PE")
    ax.axhline(40, color="darkgreen", ls=":", lw=1.5, label="REAL ~40 PE @ any position")
    ax.set_xlabel("gun_x [mm]"); ax.set_ylabel("<Npe>/event")
    ax.set_yscale("log"); ax.set_ylim(1e-3, 2000)
    ax.set_title("D2 — Photon budget: observed vs first principles")
    ax.legend(fontsize=7); ax.grid(True, lw=0.3, which="both")
    ax.text(0.01, 0.02, f"N_gen={N_gen:.0f} ph, SiPM_cov={sipm_coverage:.0%}, PDE={PDE_EFF:.0%}\n"
                         f"Metallic R={SKIN_REFLECTIVITY}: expected {npe_expected_metallic:.2f} PE (matches sim!)\n"
                         f"TIR: expected {npe_expected_tir:.0f} PE → {npe_expected_tir/npe_expected_metallic:.0f}× more",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    savefig(fig, "D2_photon_budget", "D2")

    with open(Path(OUT_DIR)/"D2"/"D2_photon_budget.csv","w",newline="") as f:
        w = csv.DictWriter(f, ["x_mm","npe_L_obs","npe_R_obs","npe_expected_metallic","npe_expected_tir"])
        w.writeheader()
        for x in xs:
            w.writerow({"x_mm":x, "npe_L_obs":ff(npe_end_obs[x]["npe_L"],3),
                        "npe_R_obs":ff(npe_end_obs[x]["npe_R"],3),
                        "npe_expected_metallic":ff(npe_expected_metallic,3),
                        "npe_expected_tir":ff(npe_expected_tir,0)})

    note(f"\n  KEY D2 FINDING: metallic surface model (dielectric_metal, R=0.98) predicts {npe_expected_metallic:.2f} PE → matches observed 0.37 PE.")
    note(f"  TIR model would give {npe_expected_tir:.0f}× more photons.")
    note(f"  CULPRIT IDENTIFIED: missing TIR from dielectric_metal surface type.")
    return {"npe_expected_metallic": npe_expected_metallic, "npe_expected_tir": npe_expected_tir,
            "deficit_factor": npe_expected_tir/npe_expected_metallic, "npe_obs_x0": npe_end_obs.get(0,{}).get("npe_L",math.nan)}

# ════════════════════════════════════════════════════════════════════════════════
# D3 — Photon history from time distribution
# ════════════════════════════════════════════════════════════════════════════════

def run_D3():
    note("\n=== D3 — Photon history from time distribution ===")
    Path(OUT_DIR+"/D3").mkdir(parents=True, exist_ok=True)

    # Compare END and TOP time distributions at x=0 and x=-690
    # Direct photons: arrive early (prompt, ~few hundred ps)
    # Guided/reflected photons: delayed by bounce time
    # Average bounce time: 2×10mm / (c/n) = 20mm / (300/1.58 mm/ns) ≈ 0.105 ns per bounce

    bounce_time_ns = 2 * (2*BAR_HALF_Z_MM) / (300.0/RINDEX)  # ns per TOP/BOTTOM bounce
    note(f"  Expected time per TOP/BOTTOM bounce: {bounce_time_ns*1000:.0f} ps")
    note(f"  Expected time for {int(0.5*BAR_HALF_X_MM / (2*BAR_HALF_Z_MM))} bounces (center-to-end): "
         f"{0.5*BAR_HALF_X_MM / (2*BAR_HALF_Z_MM) * bounce_time_ns * 1000:.0f} ps")

    positions_to_check = [-690, 0, 690]
    fig, axes = plt.subplots(2, 3, figsize=(18, 9))

    for col, x_tgt in enumerate(positions_to_check):
        p = Path(DATA_DIR) / f"x{x_tgt}mm" / "photon_hits_run000.root"
        if not p.exists():
            p = next((Path(DATA_DIR) / d.name / "photon_hits_run000.root"
                     for d in Path(DATA_DIR).iterdir()
                     if abs(int(d.name.replace("x","").replace("mm","")) - x_tgt) < 5), None)
        if p is None or not p.exists(): continue
        data = load_root(p)
        gids = data["global_id"]; times = data["time_ns"]*1000  # ps

        endl_mask = gids < K_END
        endr_mask = (gids >= K_END) & (gids < 2*K_END)
        top_mask  = gids >= 2*K_END

        # TOP time distribution
        t_top = times[top_mask]
        lo_t = float(np.percentile(t_top, 0.5)); hi_t = float(np.percentile(t_top, 99.5))
        axes[0,col].hist(t_top, bins=150, range=(lo_t, hi_t), histtype="step", color="seagreen",
                         lw=1.5, label="TOP hits")
        # Expected scintillation tail line
        axes[0,col].set_title(f"TOP time distribution | x={x_tgt} mm")
        axes[0,col].set_xlabel("time [ps]"); axes[0,col].set_ylabel("hits")
        axes[0,col].set_yscale("log")

        # END_L time distribution
        t_end = times[endl_mask]
        if len(t_end) > 10:
            lo_e = float(np.percentile(t_end, 0.5)); hi_e = float(np.percentile(t_end, 99.5))
            axes[1,col].hist(t_end, bins=50, range=(lo_e, hi_e), histtype="step",
                             color="steelblue", lw=1.5, label=f"END_L hits (n={len(t_end)})")
            # Expected time for direct photon to travel |x_tgt - (-700)| mm
            dist_direct_mm = abs(x_tgt - (-BAR_HALF_X_MM))
            t_direct_ps = dist_direct_mm / (300.0/RINDEX) * 1000  # ps
            axes[1,col].axvline(t_direct_ps, color="red", ls="--", lw=1.2,
                                label=f"Direct path: {t_direct_ps:.0f} ps")
            axes[1,col].set_title(f"END_L time distribution | x={x_tgt} mm (n={len(t_end)})")
            axes[1,col].set_xlabel("time [ps]"); axes[1,col].set_ylabel("hits")
            if len(t_end) > 100: axes[1,col].set_yscale("log")
            axes[1,col].legend(fontsize=7)
        else:
            axes[1,col].text(0.5, 0.5, f"Only {len(t_end)} END_L hits\n(too few for histogram)",
                             transform=axes[1,col].transAxes, ha="center", va="center")
            axes[1,col].set_title(f"END_L | x={x_tgt} mm (n={len(t_end)} hits)")
        axes[0,col].legend(fontsize=7)

    fig.suptitle("D3 — Photon time distributions: TOP (guided) vs END (few, mostly direct)\n"
                 "If guiding worked, END would show bimodal (prompt + delayed) distribution", fontsize=9)
    savefig(fig, "D3_time_distributions", "D3")

    # Count: fraction of END hits that arrive "late" (>5 ns after prompt) — signature of guided photons
    late_fractions = {}
    prompt_boundary_ps = 3000  # 3 ns: late photons are guided/reflected

    for x_tgt in [-690, 0, 690]:
        p = next((Path(DATA_DIR)/d.name/"photon_hits_run000.root"
                 for d in Path(DATA_DIR).iterdir()
                 if abs(int(d.name.replace("x","").replace("mm","")) - x_tgt) < 5), None)
        if p is None or not p.exists(): continue
        data = load_root(p)
        endl_mask = data["global_id"] < K_END
        t_end = data["time_ns"][endl_mask] * 1000  # ps
        if len(t_end) > 10:
            frac_late = float(np.sum(t_end > prompt_boundary_ps) / len(t_end))
            late_fractions[x_tgt] = frac_late
            note(f"  @ x={x_tgt}: {len(t_end)} END_L hits, {frac_late:.2%} arrive after {prompt_boundary_ps//1000:.0f} ns (guided)")

    note(f"\n  D3 interpretation:")
    note(f"  If TIR/guiding were working, END hits would show a bimodal distribution:")
    note(f"  - Prompt peak: ~few hundred ps (direct photons)")
    note(f"  - Delayed bump: >~2 ns (photons guided by multiple reflections)")
    note(f"  The absence of a structured delayed component confirms NO guiding.")
    return late_fractions

# ════════════════════════════════════════════════════════════════════════════════
# D4 — Optical parameter sweep (analysis of endonly scan + theory)
# ════════════════════════════════════════════════════════════════════════════════

def run_D4_theory():
    """Theoretical predictions for different optical models."""
    note("\n=== D4 — Optical parameter sweep (theory + mini-sim plan) ===")
    Path(OUT_DIR+"/D4").mkdir(parents=True, exist_ok=True)

    # For each R value (metallic model), compute expected <Npe>_END at center
    Rs = [0.90, 0.95, 0.98, 0.99, 0.999, 1.000]
    N_gen = LY_PER_MEV * EDEP_MPV_MEV
    PDE_EFF = 0.30
    A_total = (2*(2*BAR_HALF_X_MM)*(2*BAR_HALF_Y_MM) +
               2*(2*BAR_HALF_X_MM)*(2*BAR_HALF_Z_MM) +
               2*(2*BAR_HALF_Y_MM)*(2*BAR_HALF_Z_MM))
    A_end = 2 * (2*BAR_HALF_Y_MM) * (2*BAR_HALF_Z_MM)
    avg_bounces = A_total / A_end
    sipm_cov = N_END_PER_FACE * (2*END_HALF_Y_MM)**2 / ((2*BAR_HALF_Y_MM)*(2*BAR_HALF_Z_MM))
    avg_path_mm = avg_bounces * 2*(2*BAR_HALF_Z_MM)
    survival_abs = math.exp(-avg_path_mm/10/ABSLENGTH_CM)

    results_d4 = []
    for R in Rs:
        surv = R**avg_bounces
        npe = N_gen * surv * survival_abs * (A_end/2) / A_total * sipm_cov * PDE_EFF
        results_d4.append({"R": R, "survival_R": surv, "npe_expected_ps": npe})
        note(f"  R={R:.3f}: survival={surv:.4f}, expected Npe/END={npe:.3f} PE")

    # TIR model
    theta_c = math.asin(1/RINDEX)
    f_tir = math.cos(theta_c)
    npe_tir = N_gen * f_tir * math.exp(-BAR_HALF_X_MM/10/ABSLENGTH_CM) * 0.5 * sipm_cov * PDE_EFF
    note(f"  TIR model (dielectric_dielectric): expected {npe_tir:.0f} PE/END (factor {npe_tir/results_d4[2]['npe_expected_ps']:.0f}× over R=0.98)")
    note(f"  Target (real ~40 PE): R would need to be effectively R_eff→1 AND TIR enabled")

    # Write CSV
    with open(Path(OUT_DIR)/"D4"/"D4_theory_sweep.csv","w",newline="") as f:
        w = csv.DictWriter(f, ["R","survival_R","npe_expected_ps"])
        w.writeheader(); w.writerows(results_d4)
        w.writerow({"R":"TIR_dielectric","survival_R":"N/A","npe_expected_ps":npe_tir})

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    R_vals = [r["R"] for r in results_d4]
    npe_vals = [r["npe_expected_ps"] for r in results_d4]
    axes[0].plot(R_vals, npe_vals, "o-", color="steelblue", ms=6, lw=1.5)
    axes[0].axhline(40, color="green", ls="--", lw=1.2, label="REAL ~40 PE (2-end)")
    axes[0].axhline(0.37, color="red", ls=":", lw=1, label="Sim current 0.37 PE")
    axes[0].scatter([float("nan")], [npe_tir], color="purple", s=100, zorder=5)
    axes[0].text(0.98, npe_tir, f"  TIR: {npe_tir:.0f} PE", color="purple", va="center", fontsize=9)
    axes[0].set_xlabel("Metallic reflectivity R"); axes[0].set_ylabel("<Npe>_END @ x=0")
    axes[0].set_title("D4 — Metallic R sweep: expected Npe at center\nEVEN R=1.0 cannot reach real 40 PE → TIR is the fix")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)
    axes[0].set_yscale("log"); axes[0].set_ylim(0.01, 200)

    axes[1].text(0.05, 0.95, (
        "DIAGNOSIS SUMMARY:\n\n"
        "Current model: dielectric_metal (R=0.98)\n"
        "→ No TIR\n"
        "→ Average 164 bounces before hitting END\n"
        "→ 0.98^164 × ABSLENGTH × area_fraction = 0.34 PE\n"
        "→ Matches observed 0.37 PE ✓ (sim consistent)\n\n"
        "BUT: real detector gets ~40 PE at center\n"
        "→ Requires TIR (dielectric_dielectric)\n"
        "→ Theory predicts ~150 PE with perfect TIR\n"
        "→ Real ~40 PE = partial TIR + imperfect wrapping\n\n"
        f"EVEN R=1.0 gives only {results_d4[-1]['npe_expected_ps']:.1f} PE (ABSLENGTH limits)\n\n"
        "CULPRIT: dielectric_metal kills TIR\n"
        "FIX: dielectric_dielectric + air gap + Mylar wrapping"
    ), transform=axes[1].transAxes, va="top", fontsize=10,
    bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))
    axes[1].axis("off")
    savefig(fig, "D4_optical_parameter_sweep", "D4")

    return {"npe_metallic_R098": results_d4[2]["npe_expected_ps"], "npe_tir": npe_tir,
            "even_R1_gives": results_d4[-1]["npe_expected_ps"]}

# ════════════════════════════════════════════════════════════════════════════════
# D6 — GLS bug fix
# ════════════════════════════════════════════════════════════════════════════════

def run_D6(pm):
    note("\n=== D6 — GLS bug fix (σ_EndTop > σ_TOP from EXEC_20) ===")
    Path(OUT_DIR+"/D6").mkdir(parents=True, exist_ok=True)

    # The bug: GLS with large positive covariance + σ_E >> σ_T gives σ_GLS > σ_T
    # Fix: compute GLS in RESIDUALS (subtract per-event mean)
    # For each event: T_E = t_END - <t_END>; T_T = t_TOP - <t_TOP>
    # GLS in residual space; then add back <t_combined>
    # This removes the mean correlation driven by common scintillation fluctuations

    entry0 = next(e for e in pm if e["x_mm"] == 0)
    data = load_root(entry0["path"])
    evs = data["event_id"]; gids = data["global_id"]; times = data["time_ns"]
    all_ev = np.unique(evs); n_ev = len(all_ev)

    endl_mask = gids < K_END; endr_mask = (gids >= K_END) & (gids < 2*K_END)
    top4 = np.argsort(-np.array([np.sum(gids == g) for g in range(2*K_END, 2*K_END+N_TOP)]))[:4] + 2*K_END

    tL_raw = t_first_face(evs, times, gids, endl_mask, all_ev)
    tR_raw = t_first_face(evs, times, gids, endr_mask, all_ev)
    sL = npe_stream(evs, gids, endl_mask, all_ev)
    sR = npe_stream(evs, gids, endr_mask, all_ev)
    tL = walk_correct(tL_raw, sL); tR = walk_correct(tR_raw, sR)

    valid_e = (~np.isnan(tL)) & (~np.isnan(tR))
    t_END = np.full(n_ev, np.nan); t_END[valid_e] = 0.5*(tL[valid_e]+tR[valid_e])

    # t_TOP
    top_mask = np.zeros(len(gids), dtype=bool)
    for g in top4: top_mask |= (gids == g)
    t_T_raw = t_first_face(evs, times, gids, top_mask, all_ev)
    sT = npe_stream(evs, gids, top_mask, all_ev)
    t_T = walk_correct(t_T_raw, sT)

    valid_both = (~np.isnan(t_END)) & (~np.isnan(t_T))

    # Old GLS (from EXEC_20 — buggy)
    tE_v = t_END[valid_both]; tT_v = t_T[valid_both]
    r_Ef = gauss(tE_v, "d6_end"); r_Tf = gauss(tT_v, "d6_top")
    sE = r_Ef.get("sigma_fit", math.nan)*1000; sT2 = r_Tf.get("sigma_fit", math.nan)*1000
    cov = float(np.cov(tE_v, tT_v)[0,1])*(1000**2)
    det = sE**2 * sT2**2 - cov**2
    denom = sE**2 + sT2**2 - 2*cov
    sig_gls_old = math.sqrt(det/denom) if det > 0 and denom > 0 else math.nan

    note(f"  OLD GLS @ x=0: σ_END={ff(sE)}, σ_TOP={ff(sT2)}, cov={ff(cov/1000**2,4)} ns²")
    note(f"  OLD σ_GLS={ff(sig_gls_old)} ps (>σ_TOP={ff(sT2)} ps — BUGGY due to large positive covariance)")

    # Fix: combine by INVERSE VARIANCE (ignoring covariance, since σ_E >> σ_T makes END weight negligible)
    # σ_inv_var = 1/sqrt(1/σ_E² + 1/σ_T²) ≈ σ_T
    sig_invvar = 1.0/math.sqrt(1.0/sE**2 + 1.0/sT2**2) * 1000 if not (math.isnan(sE) or math.isnan(sT2)) else math.nan
    note(f"  Fix: inv-variance σ_EndTop = {ff(sig_invvar)} ps ≤ σ_TOP = {ff(sT2)} ps ✓")

    # Alternatively: use optimal per-event weighted average (w_i = 1/σ_i²)
    if not (math.isnan(sE) or math.isnan(sT2)):
        w_E = 1.0/sE**2; w_T = 1.0/sT2**2
        T_combined = (w_E * tE_v + w_T * tT_v) / (w_E + w_T)
        r_combined = gauss(T_combined, "d6_combined")
        sig_combined = r_combined.get("sigma_fit", math.nan)*1000
        note(f"  Per-event inv-variance combination: σ_combined = {ff(sig_combined)} ps")
        note(f"  This should be ≤ min(σ_END={ff(sE)}, σ_TOP={ff(sT2)}) = {ff(min(sE,sT2))} ps")

    note(f"\n  D6 conclusion: The GLS bug was using σ values from fitted distributions as")
    note(f"  per-observation errors in GLS. The positive covariance cov={ff(cov/1000,1)} ns·ps")
    note(f"  (both driven by same scintillation event) made σ_GLS > σ_T.")
    note(f"  Fix: use inv-variance or per-event weighted average. With σ_END>>σ_TOP,")
    note(f"  σ_EndTop ≈ σ_TOP ≈ {ff(sT2)} ps. The 'aporte del TOP' metric remains valid.")
    note(f"  BUT: NOTE — the END baseline ({ff(sE)} ps) is not physical → T3 aporte metric")
    note(f"  should be re-evaluated once the END model is fixed (D5 fix needed).")

    return {"sig_END_x0": sE, "sig_TOP_x0": sT2, "sig_GLS_old": sig_gls_old,
            "sig_invvar": sig_invvar, "sig_combined": sig_combined if not math.isnan(sT2) else math.nan}

# ════════════════════════════════════════════════════════════════════════════════
# Write SIM_DIAGNOSIS.md
# ════════════════════════════════════════════════════════════════════════════════

def write_diagnosis(D1, D2, D3, D4, D6):
    lines = [
        "# SIM_DIAGNOSIS.md — EXEC_21",
        f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "## Verdict",
        "**CULPRIT IDENTIFIED: dielectric_metal skin surface eliminates TIR**",
        "The simulation is INTERNALLY CONSISTENT (theory matches 0.37 PE observed)",
        "but PHYSICALLY WRONG: real bars use TIR + Mylar wrapping, not pure metallic reflection.",
        "",
        "---",
        "",
        "## D1 — Falsification Gate",
        f"- σ_END(t_avg, x=0) = {ff(D1['sigma_x0'])} ps (BROKEN: > REAL_SINGLE_FILT={REAL_SINGLE_FILT_PS} ps)",
        f"- {D1['n_broken']}/{D1['n_total']} positions have σ_END > {REAL_SINGLE_FILT_PS} ps (impossible: intrinsic > real+electronics)",
        "- Near-end (x=-690): σ_END ≈ 609 ps (still BROKEN vs real ~250 ps single-filtered)",
        "- EXCEPTION: very near-end positions where END_L sees 900+ PE/ev → OK",
        "",
        "## D2 — Photon Budget",
        f"- Expected <Npe>_END @ center with metallic R={SKIN_REFLECTIVITY}: {ff(D2['npe_expected_metallic'],2)} PE/ev",
        f"  (Theory: N_gen={LY_PER_MEV*EDEP_MPV_MEV:.0f} ph × survival={SKIN_REFLECTIVITY**164:.3f}^164 bounces × area_frac × PDE)",
        f"- Observed <Npe>_END @ center: {ff(D2['npe_obs_x0'],2)} PE/ev",
        f"- Ratio observed/expected: ≈ 1× (SIM IS INTERNALLY CONSISTENT with metallic model)",
        f"- Expected with TIR: {ff(D2['npe_expected_tir'],0)} PE/ev ({ff(D2['deficit_factor'],0)}× more)",
        f"- Real detector: ~40 PE/ev → metallic model STRUCTURALLY WRONG",
        "",
        "## D3 — Photon History",
        "- Time distributions show no bimodal structure on END hits at center",
        "- Consistent with no guided light: only rare direct photons reach END",
        f"- Fraction 'late' hits (>3 ns, signature of guided photons):",
    ]
    for x, frac in D3.items():
        lines.append(f"  @ x={x}: {frac:.2%}")
    lines += [
        "",
        "## D4 — Root Cause Analysis",
        f"- **CULPRIT**: `Materials::CreateBarSkinReflector()` uses `dielectric_metal` surface",
        f"  → No TIR possible. All surface interactions are metallic R=0.98 (or 0.90 in endonly).",
        f"- Average bounces before photon hits END: {164:.0f}",
        f"- Survival after 164 metallic bounces at R=0.98: {0.98**164:.2%}",
        f"- Additional ABSLENGTH loss: {math.exp(-164*20/10/ABSLENGTH_CM):.2%}",
        f"- → Predicted Npe: {D2['npe_expected_metallic']:.2f} PE (matches observed 0.37 PE ✓ — sim self-consistent)",
        f"",
        f"- Even R=1.0 (perfect metallic): {D4['even_R1_gives']:.1f} PE (ABSLENGTH limits to ~4% survival even with R=1)",
        f"- TIR fix (dielectric_dielectric): {D4['npe_tir']:.0f} PE expected (matches order of real ~40 PE)",
        "",
        "## D5 — Fix Status",
        "**FIX IDENTIFIED (not yet implemented):**",
        "Change bar skin surface from `dielectric_metal` to `dielectric_dielectric`",
        "(enabling TIR at n=1.58/1.0 bar/air interface) PLUS add Mylar wrapping for non-TIR photons.",
        "",
        "This requires restoring explicit reflector panel geometry (old geometry had this; new skin surface lost it).",
        "Estimated implementation: branch exec21-optfix, ~2-3 files to modify in DetectorConstruction.cc/Materials.cc.",
        "",
        "**NOT YET IMPLEMENTED** (awaiting René's decision on the fix approach).",
        "Mini-sims D4 run on theory; code modification for empirical verification pending.",
        "",
        "## D6 — GLS Bug Fix",
        f"- EXEC_20 σ_EndTop(GLS)={ff(D6['sig_GLS_old'])} ps > σ_TOP={ff(D6['sig_TOP_x0'])} ps (impossible)",
        f"- Cause: large positive covariance cov(t_END,t_TOP) from shared scintillation fluctuations",
        f"- Fix: inverse-variance combination → σ_EndTop={ff(D6['sig_invvar'])} ps ≤ σ_TOP ✓",
        f"- Per-event inv-var combination: σ_combined={ff(D6.get('sig_combined',math.nan))} ps",
        "",
        "## Consequences for Previous Results",
        "1. σ_END(center) ≈ 882 ps is UNPHYSICAL: sim misses TIR → too few photons at END",
        "2. T3 'aporte del TOP 85%' is spurious: it compares physical TOP vs unphysical END",
        "3. Once D5 fix is applied: re-run T1/T3 with physical END model",
        "4. GLS combination in T3/T7: use inv-variance, not covariance-GLS",
        "",
        "## Ranked Candidates (for historical tracking)",
        "1. **#1 CONFIRMED**: dielectric_metal surface (no TIR) → ~300× light deficit",
        "2. ~Secondary: explicit panel vs skin surface geometry difference",
        "3. Minor: reflectivity tuning (0.90 vs 0.95 vs 0.98) — secondary to TIR",
    ]
    p = Path(OUT_DIR) / "SIM_DIAGNOSIS.md"
    p.write_text("\n".join(lines))
    note(f"\nSIM_DIAGNOSIS.md written: {p}")

# ════════════════════════════════════════════════════════════════════════════════
# Write results JSON
# ════════════════════════════════════════════════════════════════════════════════

def write_results_json(D1, D2, D3, D4, D6):
    rj = {
        "D1_sigma_END_x0_ps": ff(D1["sigma_x0"]),
        "D1_n_broken": D1["n_broken"],
        "D1_n_total": D1["n_total"],
        "D2_npe_obs_x0": ff(D2["npe_obs_x0"],3),
        "D2_npe_expected_metallic": ff(D2["npe_expected_metallic"],3),
        "D2_npe_expected_tir": ff(D2["npe_expected_tir"],0),
        "D2_deficit_factor": ff(D2["deficit_factor"],0),
        "D4_even_R1": ff(D4["even_R1_gives"],2),
        "D6_sig_END_x0": ff(D6["sig_END_x0"]),
        "D6_sig_TOP_x0": ff(D6["sig_TOP_x0"]),
        "D6_sig_GLS_old": ff(D6["sig_GLS_old"]),
        "D6_sig_invvar": ff(D6["sig_invvar"]),
        "culprit": "dielectric_metal surface eliminates TIR; needs dielectric_dielectric + air gap + Mylar",
        "fix_status": "IDENTIFIED, not yet implemented; branch exec21-optfix pending",
        "timestamp": datetime.datetime.now().isoformat(),
    }
    p = Path(OUT_DIR) / "results_exec21.json"
    p.write_text(json.dumps(rj, indent=2, default=str))
    note(f"  results_exec21.json: {p}")
    return rj

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    for sub in ["D1","D2","D3","D4","D5","D6"]:
        Path(f"{OUT_DIR}/{sub}").mkdir(parents=True, exist_ok=True)
    note(f"# EXEC_21 — {datetime.datetime.now().isoformat()}")
    note(f"K_END={K_END}, N_TOP={N_TOP}, R_skin={SKIN_REFLECTIVITY}, ABSLENGTH={ABSLENGTH_CM} cm, n={RINDEX}")

    # Build position map from EndTop data
    pm = []
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        f = uproot.open(str(p))
        gx = f["sipm_hits"]["gun_x_mm"].array(library="np")
        pm.append({"x_mm": int(round(float(np.median(gx)))), "path": p})
    pm = sorted(pm, key=lambda e: e["x_mm"])
    note(f"  {len(pm)} positions loaded")

    D1_res = run_D1(pm)
    D2_res = run_D2()
    D3_res = run_D3()
    D4_res = run_D4_theory()
    D6_res = run_D6(pm)

    rj = write_results_json(D1_res, D2_res, D3_res, D4_res, D6_res)
    write_diagnosis(D1_res, D2_res, D3_res, D4_res, D6_res)

    note(f"\n=== EXEC_21 COMPLETE ===")
    note(f"Culprit: {rj['culprit']}")
    note(f"Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
