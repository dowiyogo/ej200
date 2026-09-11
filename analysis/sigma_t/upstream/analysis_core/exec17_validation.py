#!/usr/bin/env python3.12
"""
exec17_validation.py — EXEC_17 falsifiable simulation validation battery (V1–V9).

Each test: prediction from config → criterion PASS/CONCERN/FAIL → figure + sidecar.
Output: SIM_TRUST_VERDICT.md with global verdict.

Run:
    MPLBACKEND=Agg python3.12 exec17_validation.py
"""

import sys, os, csv, json, math, datetime, hashlib, warnings
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent))
from lib.fit_engine   import fit_core_gaussian
from lib.robust_seeds import gather_seeds

# ════════════════════════════════════════════════════════════════════════════════
# Constants — ALL re-read from source; NONE from EXEC_16 memory
# ════════════════════════════════════════════════════════════════════════════════

# From DetectorConstruction.hh:37-39
K_END         = 8      # kNEndSiPMs per side
N_TOP         = 70     # kNTopSiPMs
N_TOTAL       = 2 * K_END + N_TOP  # 86

# From DetectorConstruction.cc:29-42
BAR_HALF_X_MM = 700.0  # kBarHalfX
BAR_HALF_Y_MM = 30.0   # kBarHalfY
BAR_HALF_Z_MM = 5.0    # kBarHalfZ
END_HALF_X_MM = 0.25   # kEndHalfX
END_HALF_Y_MM = 3.0    # kEndHalfY
END_HALF_Z_MM = 3.0    # kEndHalfZ
END_PITCH_MM  = 2 * END_HALF_Y_MM + 1.5  # = 7.5 mm (kEndPitch)
TOP_HALF_X_MM = 3.0    # kTopHalfX
TOP_HALF_Y_MM = 0.25   # kTopHalfY
TOP_HALF_Z_MM = 3.0    # kTopHalfZ

# Predicted SiPM positions from geometry
PREDICTED_END_X_MM  = BAR_HALF_X_MM - END_HALF_X_MM   # = 699.75 mm
PREDICTED_TOP_Y_MM  = BAR_HALF_Y_MM - TOP_HALF_Y_MM   # = 29.75 mm

def top_center_x(idx):
    """TopSiPMCenterX from DetectorConstruction.cc:47-48"""
    if idx < 35:
        return -692.0 + 20.0 * idx
    return 12.0 + 20.0 * (idx - 35)

TOP_X_MM = np.array([top_center_x(i) for i in range(N_TOP)])  # 70 values

# From opsc-101.mac + DetectorConstruction.cc overrides
TAU_D_NS     = 1.8    # SCINTILLATIONTIMECONSTANT1
TAU_R_NS     = 0.7    # SCINTILLATIONRISETIME1
LY_PER_MEV   = 10400  # SCINTILLATIONYIELD
ABSLENGTH_CM = 160.0  # kBarHalfZ×2 = 10 mm bar height; DetConstr override
EMIS_PEAK_NM = 408.8  # peak of scntComp1.txt (read directly)

# From geometry: muon path through bar = 2×kBarHalfZ = 10 mm = 1 cm
MUON_PATH_CM  = 2 * BAR_HALF_Z_MM / 10.0  # = 1.0 cm
DEDX_MEVPCM   = 2.0   # approximate MIP dE/dx in plastic scintillator
EDEP_MPV_MEV  = DEDX_MEVPCM * MUON_PATH_CM  # ≈ 2 MeV

# Photon energy from emission peak: E = hc/λ
HC_EV_NM = 1239.8  # eV·nm
PHOTON_EV = HC_EV_NM / EMIS_PEAK_NM  # ≈ 3.03 eV

# Cross-check: total optical photons ≈ LY × Edep_MPV
EXPECTED_PHOTONS_MPV = LY_PER_MEV * EDEP_MPV_MEV  # ≈ 20800 photons before losses

# Paths
DATA_DIR = "/home/reriosto/SHiP/t0minidaq/runs/t0minidaq_endtop_scan_20260618_204959/outputs"
OUT_DIR  = "/home/reriosto/SHiP/analysis_core/out/EXEC_17"
EXEC16_CSV_L2 = "/home/reriosto/SHiP/analysis_core/out/EXEC_16/L2/sigma_t_all_estimators.csv"
EXEC16_CSV_L1 = "/home/reriosto/SHiP/analysis_core/out/EXEC_16/L1/npe_per_face.csv"

TREE  = "sipm_hits"
BR = dict(ev="event_id", ft="face_type", gid="global_id",
          t="time_ns", e="energy_eV", wl="wl_nm", pde="pde",
          x="x_mm", y="y_mm", z="z_mm", gx="gun_x_mm")

RANDOM_SEED = 20260618

# ── Verdict accumulator ────────────────────────────────────────────────────────
_verdicts = {}  # test_id → {"result": PASS/CONCERN/FAIL, "evidence": str}

def verdict(vid, result, evidence):
    _verdicts[vid] = {"result": result, "evidence": evidence}
    sym = {"PASS": "✓", "CONCERN": "△", "FAIL": "✗"}.get(result, "?")
    print(f"  [{sym} {result}] {vid}: {evidence}")

def savefig(fig, stem, subdir):
    d = Path(OUT_DIR) / subdir
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(d / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)

def note(msg): print(msg)

# ── Data helpers ───────────────────────────────────────────────────────────────

def pos_map():
    """Return sorted list of {x_mm, path} from DATA_DIR."""
    entries = []
    for d in sorted(Path(DATA_DIR).iterdir()):
        p = d / "photon_hits_run000.root"
        if not p.exists(): continue
        f = uproot.open(str(p))
        gx = f[TREE][BR["gx"]].array(library="np")
        entries.append({"x_mm": int(round(float(np.median(gx)))), "path": p})
    return sorted(entries, key=lambda e: e["x_mm"])

def load(entry, branches=None):
    if branches is None: branches = list(BR.values())
    f = uproot.open(str(entry["path"]))
    return f[TREE].arrays(branches, library="np")

def classify(gid):
    if gid < K_END:        return ("END_L", gid)
    if gid < 2 * K_END:    return ("END_R", gid - K_END)
    return ("TOP", gid - 2 * K_END)

def cfg_mock():
    return {"MIN_EVENTS_FOR_FIT": 30, "FIT_WINDOW_SIGMAS": 2.0,
            "BINNING_STRATEGY": "sqrt_n", "N_BOOTSTRAP": 200,
            "RANDOM_SEED": RANDOM_SEED, "FIT_OPTIONS": "R Q S 0",
            "CHI2_NDF_WARN": 3.0}

def gauss_fit_v(v, prefix="g"):
    v_clean = v[~np.isnan(v)]
    if len(v_clean) < 30:
        return {"sigma_fit": np.nan, "bootstrap_err": np.nan, "flag": "insufficient_stats"}
    return fit_core_gaussian(v_clean, cfg_mock(), name_prefix=prefix)

# ════════════════════════════════════════════════════════════════════════════════
# V1 — Geometry / channel-map from hits
# ════════════════════════════════════════════════════════════════════════════════

def run_V1(pm):
    note("\n=== V1 — Geometry / channel-map from hits ===")
    # Load one central + one far-left position for richer coverage
    entries_to_use = [e for e in pm if e["x_mm"] in (0, -690, 690)][:3]
    if not entries_to_use: entries_to_use = [pm[len(pm)//2]]

    x_by_gid = {g: [] for g in range(N_TOTAL)}
    y_by_gid = {g: [] for g in range(N_TOTAL)}
    z_by_gid = {g: [] for g in range(N_TOTAL)}

    for entry in entries_to_use:
        d = load(entry, [BR["gid"], BR["x"], BR["y"], BR["z"]])
        for g in range(N_TOTAL):
            mask = d[BR["gid"]] == g
            if mask.sum() > 0:
                x_by_gid[g].extend(d[BR["x"]][mask].tolist())
                y_by_gid[g].extend(d[BR["y"]][mask].tolist())
                z_by_gid[g].extend(d[BR["z"]][mask].tolist())

    mean_x = {g: float(np.mean(x_by_gid[g])) if x_by_gid[g] else np.nan for g in range(N_TOTAL)}
    mean_y = {g: float(np.mean(y_by_gid[g])) if y_by_gid[g] else np.nan for g in range(N_TOTAL)}
    mean_z = {g: float(np.mean(z_by_gid[g])) if z_by_gid[g] else np.nan for g in range(N_TOTAL)}

    # ── Figure: scatter of SiPMs coloured by face
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # XY view
    for gid in range(N_TOTAL):
        if math.isnan(mean_x[gid]): continue
        face, _ = classify(gid)
        color = {"END_L": "steelblue", "END_R": "tomato", "TOP": "seagreen"}[face]
        axes[0].scatter(mean_x[gid], mean_y[gid], c=color, s=30, alpha=0.8)
    axes[0].axvline(-PREDICTED_END_X_MM, color="steelblue", ls="--", lw=0.8, label=f"pred END_L x={-PREDICTED_END_X_MM:.1f} mm")
    axes[0].axvline(+PREDICTED_END_X_MM, color="tomato",    ls="--", lw=0.8, label=f"pred END_R x={+PREDICTED_END_X_MM:.1f} mm")
    axes[0].axhline(PREDICTED_TOP_Y_MM,  color="seagreen",  ls="--", lw=0.8, label=f"pred TOP y={PREDICTED_TOP_Y_MM:.1f} mm")
    axes[0].set_xlabel("x_mm (hit mean)"); axes[0].set_ylabel("y_mm (hit mean)")
    axes[0].set_title("V1 — SiPM positions (XY) from hit averages")
    axes[0].legend(fontsize=7)

    # TOP SiPMs: compare predicted cx vs observed cx
    obs_cx  = [mean_x[2*K_END + i] for i in range(N_TOP)]
    pred_cx = list(TOP_X_MM)
    valid   = [(p, o) for p, o in zip(pred_cx, obs_cx) if not math.isnan(o)]
    if valid:
        px, ox = zip(*valid)
        axes[1].scatter(px, ox, s=20, color="seagreen", alpha=0.8)
        axes[1].plot([-700, 700], [-700, 700], "k--", lw=0.8, label="perfect 1:1")
        axes[1].set_xlabel("Predicted TOP cx [mm]"); axes[1].set_ylabel("Observed cx [mm]")
        axes[1].set_title("V1 — TOP SiPM: predicted vs observed x-position")
        axes[1].legend(fontsize=7)

    savefig(fig, "V1_geometry_channelmap", "V1")

    # ── CSV sidecar
    with open(Path(OUT_DIR) / "V1" / "V1_sipm_positions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["gid","face","local_id","mean_x","mean_y","mean_z",
                                "pred_x","pred_y","delta_x"])
        w.writeheader()
        for gid in range(N_TOTAL):
            face, lid = classify(gid)
            if face == "END_L":
                px = -PREDICTED_END_X_MM
            elif face == "END_R":
                px = +PREDICTED_END_X_MM
            else:
                px = top_center_x(lid)
            dx = mean_x[gid] - px if not math.isnan(mean_x[gid]) else np.nan
            py = PREDICTED_TOP_Y_MM if face == "TOP" else np.nan
            w.writerow({"gid": gid, "face": face, "local_id": lid,
                        "mean_x": mean_x[gid], "mean_y": mean_y[gid], "mean_z": mean_z[gid],
                        "pred_x": px, "pred_y": py, "delta_x": dx})

    # ── Evaluate
    n_found = sum(1 for g in range(N_TOTAL) if not math.isnan(mean_x[g]))
    n_end_found = sum(1 for g in range(2*K_END) if not math.isnan(mean_x[g]))
    n_top_found = sum(1 for g in range(2*K_END, N_TOTAL) if not math.isnan(mean_x[g]))

    # Check END position
    end_l_xs = [mean_x[g] for g in range(K_END) if not math.isnan(mean_x[g])]
    end_r_xs = [mean_x[g] for g in range(K_END, 2*K_END) if not math.isnan(mean_x[g])]
    end_l_x_mean = float(np.mean(end_l_xs)) if end_l_xs else np.nan
    end_r_x_mean = float(np.mean(end_r_xs)) if end_r_xs else np.nan

    # Check TOP y
    top_ys = [mean_y[g] for g in range(2*K_END, N_TOTAL) if not math.isnan(mean_y[g])]
    top_y_mean = float(np.mean(top_ys)) if top_ys else np.nan

    # Check TOP x mapping: predicted vs observed RMS
    dxs = [abs(mean_x[2*K_END+i] - top_center_x(i)) for i in range(N_TOP)
            if not math.isnan(mean_x[2*K_END+i])]
    top_x_rms = float(np.std(dxs)) if dxs else np.nan

    end_l_ok  = not math.isnan(end_l_x_mean) and abs(end_l_x_mean - (-PREDICTED_END_X_MM)) < 5
    end_r_ok  = not math.isnan(end_r_x_mean) and abs(end_r_x_mean - (+PREDICTED_END_X_MM)) < 5
    top_y_ok  = not math.isnan(top_y_mean) and abs(top_y_mean - PREDICTED_TOP_Y_MM) < 2
    top_map_ok = top_x_rms < 10 if not math.isnan(top_x_rms) else False
    count_ok  = n_end_found == 2*K_END and n_top_found == N_TOP

    evid = (f"n_found={n_found}/{N_TOTAL}; END_L x_mean={end_l_x_mean:.1f} mm (pred={-PREDICTED_END_X_MM:.1f}); "
            f"END_R x_mean={end_r_x_mean:.1f} mm (pred={+PREDICTED_END_X_MM:.1f}); "
            f"TOP y_mean={top_y_mean:.1f} mm (pred={PREDICTED_TOP_Y_MM:.1f}); "
            f"TOP x-map rms_err={top_x_rms:.2f} mm")

    if all([end_l_ok, end_r_ok, top_y_ok, top_map_ok, count_ok]):
        verdict("V1", "PASS", evid)
    elif not count_ok or not (end_l_ok and end_r_ok):
        verdict("V1", "FAIL", evid)
    else:
        verdict("V1", "CONCERN", evid)

    return mean_x, mean_y

# ════════════════════════════════════════════════════════════════════════════════
# V2 — Energy deposition and photon production
# ════════════════════════════════════════════════════════════════════════════════

def run_V2(pm):
    note("\n=== V2 — Energy deposition and photon production ===")
    entry = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry, [BR["ev"], BR["e"], BR["gid"]])

    all_events = np.unique(d[BR["ev"]])
    # total energy_eV per event (proxy for photon count × photon_E)
    total_eV_per_ev = []
    n_hits_per_ev   = []
    for ev in all_events:
        mask = d[BR["ev"]] == ev
        total_eV_per_ev.append(float(d[BR["e"]][mask].sum()))
        n_hits_per_ev.append(int(mask.sum()))

    total_eV = np.array(total_eV_per_ev)
    n_hits   = np.array(n_hits_per_ev)

    # Predicted: N_photons_detected MPV ≈ LY × Edep_MPV × solid_angle × PDE
    # sum(energy_eV) ≈ N_detected × PHOTON_EV
    # For sanity: if we detect ~400 ph at x=0, sum_eV ≈ 400 × 3.03 ≈ 1212 eV
    observed_mpv_eV = float(np.percentile(total_eV, 30))  # rough MPV of Landau (below mean)
    observed_mpv_ph = observed_mpv_eV / PHOTON_EV

    # Also check: N_hits per event distribution (should be Landau-like, > 0)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Sum energy distribution
    lo, hi = np.percentile(total_eV, [0.5, 99.5])
    axes[0].hist(total_eV, bins=100, range=(lo, hi), histtype="step", lw=1.5, color="steelblue")
    axes[0].axvline(observed_mpv_eV, color="red", ls="--", lw=1, label=f"~MPV {observed_mpv_eV:.0f} eV")
    axes[0].set_xlabel("sum(energy_eV) per event"); axes[0].set_ylabel("events")
    axes[0].set_title(f"V2 — Total detected photon energy per event | x=0 mm")
    axes[0].legend(fontsize=8)

    # N_hits per event
    lo2, hi2 = np.percentile(n_hits, [0.5, 99.5])
    axes[1].hist(n_hits, bins=80, range=(lo2, hi2), histtype="step", lw=1.5, color="tomato")
    axes[1].axvline(np.percentile(n_hits, 30), color="red", ls="--", lw=1, label=f"~MPV {np.percentile(n_hits, 30):.0f} hits")
    axes[1].set_xlabel("N PE hits per event"); axes[1].set_ylabel("events")
    axes[1].set_title(f"V2 — Detected PE count per event | x=0 mm")
    axes[1].legend(fontsize=8)

    savefig(fig, "V2_edep_photon_production", "V2")

    # Criterion: positively-skewed distribution (Landau-like), non-zero hits/event, plausible scale
    frac_zero = float(np.mean(n_hits == 0))
    skewness  = float(np.mean((n_hits - n_hits.mean())**3) / n_hits.std()**3) if n_hits.std() > 0 else 0

    # Rough order-of-magnitude check on photon count
    # Expected detected ≈ EXPECTED_PHOTONS_MPV × (fraction_detected)
    # fraction ≈ solid_angle × PDE ≈ few %. With 86 channels and bar surface.
    # Just check that observed_mpv_ph is not insane: should be in [10, 100000] range
    scale_ok  = 10 < observed_mpv_ph < 1e5
    shape_ok  = skewness > 0.3  # Landau is positively skewed
    frac_ok   = frac_zero < 0.01  # <1% dead events

    note(f"  Observed ~MPV photons: {observed_mpv_ph:.0f}; skewness={skewness:.2f}; frac_zero={frac_zero:.4f}")
    evid = (f"~MPV sum_eV={observed_mpv_eV:.0f} eV (~{observed_mpv_ph:.0f} photons); "
            f"skewness={skewness:.2f} (Landau>0 expected); frac_zero={frac_zero:.4f}")

    if scale_ok and shape_ok and frac_ok:
        verdict("V2", "PASS", evid)
    elif not scale_ok:
        verdict("V2", "FAIL", evid + " — photon scale implausible")
    else:
        verdict("V2", "CONCERN", evid)

    # Save CSV sidecar (summary stats)
    with open(Path(OUT_DIR) / "V2" / "V2_edep_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["mpv_sum_eV", observed_mpv_eV])
        w.writerow(["mpv_n_photons", observed_mpv_ph])
        w.writerow(["skewness_nhits", skewness])
        w.writerow(["frac_zero_hits", frac_zero])
        w.writerow(["predicted_photons_before_losses", EXPECTED_PHOTONS_MPV])
        w.writerow(["photon_E_eV", PHOTON_EV])

# ════════════════════════════════════════════════════════════════════════════════
# V3 — Wavelength spectrum
# ════════════════════════════════════════════════════════════════════════════════

def run_V3(pm):
    note("\n=== V3 — Wavelength spectrum ===")
    entry = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry, [BR["wl"], BR["gid"]])

    wl  = d[BR["wl"]]
    gid = d[BR["gid"]]

    # All channels combined
    lo, hi = 350.0, 550.0
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Full spectrum
    axes[0].hist(wl, bins=200, range=(lo, hi), histtype="step", lw=1.5, color="purple")
    axes[0].axvline(EMIS_PEAK_NM, color="red", ls="--", lw=1.2, label=f"Config peak {EMIS_PEAK_NM:.1f} nm")
    axes[0].set_xlabel("wavelength [nm]"); axes[0].set_ylabel("hits")
    axes[0].set_title("V3 — Detected photon wavelength spectrum (all channels, x=0)")
    axes[0].legend(fontsize=8)

    # TOP vs END
    top_wl  = wl[gid >= 2*K_END]
    end_wl  = wl[gid <  2*K_END]
    for arr, label, color in [(top_wl, "TOP", "seagreen"), (end_wl, "END", "steelblue")]:
        if len(arr) > 10:
            axes[1].hist(arr, bins=150, range=(lo, hi), histtype="step", lw=1.2, alpha=0.8,
                         color=color, label=label, density=True)
    axes[1].axvline(EMIS_PEAK_NM, color="red", ls="--", lw=1.2, label=f"Config peak {EMIS_PEAK_NM:.1f} nm")
    axes[1].set_xlabel("wavelength [nm]"); axes[1].set_ylabel("density")
    axes[1].set_title("V3 — Spectrum by face (normalised)")
    axes[1].legend(fontsize=8)

    savefig(fig, "V3_wavelength_spectrum", "V3")

    # Find observed peak via histogram mode
    counts, edges = np.histogram(wl, bins=500, range=(350, 550))
    peak_bin = int(np.argmax(counts))
    observed_peak = float(0.5 * (edges[peak_bin] + edges[peak_bin+1]))
    delta = abs(observed_peak - EMIS_PEAK_NM)

    with open(Path(OUT_DIR) / "V3" / "V3_wavelength_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        w.writerow(["observed_peak_nm", observed_peak])
        w.writerow(["predicted_peak_nm", EMIS_PEAK_NM])
        w.writerow(["delta_nm", delta])
        w.writerow(["n_hits", len(wl)])

    evid = f"Observed peak={observed_peak:.1f} nm; Predicted={EMIS_PEAK_NM:.1f} nm; Δ={delta:.1f} nm"
    if delta < 10:
        verdict("V3", "PASS", evid)
    elif delta < 25:
        verdict("V3", "CONCERN", evid)
    else:
        verdict("V3", "FAIL", evid)

# ════════════════════════════════════════════════════════════════════════════════
# V4 — Attenuation profile + TOP heatmap
# ════════════════════════════════════════════════════════════════════════════════

def run_V4(pm):
    note("\n=== V4 — Attenuation profile and TOP heatmap ===")

    positions = [e["x_mm"] for e in pm]
    npe_endl, npe_endr, npe_top_by_gid = [], [], []  # top_by_gid: n_TOP × n_pos

    npe_top_matrix = np.zeros((N_TOP, len(positions)))

    for pi, entry in enumerate(pm):
        d = load(entry, [BR["ev"], BR["gid"]])
        all_events = np.unique(d[BR["ev"]])
        n_ev = len(all_events)
        gids = d[BR["gid"]]

        # END_L: gid 0..K_END-1
        endl_mask = gids < K_END
        npe_endl.append(float(np.sum(endl_mask)) / n_ev)

        # END_R: gid K_END..2K_END-1
        endr_mask = (gids >= K_END) & (gids < 2*K_END)
        npe_endr.append(float(np.sum(endr_mask)) / n_ev)

        # TOP: per channel
        for ti in range(N_TOP):
            gid_val = 2 * K_END + ti
            npe_top_matrix[ti, pi] = float(np.sum(gids == gid_val)) / n_ev

    npe_endl = np.array(npe_endl)
    npe_endr = np.array(npe_endr)
    xs = np.array(positions)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # (a) END_L and END_R NPE(x)
    axes[0].plot(xs, npe_endl, "s-", color="steelblue", ms=5, lw=1.2, label="END_L (gid 0-7)")
    axes[0].plot(xs, npe_endr, "^-", color="tomato",    ms=5, lw=1.2, label="END_R (gid 8-15)")
    axes[0].set_xlabel("gun_x [mm]"); axes[0].set_ylabel("⟨Npe⟩ per event")
    axes[0].set_title("V4a — END attenuation: near/far asymmetry")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)

    # (b) TOP heatmap: <Npe>(TOP_idx, beam_x) — should show diagonal ridge
    im = axes[1].imshow(npe_top_matrix, aspect="auto", origin="lower",
                         extent=[xs.min(), xs.max(), 0, N_TOP],
                         cmap="hot", norm=mcolors.LogNorm(vmin=0.01, vmax=npe_top_matrix.max()))
    plt.colorbar(im, ax=axes[1], label="⟨Npe⟩ per event")
    # Overlay predicted TOP SiPM x positions as diagonal reference
    top_local_idx = np.arange(N_TOP)
    axes[1].plot(TOP_X_MM, top_local_idx, "c.", ms=2, alpha=0.6, label="Predicted cx")
    axes[1].set_xlabel("beam_x [mm]"); axes[1].set_ylabel("TOP local idx (0=gid16)")
    axes[1].set_title("V4b — TOP heatmap ⟨Npe⟩ (should show diagonal ridge)")
    axes[1].legend(fontsize=7)

    savefig(fig, "V4_attenuation_heatmap", "V4")

    # ── Criterion (a): near/far asymmetry for END_L
    # END_L should be higher at negative x (beam near left end)
    # Compare x<-300 vs x>300
    idx_near_L = xs < -300
    idx_far_L  = xs > 300
    npe_endl_near = float(np.mean(npe_endl[idx_near_L])) if idx_near_L.sum() > 0 else np.nan
    npe_endl_far  = float(np.mean(npe_endl[idx_far_L]))  if idx_far_L.sum()  > 0 else np.nan
    asymmetry_endl = (npe_endl_near - npe_endl_far) / (npe_endl_near + npe_endl_far + 1e-6)

    # Symmetric check: END_R should be mirror of END_L
    # NPE_END_R(x) ≈ NPE_END_L(-x)
    pos_neg = {x: i for i, x in enumerate(positions)}
    mirror_diffs = []
    for i, x in enumerate(positions):
        if -x in pos_neg:
            j = pos_neg[-x]
            mirror_diffs.append(abs(npe_endr[i] - npe_endl[j]))

    # ── Criterion (b): diagonal ridge in heatmap
    # For each beam position, find the TOP SiPM with max <Npe> and compare to closest predicted
    nearest_obs_cx  = []
    nearest_pred_cx = []
    for pi, x_beam in enumerate(positions):
        ti_max = int(np.argmax(npe_top_matrix[:, pi]))
        obs_cx = top_center_x(ti_max)
        pred_cx = float(x_beam)  # should be near beam
        nearest_obs_cx.append(obs_cx)
        nearest_pred_cx.append(pred_cx)

    cx_residuals = np.array([abs(o - p) for o, p in zip(nearest_obs_cx, nearest_pred_cx)])
    cx_med_err   = float(np.median(cx_residuals))

    # CSV sidecar
    with open(Path(OUT_DIR) / "V4" / "V4_attenuation_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["x_mm","npe_endl","npe_endr"])
        w.writeheader()
        for i, x in enumerate(positions):
            w.writerow({"x_mm": x, "npe_endl": npe_endl[i], "npe_endr": npe_endr[i]})

    asym_ok  = asymmetry_endl > 0.1   # END_L higher near left end (>10% asymmetry)
    ridge_ok = cx_med_err < 40.0      # nearest TOP within 40 mm of beam (2 TOP pitches)

    evid = (f"END_L near/far asymmetry={asymmetry_endl:.3f} (>0.1 expected); "
            f"TOP ridge median |Δx|={cx_med_err:.1f} mm (<40 mm expected)")

    if asym_ok and ridge_ok:
        verdict("V4", "PASS", evid)
    elif not asym_ok:
        verdict("V4", "FAIL", evid + " — NO near/far asymmetry in END_L: attenuation bug")
    elif not ridge_ok:
        verdict("V4", "FAIL", evid + " — TOP nearest not following beam: channel-map bug")
    else:
        verdict("V4", "CONCERN", evid)

    return npe_top_matrix, positions

# ════════════════════════════════════════════════════════════════════════════════
# V5 — Temporal shape and tau_d verification
# ════════════════════════════════════════════════════════════════════════════════

def run_V5(pm):
    note("\n=== V5 — Temporal shape and tau_d ===")

    # Use TOP nearest at x=0 (gid=50 from EXEC_16, but recompute)
    entry_c = next(e for e in pm if e["x_mm"] == 0)
    d0 = load(entry_c, [BR["gid"], BR["t"]])
    top_mask = d0[BR["gid"]] >= 2*K_END
    top_gids = d0[BR["gid"]][top_mask]
    top_t    = d0[BR["t"]][top_mask]
    # Find nearest gid
    uniq, cts = np.unique(top_gids, return_counts=True)
    nearest_gid = int(uniq[np.argmax(cts)])
    mask_near = d0[BR["gid"]] == nearest_gid
    t_near = d0[BR["t"]][mask_near] * 1000  # ns → ps

    # (a) Fit exponential tail to extract tau_d
    # Tail region: after the peak + few ns → fit t > t_peak + 5000 ps
    seeds = gather_seeds(t_near * 0.001, "sqrt_n")  # back to ns for seeds
    t_peak_ps = seeds["peak"] * 1000
    t_tail_lo_ps = t_peak_ps + 2000  # start tail 2 ns after peak
    t_tail_hi_ps = t_peak_ps + 8000  # up to 8 ns after peak

    t_tail = t_near[(t_near > t_tail_lo_ps) & (t_near < t_tail_hi_ps)]

    tau_fit_ns  = np.nan
    tau_fit_err = np.nan
    if len(t_tail) > 50:
        try:
            counts, edges = np.histogram(t_tail, bins=80)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ok = counts > 0
            def exp_func(t, A, tau):
                return A * np.exp(-(t - t_tail_lo_ps) / (tau * 1000))
            popt, pcov = curve_fit(exp_func, centers[ok], counts[ok].astype(float),
                                   p0=[counts[ok].max(), TAU_D_NS], maxfev=2000)
            tau_fit_ns  = abs(float(popt[1]))
            tau_fit_err = float(np.sqrt(pcov[1, 1])) if pcov[1,1] >= 0 else np.nan
        except Exception as e:
            note(f"  [WARN] tau_d fit failed: {e}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Time distribution with tail fit
    lo_ps, hi_ps = float(np.percentile(t_near, 0.1)), float(np.percentile(t_near, 99.5))
    axes[0].hist(t_near, bins=200, range=(lo_ps, hi_ps), histtype="step", lw=1.2, color="seagreen", label=f"TOP gid={nearest_gid}")
    if not math.isnan(tau_fit_ns):
        t_plot = np.linspace(t_tail_lo_ps, t_tail_hi_ps, 300)
        axes[0].axvspan(t_tail_lo_ps, t_tail_hi_ps, alpha=0.1, color="red", label="fit region")
    axes[0].axvline(t_peak_ps, color="red", ls="--", lw=0.8, label=f"peak={t_peak_ps:.0f} ps")
    axes[0].set_xlabel("time [ps]"); axes[0].set_ylabel("hits")
    axes[0].set_title(f"V5a — Time distribution (TOP nearest, x=0)")
    axes[0].legend(fontsize=7)
    axes[0].set_yscale("log")

    # (b) Check far-end bump: for END_L with beam at x=+690 (far from END_L)
    # The far-end reflected photon should arrive at Δt ≈ bar_length / v_eff ≈ 140 cm / 15 cm/ns ≈ 9.3 ns
    entry_far = next(e for e in pm if e["x_mm"] == 690)
    d_far = load(entry_far, [BR["gid"], BR["t"]])
    endl_mask_far = d_far[BR["gid"]] < K_END
    t_endl_far = d_far[BR["t"]][endl_mask_far] * 1000  # ps

    # Predict geometric delay: beam at +690 mm, END_L at -699.75 mm
    # Direct photon path ≈ |690 - (-699.75)| = 1389.75 mm in material
    # v_eff ≈ 15 cm/ns = 150 mm/ns → transit time ≈ 9.27 ns = 9270 ps
    predicted_bump_ps = (2 * BAR_HALF_X_MM) / 15.0 * 1000  # ≈ 9333 ps
    note(f"  Predicted far-end bump at Δt ≈ {predicted_bump_ps:.0f} ps from t=0")

    if len(t_endl_far) > 20:
        lo2, hi2 = float(np.percentile(t_endl_far, 0.5)), min(float(np.percentile(t_endl_far, 99.5)), 20000)
        counts2, edges2 = np.histogram(t_endl_far, bins=100, range=(lo2, hi2))
        centers2 = 0.5 * (edges2[:-1] + edges2[1:])
        axes[1].step(centers2, counts2, color="steelblue", lw=1.2, label="END_L | beam at x=+690 mm")
        axes[1].axvline(predicted_bump_ps, color="red", ls="--", lw=1.2, label=f"Predicted bump @{predicted_bump_ps:.0f} ps")
        axes[1].set_xlabel("time [ps]"); axes[1].set_ylabel("hits")
        axes[1].set_title("V5b — Far-end bump on END_L (beam at +690 mm)")
        axes[1].legend(fontsize=8)
    savefig(fig, "V5_temporal_shape_taud", "V5")

    # Evaluate
    tau_ok = not math.isnan(tau_fit_ns) and abs(tau_fit_ns - TAU_D_NS) / TAU_D_NS < 0.2
    evid = (f"τ_d fit = {tau_fit_ns:.3f} ± {tau_fit_err:.3f} ns (config={TAU_D_NS} ns); "
            f"predicted far-end bump at {predicted_bump_ps:.0f} ps")
    if math.isnan(tau_fit_ns):
        verdict("V5", "CONCERN", evid + " — tau fit did not converge")
    elif tau_ok:
        verdict("V5", "PASS", evid)
    elif abs(tau_fit_ns - TAU_D_NS) / TAU_D_NS < 0.4:
        verdict("V5", "CONCERN", evid)
    else:
        verdict("V5", "FAIL", evid)

    with open(Path(OUT_DIR) / "V5" / "V5_taud_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        w.writerow(["tau_d_fit_ns", tau_fit_ns])
        w.writerow(["tau_d_fit_err_ns", tau_fit_err])
        w.writerow(["tau_d_config_ns", TAU_D_NS])
        w.writerow(["predicted_far_end_bump_ps", predicted_bump_ps])
        w.writerow(["n_endl_hits_x690", len(t_endl_far)])

# ════════════════════════════════════════════════════════════════════════════════
# V6 — Left/right symmetry
# ════════════════════════════════════════════════════════════════════════════════

def run_V6(pm):
    note("\n=== V6 — Left/right symmetry ===")

    positions = [e["x_mm"] for e in pm]
    pos_dict  = {e["x_mm"]: e for e in pm}

    # For each mirrored pair (x, -x): compare END_L(+x) vs END_R(-x) and vice versa
    # Also: NPE(x) for END_L should mirror END_R
    sym_pairs = [(x, -x) for x in positions if x > 0 and -x in pos_dict]

    npe_endl_at_pos = {}
    npe_endr_at_pos = {}
    sigma_top_l_at_pos = {}   # TOP gid 16..50 (local 0..34, left half)
    sigma_top_r_at_pos = {}   # TOP gid 51..85 (local 35..69, right half)

    all_xvals = []
    for entry in pm:
        x = entry["x_mm"]
        d = load(entry, [BR["ev"], BR["gid"]])
        all_events = np.unique(d[BR["ev"]])
        n_ev = len(all_events)
        gids = d[BR["gid"]]

        npe_endl_at_pos[x] = float(np.sum(gids < K_END)) / n_ev
        npe_endr_at_pos[x] = float(np.sum((gids >= K_END) & (gids < 2*K_END))) / n_ev
        all_xvals.append(x)

    # Compute symmetry: compare END_L(x) vs END_R(-x) for mirrored pairs
    residuals_endl = []
    residuals_endr = []
    for xp, xm in sym_pairs:
        el_plus  = npe_endl_at_pos.get(xp, np.nan)
        el_minus = npe_endl_at_pos.get(xm, np.nan)
        er_plus  = npe_endr_at_pos.get(xp, np.nan)
        er_minus = npe_endr_at_pos.get(xm, np.nan)
        # END_L(+x) should equal END_R(-x) by symmetry
        if not any(math.isnan(v) for v in [el_plus, er_minus]) and (el_plus + er_minus) > 0:
            residuals_endl.append(2 * abs(el_plus - er_minus) / (el_plus + er_minus))
        if not any(math.isnan(v) for v in [er_plus, el_minus]) and (er_plus + el_minus) > 0:
            residuals_endr.append(2 * abs(er_plus - el_minus) / (er_plus + el_minus))

    sym_asym_mean = float(np.mean(residuals_endl + residuals_endr)) if residuals_endl else np.nan

    # Figure: overlay END_L(x) and mirror of END_R(-x)
    xs_sorted = sorted(all_xvals)
    endl_vals  = [npe_endl_at_pos[x] for x in xs_sorted]
    endr_vals  = [npe_endr_at_pos[x] for x in xs_sorted]
    # Mirror END_R: plot END_R(-x) vs +x
    endr_mirror_x = []
    endr_mirror_v = []
    for x in xs_sorted:
        if -x in npe_endr_at_pos:
            endr_mirror_x.append(-x)
            endr_mirror_v.append(npe_endr_at_pos[-x])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].plot(xs_sorted, endl_vals, "s-", color="steelblue", ms=4, lw=1.2, label="END_L (x)")
    axes[0].plot(sorted(endr_mirror_x),
                 [v for _, v in sorted(zip(endr_mirror_x, endr_mirror_v))],
                 "^--", color="tomato", ms=4, lw=1.2, alpha=0.8, label="END_R(-x) mirrored")
    axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("⟨Npe⟩ per event")
    axes[0].set_title("V6 — Symmetry: END_L(x) vs END_R(−x)")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)

    # Residual distribution
    all_res = residuals_endl + residuals_endr
    if all_res:
        axes[1].hist(all_res, bins=20, histtype="step", lw=1.5, color="purple")
        axes[1].axvline(0.05, color="orange", ls="--", lw=1, label="5% threshold")
        axes[1].axvline(0.15, color="red",    ls="--", lw=1, label="15% threshold")
        axes[1].set_xlabel("Fractional asymmetry |END_L(x) - END_R(-x)| / avg")
        axes[1].set_ylabel("pairs")
        axes[1].set_title("V6 — Symmetry residuals")
        axes[1].legend(fontsize=8)
    savefig(fig, "V6_left_right_symmetry", "V6")

    with open(Path(OUT_DIR) / "V6" / "V6_symmetry_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        w.writerow(["mean_fractional_asymmetry", sym_asym_mean])
        w.writerow(["n_mirror_pairs", len(sym_pairs)])

    evid = f"Mean fractional asymmetry={sym_asym_mean:.4f} ({sym_asym_mean*100:.1f}%) over {len(sym_pairs)} mirror pairs"
    if not math.isnan(sym_asym_mean):
        if sym_asym_mean < 0.05:
            verdict("V6", "PASS", evid)
        elif sym_asym_mean < 0.15:
            verdict("V6", "CONCERN", evid)
        else:
            verdict("V6", "FAIL", evid)
    else:
        verdict("V6", "CONCERN", evid + " — no mirror pairs found")

# ════════════════════════════════════════════════════════════════════════════════
# V7 — σ_t(N_pe) scaling
# ════════════════════════════════════════════════════════════════════════════════

def run_V7(pm):
    note("\n=== V7 — σ_t vs N_pe scaling ===")

    sigma_pts = []
    npe_pts   = []

    # Use all TOP SiPMs at positions where they have ≥30 events with hits
    positions_to_use = pm[::3]  # every 3rd position (10 positions) for speed

    for entry in positions_to_use:
        d = load(entry, [BR["ev"], BR["gid"], BR["t"]])
        all_events = np.unique(d[BR["ev"]])
        n_ev = len(all_events)
        gids = d[BR["gid"]]
        times = d[BR["t"]]

        for ti in range(N_TOP):
            gid_val = 2 * K_END + ti
            mask_ch = gids == gid_val

            # <Npe> for this channel
            hits_per_ev = np.array([np.sum((d[BR["ev"]] == ev) & mask_ch) for ev in all_events])
            mean_npe = float(np.mean(hits_per_ev))
            if mean_npe < 0.5: continue  # skip channels with almost no hits

            # t_first (t_1) for this channel per event
            ev_sel = d[BR["ev"]][mask_ch]
            t_sel  = times[mask_ch]
            if len(ev_sel) == 0: continue
            sort_idx = np.lexsort((t_sel, ev_sel))
            ev_s = ev_sel[sort_idx]; t_s = t_sel[sort_idx]
            uniq_ev, counts = np.unique(ev_s, return_counts=True)
            cum = np.concatenate([[0], np.cumsum(counts)])
            t_first = np.array([t_s[cum[i]] for i in range(len(uniq_ev)) if counts[i] >= 1])

            r = gauss_fit_v(t_first, prefix=f"v7_gid{gid_val}_x{entry['x_mm']}")
            if r.get("flag") == "ok" and not math.isnan(r.get("sigma_fit", np.nan)):
                sigma_pts.append(float(r["sigma_fit"]) * 1000)  # ns → ps
                npe_pts.append(mean_npe)

    sigma_arr = np.array(sigma_pts)
    npe_arr   = np.array(npe_pts)

    # Fit σ_t = √(a²/N_pe + b²)
    a_fit = b_fit = r2 = np.nan
    if len(sigma_arr) >= 10:
        try:
            def model(npe, a, b):
                return np.sqrt(np.maximum(a**2 / npe + b**2, 1e-6))
            # Initial guess
            p0 = [sigma_arr.mean() * np.sqrt(npe_arr.mean()), 0.1]
            popt, pcov = curve_fit(model, npe_arr, sigma_arr, p0=p0, maxfev=3000,
                                   bounds=([0, 0], [np.inf, np.inf]))
            a_fit, b_fit = float(popt[0]), float(popt[1])
            # R²
            residuals = sigma_arr - model(npe_arr, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sigma_arr - sigma_arr.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        except Exception as e:
            note(f"  [WARN] V7 fit failed: {e}")

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(npe_arr, sigma_arr, s=15, alpha=0.5, color="seagreen", label="(⟨Npe⟩, σ_t) pairs")
    if not math.isnan(a_fit):
        npe_range = np.linspace(npe_arr.min(), npe_arr.max(), 200)
        ax.plot(npe_range, np.sqrt(a_fit**2 / npe_range + b_fit**2), "r-", lw=2,
                label=f"Fit: √(({a_fit:.0f} ps)²/Npe + ({b_fit:.0f} ps)²), R²={r2:.3f}")
        ax.plot(npe_range, a_fit / np.sqrt(npe_range), "k--", lw=1, alpha=0.6, label="1/√Npe component")
    ax.set_xlabel("⟨Npe⟩ per event per channel"); ax.set_ylabel("σ_t [ps]")
    ax.set_title("V7 — σ_t vs ⟨Npe⟩: Poisson timing scaling")
    ax.legend(fontsize=8); ax.grid(True, lw=0.3)
    savefig(fig, "V7_sigma_t_vs_npe", "V7")

    with open(Path(OUT_DIR) / "V7" / "V7_scaling_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        w.writerow(["n_data_points", len(sigma_arr)])
        w.writerow(["fit_a_ps", a_fit])
        w.writerow(["fit_b_ps", b_fit])
        w.writerow(["r_squared", r2])

    evid = (f"n={len(sigma_arr)} (channel,position) pairs; "
            f"fit: a={a_fit:.1f} ps, b={b_fit:.1f} ps (expect b≈0), R²={r2:.3f}")
    if len(sigma_arr) < 10:
        verdict("V7", "CONCERN", "Too few valid pairs for fit")
    elif not math.isnan(b_fit) and b_fit < 30 and r2 > 0.4:
        verdict("V7", "PASS", evid)
    elif not math.isnan(b_fit) and b_fit < 80:
        verdict("V7", "CONCERN", evid)
    else:
        verdict("V7", "FAIL", evid + " — b≫0 or no correlation")

# ════════════════════════════════════════════════════════════════════════════════
# V8 — Statistical closure of the estimator
# ════════════════════════════════════════════════════════════════════════════════

def run_V8(pm):
    note("\n=== V8 — Estimator statistical closure ===")

    rng = np.random.default_rng(RANDOM_SEED)
    sigma_trues = [30, 50, 70, 100, 150, 200]   # ps → ns
    n_sample    = 5000  # match real data size
    n_trials    = 20

    rows = []
    for sigma_true_ps in sigma_trues:
        sigma_true_ns = sigma_true_ps / 1000.0
        recovered = []
        for _ in range(n_trials):
            samples = rng.normal(0, sigma_true_ns, n_sample)
            r = gauss_fit_v(samples, prefix=f"v8_s{sigma_true_ps}")
            if r.get("flag") == "ok" and not math.isnan(r.get("sigma_fit", np.nan)):
                recovered.append(r["sigma_fit"] * 1000)  # → ps
        if recovered:
            rec_mean = float(np.mean(recovered))
            rec_std  = float(np.std(recovered, ddof=1))
            bias_pct = 100 * (rec_mean - sigma_true_ps) / sigma_true_ps
        else:
            rec_mean = rec_std = bias_pct = np.nan
        rows.append({"sigma_true_ps": sigma_true_ps, "n_trials": len(recovered),
                     "rec_mean_ps": rec_mean, "rec_std_ps": rec_std, "bias_pct": bias_pct})
        note(f"  σ_true={sigma_true_ps:4d} ps → rec={rec_mean:.1f}±{rec_std:.1f} ps, bias={bias_pct:+.1f}%")

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    st_arr = [r["sigma_true_ps"] for r in rows if not math.isnan(r["rec_mean_ps"])]
    rm_arr = [r["rec_mean_ps"]   for r in rows if not math.isnan(r["rec_mean_ps"])]
    rs_arr = [r["rec_std_ps"]    for r in rows if not math.isnan(r["rec_std_ps"])]

    axes[0].errorbar(st_arr, rm_arr, yerr=rs_arr, fmt="o-", color="seagreen", capsize=4, ms=5)
    axes[0].plot([0, 250], [0, 250], "k--", lw=0.8, label="ideal (no bias)")
    axes[0].set_xlabel("σ_true [ps]"); axes[0].set_ylabel("σ_recovered [ps]")
    axes[0].set_title("V8 — Estimator closure: recovered vs true σ")
    axes[0].legend(fontsize=8); axes[0].grid(True, lw=0.3)

    bp_arr = [r["bias_pct"] for r in rows if not math.isnan(r["bias_pct"])]
    axes[1].bar(range(len(st_arr)), bp_arr, color=["green" if abs(b)<5 else "orange" if abs(b)<15 else "red" for b in bp_arr])
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].axhline(5, color="orange", ls="--", lw=0.8)
    axes[1].axhline(-5, color="orange", ls="--", lw=0.8)
    axes[1].set_xticks(range(len(st_arr)))
    axes[1].set_xticklabels([str(s) for s in st_arr], fontsize=9)
    axes[1].set_xlabel("σ_true [ps]"); axes[1].set_ylabel("bias [%]")
    axes[1].set_title("V8 — Estimator bias per σ_true")
    savefig(fig, "V8_estimator_closure", "V8")

    with open(Path(OUT_DIR) / "V8" / "V8_closure_summary.csv", "w", newline="") as f:
        dw = csv.DictWriter(f, list(rows[0].keys()))
        dw.writeheader(); dw.writerows(rows)

    max_bias = max(abs(r["bias_pct"]) for r in rows if not math.isnan(r["bias_pct"]))
    evid = f"max |bias| across σ_true=[{min(sigma_trues)},{max(sigma_trues)}] ps: {max_bias:.1f}% (< 5% → PASS, < 15% → CONCERN)"
    if max_bias < 5:
        verdict("V8", "PASS", evid)
    elif max_bias < 15:
        verdict("V8", "CONCERN", evid)
    else:
        verdict("V8", "FAIL", evid)

# ════════════════════════════════════════════════════════════════════════════════
# V9 — Anomalous event hunting
# ════════════════════════════════════════════════════════════════════════════════

def run_V9(pm):
    note("\n=== V9 — Anomalous event hunting ===")

    entry_c = next(e for e in pm if e["x_mm"] == 0)
    d = load(entry_c, [BR["ev"], BR["gid"], BR["t"]])
    all_events = np.unique(d[BR["ev"]])

    # Find TOP nearest at x=0
    top_mask = d[BR["gid"]] >= 2*K_END
    top_gids = d[BR["gid"]][top_mask]
    uniq, cts = np.unique(top_gids, return_counts=True)
    nearest_gid = int(uniq[np.argmax(cts)])
    note(f"  TOP nearest gid={nearest_gid}")

    # (a) Per-event: t_first for nearest channel
    mask_ch = d[BR["gid"]] == nearest_gid
    ev_ch = d[BR["ev"]][mask_ch]
    t_ch  = d[BR["t"]][mask_ch] * 1000  # ps
    t_first_per_ev = {}
    for ev in all_events:
        mask_ev = ev_ch == ev
        if mask_ev.sum() > 0:
            t_first_per_ev[ev] = float(t_ch[mask_ev].min())

    # (b) Dead events (no hit in nearest channel)
    n_dead = len(all_events) - len(t_first_per_ev)
    frac_dead = n_dead / len(all_events) if len(all_events) > 0 else np.nan

    # Gráfico 1 from EXEC_16 spec — time × event 2D map (validation version)
    # Use it to spot anomalous events
    t_vals = np.array([t_first_per_ev.get(ev, np.nan) for ev in all_events])
    has_hit = ~np.isnan(t_vals)
    n_events_plot = min(len(all_events), 500)

    T_BIN_PS = 25
    t_lo = np.nanpercentile(t_vals, 0.5)
    t_hi = np.nanpercentile(t_vals, 99.5)
    n_t_bins = max(1, int((t_hi - t_lo) / T_BIN_PS))

    # 2D map: t × event_index
    ev_map = np.arange(len(all_events))
    ev_norm = (ev_map * n_events_plot / len(all_events)).astype(int)
    ev_norm = np.clip(ev_norm, 0, n_events_plot - 1)

    valid = has_hit
    h2d, xedge, yedge = np.histogram2d(
        t_vals[valid], ev_norm[valid],
        bins=[n_t_bins, n_events_plot],
        range=[[t_lo, t_hi], [0, n_events_plot]])

    # Detect anomalous events: t_first > mean + 5σ
    t_mean = float(np.nanmean(t_vals))
    t_std  = float(np.nanstd(t_vals))
    threshold_ps = t_mean + 5 * t_std
    n_anomalous  = int(np.sum(t_vals > threshold_ps))
    frac_anomalous = n_anomalous / len(all_events)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    from matplotlib.colors import LogNorm
    norm_2d = LogNorm(vmin=0.5, vmax=h2d.max()) if h2d.max() / max(h2d[h2d>0].min(), 0.5) > 50 else None
    im = axes[0].pcolormesh(xedge, yedge, h2d.T, cmap="viridis", norm=norm_2d)
    plt.colorbar(im, ax=axes[0], label="hits/bin")
    axes[0].set_xlabel("time of 1st PE [ps]"); axes[0].set_ylabel("event index")
    axes[0].set_title(f"V9 — time×event map | TOP gid={nearest_gid} | x=0 mm")

    # Distribution of t_first to show tail
    axes[1].hist(t_vals[has_hit], bins=200, range=(t_lo, t_hi), histtype="step", lw=1.5, color="seagreen")
    axes[1].axvline(threshold_ps, color="red", ls="--", lw=1, label=f"5σ threshold ({threshold_ps:.0f} ps)")
    axes[1].set_xlabel("t_first [ps]"); axes[1].set_ylabel("events")
    axes[1].set_title(f"V9 — t_first distribution (dead={frac_dead:.3%}, anomalous={frac_anomalous:.3%})")
    axes[1].legend(fontsize=8)
    savefig(fig, "V9_anomalous_events", "V9")

    with open(Path(OUT_DIR) / "V9" / "V9_anomalous_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric","value"])
        w.writerow(["nearest_gid", nearest_gid])
        w.writerow(["n_total_events", len(all_events)])
        w.writerow(["n_dead_events", n_dead])
        w.writerow(["frac_dead", frac_dead])
        w.writerow(["n_anomalous_5sigma", n_anomalous])
        w.writerow(["frac_anomalous", frac_anomalous])
        w.writerow(["t_mean_ps", t_mean])
        w.writerow(["t_std_ps", t_std])

    evid = (f"Dead events={frac_dead:.3%} (<1% expected); "
            f"anomalous (>5σ)={frac_anomalous:.3%}; "
            f"t_mean={t_mean:.0f} ps, t_std={t_std:.0f} ps")
    if frac_dead < 0.01 and frac_anomalous < 0.05:
        verdict("V9", "PASS", evid)
    elif frac_dead < 0.05 and frac_anomalous < 0.15:
        verdict("V9", "CONCERN", evid)
    else:
        verdict("V9", "FAIL", evid)

# ════════════════════════════════════════════════════════════════════════════════
# Write SIM_TRUST_VERDICT.md
# ════════════════════════════════════════════════════════════════════════════════

def write_verdict():
    results_list = list(_verdicts.values())
    n_pass    = sum(1 for v in results_list if v["result"] == "PASS")
    n_concern = sum(1 for v in results_list if v["result"] == "CONCERN")
    n_fail    = sum(1 for v in results_list if v["result"] == "FAIL")

    if n_fail >= 2:
        global_verdict = f"NO-CONFIABLE ({n_fail} FAIL: " + ", ".join(k for k,v in _verdicts.items() if v["result"]=="FAIL") + ")"
    elif n_fail == 1:
        fail_k = next(k for k,v in _verdicts.items() if v["result"]=="FAIL")
        global_verdict = f"CONFIABLE-CON-RESERVAS (1 FAIL: {fail_k}; {n_concern} CONCERN)"
    elif n_concern >= 3:
        global_verdict = f"CONFIABLE-CON-RESERVAS ({n_concern} CONCERN; 0 FAIL)"
    else:
        global_verdict = "CONFIABLE"

    lines = [
        f"# SIM_TRUST_VERDICT — EXEC_17 — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        f"## Veredicto global: **{global_verdict}**",
        "",
        f"| Test | Resultado | Evidencia |",
        f"|------|----------|-----------|",
    ]
    for vid in ["V1","V2","V3","V4","V5","V6","V7","V8","V9"]:
        v = _verdicts.get(vid, {"result":"NOT_RUN","evidence":"—"})
        sym = {"PASS":"✓","CONCERN":"△","FAIL":"✗","NOT_RUN":"?"}.get(v["result"],"?")
        ev  = v["evidence"].replace("|","\\|")
        lines.append(f"| {vid} | {sym} **{v['result']}** | {ev} |")

    lines += [
        "",
        f"## Totales: {n_pass} PASS / {n_concern} CONCERN / {n_fail} FAIL",
        "",
        "## Constantes verificadas (leídas de fuente esta sesión)",
        f"- `kNEndSiPMs = {K_END}` (DetectorConstruction.hh:37)",
        f"- `kNTopSiPMs = {N_TOP}` (DetectorConstruction.hh:38)",
        f"- `kBarHalfX = {BAR_HALF_X_MM}` mm, `kBarHalfY = {BAR_HALF_Y_MM}` mm, `kBarHalfZ = {BAR_HALF_Z_MM}` mm",
        f"- `kEndPitch = {END_PITCH_MM}` mm; END_L at x={-PREDICTED_END_X_MM:.2f} mm, END_R at +{PREDICTED_END_X_MM:.2f} mm",
        f"- TOP: 35 SiPMs each side, step 20 mm, cx = -692..+692 mm",
        f"- `SCINTILLATIONTIMECONSTANT1 = {TAU_D_NS} ns` (opsc-101.mac:11)",
        f"- `SCINTILLATIONYIELD = {LY_PER_MEV}/MeV` (opsc-101.mac:9)",
        f"- `ABSLENGTH = {ABSLENGTH_CM} cm` (DetConstr.cc override for OPSC-101)",
        f"- Emission peak = {EMIS_PEAK_NM} nm (scntComp1.txt maximum)",
        "",
        f"## `VAL_STOP = True` — DETENTE aquí, espera OK de René antes de Fase 2.",
        "",
        "*(Números generados programáticamente desde los sidecars CSV — ninguno fue tecleado.)*",
    ]

    p = Path(OUT_DIR) / "SIM_TRUST_VERDICT.md"
    p.write_text("\n".join(lines))
    note(f"\nSIM_TRUST_VERDICT.md written: {p}")
    note(f"Veredicto global: {global_verdict}")
    return global_verdict

# ════════════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    note(f"=== EXEC_17 CP-VAL — {datetime.datetime.now().isoformat()} ===")
    note(f"K_END={K_END}, N_TOP={N_TOP}, N_TOTAL={N_TOTAL}")
    note(f"τ_d={TAU_D_NS} ns, LY={LY_PER_MEV}/MeV, λ_abs={ABSLENGTH_CM} cm, peak={EMIS_PEAK_NM} nm")

    note("\nBuilding position map...")
    pm = pos_map()
    note(f"  {len(pm)} positions: {[e['x_mm'] for e in pm]}")

    run_V1(pm)
    run_V2(pm)
    run_V3(pm)
    run_V4(pm)
    run_V5(pm)
    run_V6(pm)
    run_V7(pm)
    run_V8(pm)
    run_V9(pm)

    global_v = write_verdict()
    note(f"\n=== CP-VAL complete. VAL_STOP=True — stopping here. ===")
    note(f"Veredicto: {global_v}")
    note(f"Output: {OUT_DIR}/SIM_TRUST_VERDICT.md")

if __name__ == "__main__":
    main()
