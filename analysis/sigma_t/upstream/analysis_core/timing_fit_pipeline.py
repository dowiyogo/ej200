#!/usr/bin/env python3.12
"""
timing_fit_pipeline.py — EXEC_16 main loop.

Processes 31-position ROOT scan files for one or two materials,
runs the seeded Gaussian fit on each (group, N, position) combination,
writes .root/.csv/.meta.json sidecars, and generates publication figures.

Usage:
    python3.12 timing_fit_pipeline.py --config config/exec16_config.yaml [--materials EJ-230]

See README.md for full documentation.
"""

import argparse
import math
import re
import sys
import subprocess
from pathlib import Path

import numpy as np
import yaml
import uproot
import matplotlib
matplotlib.use("Agg")   # no display — WSL2 without X server
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# ── add the package root to sys.path so `lib.*` imports work anywhere ──────
sys.path.insert(0, str(Path(__file__).parent))

from lib.fit_engine    import fit_core_gaussian
from lib.robust_seeds  import gather_seeds
from lib.gates         import (gate_0_verify_sha, gate_1_verify_schema,
                                gate_2_physics_bounds, gate_3_check_group,
                                gate_3d_write_manifest)
from lib.sidecar       import (write_root_sidecar, write_csv_sidecar,
                                write_metadata_json)

# ─── constants ────────────────────────────────────────────────────────────────
FACE_TOP = 2    # face_type value for TOP SiPM channels
FACE_EL  = 0    # face_type value for END-left channels
FACE_ER  = 1    # face_type value for END-right channels


# ─── file name helpers ────────────────────────────────────────────────────────

def x_to_tag(x_mm: int) -> str:
    """
    Convert an x position in mm to the filename segment used in the ROOT files.
    Examples:  -690 → "xm690mm",   0 → "xp0mm",   +690 → "xp690mm".
    """
    sign = "m" if x_mm < 0 else "p"      # "m" for minus, "p" for plus
    return f"x{sign}{abs(x_mm)}mm"


def build_filepath(cfg: dict, mat_cfg: dict, x_mm: int) -> Path:
    """
    Construct the absolute path to one raw ROOT file for a given position.

    Expected pattern:
        {INPUT_ROOT_BASE}/{subdir}/{file_prefix}_x{sign}{abs}mm_300ev_{sha_tag}.root
    Example:
        .../raw/ej230_endtop/ej230_endtop_xm690mm_300ev_ca2f1c3.root
    """
    base    = Path(cfg["INPUT_ROOT_BASE"])
    subdir  = mat_cfg["subdir"]          # e.g. "ej230_endtop"
    prefix  = mat_cfg["file_prefix"]     # e.g. "ej230_endtop"
    sha_tag = mat_cfg["sha_tag"]         # e.g. "ca2f1c3"
    xtag    = x_to_tag(x_mm)
    fname   = f"{prefix}_{xtag}_300ev_{sha_tag}.root"
    return base / subdir / fname


# ─── group construction ───────────────────────────────────────────────────────

def best_gids_by_count(face_mask: np.ndarray, global_id: np.ndarray, k: int) -> np.ndarray:
    """
    Return the k global_ids with the most hits in the given face.

    'Most hits' is used as a proxy for 'nearest to the gun': the top SiPM
    directly above the muon gun position collects the most photons.
    For the END SiPMs, all 8 have roughly equal hit counts at a given x
    (they're all 8 on the same face), so SUM4 takes the 4 with the most.

    Parameters
    ----------
    face_mask : np.ndarray  bool mask selecting one face from the full hit array
    global_id : np.ndarray  full global_id array (not yet masked)
    k         : int         number of SiPMs to select

    Returns
    -------
    np.ndarray of int — the k selected global_ids, sorted descending by count
    """
    gids_in_face = global_id[face_mask]        # IDs of all hits on this face
    if len(gids_in_face) == 0:
        return np.array([], dtype=int)         # no hits at all on this face
    uniq, cts = np.unique(gids_in_face, return_counts=True)
    order = np.argsort(-cts)                   # sort descending by hit count
    return uniq[order[:k]]                     # return top k global_ids


def compute_tN(event_id: np.ndarray, time_ns: np.ndarray,
               hit_mask: np.ndarray, N: int, all_events: np.ndarray) -> np.ndarray:
    """
    For each event, find the N-th smallest hit time among selected hits.

    Uses a sort + cumsum approach to avoid a Python-level event loop:
    all sorting and indexing is done in numpy, which is ~50x faster than
    per-event slicing for 300 events × ~30 000 hits.

    Parameters
    ----------
    event_id   : np.ndarray  — event_id branch (full array)
    time_ns    : np.ndarray  — time_ns branch (full array)
    hit_mask   : np.ndarray  — boolean mask for which hits belong to this group
    N          : int         — photon multiplicity threshold
    all_events : np.ndarray  — sorted unique event_ids in this file

    Returns
    -------
    np.ndarray  — t_N values (ns) for events that had >= N hits; length <= len(all_events)
    """
    ev_sel = event_id[hit_mask]     # event ids for selected hits only
    t_sel  = time_ns[hit_mask]      # times for selected hits

    if len(ev_sel) == 0:
        return np.array([], dtype=float)    # no hits at all for this group at this position

    # sort all selected hits by (event_id, time_ns) so we can find the N-th per event
    sort_idx  = np.lexsort((t_sel, ev_sel))     # primary key: event_id; secondary: time_ns
    ev_sorted = ev_sel[sort_idx]
    t_sorted  = t_sel[sort_idx]

    # get sorted unique events and how many hits each has
    uniq_ev, counts = np.unique(ev_sorted, return_counts=True)

    # cumulative sum of counts → starting index in t_sorted for each event
    cum = np.concatenate([[0], np.cumsum(counts)])

    # collect t_N for events that have >= N hits
    t_N_list = []
    for i, cnt in enumerate(counts):
        if cnt >= N:
            # the N-th hit in this event is at position cum[i] + N-1 in the sorted array
            t_N_list.append(t_sorted[cum[i] + N - 1])

    return np.array(t_N_list, dtype=float)


def extract_group_tN(arr: dict, group: str, N: int) -> tuple:
    """
    Build the t_N array for a given group and N from the raw hit arrays.

    Parameters
    ----------
    arr   : dict  — {branch_name: np.ndarray} loaded from uproot
    group : str   — one of the nine group names
    N     : int   — photon multiplicity

    Returns
    -------
    t_N       : np.ndarray   — t_N values per event (ns)
    sel_gids  : list[int]    — global_ids that were selected for this group/position
    n_events  : int          — total events in this file (denominator for efficiency)
    """
    event_id  = arr["event_id"]
    face_type = arr["face_type"]
    global_id = arr["global_id"]
    time_ns   = arr["time_ns"]

    all_events = np.unique(event_id)     # sorted unique events; used as denominator
    n_events   = len(all_events)

    # ── group → (face_mask, selected_gids) mapping ──────────────────────────
    # Each group is defined by which face(s) and which global_ids to include.

    if group == "TOP_NEAREST":
        # single highest-hit-count TOP SiPM for this position
        mask = (face_type == FACE_TOP)
        sel_gids = best_gids_by_count(mask, global_id, 1)
        if len(sel_gids) == 0:
            return np.array([]), [], n_events
        hit_mask = mask & np.isin(global_id, sel_gids)

    elif group == "TOP_SUM4":
        # 4 highest-hit-count TOP SiPMs merged
        mask = (face_type == FACE_TOP)
        sel_gids = best_gids_by_count(mask, global_id, 4)
        if len(sel_gids) == 0:
            return np.array([]), [], n_events
        hit_mask = mask & np.isin(global_id, sel_gids)

    elif group == "TOP_SUM8":
        # 8 highest-hit-count TOP SiPMs merged
        mask = (face_type == FACE_TOP)
        sel_gids = best_gids_by_count(mask, global_id, 8)
        if len(sel_gids) == 0:
            return np.array([]), [], n_events
        hit_mask = mask & np.isin(global_id, sel_gids)

    elif group == "END_L_SUM4":
        # 4 best END-left SiPMs (face_type=0, global_id 0–7)
        mask = (face_type == FACE_EL)
        sel_gids = best_gids_by_count(mask, global_id, 4)
        if len(sel_gids) == 0:
            return np.array([]), [], n_events
        hit_mask = mask & np.isin(global_id, sel_gids)

    elif group == "END_R_SUM4":
        # 4 best END-right SiPMs (face_type=1, global_id 8–15)
        mask = (face_type == FACE_ER)
        sel_gids = best_gids_by_count(mask, global_id, 4)
        if len(sel_gids) == 0:
            return np.array([]), [], n_events
        hit_mask = mask & np.isin(global_id, sel_gids)

    elif group == "END_L_SUM8":
        # all 8 END-left SiPMs (no selection needed: just use face_type=0)
        hit_mask = (face_type == FACE_EL)
        sel_gids = list(np.unique(global_id[hit_mask]))   # informational: all IDs present

    elif group == "END_R_SUM8":
        # all 8 END-right SiPMs (face_type=1)
        hit_mask = (face_type == FACE_ER)
        sel_gids = list(np.unique(global_id[hit_mask]))

    elif group == "ENDS_SUM4PAIR":
        # best 4 from END-left + best 4 from END-right merged into one stream
        mask_L  = (face_type == FACE_EL)
        mask_R  = (face_type == FACE_ER)
        gids_L  = best_gids_by_count(mask_L, global_id, 4)
        gids_R  = best_gids_by_count(mask_R, global_id, 4)
        # combine: both face masks AND selected global_ids
        hit_mask = (
            (mask_L & np.isin(global_id, gids_L)) |
            (mask_R & np.isin(global_id, gids_R))
        )
        sel_gids = list(gids_L) + list(gids_R)

    elif group == "ENDS_SUM8PAIR":
        # all 8 END-left + all 8 END-right merged
        hit_mask = (face_type == FACE_EL) | (face_type == FACE_ER)
        sel_gids = list(np.unique(global_id[hit_mask]))

    else:
        raise ValueError(f"extract_group_tN: unknown group '{group}'")

    t_N = compute_tN(event_id, time_ns, hit_mask, N, all_events)
    return t_N, list(map(int, sel_gids)), n_events


# ─── figure generation ────────────────────────────────────────────────────────

def slugify(group: str) -> str:
    """Convert group name to lowercase file-safe slug. 'TOP_SUM4' → 'top_sum4'."""
    return group.lower()


def _info_box(ax, mat_label, branch, sha, n_ev_pos, n_pos, N):
    """Add a small metadata box to the figure corner."""
    ax.text(
        0.99, 0.01,
        f"{mat_label} | {branch} @ {sha}\n"
        f"{n_ev_pos} ev/pos | {n_pos} pos | N={N} | EXEC_16",
        transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
        family="monospace", bbox=dict(boxstyle="round,pad=0.3", alpha=0.1),
    )


def plot_resolution(results_by_group: dict, group: str, mat_slug: str, mat_label: str,
                    branch: str, sha: str, n_pos: int, cfg: dict, output_dir: str) -> list:
    """
    Figure 1: σ_fit [ps] vs x_gun [mm], two panels (N=4 and N=20).

    Low-efficiency points (eff < EFFICIENCY_FLOOR) are shown as hollow
    markers with dashed error bars to indicate unreliable estimates.
    """
    n_ev_pos = 300   # events per position (hard-coded for this scan)
    eff_floor = cfg.get("EFFICIENCY_FLOOR", 0.05)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"{mat_label} — {group}  σ_fit(x)  [EndTop readout]",
        fontsize=12, fontweight="bold"
    )

    N_list   = sorted(results_by_group[group].keys())   # [4, 20]
    colors   = {4: "#2166AC", 20: "#D6604D"}            # blue for N=4, red for N=20

    for ax, N in zip(axes, N_list):
        rows = results_by_group[group][N]
        x_all  = np.array([r["x_mm"]                          for r in rows])
        sf_all = np.array([r.get("sigma_fit",     float("nan")) for r in rows]) * 1000  # ps
        # combined error: max(TF1 error, bootstrap)
        err_all = np.array([
            max(r.get("sigma_fit_err", 0)*1000, r.get("bootstrap_err", 0)*1000)
            for r in rows
        ])
        eff_all = np.array([r.get("efficiency", 0.0) for r in rows])

        # split into high-efficiency and low-efficiency points
        hi  = eff_all >= eff_floor
        lo  = (eff_all < eff_floor) & (~np.isnan(sf_all))

        col = colors.get(N, "gray")

        if hi.any():
            ax.errorbar(
                x_all[hi], sf_all[hi], yerr=err_all[hi],
                fmt="o-", color=col, ms=5, lw=1.5, capsize=3, elinewidth=1,
                label=f"N={N}  (eff ≥ {eff_floor:.0%})"
            )
        if lo.any():
            # hollow markers for low-efficiency points
            ax.errorbar(
                x_all[lo], sf_all[lo], yerr=err_all[lo],
                fmt="o--", color=col, ms=5, lw=1, capsize=3, elinewidth=0.5,
                mfc="none",   # hollow marker
                label=f"N={N}  (eff < {eff_floor:.0%}, low stat)",
                alpha=0.5,
            )

        ax.set_title(f"N = {N}", fontsize=11)
        ax.set_xlabel("x_gun [mm]", fontsize=10)
        ax.set_ylabel("σ_fit [ps]", fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        _info_box(ax, mat_label, branch, sha, n_ev_pos, n_pos, N)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    slug  = slugify(group)
    saved = []
    for ext in ("png", "pdf"):
        fpath = str(Path(output_dir) / f"{mat_slug}_{slug}_resolution_N4_N20.{ext}")
        fig.savefig(fpath, dpi=150, bbox_inches="tight")
        saved.append(fpath)
    plt.close(fig)
    return saved


def plot_chi2ndf(results_by_group: dict, group: str, mat_slug: str, mat_label: str,
                 branch: str, sha: str, n_pos: int, cfg: dict, output_dir: str) -> list:
    """
    Figure 2: χ²/NDF vs x_gun [mm], two panels (N=4 and N=20).

    A horizontal dashed line marks CHI2_NDF_WARN (= 3.0).
    Positions above the threshold are highlighted with a red background.
    """
    chi2_warn = cfg.get("CHI2_NDF_WARN", 3.0)
    n_ev_pos  = 300

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"{mat_label} — {group}  χ²/NDF(x)  [core Gaussian fit quality]",
        fontsize=12, fontweight="bold"
    )

    N_list  = sorted(results_by_group[group].keys())
    colors  = {4: "#2166AC", 20: "#D6604D"}

    for ax, N in zip(axes, N_list):
        rows    = results_by_group[group][N]
        x_all   = np.array([r["x_mm"]                              for r in rows])
        chi_all = np.array([r.get("chi2_ndf", float("nan"))        for r in rows])
        eff_all = np.array([r.get("efficiency", 0.0)               for r in rows])

        valid = ~np.isnan(chi_all)
        col   = colors.get(N, "gray")

        if valid.any():
            ax.step(x_all[valid], chi_all[valid], where="mid", color=col,
                    lw=1.5, label=f"χ²/NDF (N={N})")

            # highlight positions above the threshold with a red background
            above = valid & (chi_all > chi2_warn)
            for xi in x_all[above]:
                ax.axvspan(xi - 20, xi + 20, color="red", alpha=0.12, zorder=0)

        # reference line at χ²/NDF = 1 (ideal Gaussian core)
        ax.axhline(1.0, color="green", ls="-",  lw=1, alpha=0.5, label="χ²/NDF = 1 (ideal)")
        # warning threshold line
        ax.axhline(chi2_warn, color="red", ls="--", lw=1.2,
                   label=f"warn threshold = {chi2_warn:.1f}")

        ax.set_title(f"N = {N}", fontsize=11)
        ax.set_xlabel("x_gun [mm]", fontsize=10)
        ax.set_ylabel("χ²/NDF", fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        _info_box(ax, mat_label, branch, sha, n_ev_pos, n_pos, N)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    slug  = slugify(group)
    saved = []
    for ext in ("png", "pdf"):
        fpath = str(Path(output_dir) / f"{mat_slug}_{slug}_chi2ndf_N4_N20.{ext}")
        fig.savefig(fpath, dpi=150, bbox_inches="tight")
        saved.append(fpath)
    plt.close(fig)
    return saved


def plot_overlay(results_by_group: dict, group: str, mat_slug: str, mat_label: str,
                 branch: str, sha: str, N: int, repr_positions: list,
                 cfg: dict, output_dir: str) -> list:
    """
    Figure 3: histogram + Gaussian fit overlay for 3 representative positions.

    One subplot per position; the fit window is shown as a shaded band.
    Vertical lines show ±σ_fit (red) and ±σ_MAD (orange, dashed), so the
    reader can see how much the core fit differs from the robust seed.

    Parameters
    ----------
    N              : int  — which N to show
    repr_positions : list — e.g. [-690, 0, 690]
    """
    n_ev_pos = 300
    fit_k    = cfg.get("FIT_WINDOW_SIGMAS", 2.0)   # number of σ_MAD for the window

    rows_by_x = {r["x_mm"]: r for r in results_by_group[group][N]}

    # filter repr_positions to those actually present in the data
    rep_ok = [x for x in repr_positions if x in rows_by_x]
    if not rep_ok:
        return []

    fig, axes = plt.subplots(1, len(rep_ok), figsize=(6 * len(rep_ok), 5))
    if len(rep_ok) == 1:
        axes = [axes]   # ensure iterable even for a single panel

    fig.suptitle(
        f"{mat_label} — {group}  overlay histograms + Gaussian core fit  [N={N}]",
        fontsize=12, fontweight="bold"
    )

    for ax, x_mm in zip(axes, rep_ok):
        row = rows_by_x[x_mm]
        t_N = row.get("t_N", None)   # stored numpy array, needed for histogram

        if t_N is None or len(t_N) < 4:
            ax.text(0.5, 0.5, f"x={x_mm} mm\nn={len(t_N) if t_N is not None else 0} events",
                    ha="center", va="center", transform=ax.transAxes, color="gray")
            ax.set_title(f"x = {x_mm} mm  |  N={N}", fontsize=9)
            continue

        # ── histogram ──────────────────────────────────────────────────────
        n_bins = max(10, int(np.ceil(np.sqrt(len(t_N)))))
        t_ps   = t_N * 1000   # ns → ps for display

        counts, edges, _ = ax.hist(
            t_ps, bins=n_bins,
            color="#4878CF", alpha=0.6, edgecolor="k", lw=0.3,
            label=f"data  n={len(t_N)}",
            zorder=2,
        )

        sig_fit = row.get("sigma_fit", float("nan"))
        mu_fit  = row.get("mu_fit",    float("nan"))
        sig_mad = row.get("sigma_mad", float("nan"))
        sig_std = row.get("sigma_std", float("nan"))
        chi2    = row.get("chi2_ndf",  float("nan"))

        # ── fit window shading ─────────────────────────────────────────────
        if not math.isnan(mu_fit) and not math.isnan(sig_mad):
            win_lo = (mu_fit - fit_k * sig_mad) * 1000
            win_hi = (mu_fit + fit_k * sig_mad) * 1000
            ax.axvspan(win_lo, win_hi, alpha=0.12, color="red",
                       label=f"fit window [μ ± {fit_k:.0f}·σ_MAD]", zorder=1)

        # ── Gaussian fit curve ─────────────────────────────────────────────
        if not math.isnan(sig_fit) and not math.isnan(mu_fit):
            t_plot = np.linspace(min(t_ps), max(t_ps), 300)
            A_fit  = row.get("f_root").GetParameter(0) if row.get("f_root") else np.nan
            if not math.isnan(A_fit):
                y_fit = A_fit * np.exp(-0.5 * ((t_plot / 1000 - mu_fit) / sig_fit) ** 2)
                se    = row.get("sigma_fit_err", 0)
                be    = row.get("bootstrap_err", 0)
                ax.plot(
                    t_plot, y_fit, color="red", lw=2.0, zorder=4,
                    label=(
                        f"Gauss fit:\n"
                        f"  μ = {mu_fit*1000:.1f} ± {row.get('mu_fit_err',0)*1000:.1f} ps\n"
                        f"  σ = {sig_fit*1000:.1f} ± {max(se,be)*1000:.1f} ps\n"
                        f"  χ²/NDF = {chi2:.2f}"
                    ),
                )

            # vertical lines: ±σ_fit (red dashed), ±σ_MAD (orange), ±σ_std (green)
            ax.axvline((mu_fit + sig_fit) * 1000, color="red",    ls="--", lw=1)
            ax.axvline((mu_fit - sig_fit) * 1000, color="red",    ls="--", lw=1)
            ax.axvline((mu_fit + sig_mad) * 1000, color="orange", ls="-.", lw=1,
                       label=f"σ_MAD = {sig_mad*1000:.1f} ps")
            ax.axvline((mu_fit - sig_mad) * 1000, color="orange", ls="-.", lw=1)
            ax.axvline((mu_fit + sig_std) * 1000, color="green",  ls=":",  lw=1,
                       label=f"σ_std = {sig_std*1000:.1f} ps")
            ax.axvline((mu_fit - sig_std) * 1000, color="green",  ls=":",  lw=1)
            ax.axvline(mu_fit * 1000, color="red", ls="-", lw=0.8, alpha=0.5)

        ax.set_title(f"x = {x_mm} mm  |  N={N}  |  eff={row.get('efficiency',0):.2f}",
                     fontsize=9)
        ax.set_xlabel("t_N [ps]", fontsize=9)
        ax.set_ylabel("events / bin", fontsize=9)
        ax.legend(loc="upper left", fontsize=7, framealpha=0.85)
        ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    slug  = slugify(group)
    n_tag = f"N{N}"
    saved = []
    for ext in ("png", "pdf"):
        fpath = str(Path(output_dir) / f"{mat_slug}_{slug}_overlay_m690_0_p690_{n_tag}.{ext}")
        fig.savefig(fpath, dpi=150, bbox_inches="tight")
        saved.append(fpath)
    plt.close(fig)
    return saved


# ─── main pipeline loop ───────────────────────────────────────────────────────

def run_material(material_key: str, cfg: dict, output_dir: str) -> dict:
    """
    Run the full pipeline for ONE material (all 31 positions, all groups, all N).

    Parameters
    ----------
    material_key : str  — e.g. "EJ-230"
    cfg          : dict — loaded exec16_config.yaml
    output_dir   : str

    Returns
    -------
    dict — {group: {N: list of result dicts}}
    """
    mat_cfg   = cfg["MATERIALS"][material_key]
    mat_slug  = mat_cfg["subdir"]            # e.g. "ej230_endtop"
    mat_label = material_key                 # e.g. "EJ-230"
    branch    = mat_cfg["branch"]
    repo_path = mat_cfg["repo_path"]
    exp_sha   = mat_cfg["expected_sha"]
    groups    = cfg["GROUPS"]
    N_values  = cfg["N_VALUES"]
    positions = cfg["SCAN_POSITIONS_MM"]
    req_br    = cfg["REQUIRED_BRANCHES"]
    repr_pos  = cfg.get("REPR_POSITIONS_MM", [-690, 0, 690])

    print(f"\n{'='*64}")
    print(f"  {mat_label}  |  {branch}")
    print(f"  {len(positions)} positions × {len(groups)} groups × {len(N_values)} N values = "
          f"{len(positions)*len(groups)*len(N_values)} fits")
    print(f"{'='*64}\n")

    # QA-0: verify git SHA before touching any data
    gate_0_verify_sha(repo_path, branch, exp_sha, mat_label)

    # Collect verified runtime SHA for metadata
    proc = subprocess.run(
        ["git", "-C", repo_path, "rev-parse", "--short", f"refs/heads/{branch}"],
        capture_output=True, text=True, check=True,
    )
    runtime_sha = proc.stdout.strip()

    # Collect all input file paths and verify they exist
    input_files = []
    for x_mm in positions:
        fp = build_filepath(cfg, mat_cfg, x_mm)
        if not fp.exists():
            raise SystemExit(f"[ABORT] Missing input file: {fp}")
        input_files.append(fp)

    # ── main data structure: results[group][N] = list of dicts ──────────────
    results = {g: {N: [] for N in N_values} for g in groups}

    # ── outer loop: positions ────────────────────────────────────────────────
    for pos_idx, x_mm in enumerate(positions):
        fp = build_filepath(cfg, mat_cfg, x_mm)

        # load the TTree for this position (branch-selective to control RAM)
        with uproot.open(str(fp)) as uf:
            tree = uf[cfg["TTREE_NAME"]]

            # QA-1: verify schema (only on first file; assumption: all files same schema)
            if pos_idx == 0:
                gate_1_verify_schema(set(tree.keys()), req_br, str(fp))

            # read only the 4 branches we need for group construction + timing
            arr = tree.arrays(
                ["event_id", "face_type", "global_id", "time_ns"],
                library="np",
            )

        # ── inner loops: group × N ───────────────────────────────────────────
        for group in groups:
            for N in N_values:
                # construct the t_N array for this (group, N, position)
                t_N, sel_gids, n_total = extract_group_tN(arr, group, N)
                eff = len(t_N) / n_total if n_total > 0 else 0.0

                # build a unique name prefix for ROOT objects
                x_tag    = x_to_tag(x_mm)
                grp_slug = slugify(group)
                name_pfx = f"{mat_slug}_{grp_slug}_{x_tag}_N{N}"

                # run the Gaussian fit (or return NaN if too few events)
                context = f"{mat_label}/{group}/x={x_mm}mm/N={N}"
                res = fit_core_gaussian(t_N, cfg, name_prefix=name_pfx)

                # QA-2: abort on unphysical fit result
                gate_2_physics_bounds(res, cfg, context)

                # attach ancillary fields needed by sidecar and figures
                res["x_mm"]      = x_mm
                res["group"]     = group
                res["N"]         = N
                res["sel_gids"]  = sel_gids
                res["efficiency"] = eff
                res["t_N"]       = t_N        # store raw array for overlay plots

                results[group][N].append(res)

        # progress report every 5 positions
        if (pos_idx + 1) % 5 == 0 or pos_idx == 0 or pos_idx == len(positions) - 1:
            print(f"  [{pos_idx+1:2d}/{len(positions)}] x = {x_mm:+d} mm done")

    # ── QA-3: per-group Gaussian verdict ─────────────────────────────────────
    print("\n  QA-3 checks:")
    flagged_groups = []
    for group in groups:
        for N in N_values:
            flagged, frac, n_flag = gate_3_check_group(results[group][N], f"{group}/N={N}", cfg)
            if flagged:
                flagged_groups.append(f"{group}/N={N}")

    # ── write sidecars ────────────────────────────────────────────────────────
    print("\n  Writing sidecars …")
    all_output_files = []

    geant4_version = "11.4.0"   # from earlier system checks; not dynamically queried

    for group in groups:
        grp_slug = slugify(group)

        root_path = write_root_sidecar(
            results[group], grp_slug, mat_slug, output_dir
        )
        csv_path = write_csv_sidecar(
            results[group], grp_slug, mat_slug, output_dir
        )
        json_path = write_metadata_json(
            results[group], grp_slug, mat_slug, cfg,
            input_files, runtime_sha, geant4_version, output_dir
        )
        all_output_files.extend([root_path, csv_path, json_path])

    # ── generate figures ──────────────────────────────────────────────────────
    print("\n  Generating figures …")
    fig_files = []

    for group in groups:
        grp_slug = slugify(group)
        print(f"    {group} …", end=" ", flush=True)

        # Fig 1: resolution vs position (N=4 and N=20 panels)
        fig_files.extend(
            plot_resolution(results, group, mat_slug, mat_label,
                            branch, runtime_sha, len(positions), cfg, output_dir)
        )

        # Fig 2: χ²/NDF diagnostic (N=4 and N=20 panels)
        fig_files.extend(
            plot_chi2ndf(results, group, mat_slug, mat_label,
                         branch, runtime_sha, len(positions), cfg, output_dir)
        )

        # Fig 3: overlay histograms (one figure per N)
        for N in N_values:
            fig_files.extend(
                plot_overlay(results, group, mat_slug, mat_label,
                             branch, runtime_sha, N, repr_pos, cfg, output_dir)
            )

        print("done")

    all_output_files.extend(fig_files)

    # QA-3d: SHA-256 manifest of all outputs
    manifest_path = str(Path(output_dir) / f"{mat_slug}_manifest_sha256.txt")
    gate_3d_write_manifest(all_output_files, manifest_path)
    all_output_files.append(manifest_path)

    return results, flagged_groups, all_output_files


# ─── summary writer ───────────────────────────────────────────────────────────

def write_summary(results: dict, material_key: str, mat_cfg: dict, cfg: dict,
                  flagged_groups: list, all_output_files: list,
                  output_dir: str) -> str:
    """Write cp3_summary_{mat_slug}.txt with key metrics per group."""

    mat_slug = mat_cfg["subdir"]
    N_values = cfg["N_VALUES"]
    positions= cfg["SCAN_POSITIONS_MM"]
    eff_floor= cfg.get("EFFICIENCY_FLOOR", 0.05)

    lines = [
        "═" * 65,
        f"CP3 Summary — {material_key} EndTop, {len(positions)} positions,"
        f" {len(cfg['GROUPS'])} groups, {len(N_values)} N values",
        "═" * 65,
        "",
        f"Material:  {material_key} ({mat_cfg['sslg4_code']})",
        f"Branch:    {mat_cfg['branch']} @ {mat_cfg['expected_sha']}",
        f"Readout:   EndTop",
        f"Positions: {len(positions)} ({min(positions)}…{max(positions)} mm, "
        f"Δ={positions[1]-positions[0]} mm)",
        f"Total fits:{len(positions)} × {len(cfg['GROUPS'])} × {len(N_values)} = "
        f"{len(positions)*len(cfg['GROUPS'])*len(N_values)}",
        "",
        "Fit engine: PyROOT TF1 Gaussian seeded at histogram peak with MAD×1.4826,",
        f"            window = peak ± {cfg['FIT_WINDOW_SIGMAS']} σ_MAD",
        f"Bootstrap:  {cfg['N_BOOTSTRAP']} resamples (scipy.optimize for speed)",
        f"Binning:    {cfg['BINNING_STRATEGY']} rule per position",
        "",
    ]

    for N in N_values:
        lines.append(f"SUMMARY BY GROUP — N={N}")
        lines.append("─" * 65)
        hdr = f"{'Group':<20} {'σ_fit@x=0':>10} {'σ_std@x=0':>10}  "
        hdr += f"{'Δσ%':>6}  {'χ²/NDF range':>14}  {'eff range':>12}  flagged?"
        lines.append(hdr)
        lines.append("─" * 65)

        for group in cfg["GROUPS"]:
            row_list = results[group][N]
            # find x=0 row
            row0 = next((r for r in row_list if r["x_mm"] == 0), None)

            sf0  = row0.get("sigma_fit", float("nan")) * 1000 if row0 else float("nan")
            ss0  = row0.get("sigma_std", float("nan")) * 1000 if row0 else float("nan")
            delta = (ss0 - sf0) / ss0 * 100 if (not math.isnan(ss0) and ss0 > 0) else float("nan")

            chi2_vals = [r.get("chi2_ndf", float("nan")) for r in row_list
                         if not math.isnan(r.get("chi2_ndf", float("nan")))]
            if chi2_vals:
                chi2_range = f"({min(chi2_vals):.2f}, {max(chi2_vals):.2f})"
            else:
                chi2_range = "N/A"

            eff_vals = [r.get("efficiency", 0) for r in row_list if r.get("efficiency",0) > 0]
            if eff_vals:
                eff_range = f"{min(eff_vals)*100:.0f}%–{max(eff_vals)*100:.0f}%"
            else:
                eff_range = "< 5%"

            flagged = f"{group}/N={N}" in flagged_groups
            flag_str = "⚠ YES" if flagged else "✓ No"

            sf0_str  = f"{sf0:.1f} ps"  if not math.isnan(sf0)  else "—"
            ss0_str  = f"{ss0:.1f} ps"  if not math.isnan(ss0)  else "—"
            d_str    = f"{delta:.1f}%"  if not math.isnan(delta) else "—"

            lines.append(
                f"{group:<20} {sf0_str:>10} {ss0_str:>10}  {d_str:>6}  "
                f"{chi2_range:>14}  {eff_range:>12}  {flag_str}"
            )
        lines.append("")

    lines += [
        f"Core-not-Gaussian groups: {flagged_groups if flagged_groups else 'none'}",
        "",
        "Output files:",
        f"  Sidecars (.root): {output_dir}/{mat_slug}_*_sidecars.root  (9 files)",
        f"  CSVs:             {output_dir}/{mat_slug}_*_results.csv    (9 files)",
        f"  Metadata (.json): {output_dir}/{mat_slug}_*_metadata.json  (9 files)",
        f"  Figures:          {output_dir}/{mat_slug}_*.{{png,pdf}}",
        f"  Manifest:         {output_dir}/{mat_slug}_manifest_sha256.txt",
        f"  Total outputs:    {len(all_output_files)} files",
        "",
        "READY FOR CP4 (EJ-204) — awaiting René's review.",
        "═" * 65,
    ]

    summary_text = "\n".join(lines)
    print("\n" + summary_text)

    out_path = str(Path(output_dir) / f"cp3_summary_{mat_slug}.txt")
    Path(out_path).write_text(summary_text)
    print(f"\n  Summary written: {out_path}")
    return out_path


# ─── entry point ──────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="EXEC_16 timing-fit pipeline")
    ap.add_argument("--config",    default="config/exec16_config.yaml")
    ap.add_argument("--materials", nargs="+", default=None,
                    help="Subset of materials to run; default = all in config")
    ap.add_argument("--outdir",    default=None,
                    help="Override OUTPUT_DIR from config")
    args = ap.parse_args()

    # load config — all constants from here, nothing else hard-coded
    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise SystemExit(f"Config not found: {cfg_path}")
    cfg = yaml.safe_load(cfg_path.read_text())

    output_dir = args.outdir or cfg["OUTPUT_DIR"]
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # determine which materials to run
    all_mat = list(cfg["MATERIALS"].keys())
    materials = args.materials if args.materials else all_mat

    for mat_key in materials:
        if mat_key not in cfg["MATERIALS"]:
            raise SystemExit(f"Material '{mat_key}' not in config MATERIALS section.")

        results, flagged, all_files = run_material(mat_key, cfg, output_dir)

        write_summary(
            results, mat_key, cfg["MATERIALS"][mat_key],
            cfg, flagged, all_files, output_dir
        )

    print(f"\n{'='*64}")
    print("  Pipeline complete.")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*64}\n")


if __name__ == "__main__":
    main()
