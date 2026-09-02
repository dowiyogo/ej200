#!/usr/bin/env python3
"""
top_npe_diag.py — Phase 3 diagnostic: N_pe^TOP and N_pe^END vs muon position.

Reads:
  N_pe^TOP: results/scan_end_vik_sparse_top_v2/photon_hits_ej230_ntop20_x{x}.root
             face_type=2, raw hit count pooled over all 20 TOP SiPMs
  N_pe^END: results/scan_end_vikuiti/photon_hits_ej230_x{x}.root
             face_type 0+1, raw hit count / 2 (per-face mean, consistent with
             NapNpesimejC=701.3 in napkin_macros.tex)

Output: analysis/top_npe_diag/top_npe_diag.csv
        analysis/top_npe_diag/top_npe_diag_meta.json

Gate: N_pe^END(x=0) must reproduce 701.3 ± 1.
"""

import json
import csv
from pathlib import Path
import numpy as np
import uproot

REPO     = Path("/home/rrios/ej200")
INDIR_SP = REPO / "results/scan_end_vik_sparse_top_v2"
INDIR_EO = REPO / "results/scan_end_vikuiti"
OUTDIR   = REPO / "analysis/top_npe_diag"
OUTDIR.mkdir(parents=True, exist_ok=True)

TOP_XPOS = np.array([
    -664.3, -594.3, -524.3, -454.3, -384.4, -314.4, -244.4, -174.4,
    -104.5,  -34.5,  34.5,  104.5,  174.4,  244.4,  314.4,  384.4,
     454.3,  524.3,  594.3,  664.3
])

x_scan = list(range(-600, 601, 100))


def iqr_sigma(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 4:
        return np.nan
    q1, q3 = np.percentile(a, [25, 75])
    return (q3 - q1) / 1.349


def d_nearest_two(x, positions):
    dists = np.abs(np.asarray(positions) - x)
    idx = np.argsort(dists)
    return float(dists[idx[0]]), float(dists[idx[1]])


rows = []
header = f"{'x':>6}  {'Npe_END/face':>13}  {'Npe_TOP_mean':>13}  {'Npe_TOP_med':>12}  {'Npe_TOP_IQRσ':>13}  {'d_near':>8}  {'d_sec':>8}"
print(header)
print("-" * len(header))

for x in x_scan:
    # N_pe^TOP from hybrid scan (face_type 2, all TOP SiPMs pooled, raw count)
    path_sp = INDIR_SP / f"photon_hits_ej230_ntop20_x{x}.root"
    with uproot.open(path_sp) as f:
        a = f['sipm_hits'].arrays(['event_id', 'face_type'], library='np')
    ev_sp = a['event_id']; ft_sp = a['face_type']
    ev_uniq_sp = np.unique(ev_sp)
    n_ev_sp = len(ev_uniq_sp)
    mask_top = ft_sp == 2
    npe_top = np.bincount(
        np.searchsorted(ev_uniq_sp, ev_sp[mask_top]),
        minlength=n_ev_sp
    ).astype(float)

    # N_pe^END from HYBRID scan (face_type 0+1, raw count / 2 = per-face mean)
    mask_end_sp = (ft_sp == 0) | (ft_sp == 1)
    npe_end_hybrid_pf = np.bincount(
        np.searchsorted(ev_uniq_sp, ev_sp[mask_end_sp]),
        minlength=n_ev_sp
    ).astype(float) / 2.0

    # N_pe^END from END-only scan (face_type 0+1, raw count / 2 = per-face mean)
    path_eo = INDIR_EO / f"photon_hits_ej230_x{x}.root"
    with uproot.open(path_eo) as f:
        a = f['sipm_hits'].arrays(['event_id', 'face_type'], library='np')
    ev_eo = a['event_id']; ft_eo = a['face_type']
    ev_uniq_eo = np.unique(ev_eo)
    n_ev_eo = len(ev_uniq_eo)
    mask_end = (ft_eo == 0) | (ft_eo == 1)
    npe_end_total = np.bincount(
        np.searchsorted(ev_uniq_eo, ev_eo[mask_end]),
        minlength=n_ev_eo
    ).astype(float)
    npe_end_per_face = npe_end_total / 2.0  # per-face mean (consistent with 701.3 at x=0)

    d_near, d_sec = d_nearest_two(x, TOP_XPOS)

    eo_mean = float(np.mean(npe_end_per_face))
    hyb_mean = float(np.mean(npe_end_hybrid_pf))
    row = dict(
        x=x,
        npe_end_per_face=eo_mean,
        npe_top_mean=float(np.mean(npe_top)),
        npe_top_median=float(np.median(npe_top)),
        npe_top_iqr_sigma=float(iqr_sigma(npe_top)),
        d_nearest=d_near,
        d_second=d_sec,
        npe_end_hybrid=hyb_mean,
        npe_end_endonly=eo_mean,
        delta_pe=hyb_mean - eo_mean,
        delta_pct=(hyb_mean - eo_mean) / eo_mean * 100.0,
    )
    rows.append(row)
    print(f"{x:>6}  {row['npe_end_per_face']:>13.2f}  {row['npe_top_mean']:>13.2f}  "
          f"{row['npe_top_median']:>12.2f}  {row['npe_top_iqr_sigma']:>13.2f}  "
          f"{d_near:>8.1f}  {d_sec:>8.1f}")

# Gate check
npe_end_x0 = next(r['npe_end_per_face'] for r in rows if r['x'] == 0)
gate_ok = abs(npe_end_x0 - 701.3) <= 1.0
print(f"\nGate: N_pe^END/face(x=0) = {npe_end_x0:.2f}  (expected 701.3 ± 1)  {'PASS' if gate_ok else 'FAIL'}")
if not gate_ok:
    raise SystemExit(f"GATE FAILED: N_pe^END/face(x=0) = {npe_end_x0:.2f} outside 700.3–702.3")

# Write CSV
csv_path = OUTDIR / "top_npe_diag.csv"
fieldnames = ['x', 'npe_end_per_face', 'npe_top_mean', 'npe_top_median',
              'npe_top_iqr_sigma', 'd_nearest', 'd_second',
              'npe_end_hybrid', 'npe_end_endonly', 'delta_pe', 'delta_pct']
with open(csv_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
print(f"CSV: {csv_path}")

# Write meta.json
meta = dict(
    script="top_npe_diag.py",
    datasets=dict(
        npe_top="scan_end_vik_sparse_top_v2 (EJ-230, N_TOP=20, face_type=2)",
        npe_end_endonly="scan_end_vikuiti (EJ-230, END-only, face_type 0+1, per-face mean)",
        npe_end_hybrid="scan_end_vik_sparse_top_v2 (EJ-230, N_TOP=20, face_type 0+1, per-face mean)",
    ),
    scintillator="EJ-230",
    ntop_hybrid=20,
    x_positions=x_scan,
    n_events_sp=int(n_ev_sp),
    n_events_eo=int(n_ev_eo),
    columns=dict(
        x="muon gun x position [mm]",
        npe_end_per_face="mean hits per END face per event, from scan_end_vikuiti (raw count/2)",
        npe_top_mean="mean total TOP hits per event, all 20 SiPMs pooled (scan_end_vik_sparse_top_v2)",
        npe_top_median="median total TOP hits per event",
        npe_top_iqr_sigma="IQR-based sigma of total TOP hits per event",
        d_nearest="distance from x to nearest TOP SiPM center [mm]",
        d_second="distance from x to second-nearest TOP SiPM center [mm]",
        npe_end_hybrid="mean hits per END face per event from hybrid scan (scan_end_vik_sparse_top_v2, N_TOP=20)",
        npe_end_endonly="same as npe_end_per_face — mean hits per END face per event from END-only scan",
        delta_pe="npe_end_hybrid - npe_end_endonly [hits/face/event]: effect of TOP SiPMs on END count",
        delta_pct="delta_pe / npe_end_endonly * 100 [%]",
    ),
    gate=dict(
        quantity="npe_end_per_face at x=0 (scan_end_vikuiti)",
        expected=701.3,
        tolerance=1.0,
        actual=round(npe_end_x0, 4),
        status="PASS" if gate_ok else "FAIL",
    ),
    top_sipm_positions=TOP_XPOS.tolist(),
)
meta_path = OUTDIR / "top_npe_diag_meta.json"
with open(meta_path, 'w') as f:
    json.dump(meta, f, indent=2)
print(f"Meta: {meta_path}")
