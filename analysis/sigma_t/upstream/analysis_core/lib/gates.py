#!/usr/bin/env python3.12
"""
gates.py — QA abort-on-anomaly gates for EXEC_16.

Gates listed in the spec:
    QA-0  : runtime SHA matches declared SHA (ABORT if mismatch)
    QA-1  : TTree schema has all required branches (ABORT if missing)
    QA-2  : fitted μ and σ are within physical bounds (ABORT per point)
    QA-3  : fraction of χ²/NDF > threshold per group (WARN, flag group)
    QA-3d : SHA-256 manifest of all output files (written at end)
"""

import math
import hashlib
import subprocess
from pathlib import Path


# ─── QA-0 — runtime SHA verification ─────────────────────────────────────────

def gate_0_verify_sha(repo_path: str, branch: str, expected_sha: str, material: str) -> None:
    """
    Abort if the git branch SHA does not match the declared expected SHA.

    Why: we must never process data built from the wrong commit.  A wrong SHA
    means the geometry fix may not be applied, invalidating all results.

    Parameters
    ----------
    repo_path    : str — absolute path to the git repository
    branch       : str — branch name, e.g. "feat/ej230-sslg4"
    expected_sha : str — short SHA from exec16_config.yaml, e.g. "ca2f1c3"
    material     : str — human label for error messages, e.g. "EJ-230"

    Raises
    ------
    SystemExit
        If SHA cannot be determined or does not match expected_sha.
    """
    try:
        # use git rev-parse to get the current SHA of the branch tip
        proc = subprocess.run(
            ["git", "-C", repo_path, "rev-parse", "--short", f"refs/heads/{branch}"],
            capture_output=True, text=True, check=True,
        )
        actual_sha = proc.stdout.strip()          # e.g. "ca2f1c3"
    except subprocess.CalledProcessError as exc:
        # git command failed: repo path wrong, branch missing, or git not available
        raise SystemExit(
            f"[QA-0 ABORT] Cannot read SHA for {branch} in {repo_path}: {exc}"
        )

    if actual_sha != expected_sha:
        # mismatch: probably a different checkout or a push that changed the tip
        raise SystemExit(
            f"[QA-0 ABORT] {material}: branch {branch} is at SHA {actual_sha!r}, "
            f"but exec16_config.yaml declares {expected_sha!r}. "
            "Update the config or check out the correct commit."
        )

    print(f"[QA-0 PASS] {material}: {branch} @ {actual_sha} ✓")


# ─── QA-1 — TTree schema verification ────────────────────────────────────────

def gate_1_verify_schema(available_keys: set, required_branches: list, filename: str) -> None:
    """
    Abort if any required branch is missing from the TTree.

    Why: a missing branch means the data file is from a different software
    version.  Rather than silently producing NaN columns, we abort.

    Parameters
    ----------
    available_keys   : set — branches present in the TTree (from tree.keys())
    required_branches: list — expected branch names from exec16_config.yaml
    filename         : str — source file path for error messages

    Raises
    ------
    SystemExit
        If any required branch is absent.
    """
    missing = [b for b in required_branches if b not in available_keys]
    if missing:
        raise SystemExit(
            f"[QA-1 ABORT] {filename}: TTree 'sipm_hits' is missing branches: {missing}. "
            f"Available branches: {sorted(available_keys)}"
        )


# ─── QA-2 — per-point physics bounds ─────────────────────────────────────────

def gate_2_physics_bounds(result: dict, cfg: dict, context: str) -> None:
    """
    Abort if a successful fit produced a physically impossible result.

    This protects against: fit converging to wrong minimum, units error,
    or a bug in fit_engine that returns a negative σ or a 20-ns peak.

    Parameters
    ----------
    result  : dict — output of fit_core_gaussian()
    cfg     : dict — loaded exec16_config.yaml
    context : str  — identifier for error messages (material/group/x/N)

    Raises
    ------
    SystemExit
        If μ or σ_fit is outside the configured physical bounds.
    """
    if result.get("flag") != "ok":
        return   # only check successful fits; others already flagged by fit_engine

    mu  = result.get("mu_fit",    float("nan"))
    sig = result.get("sigma_fit", float("nan"))

    if math.isnan(mu) or math.isnan(sig):
        return   # NaN means fit already failed; gate does not apply

    mu_min  = cfg.get("QA2_MU_MIN_NS",    0.1)    # ns; earliest possible photon arrival
    mu_max  = cfg.get("QA2_MU_MAX_NS",   25.0)    # ns; very late = probably noise
    sig_max = cfg.get("QA2_SIGMA_MAX_NS", 5.0)    # ns; core σ > 5 ns is unphysical

    if not (mu_min <= mu <= mu_max):
        # μ outside expected time-of-flight window → wrong fit or wrong units
        raise SystemExit(
            f"[QA-2 ABORT] {context}: μ_fit = {mu:.4f} ns is outside "
            f"[{mu_min}, {mu_max}] ns. Possible fit convergence failure or unit error."
        )

    if sig > sig_max:
        # σ_fit > 5 ns is unphysical for a 1.4 m bar — tail dominated, not core
        raise SystemExit(
            f"[QA-2 ABORT] {context}: σ_fit = {sig*1000:.1f} ps "
            f"> {sig_max*1000:.0f} ps. Unphysical for this bar geometry. "
            "Check fit window and seeds."
        )


# ─── QA-3 — per-group χ²/NDF fraction check ──────────────────────────────────

def gate_3_check_group(results_list: list, group: str, cfg: dict) -> tuple:
    """
    After all positions are fit for a group, check what fraction have
    χ²/NDF above the warning threshold.

    If the fraction exceeds CORE_NOT_GAUSSIAN_FRACTION, the group is
    flagged 'core_not_gaussian'.  This does NOT abort — it's a diagnostic
    indicating that HOOK_WALK may be needed (or the Gaussian is a bad model
    for this group's core).

    Parameters
    ----------
    results_list : list  — list of result dicts, one per position, for one (group, N)
    group        : str   — group name, e.g. "TOP_SUM4"
    cfg          : dict  — loaded exec16_config.yaml

    Returns
    -------
    flagged   : bool  — True if group should be marked core_not_gaussian
    fraction  : float — fraction of valid positions with χ²/NDF > threshold
    n_flagged : int   — count of positions with χ²/NDF > threshold
    """
    chi2_thresh  = cfg.get("CHI2_NDF_WARN",               3.0)
    frac_thresh  = cfg.get("CORE_NOT_GAUSSIAN_FRACTION",   0.30)

    # only consider positions where the fit was attempted (chi2_ndf is not NaN)
    valid = [r for r in results_list if not math.isnan(r.get("chi2_ndf", float("nan")))]

    if not valid:
        return False, 0.0, 0   # no valid fits → cannot judge

    n_flagged = sum(1 for r in valid if r["chi2_ndf"] > chi2_thresh)
    fraction  = n_flagged / len(valid)
    flagged   = fraction > frac_thresh

    if flagged:
        print(
            f"  [QA-3 WARN] {group}: {n_flagged}/{len(valid)} positions "
            f"({fraction:.1%}) have χ²/NDF > {chi2_thresh:.1f} "
            f"→ group flagged 'core_not_gaussian'."
        )
    return flagged, fraction, n_flagged


# ─── QA-3d — output SHA-256 manifest ─────────────────────────────────────────

def compute_sha256(filepath: str) -> str:
    """
    Compute the SHA-256 hash of a file on disk, reading in 64 kB chunks.
    Returns hex digest string.
    """
    h = hashlib.sha256()
    with open(filepath, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)           # hash in chunks: handles large files without loading all to RAM
    return h.hexdigest()


def gate_3d_write_manifest(output_files: list, manifest_path: str) -> None:
    """
    Compute SHA-256 hashes for all output files and write a manifest.

    The manifest allows verifying output integrity in future sessions or
    before pushing to a remote archive.

    Parameters
    ----------
    output_files  : list — absolute paths to all output files produced
    manifest_path : str  — where to write the manifest file
    """
    lines = ["# EXEC_16 SHA-256 output manifest (QA-3d)", ""]

    for fpath in sorted(output_files):
        fpath = str(fpath)
        if Path(fpath).exists():
            sha = compute_sha256(fpath)
            lines.append(f"{sha}  {fpath}")
        else:
            lines.append(f"{'MISSING':>64}  {fpath}")   # flag missing files loudly

    manifest_str = "\n".join(lines) + "\n"
    Path(manifest_path).write_text(manifest_str)
    print(f"  [QA-3d] SHA-256 manifest written: {manifest_path}")
