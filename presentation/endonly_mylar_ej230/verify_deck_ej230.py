#!/usr/bin/env python3
"""
verify_deck_ej230.py
Verify that build_deck_ej230.py output (main.tex) satisfies all physics
and editorial constraints from the ADDENDUM FINAL.

Checks:
  1. No unresolved @PLACEHOLDER@ tokens remain
  2. Forbidden phrases not present (addendum editorial corrections)
  3. Required phrases are present
  4. sigma_eq used (not sigma_single in display labels)
  5. Bootstrap metadata present and correct
  6. deck_values.json keys present

Exit code: 0 = all pass, 1 = failures found.
"""

import json
import pathlib
import re
import sys

HERE     = pathlib.Path(__file__).parent
TEX_PATH = HERE / "main.tex"
JSON_PATH= HERE / "deck_values.json"


# ── Helpers ────────────────────────────────────────────────────────────────────

def load_tex() -> str:
    if not TEX_PATH.exists():
        print(f"FAIL: main.tex not found at {TEX_PATH}")
        sys.exit(1)
    return TEX_PATH.read_text()


def load_json() -> dict:
    if not JSON_PATH.exists():
        print(f"FAIL: deck_values.json not found at {JSON_PATH}")
        sys.exit(1)
    with open(JSON_PATH) as f:
        return json.load(f)


def check_no_placeholders(tex: str) -> list[str]:
    found = re.findall(r"@[A-Z0-9_]+@", tex)
    if found:
        return [f"Unresolved placeholder: {p}" for p in sorted(set(found))]
    return []


def check_forbidden(tex: str) -> list[str]:
    """Addendum editorial corrections: these phrases must NOT appear."""
    forbidden = [
        # Item 7: M3 offset C must not be called 'recirculated floor'
        (r"recirculated floor",
         "M3 offset C must use 'finite-range offset', not 'recirculated floor'"),
        # Item 6: EndTop comparison
        (r"due to better reflector",
         "EndTop comparison: don't say 'due to better reflector' — use 'consistent with higher-reflectivity boundaries'"),
        # Item 1: M2 description
        (r"good fit",
         "M2 must not be called a 'good fit'; chi2/ndf=101>>1"),
        # Item 8: photon budget — 'dominates' near absorption
        (r"Mylar absorption dominates",
         "Photon budget: census can't separate surface from bulk; don't say 'Mylar absorption dominates'"),
        # Item 4: FPT causal language
        (r"inflat",
         "FPT: avoid 'inflate' — use 'estimator-dependent apparent longitudinal response'"),
        # Item 10: sigma_single must not appear as a display label
        (r"\\sigma_\{\\rm single\}",
         "Use sigma_eq not sigma_single in LaTeX math labels"),
        (r"sigma.single.ps.*label",
         "sigma_single label must not appear; use sigma_eq"),
        # Item 2: covariance vs bootstrap
        (r"covariance.*bootstrap",
         "Distinguish covariance errors and bootstrap CIs explicitly — do not conflate them in a single phrase (covariance.*bootstrap)"),
        # Addendum extra: EJ-204 material name must not appear in EJ-230 report
        (r"EJ-204",
         "EJ-204 must not appear in EJ-230 report"),
        (r"OPSC-101",
         "OPSC-101 must not appear in EJ-230 report"),
    ]
    errors = []
    for pattern, message in forbidden:
        if re.search(pattern, tex, re.IGNORECASE):
            errors.append(f"FORBIDDEN: {message}")
    return errors


def check_required(tex: str) -> list[str]:
    """Required phrases that must be present."""
    required = [
        # Item 7: must use 'finite-range offset'
        (r"finite.range offset",
         "M3 offset C must use 'finite-range offset'"),
        # Item 6: EndTop comparison language
        (r"configuration.level",
         "EndTop comparison must say 'configuration-level result'"),
        # Item 1: chi2/ndf cited for M2
        (r"\\ChiMtwo|chi.*?ndf.*?\\texttt|chi.?2.*?ndf",
         "M2 chi2/ndf must be cited in deck"),
        # Item 3: bootstrap metadata
        (r"pairs.bootstrap",
         "Bootstrap must be documented as 'pairs bootstrap' or 'pairs-bootstrap'"),
        # Item 10: sigma_eq definition
        (r"sigma_\{\\rm eq\}|\\sigma_{\rm eq}|SigEq",
         "sigma_eq (as sigma_eq or SigEq macro) must appear in deck"),
        # Item 9: escape category in motivation
        (r"escape|escaped",
         "Motivation must mention photons that 'escape' the bar"),
        # Item 4: estimator-dependent language
        (r"estimator.dependent",
         "Must use 'estimator-dependent' (not causal 'inflating') for v_app"),
        # Item 5: both NPE metrics
        (r"center.point|central.avg|central-avg|central.region|NpeCentral",
         "Must distinguish center-point NPE from central-region average"),
        # Item 12: bootstrap CI in summary
        (r"bootstrap.*?CI|68.*?CI|P16.*?P84",
         "Summary must cite bootstrap CI"),
        # Item 11: M2 full params
        (r"A_s|A_l|AmpShort|AmpLong|\\AmpShort|\\AmpLong",
         "M2 full parameters (A_s, A_l) must appear in four-models table"),
    ]
    errors = []
    for pattern, message in required:
        if not re.search(pattern, tex, re.IGNORECASE):
            errors.append(f"MISSING: {message}")
    return errors


def check_json_keys(v: dict) -> list[str]:
    required_keys = [
        "scint_material",
        "lambda_bulk_cm",
        "peak_emission_nm",
        "refractive_index",
        "mylar_reflectivity",
        "mylar_specular_lobe",
        "mylar_sigma_alpha_deg",
        "pde_peak_frac",
        "att_left",
        "att_right",
        "sigma_t_at_key_positions",
        "photon_budget_pct_surface_bulk",
        "center_point_npe_left",
        "central_region_npe_left",
        "endtop_center_point_npe_left",
        "endtop_central_region_npe_left",
        "bootstrap_n_replicas",
        "bootstrap_engine",
        "bootstrap_M4_left_median_cm",
        "bootstrap_M4_left_p16_cm",
        "bootstrap_M4_left_p84_cm",
    ]
    missing = []
    for k in required_keys:
        if k not in v:
            missing.append(f"Missing deck_values.json key: {k}")
    # Check bootstrap_cpp exists in att_left
    if "att_left" in v and "bootstrap_cpp" not in v["att_left"]:
        missing.append("Missing att_left.bootstrap_cpp (C++ bootstrap results)")
    # Check scint_material is EJ-230 not EJ-204
    sm = v.get("scint_material", "")
    if "EJ-204" in sm or "OPSC-101" in sm:
        missing.append(f"scint_material is '{sm}': must be EJ-230/OPSC-106")
    return missing


def check_figures(pres_dir: pathlib.Path) -> list[str]:
    required_figs = [
        "fig_npe_vs_x.png",
        "fig_attenuation.png",
        "fig_attenuation_twocomp.png",
        "fig_endtop_comparison.png",
        "fig_far_end_yield.png",
        "fig_timing_histograms.png",
        "fig_fpt_vs_t50.png",
        "fig_mean_arrival_time.png",
        "fig_sigma_t.png",
        "fig_photon_budget.png",
    ]
    missing = []
    for fig in required_figs:
        p = pres_dir / "figures" / fig
        if not p.exists():
            missing.append(f"Missing figure: {p}")
    return missing


def main():
    tex = load_tex()
    v   = load_json()

    all_errors = []

    errors = check_no_placeholders(tex)
    all_errors.extend(errors)

    errors = check_forbidden(tex)
    all_errors.extend(errors)

    errors = check_required(tex)
    all_errors.extend(errors)

    errors = check_json_keys(v)
    all_errors.extend(errors)

    errors = check_figures(HERE)
    all_errors.extend(errors)

    if all_errors:
        print(f"VERIFY FAILED: {len(all_errors)} issue(s) found\n")
        for e in all_errors:
            print(f"  {e}")
        sys.exit(1)
    else:
        print("VERIFY PASSED: all checks OK")
        print(f"  TEX:  {TEX_PATH}")
        print(f"  JSON: {JSON_PATH}")
        print(f"  EJ-230 scint_material: {v.get('scint_material')}")
        print(f"  M4 bootstrap: {v.get('bootstrap_M4_left_median_cm', '?'):.2f} "
              f"[{v.get('bootstrap_M4_left_p16_cm', '?'):.2f}, "
              f"{v.get('bootstrap_M4_left_p84_cm', '?'):.2f}] cm")
        sys.exit(0)


if __name__ == "__main__":
    main()
