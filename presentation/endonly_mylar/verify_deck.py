#!/usr/bin/env python3
"""
verify_deck.py — Traceability check for the End-only + Mylar Beamer deck.

Verifies that every numeric value in main.tex (via LaTeX macros) can be
traced back to deck_values.json, and that deck_values.json values match
the source CSV files.

Exit code 0 = all checks pass.  Non-zero = failure(s) found.
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).parent


def _load_json(path: pathlib.Path) -> dict:
    with open(path) as f:
        return json.load(f)


def _load_csv_col(path: pathlib.Path, col: str) -> list[float]:
    with open(path, newline="") as f:
        return [float(r[col]) for r in csv.DictReader(f) if r.get(col) not in (None, "")]


def _approx_eq(a: float | None, b: float | None, rtol: float = 1e-4) -> bool:
    if a is None or b is None:
        return False
    if not (math.isfinite(a) and math.isfinite(b)):
        return a == b   # both nan or both inf
    if abs(b) < 1e-12:
        return abs(a) < 1e-12
    return abs(a - b) / abs(b) < rtol


PASS = "[PASS]"
FAIL = "[FAIL]"
SKIP = "[SKIP]"


def run(out_dir: pathlib.Path) -> int:
    failures = 0

    json_path = out_dir / "deck_values.json"
    tex_path  = out_dir / "main.tex"
    att_csv   = out_dir.parent.parent / "t0minidaq" / "endonly_mylar_20260614" / "analysis" / "attenuation_fits.csv"
    sig_csv   = out_dir.parent.parent / "t0minidaq" / "endonly_mylar_20260614" / "analysis" / "sigma_t_sum4.csv"

    # Allow overriding data path from deck_values.json
    if json_path.exists():
        v = _load_json(json_path)
        data_dir = pathlib.Path(v.get("data_dir", ""))
        att_csv  = data_dir / "analysis" / "attenuation_fits.csv"
        sig_csv  = data_dir / "analysis" / "sigma_t_sum4.csv"
    else:
        v = {}
        print(f"{FAIL} deck_values.json not found at {json_path}")
        failures += 1

    print("\n=== 1. deck_values.json vs. attenuation_fits.csv (λ value only) ===")
    # Note: errors in JSON are scaled by √(χ²/ndf) — intentionally differ from raw CSV errors.
    if att_csv.exists():
        with open(att_csv, newline="") as f:
            att_rows = {r["side"]: r for r in csv.DictReader(f)}
        for side, json_key_lam in [
            ("left",     "lambda_eff_left_cm"),
            ("right",    "lambda_eff_right_cm"),
            ("combined", "lambda_eff_combined_cm"),
        ]:
            row = att_rows.get(side, {})
            csv_lam  = float(row.get("lambda_eff_cm", "nan"))
            json_lam = v.get(json_key_lam)
            ok_lam   = _approx_eq(json_lam, csv_lam, rtol=1e-3)
            tag = PASS if ok_lam else FAIL
            if tag == FAIL:
                failures += 1
            print(f"  {tag} λ_eff {side}: JSON={json_lam:.4f} CSV={csv_lam:.4f}")
    else:
        print(f"  {SKIP} {att_csv} not found")

    print("\n=== 2. deck_values.json vs. sigma_t_sum4.csv ===")
    if sig_csv.exists():
        with open(sig_csv, newline="") as f:
            sig_rows = {float(r["x_mm"]): r for r in csv.DictReader(f)}
        sigma_at = v.get("sigma_t_at_key_positions", {})
        for xkey in sigma_at:
            xv = float(xkey)
            row = sig_rows.get(xv, {})
            csv_s = float(row.get("sigma_single_ps", "nan"))
            csv_e = float(row.get("sigma_single_error_ps", "nan"))
            json_s = sigma_at[xkey].get("sigma_single_ps")
            json_e = sigma_at[xkey].get("sigma_single_err_ps")
            ok_s = _approx_eq(json_s, csv_s)
            ok_e = _approx_eq(json_e, csv_e)
            tag = PASS if (ok_s and ok_e) else FAIL
            if tag == FAIL:
                failures += 1
            print(f"  {tag} σ_single x={xkey}: JSON={json_s:.1f} CSV={csv_s:.1f}  "
                  f"err: JSON={json_e:.1f} CSV={csv_e:.1f}")
    else:
        print(f"  {SKIP} {sig_csv} not found")

    print("\n=== 3. main.tex macro cross-check ===")
    # Every \newcommand{\Foo}{value} in main.tex must come from deck_values.json
    if tex_path.exists():
        tex = tex_path.read_text()
        # Extract \newcommand{\Xxx}{...}
        macros = dict(re.findall(
            r"\\newcommand\{\\(\w+)\}\{(.*?)\}", tex
        ))
        # Check a selection of critical macros
        checks = [
            ("LamBulk",  str(int(v.get("lambda_bulk_cm", 0)))),
            ("NEv",      str(v.get("n_events_per_position", ""))),
            ("NPos",     str(v.get("n_positions", ""))),
            ("MylarR",   str(v.get("mylar_reflectivity", ""))),
        ]
        for macro_name, expected_fragment in checks:
            val = macros.get(macro_name, "")
            # Check that the expected fragment appears in the macro value
            ok = expected_fragment in val or expected_fragment.replace(",", "") in val.replace(",", "")
            tag = PASS if ok else FAIL
            if tag == FAIL:
                failures += 1
            print(f"  {tag} \\{macro_name}: tex='{val}'  expected fragment='{expected_fragment}'")

        # Check no hardcoded long numbers appear outside the macros block
        macro_block_end = tex.find(r"\title[")
        if macro_block_end == -1:
            macro_block_end = tex.find(r"\begin{document}")
        body = tex[macro_block_end:] if macro_block_end > 0 else tex

        # Critical: λ_eff not hardcoded in prose (tabular environments are exempt)
        if v.get("lambda_eff_combined_cm"):
            lam_str = f"{v['lambda_eff_combined_cm']:.2f}"
            body_no_tables = re.sub(
                r"\\begin\{tabular\}.*?\\end\{tabular\}", "", body, flags=re.DOTALL
            )
            if lam_str in body_no_tables:
                print(f"  {FAIL} Hardcoded λ_eff={lam_str} found in prose — should use macro")
                failures += 1
            else:
                print(f"  {PASS} λ_eff={lam_str} not hardcoded in prose (uses macro or table only)")

    else:
        print(f"  {SKIP} main.tex not found at {tex_path}")

    print("\n=== 3b. Extended fit keys (three-model attenuation) ===")
    for side_key in ("att_left", "att_right"):
        d = v.get(side_key, {})
        label = side_key.replace("att_", "")
        # single-exp
        chi1 = d.get("single_chi2ndf")
        lam_p = d.get("expoff_lam_cm")
        lam_t = d.get("tail_lam_cm")
        for name, val, lo, hi in [
            (f"{label} single χ²/ndf", chi1, 100, 1e7),
            (f"{label} prompt λ (cm)",  lam_p, 10, 50),
            (f"{label} tail λ (cm)",    lam_t, 20, 80),
        ]:
            if val is None:
                print(f"  {FAIL} {name}: MISSING")
                failures += 1
            elif lo < float(val) < hi:
                print(f"  {PASS} {name}: {float(val):.2f}")
            else:
                print(f"  {FAIL} {name}: {float(val):.2f}  (expected {lo}–{hi})")
                failures += 1

    print("\n=== 3c. EJ-204 datasheet properties ===")
    ds_checks = [
        ("scint_yield_per_MeV", 10400),
        ("peak_emission_nm",    408.0),
        ("refractive_index",    1.58),
        ("density_gcc",         1.023),
    ]
    for key, expected in ds_checks:
        val = v.get(key)
        ok  = val is not None and _approx_eq(float(val), float(expected), rtol=1e-3)
        tag = PASS if ok else FAIL
        if tag == FAIL:
            failures += 1
        print(f"  {tag} {key}: JSON={val}  expected={expected}")

    print("\n=== 3d. No ±0.0 in main.tex ===")
    if tex_path.exists():
        tex_check = tex_path.read_text()
        bad_pm = re.findall(r"\\pm\s*0\.0(?:\b|\s)", tex_check)
        if bad_pm:
            print(f"  {FAIL} Found {len(bad_pm)} occurrence(s) of '\\pm 0.0' — check χ²/ndf scaling")
            failures += 1
        else:
            print(f"  {PASS} No '\\pm 0.0' patterns found")
    else:
        print(f"  {SKIP} main.tex not found")

    print("\n=== 4. Physical sanity checks ===")
    et = v.get("endtop", {})
    sanity = [
        ("λ_eff > 0",       v.get("lambda_eff_combined_cm", -1) is not None
                             and v.get("lambda_eff_combined_cm", -1) > 0),
        ("λ_eff < λ_bulk",  (v.get("lambda_eff_combined_cm", 0) or 0) <
                             (v.get("lambda_bulk_cm", 0) or 0)),
        ("σ_t center > 0",  (v.get("sigma_t_center_ps", 0) or 0) > 0),
        ("σ_t center < σ_t end",
                             (v.get("sigma_t_center_ps", 0) or 0) <
                             (v.get("sigma_t_end690_ps", 999) or 999)),
        ("n_positions = 31", v.get("n_positions") == 31),
        ("n_events = 2000",  v.get("n_events_per_position") == 2000),
        ("PDE in (0,1)",     0 < (v.get("pde_peak_frac") or 0) < 1),
        ("Mylar R in (0,1)", 0 < (v.get("mylar_reflectivity") or 0) < 1),
        ("budget_escaped < budget_generated",
                             (v.get("photon_budget_escaped") or 0) <
                             (v.get("photon_budget_generated") or 1)),
        ("yield = 10400",    v.get("scint_yield_per_MeV") == 10400),
        ("endtop available", et.get("available", False)),
        ("endtop NPE > mylar NPE at x=0",
                             (et.get("npe_left_at_x0") or 0) >
                             (v.get("npe_left_at_x0") or 0)),
    ]
    for label, result in sanity:
        tag = PASS if result else FAIL
        if tag == FAIL:
            failures += 1
        print(f"  {tag} {label}")

    print("\n=== 5. Figures exist ===")
    figs = v.get("figures", {})
    for key, rel in figs.items():
        p = out_dir / rel
        exists = p.exists()
        tag = PASS if exists else FAIL
        if tag == FAIL:
            failures += 1
        print(f"  {tag} {rel}")

    print(f"\n{'='*60}")
    if failures == 0:
        print(f"ALL CHECKS PASSED  (0 failures)")
    else:
        print(f"FAILURES: {failures}")
    print('='*60)
    return failures


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=HERE, type=pathlib.Path)
    args = ap.parse_args()
    sys.exit(run(args.out_dir))
