#!/usr/bin/env python3
"""Render the EXEC43 report and hash its analysis artifacts."""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path


MATERIALS = ["EJ-200", "EJ-204", "EJ-230"]
FACES = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
         "OPEN_+Y_OR_+/-X_UNRESOLVED"]


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def f(row, key):
    return float(row[key])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    results = json.loads((args.analysis / "results.json").read_text())
    optical = read_csv(args.analysis / "optical_properties.csv")
    absent = read_csv(args.analysis / "absent_photons.csv")
    closure = read_csv(args.analysis / "attenuation_closure.csv")
    directional = read_csv(args.analysis / "directional_path_test.csv")
    position = read_csv(args.analysis / "position_dependence.csv")
    prereg_hash = sha(args.preregistration)

    artifact_paths = sorted(path for path in args.analysis.rglob("*") if path.is_file()
                            and path.name != "artifact_manifest.json")
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "analysis_root": str(args.analysis),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": prereg_hash,
        "artifact_count_excluding_this_manifest": len(artifact_paths),
        "artifacts": {str(path): sha(path) for path in artifact_paths},
    }
    (args.analysis / "artifact_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n")

    m2 = results["M2"]
    m1 = results["M1"]
    central_failures = sum(f(row, "central_abs_residual_over_se") > 3 for row in closure)
    envelope_widths = [f(row, "envelope_width") for row in closure]
    central_z = [f(row, "central_abs_residual_over_se") for row in closure]
    lines = [
        "# EXEC_43 — mechanism of material-dependent first-encounter escape",
        "",
        "**Date:** 2026-09-14  ",
        "**Scope:** analysis of the 21 existing EXEC_42 ROOT files; no simulation was run.  ",
        "**Formal preregistered result:** **M1 CONFIRMED WITHIN PERSISTED-DATA LIMITS; M2 EXCLUDED.**  ",
        "**Scientific qualification:** the exact path distribution is not recoverable from the production schema. The preregistered WorldPV envelope closes 14/14 targets, while the central geometric assignment misses all 14 by more than 3 SE. M1 is strongly supported as the mechanism, but this dataset does not establish the precision closure needed to adopt H3′ as a golden test.",
        "",
        "## Decision table",
        "",
        "| Test | Prediction fixed before calculation | Measurement | Result |",
        "|---|---|---|---|",
        f"| M2 effective index | Exclude M2 if max−min n_eff < 1e-4 | n_eff=1.58 for all materials; spread {m2['n_effective_spread']:.3g} | **M2 EXCLUDED** |",
        f"| M1 directional path factor | median sec(theta) non-escape > escape for every evaluable nonsensor face | {m1['directional_comparisons_pass']}/{m1['directional_comparisons_total']} comparisons | **PASS** |",
        f"| M1 absent photons | EJ-200 < EJ-204 < EJ-230; both adjacent differences >3 paired SE at all x | {m1['adjacent_absent_pairs_pass']}/{m1['adjacent_absent_pairs_total']} adjacent comparisons | **PASS** |",
        f"| M1 quantitative closure | EJ-200 reweighting envelope contains EJ-204/EJ-230 at all x | {m1['closure_cells_pass']}/{m1['closure_cells_total']} targets | **FORMAL PASS** |",
        f"| Central closure diagnostic | WorldPV assigned to +Y at 30 mm | {central_failures}/14 outside 3 SE; residual magnitude {min(central_z):.2f}–{max(central_z):.2f} SE | **DOES NOT CLOSE** |",
        "| Exact absolute paths | Need creation/encounter coordinates or track length | absent for ±Z and unresolved WorldPV | **NOT EVALUABLE** |",
        f"| Overall M1 | M2 excluded plus directional, survival, and envelope closure tests pass | `{m1['status']}` | **SUPPORTED / FORMAL CONFIRMATION WITH LIMIT** |",
        "",
        "The formal and qualified statements are both required: the registered envelope criterion passes, but it is wider than the effect being tested and therefore cannot establish a precision prediction.",
        "",
        "## Provenance and method",
        "",
        "The analysis used the EXEC_42 inventory `/home/rrios/exec42_20260913/analysis/input_inventory.csv`, SHA-256 `012bb6e02f3cacdfd05c418ada06c031942032882aa4f9a886bdfadc36790ea7`. Every cell has N=10,000 and common simulation seeds 26092601 and 8349041. Uncertainties use 500 common generated-event bootstrap resamples, NumPy seed 43091401 and `ddof=1`.",
        "",
        f"The EXEC_43 preregistration SHA-256 is `{prereg_hash}`. It was committed before calculation as `3bc6ed6`; the frozen analyzer is commit `01549fe`.",
        "",
        "Exact command:",
        "",
        "```bash",
        "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/analyze_exec43.py --inventory /home/rrios/exec42_20260913/analysis/input_inventory.csv --exec42-results /home/rrios/exec42_20260913/analysis/contract_results.json --sslg4 /home/rrios/exec40_20260913/build_off/sslg4 --preregistration analysis/validation/EXEC43_PREREGISTRATION.md --out /home/rrios/exec43_20260914",
        "```",
        "",
        "The read-only pass traversed 4.244 billion stored first-encounter rows and completed in approximately 818 s (13 min 38 s). Every streamed source-1 exiting count closed exactly against the independent EXEC_42 per-event encounter counts.",
        "",
        "## Test 1 — effective refractive index",
        "",
        "Each SSLG4 RINDEX table contains the two endpoints `(200 nm, 1.58)` and `(800 nm, 1.58)`. Thus the implemented model has no optical dispersion over any of the three emission spectra.",
        "",
        "| Material | OPSC | Emission peak (nm) | Emission mean (nm) | Attenuation (mm) | n_eff | 1−cos(theta_c) | Measured ±Z conditional | Measured aggregate mean over x |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in optical:
        lines.append(
            f"| {row['material']} | {row['opsc_code']} | {f(row,'emission_peak_nm'):.3f} | "
            f"{f(row,'emission_mean_nm'):.3f} | {f(row,'attenuation_mm'):.0f} | "
            f"{f(row,'n_effective'):.12f} | {f(row,'geometric_escape'):.9f} | "
            f"{f(row,'measured_Z_conditional_escape'):.9f} | {f(row,'measured_escape_mean_x'):.9f} |")
    lines += [
        "",
        "EJ-230 has a two-point maximum plateau at 388.684 and 390.526 nm; the table reports the first maximum. The supplied nominal 391 nm is consistent with that plateau. The spectrum-weighted means differ substantially, but the weighted indices do not. M2 therefore predicts zero material shift in this implementation and is quantitatively excluded by the 1e-4 rule.",
        "",
        "The geometric cone value 0.225775763 agrees closely with the measured combined ±Z conditional fraction for EJ-200 (0.225746695). The aggregate escape is higher because -Y and unresolved open faces contribute conditional escape near 0.58–0.67; it is not expected to equal the one-face cone fraction.",
        "",
        "## Test 2 — path information and directional filtering",
        "",
        "The production `first_bar_encounters` tree contains `event_id, track_id, source, pre_copy, post_copy, pre_volume_id, post_volume_id, outcome, exiting_bar, cos_incidence`. It contains no creation position, encounter position, global time, or track length. Exact absolute first paths are therefore unavailable for the dominant ±Z faces, where the unknown production distance is 0–10 mm, and for WorldPV, whose face is unresolved.",
        "",
        "The measured proxy is `sec(theta)=1/cos_incidence`. The table gives the range of the per-cell medians across the seven x positions and total photon counts. A larger secant means a longer path for fixed perpendicular distance.",
        "",
        "| Material | Face | Non-escape median sec range | Escape median sec range | Non-escape count | Escape count | Absolute-path status |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for material in MATERIALS:
        for face in FACES:
            selected = [row for row in directional if row["material"] == material and row["face"] == face]
            non = [f(row, "non_escape_secant_median") for row in selected]
            esc = [f(row, "escape_secant_median") for row in selected
                   if row["escape_secant_median"] != "nan"]
            if face in ("-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR"):
                absolute = "fixed perpendicular distance"
            else:
                absolute = "NOT EVALUABLE"
            esc_text = f"{min(esc):.4f}–{max(esc):.4f}" if esc else "—"
            lines.append(
                f"| {material} | {face} | {min(non):.4f}–{max(non):.4f} | {esc_text} | "
                f"{sum(int(row['n_non_escape']) for row in selected):,} | "
                f"{sum(int(row['n_escape']) for row in selected):,} | {absolute} |")
    lines += [
        "",
        "All 84 evaluable comparisons on +Z, −Z, −Y and unresolved WorldPV have a larger non-escape median. On −Y, where d=30 mm is known exactly, the median absolute path is 45.7–48.9 mm for non-escape and 31.77–31.79 mm for escape. Sensor faces have no escaped population under the registered definition and cannot test the within-face prediction.",
        "",
        "The full per-cell distributions, including q10, median and q90, are stored in each `cells/<cell_id>/path_proxy_summary` sidecar; the 2,000-bin cosine distributions are in `path_proxy_histogram`.",
        "",
        "## Test 3 — photons absent before the first encounter",
        "",
        "`Absent = produced scintillation photons − recorded first physical BarLV encounters`. Fractions are estimate ± one event-bootstrap SE.",
        "",
        "| Cell | Produced | First encounters | Absent | Absent fraction |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in absent:
        lines.append(
            f"| {row['cell_id']} | {int(row['produced']):,} | {int(row['first_encounters']):,} | "
            f"{int(row['absent']):,} | {100*f(row,'absent_fraction'):.4f}% ± "
            f"{100*f(row,'absent_fraction_se'):.4f}% |")
    by_material = {material: [f(row, "absent_fraction") for row in absent
                              if row["material"] == material] for material in MATERIALS}
    lines += [
        "",
        f"The ranges are EJ-200 {100*min(by_material['EJ-200']):.4f}%–{100*max(by_material['EJ-200']):.4f}%, EJ-204 {100*min(by_material['EJ-204']):.4f}%–{100*max(by_material['EJ-204']):.4f}%, and EJ-230 {100*min(by_material['EJ-230']):.4f}%–{100*max(by_material['EJ-230']):.4f}%. The ordering follows 1/lambda_att at all positions. All 14 adjacent comparisons exceed 3 paired SE; their significances are 209–515 SE.",
        "",
        "This directly rejects the premise that the registered first-encounter population is purely geometric: material-dependent bulk survival changes its denominator before any surface is reached.",
        "",
        "## Test 4 — quantitative attenuation closure",
        "",
        "The prediction reweights the measured EJ-200 face/cosine population from 3800 mm to 1600 or 1200 mm. ±Z integrates an unrecorded uniform 0–10 mm production distance. WorldPV is evaluated at the three preregistered distance assignments; the minimum and maximum form the envelope. No attenuation length, distance, or normalization was fitted.",
        "",
        "Values below are shifts relative to the matching EJ-200 cell. The central prediction assigns WorldPV to +Y at 30 mm. `z` is the paired residual divided by its SE.",
        "",
        "| x (mm) | Target | Observed shift | Central predicted shift | Observed−central | abs(z) | Predicted-shift envelope | Formal envelope result |",
        "|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    event = json.loads(Path("/home/rrios/exec42_20260913/analysis/contract_results.json").read_text())["cells"]
    for row in closure:
        baseline = event[row["baseline_cell"]]["V2_H1"]["value"]
        observed_shift = f(row, "observed_escape") - baseline
        central_shift = f(row, "predicted_world_y") - baseline
        low_shift = f(row, "predicted_low") - baseline
        high_shift = f(row, "predicted_high") - baseline
        lines.append(
            f"| {row['x_mm']} | {row['target_material']} | {observed_shift:+.7f} | "
            f"{central_shift:+.7f} | {f(row,'central_residual'):+.7f} | "
            f"{f(row,'central_abs_residual_over_se'):.2f} | [{low_shift:+.7f}, {high_shift:+.7f}] | "
            f"**{row['status']}** |")
    lines += [
        "",
        f"The preregistered envelope accepts 14/14 targets. Its absolute width is {min(envelope_widths):.6f}–{max(envelope_widths):.6f}, much larger than the observed material shifts. The central model overpredicts every shift: EJ-204 residuals are −0.000219 to −0.000150 ({min(central_z):.2f}–{max(f(row,'central_abs_residual_over_se') for row in closure if row['target_material']=='EJ-204'):.2f} SE), and EJ-230 residuals are −0.000345 to −0.000259 ({min(f(row,'central_abs_residual_over_se') for row in closure if row['target_material']=='EJ-230'):.2f}–{max(central_z):.2f} SE).",
        "",
        "Thus the sign and scale are explained by attenuation, and the registered broad envelope closes, but the current data do not support a unique path-based prediction at the required sub-0.001 precision. The envelope pass must not be read as agreement of the central calculation.",
        "",
        "## Test 5 — position dependence",
        "",
        "The constant-difference test uses the covariance of the seven paired-bootstrap profile points. Position dependence is called resolved only at p<0.01.",
        "",
        "| Pair | Mean difference | Min–max | Peak-to-peak | Fraction of mean | chi2/dof | p | Conclusion |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in position:
        lines.append(
            f"| {row['material_a']} − {row['material_b']} | {f(row,'arithmetic_mean_difference'):+.7f} | "
            f"{f(row,'minimum_difference'):+.7f} to {f(row,'maximum_difference'):+.7f} | "
            f"{f(row,'peak_to_peak'):.7f} | {100*f(row,'fractional_peak_to_peak'):.1f}% | "
            f"{f(row,'chi_square'):.2f}/{row['dof']} | {f(row,'p_value'):.5f} | "
            f"{row['position_dependence']} |")
    lines += [
        "",
        "The absolute variation is weak (peak-to-peak 0.000167–0.000247), but only EJ-200−EJ-230 is statistically resolved (`p=0.00134`). EJ-200−EJ-204 (`p=0.227`) and EJ-204−EJ-230 (`p=0.0793`) are not resolved at the preregistered 0.01 threshold. Shared seeds mean these seven-point profiles remain correlated.",
        "",
        "## Proposed H3′ — not applied",
        "",
        "> At fixed geometry and source position, the material dependence of first-encounter escape is generated by pre-encounter bulk-absorption filtering. Starting from a declared common reference population with exact pre-encounter path length L, predict each target material by the registered ABSLENGTH(lambda) survival weight. For every material pair and position, require abs(Delta_observed − Delta_predicted) <= 3 paired event-bootstrap SE, with all populations, interpolation rules, reference direction, and any numerical systematic floor fixed before evaluation.",
        "",
        "For these SSLG4 materials ABSLENGTH is constant over wavelength, so the weight reduces to `exp[-L(1/lambda_target−1/lambda_reference)]`. A future production record must persist exact track length, or creation and first-encounter coordinates plus any scattering history. The `sec(theta)` proxy and WorldPV envelope used here are insufficient for golden acceptance at the observed precision.",
        "",
        "H3′ is a proposal only. `analysis/validation/GOLDEN_REFERENCE_20260913.json` was not edited, and H3 remains the active failed contract pending René's decision.",
        "",
        "## Artifacts and final gate",
        "",
        f"The analysis directory `/home/rrios/exec43_20260914` contains {len(artifact_paths)+1} files including `artifact_manifest.json`. All 55 CSV/ROOT/meta.json sidecar triplets were reopened and their hashes and row counts verified. `results.json`, `schema_audit.json`, and `analysis.log` preserve the decisions, limitation, and command output.",
        "",
        "No simulation, push, merge, deck edit, or golden-contract modification was performed.",
        "",
    ]
    args.out.write_text("\n".join(lines))
    print(args.out)


if __name__ == "__main__":
    main()
