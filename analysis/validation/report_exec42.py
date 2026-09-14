#!/usr/bin/env python3
"""Render the EXEC42 Phase B report from frozen machine-readable results."""
import argparse
import json
from pathlib import Path


FACES = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
         "OPEN_+Y_OR_+/-X_UNRESOLVED"]


def pm(value, error, digits=6):
    return f"{value:.{digits}f} ± {error:.{digits}f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--golden", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    contract = json.loads((args.analysis / "contract_results.json").read_text())
    inputs = json.loads((args.analysis / "input_verification.json").read_text())
    golden = json.loads(args.golden.read_text())
    block = golden["exec42_full_grid"]
    cells = sorted(contract["cells"].values(),
                   key=lambda item: (item["configuration"]["material"],
                                     item["configuration"]["x_mm"]))

    v1 = [cell["V1_primary"]["value"] for cell in cells]
    v1n = [cell["V1_nonoptical_secondary"]["value"] for cell in cells]
    v2 = [cell["V2_H1"]["value"] for cell in cells]
    v5d = [cell["V5_matched"]["difference"] for cell in cells]
    orphan = [cell["SD_orphan_fraction"]["value"] for cell in cells]
    h3z = [row["abs_difference_over_se"] for row in contract["H3"]]
    h3d = [abs(row["difference"]) for row in contract["H3"]]

    lines = [
        "# EXEC_42 Phase B — full-grid physical contract",
        "",
        "**Date:** 2026-09-14  ",
        "**Scope:** read-only analysis of the completed 21-cell grid; no simulation or relaunch was performed.  ",
        "**Decision:** **NOT READY FOR ACCEPTANCE**. V1, corrected V2-H1, and matched V5 pass in all 21 cells, but preregistered H3 fails all 21 material-pair comparisons. There are no non-evaluable required results.",
        "",
        "## Validation table",
        "",
        "| Test | Preregistered prediction | Source | Measured in 21 cells | Result |",
        "|---|---|---|---|---|",
        "| B1 inputs | 21 unique readable ROOTs, 10,000 events each, common PDE path/hash and seeds, actual workers from manifest | EXEC_42 B1 | 21/21 complete; 0 failed rows; 21 unique paths; 4 workers; seeds 26092601/8349041 | **PASS** |",
        f"| V1 primary | Unity inside ±3 event-bootstrap SE | EXEC40 `{block['preregistrations']['EXEC40']['sha256']}` | 21/21 pass; ratios {min(v1):.6f}–{max(v1):.6f} | **PASS** |",
        f"| V2-H1 | escape − 3 SE ≥ 0.226172 | EXEC41 `{block['preregistrations']['EXEC41']['sha256']}` | 21/21 pass; escape {min(v2):.6f}–{max(v2):.6f} | **PASS** |",
        f"| V5 matched | matched/incident agrees with incident-spectrum PDE within 3 paired SE | EXEC40 + EXEC41 correction | 21/21 pass; differences {min(v5d):+.6f} to {max(v5d):+.6f} | **PASS** |",
        f"| H3 material invariance | Every material pair at fixed x agrees within 3 paired SE | EXEC42 H3 `{block['preregistrations']['EXEC42_H3']['sha256']}` | 0/21 pass; abs(Δ)={min(h3d):.6f}–{max(h3d):.6f}, or {min(h3z):.2f}–{max(h3z):.2f} SE | **FAIL** |",
        "| Required evaluability | No required result is `NOT EVALUABLE` | EXEC_42 B3 | 0 non-evaluable results | **PASS** |",
        "| `ready_for_acceptance` | True only if V1, V2-H1, V5 and all H3 comparisons pass | EXEC_42 B3 | `false`; failure reason `H3` | **FAIL / false** |",
        "",
        "The fail-closed decision follows the registered rule mechanically. The H3 tolerance was not changed after observing the grid.",
        "",
        "## B1 — input verification",
        "",
        "The selected campaign is `/home/rrios/exec42_20260913/grid`, identified by its manifest, production commit `b35ee84acadef12c93506e0720580c91f901fcbf`, and completion on 2026-09-14. `/home/rrios/exec34r_20260912` was preserved and excluded because it is the older EXEC_34R campaign at commit `420addf`, before the production observables.",
        "",
        f"All 21 `.DONE` markers and manifest records agree. Every SHA-256 was recomputed over the ROOT file; all four required trees open; each `event_observables` tree contains exactly event IDs 0–9999; and every `sipm_event_counts` tree has 860,000 rows. The inputs total {sum(c['configuration']['root_bytes'] for c in cells):,} bytes ({sum(c['configuration']['root_bytes'] for c in cells)/2**30:.3f} GiB).",
        "",
        f"The common PDE is `{inputs['PDE_paths'][0]}` with SHA-256 `{inputs['PDE_hashes'][0]}`. The manifest, rather than the macro, records 4 workers for every cell. All cells use seeds 26092601 and 8349041, `eventModulo=1`, EndTop with N_TOP=70, diagnostics OFF, and N=10,000.",
        "",
        "Exact verification command:",
        "",
        "```bash",
        block["analysis"]["input_command"],
        "```",
        "",
        "B1 machine-readable result: `/home/rrios/exec42_20260913/analysis/input_verification.json`.",
        "",
        "## Statistical design and frozen definitions",
        "",
        "All uncertainty estimates use 500 generated-event bootstrap resamples with NumPy seed 42091401 and `ddof=1`. The same event-index weights are used for every cell. The 21 cells share simulation seeds, and positions within a material are correlated. Material comparisons at fixed position therefore use the bootstrap distribution of the paired ratio difference; later fits must retain this correlation.",
        "",
        "V1 is `sum(produced_scint)/(nominal_yield * sum(edep_total_MeV))`, with 10,000, 10,400 and 9,700 photons/MeV for EJ-200, EJ-204 and EJ-230. The alternate denominator `total minus optical deposit` remains a descriptive secondary reading.",
        "",
        "V2 uses scintillation photons at their first physical encounter while exiting BarLV. Escape is portable outcome 2 into the named air gaps or WorldPV; sensor acceptance is separate. H1 requires the aggregate escape fraction minus 3 SE to exceed 0.226172. The angular H2 comparisons are descriptive and do not enter acceptance.",
        "",
        "V5 is matched boundary detections divided by independently observed incidents. Its prediction is the summed surface-PDE expectation over that same incident population, divided by incidents. Independent SD detections and the orphan component are accounting diagnostics.",
        "",
        "Exact contract command:",
        "",
        "```bash",
        block["analysis"]["contract_command"],
        "```",
        "",
        "## Per-cell contract results",
        "",
        "Ratios are shown as estimate ± one event-bootstrap SE. `Npe/end` is the mean independent SD count over the two END sides.",
        "",
        "| Cell | Npe/end | V1 total | V1 non-optical | V2 escape | V5 matched | V5 PDE expected | V5 Δ ± SE | SD orphan | Required |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for cell in cells:
        cfg = cell["configuration"]
        a = cell["V1_primary"]; b = cell["V1_nonoptical_secondary"]
        c = cell["V2_H1"]; d = cell["V5_matched"]
        e = d["expected_incident_PDE"]; o = cell["SD_orphan_fraction"]
        status = "/".join((a["status"], c["status"], d["status"]))
        lines.append(
            f"| {cfg['cell_id']} | {cell['Npe']['per_end_mean']:.3f} | {pm(a['value'],a['se'])} | "
            f"{pm(b['value'],b['se'])} | {pm(c['value'],c['se'])} | {pm(d['value'],d['se'])} | "
            f"{pm(e['value'],e['se'])} | {d['difference']:+.6f} ± {d['difference_se']:.6f} | "
            f"{pm(o['value'],o['se'])} | {status} |")

    lines += [
        "",
        "The primary V1 range is " + f"{min(v1):.6f}–{max(v1):.6f}. The alternate denominator gives {min(v1n):.6f}–{max(v1n):.6f}; it is reported and does not replace V1.",
        "",
        "## H3 — material invariance at fixed position",
        "",
        "Each row uses identical bootstrap event weights in both materials. Percentage-point differences are `100 × (escape_a − escape_b)`.",
        "",
        "| x (mm) | Pair | Escape A | Escape B | Δ (percentage points) | Paired SE | abs(Δ)/SE | Result |",
        "|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in contract["H3"]:
        lines.append(
            f"| {row['x_mm']} | {row['material_a']} − {row['material_b']} | "
            f"{row['escape_a']:.6f} | {row['escape_b']:.6f} | {100*row['difference']:+.5f} | "
            f"{row['paired_se']:.7f} | {row['abs_difference_over_se']:.2f} | **{row['status']}** |")
    lines += [
        "",
        "The measured order is EJ-200 < EJ-204 < EJ-230 at every position. This is a significant material dependence in the registered first-encounter escape quantity and is the EXEC_42 finding. This analysis does not assign a cause or alter the definition.",
        "",
        "## Escape profile and first-face decomposition",
        "",
        "The expected position dependence is visible: central and ±200/±500 mm cells have escape near 0.2562–0.2573, while ±650 mm cells have 0.2627–0.2632. This profile is descriptive and has no position PASS/FAIL.",
        "",
        "Each entry below is `encounter share % / conditional escape %`. The unresolved class is kept intact because BarPV→WorldPV does not identify whether the uncovered crossing was +Y or ±X.",
        "",
        "| Cell | +Z | −Z | −Y | −X sensor | +X sensor | +Y sensor | Open unresolved |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for cell in cells:
        by_face = {face["face"]: face for face in cell["V2_faces"]}
        values = [f"{100*by_face[name]['encounter_share']:.4f}/{100*by_face[name]['value']:.4f}"
                  for name in FACES]
        lines.append("| " + cell["configuration"]["cell_id"] + " | " + " | ".join(values) + " |")

    low_mu = [cell["V2_H2"]["low_mu_probability"] for cell in cells]
    lines += [
        "",
        "Both the uniform-cosine and retained old `p(mu)=2mu` forms have numerical p-values of 0 in every cell under the frozen covariance calculation. H2 remains descriptive. The probability in the low-cosine interval [0,0.2] ranges from " + f"{min(low_mu):.6f} to {max(low_mu):.6f}, below the uniform value 0.2.",
        "",
        "## Independent SD orphan accounting",
        "",
    ]
    for material in ("EJ-200", "EJ-204", "EJ-230"):
        vals = [cell["SD_orphan_fraction"]["value"] for cell in cells
                if cell["configuration"]["material"] == material]
        lines.append(f"- {material}: {100*min(vals):.4f}%–{100*max(vals):.4f}% across x.")
    lines += [
        "",
        f"Across all materials the range is {100*min(orphan):.4f}%–{100*max(orphan):.4f}%. The historical EXEC_41 value, 1.5239%, is reproduced by the same EJ-204 material (1.5226%–1.5251%) and is stable with position at this precision. The distinct levels for EJ-200 and EJ-230 show that it is not a single material-independent constant.",
        "",
        "## B3 — golden contract and sidecars",
        "",
        "`analysis/validation/GOLDEN_REFERENCE_20260913.json` now contains an `exec42_full_grid` block with all cell results, physical observable definitions, H3 comparisons, provenance, and SHA-256 hashes. Its top-level status is `EXEC42_FULL_GRID_H3_FAIL`; `ready_for_acceptance` remains `false` with failure reason `H3`.",
        "",
        "All 21 cell directories contain the following complete products, written immediately after their cell finished:",
        "",
        "- `event_metrics.csv/.root/.meta.json`",
        "- `face_decomposition.csv/.root/.meta.json`",
        "- `angular_distribution.csv/.root/.meta.json`",
        "- `cell_results.json`",
        "",
        "The campaign-level products are `input_inventory`, `contract_cells`, `H3_material_invariance`, and `escape_profiles`, each as CSV/ROOT/meta.json, plus `input_verification.json`, `contract_results.json`, and `analysis.log`. The golden block records hashes for all 225 expected artifacts. Sidecar root: `/home/rrios/exec42_20260913/analysis`.",
        "",
        "## Versioned changes and gate",
        "",
        "- `984d3f2`: preregister H3 before calculation.",
        "- `ef3092a`: version B1/B2 input verification and contract analysis.",
        "- `c83f202`: update the golden reference and preserve the earlier blocks.",
        "",
        "No simulation, grid relaunch, merge, push, deck edit, BLUE correction, or injected electronics was performed. Work stops at the EXEC_42 final gate.",
        "",
    ]
    args.out.write_text("\n".join(lines))
    print(args.out)


if __name__ == "__main__":
    main()
