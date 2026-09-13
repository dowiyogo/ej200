#!/usr/bin/env python3
"""Prepare the EXEC42 production-observation grid without launching it."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil

import detached_grid as launcher


POSITIONS = [0, 200, -200, 500, -500, 650, -650]
MATERIALS = [("EJ200", "EJ-200", "OPSC-100"),
             ("EJ204", "EJ-204", "OPSC-101"),
             ("EJ230", "EJ-230", "OPSC-106")]


def command(raw, name):
    return launcher.macro_value(raw, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--source-campaign", type=Path, required=True)
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--pilot-results", type=Path, required=True)
    parser.add_argument("--pilot-rss-bytes", type=int, required=True)
    parser.add_argument("--handoff", type=Path, required=True)
    parser.add_argument("--scaling", type=Path, required=True)
    parser.add_argument("--reproducibility", type=Path, required=True)
    parser.add_argument("--exec40-prereg", type=Path, required=True)
    parser.add_argument("--exec41-prereg", type=Path, required=True)
    args = parser.parse_args()
    dest = args.directory.resolve()
    launcher.require(not dest.exists(), "Campaign output directory already exists")
    launcher.require(args.binary.is_file() and os.access(args.binary, os.X_OK), "Binary absent/not executable")
    cache = (args.binary.parent/"CMakeCache.txt").read_text()
    launcher.require("EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF" in cache, "Binary is not diagnostics OFF")
    launcher.require((args.binary.parent/"sslg4").is_dir(), "Binary sslg4 resources absent")
    prior_handoff = launcher.read(args.handoff)
    launcher.require(prior_handoff.get("EJ200_OPSC_CODE") == "OPSC-100", "EJ-200 OPSC handoff mismatch")
    pilot = launcher.read(args.pilot_results)
    launcher.require(pilot["G_P_1"]["status"] == "PASS" and pilot["G_P_2"]["status"] == "PASS",
                     "Pilot physics gates did not pass")
    source = launcher.read(args.source_campaign)
    source_cells = {cell["cell_id"]: cell for cell in source["cells"]}
    launcher.require(len(source_cells) == 21, "Source campaign does not contain 21 cells")

    dest.mkdir(parents=True)
    (dest/"grid_logs").mkdir()
    local_handoff = dest/"exec42_handoff.json"
    launcher.write(local_handoff, {
        "schema": "EXEC42_HANDOFF_v1", "created_utc": launcher.now(),
        "EJ200_OPSC_CODE": "OPSC-100", "source_handoff": str(args.handoff),
        "source_handoff_sha256": launcher.sha(args.handoff),
        "pilot_results": str(args.pilot_results), "pilot_results_sha256": launcher.sha(args.pilot_results),
        "G_P_1": pilot["G_P_1"], "G_P_2": pilot["G_P_2"], "G_P_3": pilot["G_P_3"],
    })

    resources = {}
    for _, _, opsc in MATERIALS:
        macro = args.binary.parent/"sslg4"/"macros"/"oscnt"/(opsc.lower()+".mac")
        data = args.binary.parent/"sslg4"/"data"/"oscnt"/opsc.lower()
        launcher.require(macro.is_file() and data.is_dir(), "Missing resource set for "+opsc)
        resources[str(macro)] = launcher.sha(macro)
        for path in sorted(data.iterdir()):
            if path.is_file():
                resources[str(path)] = launcher.sha(path)
    pde = Path(pilot["metadata"]["candidate"]).parent.parent.parent/"ej200_exec40_20260913"/"data"/"sipm"/"AFBR-S4N66P024M_pde.txt"
    # Use the source-worktree path recorded by the validated pilot log.
    matches = re.findall(r"SiPM PDE file\s*:\s*(.*)",
                         (args.pilot_root.parent/"stdout.log").read_text(errors="replace"))
    launcher.require(matches, "Pilot did not record a PDE path")
    pde = Path(matches[-1].strip()).resolve()
    launcher.require(pde.is_file(), "Pilot PDE path is missing")

    cells = []
    for prefix, material, opsc in MATERIALS:
        for x in POSITIONS:
            suffix = "xp0" if x == 0 else ("xp" if x > 0 else "xm") + str(abs(x))
            cell_id = f"{prefix}_{suffix}"
            original = source_cells[cell_id]
            src = Path(original["source_macro"])
            raw = src.read_bytes()
            adjusted = launcher.adjusted_macro(raw)
            expected = {
                "/run/numberOfThreads": "4", "/run/eventModulo": "1",
                "/random/setSeeds": "26092601 8349041", "/det/readout": "EndTop",
                "/det/scintillator": opsc, "/sipm/model": "AFBR-S4N66P024M",
                "/sipm/jitterSigma": "0 ns", "/gun/particle": "mu-",
                "/gun/energy": "1 GeV", "/muon/angle": "0",
                "/muon/gunX": f"{x} mm", "/run/beamOn": "10000",
            }
            for name, value in expected.items():
                launcher.require(command(adjusted, name) == value,
                                 f"{cell_id}: {name} is not {value}")
            cell_dir = dest/"cells"/cell_id
            launcher.require(not cell_dir.exists(), "Cell output directory already exists: "+cell_id)
            cell_dir.mkdir(parents=True)
            (cell_dir/"run.mac").write_bytes(adjusted)
            (cell_dir/"sslg4").symlink_to(args.binary.parent/"sslg4", target_is_directory=True)
            cells.append({
                "cell_id": cell_id, "material": material, "opsc": opsc, "x_mm": x,
                "output": str(cell_dir), "source_macro": str(src),
                "source_macro_sha256": launcher.sha(src),
                "macro_sha256": launcher.sha(cell_dir/"run.mac"),
            })

    evidence_paths = [args.source_campaign, args.handoff, args.scaling, args.reproducibility,
                      args.exec40_prereg, args.exec41_prereg, args.pilot_results]
    four = next(row for row in launcher.read(args.scaling) if row["workers"] == 4)
    launcher.require(four["status"] == "PASS" and four["efficiency"] > 0.95,
                     "Four-worker scaling evidence failed")
    launcher.require(all(row["different_events_from_S1"] == 0
                         for row in launcher.read(args.reproducibility)),
                     "Worker event reproducibility failed")
    memory_budget = 2*args.pilot_rss_bytes + 256*1024**2
    campaign = {
        "schema": "EXEC42 production-observation grid v1",
        "prepared_utc": launcher.now(), "output_directory": str(dest),
        "manifest": str(dest/"manifest.jsonl"),
        "binary": str(args.binary.resolve()), "binary_sha256": launcher.sha(args.binary),
        "simulation_commit": "b35ee84acadef12c93506e0720580c91f901fcbf",
        "pilot_ROOT": str(args.pilot_root.resolve()), "pilot_ROOT_size": args.pilot_root.stat().st_size,
        "pilot_N_generated": 2000,
        "PDE_path": str(pde), "PDE_sha256": launcher.sha(pde), "geant4_version": "11.4.0",
        "evidence_sha256": {str(path.resolve()): launcher.sha(path) for path in evidence_paths},
        "material_resources_sha256": resources,
        "runtime_environment": {key: value for key, value in os.environ.items()
                                if key.startswith("G4") or key == "LD_LIBRARY_PATH"},
        "handoff": str(local_handoff), "EJ200_OPSC_CODE": "OPSC-100", "excluded": [],
        "concurrency": 6, "workers": 4, "eventModulo": 1, "N_generated": 10000,
        "timeout_s": None, "cells": cells,
        "expected_rss_per_process_bytes": args.pilot_rss_bytes,
        "memory_budget_per_process_bytes": memory_budget,
        "parallelism_provenance": {
            "source": "EXEC33 S1/S2/S3/S4 and EXEC42 A2",
            "scaling": str(args.scaling), "reproducibility": str(args.reproducibility),
            "workers_tested": [1, 4, 12, 24], "exact_delta_npe_end": 0,
            "different_events": 0, "selected_workers": 4,
            "efficiency_four_workers": four["efficiency"],
            "A2_max_rss_bytes": args.pilot_rss_bytes,
            "concurrency_selection": "6 retained from validated EXEC34R; measured RAM permits it with a large reserve",
        },
    }
    launcher.write(dest/"campaign.json", campaign)
    print(f"Prepared {len(cells)} cells at {dest}; EXECUTED=0")


if __name__ == "__main__":
    main()
