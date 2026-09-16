#!/usr/bin/env python3
"""Prepare, but never launch, the 21-cell EXEC_46 production campaign."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil


SOURCE_GRID = Path("/home/rrios/exec42_20260913/grid")
DEFAULT_OUTPUT = Path("/home/rrios/exec46_20260915/full_grid")
DEFAULT_BINARY = Path("/home/rrios/exec46_20260915/build_baseline/ej200_bar_sim")
PILOT_ROOT = Path("/home/rrios/exec46_20260915/final_x0/photon_hits_run000.root")
VALIDATION_METRICS = Path("/home/rrios/exec46_20260915/validation/validation_metrics.json")
PDE_PATH = Path("/home/rrios/ej200/data/sipm/AFBR-S4N66P024M_pde.txt")
SIMULATION_COMMIT = "4967ec8"
EJ200_OPSC_CODE = "OPSC-100"
CONCURRENCY = 6
WORKERS = 4
EVENT_MODULO = 1
EVENTS = 10000
PILOT_EVENTS = 500
MEASURED_RSS_KIB = 167100


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def command_value(text, command):
    values = re.findall(r"^\s*" + re.escape(command) + r"\s+([^#\r\n]+)", text, re.M)
    require(len(values) == 1, f"missing or duplicate {command}")
    return values[0].strip()


def write_new(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def prepare(output, binary, ej200_sslg4_source=None):
    output = output.resolve()
    binary = binary.resolve()
    require(not output.exists(), f"refusing to overwrite {output}")
    require(binary.is_file() and os.access(binary, os.X_OK), "binary is not executable")
    require("EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF" in
            (binary.parent / "CMakeCache.txt").read_text(), "diagnostics are not OFF")
    require((binary.parent / "sslg4").is_dir(), "missing SSLG4 runtime directory")
    if ej200_sslg4_source is not None:
        ej200_sslg4_source = ej200_sslg4_source.resolve()
        require(ej200_sslg4_source.is_dir(), "alternate EJ-200 SSLG4 directory is missing")
        for relative in ("macros/oscnt/opsc-100.mac", "data/oscnt/opsc-100/rIndex.txt",
                         "data/oscnt/opsc-100/absLength.txt"):
            require((ej200_sslg4_source / relative).is_file(),
                    f"alternate EJ-200 SSLG4 is incomplete: {relative}")
    require(PILOT_ROOT.is_file() and VALIDATION_METRICS.is_file(), "missing validation evidence")

    source_campaign = json.loads((SOURCE_GRID / "campaign.json").read_text())
    source_cells = sorted(source_campaign["cells"], key=lambda cell: cell["cell_id"])
    require(len(source_cells) == 21, "source campaign is not a 21-cell grid")
    require({cell["opsc"] for cell in source_cells if cell["material"] == "EJ-200"} ==
            {EJ200_OPSC_CODE}, "EJ-200 OPSC code mismatch")

    output.mkdir(parents=True)
    (output / "cells").mkdir()
    (output / "grid_logs").mkdir()
    cells = []
    for source_cell in source_cells:
        source_macro = SOURCE_GRID / "cells" / source_cell["cell_id"] / "run.mac"
        text = source_macro.read_text()
        require(command_value(text, "/run/numberOfThreads") == str(WORKERS),
                f"workers changed in {source_cell['cell_id']}")
        require(command_value(text, "/run/eventModulo") == str(EVENT_MODULO),
                f"eventModulo changed in {source_cell['cell_id']}")
        require(command_value(text, "/run/beamOn") == str(EVENTS),
                f"event count changed in {source_cell['cell_id']}")
        require(command_value(text, "/det/readout") == "EndTop",
                f"readout changed in {source_cell['cell_id']}")
        require(command_value(text, "/sipm/jitterSigma") == "0 ns",
                f"jitter changed in {source_cell['cell_id']}")
        require(command_value(text, "/det/scintillator") == source_cell["opsc"],
                f"material changed in {source_cell['cell_id']}")
        target = output / "cells" / source_cell["cell_id"]
        target.mkdir()
        shutil.copyfile(source_macro, target / "run.mac")
        runtime_sslg4 = (ej200_sslg4_source if source_cell["material"] == "EJ-200"
                         and ej200_sslg4_source is not None else binary.parent / "sslg4")
        (target / "sslg4").symlink_to(runtime_sslg4, target_is_directory=True)
        cells.append({
            "cell_id": source_cell["cell_id"], "material": source_cell["material"],
            "opsc": source_cell["opsc"], "x_mm": source_cell["x_mm"],
            "output": str(target), "source_macro": str(source_macro),
            "source_macro_sha256": sha256(source_macro),
            "macro_sha256": sha256(target / "run.mac"),
            "sslg4_runtime": str(runtime_sslg4),
            "rindex_sha256": sha256(runtime_sslg4 / "data" / "oscnt"
                                    / source_cell["opsc"].lower() / "rIndex.txt"),
            "absLength_sha256": sha256(runtime_sslg4 / "data" / "oscnt"
                                       / source_cell["opsc"].lower() / "absLength.txt"),
        })

    handoff_path = output / "exec46_handoff.json"
    handoff = {
        "campaign": "EXEC_46", "EJ200_OPSC_CODE": EJ200_OPSC_CODE,
        "simulation_commit": SIMULATION_COMMIT, "binary": str(binary),
        "binary_sha256": sha256(binary), "validation_metrics": str(VALIDATION_METRICS),
        "validation_metrics_sha256": sha256(VALIDATION_METRICS),
    }
    write_new(handoff_path, handoff)
    evidence = {
        str(SOURCE_GRID / "campaign.json"): sha256(SOURCE_GRID / "campaign.json"),
        str(VALIDATION_METRICS): sha256(VALIDATION_METRICS),
        str(handoff_path): sha256(handoff_path),
    }
    expected_rss = MEASURED_RSS_KIB * 1024
    campaign = {
        "schema": "EXEC46 track-mechanism production v1",
        "output_directory": str(output), "manifest": str(output / "manifest.jsonl"),
        "binary": str(binary), "binary_sha256": sha256(binary),
        "simulation_commit": SIMULATION_COMMIT,
        "pilot_ROOT": str(PILOT_ROOT), "pilot_ROOT_size": PILOT_ROOT.stat().st_size,
        "pilot_N_generated": PILOT_EVENTS,
        "PDE_path": str(PDE_PATH), "PDE_sha256": sha256(PDE_PATH),
        "geant4_version": "11.4.0", "evidence_sha256": evidence,
        "runtime_environment": {
            "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
            "EJ200_DATA_DIR": str(PDE_PATH.parent.parent),
        },
        "handoff": str(handoff_path), "EJ200_OPSC_CODE": EJ200_OPSC_CODE,
        "EJ200_SSLG4_source": (str(ej200_sslg4_source)
                                if ej200_sslg4_source is not None else None),
        "excluded": [], "concurrency": CONCURRENCY, "workers": WORKERS,
        "eventModulo": EVENT_MODULO, "N_generated": EVENTS, "timeout_s": None,
        "cells": cells, "expected_rss_per_process_bytes": expected_rss,
        "memory_budget_per_process_bytes": 2 * expected_rss + 256 * 1024 ** 2,
        "parallelism_provenance": {
            "source": "EXEC33 reproducibility plus EXEC46 x=0 benchmark",
            "workers_tested": [1, 4, 12, 24], "exact_delta_npe_end": 0,
            "different_events": 0, "selected_workers": WORKERS,
            "EXEC46_uninstrumented_wall_s": 112.65,
            "EXEC46_instrumented_wall_s": 117.68,
            "EXEC46_slowdown": 117.68 / 112.65,
            "scope": "Parallel physics equivalence from EXEC33; instrumentation cost from EXEC46.",
        },
    }
    write_new(output / "campaign.json", campaign)
    print(f"Prepared {len(cells)} cells at {output}; EXECUTED=0")
    print("Launch manually:")
    print("python3 /home/rrios/ej200/analysis/sigma_t/orchestration/detached_grid.py "
          f"launch --directory {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--ej200-sslg4-source", type=Path,
                        help="Alternate complete SSLG4 runtime used only by EJ-200 cells")
    args = parser.parse_args()
    prepare(args.output, args.binary, args.ej200_sslg4_source)


if __name__ == "__main__":
    main()
