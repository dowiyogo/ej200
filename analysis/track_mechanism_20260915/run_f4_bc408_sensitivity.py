#!/usr/bin/env python3
"""Prepare and run the two authorized EXEC_46 F4 BC-408 sensitivity cells."""

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from datetime import datetime, timezone

import numpy as np
import uproot


BASE = Path("/home/rrios/exec46_20260915")
OUTPUT = BASE / "f4_bc408_sensitivity"
BINARY = BASE / "build_baseline" / "ej200_bar_sim"
SOURCE_SSLG4 = BASE / "build_baseline" / "sslg4"
SOURCE_MACRO = BASE / "full_grid" / "cells" / "EJ200_xm650" / "run.mac"
EXPECTED_EVENTS = 10_000
EXPECTED_WORKERS = 4
EXPECTED_SEEDS = "26092601 8349041"
FIT_A = 1.518
FIT_B = 0.640
FIT_C_PER_NM = 0.00423
RINDEX_MEASURED_MIN_NM = 370.0
RINDEX_MEASURED_MAX_NM = 660.0
RINDEX_STEP_NM = 2.0
UV_ABS_WAVELENGTH_NM = 372.0
UV_ABS_LENGTH_MM = 28.23
VISIBLE_ANCHOR_NM = 439.0
VARIANTS = {
    "visible_lower_764mm": 764.0,
    "visible_current_3800mm": 3800.0,
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def n_bc408(wavelength_nm):
    return FIT_A + FIT_B * math.exp(-FIT_C_PER_NM * wavelength_nm)


def rindex_rows():
    measured = np.arange(RINDEX_MEASURED_MIN_NM,
                         RINDEX_MEASURED_MAX_NM + 0.5 * RINDEX_STEP_NM,
                         RINDEX_STEP_NM)
    rows = [(200.0, n_bc408(RINDEX_MEASURED_MIN_NM))]
    rows.extend((float(wavelength), n_bc408(float(wavelength)))
                for wavelength in measured)
    rows.append((800.0, n_bc408(RINDEX_MEASURED_MAX_NM)))
    return rows


def write_table(path, rows, value_scale=1.0):
    with path.open("w") as stream:
        for wavelength, value in rows:
            stream.write(f"{wavelength:.2f} {value / value_scale:.12g}\n")


def prepare_variant(name, visible_abs_length_mm):
    directory = OUTPUT / name
    if directory.exists():
        meta_path = directory / "configuration.meta.json"
        require(meta_path.is_file(), f"partial directory without metadata: {directory}")
        metadata = json.loads(meta_path.read_text())
        macro_path = directory / "run.mac"
        mpt = directory / "sslg4" / "data" / "oscnt" / "opsc-100"
        require(metadata.get("variant") == name, f"variant mismatch: {directory}")
        require(metadata.get("absorption", {}).get("visible_scenario_mm")
                == visible_abs_length_mm, f"absorption mismatch: {directory}")
        require(metadata.get("binary_sha256") == sha256(BINARY),
                f"binary changed since preparation: {directory}")
        require(metadata.get("source_macro_sha256") == sha256(SOURCE_MACRO),
                f"source macro changed since preparation: {directory}")
        require(macro_path.is_file() and sha256(macro_path) == sha256(SOURCE_MACRO),
                f"prepared macro differs from source: {directory}")
        require((mpt / "rIndex.txt").is_file()
                and metadata.get("rindex", {}).get("sha256")
                == sha256(mpt / "rIndex.txt"), f"RINDEX changed: {directory}")
        require((mpt / "absLength.txt").is_file()
                and metadata.get("absorption", {}).get("sha256")
                == sha256(mpt / "absLength.txt"), f"ABSLENGTH changed: {directory}")
        require(np.allclose(np.loadtxt(mpt / "rIndex.txt"),
                            np.asarray(rindex_rows()), rtol=0.0, atol=5e-12),
                f"RINDEX no longer matches the declared model: {directory}")
        expected_absorption_cm = np.asarray([
            (200.0, UV_ABS_LENGTH_MM / 10.0),
            (UV_ABS_WAVELENGTH_NM, UV_ABS_LENGTH_MM / 10.0),
            (VISIBLE_ANCHOR_NM, visible_abs_length_mm / 10.0),
            (800.0, visible_abs_length_mm / 10.0),
        ])
        require(np.allclose(np.loadtxt(mpt / "absLength.txt"),
                            expected_absorption_cm, rtol=0.0, atol=5e-12),
                f"ABSLENGTH no longer matches the declared scenario: {directory}")
        return directory
    directory.mkdir(parents=True)
    shutil.copyfile(SOURCE_MACRO, directory / "run.mac")
    shutil.copytree(SOURCE_SSLG4, directory / "sslg4")
    macro = (directory / "run.mac").read_text()
    require(macro.count("/run/numberOfThreads 4") == 1, "workers changed")
    require(macro.count("/run/eventModulo 1") == 1, "eventModulo changed")
    require(macro.count(f"/random/setSeeds {EXPECTED_SEEDS}") == 1, "seeds changed")
    require(macro.count("/run/beamOn 10000") == 1, "event count changed")
    require(macro.count("/muon/gunX -650 mm") == 1, "gun position changed")
    require(macro.count("/det/scintillator OPSC-100") == 1, "material changed")
    mpt = directory / "sslg4" / "data" / "oscnt" / "opsc-100"
    write_table(mpt / "rIndex.txt", rindex_rows())
    # Values are written in cm because opsc-100.mac declares the table unit as cm.
    # The 372--439 nm segment is Geant4's interpolation between the measurement
    # and the visible scenario. Below 372 nm the 372 nm measurement is held
    # constant; this is an explicit hypothesis, not a UV extrapolation.
    absorption = [
        (200.0, UV_ABS_LENGTH_MM),
        (UV_ABS_WAVELENGTH_NM, UV_ABS_LENGTH_MM),
        (VISIBLE_ANCHOR_NM, visible_abs_length_mm),
        (800.0, visible_abs_length_mm),
    ]
    write_table(mpt / "absLength.txt", absorption, value_scale=10.0)
    payload = {
        "variant": name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "binary": str(BINARY),
        "binary_sha256": sha256(BINARY),
        "source_macro": str(SOURCE_MACRO),
        "source_macro_sha256": sha256(SOURCE_MACRO),
        "macro_sha256": sha256(directory / "run.mac"),
        "events": EXPECTED_EVENTS,
        "workers": EXPECTED_WORKERS,
        "seeds": [26092601, 8349041],
        "cell": "EJ200_xm650",
        "x_mm": -650,
        "opsc": "OPSC-100",
        "rindex": {
            "published_model": "n(lambda)=1.518+0.640*exp(-0.00423*lambda_nm)",
            "sample_step_nm": RINDEX_STEP_NM,
            "measured_domain_nm": [370.0, 660.0],
            "below_domain": "constant n(370 nm) from 200 to 370 nm; explicit hypothesis",
            "above_domain": "constant n(660 nm) from 660 to 800 nm",
            "sha256": sha256(mpt / "rIndex.txt"),
        },
        "absorption": {
            "uv_measurement": "28.23 +/- 2.7 mm at 372 nm",
            "visible_scenario_mm": visible_abs_length_mm,
            "visible_anchor_nm": 439.0,
            "interpretation": ("90% lower-bound scenario, not a central value"
                               if visible_abs_length_mm == 764.0
                               else "pre-existing EJ-200 visible value"),
            "below_372_nm": "constant 28.23 mm; explicit unmeasured-UV hypothesis",
            "between_372_and_439_nm": "Geant4 material-table interpolation",
            "at_and_above_439_nm": "constant visible scenario through 800 nm",
            "sha256": sha256(mpt / "absLength.txt"),
        },
    }
    (directory / "configuration.meta.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return directory


def run_variant(directory):
    root_path = directory / "photon_hits_run000.root"
    done_path = directory / ".DONE.json"
    if done_path.is_file():
        done = json.loads(done_path.read_text())
        require(root_path.is_file() and sha256(root_path) == done["root_sha256"],
                f"invalid existing result: {directory}")
        return done
    command = [str(BINARY), "-m", str(directory / "run.mac")]
    environment = os.environ.copy()
    environment.update({
        "EJ200_DATA_DIR": "/home/rrios/ej200/data",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    start = datetime.now(timezone.utc)
    monotonic_start = time.monotonic()
    with (directory / "simulation.log").open("wb") as log:
        process = subprocess.run(command, cwd=directory, env=environment,
                                 stdout=log, stderr=subprocess.STDOUT, check=False)
    wall = time.monotonic() - monotonic_start
    require(process.returncode == 0, f"simulation failed: {directory}")
    require(root_path.is_file(), f"ROOT missing: {directory}")
    with uproot.open(root_path) as root_file:
        tree = root_file["sipm_hits"]
        entries = int(tree.num_entries)
        seen = np.zeros(EXPECTED_EVENTS, dtype=bool)
        for batch in tree.iterate(["event_id"], step_size="256 MB", library="np"):
            events = batch["event_id"].astype(np.int64, copy=False)
            require(np.all((events >= 0) & (events < EXPECTED_EVENTS)), "bad event id")
            seen[events] = True
    require(np.all(seen), f"not all {EXPECTED_EVENTS} events represented")
    done = {
        "status": "COMPLETE",
        "start_utc": start.isoformat(),
        "end_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "cwd": str(directory),
        "exit_code": process.returncode,
        "wall_s": wall,
        "events_with_hits": int(seen.sum()),
        "root_entries": entries,
        "root_path": str(root_path),
        "root_size_bytes": root_path.stat().st_size,
        "root_sha256": sha256(root_path),
    }
    done_path.write_text(json.dumps(done, indent=2, sort_keys=True) + "\n")
    return done


def main():
    require(BINARY.is_file() and os.access(BINARY, os.X_OK), "binary unavailable")
    cache = (BINARY.parent / "CMakeCache.txt").read_text()
    require("EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF" in cache, "diagnostics binary required OFF")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    directories = [prepare_variant(name, value) for name, value in VARIANTS.items()]
    # The two independent 4-worker cells run concurrently. No timeout is used.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(run_variant, directories))
    campaign = {
        "status": "COMPLETE",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parallel_processes": 2,
        "workers_per_process": EXPECTED_WORKERS,
        "timeout_s": None,
        "results": results,
    }
    (OUTPUT / "campaign.meta.json").write_text(
        json.dumps(campaign, indent=2, sort_keys=True) + "\n")
    print(json.dumps(campaign, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
