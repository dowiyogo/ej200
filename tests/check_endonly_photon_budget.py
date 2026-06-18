#!/usr/bin/env python3
"""Run a short END-only job and verify that only END detectors collect photons."""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys
import tempfile


def value(log: str, label: str) -> int:
    matches = re.findall(rf"{re.escape(label)}\s*:\s*(\d+)", log)
    if not matches:
        raise SystemExit(f"missing diagnostic '{label}'")
    return int(matches[-1])


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    macro = pathlib.Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="endonly_guard_") as temporary:
        work = pathlib.Path(temporary)
        (work / "sslg4").symlink_to(executable.parent / "sslg4", target_is_directory=True)
        completed = subprocess.run(
            [str(executable), "-t", "4", "-m", str(macro)],
            cwd=work,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
    log = completed.stdout
    active_end = value(log, "Active End SiPMs")
    active_top = value(log, "Active Top SiPMs")
    top_photons = value(log, "Top SiPM  photons")
    generated = value(log, "Scint photons generated")
    detected = value(log, "Total photons detected")
    if active_end != 16 or active_top != 0:
        raise SystemExit(f"wrong active channel map: END={active_end}, TOP={active_top}")
    if top_photons != 0:
        raise SystemExit(f"+Y/TOP detection is non-zero: {top_photons}")
    if generated <= 0 or detected <= 0:
        raise SystemExit(
            f"invalid photon budget: generated={generated}, detected={detected}"
        )
    print(
        "endonly_photon_budget_guard PASSED: "
        f"END={active_end}, TOP={active_top}, top_photons={top_photons}, detected={detected}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
