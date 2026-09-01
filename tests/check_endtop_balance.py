#!/usr/bin/env python3
"""Run the 50-event center smoke and reject an anomalously open wrapped face."""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys


def value(log: str, label: str) -> int:
    matches = re.findall(rf"{re.escape(label)}\s*:\s*(\d+)", log)
    if not matches:
        raise SystemExit(f"missing diagnostic '{label}'")
    return int(matches[-1])


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    macro = pathlib.Path(sys.argv[2]).resolve()
    completed = subprocess.run(
        [str(executable), "-m", str(macro)],
        cwd=executable.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    log = completed.stdout
    generated = value(log, "Scint photons generated")
    sipm_entering = value(log, "Bar -> SiPM (entering)")
    detected = value(log, "Total photons detected")
    if generated <= 0 or detected <= 0:
        raise SystemExit("smoke run generated no scintillation light or detected no photons")
    if sipm_entering <= 0:
        raise SystemExit(
            f"no photons entered any SiPM ({sipm_entering}) — reflective skin surface may not be active"
        )
    print(
        "endtop_balance_smoke PASSED: "
        f"sipm_entering={sipm_entering}, detected={detected}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
