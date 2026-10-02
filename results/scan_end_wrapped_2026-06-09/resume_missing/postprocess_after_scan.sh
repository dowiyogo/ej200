#!/usr/bin/env bash
set -euo pipefail

repo=/mnt/d/SHiP/ej200_edge_scan
scan_dir="$repo/results/scan_end_wrapped_2026-06-09"
resume_dir="$scan_dir/resume_missing"
sim_pid=9406

while kill -0 "$sim_pid" 2>/dev/null; do
    sleep 60
done

python3 - "$resume_dir" "$scan_dir" <<'PY'
import csv
import shutil
import sys
from pathlib import Path

import numpy as np
import uproot

resume_dir = Path(sys.argv[1])
scan_dir = Path(sys.argv[2])
positions = [-150, -100, -50, 0, 50, 100, 150, 200, 250, 300, 350, 400,
             450, 500, 550, 600, 650, 670, 690]
rows = []

for local_run, expected_x in enumerate(positions):
    source = resume_dir / f"photon_hits_run{local_run:03d}.root"
    target = scan_dir / f"photon_hits_run{local_run + 12:03d}.root"
    with uproot.open(source) as root_file:
        tree = root_file["sipm_hits"]
        arrays = tree.arrays(["event_id", "gun_x_mm"], library="np")
    event_ids = np.unique(arrays["event_id"])
    x_values = np.unique(arrays["gun_x_mm"])
    complete = (
        len(event_ids) == 10000
        and event_ids[0] == 0
        and event_ids[-1] == 9999
        and x_values.tolist() == [float(expected_x)]
    )
    rows.append((local_run, local_run + 12, expected_x, tree.num_entries,
                 len(event_ids), complete))
    if not complete:
        raise RuntimeError(f"Incomplete or mismatched ROOT: {source}")
    if target.exists():
        raise RuntimeError(f"Refusing to overwrite existing target: {target}")

for local_run, _, _, _, _, _ in rows:
    source = resume_dir / f"photon_hits_run{local_run:03d}.root"
    target = scan_dir / f"photon_hits_run{local_run + 12:03d}.root"
    shutil.move(source, target)

with (scan_dir / "resume_validation.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["resume_run", "final_run", "x_mm", "entries", "events", "complete"])
    writer.writerows(rows)
PY

mkdir -p "$scan_dir/sum8_analysis"
/opt/hep/root/bin/root -l -b -q \
    "/tmp/endred_photon_budget.C(\"$scan_dir\",\"$scan_dir/sum8_analysis\",10000,8)" \
    > "$scan_dir/sum8_analysis.log" 2>&1

ctest --test-dir "$repo/build" --output-on-failure \
    > "$scan_dir/ctest_full.log" 2>&1

touch "$scan_dir/resume_complete.ok"
