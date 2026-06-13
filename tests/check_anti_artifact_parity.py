#!/usr/bin/env python3
"""Protect END/SiPM/EJ-204 bytes and timing-estimator parity with the base."""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys


INVARIANT_FILES = (
    "src/SiPMSD.cc",
    "src/PrimaryGeneratorAction.cc",
    "src/EventAction.cc",
    "analysis/congruent_sum4_timing.C",
    "data/sipm/AFBR-S4N66P024M_pde.txt",
    "src/external/SSLG4/macros/oscnt/opsc-101.mac",
    "src/external/SSLG4/data/oscnt/opsc-101/absLength.txt",
    "src/external/SSLG4/data/oscnt/opsc-101/rIndex.txt",
    "src/external/SSLG4/data/oscnt/opsc-101/scntComp1.txt",
)


def git(repo: pathlib.Path, *args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(repo), *args])


def base_ref(repo: pathlib.Path) -> str:
    for candidate in ("feat/endtop-sslg4", "origin/feat/endtop-sslg4"):
        if subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--verify", candidate],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode == 0:
            return candidate
    raise SystemExit("cannot resolve feat/endtop-sslg4 base for parity guard")


def block(text: str, start: str, stop: str) -> str:
    begin = text.index(start)
    end = text.index(stop, begin)
    return text[begin:end]


def constant(text: str, name: str) -> float:
    match = re.search(rf"{re.escape(name)}\s*=\s*([0-9.]+)", text)
    if not match:
        raise SystemExit(f"missing timing constant {name}")
    return float(match.group(1))


def main() -> int:
    repo = pathlib.Path(sys.argv[1]).resolve()
    base = base_ref(repo)
    for path in INVARIANT_FILES:
        current = (repo / path).read_bytes()
        original = git(repo, "show", f"{base}:{path}")
        if current != original:
            raise SystemExit(f"base parity violation: {path}")

    detector = (repo / "src/DetectorConstruction.cc").read_text()
    base_detector = git(repo, "show", f"{base}:src/DetectorConstruction.cc").decode()
    end_block = block(detector, "    if (IsEndInstrumented()) {", "    if (IsTopInstrumented()) {")
    base_end_block = block(
        base_detector, "    if (IsEndInstrumented()) {", "    if (IsTopInstrumented()) {"
    )
    if end_block != base_end_block:
        raise SystemExit("base parity violation: END geometry/placements changed")

    construction = block(
        detector,
        "G4VPhysicalVolume* DetectorConstruction::Construct() {",
        "    // ── World volume",
    )
    base_construction = block(
        base_detector,
        "G4VPhysicalVolume* DetectorConstruction::Construct() {",
        "    // ── World volume",
    )
    if construction != base_construction:
        raise SystemExit("base parity violation: scintillator/SiPM material construction changed")

    materials = (repo / "src/Materials.cc").read_text()
    base_materials = git(repo, "show", f"{base}:src/Materials.cc").decode()
    for label, start, stop in (
        (
            "SiPM coupling",
            "G4Material* CreateSiPMCoupling() {",
            "G4OpticalSurface* CreateBarSurface() {",
        ),
        (
            "SiPM optical surface",
            "G4OpticalSurface* CreateSiPMSurface",
            "G4Material* CreateMylar() {",
        ),
    ):
        if block(materials, start, stop) != block(base_materials, start, stop):
            raise SystemExit(f"base parity violation: {label} changed")

    cpp = (repo / "analysis/congruent_sum4_timing.C").read_text()
    python = (repo / "analysis/endonly_sum4.py").read_text()
    pairs = (
        ("kSprRiseNs", "SPR_RISE_NS"),
        ("kSprFallNs", "SPR_FALL_NS"),
        ("kThresholdPe", "LEADING_EDGE_THRESHOLD_PE"),
    )
    for cpp_name, python_name in pairs:
        if constant(cpp, cpp_name) != constant(python, python_name):
            raise SystemExit(f"timing parity violation: {cpp_name} != {python_name}")
    required_python = (
        "return current + high",
        "end_id // 4",
        "earliest(a, b)",
        "1000.0 * sigma / SQRT2",
        "time_origin=Geant4 event start; input time_ns",
    )
    for snippet in required_python:
        if snippet not in python:
            raise SystemExit(f"timing parity violation: missing '{snippet}'")

    print(f"anti_artifact_parity_guard PASSED against {base}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
