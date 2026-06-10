#!/usr/bin/env python3
"""Verify the 70 Top placements from an exported EndTop GDML file."""

from __future__ import annotations

import math
import pathlib
import sys
import xml.etree.ElementTree as ET


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def expected_positions() -> list[float]:
    return [-692.0 + 20.0 * i for i in range(35)] + [
        12.0 + 20.0 * i for i in range(35)
    ]


def main() -> int:
    path = pathlib.Path(sys.argv[1])
    root = ET.parse(path).getroot()
    positions_by_name: dict[str, tuple[float, float, float]] = {}
    for element in root.iter():
        if local_name(element.tag) != "position":
            continue
        name = element.attrib.get("name")
        if name:
            positions_by_name[name] = tuple(
                float(element.attrib.get(axis, "0")) for axis in ("x", "y", "z")
            )

    top_x: list[float] = []
    for physvol in root.iter():
        if local_name(physvol.tag) != "physvol":
            continue
        if not physvol.attrib.get("name", "").startswith("TopSiPMPV"):
            continue

        position = next(
            (child for child in physvol if local_name(child.tag) == "position"), None
        )
        if position is not None:
            xyz = tuple(float(position.attrib.get(axis, "0")) for axis in ("x", "y", "z"))
        else:
            ref = next(
                (child for child in physvol if local_name(child.tag) == "positionref"), None
            )
            if ref is None:
                raise SystemExit("TopSiPMPV has no position or positionref")
            xyz = positions_by_name[ref.attrib["ref"]]
        if abs(xyz[2]) > 0.01:
            raise SystemExit(f"Top SiPM z is not zero: {xyz[2]}")
        top_x.append(xyz[0])

    expected = expected_positions()
    if len(top_x) != 70:
        raise SystemExit(f"expected 70 Top SiPMs in GDML, found {len(top_x)}")
    for actual, target in zip(sorted(top_x), expected):
        if not math.isclose(actual, target, abs_tol=0.01):
            raise SystemExit(f"Top position mismatch: actual={actual}, expected={target}")
    if not math.isclose(expected[35] - expected[34], 24.0, abs_tol=0.01):
        raise SystemExit("central Top pair is not separated by 24 mm")

    print("endtop_gdml_check PASSED: 70 Top SiPMs at exact EXEC_07 positions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
