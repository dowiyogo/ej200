#!/usr/bin/env python3
"""Assert the intrinsic < total resolution guard and its negative control."""

from __future__ import annotations

import ast
import pathlib
import sys


def main() -> int:
    analysis_path = pathlib.Path(sys.argv[1]).resolve()
    tree = ast.parse(analysis_path.read_text())
    selected = [
        node
        for node in tree.body
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "DEFAULT_READOUT_JITTER_PS"
                for target in node.targets
            )
        )
        or (
            isinstance(node, ast.FunctionDef)
            and node.name in {"modeled_total_sigma_ps", "require_physical_ordering"}
        )
    ]
    namespace: dict[str, object] = {}
    module = ast.fix_missing_locations(
        ast.Module([ast.Import(names=[ast.alias(name="math")]), *selected], [])
    )
    exec(compile(module, str(analysis_path), "exec"), namespace)

    for intrinsic in (1.0, 20.0, 88.0, 200.0):
        total = namespace["modeled_total_sigma_ps"](intrinsic)
        namespace["require_physical_ordering"](intrinsic, total)
    try:
        namespace["require_physical_ordering"](100.0, 100.0)
    except ValueError:
        pass
    else:
        raise SystemExit("physical-ordering guard accepted total <= intrinsic")
    print(
        "physical_ordering_guard PASSED: modeled total includes positive "
        f"{namespace['DEFAULT_READOUT_JITTER_PS']:.1f} ps readout jitter"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
