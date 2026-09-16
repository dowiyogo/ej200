#!/usr/bin/env python3
"""Regression test for per-cell SSLG4 source targets."""

import importlib.util
import shutil
import tempfile
from pathlib import Path


SCRIPT = Path(__file__).with_name("prepare_campaign.py")
spec = importlib.util.spec_from_file_location("prepare_campaign", SCRIPT)
prepare_campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_campaign)


def main():
    with tempfile.TemporaryDirectory(prefix="exec46_prepare_test_") as directory:
        output = Path(directory) / "campaign"
        prepare_campaign.prepare(
            output,
            prepare_campaign.DEFAULT_BINARY,
            prepare_campaign.I2_BC404_SSLG4,
        )
        targets = {
            cell_id: (output / "cells" / cell_id / "sslg4").readlink()
            for cell_id in ("EJ200_xm200", "EJ204_xm200")
        }
        for cell_id, target in targets.items():
            assert target.name == cell_id, f"{cell_id}: unexpected target {target}"
        shutil.rmtree(output)


if __name__ == "__main__":
    main()