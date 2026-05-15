#!/usr/bin/env python3
"""Compatibility wrapper for the current liftover_by_id implementation.

This keeps the Nextflow interface stable now. A later cleanup pass can vendor the
full script here if nf-liftover must become completely independent of tc-pytools.
"""

from pathlib import Path
import runpy


SOURCE = Path("/public/scripts/tc-pytools/liftover/liftover_by_id.py")


if __name__ == "__main__":
    if not SOURCE.is_file():
        raise SystemExit(f"Missing liftover implementation: {SOURCE}")
    runpy.run_path(str(SOURCE), run_name="__main__")
