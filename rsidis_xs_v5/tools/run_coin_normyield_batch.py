#!/usr/bin/env python3
"""Run RP_get_coin_normyield for every completed Pass1 good-coin table."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "results" / "Pass1" / "Tables" / "RP_get_good_coin_ev"
MACRO = "macros/RP_get_coin_normyield.C"
EXCLUDED = {
    "RP_get_good_coin_ev_Summary.csv",
    "RP_get_good_coin_ev_Problematic.csv",
}


def main() -> int:
    inputs = sorted(
        path for path in INPUT_DIR.glob("phase*.csv") if path.name not in EXCLUDED
    )
    if not inputs:
        print(f"No per-setting good-coin tables found under {INPUT_DIR}", file=sys.stderr)
        return 1

    failures: list[tuple[str, int]] = []
    for index, path in enumerate(inputs, 1):
        relative_input = path.relative_to(PROJECT_DIR)
        invocation = f'{MACRO}("{relative_input.as_posix()}")'
        result = subprocess.run(
            ["root", "-l", "-b", "-q", invocation],
            cwd=PROJECT_DIR,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        status = "OK" if result.returncode == 0 else "FAILED"
        print(f"[{index:03d}/{len(inputs)}] {path.stem}: {status}", flush=True)
        if result.returncode:
            failures.append((relative_input.as_posix(), result.returncode))

    summary = subprocess.run(
        [sys.executable, "tools/RP_get_coin_normyield_Summary.py"],
        cwd=PROJECT_DIR,
        check=False,
    )
    if summary.returncode:
        return summary.returncode
    if failures:
        print(f"Setting failures: {len(failures)}", file=sys.stderr)
        for path, returncode in failures:
            print(f"  {path}: ROOT exit {returncode}", file=sys.stderr)
        return 1
    print(f"Completed {len(inputs)} settings and rebuilt summary tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
