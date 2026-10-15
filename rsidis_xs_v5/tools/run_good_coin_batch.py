#!/usr/bin/env python3
"""Run RP_get_good_coin_ev for every Pass1 setting, then rebuild summaries."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent.parent
MANIFEST = PROJECT_DIR / "bigtable" / "Pass1" / "filter_pass1_sidis_bigtable_summary.csv"
MACRO = "macros/RP_get_good_coin_ev.C"


def main() -> int:
    with MANIFEST.open(newline="", encoding="utf-8-sig") as stream:
        settings = list(csv.DictReader(stream))

    failures: list[tuple[str, int]] = []
    for index, setting in enumerate(settings, 1):
        relative_input = Path("bigtable") / "Pass1" / setting["output_file"]
        invocation = f'{MACRO}("{relative_input.as_posix()}")'
        result = subprocess.run(
            ["root", "-l", "-b", "-q", invocation],
            cwd=PROJECT_DIR,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        status = "OK" if result.returncode == 0 else "FAILED"
        print(
            f"[{index:03d}/{len(settings)}] {relative_input.stem}: {status}",
            flush=True,
        )
        if result.returncode:
            failures.append((relative_input.as_posix(), result.returncode))

    if failures:
        print(f"Setting failures: {len(failures)}", file=sys.stderr)
        for path, returncode in failures:
            print(f"  {path}: ROOT exit {returncode}", file=sys.stderr)
        return 1

    summary = subprocess.run(
        [sys.executable, "tools/RP_get_good_coin_ev_Summary.py"],
        cwd=PROJECT_DIR,
        check=False,
    )
    if summary.returncode:
        return summary.returncode
    print(f"Completed {len(settings)} settings and rebuilt summary tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
