#!/usr/bin/env python3
"""Combine setting-level SIDIS extraction rows into Summary reports."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import tempfile


def atomic_write(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main() -> int:
    analysis_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        default=analysis_root / "results/Tables/RP_extract_sidis_xsec")
    args = parser.parse_args()
    sources = sorted(path for path in args.input_dir.glob("phase*.csv")
                     if "Summary" not in path.name and "Problematic" not in path.name)
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for source in sources:
        with source.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if not fields:
                fields = reader.fieldnames or []
            rows.extend(reader)
    rows.sort(key=lambda row: tuple(row.get(key, "") for key in
                                   ("Phase", "Pass", "Run_type", "Target", "Setting")))
    problematic = [row for row in rows if row.get("Extraction_status") in {"WARNING", "ERROR"}]
    summary = args.input_dir / "RP_extract_sidis_xsec_Summary.csv"
    problem = args.input_dir / "RP_extract_sidis_xsec_Problematic.csv"
    atomic_write(summary, fields, rows)
    atomic_write(problem, fields, problematic)
    print(f"Wrote {len(rows)} rows to {summary}")
    print(f"Wrote {len(problematic)} rows to {problem}")
    return 1 if any(row.get("Extraction_status") == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
