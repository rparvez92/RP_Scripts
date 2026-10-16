#!/usr/bin/env python3
"""Combine Pass1 setting-level SIDIS extraction rows into summary reports."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import tempfile


def atomic_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        default=root / "results/Pass1/Tables/RP_extract_sidis_xsec")
    args = parser.parse_args()
    sources = sorted(path for path in args.input_dir.glob("phase*.csv"))
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for source in sources:
        with source.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if not fields: fields = reader.fieldnames or []
            rows.extend(reader)
    order = ("Phase", "Pass", "Run_type", "Target", "setting_id")
    rows.sort(key=lambda row: tuple(row.get(field, "") for field in order))
    problematic = [row for row in rows if row.get("Extraction_status") != "OK"]
    summary = args.input_dir / "RP_extract_sidis_xsec_Summary.csv"
    problem = args.input_dir / "RP_extract_sidis_xsec_Problematic.csv"
    atomic_csv(summary, fields, rows); atomic_csv(problem, fields, problematic)
    print(f"Wrote {len(rows)} rows to {summary}")
    print(f"Wrote {len(problematic)} rows to {problem}")
    return 1 if any(row.get("Extraction_status") == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
