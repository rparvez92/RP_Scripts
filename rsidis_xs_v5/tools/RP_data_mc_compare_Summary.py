#!/usr/bin/env python3
"""Summarize setting-wise RP data-to-MC comparison tables."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Mapping


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = PROJECT_DIR / "results" / "Tables" / "RP_data_mc_compare"
SUMMARY_NAME = "RP_data_mc_compare_Summary.csv"
PROBLEMATIC_NAME = "RP_data_mc_compare_Problematic.csv"
COLUMNS = [
    "Phase", "Pass", "Run_type", "Target", "Setting", "Status", "Reason",
    "S_dummy", "Binning_version", "Provisional_flag", "N_bin_rows", "N_variables",
    "Max_data_flow_fraction", "Max_MC_flow_fraction",
    "Delta_Data_final", "Delta_Data_final_err", "Delta_MC_total",
    "Delta_MC_total_err", "Delta_Data_by_MC",
    "Full_Data_final", "Full_Data_final_err", "Full_MC_total",
    "Full_MC_total_err", "Full_Data_by_MC",
]


def atomic_write(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def summarize(path: Path) -> Dict[str, object]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty comparison table: {path}")
    first = rows[0]
    output: Dict[str, object] = {
        key: first.get(key, "") for key in
        ("Phase", "Pass", "Run_type", "Target", "Setting", "S_dummy", "Binning_version")
    }
    statuses = [row.get("Status", "") for row in rows]
    precedence = ("ERROR", "WARNING", "PENDING", "SKIPPED", "OK")
    output["Status"] = next((status for status in precedence if status in statuses), "ERROR")
    output["Reason"] = ";".join(sorted({
        reason for row in rows for reason in row.get("Reason", "").split(";") if reason
    }))
    physics_rows = [row for row in rows if row.get("Variable")]
    output["N_bin_rows"] = len(physics_rows)
    output["N_variables"] = len({row["Variable"] for row in physics_rows})
    output["Provisional_flag"] = first.get("Provisional_flag", "")
    data_flows = [
        number(row.get("Data_flow_fraction", "")) for row in physics_rows
    ]
    mc_flows = [
        number(row.get("MC_flow_fraction", "")) for row in physics_rows
    ]
    output["Max_data_flow_fraction"] = max(
        (value for value in data_flows if math.isfinite(value)), default=math.nan
    )
    output["Max_MC_flow_fraction"] = max(
        (value for value in mc_flows if math.isfinite(value)), default=math.nan
    )
    for tier, prefix in (("delta", "Delta"), ("full", "Full")):
        # One variable is enough for integrated normalization; hdelta includes flow.
        selected = [
            row for row in physics_rows
            if row.get("Tier") == tier and row.get("Variable") == "hdelta"
        ]
        complete = bool(selected) and all(
            row.get("MC_complete") == "1" for row in selected
        )
        for source in ("Data_final", "MC_total"):
            if source == "MC_total" and not complete:
                output[f"{prefix}_{source}"] = ""
                output[f"{prefix}_{source}_err"] = ""
                continue
            values = [number(row[source]) for row in selected]
            errors = [number(row[source + "_err"]) for row in selected]
            total = sum(value for value in values if math.isfinite(value))
            variance = sum(error * error for error in errors if math.isfinite(error))
            output[f"{prefix}_{source}"] = total
            output[f"{prefix}_{source}_err"] = math.sqrt(variance)
        data, mc = output[f"{prefix}_Data_final"], output[f"{prefix}_MC_total"]
        output[f"{prefix}_Data_by_MC"] = (
            data / mc if isinstance(mc, (float, int)) and mc else math.nan
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    directory = args.input_dir.expanduser().resolve()
    paths = sorted(
        path for path in directory.glob("phase*_pass*_PI*.csv")
        if path.name not in {SUMMARY_NAME, PROBLEMATIC_NAME}
    )
    rows = [summarize(path) for path in paths]
    problematic = [
        row for row in rows if row["Status"] in {"WARNING", "ERROR"}
    ]
    atomic_write(directory / SUMMARY_NAME, rows)
    atomic_write(directory / PROBLEMATIC_NAME, problematic)
    counts = Counter(str(row["Status"]) for row in rows)
    print(f"Settings: {len(rows)}; " + ", ".join(
        f"{status}={count}" for status, count in sorted(counts.items())
    ))
    print(f"Summary: {directory / SUMMARY_NAME}")
    print(f"Problematic: {directory / PROBLEMATIC_NAME}")
    return 2 if counts["ERROR"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
