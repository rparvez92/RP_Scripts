#!/usr/bin/env python3
"""Combine per-setting RP coin-normalized-yield tables and flag bad runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = PROJECT_DIR / "results" / "Tables" / "RP_get_coin_normyield"
SUMMARY_NAME = "RP_get_coin_normyield_Summary.csv"
PROBLEMATIC_NAME = "RP_get_coin_normyield_Problematic.csv"

REQUIRED_COLUMNS = {
    "Run",
    "Run_type",
    "RP_Goodcoin_delta",
    "RP_Goodcoin_err_delta",
    "RP_Goodcoin_full",
    "RP_Goodcoin_err_full",
    "Norm_factor",
    "RP_Normyield_delta",
    "RP_Normyield_err_delta",
    "RP_Normyield_full",
    "RP_Normyield_err_full",
    "RP_Normyield_full_by_delta",
    "normyield",
    "normyield_err",
    "normyield_by_RP_Normyield_delta",
    "Norm_status",
    "Norm_reason",
}

FLAG_COLUMNS = (
    "Normalization_flag",
    "Normyield_ratio_flag",
    "Yield_flag",
    "Numeric_flag",
    "Problematic",
    "Problem_reason",
)


def finite(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def tf(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def flag_row(row: Dict[str, str]) -> Dict[str, str]:
    status = row.get("Norm_status", "").strip()
    normyield_delta = finite(row.get("RP_Normyield_delta", ""))
    normyield_err_delta = finite(row.get("RP_Normyield_err_delta", ""))
    normyield_full = finite(row.get("RP_Normyield_full", ""))
    normyield_err_full = finite(row.get("RP_Normyield_err_full", ""))
    normfactor = finite(row.get("Norm_factor", ""))
    replay_normyield = finite(row.get("normyield", ""))
    ratio = finite(row.get("normyield_by_RP_Normyield_delta", ""))

    normalization_flag = status != "OK"
    yield_flag = (
        normyield_delta is None
        or normyield_delta <= 0.0
        or normyield_full is None
        or normyield_full <= 0.0
    )
    numeric_flag = (
        normfactor is None
        or normfactor <= 0.0
        or normyield_err_delta is None
        or normyield_err_delta < 0.0
        or normyield_err_full is None
        or normyield_err_full < 0.0
    )

    # A missing legacy replay result is not a failure of the new normalization.
    # Apply the comparison threshold only when the replay normyield is usable.
    comparison_available = replay_normyield is not None and replay_normyield != -999.0
    ratio_flag = comparison_available and (
        ratio is None or abs(ratio - 1.0) > 0.20
    )

    reasons: List[str] = []
    if normalization_flag:
        reasons.append("NORMALIZATION_FAILED")
    if ratio_flag:
        reasons.append("NORMYIELD_RATIO_DEVIATION_GT_20PCT")
    if yield_flag:
        reasons.append("RP_NORMYIELD_NONPOSITIVE_OR_INVALID")
    if numeric_flag:
        reasons.append("INVALID_NORMALIZATION_NUMERIC_VALUE")

    row.update(
        {
            "Normalization_flag": tf(normalization_flag),
            "Normyield_ratio_flag": tf(ratio_flag),
            "Yield_flag": tf(yield_flag),
            "Numeric_flag": tf(numeric_flag),
            "Problematic": tf(bool(reasons)),
            "Problem_reason": ";".join(reasons),
        }
    )
    return row


def run_key(row: Dict[str, str]) -> tuple[int, str]:
    raw = row.get("Run", "")
    try:
        return int(float(raw)), raw
    except ValueError:
        return 2**63 - 1, raw


def read_tables(paths: Sequence[Path]) -> tuple[List[str], List[Dict[str, str]]]:
    header: List[str] | None = None
    rows: List[Dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ValueError(f"Input table has no header: {path}")
            current = list(reader.fieldnames)
            missing = REQUIRED_COLUMNS.difference(current)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
            if header is None:
                header = current
            elif current != header:
                raise ValueError(f"Column order differs in input table: {path}")
            rows.extend(flag_row(dict(row)) for row in reader)
    if header is None:
        raise ValueError("No per-setting input tables were found")
    return header, rows


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine RP_get_coin_normyield per-setting CSVs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Per-setting table directory (default: {DEFAULT_INPUT_DIR})",
    )
    args = parser.parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    paths = sorted(
        path
        for path in input_dir.glob("phase*.csv")
        if path.name not in {SUMMARY_NAME, PROBLEMATIC_NAME}
    )
    header, rows = read_tables(paths)
    rows.sort(key=run_key)

    run_type_index = header.index("Run_type")
    output_header = [
        *header[: run_type_index + 1],
        *FLAG_COLUMNS,
        *header[run_type_index + 1 :],
    ]
    summary = input_dir / SUMMARY_NAME
    problematic = input_dir / PROBLEMATIC_NAME
    bad_rows = [row for row in rows if row["Problematic"] == "TRUE"]
    write_csv(summary, output_header, rows)
    write_csv(problematic, output_header, bad_rows)

    print(f"Input setting tables: {len(paths)}")
    print(f"All runs: {len(rows)} -> {summary}")
    print(f"Problematic runs: {len(bad_rows)} -> {problematic}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
