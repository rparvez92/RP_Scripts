#!/usr/bin/env python3
"""Add RP coincidence-analysis results to the Pass1 SIDIS bigtable."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path
from typing import Dict, Iterable


PROJECT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_BIGTABLE = (
    PROJECT_DIR / "bigtable" / "Pass1" / "rp_rsidis_bigtable_pass1_sidis.csv"
)
DEFAULT_GOODCOIN = (
    PROJECT_DIR
    / "results"
    / "Pass1"
    / "Tables"
    / "RP_get_good_coin_ev"
    / "RP_get_good_coin_ev_Summary.csv"
)
DEFAULT_NORMYIELD = (
    PROJECT_DIR
    / "results"
    / "Pass1"
    / "Tables"
    / "RP_get_coin_normyield"
    / "RP_get_coin_normyield_Summary.csv"
)

NEW_COLUMNS = {
    "coin": ("RP_Coin",),
    "ransubcoin_err": ("RP_Goodcoin", "RP_Goodcoin_err"),
    "normyield_err": ("RP_Normyield", "RP_Normyield_err"),
    "ctsigma": ("CTmean", "CTsigma"),
}


def canonical_run(value: str) -> int:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid run number: {value!r}") from exc
    if not number.is_integer():
        raise ValueError(f"Non-integral run number: {value!r}")
    return int(number)


def one_decimal(value: str) -> str:
    try:
        return f"{float(value):.1f}"
    except ValueError as exc:
        raise ValueError(f"Invalid RP_Coin value: {value!r}") from exc


def read_csv(path: Path) -> tuple[list[str], list[Dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def index_rows(
    path: Path, required: Iterable[str]
) -> dict[int, Dict[str, str]]:
    header, rows = read_csv(path)
    missing = set(required).difference(header)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    indexed: dict[int, Dict[str, str]] = {}
    for row in rows:
        run = canonical_run(row["Run"])
        if run in indexed:
            raise ValueError(f"Duplicate run {run} in {path}")
        indexed[run] = row
    return indexed


def output_header(input_header: list[str]) -> list[str]:
    required = set(NEW_COLUMNS)
    missing = required.difference(input_header)
    if missing:
        raise ValueError(f"Bigtable is missing anchor columns: {sorted(missing)}")
    added = {name for names in NEW_COLUMNS.values() for name in names}
    header: list[str] = []
    for name in input_header:
        if name not in added:
            header.append(name)
        header.extend(NEW_COLUMNS.get(name, ()))
    return header


def add_columns(
    rows: list[Dict[str, str]],
    goodcoin: dict[int, Dict[str, str]],
    normyield: dict[int, Dict[str, str]],
    missing_value: str,
) -> tuple[int, int]:
    seen: set[int] = set()
    complete = 0
    for row in rows:
        run = canonical_run(row["run"])
        if run in seen:
            raise ValueError(f"Duplicate run {run} in bigtable")
        seen.add(run)
        good = goodcoin.get(run)
        norm = normyield.get(run)
        row["RP_Coin"] = one_decimal(good["N_coin"]) if good else missing_value
        row["RP_Goodcoin"] = good["RP_Goodcoin"] if good else missing_value
        row["RP_Goodcoin_err"] = good["RP_Goodcoin_err"] if good else missing_value
        row["RP_Normyield"] = norm["RP_Normyield"] if norm else missing_value
        row["RP_Normyield_err"] = norm["RP_Normyield_err"] if norm else missing_value
        row["CTmean"] = good["CTmean"] if good else missing_value
        row["CTsigma"] = good["CTsigma"] if good else missing_value
        complete += good is not None and norm is not None

    unknown_good = set(goodcoin).difference(seen)
    unknown_norm = set(normyield).difference(seen)
    if unknown_good or unknown_norm:
        raise ValueError(
            "RP summaries contain runs absent from the bigtable: "
            f"good-coin={sorted(unknown_good)}, normyield={sorted(unknown_norm)}"
        )
    return complete, len(rows) - complete


def write_atomic(path: Path, header: list[str], rows: list[Dict[str, str]]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=header, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Join Pass1 RP analysis values into the SIDIS bigtable."
    )
    parser.add_argument("--bigtable", type=Path, default=DEFAULT_BIGTABLE)
    parser.add_argument("--goodcoin-summary", type=Path, default=DEFAULT_GOODCOIN)
    parser.add_argument("--normyield-summary", type=Path, default=DEFAULT_NORMYIELD)
    parser.add_argument("--output", type=Path, help="Output path; defaults to in-place update")
    parser.add_argument("--missing-value", default="-999")
    args = parser.parse_args()

    bigtable = args.bigtable.expanduser().resolve()
    destination = args.output.expanduser().resolve() if args.output else bigtable
    header, rows = read_csv(bigtable)
    if "run" not in header:
        raise ValueError(f"Bigtable is missing run column: {bigtable}")
    goodcoin = index_rows(
        args.goodcoin_summary.expanduser().resolve(),
        ("Run", "N_coin", "RP_Goodcoin", "RP_Goodcoin_err", "CTmean", "CTsigma"),
    )
    normyield = index_rows(
        args.normyield_summary.expanduser().resolve(),
        ("Run", "RP_Normyield", "RP_Normyield_err"),
    )
    complete, missing = add_columns(rows, goodcoin, normyield, args.missing_value)
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(destination, output_header(header), rows)
    print(f"Bigtable rows: {len(rows)}")
    print(f"Rows with good-coin and normyield results: {complete}")
    print(f"Rows using {args.missing_value} for unavailable RP results: {missing}")
    print(f"Wrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
