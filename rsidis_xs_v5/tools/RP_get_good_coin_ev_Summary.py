#!/usr/bin/env python3
"""Combine RP_get_good_coin_ev setting tables and flag problematic runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = PROJECT_DIR / "results" / "Tables" / "RP_get_good_coin_ev"
SUMMARY_NAME = "RP_get_good_coin_ev_Summary.csv"
PROBLEMATIC_NAME = "RP_get_good_coin_ev_Problematic.csv"

REQUIRED_COLUMNS = {
    "Run",
    "CT_method",
    "CT_POS_ref_status",
    "CT_POS_ref_count",
    "CT_POS_ref_runs",
    "CT_POS_ref_ctmean_stddev",
    "CT_POS_ref_ctmean_range",
    "CTmean",
    "CTmean_err",
    "CTsigma",
    "CTsigma_err",
    "CTmean_residual",
    "N_rndm_peak",
    "RP_Goodcoin",
    "ransubcoin",
    "ransubcoin_by_RP_Goodcoin",
    "Fit_status",
}

FLAG_COLUMNS = (
    "Residual_flag",
    "Ratio_flag",
    "Fit_flag",
    "CT_POS_ref_flag",
    "Random_window_flag",
    "Yield_flag",
    "Numeric_flag",
    "Problematic",
    "Problem_reason",
)


def parse_finite(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def true_false(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def flag_row(row: Dict[str, str]) -> Dict[str, str]:
    method = row.get("CT_method", "").strip()
    is_positron_reference = method == "PRIOR_ELEC_AVERAGE"
    residual = parse_finite(row.get("CTmean_residual", ""))
    ratio = parse_finite(row.get("ransubcoin_by_RP_Goodcoin", ""))
    goodcoin = parse_finite(row.get("RP_Goodcoin", ""))
    sigma = parse_finite(row.get("CTsigma", ""))

    numerical_fields = [
        "CTmean",
        "CTmean_err",
        "CTsigma",
        "CTsigma_err",
    ]
    if not is_positron_reference:
        numerical_fields.append("CTmean_residual")
    numeric_flag = any(parse_finite(row.get(name, "")) is None for name in numerical_fields)
    if sigma is not None and sigma <= 0.0:
        numeric_flag = True

    residual_flag = (
        not is_positron_reference
        and residual is not None
        and abs(residual) > 2.0
    )
    ratio_flag = (
        not is_positron_reference
        and (ratio is None or abs(ratio - 1.0) > 0.30)
    )
    expected_fit_status = "NOT_APPLICABLE" if is_positron_reference else "OK"
    fit_flag = row.get("Fit_status", "").strip() != expected_fit_status

    ref_status = row.get("CT_POS_ref_status", "").strip()
    try:
        ref_count = int(row.get("CT_POS_ref_count", ""))
    except (TypeError, ValueError):
        ref_count = -1
    ref_status_usable = ref_status == "OK" or ref_status.startswith("WARNING_")
    positron_ref_flag = is_positron_reference and (
        ref_count != 5 or not ref_status_usable
    )

    try:
        random_window_flag = int(row.get("N_rndm_peak", "")) != 6
    except (TypeError, ValueError):
        random_window_flag = True

    yield_flag = goodcoin is None or goodcoin <= 0.0

    reasons: List[str] = []
    if residual_flag:
        reasons.append("ABS_CTMEAN_RESIDUAL_GT_2NS")
    if ratio_flag:
        reasons.append("RATIO_INVALID_OR_DEVIATION_GT_30PCT")
    if fit_flag:
        reasons.append("FIT_STATUS_UNEXPECTED_FOR_CT_METHOD")
    if positron_ref_flag:
        reasons.append("POSITRON_REFERENCE_INVALID")
    if random_window_flag:
        reasons.append("N_RNDM_PEAK_NOT_6")
    if yield_flag:
        reasons.append("RP_GOODCOIN_NONPOSITIVE_OR_INVALID")
    if numeric_flag:
        reasons.append("INVALID_FIT_NUMERIC_VALUE")

    row.update(
        {
            "Residual_flag": true_false(residual_flag),
            "Ratio_flag": true_false(ratio_flag),
            "Fit_flag": true_false(fit_flag),
            "CT_POS_ref_flag": true_false(positron_ref_flag),
            "Random_window_flag": true_false(random_window_flag),
            "Yield_flag": true_false(yield_flag),
            "Numeric_flag": true_false(numeric_flag),
            "Problematic": true_false(bool(reasons)),
            "Problem_reason": ";".join(reasons),
        }
    )
    return row


def run_sort_key(row: Dict[str, str]) -> tuple[int, str]:
    raw_run = row.get("Run", "")
    try:
        return int(raw_run), raw_run
    except ValueError:
        return 2**63 - 1, raw_run


def read_tables(paths: Sequence[Path]) -> tuple[List[str], List[Dict[str, str]]]:
    common_header: List[str] | None = None
    rows: List[Dict[str, str]] = []

    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ValueError(f"Input table has no header: {path}")
            header = list(reader.fieldnames)
            missing = REQUIRED_COLUMNS.difference(header)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
            if common_header is None:
                common_header = header
            elif header != common_header:
                raise ValueError(f"Column order differs in input table: {path}")
            rows.extend(flag_row(dict(row)) for row in reader)

    if common_header is None:
        raise ValueError("No input tables were found")
    return common_header, rows


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine RP_get_good_coin_ev per-setting CSVs and flag problematic runs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Per-setting CSV directory (default: {DEFAULT_INPUT_DIR})",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    input_paths = sorted(input_dir.glob("phase*.csv"))
    header, rows = read_tables(input_paths)
    rows.sort(key=run_sort_key)
    base_header = [column for column in header if column != "Fit_status"]
    run_type_index = base_header.index("Run_type")
    output_header = [
        *base_header[: run_type_index + 1],
        "Fit_status",
        *FLAG_COLUMNS,
        *base_header[run_type_index + 1 :],
    ]

    summary_path = input_dir / SUMMARY_NAME
    problematic_path = input_dir / PROBLEMATIC_NAME
    problematic_rows = [row for row in rows if row["Problematic"] == "TRUE"]

    write_csv(summary_path, output_header, rows)
    write_csv(problematic_path, output_header, problematic_rows)

    print(f"Input setting tables: {len(input_paths)}")
    print(f"All runs: {len(rows)} -> {summary_path}")
    print(f"Problematic runs: {len(problematic_rows)} -> {problematic_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
