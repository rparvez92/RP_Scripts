#!/usr/bin/env python3
"""
MeasuredOffsetsBySetting_4passScaled.py

Diagnostic CSV transform for 4-pass measured offsets.

The measured offsets themselves are left unchanged. Only offset_err_MeV is
recomputed by scaling the 4-pass Gaussian mean-fit errors to the per-observable
average entriesD/entriesS of the current 5-pass C/E/D measured-offset table.

Run from heep_check_v3:
  python3 macros/MeasuredOffsetsBySetting_4passScaled.py

Output:
  results/tables/Scaled/MeasuredOffsetsBySetting_4passScaled.csv
  results/tables/Scaled/MeasuredOffsetsBySetting_4passScaled_summary.csv
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


INPUT_4PASS = Path("results/tables/MeasuredOffsetsBySetting_4pass.csv")
INPUT_5PASS = Path("results/tables/MeasuredOffsetsBySetting_5pass.csv")
OUT_DIR = Path("results/tables/Scaled")
OUTPUT_CSV = OUT_DIR / "MeasuredOffsetsBySetting_4passScaled.csv"
SUMMARY_CSV = OUT_DIR / "MeasuredOffsetsBySetting_4passScaled_summary.csv"

OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
REQUIRED_COLUMNS = {
    "var",
    "entriesD",
    "entriesS",
    "muD_err_MeV",
    "muS_err_MeV",
    "offset_err_MeV",
}


def f6(x: float) -> str:
    return f"{x:.6f}" if math.isfinite(x) else "nan"


def require_columns(fieldnames: list[str] | None, path: Path) -> None:
    present = set(fieldnames or [])
    missing = sorted(REQUIRED_COLUMNS - present)
    if missing:
        raise RuntimeError(f"Missing required columns in {path}: {missing}")


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as f:
        reader = csv.DictReader(f)
        require_columns(reader.fieldnames, path)
        return list(reader.fieldnames or []), list(reader)


def average_5pass_entries(rows: list[dict[str, str]]) -> dict[str, tuple[float, float, int]]:
    entries_d: dict[str, list[float]] = defaultdict(list)
    entries_s: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        var = row["var"]
        if var not in OBSERVABLES:
            continue
        if row.get("fit_valid", "1") != "1":
            continue
        entries_d[var].append(float(row["entriesD"]))
        entries_s[var].append(float(row["entriesS"]))

    out: dict[str, tuple[float, float, int]] = {}
    for var in OBSERVABLES:
        if not entries_d[var] or not entries_s[var]:
            raise RuntimeError(f"No 5-pass reference entries found for {var}")
        avg_d = sum(entries_d[var]) / len(entries_d[var])
        avg_s = sum(entries_s[var]) / len(entries_s[var])
        out[var] = (avg_d, avg_s, len(entries_d[var]))
    return out


def main() -> None:
    header4, rows4 = read_rows(INPUT_4PASS)
    _, rows5 = read_rows(INPUT_5PASS)
    ref = average_5pass_entries(rows5)

    scaled_rows: list[dict[str, str]] = []
    for row in rows4:
        var = row["var"]
        if var not in ref:
            scaled_rows.append(dict(row))
            continue
        avg_d, avg_s, _ = ref[var]
        entries_d = float(row["entriesD"])
        entries_s = float(row["entriesS"])
        mu_d_err = float(row["muD_err_MeV"])
        mu_s_err = float(row["muS_err_MeV"])
        scaled_mu_d_err = mu_d_err * math.sqrt(entries_d / avg_d)
        scaled_mu_s_err = mu_s_err * math.sqrt(entries_s / avg_s)
        scaled_offset_err = math.sqrt(scaled_mu_d_err**2 + scaled_mu_s_err**2)

        out = dict(row)
        out["offset_err_MeV"] = f6(scaled_offset_err)
        scaled_rows.append(out)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header4)
        writer.writeheader()
        writer.writerows(scaled_rows)

    with SUMMARY_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["variable", "avg_entriesD_5pass", "avg_entriesS_5pass", "n_reference_rows"])
        for var in OBSERVABLES:
            avg_d, avg_s, nrows = ref[var]
            writer.writerow([var, f6(avg_d), f6(avg_s), nrows])

    print(f"[INFO] Wrote {OUTPUT_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")
    for var in OBSERVABLES:
        avg_d, avg_s, nrows = ref[var]
        print(f"[INFO] {var}: avg_entriesD_5pass={avg_d:.3f}, avg_entriesS_5pass={avg_s:.3f}, n={nrows}")


if __name__ == "__main__":
    main()
