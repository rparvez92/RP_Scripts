#!/usr/bin/env python3
"""
PlotSimultaneousDpeOnly.py

Python version of macros/PlotSimultaneousDpeOnly.C.

This script reads the CSV produced by SimultaneousDpeOnly.py and writes three
PNG files per setting:

  results/Python/SimultaneousDpeOnly/settingA_kinOffsets.png
  results/Python/SimultaneousDpeOnly/settingA_residuals.png
  results/Python/SimultaneousDpeOnly/settingA_chi2ndf.png

and the corresponding setting B files.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, List

import matplotlib.pyplot as plt


@dataclass
class Row:
    setting: str
    bin: str
    dthe_fit: float
    dthe_fit_err: float
    dthp_fit: float
    dthp_fit_err: float
    dpp_fit: float
    dpp_fit_err: float
    dpe_fit: float
    dpe_fit_err: float
    dW_resid: float
    dW_resid_err: float
    dEm_resid: float
    dEm_resid_err: float
    dPmz_resid: float
    dPmz_resid_err: float
    dPmy_resid: float
    dPmy_resid_err: float
    chi2_ndf: float
    fit_valid: int


def normalize_setting(text: str) -> str:
    t = text.strip().lower()
    if t in {"a", "setting a", "settinga"}:
        return "A"
    if t in {"b", "setting b", "settingb"}:
        return "B"
    return text.strip()


def to_float(text: str) -> float:
    try:
        return float(text.strip())
    except Exception:
        return math.nan


def to_int(text: str, default: int = 0) -> int:
    try:
        return int(text.strip())
    except Exception:
        return default


def bin_order(bin_label: str) -> int:
    return {"b1": 1, "b2": 2, "b3": 3, "b4": 4, "b5": 5}.get(bin_label.strip().lower(), 99)


def bin_center(bin_label: str) -> float:
    return {
        "b1": -1.5,
        "b2": -0.5,
        "b3": 0.5,
        "b4": 1.5,
        "b5": 2.5,
    }.get(bin_label.strip().lower(), math.nan)


def get_value(row: Row, key: str) -> float:
    return {
        "dthe": row.dthe_fit,
        "dpe": row.dpe_fit,
        "dthp": row.dthp_fit,
        "dpp": row.dpp_fit,
        "W": row.dW_resid,
        "Em": row.dEm_resid,
        "Pmz": row.dPmz_resid,
        "Pmy": row.dPmy_resid,
        "chi2ndf": row.chi2_ndf,
    }[key]


def get_error(row: Row, key: str) -> float:
    return {
        "dthe": row.dthe_fit_err,
        "dpe": row.dpe_fit_err,
        "dthp": row.dthp_fit_err,
        "dpp": row.dpp_fit_err,
        "W": row.dW_resid_err,
        "Em": row.dEm_resid_err,
        "Pmz": row.dPmz_resid_err,
        "Pmy": row.dPmy_resid_err,
    }.get(key, 0.0)


def sort_rows(rows: List[Row]) -> List[Row]:
    return sorted(rows, key=lambda row: bin_order(row.bin))


def compute_y_range(rows: List[Row], key: str, use_errors: bool) -> tuple[float, float]:
    ymin = math.inf
    ymax = -math.inf
    for row in rows:
        y = get_value(row, key)
        if not math.isfinite(y):
            continue
        ey = get_error(row, key) if use_errors else 0.0
        if not math.isfinite(ey):
            ey = 0.0
        ymin = min(ymin, y - ey)
        ymax = max(ymax, y + ey)

    if not (math.isfinite(ymin) and math.isfinite(ymax)):
        return -1.0, 1.0
    if abs(ymax - ymin) < 1e-12:
        pad = 0.25 * abs(ymax) if abs(ymax) > 0.0 else 1.0
        return ymin - pad, ymax + pad
    span = ymax - ymin
    pad = 0.18 * span
    return ymin - pad, ymax + pad


def read_simultaneous_csv(csv_path: str) -> List[Row]:
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Empty CSV: {csv_path}")

        required = {
            "setting",
            "bin",
            "dthe_fit",
            "dthe_fit_err",
            "dthp_fit",
            "dthp_fit_err",
            "dpp_fit",
            "dpp_fit_err",
            "dw_resid",
            "dw_resid_err",
            "dem_resid",
            "dem_resid_err",
            "dpmz_resid",
            "dpmz_resid_err",
            "dpmy_resid",
            "dpmy_resid_err",
            "chi2_ndf",
            "fit_valid",
        }
        lower_map = {name.lower(): name for name in reader.fieldnames}
        missing = sorted(required - set(lower_map))
        if missing:
            raise RuntimeError(
                f"Missing required columns in {csv_path}: {', '.join(missing)}"
            )

        rows: List[Row] = []
        for raw_row in reader:
            setting = normalize_setting(raw_row[lower_map["setting"]])
            bin_label = raw_row[lower_map["bin"]].strip()
            dpe_col = lower_map.get(f"dpe_{bin_label.lower()}_fit")
            dpe_err_col = lower_map.get(f"dpe_{bin_label.lower()}_fit_err")
            rows.append(
                Row(
                    setting=setting,
                    bin=bin_label,
                    dthe_fit=to_float(raw_row[lower_map["dthe_fit"]]),
                    dthe_fit_err=to_float(raw_row[lower_map["dthe_fit_err"]]),
                    dthp_fit=to_float(raw_row[lower_map["dthp_fit"]]),
                    dthp_fit_err=to_float(raw_row[lower_map["dthp_fit_err"]]),
                    dpp_fit=to_float(raw_row[lower_map["dpp_fit"]]),
                    dpp_fit_err=to_float(raw_row[lower_map["dpp_fit_err"]]),
                    dpe_fit=to_float(raw_row[dpe_col]) if dpe_col else math.nan,
                    dpe_fit_err=to_float(raw_row[dpe_err_col]) if dpe_err_col else math.nan,
                    dW_resid=to_float(raw_row[lower_map["dw_resid"]]),
                    dW_resid_err=to_float(raw_row[lower_map["dw_resid_err"]]),
                    dEm_resid=to_float(raw_row[lower_map["dem_resid"]]),
                    dEm_resid_err=to_float(raw_row[lower_map["dem_resid_err"]]),
                    dPmz_resid=to_float(raw_row[lower_map["dpmz_resid"]]),
                    dPmz_resid_err=to_float(raw_row[lower_map["dpmz_resid_err"]]),
                    dPmy_resid=to_float(raw_row[lower_map["dpmy_resid"]]),
                    dPmy_resid_err=to_float(raw_row[lower_map["dpmy_resid_err"]]),
                    chi2_ndf=to_float(raw_row[lower_map["chi2_ndf"]]),
                    fit_valid=to_int(raw_row[lower_map["fit_valid"]], 0),
                )
            )
    return rows


def draw_series(ax, rows: List[Row], key: str, pad_title: str, y_title: str, use_errors: bool) -> None:
    xvals: List[float] = []
    yvals: List[float] = []
    yerrs: List[float] = []
    for row in sort_rows(rows):
        x = bin_center(row.bin)
        y = get_value(row, key)
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        xvals.append(x)
        yvals.append(y)
        err = get_error(row, key) if use_errors else 0.0
        yerrs.append(err if math.isfinite(err) else 0.0)

    ymin, ymax = compute_y_range(rows, key, use_errors)
    ax.errorbar(
        xvals,
        yvals,
        yerr=yerrs if use_errors else None,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.5,
        capsize=3,
        markersize=5,
    )
    ax.set_xlim(-2.2, 3.2)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("#delta bin center (%)")
    ax.set_ylabel(y_title)
    ax.set_title(pad_title)
    ax.grid(True)


def make_offsets_figure(setting: str, rows: List[Row], out_dir: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    defs = [
        ("dthe", "dthe", "Offset (mrad)"),
        ("dpe", "dpe", "Offset (0.1%)"),
        ("dthp", "dthp", "Offset (mrad)"),
        ("dpp", "dpp", "Offset (0.1%)"),
    ]
    for ax, (key, title, ylabel) in zip(axes.flat, defs):
        draw_series(ax, rows, key, title, ylabel, True)
    fig.savefig(os.path.join(out_dir, f"setting{setting}_kinOffsets.png"), dpi=150)
    plt.close(fig)


def make_residuals_figure(setting: str, rows: List[Row], out_dir: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    defs = [
        ("W", "W", "Residuals (MeV)"),
        ("Em", "Em", "Residuals (MeV)"),
        ("Pmz", "Pmz", "Residuals (MeV)"),
        ("Pmy", "Pmy", "Residuals (MeV)"),
    ]
    for ax, (key, title, ylabel) in zip(axes.flat, defs):
        draw_series(ax, rows, key, title, ylabel, True)
    fig.savefig(os.path.join(out_dir, f"setting{setting}_residuals.png"), dpi=150)
    plt.close(fig)


def make_chi2_figure(setting: str, rows: List[Row], out_dir: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    draw_series(ax, rows, "chi2ndf", "chi2/ndf", "chi2/ndf", False)
    fig.savefig(os.path.join(out_dir, f"setting{setting}_chi2ndf.png"), dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Python SimultaneousDpeOnly fit results.")
    parser.add_argument(
        "--in",
        dest="in_csv",
        default="results/tables/SimultaneousDpeOnly_python.csv",
        help="Input CSV from SimultaneousDpeOnly.py.",
    )
    parser.add_argument(
        "--outdir",
        dest="out_dir",
        default="results/Python/SimultaneousDpeOnly",
        help="Output directory for PNGs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        all_rows = read_simultaneous_csv(args.in_csv)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    if not all_rows:
        print(f"[ERROR] No rows read from {args.in_csv}", file=sys.stderr)
        return 1

    os.makedirs(args.out_dir, exist_ok=True)

    by_setting: Dict[str, List[Row]] = {}
    skipped_invalid = 0
    for row in all_rows:
        if row.fit_valid != 1:
            skipped_invalid += 1
            continue
        by_setting.setdefault(row.setting, []).append(row)

    if not by_setting:
        print(f"[ERROR] No valid rows found in {args.in_csv}", file=sys.stderr)
        return 1

    for setting in ["A", "B"]:
        if setting not in by_setting:
            continue
        make_offsets_figure(setting, by_setting[setting], args.out_dir)
        make_residuals_figure(setting, by_setting[setting], args.out_dir)
        make_chi2_figure(setting, by_setting[setting], args.out_dir)

    for setting, rows in by_setting.items():
        if setting in {"A", "B"}:
            continue
        make_offsets_figure(setting, rows, args.out_dir)
        make_residuals_figure(setting, rows, args.out_dir)
        make_chi2_figure(setting, rows, args.out_dir)

    print(
        f"[INFO] Wrote PlotSimultaneousDpeOnly PNGs to {args.out_dir} "
        f"(skipped invalid rows: {skipped_invalid})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
