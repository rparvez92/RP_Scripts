#!/usr/bin/env python3
"""
Combined4and5PassErrors.py

Toy-MC uncertainties for the combined 4+5-pass SciPy/Nelder-Mead fit.
This uses the exact same model as Combined4and5PassFitParams.py:
all setting-dependent dthe/dpe/dpp offsets float, while dthp is global.

Inputs:
  results/tables/MeasuredOffsetsBySetting_4pass.csv
  results/tables/MeasuredOffsetsBySetting_5pass.csv
  results/tables/Combined4and5Pass/Combined4and5PassFitParams.csv

Outputs:
  results/tables/Combined4and5Pass/Combined4and5PassErrors_toys.csv
  results/tables/Combined4and5Pass/Combined4and5PassErrors_params.csv
  results/tables/Combined4and5Pass/Combined4and5PassErrors_summary.csv
  results/tables/Combined4and5Pass/Combined4and5PassResiduals_WithToyErrors.csv

Run from heep_check_v3:
  python3 macros/Combined4and5PassErrors.py
  python3 macros/Combined4and5PassErrors.py --n-toys 5000 --seed 314159
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Combined4and5PassFitParams as fitmod  # noqa: E402


BAD_BRANCH_ABS_PARAM_THRESHOLD = 5.0
OUT_DIR = Path("results/tables/Combined4and5Pass")
FITPARAMS_CSV = OUT_DIR / "Combined4and5PassFitParams.csv"
PARAMS_CSV = OUT_DIR / "Combined4and5PassErrors_params.csv"
TOYS_CSV = OUT_DIR / "Combined4and5PassErrors_toys.csv"
SUMMARY_CSV = OUT_DIR / "Combined4and5PassErrors_summary.csv"
RESIDUALS_CSV = OUT_DIR / "Combined4and5PassResiduals_WithToyErrors.csv"


def f4(x: float) -> str:
    return f"{x:.4f}" if np.isfinite(x) else "nan"


def read_central_params(path: Path) -> np.ndarray:
    if not path.exists():
        raise RuntimeError(f"Missing central fit CSV: {path}")
    with path.open() as f:
        row = next(csv.DictReader(f))
    return np.array([float(row[name]) for name in fitmod.PARAM_NAMES], dtype=float)


def predictions(params: np.ndarray, rows: list[fitmod.ObsRow]) -> np.ndarray:
    return np.array([fitmod.predict_one(params, row) for row in rows], dtype=float)


def bad_branch(params: np.ndarray) -> bool:
    return bool(np.max(np.abs(params)) > BAD_BRANCH_ABS_PARAM_THRESHOLD)


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q)) if values.size else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input4", type=Path, default=fitmod.INPUT_4PASS)
    parser.add_argument("--input5", type=Path, default=fitmod.INPUT_5PASS)
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    parser.add_argument("--fitparams", type=Path, default=FITPARAMS_CSV)
    parser.add_argument("--n-toys", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=314159)
    args = parser.parse_args()

    out_dir = args.outdir
    out_dir.mkdir(parents=True, exist_ok=True)
    params_csv = out_dir / PARAMS_CSV.name
    toys_csv = out_dir / TOYS_CSV.name
    summary_csv = out_dir / SUMMARY_CSV.name
    residuals_csv = out_dir / RESIDUALS_CSV.name

    rows = fitmod.read_inputs(args.input4, args.input5)
    central_params = read_central_params(args.fitparams)
    central_measured = np.array([row.measured for row in rows], dtype=float)
    measured_err = np.array([row.measured_err for row in rows], dtype=float)
    central_pred = predictions(central_params, rows)
    central_chi2 = fitmod.chi2(central_params, rows, central_measured)

    rng = np.random.default_rng(args.seed)
    toy_rows: list[list[object]] = []
    toy_params: list[np.ndarray] = []
    toy_preds: list[np.ndarray] = []
    success_flags: list[bool] = []
    bad_flags: list[bool] = []

    for itoy in range(args.n_toys):
        toy_measured = rng.normal(central_measured, measured_err)
        result = fitmod.fit_combined(rows, toy_measured)
        params = result.x
        fit_chi2 = fitmod.chi2(params, rows, toy_measured)
        is_bad = bad_branch(params)
        physical = bool(result.success and not is_bad)

        toy_params.append(params)
        toy_preds.append(predictions(params, rows))
        success_flags.append(bool(result.success))
        bad_flags.append(is_bad)
        toy_rows.append(
            [
                itoy,
                bool(result.success),
                str(result.message),
                int(result.nit),
                int(result.nfev),
                f4(fit_chi2),
                int(is_bad),
                int(physical),
                *[f4(x) for x in params],
            ]
        )

    toy_param_array = np.array(toy_params, dtype=float)
    toy_pred_array = np.array(toy_preds, dtype=float)
    success_mask = np.array(success_flags, dtype=bool)
    bad_mask = np.array(bad_flags, dtype=bool)
    physical_mask = success_mask & (~bad_mask)
    physical_params = toy_param_array[physical_mask]
    physical_preds = toy_pred_array[physical_mask]

    with toys_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "toy",
                "success",
                "message",
                "nit",
                "nfev",
                "chi2",
                "bad_branch_flag_abs_param_gt_5",
                "physical_branch_flag",
                *fitmod.PARAM_NAMES,
            ]
        )
        writer.writerows(toy_rows)

    with params_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "parameter",
                "central_value",
                "toy_mean",
                "toy_std",
                "p16",
                "p50",
                "p84",
                "minus_err",
                "plus_err",
                "n_toys_physical",
            ]
        )
        for ipar, name in enumerate(fitmod.PARAM_NAMES):
            vals = physical_params[:, ipar] if physical_params.size else np.array([], dtype=float)
            p16 = percentile(vals, 16)
            p50 = percentile(vals, 50)
            p84 = percentile(vals, 84)
            writer.writerow(
                [
                    name,
                    f4(central_params[ipar]),
                    f4(float(np.mean(vals))) if vals.size else "nan",
                    f4(float(np.std(vals, ddof=1))) if vals.size > 1 else "nan",
                    f4(p16),
                    f4(p50),
                    f4(p84),
                    f4(central_params[ipar] - p16),
                    f4(p84 - central_params[ipar]),
                    int(vals.size),
                ]
            )

    pred_errs = (
        np.std(physical_preds, axis=0, ddof=1)
        if len(physical_preds) > 1
        else np.full(len(rows), np.nan)
    )
    with residuals_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "pass",
                "setting",
                "shms_p_GeV",
                "dp_label",
                "delta_center",
                "dp_lo",
                "dp_hi",
                "variable",
                "measured",
                "measured_err",
                "predicted",
                "predicted_err_FromToy",
                "residual",
                "residual_err_Total",
                "pull_Total",
            ]
        )
        for i, row in enumerate(rows):
            resid = row.measured - central_pred[i]
            residual_err = np.sqrt(row.measured_err**2 + pred_errs[i] ** 2)
            pull = resid / residual_err
            writer.writerow(
                [
                    row.pass_tag,
                    row.setting,
                    f4(row.shms_p_GeV),
                    row.dp_label,
                    f4(row.delta_center),
                    f4(row.dp_lo),
                    f4(row.dp_hi),
                    row.variable,
                    f4(row.measured),
                    f4(row.measured_err),
                    f4(central_pred[i]),
                    f4(pred_errs[i]),
                    f4(resid),
                    f4(residual_err),
                    f4(pull),
                ]
            )

    n_obs = len(rows)
    n_params = len(fitmod.PARAM_NAMES)
    ndf = n_obs - n_params
    with summary_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        writer.writerow(["method", "combined_4and5_scipy_nelder_mead_toy_errors"])
        writer.writerow(["input4_csv", str(args.input4)])
        writer.writerow(["input5_csv", str(args.input5)])
        writer.writerow(["fitparams_csv", str(args.fitparams)])
        writer.writerow(["n_toys_requested", args.n_toys])
        writer.writerow(["seed", args.seed])
        writer.writerow(["bad_branch_abs_param_threshold", f4(BAD_BRANCH_ABS_PARAM_THRESHOLD)])
        writer.writerow(["central_chi2", f4(central_chi2)])
        writer.writerow(["central_n_obs", n_obs])
        writer.writerow(["central_n_params", n_params])
        writer.writerow(["central_ndf", ndf])
        writer.writerow(["central_chi2_ndf", f4(central_chi2 / ndf)])
        writer.writerow(["central_bad_branch_flag_abs_param_gt_5", int(bad_branch(central_params))])
        writer.writerow(["toy_success_count", int(np.sum(success_mask))])
        writer.writerow(["toy_bad_branch_count_abs_param_gt_5", int(np.sum(bad_mask))])
        writer.writerow(["toy_physical_branch_count", int(np.sum(physical_mask))])
        writer.writerow(["toy_physical_branch_fraction", f4(float(np.mean(physical_mask)))])

    print(f"[INFO] Wrote {params_csv}")
    print(f"[INFO] Wrote {toys_csv}")
    print(f"[INFO] Wrote {summary_csv}")
    print(f"[INFO] Wrote {residuals_csv}")
    print(f"[INFO] combined 4+5 physical toy fraction: {float(np.mean(physical_mask)):.4f}")


if __name__ == "__main__":
    main()
