#!/usr/bin/env python3
"""
SimultaneousDpeDtheErrors.py

Toy-MC uncertainties for the combined SciPy/Nelder-Mead dpe/dthe fits.

Supported modes:
  4-pass:
    input central fit:
      results/tables/SimultaneousDpeDtheFitParams_4pass.csv
    outputs:
      results/tables/SimultaneousDpeDtheErrors_4pass_toys.csv
      results/tables/SimultaneousDpeDtheErrors_4pass_params.csv
      results/tables/SimultaneousDpeDtheErrors_4pass_summary.csv
      results/tables/SimultaneousDpeDtheResiduals_4pass.csv

  5-pass:
    input central fit:
      results/tables/SimultaneousDpeDtheFitParams_5pass.csv
    outputs:
      results/tables/SimultaneousDpeDtheErrors_5pass_toys.csv
      results/tables/SimultaneousDpeDtheErrors_5pass_params.csv
      results/tables/SimultaneousDpeDtheErrors_5pass_summary.csv
      results/tables/SimultaneousDpeDtheResiduals_5pass.csv

Run from heep_check_v3:
  python3 macros/SimultaneousDpeDtheErrors.py --pass 4pass
  python3 macros/SimultaneousDpeDtheErrors.py --pass 5pass

Optional:
  python3 macros/SimultaneousDpeDtheErrors.py --pass 4pass --n-toys 5000 --seed 314159
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import minimize


OBSERVABLES = ("W", "Em", "Pmz", "Pmx")
REQUIRED_INPUT_COLUMNS = {
    "setting",
    "shms_p_GeV",
    "dp_label",
    "dp_lo",
    "dp_hi",
    "bin_center",
    "var",
    "offset_MeV",
    "offset_err_MeV",
    "sigD_MeV",
    "sigD_err_MeV",
    "fit_valid",
}
OPT_OPTIONS = {"maxiter": 30000, "maxfev": 50000, "xatol": 1e-10, "fatol": 1e-10}
BAD_BRANCH_ABS_PARAM_THRESHOLD = 5.0

# Rows: W, Em, Pmz, Pmx
# Columns: dthe, dpe, dthp, dpp
RESPONSE_5PASS = np.array(
    [
        [-24.63, -10.73, 0.00, 0.00],
        [0.00, -6.71, 0.00, -4.73],
        [-4.73, 4.75, 0.00, 4.81],
        [4.75, 4.73, -4.81, 0.00],
    ],
    dtype=float,
)

RESPONSE_4PASS = {
    "A": np.array(
        [
            [-14.08, -8.62, 0.00, 0.00],
            [0.00, -7.06, 0.00, -2.10],
            [-5.75, 4.10, 0.00, 2.27],
            [4.10, 5.75, -2.27, 0.00],
        ],
        dtype=float,
    ),
    "B": np.array(
        [
            [-17.33, -8.62, 0.00, 0.00],
            [0.00, -5.66, 0.00, -3.63],
            [-4.30, 3.69, 0.00, 3.74],
            [3.69, 4.30, -3.74, 0.00],
        ],
        dtype=float,
    ),
}


@dataclass(frozen=True)
class ObsRow:
    setting: str
    shms_p_GeV: float
    dp_label: str
    delta_center: float
    dp_lo: float
    dp_hi: float
    variable: str
    measured: float
    measured_err: float
    sigma: float
    sigma_err: float


@dataclass(frozen=True)
class FitConfig:
    pass_tag: str
    input_csv: Path
    fitparams_csv: Path
    params_csv: Path
    toys_csv: Path
    summary_csv: Path
    residuals_csv: Path
    settings: tuple[str, ...]
    dp_labels: tuple[str, ...]
    param_names: tuple[str, ...]
    param_vector: Callable[[np.ndarray, ObsRow], np.ndarray]
    response_matrix: Callable[[ObsRow], np.ndarray]


def f4(x: float) -> str:
    return f"{x:.4f}" if np.isfinite(x) else "nan"


def normalize_pass_tag(pass_arg: str) -> str:
    tag = pass_arg.strip().lower()
    if tag in {"4", "4pass", "pass4"}:
        return "4pass"
    if tag in {"5", "5pass", "pass5"}:
        return "5pass"
    raise ValueError(f"Unknown pass '{pass_arg}'. Use 4pass or 5pass.")


def make_config(pass_arg: str) -> FitConfig:
    pass_tag = normalize_pass_tag(pass_arg)
    if pass_tag == "4pass":
        dp_labels = ("b1",)
        param_names = (
            "dthe_b1_A",
            "dthe_b1_B",
            "dpe_b1_A",
            "dpe_b1_B",
            "dthp",
            "dpp_b1_A",
            "dpp_b1_B",
        )

        def param_vector(params: np.ndarray, row: ObsRow) -> np.ndarray:
            if row.setting == "A":
                return np.array([params[0], params[2], params[4], params[5]], dtype=float)
            if row.setting == "B":
                return np.array([params[1], params[3], params[4], params[6]], dtype=float)
            raise RuntimeError(f"Unexpected 4-pass setting: {row.setting}")

        def response_matrix(row: ObsRow) -> np.ndarray:
            return RESPONSE_4PASS[row.setting]

        return FitConfig(
            pass_tag=pass_tag,
            input_csv=Path("results/tables/MeasuredOffsetsBySetting_4pass.csv"),
            fitparams_csv=Path("results/tables/SimultaneousDpeDtheFitParams_4pass.csv"),
            params_csv=Path("results/tables/SimultaneousDpeDtheErrors_4pass_params.csv"),
            toys_csv=Path("results/tables/SimultaneousDpeDtheErrors_4pass_toys.csv"),
            summary_csv=Path("results/tables/SimultaneousDpeDtheErrors_4pass_summary.csv"),
            residuals_csv=Path("results/tables/SimultaneousDpeDtheResiduals_4pass.csv"),
            settings=("A", "B"),
            dp_labels=dp_labels,
            param_names=param_names,
            param_vector=param_vector,
            response_matrix=response_matrix,
        )

    settings = ("C", "E", "D")
    param_names = ("dthe_C", "dthe_E", "dthe_D", "dpe_C", "dpe_E", "dpe_D", "dthp", "dpp")

    def param_vector(params: np.ndarray, row: ObsRow) -> np.ndarray:
        iset = settings.index(row.setting)
        return np.array([params[iset], params[3 + iset], params[6], params[7]], dtype=float)

    def response_matrix(row: ObsRow) -> np.ndarray:
        return RESPONSE_5PASS

    return FitConfig(
        pass_tag=pass_tag,
        input_csv=Path("results/tables/MeasuredOffsetsBySetting_5pass.csv"),
        fitparams_csv=Path("results/tables/SimultaneousDpeDtheFitParams_5pass.csv"),
        params_csv=Path("results/tables/SimultaneousDpeDtheErrors_5pass_params.csv"),
        toys_csv=Path("results/tables/SimultaneousDpeDtheErrors_5pass_toys.csv"),
        summary_csv=Path("results/tables/SimultaneousDpeDtheErrors_5pass_summary.csv"),
        residuals_csv=Path("results/tables/SimultaneousDpeDtheResiduals_5pass.csv"),
        settings=settings,
        dp_labels=("b1",),
        param_names=param_names,
        param_vector=param_vector,
        response_matrix=response_matrix,
    )


def require_input_columns(reader: csv.DictReader, input_csv: Path) -> None:
    present = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - present)
    if missing:
        raise RuntimeError(f"Missing required columns in {input_csv}: {missing}")


def read_input(config: FitConfig) -> list[ObsRow]:
    by_key: dict[tuple[str, str, str], ObsRow] = {}
    with config.input_csv.open() as f:
        reader = csv.DictReader(f)
        require_input_columns(reader, config.input_csv)
        for row in reader:
            setting = row["setting"]
            dp_label = row["dp_label"]
            variable = row["var"]
            if setting not in config.settings or dp_label not in config.dp_labels or variable not in OBSERVABLES:
                continue
            if row["fit_valid"] != "1":
                continue
            by_key[(setting, dp_label, variable)] = ObsRow(
                setting=setting,
                shms_p_GeV=float(row["shms_p_GeV"]),
                dp_label=dp_label,
                delta_center=float(row["bin_center"]),
                dp_lo=float(row["dp_lo"]),
                dp_hi=float(row["dp_hi"]),
                variable=variable,
                measured=float(row["offset_MeV"]),
                measured_err=float(row["offset_err_MeV"]),
                sigma=float(row["sigD_MeV"]),
                sigma_err=float(row["sigD_err_MeV"]),
            )

    rows: list[ObsRow] = []
    missing: list[tuple[str, str, str]] = []
    for setting in config.settings:
        for dp_label in config.dp_labels:
            for variable in OBSERVABLES:
                key = (setting, dp_label, variable)
                if key not in by_key:
                    missing.append(key)
                else:
                    rows.append(by_key[key])
    if missing:
        raise RuntimeError(f"Missing required rows in {config.input_csv}: {missing}")
    return rows


def read_central_params(config: FitConfig) -> np.ndarray:
    if not config.fitparams_csv.exists():
        raise RuntimeError(f"Missing central fit CSV: {config.fitparams_csv}")
    with config.fitparams_csv.open() as f:
        first = next(csv.DictReader(f))
    return np.array([float(first[name]) for name in config.param_names], dtype=float)


def predict_one(config: FitConfig, params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(config.response_matrix(row)[ivar] @ config.param_vector(params, row))


def predictions(config: FitConfig, params: np.ndarray, rows: list[ObsRow]) -> np.ndarray:
    return np.array([predict_one(config, params, row) for row in rows], dtype=float)


def chi2(config: FitConfig, params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    pred = predictions(config, params, rows)
    err = np.array([row.measured_err for row in rows], dtype=float)
    return float(np.sum(((measured - pred) / err) ** 2))


def fit_combined(config: FitConfig, rows: list[ObsRow], measured: np.ndarray):
    return minimize(
        lambda p: chi2(config, p, rows, measured),
        np.zeros(len(config.param_names), dtype=float),
        method="Nelder-Mead",
        options=OPT_OPTIONS,
    )


def bad_branch(params: np.ndarray) -> bool:
    return bool(np.max(np.abs(params)) > BAD_BRANCH_ABS_PARAM_THRESHOLD)


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q)) if values.size else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pass", dest="pass_arg", required=True, help="4pass or 5pass")
    parser.add_argument("--n-toys", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=314159)
    args = parser.parse_args()

    config = make_config(args.pass_arg)
    rows = read_input(config)
    central_params = read_central_params(config)
    central_measured = np.array([row.measured for row in rows], dtype=float)
    measured_err = np.array([row.measured_err for row in rows], dtype=float)
    central_pred = predictions(config, central_params, rows)
    central_chi2 = chi2(config, central_params, rows, central_measured)

    for path in (config.params_csv, config.toys_csv, config.summary_csv, config.residuals_csv):
        path.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    toy_param_rows: list[list[object]] = []
    toy_params: list[np.ndarray] = []
    toy_preds: list[np.ndarray] = []
    success_flags: list[bool] = []
    bad_flags: list[bool] = []

    for itoy in range(args.n_toys):
        toy_measured = rng.normal(central_measured, measured_err)
        result = fit_combined(config, rows, toy_measured)
        params = result.x
        fit_chi2 = chi2(config, params, rows, toy_measured)
        is_bad = bad_branch(params)
        physical = bool(result.success and not is_bad)
        toy_params.append(params)
        toy_preds.append(predictions(config, params, rows))
        success_flags.append(bool(result.success))
        bad_flags.append(is_bad)
        toy_param_rows.append(
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

    with config.toys_csv.open("w", newline="") as f:
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
                *config.param_names,
            ]
        )
        writer.writerows(toy_param_rows)

    with config.params_csv.open("w", newline="") as f:
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
        for ipar, name in enumerate(config.param_names):
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

    pred_errs = np.std(physical_preds, axis=0, ddof=1) if len(physical_preds) > 1 else np.full(len(rows), np.nan)
    with config.residuals_csv.open("w", newline="") as f:
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
                    config.pass_tag,
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

    with config.summary_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        writer.writerow(["method", "scipy_nelder_mead_toy_errors"])
        writer.writerow(["pass", config.pass_tag])
        writer.writerow(["input_csv", str(config.input_csv)])
        writer.writerow(["fitparams_csv", str(config.fitparams_csv)])
        writer.writerow(["n_toys_requested", args.n_toys])
        writer.writerow(["seed", args.seed])
        writer.writerow(["bad_branch_abs_param_threshold", f4(BAD_BRANCH_ABS_PARAM_THRESHOLD)])
        writer.writerow(["central_chi2", f4(central_chi2)])
        writer.writerow(["central_n_obs", len(rows)])
        writer.writerow(["central_n_params", len(config.param_names)])
        writer.writerow(["central_ndf", len(rows) - len(config.param_names)])
        writer.writerow(["central_chi2_ndf", f4(central_chi2 / (len(rows) - len(config.param_names)))])
        writer.writerow(["central_bad_branch_flag_abs_param_gt_5", int(bad_branch(central_params))])
        writer.writerow(["toy_success_count", int(np.sum(success_mask))])
        writer.writerow(["toy_bad_branch_count_abs_param_gt_5", int(np.sum(bad_mask))])
        writer.writerow(["toy_physical_branch_count", int(np.sum(physical_mask))])
        writer.writerow(["toy_physical_branch_fraction", f4(float(np.mean(physical_mask)))])

    print(f"[INFO] Wrote {config.params_csv}")
    print(f"[INFO] Wrote {config.toys_csv}")
    print(f"[INFO] Wrote {config.summary_csv}")
    print(f"[INFO] Wrote {config.residuals_csv}")
    print(f"[INFO] {config.pass_tag} physical toy fraction: {float(np.mean(physical_mask)):.4f}")


if __name__ == "__main__":
    main()
