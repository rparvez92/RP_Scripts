#!/usr/bin/env python3
"""
ScipyWithError.py

Bootstrap / toy-Monte-Carlo uncertainty estimate for the 5-pass combined
SciPy fit using Radwan's measured offsets.

Run from heep_check_v3:
  python3 macros/ScipyWithError.py

Optional:
  python3 macros/ScipyWithError.py --n-toys 5000 --seed 12345
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


INPUT_CSV = Path("results/tables/Julio/TestingRadwan_input.csv")
OUT_DIR = Path("results/tables/ScipyWithError")

PARAMS_CSV = OUT_DIR / "ScipyWithError_params.csv"
TOYS_CSV = OUT_DIR / "ScipyWithError_toys.csv"
RESIDUALS_CSV = OUT_DIR / "ScipyWithError_residuals.csv"
SUMMARY_CSV = OUT_DIR / "ScipyWithError_summary.csv"

SETTINGS = ("C", "E", "D")
OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
PARAM_NAMES = (
    "dpe_C",
    "dpe_E",
    "dpe_D",
    "dtheta_e_C",
    "dtheta_e_E",
    "dtheta_e_D",
    "dtheta_p_common",
    "dpp_common",
)

# Rows: W, Em, Pmz, Pmy
# Columns: dtheta_e, dpe, dtheta_p, dpp
RESPONSE = np.array(
    [
        [-24.63, -10.73, 0.00, 0.00],
        [0.00, -6.71, 0.00, -4.73],
        [-4.73, 4.75, 0.00, 4.81],
        [4.75, 4.73, -4.81, 0.00],
    ],
    dtype=float,
)

ZERO_INIT = np.zeros(8, dtype=float)


@dataclass(frozen=True)
class ObsRow:
    setting: str
    shms_p_GeV: float
    delta_center: float
    dp_lo: float
    dp_hi: float
    variable: str
    mu: float
    mu_err: float
    sigma: float
    sigma_err: float


def parameter_vector_for_setting(params: np.ndarray, setting: str) -> np.ndarray:
    iset = SETTINGS.index(setting)
    return np.array(
        [
            params[3 + iset],  # dtheta_e_setting
            params[iset],      # dpe_setting
            params[6],         # dtheta_p_common
            params[7],         # dpp_common
        ],
        dtype=float,
    )


def predict_one(params: np.ndarray, setting: str, variable: str) -> float:
    ivar = OBSERVABLES.index(variable)
    return float(RESPONSE[ivar] @ parameter_vector_for_setting(params, setting))


def chi2(params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    total = 0.0
    for i, row in enumerate(rows):
        pred = predict_one(params, row.setting, row.variable)
        total += ((measured[i] - pred) / row.mu_err) ** 2
    return float(total)


def exact_inversion_like(params: np.ndarray, fit_chi2: float) -> bool:
    return bool(np.max(np.abs(params)) > 20.0 or (fit_chi2 < 1.0e-6 and np.max(np.abs(params)) > 5.0))


def fit_combined(rows: list[ObsRow], measured: np.ndarray):
    return minimize(
        lambda p: chi2(p, rows, measured),
        ZERO_INIT,
        method="Nelder-Mead",
        options={"maxiter": 20000, "maxfev": 30000, "xatol": 1e-10, "fatol": 1e-10},
    )


def read_input() -> list[ObsRow]:
    rows: list[ObsRow] = []
    with INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["setting"] not in SETTINGS or row["variable"] not in OBSERVABLES:
                continue
            rows.append(
                ObsRow(
                    setting=row["setting"],
                    shms_p_GeV=float(row["shms_p_GeV"]),
                    delta_center=float(row["delta_center"]),
                    dp_lo=float(row["dp_lo"]),
                    dp_hi=float(row["dp_hi"]),
                    variable=row["variable"],
                    mu=float(row["mu"]),
                    mu_err=float(row["mu_err"]),
                    sigma=float(row["sigma"]),
                    sigma_err=float(row["sigma_err"]),
                )
            )

    rows_by_key = {(r.setting, r.variable): r for r in rows}
    ordered_rows = []
    missing = []
    for setting in SETTINGS:
      for variable in OBSERVABLES:
        key = (setting, variable)
        if key not in rows_by_key:
            missing.append(key)
        else:
            ordered_rows.append(rows_by_key[key])
    if missing:
        raise RuntimeError(f"Missing required rows in {INPUT_CSV}: {missing}")
    return ordered_rows


def vector_mu(rows: list[ObsRow]) -> np.ndarray:
    return np.array([row.mu for row in rows], dtype=float)


def vector_err(rows: list[ObsRow]) -> np.ndarray:
    return np.array([row.mu_err for row in rows], dtype=float)


def max_abs_pull(params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    pulls = []
    for i, row in enumerate(rows):
        pred = predict_one(params, row.setting, row.variable)
        pulls.append((measured[i] - pred) / row.mu_err)
    return float(np.max(np.abs(pulls)))


def write_central_residuals(rows: list[ObsRow], params: np.ndarray, measured: np.ndarray) -> None:
    with RESIDUALS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "dataset",
                "method",
                "init",
                "setting",
                "shms_p_GeV",
                "delta_center",
                "dp_lo",
                "dp_hi",
                "variable",
                "mu",
                "mu_err",
                "sigma",
                "sigma_err",
                "predicted",
                "residual",
                "pull",
            ]
        )
        for i, row in enumerate(rows):
            pred = predict_one(params, row.setting, row.variable)
            resid = measured[i] - pred
            pull = resid / row.mu_err
            writer.writerow(
                [
                    "radwan_input",
                    "scipy_nelder_mead_bootstrap_central",
                    "zero_init",
                    row.setting,
                    f"{row.shms_p_GeV:.6f}",
                    f"{row.delta_center:.6f}",
                    f"{row.dp_lo:.6f}",
                    f"{row.dp_hi:.6f}",
                    row.variable,
                    f"{measured[i]:.6f}",
                    f"{row.mu_err:.6f}",
                    f"{row.sigma:.6f}",
                    f"{row.sigma_err:.6f}",
                    f"{pred:.6f}",
                    f"{resid:.6f}",
                    f"{pull:.6f}",
                ]
            )


def percentile(x: np.ndarray, q: float) -> float:
    return float(np.percentile(x, q)) if x.size else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-toys", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=314159)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_input()
    measured = vector_mu(rows)
    errors = vector_err(rows)

    central_result = fit_combined(rows, measured)
    central_params = central_result.x
    central_chi2 = chi2(central_params, rows, measured)
    central_exact = exact_inversion_like(central_params, central_chi2)

    rng = np.random.default_rng(args.seed)
    toy_params: list[np.ndarray] = []
    toy_rows: list[list[object]] = []

    for itoy in range(args.n_toys):
        toy_measured = rng.normal(measured, errors)
        result = fit_combined(rows, toy_measured)
        params = result.x
        fit_chi2 = chi2(params, rows, toy_measured)
        exact = exact_inversion_like(params, fit_chi2)
        good_physical = bool(result.success and not exact)
        toy_params.append(params)
        toy_rows.append(
            [
                itoy,
                bool(result.success),
                str(result.message),
                int(result.nit),
                int(result.nfev),
                f"{fit_chi2:.6f}",
                f"{max_abs_pull(params, rows, toy_measured):.6f}",
                int(exact),
                int(good_physical),
                *[f"{x:.6f}" for x in params],
            ]
        )

    toy_array = np.array(toy_params, dtype=float)
    exact_flags = np.array([int(r[7]) for r in toy_rows], dtype=int)
    success_flags = np.array([bool(r[1]) for r in toy_rows], dtype=bool)
    physical_mask = (success_flags) & (exact_flags == 0)

    with TOYS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "toy",
                "success",
                "message",
                "nit",
                "nfev",
                "chi2",
                "max_abs_pull",
                "exact_inversion_flag",
                "physical_branch_flag",
                *PARAM_NAMES,
            ]
        )
        writer.writerows(toy_rows)

    with PARAMS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "dataset",
                "method",
                "init",
                "parameter",
                "central_value",
                "toy_mean_physical",
                "toy_std_physical",
                "toy_p16_physical",
                "toy_p50_physical",
                "toy_p84_physical",
                "n_toys_physical",
            ]
        )
        physical = toy_array[physical_mask]
        for ipar, name in enumerate(PARAM_NAMES):
            vals = physical[:, ipar] if physical.size else np.array([], dtype=float)
            writer.writerow(
                [
                    "radwan_input",
                    "scipy_nelder_mead_bootstrap",
                    "zero_init",
                    name,
                    f"{central_params[ipar]:.6f}",
                    f"{float(np.mean(vals)):.6f}" if vals.size else "nan",
                    f"{float(np.std(vals, ddof=1)):.6f}" if vals.size > 1 else "nan",
                    f"{percentile(vals, 16):.6f}",
                    f"{percentile(vals, 50):.6f}",
                    f"{percentile(vals, 84):.6f}",
                    int(vals.size),
                ]
            )

    with SUMMARY_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        writer.writerow(["dataset", "radwan_input"])
        writer.writerow(["method", "scipy_nelder_mead_bootstrap"])
        writer.writerow(["init", "zero_init"])
        writer.writerow(["n_toys_requested", args.n_toys])
        writer.writerow(["seed", args.seed])
        writer.writerow(["central_success", bool(central_result.success)])
        writer.writerow(["central_message", str(central_result.message)])
        writer.writerow(["central_chi2", f"{central_chi2:.6f}"])
        writer.writerow(["central_n_obs", len(rows)])
        writer.writerow(["central_n_params", len(PARAM_NAMES)])
        writer.writerow(["central_ndf", len(rows) - len(PARAM_NAMES)])
        writer.writerow(["central_chi2_ndf", f"{central_chi2 / (len(rows) - len(PARAM_NAMES)):.6f}"])
        writer.writerow(["central_max_abs_pull", f"{max_abs_pull(central_params, rows, measured):.6f}"])
        writer.writerow(["central_exact_inversion_flag", int(central_exact)])
        writer.writerow(["toy_success_count", int(np.sum(success_flags))])
        writer.writerow(["toy_exact_inversion_count", int(np.sum(exact_flags))])
        writer.writerow(["toy_physical_branch_count", int(np.sum(physical_mask))])
        writer.writerow(["toy_physical_branch_fraction", f"{float(np.mean(physical_mask)):.6f}"])

    write_central_residuals(rows, central_params, measured)

    print(f"[INFO] Wrote {PARAMS_CSV}")
    print(f"[INFO] Wrote {TOYS_CSV}")
    print(f"[INFO] Wrote {RESIDUALS_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")
    print(f"[INFO] Physical toy fraction: {float(np.mean(physical_mask)):.4f}")


if __name__ == "__main__":
    main()
