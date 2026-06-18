#!/usr/bin/env python3
"""
ScipyWithLeastSquares.py

Formal least-squares covariance diagnostic for the 5-pass combined fit using
Radwan's measured offsets.

Run from heep_check_v3:
  python3 macros/ScipyWithLeastSquares.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


INPUT_CSV = Path("results/tables/Julio/TestingRadwan_input.csv")
OUT_DIR = Path("results/tables/ScipyWithLeastSquares")

PARAMS_CSV = OUT_DIR / "ScipyWithLeastSquares_params.csv"
RESIDUALS_CSV = OUT_DIR / "ScipyWithLeastSquares_residuals.csv"
COVARIANCE_CSV = OUT_DIR / "ScipyWithLeastSquares_covariance.csv"
SUMMARY_CSV = OUT_DIR / "ScipyWithLeastSquares_summary.csv"

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


def read_input() -> list[ObsRow]:
    raw: dict[tuple[str, str], ObsRow] = {}
    with INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            setting = row["setting"]
            variable = row["variable"]
            if setting not in SETTINGS or variable not in OBSERVABLES:
                continue
            raw[(setting, variable)] = ObsRow(
                setting=setting,
                shms_p_GeV=float(row["shms_p_GeV"]),
                delta_center=float(row["delta_center"]),
                dp_lo=float(row["dp_lo"]),
                dp_hi=float(row["dp_hi"]),
                variable=variable,
                mu=float(row["mu"]),
                mu_err=float(row["mu_err"]),
                sigma=float(row["sigma"]),
                sigma_err=float(row["sigma_err"]),
            )

    rows: list[ObsRow] = []
    missing = []
    for setting in SETTINGS:
        for variable in OBSERVABLES:
            key = (setting, variable)
            if key not in raw:
                missing.append(key)
            else:
                rows.append(raw[key])
    if missing:
        raise RuntimeError(f"Missing rows in {INPUT_CSV}: {missing}")
    return rows


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


def predict_one(params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(RESPONSE[ivar] @ parameter_vector_for_setting(params, row.setting))


def residuals(params: np.ndarray, rows: list[ObsRow]) -> np.ndarray:
    return np.array([(row.mu - predict_one(params, row)) / row.mu_err for row in rows], dtype=float)


def exact_inversion_like(params: np.ndarray, chi2: float) -> bool:
    return bool(np.max(np.abs(params)) > 20.0 or (chi2 < 1.0e-6 and np.max(np.abs(params)) > 5.0))


def covariance_from_jacobian(jac: np.ndarray, chi2: float, ndf: int) -> tuple[np.ndarray, str, float, int]:
    jtj = jac.T @ jac
    singular_values = np.linalg.svd(jtj, compute_uv=False)
    rank = int(np.linalg.matrix_rank(jtj))
    condition = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else float("inf")
    scale = chi2 / ndf if ndf > 0 else 1.0

    try:
        cov = np.linalg.inv(jtj) * scale
        method = "inverse"
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(jtj) * scale
        method = "pinv"

    if not np.all(np.isfinite(cov)):
        cov = np.linalg.pinv(jtj) * scale
        method = "pinv_nonfinite_fallback"

    return cov, method, condition, rank


def write_outputs(rows: list[ObsRow], result) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    params = result.x
    pulls = residuals(params, rows)
    chi2 = float(np.sum(pulls * pulls))
    n_obs = len(rows)
    n_params = len(params)
    ndf = n_obs - n_params
    cov, cov_method, jtj_condition, jtj_rank = covariance_from_jacobian(result.jac, chi2, ndf)
    errs = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
    exact_flag = exact_inversion_like(params, chi2)

    with PARAMS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "dataset",
                "method",
                "init",
                "parameter",
                "value",
                "formal_error",
                "covariance_method",
            ]
        )
        for name, value, err in zip(PARAM_NAMES, params, errs):
            writer.writerow(
                [
                    "radwan_input",
                    "scipy_least_squares",
                    "zero_init",
                    name,
                    f"{value:.6f}",
                    f"{err:.6f}",
                    cov_method,
                ]
            )

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
        for row, pull in zip(rows, pulls):
            pred = predict_one(params, row)
            resid = row.mu - pred
            writer.writerow(
                [
                    "radwan_input",
                    "scipy_least_squares",
                    "zero_init",
                    row.setting,
                    f"{row.shms_p_GeV:.6f}",
                    f"{row.delta_center:.6f}",
                    f"{row.dp_lo:.6f}",
                    f"{row.dp_hi:.6f}",
                    row.variable,
                    f"{row.mu:.6f}",
                    f"{row.mu_err:.6f}",
                    f"{row.sigma:.6f}",
                    f"{row.sigma_err:.6f}",
                    f"{pred:.6f}",
                    f"{resid:.6f}",
                    f"{pull:.6f}",
                ]
            )

    with COVARIANCE_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["row_parameter", *PARAM_NAMES])
        for name, cov_row in zip(PARAM_NAMES, cov):
            writer.writerow([name, *[f"{x:.10g}" for x in cov_row]])

    with SUMMARY_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "value"])
        writer.writerow(["dataset", "radwan_input"])
        writer.writerow(["method", "scipy_least_squares"])
        writer.writerow(["init", "zero_init"])
        writer.writerow(["success", bool(result.success)])
        writer.writerow(["status", int(result.status)])
        writer.writerow(["message", str(result.message)])
        writer.writerow(["nfev", int(result.nfev)])
        writer.writerow(["njev", int(result.njev) if result.njev is not None else ""])
        writer.writerow(["chi2", f"{chi2:.6f}"])
        writer.writerow(["n_obs", n_obs])
        writer.writerow(["n_params", n_params])
        writer.writerow(["ndf", ndf])
        writer.writerow(["chi2_ndf", f"{chi2 / ndf:.6f}" if ndf > 0 else "nan"])
        writer.writerow(["max_abs_pull", f"{float(np.max(np.abs(pulls))):.6f}"])
        writer.writerow(["exact_inversion_flag", int(exact_flag)])
        writer.writerow(["covariance_method", cov_method])
        writer.writerow(["jtj_rank", jtj_rank])
        writer.writerow(["jtj_condition_number", f"{jtj_condition:.6g}"])


def main() -> None:
    rows = read_input()
    result = least_squares(lambda p: residuals(p, rows), ZERO_INIT, method="trf")
    write_outputs(rows, result)
    print(f"[INFO] Wrote {PARAMS_CSV}")
    print(f"[INFO] Wrote {RESIDUALS_CSV}")
    print(f"[INFO] Wrote {COVARIANCE_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")


if __name__ == "__main__":
    main()
