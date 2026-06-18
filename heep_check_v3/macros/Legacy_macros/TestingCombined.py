#!/usr/bin/env python3
"""
TestingCombined.py

SciPy combined-mode response-matrix test for 5-pass C/E/D.

Combined mode follows Julio's notebook structure:
  - dpe is setting dependent:        dpe_C, dpe_E, dpe_D
  - dtheta_e is setting dependent:   dtheta_e_C, dtheta_e_E, dtheta_e_D
  - dtheta_p is common across C/E/D
  - dpp is common across C/E/D

Run from heep_check_v3:
  python3 macros/TestingCombined.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


OUT_DIR = Path("results/tables/CombinedTest")
PARAMS_CSV = OUT_DIR / "CombinedTest_params.csv"
RESIDUALS_CSV = OUT_DIR / "CombinedTest_residuals.csv"
SUMMARY_CSV = OUT_DIR / "CombinedTest_summary.csv"

JULIO_INPUT = Path("results/tables/Julio/TestingJulio_input.csv")
RADWAN_INPUT = Path("results/tables/Julio/TestingRadwan_input.csv")
RANGESCAN_INPUT = Path("results/tables/RangeScan/TestingRadwan_RangeScan_input.csv")

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
JULIO_LIKE_INIT = np.array(
    [
        0.9869,
        1.0048,
        1.4619,
        -0.8619,
        -0.6471,
        -1.1178,
        0.0,
        0.0,
    ],
    dtype=float,
)


@dataclass(frozen=True)
class ObsRow:
    dataset: str
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


def param_vector_for_setting(params: np.ndarray, setting: str) -> np.ndarray:
    iset = SETTINGS.index(setting)
    dpe = params[iset]
    dtheta_e = params[3 + iset]
    dtheta_p = params[6]
    dpp = params[7]
    return np.array([dtheta_e, dpe, dtheta_p, dpp], dtype=float)


def prediction(params: np.ndarray, setting: str) -> np.ndarray:
    return RESPONSE @ param_vector_for_setting(params, setting)


def chi2(params: np.ndarray, rows: list[ObsRow]) -> float:
    total = 0.0
    for row in rows:
        ivar = OBSERVABLES.index(row.variable)
        pred = prediction(params, row.setting)[ivar]
        total += ((row.mu - pred) / row.mu_err) ** 2
    return float(total)


def read_standard_input(path: Path, dataset: str) -> list[ObsRow]:
    rows: list[ObsRow] = []
    with path.open() as f:
        for row in csv.DictReader(f):
            if row["setting"] not in SETTINGS or row["variable"] not in OBSERVABLES:
                continue
            rows.append(
                ObsRow(
                    dataset=dataset,
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
    validate_rows(dataset, rows)
    return rows


def read_rangescan_input(path: Path, k_sigma: float, dataset: str) -> list[ObsRow] | None:
    rows: list[ObsRow] = []
    with path.open() as f:
        for row in csv.DictReader(f):
            if abs(float(row["fit_k_sigma"]) - k_sigma) > 1e-9:
                continue
            if row.get("fit_valid", "0") != "1":
                continue
            if row["setting"] not in SETTINGS or row["variable"] not in OBSERVABLES:
                continue
            rows.append(
                ObsRow(
                    dataset=dataset,
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

    try:
        validate_rows(dataset, rows)
    except RuntimeError as exc:
        print(f"[WARN] Skipping {dataset}: {exc}")
        return None
    return rows


def validate_rows(dataset: str, rows: list[ObsRow]) -> None:
    present = {(r.setting, r.variable) for r in rows}
    missing = [
        (setting, variable)
        for setting in SETTINGS
        for variable in OBSERVABLES
        if (setting, variable) not in present
    ]
    if missing:
        raise RuntimeError(f"{dataset} missing rows: {missing}")
    if len(rows) != len(SETTINGS) * len(OBSERVABLES):
        raise RuntimeError(f"{dataset} has {len(rows)} rows, expected 12")


def fit_combined(rows: list[ObsRow], init: np.ndarray):
    return minimize(
        lambda p: chi2(p, rows),
        init,
        method="Nelder-Mead",
        options={"maxiter": 20000, "maxfev": 30000, "xatol": 1e-10, "fatol": 1e-10},
    )


def exact_inversion_like(params: np.ndarray, fit_chi2: float) -> bool:
    return bool(np.max(np.abs(params)) > 20.0 or (fit_chi2 < 1.0e-6 and np.max(np.abs(params)) > 5.0))


def write_result(
    params_writer: csv.writer,
    residuals_writer: csv.writer,
    summary_writer: csv.writer,
    dataset: str,
    init_name: str,
    rows: list[ObsRow],
    result,
) -> None:
    params = result.x
    fit_chi2 = chi2(params, rows)
    n_obs = len(rows)
    n_params = len(params)
    ndf = n_obs - n_params

    pulls: list[float] = []
    for row in rows:
        ivar = OBSERVABLES.index(row.variable)
        pred = prediction(params, row.setting)[ivar]
        resid = row.mu - pred
        pull = resid / row.mu_err
        pulls.append(float(pull))
        residuals_writer.writerow(
            [
                dataset,
                "scipy_nelder_mead",
                init_name,
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

    exact_flag = exact_inversion_like(params, fit_chi2)
    max_abs_pull = max(abs(p) for p in pulls)

    for name, value in zip(PARAM_NAMES, params):
        params_writer.writerow(
            [
                dataset,
                "scipy_nelder_mead",
                init_name,
                name,
                f"{value:.6f}",
                "nan",
            ]
        )

    summary_writer.writerow(
        [
            dataset,
            "scipy_nelder_mead",
            init_name,
            f"{fit_chi2:.6f}",
            n_obs,
            n_params,
            ndf,
            f"{fit_chi2 / ndf:.6f}" if ndf > 0 else "nan",
            f"{max_abs_pull:.6f}",
            int(exact_flag),
            bool(result.success),
            str(result.message),
            int(result.nit),
            int(result.nfev),
        ]
    )


def build_datasets() -> dict[str, list[ObsRow]]:
    datasets = {
        "julio_input": read_standard_input(JULIO_INPUT, "julio_input"),
        "radwan_input": read_standard_input(RADWAN_INPUT, "radwan_input"),
    }
    for k_sigma in (1.5, 2.0):
        name = f"radwan_rangescan_k{str(k_sigma).replace('.', 'p')}"
        rows = read_rangescan_input(RANGESCAN_INPUT, k_sigma, name)
        if rows is not None:
            datasets[name] = rows
    return datasets


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    datasets = build_datasets()

    with PARAMS_CSV.open("w", newline="") as f_params, RESIDUALS_CSV.open("w", newline="") as f_res, SUMMARY_CSV.open("w", newline="") as f_sum:
        params_writer = csv.writer(f_params)
        residuals_writer = csv.writer(f_res)
        summary_writer = csv.writer(f_sum)

        params_writer.writerow(["dataset", "method", "init", "parameter", "value", "error_or_nan"])
        residuals_writer.writerow(
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
        summary_writer.writerow(
            [
                "dataset",
                "method",
                "init",
                "chi2",
                "n_obs",
                "n_params",
                "ndf",
                "chi2_ndf",
                "max_abs_pull",
                "exact_inversion_flag",
                "success",
                "message",
                "nit",
                "nfev",
            ]
        )

        for dataset, rows in datasets.items():
            for init_name, init in (("zero_init", ZERO_INIT), ("julio_like_init", JULIO_LIKE_INIT)):
                result = fit_combined(rows, init)
                write_result(params_writer, residuals_writer, summary_writer, dataset, init_name, rows, result)

    svals = np.linalg.svd(RESPONSE, compute_uv=False)
    print(f"[INFO] Wrote {PARAMS_CSV}")
    print(f"[INFO] Wrote {RESIDUALS_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")
    print(f"[INFO] Response condition number: {svals[0] / svals[-1]:.6g}")


if __name__ == "__main__":
    main()
