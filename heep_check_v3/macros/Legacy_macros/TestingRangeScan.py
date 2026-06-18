#!/usr/bin/env python3
"""
TestingRangeScan.py

SciPy response-matrix fits for the RangeScan measured offsets.

Run from heep_check_v3:
  python3 macros/TestingRangeScan.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


INPUT_CSV = Path("results/tables/RangeScan/TestingRadwan_RangeScan_input.csv")
OUTPUT_CSV = Path("results/tables/RangeScan/TestingRadwan_RangeScan_SciPy_output.csv")
SUMMARY_CSV = Path("results/tables/RangeScan/TestingRadwan_RangeScan_Summary.csv")

OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
SETTINGS = ("C", "E", "D")

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


def chi2(params: np.ndarray, measured: np.ndarray, errors: np.ndarray) -> float:
    predicted = RESPONSE @ params
    return float(np.sum(((measured - predicted) / errors) ** 2))


def read_input() -> dict[tuple[float, str], dict[str, object]]:
    grouped: dict[tuple[float, str], dict[str, object]] = defaultdict(
        lambda: {
            "rows": {},
            "shms_p_GeV": np.nan,
            "delta_center": np.nan,
            "dp_lo": np.nan,
            "dp_hi": np.nan,
        }
    )
    with INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            if row.get("fit_valid", "0") != "1":
                continue
            k_sigma = float(row["fit_k_sigma"])
            setting = row["setting"]
            variable = row["variable"]
            if setting not in SETTINGS or variable not in OBSERVABLES:
                continue
            key = (k_sigma, setting)
            grouped[key]["shms_p_GeV"] = float(row["shms_p_GeV"])
            grouped[key]["delta_center"] = float(row.get("bin_center", row.get("delta_center", "nan")))
            grouped[key]["dp_lo"] = float(row["dp_lo"])
            grouped[key]["dp_hi"] = float(row["dp_hi"])
            grouped[key]["rows"][variable] = row  # type: ignore[index]

    incomplete = []
    for key, data in list(grouped.items()):
        rows = data["rows"]  # type: ignore[index]
        missing = [variable for variable in OBSERVABLES if variable not in rows]
        if missing:
            incomplete.append((key[0], key[1], missing))
            del grouped[key]
    if incomplete:
        print("[WARN] Skipping incomplete RangeScan groups:")
        for k_sigma, setting, missing in incomplete:
            print(f"  k={k_sigma:.2f} setting={setting}: missing {missing}")
    return dict(grouped)


def vector_for(data: dict[str, object], column: str) -> np.ndarray:
    rows = data["rows"]  # type: ignore[index]
    return np.array([float(rows[var][column]) for var in OBSERVABLES], dtype=float)  # type: ignore[index]


def fit_nelder_mead(measured: np.ndarray, errors: np.ndarray):
    return minimize(
        lambda p: chi2(p, measured, errors),
        np.zeros(4, dtype=float),
        method="Nelder-Mead",
    )


def is_exact_inversion_like(params: np.ndarray, fit_chi2: float) -> bool:
    return bool(np.max(np.abs(params)) > 20.0 or (fit_chi2 < 1.0e-6 and np.max(np.abs(params)) > 5.0))


def main() -> None:
    grouped = read_input()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_CSV.open("w", newline="") as f_out, SUMMARY_CSV.open("w", newline="") as f_sum:
        writer = csv.writer(f_out)
        summary = csv.writer(f_sum)

        writer.writerow(
            [
                "fit_type",
                "fit_k_sigma",
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
                "dtheta_e",
                "dpe",
                "dtheta_p",
                "dpp",
                "chi2",
                "n_obs",
                "n_params",
                "ndf",
                "success",
                "message",
                "nit",
                "nfev",
                "exact_inversion_flag",
            ]
        )
        summary.writerow(
            [
                "method",
                "fit_k_sigma",
                "setting",
                "dtheta_e",
                "dpe",
                "dtheta_p",
                "dpp",
                "chi2",
                "max_abs_pull",
                "exact_inversion_flag",
                "success",
            ]
        )

        for k_sigma in sorted({key[0] for key in grouped}):
            for setting in SETTINGS:
                key = (k_sigma, setting)
                if key not in grouped:
                    continue
                data = grouped[key]
                measured = vector_for(data, "mu")
                errors = vector_for(data, "mu_err")
                sigmas = vector_for(data, "sigma")
                sigma_errs = vector_for(data, "sigma_err")

                result = fit_nelder_mead(measured, errors)
                params = result.x
                predicted = RESPONSE @ params
                residual = measured - predicted
                pulls = residual / errors
                fit_chi2 = chi2(params, measured, errors)
                exact_flag = is_exact_inversion_like(params, fit_chi2)

                summary.writerow(
                    [
                        "scipy_nelder_mead",
                        f"{k_sigma:.2f}",
                        setting,
                        f"{params[0]:.6f}",
                        f"{params[1]:.6f}",
                        f"{params[2]:.6f}",
                        f"{params[3]:.6f}",
                        f"{fit_chi2:.6f}",
                        f"{np.max(np.abs(pulls)):.6f}",
                        int(exact_flag),
                        bool(result.success),
                    ]
                )

                for i, variable in enumerate(OBSERVABLES):
                    writer.writerow(
                        [
                            "scipy_nelder_mead",
                            f"{k_sigma:.2f}",
                            setting,
                            data["shms_p_GeV"],
                            data["delta_center"],
                            data["dp_lo"],
                            data["dp_hi"],
                            variable,
                            f"{measured[i]:.6f}",
                            f"{errors[i]:.6f}",
                            f"{sigmas[i]:.6f}",
                            f"{sigma_errs[i]:.6f}",
                            f"{predicted[i]:.6f}",
                            f"{residual[i]:.6f}",
                            f"{pulls[i]:.6f}",
                            f"{params[0]:.6f}",
                            f"{params[1]:.6f}",
                            f"{params[2]:.6f}",
                            f"{params[3]:.6f}",
                            f"{fit_chi2:.6f}",
                            4,
                            4,
                            0,
                            bool(result.success),
                            str(result.message),
                            int(result.nit),
                            int(result.nfev),
                            int(exact_flag),
                        ]
                    )

    svals = np.linalg.svd(RESPONSE, compute_uv=False)
    print(f"[INFO] Wrote {OUTPUT_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")
    print(f"[INFO] Response condition number: {svals[0] / svals[-1]:.6g}")


if __name__ == "__main__":
    main()
