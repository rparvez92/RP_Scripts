#!/usr/bin/env python3
"""
TestingJulio.py

Single-setting 5-pass response-matrix test using Julio's quoted measured
observable offsets. This script writes a small input CSV, runs scipy
Nelder-Mead fits independently for settings C, D, and E, and compares those
fits to Julio's reported notebook parameters.

It also builds the same input/output schema for Radwan's measured offsets from:
  results/tables/MeasuredOffsetsBySetting_5pass.csv

Run from RSIDIS/heep_check_v3:
  python3 macros/TestingJulio.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


INPUT_CSV = Path("results/tables/Julio/TestingJulio_input.csv")
OUTPUT_CSV = Path("results/tables/Julio/TestingJulio_output.csv")
RADWAN_SOURCE_CSV = Path("results/tables/MeasuredOffsetsBySetting_5pass.csv")
RADWAN_INPUT_CSV = Path("results/tables/Julio/TestingRadwan_input.csv")
RADWAN_OUTPUT_CSV = Path("results/tables/Julio/TestingRadwan_output.csv")

OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
PARAM_NAMES = ("dtheta_e", "dpe", "dtheta_p", "dpp")

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

JULIO_INPUT_ROWS = [
    # setting, shms_p_GeV, delta_center, dp_lo, dp_hi, variable, mu, mu_err, sigma, sigma_err
    ("C", -6.7200, 1.0, -0.5, 2.5, "W", 10.468058, 1.002105, 37.972155, 0.893049),
    ("C", -6.7200, 1.0, -0.5, 2.5, "Em", -8.474547, 0.186758, 7.133772, 0.172972),
    ("C", -6.7200, 1.0, -0.5, 2.5, "Pmz", 11.293317, 0.559726, 9.323317, 0.662723),
    ("C", -6.7200, 1.0, -0.5, 2.5, "Pmy", 12.310709, 0.434572, 11.728930, 0.424420),
    ("E", -6.6048, 3.0, 1.5, 4.5, "W", 6.229827, 1.096396, 35.317492, 0.958426),
    ("E", -6.6048, 3.0, 1.5, 4.5, "Em", -15.954901, 0.261251, 7.348766, 0.247712),
    ("E", -6.6048, 3.0, 1.5, 4.5, "Pmz", 16.510733, 0.513753, 9.304224, 0.554649),
    ("E", -6.6048, 3.0, 1.5, 4.5, "Pmy", 9.259874, 0.547192, 13.274062, 0.534578),
    ("D", -6.3840, 6.0, 4.5, 7.5, "W", 13.607670, 1.312226, 34.572211, 1.243736),
    ("D", -6.3840, 6.0, 4.5, 7.5, "Em", -14.401973, 0.400741, 9.055663, 0.373045),
    ("D", -6.3840, 6.0, 4.5, 7.5, "Pmz", 16.704784, 0.298057, 8.313920, 0.278653),
    ("D", -6.3840, 6.0, 4.5, 7.5, "Pmy", 8.709048, 0.431981, 11.820399, 0.388189),
]

JULIO_REPORTED_PARAMS = {
    "C": np.array([-0.8619, 0.9869, -2.4567, 0.4436], dtype=float),
    "E": np.array([-0.6471, 1.0048, -1.5902, 1.9413], dtype=float),
    "D": np.array([-1.1178, 1.4619, -1.3684, 0.9741], dtype=float),
}


def write_julio_input_csv() -> None:
    INPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with INPUT_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
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
            ]
        )
        writer.writerows(JULIO_INPUT_ROWS)


def write_radwan_input_csv() -> None:
    expected_ranges = {
        "C": (-0.5, 2.5),
        "E": (1.5, 4.5),
        "D": (4.5, 7.5),
    }

    rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    with RADWAN_SOURCE_CSV.open() as f:
        for row in csv.DictReader(f):
            setting = row["setting"]
            variable = row["var"]
            if setting not in expected_ranges:
                continue
            if row["dp_label"] != "b1":
                continue
            if row["fit_valid"] != "1":
                continue
            if variable not in OBSERVABLES:
                continue

            dp_lo = float(row["dp_lo"])
            dp_hi = float(row["dp_hi"])
            exp_lo, exp_hi = expected_ranges[setting]
            if abs(dp_lo - exp_lo) > 1e-9 or abs(dp_hi - exp_hi) > 1e-9:
                raise RuntimeError(
                    f"Unexpected {setting} b1 range: found [{dp_lo}, {dp_hi}], "
                    f"expected [{exp_lo}, {exp_hi}]"
                )

            key = (setting, variable)
            if key in rows_by_key:
                raise RuntimeError(f"Duplicate measured row for {setting} {variable}")
            rows_by_key[key] = row

    missing = [
        (setting, variable)
        for setting in ("C", "E", "D")
        for variable in OBSERVABLES
        if (setting, variable) not in rows_by_key
    ]
    if missing:
        raise RuntimeError(f"Missing Radwan measured rows: {missing}")

    RADWAN_INPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with RADWAN_INPUT_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
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
            ]
        )
        for setting in ("C", "E", "D"):
            for variable in OBSERVABLES:
                row = rows_by_key[(setting, variable)]
                dp_lo = float(row["dp_lo"])
                dp_hi = float(row["dp_hi"])
                writer.writerow(
                    [
                        setting,
                        float(row["shms_p_GeV"]),
                        0.5 * (dp_lo + dp_hi),
                        dp_lo,
                        dp_hi,
                        variable,
                        float(row["offset_MeV"]),
                        float(row["offset_err_MeV"]),
                        float(row["sigD_MeV"]),
                        float(row["sigD_err_MeV"]),
                    ]
                )


def read_input_by_setting(input_csv: Path) -> dict[str, dict[str, object]]:
    by_setting: dict[str, dict[str, object]] = {}
    with input_csv.open() as f:
        for row in csv.DictReader(f):
            setting = row["setting"]
            if setting not in by_setting:
                by_setting[setting] = {
                    "shms_p_GeV": float(row["shms_p_GeV"]),
                    "delta_center": float(row["delta_center"]),
                    "dp_lo": float(row["dp_lo"]),
                    "dp_hi": float(row["dp_hi"]),
                    "mu": {},
                    "mu_err": {},
                    "sigma": {},
                    "sigma_err": {},
                }
            variable = row["variable"]
            by_setting[setting]["mu"][variable] = float(row["mu"])  # type: ignore[index]
            by_setting[setting]["mu_err"][variable] = float(row["mu_err"])  # type: ignore[index]
            by_setting[setting]["sigma"][variable] = float(row["sigma"])  # type: ignore[index]
            by_setting[setting]["sigma_err"][variable] = float(row["sigma_err"])  # type: ignore[index]
    return by_setting


def vector_from_setting(setting_data: dict[str, object], key: str) -> np.ndarray:
    values = setting_data[key]  # type: ignore[index]
    return np.array([values[var] for var in OBSERVABLES], dtype=float)  # type: ignore[index]


def chi2(params: np.ndarray, measured: np.ndarray, errors: np.ndarray) -> float:
    predicted = RESPONSE @ params
    return float(np.sum(((measured - predicted) / errors) ** 2))


def fit_scipy_nelder_mead(measured: np.ndarray, errors: np.ndarray):
    return minimize(
        lambda p: chi2(p, measured, errors),
        np.zeros(4, dtype=float),
        method="Nelder-Mead",
    )


def append_output_rows(
    writer: csv.writer,
    fit_type: str,
    setting: str,
    setting_data: dict[str, object],
    params: np.ndarray,
    success: bool | str,
    message: str,
    nit: int | str,
    nfev: int | str,
) -> None:
    measured = vector_from_setting(setting_data, "mu")
    errors = vector_from_setting(setting_data, "mu_err")
    sigmas = vector_from_setting(setting_data, "sigma")
    sigma_errs = vector_from_setting(setting_data, "sigma_err")
    predicted = RESPONSE @ params
    residual = measured - predicted
    fit_chi2 = chi2(params, measured, errors)

    for i, variable in enumerate(OBSERVABLES):
        writer.writerow(
            [
                fit_type,
                setting,
                setting_data["shms_p_GeV"],
                setting_data["delta_center"],
                setting_data["dp_lo"],
                setting_data["dp_hi"],
                variable,
                measured[i],
                errors[i],
                sigmas[i],
                sigma_errs[i],
                predicted[i],
                residual[i],
                residual[i] / errors[i],
                params[0],
                params[1],
                params[2],
                params[3],
                fit_chi2,
                4,
                4,
                0,
                success,
                message,
                nit,
                nfev,
            ]
        )


def write_output_csv(
    by_setting: dict[str, dict[str, object]],
    output_csv: Path,
    include_julio_reference: bool,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "fit_type",
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
            ]
        )

        for setting in ("C", "E", "D"):
            setting_data = by_setting[setting]
            measured = vector_from_setting(setting_data, "mu")
            errors = vector_from_setting(setting_data, "mu_err")
            result = fit_scipy_nelder_mead(measured, errors)

            append_output_rows(
                writer,
                "scipy_best_fit",
                setting,
                setting_data,
                result.x,
                bool(result.success),
                str(result.message),
                int(result.nit),
                int(result.nfev),
            )

            if include_julio_reference:
                append_output_rows(
                    writer,
                    "julio_reported_reference",
                    setting,
                    setting_data,
                    JULIO_REPORTED_PARAMS[setting],
                    "reference",
                    "copied from Julio notebook; not refit here",
                    "",
                    "",
                )


def print_summary(output_csv: Path, label: str) -> None:
    singular_values = np.linalg.svd(RESPONSE, compute_uv=False)
    print(f"\n[INFO] Summary for {label}")
    print(f"[INFO] Wrote output: {output_csv}")
    print(f"[INFO] Response condition number: {singular_values[0] / singular_values[-1]:.6g}")

    with output_csv.open() as f:
        rows = list(csv.DictReader(f))

    fit_types = []
    for fit_type in ("scipy_best_fit", "julio_reported_reference"):
        if any(row["fit_type"] == fit_type for row in rows):
            fit_types.append(fit_type)

    for fit_type in fit_types:
        print(f"\n{fit_type}")
        for setting in ("C", "E", "D"):
            setting_rows = [r for r in rows if r["fit_type"] == fit_type and r["setting"] == setting]
            p = setting_rows[0]
            max_abs_pull = max(abs(float(r["pull"])) for r in setting_rows)
            print(
                f"  {setting}: "
                f"dtheta_e={float(p['dtheta_e']): .6f}, "
                f"dpe={float(p['dpe']): .6f}, "
                f"dtheta_p={float(p['dtheta_p']): .6f}, "
                f"dpp={float(p['dpp']): .6f}, "
                f"chi2={float(p['chi2']):.6f}, "
                f"max|pull|={max_abs_pull:.3f}"
            )


def main() -> None:
    write_julio_input_csv()
    julio_by_setting = read_input_by_setting(INPUT_CSV)
    write_output_csv(julio_by_setting, OUTPUT_CSV, include_julio_reference=True)
    print(f"[INFO] Wrote input:  {INPUT_CSV}")
    print_summary(OUTPUT_CSV, "Julio")

    write_radwan_input_csv()
    radwan_by_setting = read_input_by_setting(RADWAN_INPUT_CSV)
    write_output_csv(radwan_by_setting, RADWAN_OUTPUT_CSV, include_julio_reference=False)
    print(f"\n[INFO] Wrote input:  {RADWAN_INPUT_CSV}")
    print_summary(RADWAN_OUTPUT_CSV, "Radwan")


if __name__ == "__main__":
    main()
