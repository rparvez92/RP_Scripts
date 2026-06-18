#!/usr/bin/env python3
"""
Dry4PassTest_Fitting.py

Dry 4-pass global fit:
  parameters = dthe, dpe, dthp, dpp
  settings A/B use setting-dependent response matrices

Run from heep_check_v3:
  python3 macros/Dry4PassTest_Fitting.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


INPUT_CSV = Path("results/tables/Dry4PassTest/MeasuredOffsets.csv")
OUTPUT_CSV = Path("results/tables/Dry4PassTest/FitParams.csv")
FIVEPASS_INPUT_CSV = Path("results/tables/MeasuredOffsetsBySetting_5pass.csv")
OUTPUT_CSV_5PASS_ERRORS = Path("results/tables/Dry4PassTest/FitParams_Using5passErrors.csv")
OUTPUT_CSV_DOWNSCALE_A = Path("results/tables/Dry4PassTest/FitParams_DownscaleAErrors.csv")
OUTPUT_CSV_SCALED_TO_5PASS = Path("results/tables/Dry4PassTest/FitParams_ScaledTo5passEntries.csv")

SETTINGS = ("A", "B")
OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
PARAM_NAMES = ("dthe", "dpe", "dthp", "dpp")

RESPONSE = {
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

ZERO_INIT = np.zeros(4, dtype=float)
OPT_OPTIONS = {"maxiter": 20000, "maxfev": 30000, "xatol": 1e-10, "fatol": 1e-10}


@dataclass(frozen=True)
class ObsRow:
    setting: str
    delta_center: float
    dp_lo: float
    dp_hi: float
    variable: str
    measured: float
    measured_err: float
    mu_d_err: float
    mu_s_err: float
    entries_d: float
    entries_s: float
    sigma: float
    sigma_err: float


def with_measured_err(row: ObsRow, measured_err: float) -> ObsRow:
    return ObsRow(
        setting=row.setting,
        delta_center=row.delta_center,
        dp_lo=row.dp_lo,
        dp_hi=row.dp_hi,
        variable=row.variable,
        measured=row.measured,
        measured_err=measured_err,
        mu_d_err=row.mu_d_err,
        mu_s_err=row.mu_s_err,
        entries_d=row.entries_d,
        entries_s=row.entries_s,
        sigma=row.sigma,
        sigma_err=row.sigma_err,
    )


def f4(x: float) -> str:
    return f"{x:.4f}" if np.isfinite(x) else "nan"


def read_input() -> list[ObsRow]:
    rows_by_key: dict[tuple[str, str], ObsRow] = {}
    with INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["setting"] not in SETTINGS or row["var"] not in OBSERVABLES:
                continue
            if row["dp_label"] != "b1" or row["fit_valid"] != "1":
                continue
            rows_by_key[(row["setting"], row["var"])] = ObsRow(
                setting=row["setting"],
                delta_center=float(row["delta_center"]),
                dp_lo=float(row["dp_lo"]),
                dp_hi=float(row["dp_hi"]),
                variable=row["var"],
                measured=float(row["offset_MeV"]),
                measured_err=float(row["offset_err_MeV"]),
                mu_d_err=float(row["muD_err_MeV"]),
                mu_s_err=float(row["muS_err_MeV"]),
                entries_d=float(row["entriesD"]),
                entries_s=float(row["entriesS"]),
                sigma=float(row["sigD_MeV"]),
                sigma_err=float(row["sigD_err_MeV"]),
            )

    rows: list[ObsRow] = []
    missing = []
    for setting in SETTINGS:
        for variable in OBSERVABLES:
            key = (setting, variable)
            if key not in rows_by_key:
                missing.append(key)
            else:
                rows.append(rows_by_key[key])
    if missing:
        raise RuntimeError(f"Missing rows in {INPUT_CSV}: {missing}")
    return rows


def read_5pass_reference_errors() -> dict[tuple[str, str], float]:
    """Diagnostic only: use 5-pass C errors for A, and D errors for B."""
    source_to_target = {"C": "A", "D": "B"}
    errors: dict[tuple[str, str], float] = {}
    with FIVEPASS_INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            source_setting = row["setting"]
            if source_setting not in source_to_target:
                continue
            if row["dp_label"] != "b1" or row["fit_valid"] != "1":
                continue
            variable = row["var"]
            if variable not in OBSERVABLES:
                continue
            target_setting = source_to_target[source_setting]
            errors[(target_setting, variable)] = float(row["offset_err_MeV"])

    missing = [
        (setting, variable)
        for setting in SETTINGS
        for variable in OBSERVABLES
        if (setting, variable) not in errors
    ]
    if missing:
        raise RuntimeError(f"Missing 5-pass reference errors in {FIVEPASS_INPUT_CSV}: {missing}")
    return errors


def read_5pass_reference_entries() -> tuple[float, float]:
    entries_d: list[float] = []
    entries_s: list[float] = []
    with FIVEPASS_INPUT_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["setting"] not in {"C", "E", "D"}:
                continue
            if row["dp_label"] != "b1" or row["fit_valid"] != "1":
                continue
            if row["var"] not in OBSERVABLES:
                continue
            entries_d.append(float(row["entriesD"]))
            entries_s.append(float(row["entriesS"]))
    if not entries_d or not entries_s:
        raise RuntimeError(f"No 5-pass b1 reference entries found in {FIVEPASS_INPUT_CSV}")
    return float(np.mean(entries_d)), float(np.mean(entries_s))


def predict_one(params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(RESPONSE[row.setting][ivar] @ params)


def chi2(params: np.ndarray, rows: list[ObsRow]) -> float:
    return float(
        np.sum(
            [
                ((row.measured - predict_one(params, row)) / row.measured_err) ** 2
                for row in rows
            ]
        )
    )


def write_fit(output_csv: Path, rows: list[ObsRow], error_source: str) -> tuple[np.ndarray, float]:
    result = minimize(lambda p: chi2(p, rows), ZERO_INIT, method="Nelder-Mead", options=OPT_OPTIONS)
    params = result.x
    fit_chi2 = chi2(params, rows)
    n_obs = len(rows)
    n_params = len(params)
    ndf = n_obs - n_params
    chi2_ndf = fit_chi2 / ndf if ndf > 0 else np.nan

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "setting",
                "delta_center",
                "dp_lo",
                "dp_hi",
                "variable",
                "measured",
                "measured_err",
                "error_source",
                "sigma",
                "sigma_err",
                "predicted",
                "residual",
                "pull",
                *PARAM_NAMES,
                "chi2",
                "n_obs",
                "n_params",
                "ndf",
                "chi2_ndf",
                "success",
                "message",
                "nit",
                "nfev",
            ]
        )
        for row in rows:
            pred = predict_one(params, row)
            resid = row.measured - pred
            pull = resid / row.measured_err
            writer.writerow(
                [
                    row.setting,
                    f4(row.delta_center),
                    f4(row.dp_lo),
                    f4(row.dp_hi),
                    row.variable,
                    f4(row.measured),
                    f4(row.measured_err),
                    error_source,
                    f4(row.sigma),
                    f4(row.sigma_err),
                    f4(pred),
                    f4(resid),
                    f4(pull),
                    *[f4(x) for x in params],
                    f4(fit_chi2),
                    n_obs,
                    n_params,
                    ndf,
                    f4(chi2_ndf),
                    bool(result.success),
                    str(result.message),
                    int(result.nit),
                    int(result.nfev),
                ]
            )

    print(f"[INFO] Wrote {output_csv}")
    print(f"[INFO] {error_source} params:", dict(zip(PARAM_NAMES, params)))
    print(f"[INFO] {error_source} chi2/ndf = {chi2_ndf:.4f}")
    return params, chi2_ndf


def main() -> None:
    rows = read_input()
    write_fit(OUTPUT_CSV, rows, "dry4pass_measured_errors")

    reference_errors = read_5pass_reference_errors()
    rows_with_5pass_errors = [
        with_measured_err(row, reference_errors[(row.setting, row.variable)]) for row in rows
    ]
    write_fit(OUTPUT_CSV_5PASS_ERRORS, rows_with_5pass_errors, "A_uses_5pass_C_errors_B_uses_5pass_D_errors")

    data_scale = np.sqrt(27.0)
    sim_scale = np.sqrt(2.8)
    rows_downscale_a: list[ObsRow] = []
    for row in rows:
        if row.setting != "A":
            rows_downscale_a.append(row)
            continue
        mu_d_err = row.mu_d_err * data_scale
        mu_s_err = row.mu_s_err * sim_scale
        scaled_err = float(np.hypot(mu_d_err, mu_s_err))
        rows_downscale_a.append(with_measured_err(row, scaled_err))
    write_fit(
        OUTPUT_CSV_DOWNSCALE_A,
        rows_downscale_a,
        "A_offset_errors_scaled_by_sqrt27_B_unchanged",
    )

    avg_entries_d_5pass, avg_entries_s_5pass = read_5pass_reference_entries()
    rows_scaled_to_5pass: list[ObsRow] = []
    for row in rows:
        scale_d = row.entries_d / avg_entries_d_5pass
        scale_s = row.entries_s / avg_entries_s_5pass
        mu_d_err = row.mu_d_err * np.sqrt(scale_d)
        mu_s_err = row.mu_s_err * np.sqrt(scale_s)
        scaled_err = float(np.hypot(mu_d_err, mu_s_err))
        rows_scaled_to_5pass.append(with_measured_err(row, scaled_err))
    print(
        f"[INFO] 5-pass b1 reference entries: "
        f"avg_entriesD={avg_entries_d_5pass:.4f}, avg_entriesS={avg_entries_s_5pass:.4f}"
    )
    write_fit(
        OUTPUT_CSV_SCALED_TO_5PASS,
        rows_scaled_to_5pass,
        "A_B_errors_scaled_to_5pass_average_b1_entries",
    )


if __name__ == "__main__":
    main()
