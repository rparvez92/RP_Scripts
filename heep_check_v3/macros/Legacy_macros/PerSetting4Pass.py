#!/usr/bin/env python3
"""
PerSetting4Pass.py

4-pass per-setting SciPy/Nelder-Mead fit using the measured offsets snapshot in
results/tables/PerSetting4Pass/MeasuredOffsetsBySetting_4pass.csv.

For each setting independently:
  - dthe and dpe are fitted separately for b1..b4
  - dthp and dpp are common across b1..b4 within that setting

Run from heep_check_v3:
  python3 macros/PerSetting4Pass.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


OUT_DIR = Path("results/tables/PerSetting4Pass")
INPUT_CSV = OUT_DIR / "MeasuredOffsetsBySetting_4pass.csv"
PARAMS_CSV = OUT_DIR / "PerSetting4PassFitParams.csv"
RESIDUALS_CSV = OUT_DIR / "PerSetting4PassResiduals.csv"
SUMMARY_CSV = OUT_DIR / "PerSetting4PassSummary.csv"

SETTINGS = ("A", "B")
DP_LABELS = ("b1", "b2", "b3", "b4")
OBSERVABLES = ("W", "Em", "Pmz", "Pmy")
REQUIRED_INPUT_COLUMNS = {
    "setting",
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

# Rows: W, Em, Pmz, Pmy
# Columns: dthe, dpe, dthp, dpp
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

OPT_OPTIONS = {"maxiter": 30000, "maxfev": 50000, "xatol": 1e-10, "fatol": 1e-10}


@dataclass(frozen=True)
class ObsRow:
    setting: str
    dp_label: str
    delta_center: float
    dp_lo: float
    dp_hi: float
    variable: str
    measured: float
    measured_err: float
    sigma: float
    sigma_err: float


def f4(x: float) -> str:
    return f"{x:.4f}" if np.isfinite(x) else "nan"


def param_names(setting: str) -> tuple[str, ...]:
    return (
        *(f"dthe_{dp}_{setting}" for dp in DP_LABELS),
        *(f"dpe_{dp}_{setting}" for dp in DP_LABELS),
        f"dthp_{setting}",
        f"dpp_{setting}",
    )


def require_input_columns(reader: csv.DictReader) -> None:
    present = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - present)
    if missing:
        raise RuntimeError(f"Missing required columns in {INPUT_CSV}: {missing}")


def read_input() -> dict[str, list[ObsRow]]:
    rows_by_key: dict[tuple[str, str, str], ObsRow] = {}
    with INPUT_CSV.open() as f:
        reader = csv.DictReader(f)
        require_input_columns(reader)
        for row in reader:
            setting = row["setting"]
            dp_label = row["dp_label"]
            variable = row["var"]
            if setting not in SETTINGS or dp_label not in DP_LABELS or variable not in OBSERVABLES:
                continue
            if row["fit_valid"] != "1":
                continue
            rows_by_key[(setting, dp_label, variable)] = ObsRow(
                setting=setting,
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

    grouped: dict[str, list[ObsRow]] = {}
    missing: list[tuple[str, str, str]] = []
    for setting in SETTINGS:
        rows: list[ObsRow] = []
        for dp_label in DP_LABELS:
            for variable in OBSERVABLES:
                key = (setting, dp_label, variable)
                if key not in rows_by_key:
                    missing.append(key)
                else:
                    rows.append(rows_by_key[key])
        grouped[setting] = rows

    if missing:
        raise RuntimeError(f"Missing required rows in {INPUT_CSV}: {missing}")
    return grouped


def parameter_vector_for_row(params: np.ndarray, row: ObsRow) -> np.ndarray:
    ibin = DP_LABELS.index(row.dp_label)
    return np.array(
        [
            params[ibin],      # dthe_bin_setting
            params[4 + ibin],  # dpe_bin_setting
            params[8],         # dthp_setting
            params[9],         # dpp_setting
        ],
        dtype=float,
    )


def predict_one(params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(RESPONSE[row.setting][ivar] @ parameter_vector_for_row(params, row))


def chi2(params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    total = 0.0
    for i, row in enumerate(rows):
        pred = predict_one(params, row)
        total += ((measured[i] - pred) / row.measured_err) ** 2
    return float(total)


def fit_setting(rows: list[ObsRow]):
    measured = np.array([row.measured for row in rows], dtype=float)
    return minimize(
        lambda p: chi2(p, rows, measured),
        np.zeros(10, dtype=float),
        method="Nelder-Mead",
        options=OPT_OPTIONS,
    )


def main() -> None:
    grouped = read_input()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    for setting, rows in grouped.items():
        result = fit_setting(rows)
        measured = np.array([row.measured for row in rows], dtype=float)
        fit_chi2 = chi2(result.x, rows, measured)
        n_obs = len(rows)
        n_params = len(result.x)
        ndf = n_obs - n_params
        chi2_ndf = fit_chi2 / ndf if ndf > 0 else float("nan")
        results[setting] = {
            "rows": rows,
            "result": result,
            "chi2": fit_chi2,
            "n_obs": n_obs,
            "n_params": n_params,
            "ndf": ndf,
            "chi2_ndf": chi2_ndf,
        }

    all_param_names = tuple(name for setting in SETTINGS for name in param_names(setting))

    with PARAMS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "setting",
                *all_param_names,
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
        for setting in SETTINGS:
            record = results[setting]
            params = record["result"].x
            values_by_name = {name: "" for name in all_param_names}
            for name, value in zip(param_names(setting), params):
                values_by_name[name] = f4(value)
            writer.writerow(
                [
                    setting,
                    *[values_by_name[name] for name in all_param_names],
                    f4(record["chi2"]),
                    record["n_obs"],
                    record["n_params"],
                    record["ndf"],
                    f4(record["chi2_ndf"]),
                    bool(record["result"].success),
                    str(record["result"].message),
                    int(record["result"].nit),
                    int(record["result"].nfev),
                ]
            )

    with RESIDUALS_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "setting",
                "dp_label",
                "delta_center",
                "dp_lo",
                "dp_hi",
                "variable",
                "measured",
                "measured_err",
                "sigma",
                "sigma_err",
                "predicted",
                "residual",
                "pull",
                *all_param_names,
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
        for setting in SETTINGS:
            record = results[setting]
            params = record["result"].x
            values_by_name = {name: "" for name in all_param_names}
            for name, value in zip(param_names(setting), params):
                values_by_name[name] = f4(value)
            for row in record["rows"]:
                pred = predict_one(params, row)
                resid = row.measured - pred
                pull = resid / row.measured_err
                writer.writerow(
                    [
                        row.setting,
                        row.dp_label,
                        f4(row.delta_center),
                        f4(row.dp_lo),
                        f4(row.dp_hi),
                        row.variable,
                        f4(row.measured),
                        f4(row.measured_err),
                        f4(row.sigma),
                        f4(row.sigma_err),
                        f4(pred),
                        f4(resid),
                        f4(pull),
                        *[values_by_name[name] for name in all_param_names],
                        f4(record["chi2"]),
                        record["n_obs"],
                        record["n_params"],
                        record["ndf"],
                        f4(record["chi2_ndf"]),
                        bool(record["result"].success),
                        str(record["result"].message),
                        int(record["result"].nit),
                        int(record["result"].nfev),
                    ]
                )

    with SUMMARY_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["setting", "chi2", "n_obs", "n_params", "ndf", "chi2_ndf", "success", "nit", "nfev"])
        for setting in SETTINGS:
            record = results[setting]
            writer.writerow(
                [
                    setting,
                    f4(record["chi2"]),
                    record["n_obs"],
                    record["n_params"],
                    record["ndf"],
                    f4(record["chi2_ndf"]),
                    bool(record["result"].success),
                    int(record["result"].nit),
                    int(record["result"].nfev),
                ]
            )

    print(f"[INFO] Wrote {PARAMS_CSV}")
    print(f"[INFO] Wrote {RESIDUALS_CSV}")
    print(f"[INFO] Wrote {SUMMARY_CSV}")
    for setting in SETTINGS:
        record = results[setting]
        print(
            f"[INFO] {setting}: chi2/ndf={record['chi2_ndf']:.4f}, "
            f"success={bool(record['result'].success)}, nit={int(record['result'].nit)}, "
            f"nfev={int(record['result'].nfev)}"
        )


if __name__ == "__main__":
    main()
