#!/usr/bin/env python3
"""
Combined4and5PassFitParams.py

Central SciPy/Nelder-Mead fit using both 4-pass and 5-pass measured HEEP
offsets. This diagnostic model floats all setting-dependent offsets except
for one global HMS angle offset dthp.

Inputs:
  results/tables/MeasuredOffsetsBySetting_4pass.csv
  results/tables/MeasuredOffsetsBySetting_5pass.csv

Outputs:
  results/tables/Combined4and5Pass/Combined4and5PassFitParams.csv
  results/tables/Combined4and5Pass/Combined4and5PassResiduals.csv

Model:
  dthe_A, dthe_B, dthe_C, dthe_E, dthe_D
  dpe_A,  dpe_B,  dpe_C,  dpe_E,  dpe_D
  dpp_A,  dpp_B,  dpp_C,  dpp_E,  dpp_D
  dthp

The observable order is W, Em, Pmz, Pmx. The Pmx row uses the same response
matrix row previously mislabeled as Pmy.

Run from heep_check_v3:
  python3 macros/Combined4and5PassFitParams.py
"""

from __future__ import annotations

import csv
import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


OBSERVABLES = ("W", "Em", "Pmz", "Pmx")
SETTING_ORDER = ("A", "B", "C", "E", "D")
OPT_OPTIONS = {"maxiter": 50000, "maxfev": 80000, "xatol": 1e-10, "fatol": 1e-10}

INPUT_4PASS = Path("results/tables/MeasuredOffsetsBySetting_4pass.csv")
INPUT_5PASS = Path("results/tables/MeasuredOffsetsBySetting_5pass.csv")
OUT_DIR = Path("results/tables/Combined4and5Pass")
FITPARAMS_OUT = OUT_DIR / "Combined4and5PassFitParams.csv"
RESIDUALS_OUT = OUT_DIR / "Combined4and5PassResiduals.csv"

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

PARAM_NAMES = (
    "dthe_A",
    "dthe_B",
    "dthe_C",
    "dthe_E",
    "dthe_D",
    "dpe_A",
    "dpe_B",
    "dpe_C",
    "dpe_E",
    "dpe_D",
    "dpp_A",
    "dpp_B",
    "dpp_C",
    "dpp_E",
    "dpp_D",
    "dthp",
)

PARAM_INDEX = {name: i for i, name in enumerate(PARAM_NAMES)}


@dataclass(frozen=True)
class ObsRow:
    pass_tag: str
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


def f4(x: float) -> str:
    return f"{x:.4f}" if np.isfinite(x) else "nan"


def require_input_columns(reader: csv.DictReader, input_csv: Path) -> None:
    present = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - present)
    if missing:
        raise RuntimeError(f"Missing required columns in {input_csv}: {missing}")


def read_one_input(path: Path, pass_tag: str, settings: tuple[str, ...]) -> list[ObsRow]:
    by_key: dict[tuple[str, str, str], ObsRow] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        require_input_columns(reader, path)
        for row in reader:
            setting = row["setting"]
            dp_label = row["dp_label"]
            variable = row["var"]
            if setting not in settings or dp_label != "b1" or variable not in OBSERVABLES:
                continue
            if row["fit_valid"] != "1":
                continue
            measured_err = float(row["offset_err_MeV"])
            if not np.isfinite(measured_err) or measured_err <= 0.0:
                continue
            by_key[(setting, dp_label, variable)] = ObsRow(
                pass_tag=pass_tag,
                setting=setting,
                shms_p_GeV=float(row["shms_p_GeV"]),
                dp_label=dp_label,
                delta_center=float(row["bin_center"]),
                dp_lo=float(row["dp_lo"]),
                dp_hi=float(row["dp_hi"]),
                variable=variable,
                measured=float(row["offset_MeV"]),
                measured_err=measured_err,
                sigma=float(row["sigD_MeV"]),
                sigma_err=float(row["sigD_err_MeV"]),
            )

    rows: list[ObsRow] = []
    missing: list[tuple[str, str, str]] = []
    for setting in settings:
        for variable in OBSERVABLES:
            key = (setting, "b1", variable)
            if key not in by_key:
                missing.append(key)
            else:
                rows.append(by_key[key])
    if missing:
        raise RuntimeError(f"Missing required rows in {path}: {missing}")
    return rows


def read_inputs(input4: Path, input5: Path) -> list[ObsRow]:
    rows: list[ObsRow] = []
    rows.extend(read_one_input(input4, "4pass", ("A", "B")))
    rows.extend(read_one_input(input5, "5pass", ("C", "E", "D")))
    return rows


def response_matrix(row: ObsRow) -> np.ndarray:
    if row.pass_tag == "4pass":
        return RESPONSE_4PASS[row.setting]
    if row.pass_tag == "5pass":
        return RESPONSE_5PASS
    raise RuntimeError(f"Unexpected pass tag: {row.pass_tag}")


def param_vector(params: np.ndarray, row: ObsRow) -> np.ndarray:
    return np.array(
        [
            params[PARAM_INDEX[f"dthe_{row.setting}"]],
            params[PARAM_INDEX[f"dpe_{row.setting}"]],
            params[PARAM_INDEX["dthp"]],
            params[PARAM_INDEX[f"dpp_{row.setting}"]],
        ],
        dtype=float,
    )


def predict_one(params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(response_matrix(row)[ivar] @ param_vector(params, row))


def chi2(params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    total = 0.0
    for i, row in enumerate(rows):
        pred = predict_one(params, row)
        total += ((measured[i] - pred) / row.measured_err) ** 2
    return float(total)


def fit_combined(rows: list[ObsRow], measured: np.ndarray):
    return minimize(
        lambda p: chi2(p, rows, measured),
        np.zeros(len(PARAM_NAMES), dtype=float),
        method="Nelder-Mead",
        options=OPT_OPTIONS,
    )


def write_outputs(rows: list[ObsRow], result, out_dir: Path, fitparams_out: Path, residuals_out: Path) -> None:
    params = result.x
    measured = np.array([row.measured for row in rows], dtype=float)
    fit_chi2 = chi2(params, rows, measured)
    n_obs = len(rows)
    n_params = len(params)
    ndf = n_obs - n_params
    chi2_ndf = fit_chi2 / ndf if ndf > 0 else float("nan")

    out_dir.mkdir(parents=True, exist_ok=True)

    with fitparams_out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "fit_type",
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
        writer.writerow(
            [
                "combined_4and5_float_all_except_dthp",
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

    with residuals_out.open("w", newline="") as f:
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
                "sigma",
                "sigma_err",
                "predicted",
                "residual",
                "residual_err_FromMeasuredOnly",
                "pull_FromMeasuredOnly",
            ]
        )
        for row in rows:
            pred = predict_one(params, row)
            resid = row.measured - pred
            pull = resid / row.measured_err
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
                    f4(row.sigma),
                    f4(row.sigma_err),
                    f4(pred),
                    f4(resid),
                    f4(row.measured_err),
                    f4(pull),
                ]
            )

    print(f"[INFO] Wrote {fitparams_out}")
    print(f"[INFO] Wrote {residuals_out}")
    print(f"[INFO] combined 4+5 chi2/ndf = {chi2_ndf:.4f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input4", type=Path, default=INPUT_4PASS, help="4-pass measured-offset CSV")
    parser.add_argument("--input5", type=Path, default=INPUT_5PASS, help="5-pass measured-offset CSV")
    parser.add_argument("--outdir", type=Path, default=OUT_DIR, help="Output directory")
    args = parser.parse_args()

    fitparams_out = args.outdir / "Combined4and5PassFitParams.csv"
    residuals_out = args.outdir / "Combined4and5PassResiduals.csv"

    rows = read_inputs(args.input4, args.input5)
    measured = np.array([row.measured for row in rows], dtype=float)
    result = fit_combined(rows, measured)
    write_outputs(rows, result, args.outdir, fitparams_out, residuals_out)


if __name__ == "__main__":
    main()
