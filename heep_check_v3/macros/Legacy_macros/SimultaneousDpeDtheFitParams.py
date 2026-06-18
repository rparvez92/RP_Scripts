#!/usr/bin/env python3
"""
SimultaneousDpeDtheFitParams.py

Combined SciPy/Nelder-Mead fit for measured HEEP offsets.

Supported modes:
  4-pass:
    input  results/tables/MeasuredOffsetsBySetting_4pass.csv
    output results/tables/SimultaneousDpeDtheFitParams_4pass.csv
    model  dthe_b1_A/B, dpe_b1_A/B, common dthp, dpp_b1_A/B,
           using A/B-specific response matrices row by row.

  5-pass:
    input  results/tables/MeasuredOffsetsBySetting_5pass.csv
    output results/tables/SimultaneousDpeDtheFitParams_5pass.csv
    model  dthe_C/E/D, dpe_C/E/D, common dthp, common dpp,
           using the 5-pass response matrix.

Run from heep_check_v3:
  python3 macros/SimultaneousDpeDtheFitParams.py --pass 4pass
  python3 macros/SimultaneousDpeDtheFitParams.py --pass 5pass

Diagnostic input/output override example:
  python3 macros/SimultaneousDpeDtheFitParams.py --pass 4pass \
    --input results/tables/Scaled/MeasuredOffsetsBySetting_4passScaled.csv \
    --output results/tables/Scaled/SimultaneousDpeDtheFitParams_4passScaled.csv
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
    output_csv: Path
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
            output_csv=Path("results/tables/SimultaneousDpeDtheFitParams_4pass.csv"),
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
        output_csv=Path("results/tables/SimultaneousDpeDtheFitParams_5pass.csv"),
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


def predict_one(config: FitConfig, params: np.ndarray, row: ObsRow) -> float:
    ivar = OBSERVABLES.index(row.variable)
    return float(config.response_matrix(row)[ivar] @ config.param_vector(params, row))


def chi2(config: FitConfig, params: np.ndarray, rows: list[ObsRow], measured: np.ndarray) -> float:
    total = 0.0
    for i, row in enumerate(rows):
        pred = predict_one(config, params, row)
        total += ((measured[i] - pred) / row.measured_err) ** 2
    return float(total)


def fit_combined(config: FitConfig, rows: list[ObsRow], measured: np.ndarray):
    return minimize(
        lambda p: chi2(config, p, rows, measured),
        np.zeros(len(config.param_names), dtype=float),
        method="Nelder-Mead",
        options=OPT_OPTIONS,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pass", dest="pass_arg", required=True, help="4pass or 5pass")
    parser.add_argument("--input", type=Path, default=None, help="Optional measured-offset CSV override")
    parser.add_argument("--output", type=Path, default=None, help="Optional fit-parameter CSV override")
    args = parser.parse_args()

    config = make_config(args.pass_arg)
    if args.input is not None or args.output is not None:
        config = FitConfig(
            pass_tag=config.pass_tag,
            input_csv=args.input if args.input is not None else config.input_csv,
            output_csv=args.output if args.output is not None else config.output_csv,
            settings=config.settings,
            dp_labels=config.dp_labels,
            param_names=config.param_names,
            param_vector=config.param_vector,
            response_matrix=config.response_matrix,
        )
    rows = read_input(config)
    measured = np.array([row.measured for row in rows], dtype=float)
    result = fit_combined(config, rows, measured)
    params = result.x
    fit_chi2 = chi2(config, params, rows, measured)
    n_obs = len(rows)
    n_params = len(params)
    ndf = n_obs - n_params
    chi2_ndf = fit_chi2 / ndf if ndf > 0 else float("nan")

    config.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with config.output_csv.open("w", newline="") as f:
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
                *config.param_names,
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
            pred = predict_one(config, params, row)
            resid = row.measured - pred
            pull = resid / row.measured_err
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
                    f4(row.sigma),
                    f4(row.sigma_err),
                    f4(pred),
                    f4(resid),
                    f4(row.measured_err),
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

    print(f"[INFO] Wrote {config.output_csv}")
    print(f"[INFO] {config.pass_tag} chi2/ndf = {chi2_ndf:.4f}")


if __name__ == "__main__":
    main()
