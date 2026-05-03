#!/usr/bin/env python3
"""
SimultaneousDpeOnly.py

Python version of macros/SimultaneousDpeOnly.C.

This script reads setting-wise measured observable offsets from
results/tables/MeasuredOffsetsBySetting.csv and performs one simultaneous fit
per setting with:

  - one shared dthe
  - one shared dthp
  - one shared dpp
  - one bin-dependent dpe for each narrow bin

The objective function matches the C++ macro:

  chi2 = sum_i [ (measured_i - predicted_i)^2 / sigma_i^2 ]

where measured_i is offset_MeV and sigma_i is offset_err_MeV.

The fit is carried out with scipy.optimize.minimize(method="Nelder-Mead").
After convergence, a numerical Hessian is estimated from chi2 in order to
derive a covariance matrix and parameter uncertainties, following the same
strategy used in the C++ version.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

try:
    from scipy.optimize import minimize
    from scipy.stats import chi2 as scipy_chi2
except ImportError as exc:  # pragma: no cover - runtime dependency check
    print(
        "[ERROR] This script requires scipy. Install it first, for example with:\n"
        "        python3 -m pip install scipy",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc


@dataclass
class ObsSummary:
    value: float = math.nan
    err: float = math.nan
    ok: bool = False


@dataclass
class BinData:
    setting: str = ""
    dp_idx: int = -1
    bin: str = ""
    dp_lo: float = math.nan
    dp_hi: float = math.nan
    fit_valid: int = 0
    W: ObsSummary = field(default_factory=ObsSummary)
    Em: ObsSummary = field(default_factory=ObsSummary)
    Pmz: ObsSummary = field(default_factory=ObsSummary)
    Pmy: ObsSummary = field(default_factory=ObsSummary)


GLOBAL_PARAM_NAMES = ["dthe", "dpe_b1", "dpe_b2", "dpe_b3", "dpe_b4", "dpe_b5", "dthp", "dpp"]


def normalize_setting(text: str) -> str:
    t = text.strip().lower()
    if t in {"a", "setting a", "settinga"}:
        return "A"
    if t in {"b", "setting b", "settingb"}:
        return "B"
    return text.strip()


def is_allowed_bin_for_setting(setting: str, bin_label: str) -> bool:
    b = bin_label.strip().lower()
    if setting == "A":
        return b in {"b1", "b2", "b3", "b4"}
    if setting == "B":
        return b in {"b1", "b2", "b3", "b4", "b5"}
    return False


def bin_order_for_setting(setting: str) -> List[str]:
    if setting == "A":
        return ["b1", "b2", "b3", "b4"]
    if setting == "B":
        return ["b1", "b2", "b3", "b4", "b5"]
    return []


def to_float(text: str) -> float:
    try:
        return float(text.strip())
    except Exception:
        return math.nan


def to_int(text: str, default: int = 0) -> int:
    try:
        return int(text.strip())
    except Exception:
        return default


def read_measured_offsets(csv_path: str) -> Dict[str, Dict[str, BinData]]:
    required = {
        "setting",
        "dp_idx",
        "dp_label",
        "dp_lo",
        "dp_hi",
        "var",
        "offset_mev",
        "offset_err_mev",
        "fit_valid",
    }

    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Empty CSV: {csv_path}")
        lower_map = {name.lower(): name for name in reader.fieldnames}
        missing = sorted(required - set(lower_map))
        if missing:
            raise RuntimeError(
                f"Missing required columns in {csv_path}: {', '.join(missing)}"
            )

        out: Dict[str, Dict[str, BinData]] = {}
        for raw_row in reader:
            setting = normalize_setting(raw_row[lower_map["setting"]])
            bin_label = raw_row[lower_map["dp_label"]].strip()
            var = raw_row[lower_map["var"]].strip()
            fit_valid = to_int(raw_row[lower_map["fit_valid"]], 0)

            if bin_label.lower() == "full":
                continue
            if not is_allowed_bin_for_setting(setting, bin_label):
                continue
            if fit_valid != 1:
                continue

            row_bin = out.setdefault(setting, {}).setdefault(bin_label, BinData())
            row_bin.setting = setting
            row_bin.bin = bin_label
            row_bin.fit_valid = fit_valid
            row_bin.dp_idx = to_int(raw_row[lower_map["dp_idx"]], -1)
            row_bin.dp_lo = to_float(raw_row[lower_map["dp_lo"]])
            row_bin.dp_hi = to_float(raw_row[lower_map["dp_hi"]])

            obs = ObsSummary(
                value=to_float(raw_row[lower_map["offset_mev"]]),
                err=to_float(raw_row[lower_map["offset_err_mev"]]),
            )
            obs.ok = math.isfinite(obs.value) and math.isfinite(obs.err) and obs.err > 0.0

            if var == "W":
                row_bin.W = obs
            elif var == "Em":
                row_bin.Em = obs
            elif var == "Pmz":
                row_bin.Pmz = obs
            elif var == "Pmy":
                row_bin.Pmy = obs

    return out


def build_setting_equations(
    setting: str, by_bin: Dict[str, BinData]
) -> Tuple[List[str], List[dict], str]:
    bins: List[str] = []
    equations: List[dict] = []

    desired = bin_order_for_setting(setting)
    if not desired:
        return bins, equations, "unsupported_setting"

    for bin_label in desired:
        if bin_label not in by_bin:
            return bins, equations, f"missing_required_bin_{bin_label}"
        b = by_bin[bin_label]
        if not (b.W.ok and b.Em.ok and b.Pmz.ok and b.Pmy.ok):
            return bins, equations, f"missing_one_or_more_valid_observables_in_{bin_label}"
        bins.append(bin_label)

    for ibin, bin_label in enumerate(bins):
        b = by_bin[bin_label]
        if setting == "A":
            equations.extend(
                [
                    {"bin_index": ibin, "c_dthe": -14.08, "c_dpe": -8.62, "c_dthp": 0.0, "c_dpp": 0.0, "measured": b.W.value, "sigma": b.W.err},
                    {"bin_index": ibin, "c_dthe": 0.0, "c_dpe": -7.06, "c_dthp": 0.0, "c_dpp": -2.10, "measured": b.Em.value, "sigma": b.Em.err},
                    {"bin_index": ibin, "c_dthe": -5.75, "c_dpe": 4.10, "c_dthp": 0.0, "c_dpp": 2.27, "measured": b.Pmz.value, "sigma": b.Pmz.err},
                    {"bin_index": ibin, "c_dthe": 4.10, "c_dpe": 5.75, "c_dthp": -2.27, "c_dpp": 0.0, "measured": b.Pmy.value, "sigma": b.Pmy.err},
                ]
            )
        elif setting == "B":
            equations.extend(
                [
                    {"bin_index": ibin, "c_dthe": -17.33, "c_dpe": -8.62, "c_dthp": 0.0, "c_dpp": 0.0, "measured": b.W.value, "sigma": b.W.err},
                    {"bin_index": ibin, "c_dthe": 0.0, "c_dpe": -5.66, "c_dthp": 0.0, "c_dpp": -3.63, "measured": b.Em.value, "sigma": b.Em.err},
                    {"bin_index": ibin, "c_dthe": -4.30, "c_dpe": 3.69, "c_dthp": 0.0, "c_dpp": 3.74, "measured": b.Pmz.value, "sigma": b.Pmz.err},
                    {"bin_index": ibin, "c_dthe": 3.69, "c_dpe": 4.30, "c_dthp": -3.74, "c_dpp": 0.0, "measured": b.Pmy.value, "sigma": b.Pmy.err},
                ]
            )
        else:
            return [], [], "unsupported_setting"

    return bins, equations, ""


def build_param_names(bins: List[str]) -> List[str]:
    return ["dthe"] + [f"dpe_{bin_label}" for bin_label in bins] + ["dthp", "dpp"]


def evaluate_chi2(params: np.ndarray, equations: List[dict], dpe_index_by_bin: List[int]) -> float:
    idx_dthe = 0
    idx_dthp = len(dpe_index_by_bin) + 1
    idx_dpp = idx_dthp + 1

    chi2_value = 0.0
    for eq in equations:
        idx_dpe = dpe_index_by_bin[eq["bin_index"]]
        pred = (
            eq["c_dthe"] * params[idx_dthe]
            + eq["c_dpe"] * params[idx_dpe]
            + eq["c_dthp"] * params[idx_dthp]
            + eq["c_dpp"] * params[idx_dpp]
        )
        resid = eq["measured"] - pred
        chi2_value += resid * resid / (eq["sigma"] * eq["sigma"])
    return float(chi2_value)


def estimate_covariance_from_chi2(
    params: np.ndarray, equations: List[dict], dpe_index_by_bin: List[int]
) -> Tuple[np.ndarray, np.ndarray]:
    npar = len(params)
    f0 = evaluate_chi2(params, equations, dpe_index_by_bin)
    if not math.isfinite(f0):
        raise RuntimeError("non_finite_chi2_at_minimum")

    steps = np.maximum(1e-3, 1e-3 * np.maximum(1.0, np.abs(params)))
    hess = np.zeros((npar, npar), dtype=float)
    x = params.copy()

    for i in range(npar):
        hi = steps[i]
        x[i] = params[i] + hi
        fp = evaluate_chi2(x, equations, dpe_index_by_bin)
        x[i] = params[i] - hi
        fm = evaluate_chi2(x, equations, dpe_index_by_bin)
        x[i] = params[i]
        if not (math.isfinite(fp) and math.isfinite(fm)):
            raise RuntimeError("non_finite_hessian_diagonal")
        hess[i, i] = (fp - 2.0 * f0 + fm) / (hi * hi)

        for j in range(i + 1, npar):
            hj = steps[j]
            x[i], x[j] = params[i] + hi, params[j] + hj
            fpp = evaluate_chi2(x, equations, dpe_index_by_bin)
            x[i], x[j] = params[i] + hi, params[j] - hj
            fpm = evaluate_chi2(x, equations, dpe_index_by_bin)
            x[i], x[j] = params[i] - hi, params[j] + hj
            fmp = evaluate_chi2(x, equations, dpe_index_by_bin)
            x[i], x[j] = params[i] - hi, params[j] - hj
            fmm = evaluate_chi2(x, equations, dpe_index_by_bin)
            x[i], x[j] = params[i], params[j]

            if not all(math.isfinite(v) for v in (fpp, fpm, fmp, fmm)):
                raise RuntimeError("non_finite_hessian_off_diagonal")
            hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj)
            hess[i, j] = hij
            hess[j, i] = hij

    inv_hess = np.linalg.inv(hess)
    covariance = 2.0 * inv_hess
    variances = np.diag(covariance)
    if np.any(~np.isfinite(variances)) or np.any(variances < 0.0):
        raise RuntimeError("invalid_covariance_variances")

    return covariance, np.sqrt(variances)


def pred_variance_for_bin(
    covariance: np.ndarray,
    idx_dthe: int,
    idx_dpe: int,
    idx_dthp: int,
    idx_dpp: int,
    c_dthe: float,
    c_dpe: float,
    c_dthp: float,
    c_dpp: float,
) -> float:
    coeffs = np.array([c_dthe, c_dpe, c_dthp, c_dpp], dtype=float)
    idxs = np.array([idx_dthe, idx_dpe, idx_dthp, idx_dpp], dtype=int)
    subcov = covariance[np.ix_(idxs, idxs)]
    var = float(coeffs @ subcov @ coeffs)
    if not math.isfinite(var):
        return math.nan
    if var < 0.0 and abs(var) < 1e-12:
        return 0.0
    return var if var >= 0.0 else math.nan


def fit_prob_from_chi2(chi2_value: float, ndf: int) -> float:
    if ndf <= 0 or not math.isfinite(chi2_value):
        return math.nan
    return float(scipy_chi2.sf(chi2_value, ndf))


def param_value_or_nan(result: dict, name: str, want_error: bool) -> float:
    try:
        idx = result["param_names"].index(name)
    except ValueError:
        return math.nan
    source = result["param_errs"] if want_error else result["params"]
    return float(source[idx])


def covariance_or_nan(result: dict, pi: str, pj: str) -> float:
    try:
        i = result["param_names"].index(pi)
        j = result["param_names"].index(pj)
    except ValueError:
        return math.nan
    return float(result["covariance"][i, j])


def run_setting_fit(setting: str, by_bin: Dict[str, BinData]) -> dict:
    bins, equations, note = build_setting_equations(setting, by_bin)
    result = {
        "setting": setting,
        "fit_valid": 0,
        "fit_note": "",
        "param_names": [],
        "params": np.array([], dtype=float),
        "param_errs": np.array([], dtype=float),
        "covariance": np.empty((0, 0), dtype=float),
        "bins": [],
        "dthe": math.nan,
        "dthe_err": math.nan,
        "dthp": math.nan,
        "dthp_err": math.nan,
        "dpp": math.nan,
        "dpp_err": math.nan,
        "chi2_min": math.nan,
        "n_obs_used": 0,
        "n_params": 0,
        "ndf": -1,
        "chi2_ndf": math.nan,
        "fit_prob": math.nan,
    }
    if note:
        result["fit_note"] = note
        return result

    param_names = build_param_names(bins)
    npar = len(param_names)
    dpe_index_by_bin = [1 + i for i in range(len(bins))]
    x0 = np.zeros(npar, dtype=float)

    opt = minimize(
        lambda x: evaluate_chi2(np.asarray(x, dtype=float), equations, dpe_index_by_bin),
        x0,
        method="Nelder-Mead",
        options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000, "maxfev": 20000},
    )
    if not opt.success or not math.isfinite(float(opt.fun)):
        result["fit_note"] = "scipy_minimize_failed"
        return result

    params = np.asarray(opt.x, dtype=float)
    try:
        covariance, param_errs = estimate_covariance_from_chi2(params, equations, dpe_index_by_bin)
    except Exception:
        result["fit_note"] = "scipy_covariance_estimation_failed"
        return result

    idx_dthe = 0
    idx_dthp = npar - 2
    idx_dpp = npar - 1

    result.update(
        {
            "fit_valid": 1,
            "fit_note": "ok_scipy_nelder_mead",
            "param_names": param_names,
            "params": params,
            "param_errs": param_errs,
            "covariance": covariance,
            "dthe": float(params[idx_dthe]),
            "dthe_err": float(param_errs[idx_dthe]),
            "dthp": float(params[idx_dthp]),
            "dthp_err": float(param_errs[idx_dthp]),
            "dpp": float(params[idx_dpp]),
            "dpp_err": float(param_errs[idx_dpp]),
            "chi2_min": float(opt.fun),
            "n_obs_used": len(equations),
            "n_params": npar,
            "ndf": len(equations) - npar,
        }
    )
    if result["ndf"] > 0:
        result["chi2_ndf"] = result["chi2_min"] / result["ndf"]
        result["fit_prob"] = fit_prob_from_chi2(result["chi2_min"], result["ndf"])

    for ibin, bin_label in enumerate(bins):
        b = by_bin[bin_label]
        idx_dpe = 1 + ibin
        fit_bin = {
            "bin": bin_label,
            "dpe": float(params[idx_dpe]),
            "dpe_err": float(param_errs[idx_dpe]),
            "dW_measured": b.W.value,
            "dW_err": b.W.err,
            "dEm_measured": b.Em.value,
            "dEm_err": b.Em.err,
            "dPmz_measured": b.Pmz.value,
            "dPmz_err": b.Pmz.err,
            "dPmy_measured": b.Pmy.value,
            "dPmy_err": b.Pmy.err,
        }

        if setting == "A":
            fit_bin["dW_pred"] = -14.08 * result["dthe"] + -8.62 * fit_bin["dpe"]
            fit_bin["dEm_pred"] = -7.06 * fit_bin["dpe"] + -2.10 * result["dpp"]
            fit_bin["dPmz_pred"] = -5.75 * result["dthe"] + 4.10 * fit_bin["dpe"] + 2.27 * result["dpp"]
            fit_bin["dPmy_pred"] = 4.10 * result["dthe"] + 5.75 * fit_bin["dpe"] - 2.27 * result["dthp"]

            v_w = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -14.08, -8.62, 0.0, 0.0)
            v_em = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, 0.0, -7.06, 0.0, -2.10)
            v_pmz = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -5.75, 4.10, 0.0, 2.27)
            v_pmy = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, 4.10, 5.75, -2.27, 0.0)
        else:
            fit_bin["dW_pred"] = -17.33 * result["dthe"] + -8.62 * fit_bin["dpe"]
            fit_bin["dEm_pred"] = -5.66 * fit_bin["dpe"] + -3.63 * result["dpp"]
            fit_bin["dPmz_pred"] = -4.30 * result["dthe"] + 3.69 * fit_bin["dpe"] + 3.74 * result["dpp"]
            fit_bin["dPmy_pred"] = 3.69 * result["dthe"] + 4.30 * fit_bin["dpe"] - 3.74 * result["dthp"]

            v_w = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -17.33, -8.62, 0.0, 0.0)
            v_em = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, 0.0, -5.66, 0.0, -3.63)
            v_pmz = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -4.30, 3.69, 0.0, 3.74)
            v_pmy = pred_variance_for_bin(covariance, idx_dthe, idx_dpe, idx_dthp, idx_dpp, 3.69, 4.30, -3.74, 0.0)

        fit_bin["dW_pred_err"] = math.sqrt(v_w) if math.isfinite(v_w) else math.nan
        fit_bin["dEm_pred_err"] = math.sqrt(v_em) if math.isfinite(v_em) else math.nan
        fit_bin["dPmz_pred_err"] = math.sqrt(v_pmz) if math.isfinite(v_pmz) else math.nan
        fit_bin["dPmy_pred_err"] = math.sqrt(v_pmy) if math.isfinite(v_pmy) else math.nan

        fit_bin["dW_resid"] = fit_bin["dW_pred"] - fit_bin["dW_measured"]
        fit_bin["dEm_resid"] = fit_bin["dEm_pred"] - fit_bin["dEm_measured"]
        fit_bin["dPmz_resid"] = fit_bin["dPmz_pred"] - fit_bin["dPmz_measured"]
        fit_bin["dPmy_resid"] = fit_bin["dPmy_pred"] - fit_bin["dPmy_measured"]

        fit_bin["dW_resid_err"] = math.sqrt(fit_bin["dW_pred_err"] ** 2 + fit_bin["dW_err"] ** 2) if math.isfinite(fit_bin["dW_pred_err"]) and math.isfinite(fit_bin["dW_err"]) else math.nan
        fit_bin["dEm_resid_err"] = math.sqrt(fit_bin["dEm_pred_err"] ** 2 + fit_bin["dEm_err"] ** 2) if math.isfinite(fit_bin["dEm_pred_err"]) and math.isfinite(fit_bin["dEm_err"]) else math.nan
        fit_bin["dPmz_resid_err"] = math.sqrt(fit_bin["dPmz_pred_err"] ** 2 + fit_bin["dPmz_err"] ** 2) if math.isfinite(fit_bin["dPmz_pred_err"]) and math.isfinite(fit_bin["dPmz_err"]) else math.nan
        fit_bin["dPmy_resid_err"] = math.sqrt(fit_bin["dPmy_pred_err"] ** 2 + fit_bin["dPmy_err"] ** 2) if math.isfinite(fit_bin["dPmy_pred_err"]) and math.isfinite(fit_bin["dPmy_err"]) else math.nan

        result["bins"].append(fit_bin)

    return result


def format_float(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.6f}"


def write_output_csv(out_path: str, fit_results: List[dict], by_setting_bin: Dict[str, Dict[str, BinData]]) -> Tuple[int, int, int]:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    rows_written = 0
    valid_settings = 0
    invalid_settings = 0

    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        header = [
            "setting",
            "bin",
            "dp_idx",
            "dp_lo",
            "dp_hi",
            "dthe_fit",
            "dthe_fit_err",
        ]
        for name in GLOBAL_PARAM_NAMES:
            if name in {"dthe", "dthp", "dpp"}:
                continue
            header.extend([f"{name}_fit", f"{name}_fit_err"])
        header.extend(
            [
                "dthp_fit",
                "dthp_fit_err",
                "dpp_fit",
                "dpp_fit_err",
                "dW_measured",
                "dW_err",
                "dEm_measured",
                "dEm_err",
                "dPmz_measured",
                "dPmz_err",
                "dPmy_measured",
                "dPmy_err",
                "dW_pred",
                "dW_pred_err",
                "dEm_pred",
                "dEm_pred_err",
                "dPmz_pred",
                "dPmz_pred_err",
                "dPmy_pred",
                "dPmy_pred_err",
                "dW_resid",
                "dW_resid_err",
                "dEm_resid",
                "dEm_resid_err",
                "dPmz_resid",
                "dPmz_resid_err",
                "dPmy_resid",
                "dPmy_resid_err",
                "chi2_min",
                "n_obs_used",
                "n_params",
                "ndf",
                "chi2_ndf",
                "fit_prob",
                "fit_valid",
                "fit_note",
            ]
        )
        for pi in GLOBAL_PARAM_NAMES:
            for pj in GLOBAL_PARAM_NAMES:
                header.append(f"cov_{pi}_{pj}")
        header.append("end_marker")
        writer.writerow(header)

        for result in fit_results:
            if result["fit_valid"] != 1:
                invalid_settings += 1
                print(
                    f"[WARN] Setting {result['setting']} fit invalid: {result['fit_note']}",
                    file=sys.stderr,
                )
                continue

            valid_settings += 1
            by_bin = by_setting_bin[result["setting"]]
            for fit_bin in result["bins"]:
                b = by_bin[fit_bin["bin"]]
                row = [
                    result["setting"],
                    fit_bin["bin"],
                    str(b.dp_idx),
                    format_float(b.dp_lo),
                    format_float(b.dp_hi),
                    format_float(result["dthe"]),
                    format_float(result["dthe_err"]),
                ]
                for name in GLOBAL_PARAM_NAMES:
                    if name in {"dthe", "dthp", "dpp"}:
                        continue
                    row.extend(
                        [
                            format_float(param_value_or_nan(result, name, False)),
                            format_float(param_value_or_nan(result, name, True)),
                        ]
                    )
                row.extend(
                    [
                        format_float(result["dthp"]),
                        format_float(result["dthp_err"]),
                        format_float(result["dpp"]),
                        format_float(result["dpp_err"]),
                        format_float(fit_bin["dW_measured"]),
                        format_float(fit_bin["dW_err"]),
                        format_float(fit_bin["dEm_measured"]),
                        format_float(fit_bin["dEm_err"]),
                        format_float(fit_bin["dPmz_measured"]),
                        format_float(fit_bin["dPmz_err"]),
                        format_float(fit_bin["dPmy_measured"]),
                        format_float(fit_bin["dPmy_err"]),
                        format_float(fit_bin["dW_pred"]),
                        format_float(fit_bin["dW_pred_err"]),
                        format_float(fit_bin["dEm_pred"]),
                        format_float(fit_bin["dEm_pred_err"]),
                        format_float(fit_bin["dPmz_pred"]),
                        format_float(fit_bin["dPmz_pred_err"]),
                        format_float(fit_bin["dPmy_pred"]),
                        format_float(fit_bin["dPmy_pred_err"]),
                        format_float(fit_bin["dW_resid"]),
                        format_float(fit_bin["dW_resid_err"]),
                        format_float(fit_bin["dEm_resid"]),
                        format_float(fit_bin["dEm_resid_err"]),
                        format_float(fit_bin["dPmz_resid"]),
                        format_float(fit_bin["dPmz_resid_err"]),
                        format_float(fit_bin["dPmy_resid"]),
                        format_float(fit_bin["dPmy_resid_err"]),
                        format_float(result["chi2_min"]),
                        str(result["n_obs_used"]),
                        str(result["n_params"]),
                        str(result["ndf"]),
                        format_float(result["chi2_ndf"]),
                        format_float(result["fit_prob"]),
                        str(result["fit_valid"]),
                        result["fit_note"],
                    ]
                )
                for pi in GLOBAL_PARAM_NAMES:
                    for pj in GLOBAL_PARAM_NAMES:
                        row.append(format_float(covariance_or_nan(result, pi, pj)))
                row.append("done")
                writer.writerow(row)
                rows_written += 1

    return rows_written, valid_settings, invalid_settings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Python simultaneous Dpe-only fit.")
    parser.add_argument(
        "--in",
        dest="in_csv",
        default="results/tables/MeasuredOffsetsBySetting.csv",
        help="Input MeasuredOffsetsBySetting CSV.",
    )
    parser.add_argument(
        "--out",
        dest="out_csv",
        default="results/tables/SimultaneousDpeOnly_python.csv",
        help="Output CSV path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        by_setting_bin = read_measured_offsets(args.in_csv)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    if not by_setting_bin:
        print(f"[ERROR] No valid rows read from {args.in_csv}", file=sys.stderr)
        return 1

    fit_results: List[dict] = []
    for setting in ["A", "B"]:
        if setting not in by_setting_bin:
            continue
        fit_results.append(run_setting_fit(setting, by_setting_bin[setting]))

    if not fit_results:
        print(f"[ERROR] No supported settings found in {args.in_csv}", file=sys.stderr)
        return 1

    rows_written, valid_settings, invalid_settings = write_output_csv(
        args.out_csv, fit_results, by_setting_bin
    )
    print(f"[INFO] Wrote: {args.out_csv}")
    print(f"[INFO] Rows written: {rows_written}")
    print(f"[INFO] Valid setting fits: {valid_settings}")
    print(f"[INFO] Invalid setting fits: {invalid_settings}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

