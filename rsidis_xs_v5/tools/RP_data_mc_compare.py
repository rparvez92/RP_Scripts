#!/usr/bin/env python3
"""Build setting-wise, random/positron/dummy-corrected data-to-SIMC comparisons."""

from __future__ import annotations

import argparse
import csv
import ctypes
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_NORM_DIR = PROJECT_DIR / "results" / "Tables" / "RP_get_coin_normyield"
DEFAULT_BIGTABLE_DIR = PROJECT_DIR / "bigtable"
DEFAULT_SIMC_CATALOG = (
    PROJECT_DIR / "results" / "Tables" / "RP_extract_simc_info"
    / "RP_extract_simc_info.csv"
)
DEFAULT_TABLE_DIR = PROJECT_DIR / "results" / "Tables" / "RP_data_mc_compare"
DEFAULT_PDF_DIR = PROJECT_DIR / "results" / "PDFs" / "RP_data_mc_compare"
DEFAULT_T7_ROOT = Path("/Volumes/T7/RSIDIS")

REACTIONS = ("sidis", "rho", "delta", "exclusive")
TIERS = ("delta", "full")
CANVAS_WIDTH = 1800
CANVAS_HEIGHT = 1200
METADATA_LINE_CHARS = 90
METADATA_RUN_LINES_FIRST_PAGE = 8
METADATA_RUN_LINES_CONTINUATION = 18
S_DUMMY = {
    ("phase1", "LH2"): 1.0 / 3.5274,
    ("phase1", "LD2"): 1.0 / 3.7825,
    ("phase2", "LH2"): 1.0 / 3.9031,
    ("phase2", "LD2"): 1.0 / 4.1854,
}

PID_CUT = (
    "(P_gtr_p<=2.7 || P_hgcer_npeSum>1) && P_aero_npeSum>2 && "
    "P_cal_etottracknorm<0.8 && H_cer_npeSum>2 && H_cal_etottracknorm>0.8"
)
DATA_DELTA = "H_gtr_dp>-8 && H_gtr_dp<8 && P_gtr_dp>-10 && P_gtr_dp<22"
DATA_ANGLE_CUTS = {
    "hxp": ("H_gtr_th", -0.15, 0.15),
    "hyp": ("H_gtr_ph", -0.10, 0.10),
    "pxp": ("P_gtr_th", -0.15, 0.15),
    "pyp": ("P_gtr_ph", -0.10, 0.10),
}
SIMC_DELTA = "hsdelta>-8 && hsdelta<8 && ssdelta>-10 && ssdelta<22"
SIMC_ANGLE_CUTS = {
    "hxp": ("hsxptar", -0.15, 0.15),
    "hyp": ("hsyptar", -0.10, 0.10),
    "pxp": ("ssxptar", -0.15, 0.15),
    "pyp": ("ssyptar", -0.10, 0.10),
}

# Versioned fixed comparison bins. They intentionally extend beyond acceptance cuts.
BINNING_VERSION = "RP_DMC_V3"
BINNING = {
    "hdelta": (20, -12.0, 12.0, "HMS #delta", "%"),
    "pdelta": (20, -15.0, 25.0, "SHMS #delta", "%"),
    "hxp": (20, -0.20, 0.20, "HMS x'_{tar}", "rad"),
    "hyp": (20, -0.15, 0.15, "HMS y'_{tar}", "rad"),
    "pxp": (20, -0.20, 0.20, "SHMS x'_{tar}", "rad"),
    "pyp": (20, -0.15, 0.15, "SHMS y'_{tar}", "rad"),
    "hytar": (20, -5.0, 5.0, "HMS y_{tar}", "cm"),
    "pytar": (20, -5.0, 5.0, "SHMS y_{tar}", "cm"),
    "xbj": (20, 0.0, 0.6, "x_{bj}", ""),
    "Q2": (20, 0.0, 6.0, "Q^{2}", "GeV^{2}"),
    "W": (20, 0.0, 4.0, "W", "GeV"),
    "z": (20, 0.0, 1.0, "z", ""),
    "thetapq": (20, 0.0, 0.20, "#theta_{pq}", "rad"),
    "phipq": (20, 0.0, 2.0 * math.pi, "#phi_{pq}", "rad"),
    "pt2": (20, 0.0, 0.6, "p_{T}^{2}", "GeV^{2}"),
    "missmass": (20, 0.0, 3.0, "Missing mass", "GeV"),
}
ACCEPTANCE_VARIABLES = tuple(list(BINNING)[:8])
PHYSICS_VARIABLES = tuple(list(BINNING)[8:])

DATA_EXPRESSIONS = {
    "hdelta": "H_gtr_dp", "pdelta": "P_gtr_dp",
    "hxp": "H_gtr_th", "hyp": "H_gtr_ph",
    "pxp": "P_gtr_th", "pyp": "P_gtr_ph",
    "hytar": "H_gtr_y", "pytar": "P_gtr_y",
    "xbj": "H_kin_primary_x_bj", "Q2": "H_kin_primary_Q2",
    "W": "H_kin_primary_W",
    "z": "P_gtr_p/H_kin_primary_nu",
    "thetapq": "P_kin_secondary_th_xq",
    "phipq": "P_kin_secondary_ph_xq<0 ? P_kin_secondary_ph_xq+2*TMath::Pi() : P_kin_secondary_ph_xq",
    "pt2": "pt*pt", "missmass": "mmass",
}
SIMC_EXPRESSIONS = {
    "hdelta": "hsdelta", "pdelta": "ssdelta",
    "hxp": "hsxptar", "hyp": "hsyptar",
    "pxp": "ssxptar", "pyp": "ssyptar",
    "hytar": "hsytar", "pytar": "ssytar",
    "xbj": "xbj_recon", "Q2": "Q2_recon", "W": "W_recon",
    "z": "z_recon", "thetapq": "thetapq_recon",
    "phipq": "phipq_recon<0 ? phipq_recon+2*TMath::Pi() : phipq_recon",
    "pt2": "pt2_recon", "missmass": "missmass_recon",
}

CSV_COLUMNS = [
    "Phase", "Pass", "Run_type", "Target", "Setting", "Tier", "Variable",
    "Bin_index", "Bin_low", "Bin_high", "Binning_version", "S_dummy",
    "Target_Elec", "Target_Elec_err", "Target_Pos", "Target_Pos_err",
    "Target_pos_sub", "Target_pos_sub_err", "Dummy_Elec", "Dummy_Elec_err",
    "Dummy_Pos", "Dummy_Pos_err", "Dummy_pos_sub", "Dummy_pos_sub_err",
    "Dummy_scaled", "Dummy_scaled_err", "Data_final", "Data_final_err",
    "MC_sidis", "MC_sidis_err", "MC_rho", "MC_rho_err",
    "MC_delta", "MC_delta_err", "MC_exclusive", "MC_exclusive_err",
    "MC_total", "MC_total_err", "Data_by_MC", "Data_by_MC_err", "Pull",
    "Data_underflow", "Data_overflow", "MC_underflow", "MC_overflow",
    "Data_flow_fraction", "MC_flow_fraction", "Dummy_fraction",
    "Data_closure_rel", "MC_complete", "Provisional_flag", "Status", "Reason",
]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def atomic_csv(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def parse_norm_stem(path: Path) -> Dict[str, str]:
    parts = path.stem.split("_", 5)
    if len(parts) != 6:
        raise ValueError(f"unexpected normalized-yield filename: {path.name}")
    return dict(zip(("Phase", "Pass", "Run_type", "Target", "Charge", "Setting"), parts))


def setting_key(identity: Mapping[str, str]) -> Tuple[str, str, str, str]:
    return tuple(identity[key] for key in ("Phase", "Pass", "Run_type", "Setting"))


def discover_norm_files(directory: Path) -> Dict[Tuple[str, str, str, str], Dict[Tuple[str, str], Path]]:
    groups: Dict[Tuple[str, str, str, str], Dict[Tuple[str, str], Path]] = defaultdict(dict)
    for path in sorted(directory.glob("phase*_pass*_PI*_*.csv")):
        if path.name.endswith(("_Summary.csv", "_Problematic.csv")):
            continue
        identity = parse_norm_stem(path)
        groups[setting_key(identity)][(identity["Target"], identity["Charge"])] = path
    return groups


def leaf_bigtable_path(norm_path: Path, bigtable_dir: Path) -> Path:
    identity = parse_norm_stem(norm_path)
    return (
        bigtable_dir / identity["Phase"] / identity["Pass"]
        / identity["Run_type"] / identity["Target"] / identity["Charge"]
        / identity["Setting"] / f"{norm_path.stem}.csv"
    )


def rows_by_run(path: Path, run_field: str) -> Dict[int, Dict[str, str]]:
    output = {}
    for row in read_csv(path):
        if not finite(row.get(run_field)):
            continue
        output[int(float(row[run_field]))] = row
    return output


def signal_stability_points(
    norm_path: Path, bigtable_dir: Path
) -> Tuple[Dict[str, List[Tuple[float, float, float, float]]], List[str]]:
    """Return tier -> (run, current, yield, run) tuples and diagnostics."""
    leaf_path = leaf_bigtable_path(norm_path, bigtable_dir)
    if not leaf_path.is_file():
        return {tier: [] for tier in TIERS}, [f"BIGTABLE_LEAF_MISSING={leaf_path}"]
    leaf_rows = rows_by_run(leaf_path, "run")
    points = {tier: [] for tier in TIERS}
    reasons = []
    for row in read_csv(norm_path):
        if row.get("Norm_status") != "OK" or not finite(row.get("Run")):
            continue
        run = int(float(row["Run"]))
        leaf = leaf_rows.get(run)
        current = leaf.get("BCM2_I", "") if leaf else ""
        if not finite(current) or float(current) == -999:
            reasons.append(f"BCM2_I_UNAVAILABLE_RUN_{run}")
            continue
        for tier in TIERS:
            if not (
                finite(row.get(f"RP_Normyield_{tier}"))
                and finite(row.get(f"RP_Normyield_err_{tier}"))
            ):
                continue
            points[tier].append((
                float(row["Run"]), float(current),
                float(row[f"RP_Normyield_{tier}"]),
                float(row[f"RP_Normyield_err_{tier}"]),
            ))
    return points, sorted(set(reasons))


def filtered_simc_warning(reason: str) -> str:
    retained = [
        item for item in reason.split(";")
        if item and "_NEFF=" not in item.upper()
    ]
    return ";".join(retained)


def data_root_path(t7: Path, phase: str, run: int) -> Path:
    roots = (
        t7 / ("Phase1" if phase == "phase1" else "Phase2")
        / ("Pass0p1" if phase == "phase1" else "Passm1")
        / "ROOTfiles" / "SkimmedData"
    )
    candidates = (
        roots / f"skimmed_coin_replay_production_{run}_-1.root",
        roots / f"coin_replay_production_{run}_-1.root",
    )
    return next((path for path in candidates if path.is_file()), candidates[0])


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def make_hist(name: str, variable: str) -> ROOT.TH1D:
    bins, low, high, title, unit = BINNING[variable]
    axis = title + (f" [{unit}]" if unit else "")
    hist = ROOT.TH1D(name, f";{axis};Normalized yield", bins, low, high)
    hist.Sumw2()
    hist.SetDirectory(0)
    return hist


def add_scaled(destination: ROOT.TH1, source: ROOT.TH1, scale: float) -> None:
    destination.Add(source, scale)


def full_cut(cuts: Mapping[str, Tuple[str, float, float]], delta: str) -> str:
    return delta + " && " + " && ".join(
        f"{expression}>{low:.12g} && {expression}<{high:.12g}"
        for expression, low, high in cuts.values()
    )


def nminus1_cut(variable: str, simc: bool, tier: str = "full") -> str:
    delta = SIMC_DELTA if simc else DATA_DELTA
    angle_cuts = SIMC_ANGLE_CUTS if simc else DATA_ANGLE_CUTS
    delta_parts = {
        "hdelta": ("hsdelta>-8 && hsdelta<8" if simc else "H_gtr_dp>-8 && H_gtr_dp<8"),
        "pdelta": ("ssdelta>-10 && ssdelta<22" if simc else "P_gtr_dp>-10 && P_gtr_dp<22"),
    }
    parts = []
    for name, expression in delta_parts.items():
        if name != variable:
            parts.append(expression)
    if tier == "full":
        for name, (expression, low, high) in angle_cuts.items():
            if name != variable:
                parts.append(f"{expression}>{low:.12g} && {expression}<{high:.12g}")
    return " && ".join(parts) or "1"


def dataframe_histograms(
    path: Path,
    tree_name: str,
    variables: Sequence[str],
    expressions: Mapping[str, str],
    selection: str,
    weight_expression: str,
    prefix: str,
) -> Dict[str, ROOT.TH1D]:
    frame = ROOT.RDataFrame(tree_name, path.as_posix()).Filter(selection).Define(
        "_rp_dmc_weight", weight_expression
    )
    actions = {}
    for index, variable in enumerate(variables):
        column = f"_rp_dmc_v_{index}"
        bins, low, high, _, _ = BINNING[variable]
        selected = frame.Define(column, expressions[variable]).Filter(
            f"std::isfinite(static_cast<double>({column})) && {column}>-998"
        )
        actions[variable] = selected.Histo1D(
            (f"{prefix}_{variable}", "", bins, low, high), column, "_rp_dmc_weight"
        )
    ROOT.RDF.RunGraphs(list(actions.values()))
    result = {}
    for variable, action in actions.items():
        hist = action.GetValue().Clone(f"{prefix}_{variable}_owned")
        hist.SetDirectory(0)
        hist.Sumw2()
        result[variable] = hist
    return result


def class_histograms(
    csv_path: Path,
    t7_root: Path,
    tier: str,
    variables: Sequence[str],
    prefix: str,
    selection_override: str | None = None,
    check_closure: bool = True,
) -> Tuple[Dict[str, ROOT.TH1D], List[str], float]:
    rows = read_csv(csv_path)
    valid = [
        row for row in rows
        if row.get("Norm_status") == "OK"
        and finite(row.get("BCM2_Q")) and float(row["BCM2_Q"]) > 0
        and finite(row.get(f"RP_Normyield_{tier}"))
    ]
    output = {variable: make_hist(f"{prefix}_{tier}_{variable}", variable) for variable in variables}
    excluded = [
        row.get("Run", "?") for row in rows if row.get("Norm_status") != "OK"
    ]
    initial_reasons = (
        ["EXCLUDED_INVALID_NORM_RUNS=" + "|".join(excluded)] if excluded else []
    )
    if not valid:
        return output, initial_reasons + ["NO_VALID_NORMALIZED_RUNS"], math.nan
    total_charge = sum(float(row["BCM2_Q"]) for row in valid)
    reasons: List[str] = initial_reasons
    expected = 0.0
    for row in valid:
        run = int(float(row["Run"]))
        path = data_root_path(t7_root, parse_norm_stem(csv_path)["Phase"], run)
        if not path.is_file():
            reasons.append(f"DATA_ROOT_MISSING_RUN_{run}")
            continue
        rf = float(row["RFperiodNs"])
        mean = float(row["CTmean"])
        # RP_get_good_coin_ev integrates a 400-bin, 0--100 ns histogram.
        # Reproduce its FindFixBin(low), FindFixBin(high)-1 semantics exactly.
        align = lambda value: math.floor(value / 0.25) * 0.25
        central_low, central_high = (
            align(mean - rf / 2.0), align(mean + rf / 2.0)
        )
        random_ranges = [
            (
                align(mean + offset * rf - rf / 2.0),
                align(mean + offset * rf + rf / 2.0),
            )
            for offset in (-4, -3, -2, 2, 3, 4)
        ]
        random_test = " || ".join(
            f"(CTime_ePiCoinTime_ROC2>={low:.12g} && CTime_ePiCoinTime_ROC2<{high:.12g})"
            for low, high in random_ranges
        )
        timing_weight = (
            f"(CTime_ePiCoinTime_ROC2>={central_low:.12g} && "
            f"CTime_ePiCoinTime_ROC2<{central_high:.12g}) ? 1.0 : "
            f"(({random_test}) ? -1.0/6.0 : 0.0)"
        )
        tier_cut = (
            selection_override if selection_override is not None
            else DATA_DELTA if tier == "delta"
            else full_cut(DATA_ANGLE_CUTS, DATA_DELTA)
        )
        selection = f"({PID_CUT}) && ({tier_cut}) && ({timing_weight})!=0"
        run_hists = dataframe_histograms(
            path, "T", variables, DATA_EXPRESSIONS, selection,
            f"({float(row['Norm_factor']):.17g})*({timing_weight})",
            f"{prefix}_{tier}_{run}",
        )
        alpha = float(row["BCM2_Q"]) / total_charge
        for variable in variables:
            add_scaled(output[variable], run_hists[variable], alpha)
        expected += alpha * float(row[f"RP_Normyield_{tier}"])
    closure = output[variables[0]].Integral(
        0, output[variables[0]].GetNbinsX() + 1
    )
    relative = abs(closure - expected) / max(abs(expected), 1e-30)
    # Published normalization inputs are formatted to four decimal places.
    if check_closure and relative > 5e-3:
        reasons.append(f"DATA_CLOSURE_FAILED={relative:.6g}")
    return output, reasons, relative


def simc_expressions(reaction: str) -> Dict[str, str]:
    expressions = dict(SIMC_EXPRESSIONS)
    if reaction in {"delta", "exclusive"}:
        expressions.update({
            "xbj": (
                "(std::isfinite(static_cast<double>(Q2_recon)) && Q2_recon>-998 && "
                "std::isfinite(static_cast<double>(nu_recon)) && nu_recon>-998 && "
                "abs(nu_recon)>1e-12) ? Q2_recon/(2.0*0.9382720813*nu_recon) : "
                "std::numeric_limits<double>::quiet_NaN()"
            ),
            "z": (
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(nu_recon)) && nu_recon>-998 && "
                "abs(nu_recon)>1e-12) ? phad/nu_recon : "
                "std::numeric_limits<double>::quiet_NaN()"
            ),
            "pt2": (
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(thetapq_recon)) && thetapq_recon>-998) ? "
                "pow(phad*sin(thetapq_recon),2) : "
                "std::numeric_limits<double>::quiet_NaN()"
            ),
        })
    return expressions


def simc_histograms(
    catalog_row: Mapping[str, str],
    tier: str,
    variables: Sequence[str],
    prefix: str,
    selection_override: str | None = None,
) -> Tuple[Dict[str, ROOT.TH1D], List[str]]:
    output = {variable: make_hist(f"{prefix}_{tier}_{variable}", variable) for variable in variables}
    path = Path(catalog_row.get("Recon_root_path", ""))
    if catalog_row.get("QA_status") not in {"OK", "WARNING"}:
        return output, [f"SIMC_QA_{catalog_row.get('QA_status','MISSING')}"]
    if not path.is_file():
        return output, ["SIMC_ROOT_MISSING"]
    selection = (
        selection_override if selection_override is not None
        else SIMC_DELTA if tier == "delta"
        else full_cut(SIMC_ANGLE_CUTS, SIMC_DELTA)
    )
    try:
        histograms = dataframe_histograms(
            path, "h10", variables, simc_expressions(catalog_row["Reaction"]),
            selection, "static_cast<double>(fWeight)", prefix,
        )
        expected = number_or_nan(catalog_row.get(f"SimYield_{tier}", ""))
        observed = histograms[variables[0]].Integral(
            0, histograms[variables[0]].GetNbinsX() + 1
        )
        reasons = []
        relative = abs(observed - expected) / max(abs(expected), 1e-30)
        if selection_override is None and (
            not math.isfinite(expected) or relative > 1e-5
        ):
            reasons.append(f"SIMC_CLOSURE_FAILED={relative:.6g}")
        return histograms, reasons
    except Exception as error:
        return output, [f"SIMC_HIST_FAILED={type(error).__name__}:{error}"]


def number_or_nan(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def combine_data(classes: Mapping[str, Dict[str, ROOT.TH1D]], variable: str,
                 s_dummy: float) -> Dict[str, ROOT.TH1D]:
    def clone(key: str, name: str) -> ROOT.TH1D:
        hist = classes[key][variable].Clone(name)
        hist.SetDirectory(0)
        return hist
    target_e = clone("target_e", f"{variable}_target_e")
    target_p = clone("target_p", f"{variable}_target_p")
    dummy_e = clone("dummy_e", f"{variable}_dummy_e")
    dummy_p = clone("dummy_p", f"{variable}_dummy_p")
    target_sub = target_e.Clone(f"{variable}_target_sub"); target_sub.Add(target_p, -1)
    dummy_sub = dummy_e.Clone(f"{variable}_dummy_sub"); dummy_sub.Add(dummy_p, -1)
    dummy_scaled = dummy_sub.Clone(f"{variable}_dummy_scaled"); dummy_scaled.Scale(s_dummy)
    final = target_sub.Clone(f"{variable}_final"); final.Add(dummy_scaled, -1)
    return {
        "target_e": target_e, "target_p": target_p, "target_sub": target_sub,
        "dummy_e": dummy_e, "dummy_p": dummy_p, "dummy_sub": dummy_sub,
        "dummy_scaled": dummy_scaled, "final": final,
    }


def add_mc_total(components: Mapping[str, Dict[str, ROOT.TH1D]], variable: str) -> ROOT.TH1D:
    total = make_hist(f"{variable}_mc_total", variable)
    for reaction in REACTIONS:
        total.Add(components[reaction][variable])
    return total


def bin_value(hist: ROOT.TH1, index: int) -> Tuple[float, float]:
    return hist.GetBinContent(index), hist.GetBinError(index)


def rows_for_variable(
    identity: Mapping[str, str], tier: str, variable: str,
    data: Mapping[str, ROOT.TH1D], mc: Mapping[str, Dict[str, ROOT.TH1D]],
    mc_total: ROOT.TH1D, s_dummy: float, status: str, reasons: Sequence[str],
    closure: float, mc_complete: bool, available_reactions: set[str],
) -> List[Dict[str, object]]:
    final = data["final"]
    data_under, data_over = final.GetBinContent(0), final.GetBinContent(final.GetNbinsX() + 1)
    mc_under, mc_over = mc_total.GetBinContent(0), mc_total.GetBinContent(mc_total.GetNbinsX() + 1)
    data_regular_abs = sum(abs(final.GetBinContent(i)) for i in range(1, final.GetNbinsX() + 1))
    mc_regular_abs = sum(abs(mc_total.GetBinContent(i)) for i in range(1, mc_total.GetNbinsX() + 1))
    data_flow_fraction = (
        (abs(data_under) + abs(data_over))
        / max(data_regular_abs + abs(data_under) + abs(data_over), 1e-30)
    )
    mc_flow_fraction = (
        (abs(mc_under) + abs(mc_over))
        / max(mc_regular_abs + abs(mc_under) + abs(mc_over), 1e-30)
    )
    local_reasons = list(reasons)
    if data_flow_fraction > 0.01:
        local_reasons.append(f"DATA_FLOW_GT_1PCT_{variable}={data_flow_fraction:.6g}")
    if mc_complete and mc_flow_fraction > 0.01:
        local_reasons.append(f"MC_FLOW_GT_1PCT_{variable}={mc_flow_fraction:.6g}")
    local_status = status
    if local_reasons and status == "OK":
        local_status = "WARNING"
    output = []
    for index in range(1, final.GetNbinsX() + 1):
        values = {
            key: bin_value(hist, index) for key, hist in data.items()
        }
        mc_values = {
            reaction: bin_value(mc[reaction][variable], index) for reaction in REACTIONS
        }
        d, de = values["final"]
        m, me = bin_value(mc_total, index)
        ratio = d / m if mc_complete and finite(m) and m != 0 else math.nan
        ratio_error = (
            math.sqrt((de / m) ** 2 + (d * me / (m * m)) ** 2)
            if finite(ratio) else math.nan
        )
        pull_denominator = math.hypot(de, me)
        pull = (
            (d - m) / pull_denominator
            if mc_complete and m != 0 and pull_denominator > 0 else math.nan
        )
        dummy_fraction = (
            values["dummy_scaled"][0] / values["target_sub"][0]
            if values["target_sub"][0] != 0 else math.nan
        )
        row: Dict[str, object] = {
            **identity, "Tier": tier, "Variable": variable, "Bin_index": index,
            "Bin_low": final.GetXaxis().GetBinLowEdge(index),
            "Bin_high": final.GetXaxis().GetBinUpEdge(index),
            "Binning_version": BINNING_VERSION, "S_dummy": s_dummy,
            "Target_Elec": values["target_e"][0], "Target_Elec_err": values["target_e"][1],
            "Target_Pos": values["target_p"][0], "Target_Pos_err": values["target_p"][1],
            "Target_pos_sub": values["target_sub"][0], "Target_pos_sub_err": values["target_sub"][1],
            "Dummy_Elec": values["dummy_e"][0], "Dummy_Elec_err": values["dummy_e"][1],
            "Dummy_Pos": values["dummy_p"][0], "Dummy_Pos_err": values["dummy_p"][1],
            "Dummy_pos_sub": values["dummy_sub"][0], "Dummy_pos_sub_err": values["dummy_sub"][1],
            "Dummy_scaled": values["dummy_scaled"][0], "Dummy_scaled_err": values["dummy_scaled"][1],
            "Data_final": d, "Data_final_err": de,
            "MC_total": m if mc_complete else "", "MC_total_err": me if mc_complete else "",
            "Data_by_MC": ratio, "Data_by_MC_err": ratio_error, "Pull": pull,
            "Data_underflow": data_under, "Data_overflow": data_over,
            "MC_underflow": mc_under, "MC_overflow": mc_over,
            "Data_flow_fraction": data_flow_fraction,
            "MC_flow_fraction": mc_flow_fraction if mc_complete else "",
            "Dummy_fraction": dummy_fraction,
            "Data_closure_rel": closure, "MC_complete": int(mc_complete),
            "Provisional_flag": 1,
            "Status": local_status,
            "Reason": ";".join(sorted(set(local_reasons))),
        }
        for reaction in REACTIONS:
            if reaction in available_reactions:
                row[f"MC_{reaction}"], row[f"MC_{reaction}_err"] = mc_values[reaction]
            else:
                row[f"MC_{reaction}"], row[f"MC_{reaction}_err"] = "", ""
        output.append(row)
    return output


def metadata_row(identity: Mapping[str, str], status: str, reason: str,
                 s_dummy: float | str = "") -> Dict[str, object]:
    return {
        **identity, "Tier": "", "Variable": "", "Binning_version": BINNING_VERSION,
        "S_dummy": s_dummy, "MC_complete": 0, "Provisional_flag": 1,
        "Status": status, "Reason": reason,
    }


def catalog_index(path: Path) -> Dict[Tuple[str, str, str, str, str, str], Dict[str, str]]:
    output = {}
    for row in read_csv(path):
        key = tuple(row[field] for field in (
            "Phase", "Pass", "Run_type", "Target", "Setting", "Reaction"
        ))
        output[key] = row
    return output


def draw_text_page(
    canvas: ROOT.TCanvas, pdf: Path, title: str, lines: Sequence[str],
    text_size: float = 0.026,
) -> None:
    canvas.Clear()
    canvas.cd()
    canvas.SetMargin(0.025, 0.025, 0.025, 0.025)
    box = ROOT.TPaveText(0.035, 0.035, 0.965, 0.965, "NDC")
    box.SetFillColor(0)
    box.SetFillStyle(0)
    box.SetBorderSize(0)
    box.SetTextAlign(12)
    box.SetTextFont(42)
    box.SetTextSize(text_size)
    box.AddText(title)
    for line in lines:
        box.AddText(line)
    box.Draw()
    canvas._rp_text_page = box
    canvas.Print(pdf.as_posix())


def run_list_lines(
    files: Mapping[Tuple[str, str], Path], target: str
) -> List[str]:
    categories = (
        ("Signal e-", (target, "Elec")),
        ("Signal e+", (target, "Pos")),
        ("Dummy e-", ("DUMMY", "Elec")),
        ("Dummy e+", ("DUMMY", "Pos")),
    )
    output = []
    for label, category in categories:
        path = files.get(category)
        runs = []
        if path and path.is_file():
            runs = sorted({
                int(float(row["Run"])) for row in read_csv(path)
                if finite(row.get("Run"))
            })
        prefix = f"Runs - {label}: "
        if not runs:
            output.append(prefix + "none")
            continue
        current = prefix
        continuation_prefix = f"Runs - {label} (cont.): "
        for run in map(str, runs):
            separator = "" if current in (prefix, continuation_prefix) else ", "
            if len(current) + len(separator) + len(run) > METADATA_LINE_CHARS:
                output.append(current)
                current = continuation_prefix + run
            else:
                current += separator + run
        output.append(current)
    return output


def draw_metadata_pages(
    canvas: ROOT.TCanvas, pdf: Path, identity: Mapping[str, str],
    status: str, s_dummy: float,
    files: Mapping[Tuple[str, str], Path],
) -> None:
    target = identity["Target"]
    run_lines = run_list_lines(files, target)
    page1_prefix = [
        ", ".join(f"{key}={value}" for key, value in identity.items()),
        f"Status: {status}",
        f"S_dummy = {s_dummy:.4f} (cell-wall thickness / dummy-slab thickness)",
    ]
    page1_suffix = [
        "Data: Per-run random subtraction, charge-weighted run combination, then positron and scaled-dummy subtraction.",
        "SIMC: Event weight = fWeight, Normalized Yield = sum(fWeight), Total = sidis + rho + delta + exclusive.",
        "PROVISIONAL: comp_livetime=1, trigger efficiency=1, PID efficiency=1.",
    ]
    page1 = [
        *page1_prefix,
        *run_lines[:METADATA_RUN_LINES_FIRST_PAGE],
        *page1_suffix,
    ]
    draw_text_page(canvas, pdf, "RP Data-to-MC Comparison - Identity", page1)
    remaining_run_lines = run_lines[METADATA_RUN_LINES_FIRST_PAGE:]
    continuation = 1
    while remaining_run_lines:
        page_lines = remaining_run_lines[:METADATA_RUN_LINES_CONTINUATION]
        remaining_run_lines = remaining_run_lines[METADATA_RUN_LINES_CONTINUATION:]
        draw_text_page(
            canvas, pdf,
            f"RP Data-to-MC Comparison - Run inventory (continued {continuation})",
            page_lines,
        )
        continuation += 1

    page2 = [
        "Data random subtraction: central coincidence window minus the average of six random windows at RF offsets -4, -3, -2, +2, +3, +4.",
        "Run combination: each normalized run distribution is weighted by its fraction of the total setting charge.",
        "Corrected data = (Signal e- - Signal e+) - S_dummy x (Dummy e- - Dummy e+).",
        "Data Cuts: Delta-only = PID + Delta Acceptance",
        "Data Cuts: Full-cut = PID + Delta Acceptance + Angle Acceptance",
        "MC Cuts: Delta-only = Delta Acceptance",
        "MC Cuts: Full-cut = Delta Acceptance + Angle Acceptance",
        "PID: (P_gtr_p<=2.7 || P_hgcer_npeSum>1) && P_aero_npeSum>2 && P_cal_etottracknorm<0.8",
        "     && H_cer_npeSum>2 && H_cal_etottracknorm>0.8",
        "Delta Acceptance: -8 < HMS delta < 8 && -10 < SHMS delta < 22",
        "Angle Acceptance: -0.15 < HMS xptar < 0.15 && -0.10 < HMS yptar < 0.10",
        "                  && -0.15 < SHMS xptar < 0.15 && -0.10 < SHMS yptar < 0.10",
        "Cut-removed (N-1) Plots: Every cut applied except the plotted-variable cut.",
    ]
    draw_text_page(canvas, pdf, "RP Data-to-MC Comparison - Calculation and cuts", page2)

    bin_lines = []
    for variable, (bins, low, high, title, unit) in BINNING.items():
        suffix = f" {unit}" if unit else ""
        high_text = "2*pi" if variable == "phipq" else f"{high:g}"
        bin_lines.append(
            f"{title}: ({bins}, {low:g}, {high_text}){suffix}"
        )
    draw_text_page(
        canvas, pdf, "RP Data-to-MC Comparison - Fixed binning",
        bin_lines, text_size=0.025,
    )


def draw_run_stability(
    canvas: ROOT.TCanvas,
    pdf: Path,
    points: Mapping[str, Sequence[Tuple[float, float, float, float]]],
) -> None:
    canvas.Clear()
    canvas.cd()
    pads = (
        ROOT.TPad("stability_run", "", 0.0, 0.5, 1.0, 1.0),
        ROOT.TPad("stability_current", "", 0.0, 0.0, 1.0, 0.5),
    )
    for pad in pads:
        pad.Draw()
    series = (
        ("delta", ROOT.kBlue + 1, "Delta-only"),
        ("full", ROOT.kRed + 1, "Full-cut"),
    )
    keepalive = []
    for pad_index, (x_index, title, x_title) in enumerate((
        (0, "Signal Normalized Yield Vs Run", "Run"),
        (1, "Signal Normalized Yield Vs Current", "BCM2_I (uA)"),
    ), start=1):
        pad = pads[pad_index - 1]
        pad.cd(); pad.SetGrid(); pad.SetLeftMargin(0.10); pad.SetRightMargin(0.035)
        multigraph = ROOT.TMultiGraph()
        legend = ROOT.TLegend(0.77, 0.76, 0.95, 0.92)
        legend.SetBorderSize(0); legend.SetFillStyle(0)
        for tier, color, label in series:
            graph = ROOT.TGraphErrors()
            for run, current, value, error in points[tier]:
                point = graph.GetN()
                graph.SetPoint(point, (run, current)[x_index], value)
                graph.SetPointError(point, 0, error)
            graph.SetMarkerStyle(20); graph.SetMarkerColor(color)
            graph.SetLineColor(color); graph.SetMarkerSize(1.25)
            multigraph.Add(graph, "P"); legend.AddEntry(graph, label, "lep")
            keepalive.append(graph)
        multigraph.SetTitle(f"{title};{x_title};Normalized Yield")
        multigraph.Draw("A P"); legend.Draw()
        keepalive.extend((multigraph, legend))
    canvas._rp_run_stability = [*pads, *keepalive]
    canvas.Print(pdf.as_posix())


def draw_diagnostic_text(
    canvas: ROOT.TCanvas, pdf: Path, title: str, lines: Sequence[str]
) -> None:
    draw_text_page(canvas, pdf, title, lines)


def transparent_legend(x1: float = 0.70, y1: float = 0.64,
                       x2: float = 0.94, y2: float = 0.91) -> ROOT.TLegend:
    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    return legend


def comparison_y_range(
    payloads: Mapping[str, Tuple[Mapping[str, ROOT.TH1D],
                                 Mapping[str, Dict[str, ROOT.TH1D]], ROOT.TH1D]],
    variable: str,
) -> Tuple[float, float]:
    upper, lower = [0.0], [0.0]
    for data, mc, total in payloads.values():
        histograms = [data["final"], total] + [mc[r][variable] for r in REACTIONS]
        for hist in histograms:
            for index in range(1, hist.GetNbinsX() + 1):
                upper.append(hist.GetBinContent(index) + hist.GetBinError(index))
                lower.append(hist.GetBinContent(index) - hist.GetBinError(index))
    ymin, ymax = min(lower), max(upper)
    span = max(ymax - ymin, 1e-12)
    return ymin - 0.08 * span, ymax + 0.15 * span


def draw_not_applicable(pad: ROOT.TPad, title: str) -> None:
    pad.cd()
    pad.SetMargin(0.08, 0.04, 0.08, 0.08)
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextAlign(22)
    text.SetTextSize(0.045)
    text.DrawLatex(0.5, 0.55, title)
    text.SetTextSize(0.075)
    text.DrawLatex(0.5, 0.45, "Not Applicable")
    pad._rp_na_text = text


def draw_comparison_pair(
    canvas: ROOT.TCanvas, pdf: Path, variable: str,
    payloads: Mapping[str, Tuple[Mapping[str, ROOT.TH1D],
                                 Mapping[str, Dict[str, ROOT.TH1D]], ROOT.TH1D]],
    cut_removed: bool = False, delta_applicable: bool = True,
) -> None:
    canvas.Clear()
    canvas.cd()
    pads = {}
    for column, tier in enumerate(TIERS):
        x1, x2 = 0.5 * column, 0.5 * (column + 1)
        pads[(tier, "top")] = ROOT.TPad(
            f"pad_{tier}_{variable}_top", "", x1, 0.30, x2, 1.0
        )
        pads[(tier, "bottom")] = ROOT.TPad(
            f"pad_{tier}_{variable}_bottom", "", x1, 0.0, x2, 0.30
        )
        for pad in (pads[(tier, "top")], pads[(tier, "bottom")]):
            pad.Draw()
    if cut_removed and not delta_applicable:
        draw_not_applicable(
            pads[("delta", "top")],
            f"Delta-only, Cut-removed (N-1): {BINNING[variable][3]}",
        )

    ymin, ymax = comparison_y_range(payloads, variable)
    colors = {
        "sidis": ROOT.kGreen + 2, "rho": ROOT.kMagenta + 1,
        "delta": ROOT.kOrange + 7, "exclusive": ROOT.kCyan + 2,
    }
    keepalive = list(pads.values())
    shared_legend = None
    for tier in TIERS:
        if tier == "delta" and cut_removed and not delta_applicable:
            continue
        data, mc, total = payloads[tier]
        final = data["final"]
        axis_title = final.GetXaxis().GetTitle()
        tier_label = "Delta-only" if tier == "delta" else "Full-cut"
        heading = (
            f"{tier_label}, Cut-removed (N-1): {BINNING[variable][3]}"
            if cut_removed else f"{tier_label}: {BINNING[variable][3]}"
        )
        top = pads[(tier, "top")]
        top.cd(); top.SetGrid()
        top.SetLeftMargin(0.15); top.SetRightMargin(0.04)
        top.SetBottomMargin(0.02); top.SetTopMargin(0.10)
        final.SetTitle(f"{heading};{axis_title};Normalized Yield")
        final.SetMinimum(ymin); final.SetMaximum(ymax)
        final.GetXaxis().SetLabelSize(0)
        final.GetXaxis().SetTitleSize(0)
        final.GetYaxis().SetTitleSize(0.052)
        final.GetYaxis().SetLabelSize(0.043)
        final.GetYaxis().SetTitleOffset(1.35)
        final.SetMarkerStyle(20)
        final.SetMarkerColor(ROOT.kBlue + 1)
        final.SetLineColor(ROOT.kBlue + 1)
        final.Draw("E1")
        legend = transparent_legend()
        legend.AddEntry(final, "Corrected data", "lep")
        for reaction in REACTIONS:
            hist = mc[reaction][variable]
            hist.SetLineColor(colors[reaction]); hist.SetLineWidth(2)
            hist.SetFillColorAlpha(colors[reaction], 0.10)
            hist.Draw("E2 SAME"); hist.Draw("HIST SAME")
            legend.AddEntry(hist, reaction, "l")
        total.SetLineColor(ROOT.kBlack); total.SetLineWidth(3)
        total.SetFillColorAlpha(ROOT.kGray + 1, 0.25)
        total.Draw("E2 SAME"); total.Draw("HIST SAME")
        legend.AddEntry(total, "MC total", "l")
        final.Draw("E1 SAME")
        if tier == "delta" and delta_applicable:
            legend.Draw()
            shared_legend = legend
        keepalive.append(legend)

        if cut_removed:
            boundaries = {
                "hdelta": (-8.0, 8.0), "pdelta": (-10.0, 22.0),
                **{name: (values[1], values[2])
                   for name, values in DATA_ANGLE_CUTS.items()},
            }[variable]
            top.Update()
            for boundary in boundaries:
                line = ROOT.TLine(
                    boundary, ROOT.gPad.GetUymin(),
                    boundary, ROOT.gPad.GetUymax(),
                )
                line.SetLineColor(ROOT.kRed + 1)
                line.SetLineStyle(2); line.SetLineWidth(2)
                line.Draw()
                keepalive.append(line)

        bottom = pads[(tier, "bottom")]
        bottom.cd(); bottom.SetGrid()
        bottom.SetLeftMargin(0.15); bottom.SetRightMargin(0.04)
        bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.30)
        ratio_frame = final.Clone(f"ratio_frame_{tier}_{variable}")
        ratio_frame.Reset("ICES")
        ratio_frame.SetTitle(";"+final.GetXaxis().GetTitle()+";Data / MC")
        ratio_frame.SetMinimum(0.0); ratio_frame.SetMaximum(2.0)
        ratio_frame.GetXaxis().SetTitleSize(0.12)
        ratio_frame.GetXaxis().SetLabelSize(0.095)
        ratio_frame.GetYaxis().SetTitleSize(0.105)
        ratio_frame.GetYaxis().SetLabelSize(0.085)
        ratio_frame.GetYaxis().SetTitleOffset(0.62)
        ratio_frame.Draw("")
        ratio = ROOT.TGraphErrors()
        ratio.SetName(f"ratio_{tier}_{variable}")
        for index in range(1, final.GetNbinsX() + 1):
            data_value, data_error = bin_value(final, index)
            mc_value, mc_error = bin_value(total, index)
            if not (finite(mc_value) and mc_value != 0):
                continue
            value = data_value / mc_value
            error = math.sqrt(
                (data_error / mc_value) ** 2
                + (data_value * mc_error / (mc_value * mc_value)) ** 2
            )
            point = ratio.GetN()
            ratio.SetPoint(point, final.GetXaxis().GetBinCenter(index), value)
            ratio.SetPointError(point, 0.0, error)
        ratio.SetMarkerStyle(20); ratio.SetMarkerSize(0.9)
        ratio.SetMarkerColor(ROOT.kBlack); ratio.SetLineColor(ROOT.kBlack)
        ratio.Draw("P SAME")
        line = ROOT.TLine(
            final.GetXaxis().GetXmin(), 1,
            final.GetXaxis().GetXmax(), 1,
        )
        line.SetLineStyle(2); line.Draw()
        keepalive.extend((ratio_frame, ratio, line))

    if cut_removed and not delta_applicable:
        pads[("delta", "top")].cd()
        data, mc, total = payloads["full"]
        legend = transparent_legend()
        legend.AddEntry(data["final"], "Corrected data", "lep")
        for reaction in REACTIONS:
            legend.AddEntry(mc[reaction][variable], reaction, "l")
        legend.AddEntry(total, "MC total", "l")
        legend.Draw()
        keepalive.append(legend)
    canvas._rp_comparison_objects = keepalive
    canvas.Print(pdf.as_posix())


def integrated_value(hist: ROOT.TH1D) -> Tuple[float, float]:
    error = ctypes.c_double(0.0)
    value = hist.IntegralAndError(0, hist.GetNbinsX() + 1, error)
    return value, error.value


def draw_waterfall(
    canvas: ROOT.TCanvas, pdf: Path,
    comparisons: Mapping[str, Mapping[str, ROOT.TH1D]],
) -> None:
    canvas.Clear(); canvas.cd()
    pad = ROOT.TPad("waterfall_pad", "", 0.0, 0.0, 1.0, 1.0)
    pad.Draw(); pad.cd(); pad.SetGrid()
    pad.SetLeftMargin(0.10); pad.SetRightMargin(0.03)
    pad.SetTopMargin(0.08); pad.SetBottomMargin(0.16)
    labels = (
        "Signal e-", "Signal After e+", "Dummy After e+",
        "Scaled Dummy", "Final",
    )
    keys = ("target_e", "target_sub", "dummy_sub", "dummy_scaled", "final")
    frame = ROOT.TH1D("waterfall_frame", "Normalized Yield Correction;;Normalized Yield",
                      5, 0.5, 5.5)
    frame.SetDirectory(0)
    frame.SetLineColor(0)
    values = []
    for index, label in enumerate(labels, start=1):
        frame.GetXaxis().SetBinLabel(index, label)
    graphs = []
    for tier, color, offset in (
        ("delta", ROOT.kBlue + 1, -0.08),
        ("full", ROOT.kRed + 1, 0.08),
    ):
        graph = ROOT.TGraphErrors()
        for index, key in enumerate(keys, start=1):
            value, error = integrated_value(comparisons[tier][key])
            graph.SetPoint(index - 1, index + offset, value)
            graph.SetPointError(index - 1, 0.0, error)
            values.extend((value - error, value + error))
        graph.SetMarkerStyle(20); graph.SetMarkerSize(1.5)
        graph.SetMarkerColor(color); graph.SetLineColor(color)
        graphs.append((tier, graph))
    ymin, ymax = min([0.0, *values]), max([0.0, *values])
    span = max(ymax - ymin, 1e-12)
    frame.SetMinimum(ymin - 0.08 * span); frame.SetMaximum(ymax + 0.15 * span)
    frame.LabelsOption("h", "X"); frame.Draw()
    legend = transparent_legend(0.78, 0.78, 0.95, 0.92)
    for tier, graph in graphs:
        graph.Draw("P SAME")
        legend.AddEntry(
            graph, "Delta-only" if tier == "delta" else "Full-cut", "lep"
        )
    legend.Draw()
    canvas._rp_waterfall_objects = (pad, frame, graphs, legend)
    canvas.Print(pdf.as_posix())


def process_setting(
    key: Tuple[str, str, str, str, str],
    files: Mapping[Tuple[str, str], Path],
    catalog: Mapping[Tuple[str, str, str, str, str, str], Dict[str, str]],
    t7_root: Path,
    bigtable_dir: Path,
    table_dir: Path,
    pdf_dir: Path,
    no_pdf: bool,
) -> str:
    phase, pass_name, run_type, target, setting = key
    identity = {
        "Phase": phase, "Pass": pass_name, "Run_type": run_type,
        "Target": target, "Setting": setting,
    }
    stem = "_".join(identity.values())
    output_csv = table_dir / f"{stem}.csv"
    output_pdf = pdf_dir / f"{stem}.pdf"
    s_dummy = S_DUMMY[(phase, target)]
    setting_values = []
    for source in files.values():
        source_rows = read_csv(source)
        if source_rows:
            setting_values.extend(float(source_rows[0].get(field, "-999"))
                                  for field in ("x", "Q2", "z", "thpq"))
    if any(value == -999 for value in setting_values):
        atomic_csv(output_csv, [
            metadata_row(identity, "SKIPPED", "SENTINEL_SETTING_IDENTITY", s_dummy)
        ])
        return "SKIPPED"
    required_classes = {
        "target_e": (target, "Elec"), "target_p": (target, "Pos"),
        "dummy_e": ("DUMMY", "Elec"), "dummy_p": ("DUMMY", "Pos"),
    }
    missing_classes = [name for name, category in required_classes.items() if category not in files]
    if missing_classes:
        reason = "MISSING_DATA_CLASSES=" + "|".join(missing_classes)
        atomic_csv(output_csv, [metadata_row(identity, "PENDING", reason, s_dummy)])
        if not no_pdf:
            pdf_dir.mkdir(parents=True, exist_ok=True)
            canvas = ROOT.TCanvas(
                f"c_rp_dmc_pending_{abs(hash(stem))}", "RP DMC",
                CANVAS_WIDTH, CANVAS_HEIGHT,
            )
            canvas.Print(output_pdf.as_posix() + "[")
            draw_metadata_pages(
                canvas, output_pdf, identity, "PENDING", s_dummy, files
            )
            draw_diagnostic_text(
                canvas, output_pdf, "Pending data-to-MC comparison", [reason]
            )
            canvas.Print(output_pdf.as_posix() + "]")
        return "PENDING"
    simc_rows = {}
    missing_mc = []
    structural_mc = []
    catalog_warnings = []
    for reaction in REACTIONS:
        row = catalog.get((phase, pass_name, run_type, target, setting, reaction))
        if row is None or row.get("QA_status") in {"PENDING", "SKIPPED"}:
            missing_mc.append(reaction)
        elif row.get("QA_status") == "ERROR":
            structural_mc.append(reaction)
        else:
            simc_rows[reaction] = row
            if row.get("QA_status") == "WARNING":
                warning = filtered_simc_warning(row.get("QA_reason", ""))
                if warning:
                    catalog_warnings.append(
                        f"SIMC_WARNING_{reaction}={warning}"
                    )

    all_rows: List[Dict[str, object]] = []
    reasons: List[str] = []
    status = "ERROR" if structural_mc else "PENDING" if missing_mc else "OK"
    if missing_mc:
        reasons.append("INCOMPLETE_MC=" + "|".join(missing_mc))
    if structural_mc:
        reasons.append("STRUCTURAL_SIMC_ERROR=" + "|".join(structural_mc))
    reasons.extend(catalog_warnings)
    stability_points, stability_reasons = signal_stability_points(
        files[(target, "Elec")], bigtable_dir
    )
    reasons.extend(stability_reasons)
    canvas = ROOT.TCanvas(
        f"c_rp_dmc_{abs(hash(stem))}", "RP DMC",
        CANVAS_WIDTH, CANVAS_HEIGHT,
    )
    pdf_dir.mkdir(parents=True, exist_ok=True)
    if not no_pdf:
        canvas.Print(output_pdf.as_posix() + "[")
        draw_metadata_pages(canvas, output_pdf, identity, status, s_dummy, files)
        draw_run_stability(canvas, output_pdf, stability_points)
    variables = tuple(BINNING)
    closure_report: Dict[str, Dict[str, float]] = {}
    tier_comparisons = {}
    tier_components = {}
    tier_totals = {}
    tier_complete = {}
    for tier in TIERS:
        classes = {}
        closures = []
        closure_report[tier] = {}
        for name, category in required_classes.items():
            histograms, class_reasons, closure = class_histograms(
                files[category], t7_root, tier, variables, f"{stem}_{name}"
            )
            classes[name] = histograms
            reasons.extend(class_reasons)
            if finite(closure):
                closures.append(closure)
                closure_report[tier][name] = closure
        components = {}
        tier_mc_complete = not missing_mc and not structural_mc
        tier_available_reactions: set[str] = set()
        for reaction in REACTIONS:
            if reaction in simc_rows:
                components[reaction], component_reasons = simc_histograms(
                    simc_rows[reaction], tier, variables, f"{stem}_{reaction}_{tier}"
                )
                reasons.extend(component_reasons)
                if component_reasons:
                    tier_mc_complete = False
                else:
                    tier_available_reactions.add(reaction)
            else:
                components[reaction] = {
                    variable: make_hist(f"{stem}_{reaction}_{tier}_{variable}", variable)
                    for variable in variables
                }
        comparison = {}
        totals = {}
        for variable in variables:
            comparison[variable] = combine_data(classes, variable, s_dummy)
            totals[variable] = add_mc_total(components, variable)
            variable_status = "WARNING" if reasons and status == "OK" else status
            all_rows.extend(rows_for_variable(
                identity, tier, variable, comparison[variable], components,
                totals[variable], s_dummy, variable_status, reasons,
                max(closures) if closures else math.nan,
                tier_mc_complete, tier_available_reactions,
            ))
        tier_comparisons[tier] = comparison
        tier_components[tier] = components
        tier_totals[tier] = totals
        tier_complete[tier] = tier_mc_complete
    if not no_pdf:
        draw_waterfall(
            canvas, output_pdf,
            {tier: tier_comparisons[tier]["hdelta"] for tier in TIERS},
        )
        if all(tier_complete.values()):
            for variable in variables:
                draw_comparison_pair(
                    canvas, output_pdf, variable,
                    {
                        tier: (
                            tier_comparisons[tier][variable],
                            tier_components[tier],
                            tier_totals[tier][variable],
                        )
                        for tier in TIERS
                    },
                )
    # Tier-aware acceptance diagnostics with only the plotted cut removed.
    if not missing_mc and not structural_mc:
        for variable in ("hdelta", "pdelta", "hxp", "hyp", "pxp", "pyp"):
            nminus_payloads = {}
            applicable_tiers = TIERS if variable in {"hdelta", "pdelta"} else ("full",)
            for tier in applicable_tiers:
                nminus_classes = {}
                for name, category in required_classes.items():
                    nminus_classes[name], class_reasons, _ = class_histograms(
                        files[category], t7_root, tier, (variable,),
                        f"{stem}_{name}_{tier}_nminus1_{variable}",
                        selection_override=nminus1_cut(variable, False, tier),
                        check_closure=False,
                    )
                    reasons.extend(class_reasons)
                nminus_components = {}
                for reaction in REACTIONS:
                    nminus_components[reaction], component_reasons = simc_histograms(
                        simc_rows[reaction], tier, (variable,),
                        f"{stem}_{reaction}_{tier}_nminus1_{variable}",
                        selection_override=nminus1_cut(variable, True, tier),
                    )
                    reasons.extend(component_reasons)
                nminus_data = combine_data(nminus_classes, variable, s_dummy)
                nminus_total = add_mc_total(nminus_components, variable)
                nminus_payloads[tier] = (
                    nminus_data, nminus_components, nminus_total
                )
            if "delta" not in nminus_payloads:
                nminus_payloads["delta"] = nminus_payloads["full"]
            if not no_pdf:
                draw_comparison_pair(
                    canvas, output_pdf, variable, nminus_payloads,
                    cut_removed=True,
                    delta_applicable=variable in {"hdelta", "pdelta"},
                )
    final_status = (
        "WARNING"
        if status == "OK" and (
            reasons or any(row.get("Status") == "WARNING" for row in all_rows)
        )
        else status
    )
    final_reason = ";".join(sorted(set(reasons)))
    for row in all_rows:
        if row["Status"] == "OK" and final_status != "OK":
            row["Status"] = final_status
        existing = str(row.get("Reason", ""))
        row["Reason"] = ";".join(
            sorted({item for item in (existing + ";" + final_reason).split(";") if item})
        )
    if not no_pdf:
        diagnostic_lines = []
        for tier in TIERS:
            label = "Delta-only" if tier == "delta" else "Full-cut"
            for name, value in closure_report.get(tier, {}).items():
                diagnostic_lines.append(
                    f"{label} {name} integrated closure relative difference = {value:.6g}"
                )
        diagnostic_lines.extend(sorted({
            item for item in reasons if "CLOSURE" in item
        }))
        diagnostic_lines.append("Final status: " + final_status)
        draw_diagnostic_text(
            canvas, output_pdf, "Integrated closure and final status",
            [line for line in diagnostic_lines if line],
        )
        canvas.Print(output_pdf.as_posix() + "]")
    atomic_csv(output_csv, all_rows)
    return final_status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="one normalized-yield setting CSV")
    parser.add_argument("--all", action="store_true", help="process every discovered setting")
    parser.add_argument("--norm-dir", type=Path, default=DEFAULT_NORM_DIR)
    parser.add_argument("--bigtable-dir", type=Path, default=DEFAULT_BIGTABLE_DIR)
    parser.add_argument("--simc-catalog", type=Path, default=DEFAULT_SIMC_CATALOG)
    parser.add_argument("--t7-root", type=Path, default=DEFAULT_T7_ROOT)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--pdf-dir", type=Path, default=DEFAULT_PDF_DIR)
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    if bool(args.input) == bool(args.all):
        parser.error("provide exactly one input CSV or --all")
    groups = discover_norm_files(args.norm_dir)
    if args.input:
        identity = parse_norm_stem(args.input)
        keys = [(
            identity["Phase"], identity["Pass"], identity["Run_type"],
            identity["Target"], identity["Setting"],
        )]
    else:
        keys = sorted(
            (*base[:3], target, base[3])
            for base, files in groups.items()
            for target in ("LH2", "LD2")
            if (target, "Elec") in files
        )
    catalog = catalog_index(args.simc_catalog)
    counts = defaultdict(int)
    for index, key in enumerate(keys, start=1):
        try:
            status = process_setting(
                key, groups.get((key[0], key[1], key[2], key[4]), {}),
                catalog, args.t7_root, args.bigtable_dir,
                args.table_dir, args.pdf_dir, args.no_pdf,
            )
        except Exception as error:
            phase, pass_name, run_type, target, setting = key
            identity = {
                "Phase": phase, "Pass": pass_name, "Run_type": run_type,
                "Target": target, "Setting": setting,
            }
            reason = f"SETTING_PROCESSING_FAILED={type(error).__name__}:{error}"
            atomic_csv(
                args.table_dir / ("_".join(identity.values()) + ".csv"),
                [metadata_row(identity, "ERROR", reason, S_DUMMY[(phase, target)])],
            )
            status = "ERROR"
            print(f"ERROR detail: {reason}")
        counts[status] += 1
        print(f"[{index}/{len(keys)}] {status}: {'_'.join(key)}")
    print("Status:", ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 2 if counts["ERROR"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
