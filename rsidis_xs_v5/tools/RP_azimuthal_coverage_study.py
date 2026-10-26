#!/usr/bin/env python3
"""Study SIDIS azimuthal coverage after combining nominal-thpq settings."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from array import array
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import ROOT

import RP_data_mc_compare as dmc
from RP_sidis_binning import (
    AZIMUTHAL_BINNING_VERSION,
    AZIMUTHAL_PHI_EDGES as PHI_EDGES,
    AZIMUTHAL_PT_EDGES as PT_EDGES,
)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

IDENTITY = ("Phase", "Pass", "Run_type", "Target", "setting_id")
POLICY_VERSION = "RP_AZIMUTHAL_COVERAGE_V6_THPQ_COMBINED"
DEFAULT_MODEL_CATALOG = (
    Path(__file__).resolve().parents[1]
    / "results/Pass1/Tables/RP_sidis_model/RP_sidis_model.csv"
)
DEFAULT_TABLE_DIR = (
    Path(__file__).resolve().parents[1]
    / "results/Pass1/Tables/RP_azimuthal_coverage_study"
)
DEFAULT_PDF_DIR = (
    Path(__file__).resolve().parents[1]
    / "results/Pass1/PDFs/RP_azimuthal_coverage_study"
)
DEFAULT_COMBINED_YIELDS = (
    Path(__file__).resolve().parents[1]
    / "results/Pass1/Tables/RP_combine_thpq_yields/RP_combined_thpq_yields.csv"
)

BIN_FIELDS = [
    *IDENTITY, "x", "Q2", "thpq", "z", "Pt_bin", "Pt_low", "Pt_high",
    "Pt_high_closed", "phipq_bin", "phipq_low", "phipq_high", "phipq_center",
    "Binning_version", "Policy_version", "Max_MC_rel_err", "Max_Ratio_rel_err",
    "Ratio_warning_low", "Ratio_warning_high", "Y_Data", "Y_Data_err",
    "Data_rel_err", *(item for reaction in dmc.REACTIONS for item in
                      (f"Y_MC_{reaction}", f"Y_MC_{reaction}_err")),
    "Y_MC", "Y_MC_err", "MC_rel_err", "Data_by_MC", "Data_by_MC_err",
    "Ratio_rel_err", "MC_support_pass", "Ratio_precision_pass",
    "Independent_pass", "Selected_contiguous", "Bin_status", "Bin_reason",
]

CURVE_FIELDS = [
    *IDENTITY, "x", "Q2", "thpq", "z", "Pt_bin", "Pt_low", "Pt_high",
    "Pt_high_closed", "Binning_version", "Policy_version", "N_phi_bins_total",
    "Independent_mask", "Selected_mask", "Selected_start_bin", "Selected_end_bin",
    "Selected_wraps_zero", "N_independent_pass", "N_selected", "Selected_span_rad",
    "Selected_span_fraction", "Design_rank", "Design_condition", "Weighted_design_condition",
    "Coeff_sigma_0", "Coeff_sigma_1", "Coeff_sigma_2",
    "Coeff_correlation_01", "Coeff_correlation_02", "Coeff_correlation_12",
    "Min_selected_span_rad", "Curve_status", "Curve_reason",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def atomic_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def make_hist(name: str) -> ROOT.TH2D:
    hist = ROOT.TH2D(
        name, ";#phi_{pq} [rad];P_{t} [GeV];Normalized yield",
        len(PHI_EDGES) - 1, array("d", PHI_EDGES),
        len(PT_EDGES) - 1, array("d", PT_EDGES),
    )
    hist.Sumw2()
    hist.SetDirectory(0)
    return hist


def dataframe_hist(path: Path, tree: str, phi: str, pt: str, selection: str,
                   weight: str, name: str) -> ROOT.TH2D:
    frame = (
        ROOT.RDataFrame(tree, str(path)).Filter(selection)
        .Define("_rp_cov_phi", f"fmod(fmod((double)({phi}),2*TMath::Pi())+2*TMath::Pi(),2*TMath::Pi())")
        .Define("_rp_cov_pt", pt).Define("_rp_cov_weight", weight)
        .Filter("std::isfinite((double)_rp_cov_phi) && std::isfinite((double)_rp_cov_pt) && _rp_cov_pt>=0 && _rp_cov_pt<=0.5")
    )
    action = frame.Histo2D(
        (name, "", len(PHI_EDGES) - 1, array("d", PHI_EDGES),
         len(PT_EDGES) - 1, array("d", PT_EDGES)),
        "_rp_cov_phi", "_rp_cov_pt", "_rp_cov_weight",
    )
    hist = action.GetValue().Clone(name + "_owned")
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def data_class(path: Path, root_dir: Path, prefix: str) -> tuple[ROOT.TH2D, list[str]]:
    rows = [
        row for row in dmc.read_csv(path)
        if row.get("Norm_status") == "OK" and finite(row.get("BCM2_Q"))
        and float(row["BCM2_Q"]) > 0 and finite(row.get("RP_Normyield"))
    ]
    output = make_hist(prefix)
    if not rows:
        return output, ["NO_VALID_NORMALIZED_RUNS"]
    total_charge = sum(float(row["BCM2_Q"]) for row in rows)
    reasons: list[str] = []
    for row in rows:
        run = int(float(row["Run"]))
        if dmc.is_floored_positron_run(path, row):
            continue
        root = dmc.data_root_path(root_dir, dmc.parse_norm_stem(path)["Phase"], run)
        if not root.is_file():
            reasons.append(f"DATA_ROOT_MISSING_RUN_{run}")
            continue
        timing = dmc.timing_weight_expression(row)
        selection = f"({dmc.PID_CUT}) && ({dmc.DATA_DELTA}) && ({timing})!=0"
        hist = dataframe_hist(
            root, "T", "P_kin_secondary_ph_xq", "abs(pt)", selection,
            f"({float(row['Norm_factor']):.17g})*({timing})", f"{prefix}_{run}",
        )
        output.Add(hist, float(row["BCM2_Q"]) / total_charge)
    return output, reasons


def simc_class(row: Mapping[str, str], prefix: str) -> tuple[ROOT.TH2D, list[str]]:
    output = make_hist(prefix)
    path = Path(row.get("Recon_root_path", ""))
    if row.get("QA_status") not in {"OK", "WARNING"}:
        return output, [f"SIMC_QA_{row.get('QA_status', 'MISSING')}"]
    if not path.is_file():
        return output, ["SIMC_ROOT_MISSING"]
    expressions = dmc.simc_expressions(row["Reaction"])
    try:
        return dataframe_hist(
            path, "h10", expressions["phipq"], expressions["pt"], dmc.SIMC_DELTA,
            "static_cast<double>(fWeight)", prefix,
        ), []
    except Exception as error:
        return output, [f"SIMC_HIST_FAILED={type(error).__name__}:{error}"]


def combine_data(classes: Mapping[str, ROOT.TH2D], scale: float) -> ROOT.TH2D:
    result = classes["target_e"].Clone("coverage_data")
    result.Add(classes["target_p"], -1)
    dummy = classes["dummy_e"].Clone("coverage_dummy")
    dummy.Add(classes["dummy_p"], -1)
    result.Add(dummy, -scale)
    result.SetDirectory(0)
    return result


def longest_circular_run(mask: Sequence[bool]) -> list[int]:
    """Return zero-based indices in the longest passing circular run."""
    size = len(mask)
    if not size or not any(mask):
        return []
    if all(mask):
        return list(range(size))
    best: list[int] = []
    for start in range(size):
        if not mask[start] or mask[(start - 1) % size]:
            continue
        run = []
        index = start
        while mask[index] and len(run) < size:
            run.append(index)
            index = (index + 1) % size
        if len(run) > len(best):
            best = run
    return best


def design_diagnostics(indices: Sequence[int]) -> tuple[int, float]:
    if not indices:
        return 0, math.inf
    phi = np.asarray([(PHI_EDGES[index] + PHI_EDGES[index + 1]) / 2 for index in indices])
    design = np.column_stack((np.ones(len(phi)), np.cos(phi), np.cos(2 * phi)))
    rank = int(np.linalg.matrix_rank(design))
    condition = float(np.linalg.cond(design)) if rank == 3 else math.inf
    return rank, condition


def classify_bin(data: float, data_error: float, mc: float, mc_error: float,
                 structural_reasons: Sequence[str], advisory_reasons: Sequence[str], max_mc_rel: float,
                 max_ratio_rel: float, ratio_low: float, ratio_high: float) -> dict[str, object]:
    reasons = list(structural_reasons) + list(advisory_reasons)
    mc_rel = mc_error / mc if finite(mc) and finite(mc_error) and mc > 0 else math.nan
    data_rel = data_error / abs(data) if finite(data) and finite(data_error) and data != 0 else math.inf
    ratio = data / mc if finite(data) and finite(mc) and mc > 0 else math.nan
    ratio_error = (
        math.hypot(data_error / mc, data * mc_error / (mc * mc))
        if finite(ratio) and finite(data_error) and finite(mc_error) else math.nan
    )
    ratio_rel = ratio_error / abs(ratio) if finite(ratio) and ratio != 0 else math.inf
    structural_ok = not structural_reasons
    mc_pass = structural_ok and finite(mc) and mc > 0 and finite(mc_rel) and mc_rel <= max_mc_rel
    ratio_pass = (finite(data) and data > 0 and finite(data_error) and data_error >= 0
                  and finite(ratio_rel) and ratio_rel <= max_ratio_rel)
    if not finite(mc) or mc <= 0:
        reasons.append("NONPOSITIVE_MC")
    elif not finite(mc_rel) or mc_rel > max_mc_rel:
        reasons.append("MC_REL_ERR_ABOVE_LIMIT")
    if not finite(data):
        reasons.append("NONFINITE_CORRECTED_DATA")
    elif data <= 0:
        reasons.append("NONPOSITIVE_CORRECTED_DATA")
    if not finite(data_error) or data_error < 0:
        reasons.append("INVALID_DATA_UNCERTAINTY")
    elif not finite(ratio_rel) or ratio_rel > max_ratio_rel:
        reasons.append("RATIO_REL_ERR_ABOVE_LIMIT")
    warning = finite(ratio) and not (ratio_low <= ratio <= ratio_high)
    if warning:
        reasons.append("RATIO_OUTSIDE_WARNING_RANGE")
    independent = mc_pass and ratio_pass
    has_warning = warning or bool(advisory_reasons)
    status = "PASS" if independent and not has_warning else "WARNING" if independent else "FAIL"
    return {
        "Data_rel_err": data_rel, "MC_rel_err": mc_rel, "Data_by_MC": ratio,
        "Data_by_MC_err": ratio_error, "Ratio_rel_err": ratio_rel,
        "MC_support_pass": int(mc_pass), "Ratio_precision_pass": int(ratio_pass),
        "Independent_pass": int(independent), "Bin_status": status,
        "Bin_reason": ";".join(sorted(set(reasons))),
    }


def finalize_curve(rows: list[dict[str, object]], min_span: float,
                   max_condition: float) -> dict[str, object]:
    mask = [bool(int(row["Independent_pass"])) for row in rows]
    selected = longest_circular_run(mask)
    selected_set = set(selected)
    for index, row in enumerate(rows):
        row["Selected_contiguous"] = int(index in selected_set)
    rank, condition = design_diagnostics(selected)
    weighted_condition=math.inf; coefficient_sigma=[math.nan]*3; coefficient_correlation=[math.nan]*3
    if rank==3 and selected:
        phi=np.asarray([float(rows[index]["phipq_center"]) for index in selected])
        errors=np.asarray([float(rows[index]["Data_by_MC_err"]) for index in selected])
        design=np.column_stack((np.ones(len(phi)),np.cos(phi),np.cos(2*phi)))
        weighted=design/errors[:,None]; covariance=np.linalg.inv(weighted.T@weighted)
        weighted_condition=float(np.linalg.cond(weighted))
        coefficient_sigma=np.sqrt(np.maximum(0,np.diag(covariance))).tolist()
        coefficient_correlation=[covariance[i,j]/math.sqrt(covariance[i,i]*covariance[j,j]) for i,j in ((0,1),(0,2),(1,2))]
    span = len(selected) * 2 * math.pi / len(mask)
    reasons = []
    if span < min_span:
        reasons.append("INSUFFICIENT_ANGULAR_SPAN")
    if rank < 3:
        reasons.append("DESIGN_RANK_BELOW_3")
    elif not finite(weighted_condition) or weighted_condition > max_condition:
        reasons.append("WEIGHTED_DESIGN_CONDITION_ABOVE_LIMIT")
    first = rows[0]
    curve = {field: "" for field in CURVE_FIELDS}
    curve.update({field: first[field] for field in (*IDENTITY, "x", "Q2", "thpq", "z", "Pt_bin", "Pt_low", "Pt_high", "Pt_high_closed", "Binning_version", "Policy_version")})
    curve.update({
        "N_phi_bins_total": len(mask),
        "Independent_mask": "".join("1" if value else "0" for value in mask),
        "Selected_mask": "".join("1" if index in selected_set else "0" for index in range(len(mask))),
        "Selected_start_bin": selected[0] + 1 if selected else "",
        "Selected_end_bin": selected[-1] + 1 if selected else "",
        "Selected_wraps_zero": int(bool(selected) and selected[-1] < selected[0]),
        "N_independent_pass": sum(mask), "N_selected": len(selected),
        "Selected_span_rad": span, "Selected_span_fraction": len(selected) / len(mask),
        "Design_rank": rank, "Design_condition": condition,
        "Weighted_design_condition": weighted_condition,
        "Coeff_sigma_0": coefficient_sigma[0], "Coeff_sigma_1": coefficient_sigma[1], "Coeff_sigma_2": coefficient_sigma[2],
        "Coeff_correlation_01": coefficient_correlation[0], "Coeff_correlation_02": coefficient_correlation[1], "Coeff_correlation_12": coefficient_correlation[2],
        "Min_selected_span_rad": min_span,
        "Curve_status": "APPROVED" if not reasons else "REJECTED",
        "Curve_reason": ";".join(reasons),
    })
    return curve


def extract_setting(key, files, catalog, model, root_dir, options):
    phase, pass_name, run_type, target, setting = key
    categories = {
        "target_e": (target, "Elec"), "target_p": (target, "Pos"),
        "dummy_e": ("DUMMY", "Elec"), "dummy_p": ("DUMMY", "Pos"),
    }
    structural: list[str] = []
    advisory: list[str] = []
    classes = {}
    for name, category in categories.items():
        if category in files:
            classes[name], reasons = data_class(files[category], root_dir, f"{setting}_{name}")
            structural.extend(reasons)
        else:
            classes[name] = make_hist(f"{setting}_{name}_zero")
            advisory.append(f"ZERO_CORRECTION_CLASS={name}")
    data = combine_data(classes, dmc.S_DUMMY[(phase, target)])
    components = {}
    physical = dmc.applicable_reactions(run_type)
    for reaction in dmc.REACTIONS:
        if reaction not in physical:
            components[reaction] = make_hist(f"{setting}_{reaction}_zero")
            continue
        source = catalog.get((*key, reaction))
        if source is None:
            components[reaction] = make_hist(f"{setting}_{reaction}_missing")
            structural.append(f"MISSING_MC={reaction}")
        else:
            components[reaction], reasons = simc_class(source, f"{setting}_{reaction}")
            structural.extend(reasons)
    rows = []
    curves = []
    for ipt, (pt_low, pt_high) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:]), 1):
        curve_rows = []
        for iphi, (phi_low, phi_high) in enumerate(zip(PHI_EDGES[:-1], PHI_EDGES[1:]), 1):
            data_value = data.GetBinContent(iphi, ipt)
            data_error = data.GetBinError(iphi, ipt)
            mc_values = [components[r].GetBinContent(iphi, ipt) for r in dmc.REACTIONS]
            mc_errors = [components[r].GetBinError(iphi, ipt) for r in dmc.REACTIONS]
            mc = sum(mc_values)
            mc_error = math.sqrt(sum(value * value for value in mc_errors))
            decision = classify_bin(
                data_value, data_error, mc, mc_error, structural, advisory,
                options.max_mc_rel_err, options.max_ratio_rel_err,
                options.ratio_warning_low, options.ratio_warning_high,
            )
            row = {field: "" for field in BIN_FIELDS}
            row.update({
                "Phase": phase, "Pass": pass_name, "Run_type": run_type,
                "Target": target, "setting_id": setting, "x": model["x"],
                "Q2": model["Q2_GeV2"], "thpq": model["theta_pq_deg"], "z": model["z"],
                "Pt_bin": ipt, "Pt_low": pt_low, "Pt_high": pt_high,
                "Pt_high_closed": int(ipt == len(PT_EDGES) - 1),
                "phipq_bin": iphi, "phipq_low": phi_low, "phipq_high": phi_high,
                "phipq_center": (phi_low + phi_high) / 2,
                "Binning_version": AZIMUTHAL_BINNING_VERSION,
                "Policy_version": POLICY_VERSION,
                "Max_MC_rel_err": options.max_mc_rel_err,
                "Max_Ratio_rel_err": options.max_ratio_rel_err,
                "Ratio_warning_low": options.ratio_warning_low,
                "Ratio_warning_high": options.ratio_warning_high,
                "Y_Data": data_value, "Y_Data_err": data_error,
                "Y_MC": mc, "Y_MC_err": mc_error, **decision,
            })
            for reaction, value, error in zip(dmc.REACTIONS, mc_values, mc_errors):
                row[f"Y_MC_{reaction}"] = value
                row[f"Y_MC_{reaction}_err"] = error
            curve_rows.append(row)
        curve = finalize_curve(curve_rows, options.min_selected_span,
                               options.max_design_condition)
        rows.extend(curve_rows)
        curves.append(curve)
    return rows, curves, data, components


def z_directory(value: object) -> str:
    return "z_" + f"{float(value):.2f}".replace(".", "p")


def draw_setting_pdf(path: Path, rows: Sequence[Mapping[str, object]], data: ROOT.TH2D,
                     components: Mapping[str, ROOT.TH2D]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    total = components[dmc.REACTIONS[0]].Clone("coverage_mc_total_draw")
    total.Reset("ICES")
    for reaction in dmc.REACTIONS:
        total.Add(components[reaction])
    ratio = data.Clone("coverage_ratio_draw")
    ratio.Reset("ICES")
    for ix in range(1, data.GetNbinsX() + 1):
        for iy in range(1, data.GetNbinsY() + 1):
            denominator = total.GetBinContent(ix, iy)
            if denominator > 0:
                ratio.SetBinContent(ix, iy, data.GetBinContent(ix, iy) / denominator)
    canvas = ROOT.TCanvas("c_coverage", "coverage", 1800, 1200)
    canvas.Print(str(path) + "[")
    map_objects=[]
    for index, (hist, title) in enumerate(((data, "Corrected data"), (total, "MC total"), (ratio, "Data / MC")), 1):
        x1=.015+(index-1)*.33; pad=ROOT.TPad(f"coverage_map_pad_{index}","",x1,.12,x1+.32,.90)
        pad.SetLeftMargin(.14); pad.SetRightMargin(.18); pad.SetBottomMargin(.13); pad.SetTopMargin(.08); pad.Draw(); pad.cd()
        if index == 3:
            hist.SetMinimum(0.0); hist.SetMaximum(2.0)
            hist.GetZaxis().SetTitle("Data / MC")
        hist.GetXaxis().SetTitleSize(.05); hist.GetYaxis().SetTitleSize(.05); hist.GetZaxis().SetTitleSize(.045)
        hist.GetXaxis().SetLabelSize(.04); hist.GetYaxis().SetLabelSize(.04); hist.GetZaxis().SetLabelSize(.035)
        hist.GetXaxis().SetTitleOffset(1.05); hist.GetYaxis().SetTitleOffset(1.15)
        hist.Draw("COLZ")
        label = ROOT.TLatex(); label.SetNDC(); label.SetTextSize(.045)
        label.DrawLatex(.16, .94, title); map_objects.extend((pad,label))
        canvas.cd()
    canvas._coverage_map_labels=map_objects
    canvas.Print(str(path))
    canvas.Clear()
    canvas.Divide(3, 2)
    owned = [total, ratio]
    for ipt in range(1, len(PT_EDGES)):
        pad = canvas.cd(ipt)
        pad.SetGridx(); pad.SetGridy(); pad.SetLeftMargin(.15); pad.SetRightMargin(.05); pad.SetBottomMargin(.15); pad.SetTopMargin(.10)
        selected_rows = [row for row in rows if int(row["Pt_bin"]) == ipt]
        frame = ROOT.TH1D(
            f"coverage_frame_{ipt}", "", len(PHI_EDGES) - 1, 0, 2 * math.pi
        )
        values = [float(row["Data_by_MC"]) for row in selected_rows if finite(row["Data_by_MC"])]
        errors = [float(row["Data_by_MC_err"]) for row in selected_rows if finite(row["Data_by_MC_err"])]
        ymax = max([2.1] + [value + error for value, error in zip(values, errors)])
        frame.SetMinimum(0); frame.SetMaximum(ymax * 1.1)
        frame.GetXaxis().SetTitle("#phi_{pq} [rad]")
        frame.GetYaxis().SetTitle("Data / MC")
        frame.GetXaxis().SetTitleSize(.05); frame.GetYaxis().SetTitleSize(.05); frame.GetXaxis().SetLabelSize(.04); frame.GetYaxis().SetLabelSize(.04); frame.GetYaxis().SetTitleOffset(1.25)
        frame.Draw()
        accepted = ROOT.TGraphErrors(); rejected = ROOT.TGraphErrors()
        for row in selected_rows:
            if not finite(row["Data_by_MC"]):
                continue
            graph = accepted if int(row["Selected_contiguous"]) else rejected
            point = graph.GetN()
            graph.SetPoint(point, float(row["phipq_center"]), float(row["Data_by_MC"]))
            graph.SetPointError(point, 0, float(row["Data_by_MC_err"]) if finite(row["Data_by_MC_err"]) else 0)
        accepted.SetMarkerStyle(20); accepted.SetMarkerColor(ROOT.kGreen + 2)
        rejected.SetMarkerStyle(24); rejected.SetMarkerColor(ROOT.kRed + 1)
        accepted.Draw("P SAME"); rejected.Draw("P SAME")
        label = ROOT.TLatex(); label.SetNDC(); label.SetTextSize(.045)
        label.DrawLatex(.16, .92, f"{PT_EDGES[ipt-1]:.1f} #leq P_{{t}} {'#leq' if ipt == len(PT_EDGES)-1 else '<'} {PT_EDGES[ipt]:.1f} GeV")
        owned.extend((frame, accepted, rejected, label))
    pad = canvas.cd(6)
    pad.Clear()
    text = ROOT.TLatex(); text.SetNDC(); text.SetTextSize(.050)
    first = rows[0]
    text.DrawLatex(.08, .88, f"{first['Phase']}, {first['Pass']}, {first['Run_type']}, {first['Target']}")
    text.DrawLatex(.08, .78, f"setting={first['setting_id']}, z={float(first['z']):g}")
    text.DrawLatex(.08, .65, "green: selected contiguous range")
    text.DrawLatex(.08, .57, "red open: not selected")
    owned.append(text)
    canvas._coverage_owned = owned
    canvas.Print(str(path))
    canvas.Print(str(path) + "]")
    canvas.Close()


COMBINED_IDENTITY = (
    "Combination_group_id", "Phase", "Pass", "Run_type", "Target", "x", "Q2", "z"
)
COMBINED_BIN_FIELDS = [
    *COMBINED_IDENTITY, "Ebeam_GeV", "N_contributors", "Contributing_setting_ids",
    "Contributing_thpq_deg", "Combination_policy_version", "Pt_bin", "Pt_low",
    "Pt_high", "Pt_high_closed", "phipq_bin", "phipq_low", "phipq_high",
    "phipq_center", "Binning_version", "Policy_version", "Max_MC_rel_err",
    "Max_Ratio_rel_err", "Ratio_warning_low", "Ratio_warning_high", "Y_Data",
    "Y_Data_err", "Data_rel_err", *(item for reaction in dmc.REACTIONS for item in
    (f"Y_MC_{reaction}", f"Y_MC_{reaction}_err")), "Y_MC", "Y_MC_err",
    "MC_rel_err", "Data_by_MC", "Data_by_MC_err", "Ratio_rel_err",
    "MC_support_pass", "Ratio_precision_pass", "Independent_pass",
    "Selected_contiguous", "Survives_all_requirements", "Bin_status", "Bin_reason",
]
COMBINED_CURVE_FIELDS = [
    *COMBINED_IDENTITY, "Ebeam_GeV", "N_contributors", "Contributing_setting_ids",
    "Contributing_thpq_deg", "Combination_policy_version", "Pt_bin", "Pt_low",
    "Pt_high", "Pt_high_closed", "Binning_version", "Policy_version",
    "N_phi_bins_total", "Independent_mask", "Selected_mask", "Surviving_mask",
    "Selected_start_bin", "Selected_end_bin", "Selected_wraps_zero",
    "N_independent_pass", "N_selected", "N_surviving", "Selected_span_rad",
    "Selected_span_fraction", "Design_rank", "Design_condition",
    "Weighted_design_condition", "Coeff_sigma_0", "Coeff_sigma_1", "Coeff_sigma_2",
    "Coeff_correlation_01", "Coeff_correlation_02", "Coeff_correlation_12",
    "Min_selected_span_rad", "Curve_status", "Curve_reason",
]


def finalize_combined_curve(rows, min_span, max_condition):
    mask = [bool(int(row["Independent_pass"])) for row in rows]
    selected = longest_circular_run(mask); selected_set = set(selected)
    phi = np.asarray([float(rows[index]["phipq_center"]) for index in selected])
    rank, condition = design_diagnostics(selected)
    weighted_condition = math.inf; sigma = [math.nan] * 3; correlation = [math.nan] * 3
    if rank == 3 and selected:
        errors = np.asarray([float(rows[index]["Data_by_MC_err"]) for index in selected])
        design = np.column_stack((np.ones(len(phi)), np.cos(phi), np.cos(2 * phi)))
        weighted = design / errors[:, None]; covariance = np.linalg.inv(weighted.T @ weighted)
        weighted_condition = float(np.linalg.cond(weighted))
        sigma = np.sqrt(np.maximum(0, np.diag(covariance))).tolist()
        correlation = [covariance[i,j]/math.sqrt(covariance[i,i]*covariance[j,j])
                       for i,j in ((0,1),(0,2),(1,2))]
    span = len(selected) * 2 * math.pi / len(mask); reasons = []
    if span < min_span: reasons.append("INSUFFICIENT_ANGULAR_SPAN")
    if rank < 3: reasons.append("DESIGN_RANK_BELOW_3")
    elif not finite(weighted_condition) or weighted_condition > max_condition:
        reasons.append("WEIGHTED_DESIGN_CONDITION_ABOVE_LIMIT")
    approved = not reasons
    for index, row in enumerate(rows):
        row["Selected_contiguous"] = int(index in selected_set)
        row["Survives_all_requirements"] = int(approved and index in selected_set)
    first = rows[0]; curve = {field: "" for field in COMBINED_CURVE_FIELDS}
    copy = (*COMBINED_IDENTITY, "Ebeam_GeV", "N_contributors",
            "Contributing_setting_ids", "Contributing_thpq_deg",
            "Combination_policy_version", "Pt_bin", "Pt_low", "Pt_high",
            "Pt_high_closed", "Binning_version", "Policy_version")
    curve.update({field: first[field] for field in copy})
    curve.update({
        "N_phi_bins_total": len(mask), "Independent_mask": "".join("1" if v else "0" for v in mask),
        "Selected_mask": "".join("1" if i in selected_set else "0" for i in range(len(mask))),
        "Surviving_mask": "".join("1" if approved and i in selected_set else "0" for i in range(len(mask))),
        "Selected_start_bin": selected[0]+1 if selected else "",
        "Selected_end_bin": selected[-1]+1 if selected else "",
        "Selected_wraps_zero": int(bool(selected) and selected[-1] < selected[0]),
        "N_independent_pass": sum(mask), "N_selected": len(selected),
        "N_surviving": len(selected) if approved else 0, "Selected_span_rad": span,
        "Selected_span_fraction": len(selected)/len(mask), "Design_rank": rank,
        "Design_condition": condition, "Weighted_design_condition": weighted_condition,
        "Coeff_sigma_0": sigma[0], "Coeff_sigma_1": sigma[1], "Coeff_sigma_2": sigma[2],
        "Coeff_correlation_01": correlation[0], "Coeff_correlation_02": correlation[1],
        "Coeff_correlation_12": correlation[2], "Min_selected_span_rad": min_span,
        "Curve_status": "APPROVED" if approved else "REJECTED",
        "Curve_reason": ";".join(reasons),
    })
    return curve


def evaluate_combined_group(source_rows, options):
    rows, curves = [], []
    ordered = sorted(source_rows, key=lambda row: (int(row["Pt_bin"]), int(row["phipq_bin"])))
    for ipt in range(1, len(PT_EDGES)):
        curve_rows = []
        for source in [row for row in ordered if int(row["Pt_bin"]) == ipt]:
            structural = []
            advisory = []
            reason = source.get("Combination_reason", "")
            if source.get("Combination_status") == "ERROR": structural = list(filter(None, reason.split(";")))
            elif source.get("Combination_status") == "WARNING": advisory = list(filter(None, reason.split(";")))
            decision = classify_bin(
                float(source["Y_Data"]), float(source["Y_Data_err"]),
                float(source["Y_MC"]), float(source["Y_MC_err"]), structural, advisory,
                options.max_mc_rel_err, options.max_ratio_rel_err,
                options.ratio_warning_low, options.ratio_warning_high,
            )
            row = {field: "" for field in COMBINED_BIN_FIELDS}
            for field in row:
                if field in source: row[field] = source[field]
            row.update(decision, Policy_version=POLICY_VERSION,
                       Max_MC_rel_err=options.max_mc_rel_err,
                       Max_Ratio_rel_err=options.max_ratio_rel_err,
                       Ratio_warning_low=options.ratio_warning_low,
                       Ratio_warning_high=options.ratio_warning_high,
                       Selected_contiguous=0, Survives_all_requirements=0)
            curve_rows.append(row)
        if len(curve_rows) != len(PHI_EDGES)-1:
            raise ValueError(f"{ordered[0]['Combination_group_id']} Pt bin {ipt} has {len(curve_rows)} phi rows")
        curve = finalize_combined_curve(curve_rows, options.min_selected_span,
                                        options.max_design_condition)
        rows.extend(curve_rows); curves.append(curve)
    return rows, curves


def blank_ratio_hist(name):
    hist = make_hist(name)
    for ix in range(1, hist.GetNbinsX()+1):
        # The polar painter colors NaN cells at the palette minimum. Exact
        # zeros remain empty, so rejected or unsupported bins are not drawn.
        for iy in range(1, hist.GetNbinsY()+1): hist.SetBinContent(ix, iy, 0.0)
    hist.GetZaxis().SetTitle("Data / MC"); hist.SetMinimum(0.0); hist.SetMaximum(2.0)
    # POL otherwise retains Cartesian phi/Pt labels along the bottom and left,
    # which are misleading for the polar representation.
    for axis in (hist.GetXaxis(),hist.GetYaxis()):
        axis.SetTitle(""); axis.SetLabelSize(0); axis.SetTickLength(0)
    return hist


def wrapped_items(prefix, values, width=95):
    lines=[]; current=prefix
    for value in values:
        piece=(", " if current != prefix else "")+value
        if len(current)+len(piece)>width:
            lines.append(current); current="  "+value
        else: current+=piece
    if current.strip(): lines.append(current)
    return lines


def phi_projection(rows, accepted_only, prefix):
    """Project the five Pt bins onto phi, optionally retaining survivors only."""
    data=ROOT.TH1D(prefix+"_data","",len(PHI_EDGES)-1,array("d",PHI_EDGES))
    mc=ROOT.TH1D(prefix+"_mc","",len(PHI_EDGES)-1,array("d",PHI_EDGES))
    for hist in (data,mc): hist.Sumw2(); hist.SetDirectory(0)
    for iphi in range(1,len(PHI_EDGES)):
        selected=[row for row in rows if int(row["phipq_bin"])==iphi
                  and (not accepted_only or int(row["Survives_all_requirements"]))]
        if accepted_only and not selected: continue
        for hist,value,error in (
            (data,"Y_Data","Y_Data_err"),(mc,"Y_MC","Y_MC_err")):
            valid=[row for row in selected if finite(row[value]) and finite(row[error])]
            if not valid: continue
            hist.SetBinContent(iphi,sum(float(row[value]) for row in valid))
            hist.SetBinError(iphi,math.sqrt(sum(float(row[error])**2 for row in valid)))
    return data,mc


def draw_phi_comparison_page(canvas,path,rows):
    """Draw before/after combined-thpq phi projections with ratio subpads."""
    canvas.Clear(); canvas.cd(); owned=[]
    for column,(title,accepted_only) in enumerate((("Before coverage selection",False),
                                                   ("After all coverage requirements",True))):
        canvas.cd()
        x1=.01+column*.50; x2=.49+column*.50
        top=ROOT.TPad(f"phi_projection_top_{column}","",x1,.32,x2,.97)
        bottom=ROOT.TPad(f"phi_projection_bottom_{column}","",x1,.03,x2,.32)
        for pad in (top,bottom): pad.Draw(); owned.append(pad)
        data,mc=phi_projection(rows,accepted_only,f"phi_projection_{column}")
        values=[]
        for hist in (data,mc):
            values.extend(hist.GetBinContent(i)+hist.GetBinError(i)
                          for i in range(1,hist.GetNbinsX()+1)
                          if finite(hist.GetBinContent(i)))
        data_lows=[data.GetBinContent(i)-data.GetBinError(i)
                   for i in range(1,data.GetNbinsX()+1) if finite(data.GetBinContent(i))]
        ymax=max([1e-12]+values); ymin=min([0.0]+data_lows)
        span=max(ymax-ymin,1e-12)
        top.cd(); top.SetGrid(); top.SetLeftMargin(.15); top.SetRightMargin(.04)
        top.SetBottomMargin(.02); top.SetTopMargin(.12)
        data.SetTitle(f"{title};#phi_{{pq}} [rad];Normalized yield")
        data.SetMinimum(min(0.0,ymin-.08*span)); data.SetMaximum(ymax+.18*span)
        data.GetXaxis().SetLabelSize(0); data.GetXaxis().SetTitleSize(0)
        data.GetYaxis().SetTitleSize(.052); data.GetYaxis().SetLabelSize(.043)
        data.GetYaxis().SetTitleOffset(1.35)
        data.SetMarkerStyle(20); data.SetMarkerColor(ROOT.kBlue+1); data.SetLineColor(ROOT.kBlue+1)
        data.Draw("E1")
        mc.SetLineColor(ROOT.kBlack); mc.SetLineWidth(3); mc.SetFillColorAlpha(ROOT.kGray+1,.25)
        mc.Draw("E2 SAME"); mc.Draw("HIST SAME"); data.Draw("E1 SAME")
        legend=ROOT.TLegend(.58,.72,.92,.88); legend.SetBorderSize(0); legend.SetFillStyle(0)
        legend.AddEntry(data,"Corrected Data","lep"); legend.AddEntry(mc,"MC Total","lf"); legend.Draw()
        bottom.cd(); bottom.SetGrid(); bottom.SetLeftMargin(.15); bottom.SetRightMargin(.04)
        bottom.SetTopMargin(.03); bottom.SetBottomMargin(.30)
        frame=data.Clone(f"phi_projection_ratio_frame_{column}"); frame.Reset("ICES")
        frame.SetTitle(";#phi_{pq} [rad];Data / MC"); frame.SetMinimum(0); frame.SetMaximum(2)
        frame.GetXaxis().SetTitleSize(.12); frame.GetXaxis().SetLabelSize(.095)
        frame.GetYaxis().SetTitleSize(.105); frame.GetYaxis().SetLabelSize(.085)
        frame.GetYaxis().SetTitleOffset(.62); frame.Draw()
        ratio=ROOT.TGraphErrors()
        for ibin in range(1,data.GetNbinsX()+1):
            dv,de=data.GetBinContent(ibin),data.GetBinError(ibin)
            mv,me=mc.GetBinContent(ibin),mc.GetBinError(ibin)
            if not finite(mv) or mv<=0: continue
            point=ratio.GetN(); ratio.SetPoint(point,data.GetBinCenter(ibin),dv/mv)
            ratio.SetPointError(point,0,math.hypot(de/mv,dv*me/(mv*mv)))
        ratio.SetMarkerStyle(20); ratio.SetMarkerColor(ROOT.kBlack); ratio.SetLineColor(ROOT.kBlack)
        ratio.Draw("P SAME")
        unity=ROOT.TLine(PHI_EDGES[0],1,PHI_EDGES[-1],1)
        unity.SetLineColor(ROOT.kRed+1); unity.SetLineStyle(2); unity.SetLineWidth(2); unity.Draw("SAME")
        owned.extend((data,mc,legend,frame,ratio,unity))
    canvas._phi_projection_owned=owned; canvas.Print(str(path))


def draw_combined_pdf(path, rows, curves, options):
    path.parent.mkdir(parents=True, exist_ok=True)
    before, after = blank_ratio_hist("coverage_before"), blank_ratio_hist("coverage_after")
    for row in rows:
        ix, iy = int(row["phipq_bin"]), int(row["Pt_bin"])
        if finite(row["Data_by_MC"]) and float(row["Y_MC"]) > 0:
            before.SetBinContent(ix, iy, float(row["Data_by_MC"]))
        if int(row["Survives_all_requirements"]):
            after.SetBinContent(ix, iy, float(row["Data_by_MC"]))
    first=rows[0]; canvas=ROOT.TCanvas("c_combined_coverage","combined coverage",1800,1200)
    canvas.Print(str(path)+"["); canvas.Clear(); owned=[before,after]
    for index,(hist,title) in enumerate(((before,"Before coverage selection"),(after,"After all coverage requirements"))):
        x1=.02+index*.49; pad=ROOT.TPad(f"combined_polar_{index}","",x1,.23,x1+.465,.90)
        pad.SetLeftMargin(.06); pad.SetRightMargin(.16); pad.SetTopMargin(.05); pad.SetBottomMargin(.06)
        pad.Draw(); pad.cd(); hist.Draw("COLZ POL")
        angular=ROOT.TLatex(); angular.SetNDC(); angular.SetTextFont(42); angular.SetTextSize(.038); angular.SetTextAlign(22)
        angular.DrawLatex(.80,.50,"0#circ"); angular.DrawLatex(.45,.86,"90#circ")
        angular.DrawLatex(.12,.50,"180#circ"); angular.DrawLatex(.45,.14,"270#circ")
        owned.append(angular)
        if index==0:
            radial=ROOT.TLatex(); radial.SetNDC(); radial.SetTextFont(42); radial.SetTextSize(.026); radial.SetTextAlign(12)
            # Intersections of the five rings with a common 135-degree radial
            # line. The separate x/y radii account for the pad aspect ratio.
            for ring,value in enumerate((.1,.2,.3,.4,.5),1):
                fraction=ring/5
                radial.DrawLatex(.45-.262*fraction+.006,.50+.304*fraction,f"{value:.1f}")
            radial.DrawLatex(.12,.85,"P_{t} [GeV]")
            owned.append(radial)
        canvas.cd()
        label=ROOT.TLatex(); label.SetNDC(); label.SetTextFont(62); label.SetTextAlign(22); label.SetTextSize(.034)
        label.DrawLatex(x1+.2325,.94,title); owned.extend((pad,label))
    text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.025)
    text.DrawLatex(.04,.185,f"Bin pass: Data > 0, MC > 0, #deltaMC/MC #leq {100*options.max_mc_rel_err:g}%, #delta(Data/MC)/|Data/MC| #leq {100*options.max_ratio_rel_err:g}%.")
    span_factor=options.min_selected_span/math.pi
    span_label="#pi" if math.isclose(span_factor,1.0) else f"{span_factor:g}#pi"
    text.DrawLatex(.04,.140,f"Curve approval: accepted #Delta#phi #geq {span_label}, harmonic design rank = 3, weighted condition #leq {options.max_design_condition:g}.")
    text.DrawLatex(.04,.085,f"Shared kinematics: {first['Phase']}, {first['Pass']}, {first['Run_type']}, {first['Target']}; x = {float(first['x']):g}; Q^{{2}} = {float(first['Q2']):g} GeV^{{2}}; z = {float(first['z']):g}.")
    owned.append(text); canvas._combined_polar_owned=owned; canvas.Print(str(path))
    draw_phi_comparison_page(canvas,path,rows)
    meta=[f"{first['Combination_group_id']}: {first['Phase']}, {first['Pass']}, {first['Run_type']}, {first['Target']}",
          f"x = {float(first['x']):g}, Q2 = {float(first['Q2']):g} GeV^{{2}}, shared z = {float(first['z']):g}"]
    meta += wrapped_items("settings: ", first["Contributing_setting_ids"].split("|"))
    meta += wrapped_items("nominal #theta_{pq} [deg]: ", first["Contributing_thpq_deg"].split("|"))
    canvas.Clear(); canvas.Divide(3,2); owned=[before,after]
    curve_by_pt={int(curve["Pt_bin"]):curve for curve in curves}
    for ipt in range(1,len(PT_EDGES)):
        pad=canvas.cd(ipt); pad.SetGridx(); pad.SetGridy(); pad.SetMargin(.15,.04,.15,.12)
        selected=[row for row in rows if int(row["Pt_bin"])==ipt]
        values=[float(row["Data_by_MC"]) for row in selected if finite(row["Data_by_MC"])]
        errors=[float(row["Data_by_MC_err"]) for row in selected if finite(row["Data_by_MC_err"])]
        ymax=max([2.1]+[v+e for v,e in zip(values,errors)])
        frame=ROOT.TH1D(f"combined_frame_{ipt}","",len(PHI_EDGES)-1,0,2*math.pi)
        frame.SetMinimum(0); frame.SetMaximum(ymax*1.1); frame.GetXaxis().SetTitle("#phi_{pq} [rad]"); frame.GetYaxis().SetTitle("Data / MC"); frame.Draw()
        accepted=ROOT.TGraphErrors(); rejected=ROOT.TGraphErrors()
        for row in selected:
            if not finite(row["Data_by_MC"]): continue
            graph=accepted if int(row["Survives_all_requirements"]) else rejected; point=graph.GetN()
            graph.SetPoint(point,float(row["phipq_center"]),float(row["Data_by_MC"])); graph.SetPointError(point,0,float(row["Data_by_MC_err"]))
        accepted.SetMarkerStyle(20); accepted.SetMarkerColor(ROOT.kGreen+2)
        rejected.SetMarkerStyle(24); rejected.SetMarkerColor(ROOT.kRed+1)
        if accepted.GetN(): accepted.Draw("P SAME")
        if rejected.GetN(): rejected.Draw("P SAME")
        label=ROOT.TLatex(); label.SetNDC(); label.SetTextSize(.042)
        label.DrawLatex(.17,.92,f"{PT_EDGES[ipt-1]:.1f} #leq P_{{t}} {'#leq' if ipt==len(PT_EDGES)-1 else '<'} {PT_EDGES[ipt]:.1f} GeV")
        curve=curve_by_pt[ipt]; status=ROOT.TLatex(); status.SetNDC(); status.SetTextSize(.035)
        status.DrawLatex(.17,.84,f"{curve['Curve_status']}: {curve['N_surviving']}/{curve['N_phi_bins_total']} survive")
        owned.extend((frame,accepted,rejected,label,status))
    pad=canvas.cd(6); pad.Clear(); text=ROOT.TLatex(); text.SetNDC(); text.SetTextSize(.042); y=.88
    for line in meta:
        text.DrawLatex(.07,y,line); y-=.07
    text.DrawLatex(.07,max(.10,y-.02),"green: survives all requirements; red: excluded")
    owned.append(text); canvas._combined_pad_owned=owned; canvas.Print(str(path)); canvas.Print(str(path)+"]"); canvas.Close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combined-yields", type=Path, default=DEFAULT_COMBINED_YIELDS)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--pdf-dir", type=Path, default=DEFAULT_PDF_DIR)
    parser.add_argument("--group", action="append", default=[])
    parser.add_argument("--max-mc-rel-err", type=float, default=.30)
    parser.add_argument("--max-ratio-rel-err", type=float, default=.30)
    parser.add_argument("--ratio-warning-low", type=float, default=.5)
    parser.add_argument("--ratio-warning-high", type=float, default=1.5)
    parser.add_argument("--min-selected-span", type=float, default=math.pi)
    parser.add_argument("--max-design-condition", type=float, default=100.0)
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    grouped=defaultdict(list)
    for row in dmc.read_csv(args.combined_yields): grouped[row["Combination_group_id"]].append(row)
    if args.group: grouped={key:value for key,value in grouped.items() if key in set(args.group)}
    all_rows=[]; all_curves=[]
    for index,(group_id,source) in enumerate(sorted(grouped.items()),1):
        rows,curves=evaluate_combined_group(source,args); all_rows.extend(rows); all_curves.extend(curves)
        if not args.no_pdf:
            draw_combined_pdf(args.pdf_dir/z_directory(rows[0]["z"])/f"{group_id}.pdf",rows,curves,args)
        print(f"[{index}/{len(grouped)}] OK: {group_id}")
    atomic_csv(args.table_dir / "RP_azimuthal_coverage_bins.csv", COMBINED_BIN_FIELDS, all_rows)
    atomic_csv(args.table_dir / "RP_azimuthal_coverage_curves.csv", COMBINED_CURVE_FIELDS, all_curves)
    atomic_csv(args.table_dir / "RP_azimuthal_coverage_bins_Problematic.csv", COMBINED_BIN_FIELDS,
               [row for row in all_rows if row["Bin_status"] != "PASS"])
    atomic_csv(args.table_dir / "RP_azimuthal_coverage_curves_Problematic.csv", COMBINED_CURVE_FIELDS,
               [row for row in all_curves if row["Curve_status"] != "APPROVED"])
    print("Bins:", Counter(row["Bin_status"] for row in all_rows))
    print("Curves:", Counter(row["Curve_status"] for row in all_curves))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
