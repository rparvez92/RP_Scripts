#!/usr/bin/env python3
"""Fit coverage-approved Pass-5 reduced SIDIS cross sections versus phi."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import ROOT

from RP_sidis_binning import AZIMUTHAL_PHI_EDGES as PHI_EDGES, AZIMUTHAL_PT_EDGES as PT_EDGES

ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray + 1); ROOT.gStyle.SetGridStyle(3)
ROOT.TGaxis.SetMaxDigits(3)

POLICY_VERSION = "RP_FIT_REDUCED_XSEC_PASS5_V1"
FAMILY_FIELDS = (
    "Epsilon_family_id", "Epsilon_family_policy_version", "Family_status",
    "Family_reason", "Run_type", "Target", "x", "Q2", "z", "N_members",
    "Pass_pattern", "Has_pass3", "Has_pass4", "Has_pass5",
    "N_distinct_epsilon", "Epsilon_min", "Epsilon_max", "Epsilon_span",
    "Member_index", "Member_phase", "Member_pass",
)
PROVENANCE_FIELDS = (
    "Combination_group_id", "Combination_policy_version", "Ebeam_GeV",
    "N_contributors", "Contributing_setting_ids", "Contributing_thpq_deg",
)
BIN_FIELDS = (
    "Pt_bin", "Pt_low", "Pt_high", "Pt_high_closed", "phipq_bin",
    "phipq_low", "phipq_high", "phipq_center", "Binning_version",
)
POINT_EXTRA = (
    "Coverage_policy_version", "Coverage_bin_status", "Coverage_bin_reason",
    "Coverage_survives_all_requirements", "Coverage_curve_status",
    "Coverage_curve_reason", "Coverage_selected_mask", "Used_in_fit",
    "Fit_point_status", "Fit_point_reason",
)
FIT_FIELDS = [
    *FAMILY_FIELDS, *PROVENANCE_FIELDS, "epsilon", "Pt_bin", "Pt_low",
    "Pt_high", "Pt_high_closed", "Coverage_policy_version",
    "Coverage_curve_status", "Coverage_curve_reason", "Coverage_selected_mask",
    "N_phi_bins_total", "N_points_fitted", "Phi_coverage_rad",
    "Weighted_design_rank", "Weighted_design_condition", "C_0_pass5",
    "C_0_pass5_err", "C_1_pass5", "C_1_pass5_err", "C_2_pass5",
    "C_2_pass5_err", "Cov_00", "Cov_01", "Cov_02", "Cov_11", "Cov_12",
    "Cov_22", "Correlation_01", "Correlation_02", "Correlation_12",
    "Chi2", "NDF", "Chi2_NDF", "Fit_probability", "Fit_status",
    "Fit_reason", "Fit_line_style", "Fit_policy_version",
]


def finite(value: object) -> bool:
    try: return math.isfinite(float(value))
    except (TypeError, ValueError): return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def atomic_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def bin_key(row: Mapping[str, object]) -> tuple[str, int, int]:
    return str(row["Combination_group_id"]), int(row["Pt_bin"]), int(row["phipq_bin"])


def curve_key(row: Mapping[str, object]) -> tuple[str, int]:
    return str(row["Combination_group_id"]), int(row["Pt_bin"])


def full_2pi_curve(curve: Mapping[str, object] | None) -> bool:
    total = len(PHI_EDGES) - 1
    return bool(curve and curve.get("Curve_status") == "APPROVED"
                and int(curve.get("N_selected", 0)) == total
                and str(curve.get("Surviving_mask", "")) == "1" * total
                and math.isclose(number(curve.get("Min_selected_span_rad")),
                                 2 * math.pi, rel_tol=0, abs_tol=1e-12))


def annotate_points(points: Sequence[Mapping[str, object]],
                    coverage_bins: Mapping, coverage_curves: Mapping
                    ) -> list[dict[str, object]]:
    output = []
    for source in points:
        row = dict(source); bin_decision = coverage_bins.get(bin_key(source))
        curve = coverage_curves.get(curve_key(source)); reasons = []
        if bin_decision is None: reasons.append("COVERAGE_BIN_MISSING")
        if curve is None: reasons.append("COVERAGE_CURVE_MISSING")
        survives = bool(bin_decision and bin_decision.get("Survives_all_requirements") == "1")
        approved = full_2pi_curve(curve)
        if bin_decision is not None and not survives: reasons.append("POINT_FAILS_COVERAGE")
        if curve is not None and not approved: reasons.append("COVERAGE_CURVE_NOT_FULL_2PI_APPROVED")
        reduction_ok = source.get("Reduction_status") in {"OK", "WARNING"}
        if not reduction_ok: reasons.append("REDUCTION_POINT_NOT_USABLE")
        value = number(source.get("Reduced_data_xsec")); error = number(source.get("Reduced_data_xsec_err"))
        if not finite(value): reasons.append("NONFINITE_REDUCED_VALUE")
        if not finite(error) or error <= 0: reasons.append("NONPOSITIVE_REDUCED_UNCERTAINTY")
        usable = reduction_ok and finite(value) and finite(error) and error > 0 and survives and approved
        row.update({
            "Coverage_policy_version": bin_decision.get("Policy_version", "") if bin_decision else "",
            "Coverage_bin_status": bin_decision.get("Bin_status", "MISSING") if bin_decision else "MISSING",
            "Coverage_bin_reason": bin_decision.get("Bin_reason", "COVERAGE_BIN_MISSING") if bin_decision else "COVERAGE_BIN_MISSING",
            "Coverage_survives_all_requirements": int(survives),
            "Coverage_curve_status": curve.get("Curve_status", "MISSING") if curve else "MISSING",
            "Coverage_curve_reason": curve.get("Curve_reason", "COVERAGE_CURVE_MISSING") if curve else "COVERAGE_CURVE_MISSING",
            "Coverage_selected_mask": curve.get("Selected_mask", "") if curve else "",
            "Used_in_fit": int(usable),
            "Fit_point_status": "USABLE" if usable else "EXCLUDED",
            "Fit_point_reason": ";".join(sorted(set(reasons))),
        })
        output.append(row)
    return output


def fit_curve(points: Sequence[Mapping[str, object]]) -> dict[str, object]:
    first = points[0]; curve = [row for row in points if int(row["Used_in_fit"]) == 1]
    result = {field: "" for field in FIT_FIELDS}
    for field in (*FAMILY_FIELDS, *PROVENANCE_FIELDS, "epsilon", "Pt_bin",
                  "Pt_low", "Pt_high", "Pt_high_closed",
                  "Coverage_policy_version", "Coverage_curve_status",
                  "Coverage_curve_reason", "Coverage_selected_mask"):
        result[field] = first.get(field, "")
    n = len(curve); result.update(
        N_phi_bins_total=len(PHI_EDGES)-1, N_points_fitted=n,
        Phi_coverage_rad=n * 2 * math.pi / (len(PHI_EDGES)-1),
        Fit_policy_version=POLICY_VERSION,
    )
    if n != len(PHI_EDGES) - 1:
        result.update(Fit_status="NOT_FITTED", Fit_reason=(
            first.get("Coverage_curve_reason") or "FULL_2PI_COVERAGE_REQUIRED"
        ), Fit_line_style="none")
        return result
    phi = np.asarray([number(row["phipq_center"]) for row in curve])
    values = np.asarray([number(row["Reduced_data_xsec"]) for row in curve])
    errors = np.asarray([number(row["Reduced_data_xsec_err"]) for row in curve])
    design = np.column_stack((np.ones(n), np.cos(phi), np.cos(2 * phi)))
    weighted = design / errors[:, None]
    rank = int(np.linalg.matrix_rank(weighted)); condition = float(np.linalg.cond(weighted))
    result.update(Weighted_design_rank=rank, Weighted_design_condition=condition)
    if rank < 3:
        result.update(Fit_status="ERROR", Fit_reason="WEIGHTED_DESIGN_RANK_BELOW_3",
                      Fit_line_style="none"); return result
    normal = weighted.T @ weighted; covariance = np.linalg.inv(normal)
    parameters = covariance @ (design.T @ (values / errors**2))
    prediction = design @ parameters; chi2 = float(np.sum(((values - prediction) / errors)**2))
    ndf = n - 3; probability = float(ROOT.TMath.Prob(chi2, ndf))
    for index, name in enumerate(("C_0_pass5", "C_1_pass5", "C_2_pass5")):
        result[name] = parameters[index]; result[name + "_err"] = math.sqrt(max(0, covariance[index,index]))
    for name, i, j in (("Cov_00",0,0),("Cov_01",0,1),("Cov_02",0,2),
                       ("Cov_11",1,1),("Cov_12",1,2),("Cov_22",2,2)):
        result[name] = covariance[i,j]
    for name, i, j in (("Correlation_01",0,1),("Correlation_02",0,2),("Correlation_12",1,2)):
        result[name] = covariance[i,j] / math.sqrt(covariance[i,i] * covariance[j,j])
    reasons = []
    if probability < 0.01: reasons.append("POOR_FIT_PROBABILITY")
    if parameters[0] <= 0: reasons.append("NONPOSITIVE_C0")
    result.update(Chi2=chi2, NDF=ndf, Chi2_NDF=chi2/ndf,
                  Fit_probability=probability, Fit_status="WARNING" if reasons else "OK",
                  Fit_reason=";".join(reasons), Fit_line_style="solid")
    return result


def graph(rows, usable, marker, color):
    result = ROOT.TGraphErrors(); index = 0
    for row in rows:
        if bool(int(row["Used_in_fit"])) != usable or not finite(row["Reduced_data_xsec"]): continue
        result.SetPoint(index, number(row["phipq_center"]), number(row["Reduced_data_xsec"]))
        result.SetPointError(index, 0, number(row["Reduced_data_xsec_err"], 0)); index += 1
    result.SetMarkerStyle(marker); result.SetMarkerColor(color); result.SetLineColor(color)
    return result


def draw_pdf(path: Path, groups: Sequence[tuple[list[dict], list[dict]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas = ROOT.TCanvas("c_reduced_pass5", "reduced pass5", 1800, 1200)
    canvas.SetCanvasSize(1800,1200); canvas.Print(str(path)+"[")
    for points, fits in groups:
        first = points[0]; by_pt = defaultdict(list)
        for row in points: by_pt[int(row["Pt_bin"])].append(row)
        fit_by_pt = {int(row["Pt_bin"]): row for row in fits}
        canvas.Clear(); canvas.Divide(3,2,.002,.002); owned=[]
        for pt_bin in range(1,6):
            pad=canvas.cd(pt_bin); pad.SetGridx(); pad.SetGridy(); pad.SetMargin(.15,.04,.15,.11)
            rows=by_pt[pt_bin]; usable=[r for r in rows if int(r["Used_in_fit"]) and finite(r["Reduced_data_xsec"])]
            scale=usable or [r for r in rows if finite(r["Reduced_data_xsec"])]
            lows=[number(r["Reduced_data_xsec"])-number(r["Reduced_data_xsec_err"],0) for r in scale]
            highs=[number(r["Reduced_data_xsec"])+number(r["Reduced_data_xsec_err"],0) for r in scale]
            low=min(lows+[0.0]); high=max(highs+[1e-12]); span=max(high-low,1e-12)
            frame=ROOT.TH1D(f"rf_{abs(hash((first['Combination_group_id'],pt_bin)))}","",60,0,2*math.pi)
            frame.SetMinimum(min(0,low-.1*span)); frame.SetMaximum(high+.18*span)
            frame.SetLineColorAlpha(ROOT.kWhite,0); frame.GetXaxis().SetTitle("#phi_{pq} [rad]")
            frame.GetYaxis().SetTitle("#sigma_{red} [GeV^{-2}]"); frame.Draw()
            label=ROOT.TLatex(); label.SetNDC(); label.SetTextSize(.043)
            relation="#leq" if pt_bin==5 else "<"
            label.DrawLatex(.20,.92,f"{PT_EDGES[pt_bin-1]:.1f} #leq P_{{T}} {relation} {PT_EDGES[pt_bin]:.1f} GeV")
            zero=ROOT.TLine(0,0,2*math.pi,0); zero.SetLineColor(ROOT.kViolet+1); zero.SetLineStyle(2); zero.Draw()
            good=graph(rows,True,20,ROOT.kBlue+1); bad=graph(rows,False,24,ROOT.kGray+2)
            if bad.GetN(): bad.Draw("P SAME")
            if good.GetN(): good.Draw("P SAME")
            fit=fit_by_pt[pt_bin]
            if fit["Fit_line_style"]=="solid":
                func=ROOT.TF1(f"rfit_{abs(hash((first['Combination_group_id'],pt_bin)))}","[0]+[1]*cos(x)+[2]*cos(2*x)",0,2*math.pi)
                for i,name in enumerate(("C_0_pass5","C_1_pass5","C_2_pass5")): func.SetParameter(i,number(fit[name]))
                func.SetLineColor(ROOT.kRed+1); func.SetLineWidth(2); func.Draw("SAME"); owned.append(func)
            elif fit["Fit_status"]=="NOT_FITTED":
                notice=ROOT.TLatex(); notice.SetNDC(); notice.SetTextAlign(22); notice.SetTextSize(.042); notice.SetTextColorAlpha(ROOT.kGray+2,.65)
                notice.DrawLatex(.55,.82,"Not fitted - coverage rejected"); owned.append(notice)
            owned.extend((frame,label,zero,good,bad))
        pad=canvas.cd(6); pad.Clear(); text=ROOT.TLatex(); text.SetNDC(); text.SetTextSize(.038)
        lines=["Pass-5 anchor: reduced SIDIS xsec",f"{first['Run_type']} / {first['Target']}",
               f"x={number(first['x']):.3g}, Q^{{2}}={number(first['Q2']):.3g} GeV^{{2}}, z={number(first['z']):.3g}",
               f"E={number(first['Ebeam_GeV']):.5g} GeV, #epsilon={number(first['epsilon']):.5f}",
               f"Family: {first['Epsilon_family_id']}","Fit: C_{0}^{(5)} + C_{1}^{(5)} cos#phi + C_{2}^{(5)} cos2#phi"]
        for i,line in enumerate(lines): text.DrawLatex(.05,.93-.075*i,line)
        valid=ROOT.TMarker(.09,.38,20); valid.SetNDC(True); valid.SetMarkerColor(ROOT.kBlue+1); valid.Draw()
        invalid=ROOT.TMarker(.09,.30,24); invalid.SetNDC(True); invalid.SetMarkerColor(ROOT.kGray+2); invalid.Draw()
        line=ROOT.TLine(.07,.22,.12,.22); line.SetNDC(True); line.SetLineColor(ROOT.kRed+1); line.SetLineWidth(2); line.Draw()
        text.DrawLatex(.14,.365,"Used in fit"); text.DrawLatex(.14,.285,"Excluded diagnostic point"); text.DrawLatex(.14,.205,"Pass-5 harmonic fit")
        text.SetTextSize(.030); text.DrawLatex(.05,.10,str(first["Coverage_policy_version"])); text.DrawLatex(.05,.05,POLICY_VERSION)
        owned.extend((text,valid,invalid,line)); canvas._reduced_fit_owned=owned; canvas.Print(str(path))
    canvas.Print(str(path)+"]"); canvas.Close()


def main() -> int:
    project=Path(__file__).resolve().parents[1]; parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced-points",type=Path,default=project/"results/Pass1/Tables/RP_extract_sidis_reduced_xsec/RP_sidis_reduced_xsec_points.csv")
    parser.add_argument("--coverage-bins",type=Path,default=project/"results/Pass1/Tables/RP_azimuthal_coverage_study/RP_azimuthal_coverage_bins.csv")
    parser.add_argument("--coverage-curves",type=Path,default=project/"results/Pass1/Tables/RP_azimuthal_coverage_study/RP_azimuthal_coverage_curves.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_fit_sidis_reduced_xsec_pass5")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_fit_sidis_reduced_xsec_pass5")
    parser.add_argument("--no-pdf",action="store_true"); args=parser.parse_args()
    raw=[row for row in read_csv(args.reduced_points) if row["Member_pass"]=="pass5" and row["Family_status"].startswith("LT_READY")]
    bins={bin_key(row):row for row in read_csv(args.coverage_bins)}; curves={curve_key(row):row for row in read_csv(args.coverage_curves)}
    points=annotate_points(raw,bins,curves); grouped=defaultdict(list)
    for row in points: grouped[curve_key(row)].append(row)
    fits=[fit_curve(sorted(rows,key=lambda r:int(r["phipq_bin"]))) for _,rows in sorted(grouped.items())]
    point_fields=list(points[0]); atomic_csv(args.table_dir/"RP_sidis_reduced_xsec_pass5_points.csv",point_fields,points)
    atomic_csv(args.table_dir/"RP_sidis_reduced_xsec_pass5_points_Problematic.csv",point_fields,[r for r in points if r["Fit_point_status"]!="USABLE"])
    atomic_csv(args.table_dir/"RP_sidis_reduced_xsec_pass5_fit_parameters.csv",FIT_FIELDS,fits)
    atomic_csv(args.table_dir/"RP_sidis_reduced_xsec_pass5_fit_parameters_Problematic.csv",FIT_FIELDS,[r for r in fits if r["Fit_status"]!="OK"])
    if not args.no_pdf:
        by_group=defaultdict(list); fit_group=defaultdict(list)
        for row in points: by_group[row["Combination_group_id"]].append(row)
        for row in fits: fit_group[row["Combination_group_id"]].append(row)
        for z in sorted({row["z"] for row in points},key=float):
            groups=[(by_group[key],fit_group[key]) for key in sorted(by_group) if by_group[key][0]["z"]==z]
            zdir="z_"+f"{float(z):.2f}".replace(".","p")
            draw_pdf(args.pdf_dir/zdir/"RP_sidis_reduced_xsec_pass5_vs_phipq.pdf",groups)
    print("Points:",len(points),Counter(r["Fit_point_status"] for r in points))
    print("Fits:",len(fits),Counter(r["Fit_status"] for r in fits)); return 0


if __name__=="__main__": raise SystemExit(main())
