#!/usr/bin/env python3
"""Extract and fit Pass1 SIDIS sigma/multiplicity versus phipq."""

from __future__ import annotations

import argparse, csv, math, os, tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import ROOT

import RP_data_mc_compare as dmc
from RP_build_sidis_model import calculator_identity, run_calculator_rows
from RP_sidis_binning import (
    AZIMUTHAL_BINNING_VERSION,
    AZIMUTHAL_PHI_EDGES as PHI_EDGES,
    AZIMUTHAL_PT_EDGES as PT_EDGES,
    bin_index,
    canonical_phi,
)

ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray + 1); ROOT.gStyle.SetGridStyle(3)
ROOT.TGaxis.SetMaxDigits(3)

IDENTITY = ("Combination_group_id",)
GROUP_FIELDS = ("Phase", "Pass", "Run_type", "Target", "x", "Q2")
PROVENANCE_FIELDS = (
    "Combination_group_id", "Combination_policy_version", "Ebeam_GeV",
    "N_contributors", "Contributing_setting_ids", "Contributing_thpq_deg",
)
REACTIONS = dmc.REACTIONS
MODEL_POLICY = "CALC_SEMI_XSEC_AT_15_PHI_BIN_CENTERS"
Z_COLORS = {0.36: ROOT.kRed + 1, 0.50: ROOT.kGreen + 2,
            0.52: ROOT.kGreen + 2, 0.67: ROOT.kBlue + 1, 0.90: ROOT.kBlack}
REJECTION_CODES = {
    "INSUFFICIENT_ANGULAR_SPAN": 1,
    "DESIGN_RANK_BELOW_3": 2,
    "WEIGHTED_DESIGN_CONDITION_ABOVE_LIMIT": 3,
}

POINT_FIELDS = [
    *PROVENANCE_FIELDS, *GROUP_FIELDS, "z", "Pt_bin", "Pt_low", "Pt_high",
    "Pt_high_closed", "phipq_bin", "phipq_low", "phipq_high",
    "phipq_center", "Y_Data", "Y_Data_err",
    *(item for reaction in REACTIONS for item in (f"Y_MC_{reaction}", f"Y_MC_{reaction}_err")),
    "Y_MC", "Y_MC_err", "C_Y", "C_Y_err", "Model_value", "Model_units",
    "Model_Pt_GeV", "Model_phi_rad", "Model_phi_policy", "Extracted_value", "Extracted_err",
    "Observable_units", "Coverage_binning_version", "Coverage_policy_version",
    "Coverage_bin_status", "Coverage_bin_reason", "Coverage_independent_pass",
    "Coverage_selected_contiguous", "Coverage_curve_status", "Coverage_curve_reason",
    "Used_in_fit", "Point_status", "Point_reason",
]
FIT_COMMON = [
    *PROVENANCE_FIELDS, *GROUP_FIELDS, "z", "Pt_bin", "Pt_low", "Pt_high",
    "Pt_high_closed", "N_phi_bins_total", "N_points_extracted", "N_points_usable",
    "N_points_fitted", "Phi_coverage_rad", "Phi_coverage_fraction", "Phi_bin_mask",
    "Model_value", "Model_units", "Model_Pt_GeV", "Model_phi_policy",
    "Coverage_binning_version", "Coverage_policy_version",
    "Coverage_curve_status", "Coverage_curve_reason", "Coverage_selected_mask",
]
FIT_DIAGNOSTICS = [
    "Cov_00", "Cov_01", "Cov_02", "Cov_11", "Cov_12", "Cov_22",
    "Correlation_01", "Correlation_02", "Correlation_12", "Chi2", "NDF",
    "Chi2_NDF", "Fit_probability", "Fit_status", "Fit_reason", "Fit_line_style",
]


def finite(value) -> bool:
    try: return math.isfinite(float(value))
    except (TypeError, ValueError): return False


def number(value, default=math.nan): return float(value) if finite(value) else default


def coverage_bin_key(row):
    return (str(row["Combination_group_id"]),int(row["Pt_bin"]),int(row["phipq_bin"]))


def coverage_curve_key(row):
    return (str(row["Combination_group_id"]),int(row["Pt_bin"]))


def load_coverage(bin_path: Path, curve_path: Path):
    bin_rows=dmc.read_csv(bin_path); curve_rows=dmc.read_csv(curve_path)
    versions={row.get("Binning_version") for row in (*bin_rows,*curve_rows)}
    if versions!={AZIMUTHAL_BINNING_VERSION}:
        raise ValueError(f"coverage binning mismatch: expected {AZIMUTHAL_BINNING_VERSION}, found {sorted(versions)}")
    bins={coverage_bin_key(row):row for row in bin_rows}
    curves={coverage_curve_key(row):row for row in curve_rows}
    if len(bins)!=len(bin_rows) or len(curves)!=len(curve_rows):
        raise ValueError("duplicate coverage decision key")
    return bins,curves


def atomic_csv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def make_hist(name: str) -> ROOT.TH2D:
    hist=ROOT.TH2D(name,"",len(PHI_EDGES)-1,np.asarray(PHI_EDGES,"d"),len(PT_EDGES)-1,np.asarray(PT_EDGES,"d"))
    hist.Sumw2(); hist.SetDirectory(0); return hist


def rdf_hist(path, tree, phi_expr, pt_expr, selection, weight, name):
    frame=(ROOT.RDataFrame(tree,str(path)).Filter(selection).Define("_az_phi",phi_expr)
           .Define("_az_pt",pt_expr).Define("_az_w",weight)
           .Filter("std::isfinite((double)_az_phi) && std::isfinite((double)_az_pt) && _az_pt>=0 && _az_pt<=0.5"))
    action=frame.Histo2D((name,"",len(PHI_EDGES)-1,np.asarray(PHI_EDGES,"d"),len(PT_EDGES)-1,np.asarray(PT_EDGES,"d")),
                         "_az_phi","_az_pt","_az_w")
    hist=action.GetValue().Clone(name+"_owned"); hist.SetDirectory(0); hist.Sumw2(); return hist


def wrap_phi(expression: str) -> str:
    return f"fmod(fmod((double)({expression}),2*TMath::Pi())+2*TMath::Pi(),2*TMath::Pi())"


def data_class(path: Path, root_dir: Path, prefix: str):
    rows=[row for row in dmc.read_csv(path) if row.get("Norm_status")=="OK"
          and finite(row.get("BCM2_Q")) and float(row["BCM2_Q"])>0]
    total_charge=sum(float(row["BCM2_Q"]) for row in rows); output=make_hist(prefix)
    reasons=[]
    for row in rows:
        run=int(float(row["Run"])); root=dmc.data_root_path(root_dir,dmc.parse_norm_stem(path)["Phase"],run)
        if dmc.is_floored_positron_run(path,row): continue
        if not root.is_file(): reasons.append(f"DATA_ROOT_MISSING_RUN_{run}"); continue
        timing=dmc.timing_weight_expression(row)
        selection=f"({dmc.PID_CUT}) && ({dmc.DATA_DELTA}) && ({timing})!=0"
        hist=rdf_hist(root,"T",wrap_phi("P_kin_secondary_ph_xq"),"abs(pt)",selection,
                      f"({float(row['Norm_factor']):.17g})*({timing})",f"{prefix}_{run}")
        output.Add(hist,float(row["BCM2_Q"])/total_charge)
    return output,reasons


def simc_class(row, prefix):
    path=Path(row.get("Recon_root_path","")); output=make_hist(prefix)
    if not path.is_file(): return output,["SIMC_ROOT_MISSING"]
    expressions=dmc.simc_expressions(row["Reaction"])
    try:
        return rdf_hist(path,"h10",wrap_phi(expressions["phipq"]),expressions["pt"],
                        dmc.SIMC_DELTA,"static_cast<double>(fWeight)",prefix),[]
    except Exception as error: return output,[f"SIMC_HIST_FAILED={type(error).__name__}"]


def combine_data(classes, scale):
    result=classes["target_e"].Clone("az_data"); result.Add(classes["target_p"],-1)
    dummy=classes["dummy_e"].Clone("az_dummy"); dummy.Add(classes["dummy_p"],-1)
    result.Add(dummy,-scale); result.SetDirectory(0); return result


def fit_curve(points, observable):
    usable=[row for row in points if row["Used_in_fit"]==1]
    result={field:"" for field in FIT_DIAGNOSTICS}; n=len(usable)
    extracted=sum(finite(row.get("Extracted_value")) for row in points)
    result.update(N_points_extracted=extracted,N_points_usable=n,N_points_fitted=n,
                  N_phi_bins_total=len(PHI_EDGES)-1,Phi_bin_mask="".join("1" if any(p["phipq_bin"]==i and p["Used_in_fit"]==1 for p in points) else "0" for i in range(1,len(PHI_EDGES))),
                  Phi_coverage_rad=n*2*math.pi/(len(PHI_EDGES)-1),Phi_coverage_fraction=n/(len(PHI_EDGES)-1))
    if n<3:
        result.update(Fit_status="NOT_FITTED",Fit_reason="NO_USABLE_PHI_POINTS" if n==0 else "INSUFFICIENT_PHI_POINTS",Fit_line_style="none")
        return result
    phi=np.array([number(row["phipq_center"]) for row in usable],dtype=float)
    y=np.array([number(row["Extracted_value"]) for row in usable],dtype=float)
    error=np.array([number(row["Extracted_err"]) for row in usable],dtype=float)
    design=np.column_stack((np.ones(n),np.cos(phi),np.cos(2*phi))); weighted=design/error[:,None]
    normal=weighted.T@weighted
    if np.linalg.matrix_rank(normal)<3:
        result.update(Fit_status="ERROR",Fit_reason="SINGULAR_FIT",Fit_line_style="none"); return result
    covariance=np.linalg.inv(normal); parameters=covariance@(design.T@(y/error**2))
    if not all(finite(value) for value in parameters):
        result.update(Fit_status="ERROR",Fit_reason="NONFINITE_FIT_COEFFICIENT",Fit_line_style="none"); return result
    prediction=design@parameters; chi2=float(np.sum(((y-prediction)/error)**2)); ndf=n-3
    probability=float(ROOT.TMath.Prob(chi2,ndf)) if ndf>0 else math.nan
    names=("C_0","C_1","C_2") if observable=="sigma" else ("M_0","M_1","M_2")
    for index,name in enumerate(names): result[name]=parameters[index]; result[name+"_err"]=math.sqrt(max(0,covariance[index,index]))
    for (name,i,j) in (("Cov_00",0,0),("Cov_01",0,1),("Cov_02",0,2),("Cov_11",1,1),("Cov_12",1,2),("Cov_22",2,2)):
        result[name]=covariance[i,j]
    for name,i,j in (("Correlation_01",0,1),("Correlation_02",0,2),("Correlation_12",1,2)):
        result[name]=covariance[i,j]/math.sqrt(covariance[i,i]*covariance[j,j]) if covariance[i,i]>0 and covariance[j,j]>0 else math.nan
    reasons=[]
    if ndf==0: reasons.append("EXACTLY_DETERMINED_NDF0")
    if ndf>0 and probability<0.01: reasons.append("POOR_FIT_PROBABILITY")
    if parameters[0]<=0: reasons.append("NONPOSITIVE_NORMALIZATION")
    result.update(Chi2=chi2,NDF=ndf,Chi2_NDF=chi2/ndf if ndf else "",
                  Fit_probability=probability if ndf else "",Fit_status="WARNING" if reasons else "OK",
                  Fit_reason=";".join(reasons),Fit_line_style="dashed" if ndf==0 else "solid")
    return result


def model_values(calculator, model, pt):
    a,z,charge=calculator_identity(model["Target"],model["Run_type"])
    rows=run_calculator_rows(
        calculator,
        (number(model["Ebeam_GeV"]),number(model["Q2_GeV2"]),number(model["x"]),
         number(model["z"]),pt,a,z,charge,float(len(PHI_EDGES)-1)),
    )
    if len(rows) != len(PHI_EDGES)-1:
        raise ValueError(
            f"calc_semi_xsec returned {len(rows)} phi rows; expected {len(PHI_EDGES)-1}"
        )
    values=[]
    for index,(phi,multiplicity,sigma) in enumerate(rows):
        expected=(PHI_EDGES[index]+PHI_EDGES[index+1])/2
        if not math.isclose(phi,expected,rel_tol=0.0,abs_tol=1e-10):
            raise ValueError(
                f"calc_semi_xsec phi center mismatch at bin {index+1}: {phi} != {expected}"
            )
        if multiplicity <= 0 or sigma <= 0:
            raise ValueError(f"calc_semi_xsec returned nonpositive model at phi bin {index+1}")
        values.append((phi,sigma,multiplicity))
    return values


def extract_setting(key,files,catalog,model,calculator,root_dir,coverage_bins,coverage_curves):
    phase,pas,run_type,target,setting=key; reasons=[]
    categories={"target_e":(target,"Elec"),"target_p":(target,"Pos"),"dummy_e":("DUMMY","Elec"),"dummy_p":("DUMMY","Pos")}
    classes={}
    for name,category in categories.items():
        if category in files: classes[name],found=data_class(files[category],root_dir,f"{setting}_{name}"); reasons+=found
        else: classes[name]=make_hist(f"{setting}_{name}_zero"); reasons.append(f"ZERO_CORRECTION_CLASS={name}")
    data=combine_data(classes,dmc.S_DUMMY[(phase,target)])
    components={}; physical=dmc.applicable_reactions(run_type)
    for reaction in REACTIONS:
        if reaction not in physical: components[reaction]=make_hist(f"{setting}_{reaction}_zero"); continue
        source=catalog.get((*key,reaction))
        if source is None: components[reaction]=make_hist(f"{setting}_{reaction}_missing"); reasons.append(f"MISSING_MC={reaction}")
        else: components[reaction],found=simc_class(source,f"{setting}_{reaction}"); reasons+=found
    points={"sigma":[],"multiplicity":[]}
    for ipt,(pt_low,pt_high) in enumerate(zip(PT_EDGES[:-1],PT_EDGES[1:]),1):
        phi_models=model_values(calculator,model,(pt_low+pt_high)/2)
        curve_coverage=coverage_curves.get(tuple(key)+(ipt,))
        for iphi,(phi_low,phi_high) in enumerate(zip(PHI_EDGES[:-1],PHI_EDGES[1:]),1):
            model_phi,sigma_model,mult_model=phi_models[iphi-1]
            bin_coverage=coverage_bins.get(tuple(key)+(ipt,iphi))
            dy,de=data.GetBinContent(iphi,ipt),data.GetBinError(iphi,ipt); values=[]; errors=[]
            for reaction in REACTIONS: values.append(components[reaction].GetBinContent(iphi,ipt)); errors.append(components[reaction].GetBinError(iphi,ipt))
            mc=sum(values); me=math.sqrt(sum(item*item for item in errors)); cy=dy/mc if mc>0 else math.nan
            ce=math.hypot(de/mc,dy*me/(mc*mc)) if mc>0 else math.nan
            for observable,mod,units in (("sigma",sigma_model,"ub/GeV^2/sr^2"),("multiplicity",mult_model,"GeV^-2")):
                value=cy*mod if finite(cy) else math.nan; error=abs(mod)*ce if finite(ce) else math.nan
                point_reasons=list(reasons)
                if not finite(mc) or mc<=0: point_reasons.append("NONPOSITIVE_MC")
                if not finite(value): point_reasons.append("NONFINITE_EXTRACTED_VALUE")
                if not finite(error) or error<=0: point_reasons.append("NONPOSITIVE_UNCERTAINTY")
                if finite(value) and value<0: point_reasons.append("NEGATIVE_EXTRACTED_VALUE")
                if bin_coverage is None: point_reasons.append("COVERAGE_BIN_MISSING")
                if curve_coverage is None: point_reasons.append("COVERAGE_CURVE_MISSING")
                coverage_selected=bin_coverage is not None and bin_coverage.get("Selected_contiguous")=="1"
                coverage_approved=curve_coverage is not None and curve_coverage.get("Curve_status")=="APPROVED"
                if bin_coverage is not None and not coverage_selected: point_reasons.append("OUTSIDE_APPROVED_CONTIGUOUS_RANGE")
                if curve_coverage is not None and not coverage_approved: point_reasons.append("COVERAGE_CURVE_REJECTED")
                usable=finite(value) and finite(error) and error>0 and mc>0 and mod>0 and coverage_selected and coverage_approved
                row={field:"" for field in POINT_FIELDS}; row.update({
                    "Phase":phase,"Pass":pas,"Run_type":run_type,"Target":target,"x":model["x"],"Q2":model["Q2_GeV2"],"thpq":model["theta_pq_deg"],"z":model["z"],"setting_id":setting,
                    "Pt_bin":ipt,"Pt_low":pt_low,"Pt_high":pt_high,"Pt_high_closed":int(ipt==5),"phipq_bin":iphi,"phipq_low":phi_low,"phipq_high":phi_high,"phipq_center":(phi_low+phi_high)/2,
                    "Y_Data":dy,"Y_Data_err":de,"Y_MC":mc,"Y_MC_err":me,"C_Y":cy,"C_Y_err":ce,"Model_value":mod,"Model_units":units,"Model_Pt_GeV":(pt_low+pt_high)/2,"Model_phi_rad":model_phi,"Model_phi_policy":MODEL_POLICY,"Extracted_value":value,"Extracted_err":error,"Observable_units":units,
                    "Coverage_binning_version":bin_coverage.get("Binning_version","") if bin_coverage else "",
                    "Coverage_policy_version":bin_coverage.get("Policy_version","") if bin_coverage else "",
                    "Coverage_bin_status":bin_coverage.get("Bin_status","") if bin_coverage else "MISSING",
                    "Coverage_bin_reason":bin_coverage.get("Bin_reason","") if bin_coverage else "COVERAGE_BIN_MISSING",
                    "Coverage_independent_pass":bin_coverage.get("Independent_pass","") if bin_coverage else "",
                    "Coverage_selected_contiguous":int(coverage_selected),
                    "Coverage_curve_status":curve_coverage.get("Curve_status","") if curve_coverage else "MISSING",
                    "Coverage_curve_reason":curve_coverage.get("Curve_reason","") if curve_coverage else "COVERAGE_CURVE_MISSING",
                    "Used_in_fit":int(usable),"Point_status":"WARNING" if point_reasons else "OK","Point_reason":";".join(sorted(set(point_reasons)))})
                for reaction,v,e in zip(REACTIONS,values,errors): row[f"Y_MC_{reaction}"]=v; row[f"Y_MC_{reaction}_err"]=e
                points[observable].append(row)
    return points


def extract_group(group_rows,calculator,coverage_bins,coverage_curves):
    """Extract observables from already-combined corrected-data and MC yields."""
    ordered=sorted(group_rows,key=lambda row:(int(row["Pt_bin"]),int(row["phipq_bin"])))
    first=ordered[0]; group_id=first["Combination_group_id"]
    expected=(len(PT_EDGES)-1)*(len(PHI_EDGES)-1)
    if len(ordered)!=expected:
        raise ValueError(f"{group_id} has {len(ordered)} combined bins; expected {expected}")
    model={"Target":first["Target"],"Run_type":first["Run_type"],
           "Ebeam_GeV":first["Ebeam_GeV"],"Q2_GeV2":first["Q2"],
           "x":first["x"],"z":first["z"]}
    models={ipt:model_values(calculator,model,(low+high)/2)
            for ipt,(low,high) in enumerate(zip(PT_EDGES[:-1],PT_EDGES[1:]),1)}
    points={"sigma":[],"multiplicity":[]}
    for source in ordered:
        ipt,iphi=int(source["Pt_bin"]),int(source["phipq_bin"])
        model_phi,sigma_model,mult_model=models[ipt][iphi-1]
        bin_coverage=coverage_bins.get((group_id,ipt,iphi))
        curve_coverage=coverage_curves.get((group_id,ipt))
        dy,de=number(source["Y_Data"]),number(source["Y_Data_err"])
        mc,me=number(source["Y_MC"]),number(source["Y_MC_err"])
        cy,ce=number(source.get("Data_by_MC")),number(source.get("Data_by_MC_err"))
        for observable,mod,units in (("sigma",sigma_model,"ub/GeV^2/sr^2"),
                                     ("multiplicity",mult_model,"GeV^-2")):
            value=cy*mod if finite(cy) else math.nan
            error=abs(mod)*ce if finite(ce) else math.nan
            reasons=list(filter(None,str(source.get("Combination_reason","")).split(";")))
            combination_ok=source.get("Combination_status")!="ERROR"
            if not combination_ok: reasons.append("COMBINATION_ERROR")
            if not finite(mc) or mc<=0: reasons.append("NONPOSITIVE_MC")
            if not finite(value): reasons.append("NONFINITE_EXTRACTED_VALUE")
            if not finite(error) or error<=0: reasons.append("NONPOSITIVE_UNCERTAINTY")
            if finite(value) and value<0: reasons.append("NEGATIVE_EXTRACTED_VALUE")
            if bin_coverage is None: reasons.append("COVERAGE_BIN_MISSING")
            if curve_coverage is None: reasons.append("COVERAGE_CURVE_MISSING")
            selected=bin_coverage is not None and bin_coverage.get("Survives_all_requirements")=="1"
            approved=curve_coverage is not None and curve_coverage.get("Curve_status")=="APPROVED"
            if bin_coverage is not None and not selected: reasons.append("OUTSIDE_APPROVED_CONTIGUOUS_RANGE")
            if curve_coverage is not None and not approved: reasons.append("COVERAGE_CURVE_REJECTED")
            usable=(combination_ok and finite(value) and finite(error) and error>0
                    and mc>0 and mod>0 and selected and approved)
            row={field:"" for field in POINT_FIELDS}
            row.update({field:first.get(field,"") for field in (*PROVENANCE_FIELDS,*GROUP_FIELDS,"z")})
            row.update({
                "Pt_bin":ipt,"Pt_low":source["Pt_low"],"Pt_high":source["Pt_high"],
                "Pt_high_closed":source["Pt_high_closed"],"phipq_bin":iphi,
                "phipq_low":source["phipq_low"],"phipq_high":source["phipq_high"],
                "phipq_center":number(source["phipq_center"]),"Y_Data":dy,"Y_Data_err":de,
                "Y_MC":mc,"Y_MC_err":me,"C_Y":cy,"C_Y_err":ce,"Model_value":mod,
                "Model_units":units,"Model_Pt_GeV":(number(source["Pt_low"])+number(source["Pt_high"]))/2,
                "Model_phi_rad":model_phi,"Model_phi_policy":MODEL_POLICY,
                "Extracted_value":value,"Extracted_err":error,"Observable_units":units,
                "Coverage_binning_version":bin_coverage.get("Binning_version","") if bin_coverage else "",
                "Coverage_policy_version":bin_coverage.get("Policy_version","") if bin_coverage else "",
                "Coverage_bin_status":bin_coverage.get("Bin_status","") if bin_coverage else "MISSING",
                "Coverage_bin_reason":bin_coverage.get("Bin_reason","") if bin_coverage else "COVERAGE_BIN_MISSING",
                "Coverage_independent_pass":bin_coverage.get("Independent_pass","") if bin_coverage else "",
                "Coverage_selected_contiguous":int(selected),
                "Coverage_curve_status":curve_coverage.get("Curve_status","") if curve_coverage else "MISSING",
                "Coverage_curve_reason":curve_coverage.get("Curve_reason","") if curve_coverage else "COVERAGE_CURVE_MISSING",
                "Used_in_fit":int(usable),"Point_status":"WARNING" if reasons else "OK",
                "Point_reason":";".join(sorted(set(reasons))),
            })
            for reaction in REACTIONS:
                row[f"Y_MC_{reaction}"]=source[f"Y_MC_{reaction}"]
                row[f"Y_MC_{reaction}_err"]=source[f"Y_MC_{reaction}_err"]
            points[observable].append(row)
    return points


def z_color(z):
    return next((color for value,color in Z_COLORS.items() if abs(float(z)-value)<1e-6),ROOT.kOrange+7)


def rejection_codes(reason):
    reasons=set(filter(None,str(reason).split(";")))
    return [code for name,code in REJECTION_CODES.items() if name in reasons]


def draw_pdf(path,observable,point_rows,fit_rows):
    groups=defaultdict(list)
    for row in point_rows: groups[tuple(str(row[field]) for field in GROUP_FIELDS)].append(row)
    fits={(tuple(str(row[field]) for field in GROUP_FIELDS),str(row["z"]),int(row["Pt_bin"])):row for row in fit_rows}
    canvas=ROOT.TCanvas(f"c_{observable}","azimuth",1800,1200); canvas.SetCanvasSize(1800,1200); canvas.Print(str(path)+"[")
    text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.034)
    title="Cross Section" if observable=="sigma" else "Multiplicity"
    canvas.Clear(); text.DrawLatex(.05,.93,f"RP SIDIS {title} vs #phi_{{pq}}"); text.DrawLatex(.05,.86,f"Baseline PID + Delta-only acceptance; {len(PHI_EDGES)-1} #phi bins; five P_{{t}} bins")
    text.DrawLatex(.05,.79,"Filled points: approved contiguous coverage; open gray: excluded diagnostic points")
    text.DrawLatex(.05,.72,"Fits use only curves approved by the versioned azimuthal coverage study")
    formula="C_{0} + C_{1} cos(#phi) + C_{2} cos(2#phi)" if observable=="sigma" else "M_{0} + M_{1} cos(#phi) + M_{2} cos(2#phi)"
    text.DrawLatex(.05,.65,"Fit: "+formula)
    text.DrawLatex(.05,.58,f"Model phi policy: {MODEL_POLICY}"); canvas.Print(str(path))
    for group,members in sorted(groups.items()):
        canvas.Clear(); canvas.Divide(3,2,.002,.002); byz=defaultdict(list)
        for row in members: byz[str(row["z"])].append(row)
        owned=[]
        for ipt in range(1,6):
            pad=canvas.cd(ipt); pad.SetGridx(); pad.SetGridy(); pad.SetMargin(.15,.04,.15,.12)
            selected=[row for rows in byz.values() for row in rows if int(row["Pt_bin"])==ipt]
            reference=[row for row in selected if row["Used_in_fit"]==1 and finite(row["Extracted_value"])]
            if not reference: reference=[row for row in selected if finite(row["Extracted_value"])]
            vals=[number(row["Extracted_value"]) for row in reference]
            errs=[number(row["Extracted_err"],0) for row in reference]
            ymax=max([v+e for v,e in zip(vals,errs)]+[1e-12]); ymin=min([v-e for v,e in zip(vals,errs)]+[0.0]); span=max(ymax-ymin,1e-12)
            frame=ROOT.TH1D(f"f_{observable}_{abs(hash((group,ipt)))}",f"{PT_EDGES[ipt-1]:.1f} #leq P_{{t}} {'#leq' if ipt==5 else '<'} {PT_EDGES[ipt]:.1f} GeV;#phi_{{pq}};{title}",60,0,2*math.pi)
            frame_low=min(0,ymin-.1*span); frame_high=ymax+.2*span
            frame.SetMinimum(frame_low); frame.SetMaximum(frame_high); frame.GetXaxis().SetNdivisions(506)
            frame.SetLineColorAlpha(ROOT.kWhite,0.0); frame.Draw(); owned.append(frame)
            zero=ROOT.TLine(0,0,2*math.pi,0); zero.SetLineColor(ROOT.kViolet+1); zero.SetLineWidth(2); zero.SetLineStyle(2); zero.Draw("SAME"); owned.append(zero)
            panel_label=ROOT.TLatex(); panel_label.SetNDC(); panel_label.SetTextFont(42); panel_label.SetTextSize(.045)
            panel_label.DrawLatex(.19,.92,f"{PT_EDGES[ipt-1]:.1f} #leq P_{{t}} {'#leq' if ipt==5 else '<'} {PT_EDGES[ipt]:.1f} GeV")
            owned.append(panel_label)
            for z,rows in sorted(byz.items(),key=lambda item:float(item[0])):
                pts=sorted([row for row in rows if int(row["Pt_bin"])==ipt and finite(row["Extracted_value"])],key=lambda row:row["phipq_center"])
                if not pts: continue
                accepted=ROOT.TGraphErrors(); rejected=ROOT.TGraphErrors(); rejected_high=ROOT.TGraph(); rejected_low=ROOT.TGraph(); color=z_color(z)
                for row in pts:
                    value=number(row["Extracted_value"]); phi=number(row["phipq_center"])
                    if row["Used_in_fit"]==1:
                        i=accepted.GetN(); accepted.SetPoint(i,phi,value); accepted.SetPointError(i,0,number(row["Extracted_err"],0))
                    elif value>frame_high:
                        rejected_high.SetPoint(rejected_high.GetN(),phi,frame_high-.03*(frame_high-frame_low))
                    elif value<frame_low:
                        rejected_low.SetPoint(rejected_low.GetN(),phi,frame_low+.03*(frame_high-frame_low))
                    else:
                        i=rejected.GetN(); rejected.SetPoint(i,phi,value); rejected.SetPointError(i,0,number(row["Extracted_err"],0))
                accepted.SetMarkerStyle(20); accepted.SetMarkerColor(color); accepted.SetLineColor(color)
                if accepted.GetN(): accepted.Draw("P SAME")
                rejected.SetMarkerStyle(24); rejected.SetMarkerColor(ROOT.kGray+1); rejected.SetLineColor(ROOT.kGray+1)
                if rejected.GetN(): rejected.Draw("P SAME")
                rejected_high.SetMarkerStyle(26); rejected_low.SetMarkerStyle(32)
                for graph in (rejected_high,rejected_low):
                    graph.SetMarkerColor(ROOT.kGray+1)
                    if graph.GetN(): graph.Draw("P SAME")
                owned.extend((accepted,rejected,rejected_high,rejected_low))
                fit=fits.get((group,z,ipt),{})
                if fit.get("Fit_line_style") in {"solid","dashed"}:
                    names=("C_0","C_1","C_2") if observable=="sigma" else ("M_0","M_1","M_2")
                    func=ROOT.TF1(f"fit_{observable}_{abs(hash((group,z,ipt)))}","[0]+[1]*cos(x)+[2]*cos(2*x)",0,2*math.pi)
                    for j,name in enumerate(names): func.SetParameter(j,float(fit[name]))
                    func.SetLineColor(color); func.SetLineWidth(2); func.SetLineStyle(2 if fit["Fit_line_style"]=="dashed" else 1); func.Draw("SAME"); owned.append(func)
                elif fit.get("Coverage_curve_status")=="REJECTED":
                    codes=rejection_codes(fit.get("Coverage_curve_reason",""))
                    if codes:
                        notice=ROOT.TLatex(); notice.SetNDC(); notice.SetTextAlign(22); notice.SetTextSize(.042); notice.SetTextColorAlpha(ROOT.kGray+2,.62)
                        notice.DrawLatex(.55,.82,"Not fitted - codes: "+", ".join(map(str,codes))); owned.append(notice)
        pad=canvas.cd(6); pad.Clear(); legend=ROOT.TLegend(.12,.35,.88,.88); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.055)
        for z in sorted(byz,key=float):
            if not any(finite(row["Extracted_value"]) for row in byz[z]): continue
            marker=ROOT.TGraph(); marker.SetMarkerStyle(20); marker.SetMarkerColor(z_color(z)); legend.AddEntry(marker,f"z = {float(z):g}","p"); owned.append(marker)
        invalid=ROOT.TGraph(); invalid.SetMarkerStyle(24); invalid.SetMarkerColor(ROOT.kGray+1); legend.AddEntry(invalid,"Excluded from fit","p")
        offscale=ROOT.TGraph(); offscale.SetMarkerStyle(26); offscale.SetMarkerColor(ROOT.kGray+1); legend.AddEntry(offscale,"Excluded, off scale","p")
        legend.SetY1(.57); legend.SetY2(.90); legend.Draw()
        text2=ROOT.TLatex(); text2.SetNDC(); text2.SetTextFont(42); text2.SetTextSize(.038)
        text2.DrawLatex(.10,.51,"Fit rejection codes:"); text2.DrawLatex(.10,.44,"1  accepted #Delta#phi < #pi")
        text2.DrawLatex(.10,.37,"2  harmonic design rank < 3"); text2.DrawLatex(.10,.30,"3  weighted condition > 100")
        first=members[0]
        text2.DrawLatex(.10,.21,", ".join(group[:4])); text2.DrawLatex(.10,.15,f"x={group[4]}, Q2={group[5]}")
        settings=str(first['Contributing_setting_ids']).replace('|',', ')
        theta=str(first['Contributing_thpq_deg']).replace('|',', ')
        text2.DrawLatex(.10,.09,f"settings: {settings}")
        text2.DrawLatex(.10,.03,f"thpq [deg]: {theta}")
        owned.extend((legend,invalid,offscale,text2)); canvas.Print(str(path))
    canvas.Print(str(path)+"]"); canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1]; jlab=project.parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combined-yields",type=Path,default=project/"results/Pass1/Tables/RP_combine_thpq_yields/RP_combined_thpq_yields.csv")
    parser.add_argument("--coverage-bins",type=Path,default=project/"results/Pass1/Tables/RP_azimuthal_coverage_study/RP_azimuthal_coverage_bins.csv")
    parser.add_argument("--coverage-curves",type=Path,default=project/"results/Pass1/Tables/RP_azimuthal_coverage_study/RP_azimuthal_coverage_curves.csv")
    parser.add_argument("--calculator",type=Path,default=jlab/"simc_gfortran/util/sidisxsec/calc_semi_xsec")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_azimuthal_analysis"); parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_azimuthal_analysis")
    parser.add_argument("--group",action="append",default=[],help="Process only this combination_group_id; repeat as needed")
    parser.add_argument("--no-pdf",action="store_true"); args=parser.parse_args()
    coverage_bins,coverage_curves=load_coverage(args.coverage_bins,args.coverage_curves)
    grouped=defaultdict(list)
    for row in dmc.read_csv(args.combined_yields): grouped[row["Combination_group_id"]].append(row)
    keys=sorted(grouped)
    if args.group: keys=[key for key in keys if key in set(args.group)]
    all_points={"sigma":[],"multiplicity":[]}
    for index,key in enumerate(keys,1):
        points=extract_group(grouped[key],args.calculator,coverage_bins,coverage_curves)
        for observable in all_points: all_points[observable]+=points[observable]
        print(f"[{index}/{len(keys)}] OK: {key}")
    args.table_dir.mkdir(parents=True,exist_ok=True); args.pdf_dir.mkdir(parents=True,exist_ok=True)
    for observable,points in all_points.items():
        fit_fields=FIT_COMMON+((["C_0","C_0_err","C_1","C_1_err","C_2","C_2_err"] if observable=="sigma" else ["M_0","M_0_err","M_1","M_1_err","M_2","M_2_err"])+FIT_DIAGNOSTICS)
        curves=defaultdict(list)
        for row in points: curves[(str(row["Combination_group_id"]),int(row["Pt_bin"]))].append(row)
        fits=[]
        for key,curve in sorted(curves.items(),key=lambda item:str(item[0])):
            first=curve[0]; fit={field:"" for field in fit_fields}; fit.update({field:first[field] for field in (*PROVENANCE_FIELDS,*GROUP_FIELDS,"z","Pt_bin","Pt_low","Pt_high","Pt_high_closed","Model_value","Model_units","Model_Pt_GeV","Model_phi_policy","Coverage_binning_version","Coverage_policy_version","Coverage_curve_status","Coverage_curve_reason")}); fit["Coverage_selected_mask"]="".join("1" if row["Coverage_selected_contiguous"]==1 else "0" for row in sorted(curve,key=lambda row:int(row["phipq_bin"]))); fit.update(fit_curve(curve,observable))
            if first["Coverage_curve_status"]!="APPROVED":
                reasons=set(filter(None,str(fit["Fit_reason"]).split(";"))); reasons.add("COVERAGE_CURVE_NOT_APPROVED")
                if first["Coverage_curve_reason"]: reasons.update(filter(None,str(first["Coverage_curve_reason"]).split(";")))
                fit["Fit_reason"]=";".join(sorted(reasons)); fit["Fit_status"]="NOT_FITTED"; fit["Fit_line_style"]="none"
            fits.append(fit)
        prefix=f"RP_sidis_{observable}_vs_phipq"; atomic_csv(args.table_dir/f"{prefix}.csv",POINT_FIELDS,points); atomic_csv(args.table_dir/f"{prefix}_Problematic.csv",POINT_FIELDS,[row for row in points if row["Point_status"]!="OK"])
        atomic_csv(args.table_dir/f"{prefix}_FitParameters.csv",fit_fields,fits); atomic_csv(args.table_dir/f"{prefix}_FitParameters_Problematic.csv",fit_fields,[row for row in fits if row["Fit_status"]!="OK"])
        if not args.no_pdf:
            for z in sorted({str(row["z"]) for row in points},key=float):
                zdir="z_"+f"{float(z):.2f}".replace(".","p")
                selected_points=[row for row in points if str(row["z"])==z]
                selected_fits=[row for row in fits if str(row["z"])==z]
                output=args.pdf_dir/zdir/f"{prefix}.pdf"; output.parent.mkdir(parents=True,exist_ok=True)
                draw_pdf(output,observable,selected_points,selected_fits)
        print(observable,"points",len(points),"fits",Counter(row["Fit_status"] for row in fits))
    return 0


if __name__=="__main__": raise SystemExit(main())
