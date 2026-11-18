#!/usr/bin/env python3
"""Fit lower-pass virtual-photon C0 with Pass-5-anchored harmonics."""

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
ROOT.gStyle.SetGridColor(ROOT.kGray+1); ROOT.gStyle.SetGridStyle(3); ROOT.TGaxis.SetMaxDigits(3)

POLICY_VERSION="RP_FIT_VIRTUAL_PHOTON_XSEC_LOWER_PASS_V2"
MAX_MC_REL_ERR=.30; MAX_RATIO_REL_ERR=.30; RATIO_WARNING_LOW=.50; RATIO_WARNING_HIGH=1.50
MIN_POINTS=2
PASS_ORDER={"pass3":3,"pass4":4,"pass5":5}
PASS_COLORS={"pass3":ROOT.kOrange+7,"pass4":ROOT.kGreen+2,"pass5":ROOT.kBlue+1}
IDENTITY=("Epsilon_family_id","Epsilon_family_policy_version","Family_status","Family_reason",
          "Run_type","Target","x","Q2","z","N_members","Pass_pattern","Has_pass3",
          "Has_pass4","Has_pass5","N_distinct_epsilon","Epsilon_min","Epsilon_max","Epsilon_span")
PROVENANCE=("Combination_group_id","Combination_policy_version","Ebeam_GeV","N_contributors",
            "Contributing_setting_ids","Contributing_thpq_deg")
POINT_EXTRA=("Lower_fit_policy_version","MC_rel_err","Ratio_rel_err","MC_support_pass",
             "Ratio_precision_pass","Ratio_warning","Anchor_fit_status","Anchor_fit_reason",
             "Used_in_lower_fit","Lower_point_status","Lower_point_reason")
FIT_FIELDS=[*IDENTITY,"Member_phase","Member_pass",*PROVENANCE,"epsilon","Pt_bin","Pt_low",
            "Pt_high","Pt_high_closed","Anchor_group_id","Anchor_epsilon","Anchor_fit_status",
            "Anchor_fit_reason","Convention_version","Cross_section_level",
            "Virtual_photon_xsec_units","Coefficient_units","C_0_interpretation",
            "C_1_interpretation","C_2_interpretation","C_1_pass5","C_2_pass5","Scale_C1","Scale_C2","C_1_fixed",
            "C_1_fixed_err","C_2_fixed","C_2_fixed_err","Cov_fixed_12","N_points_available",
            "N_points_fitted","Phi_span_observed_rad","C_0","C_0_err","C_0_lower_stat_err",
            "C_0_anchor_err","dC0_dC1_pass5","dC0_dC2_pass5","Covariance_rank",
            "Covariance_condition","Chi2","NDF","Chi2_NDF","Fit_probability","Fit_status",
            "Fit_reason","Fit_policy_version"]


def finite(value):
    try:return math.isfinite(float(value))
    except (TypeError,ValueError):return False


def number(value,default=math.nan):return float(value) if finite(value) else default


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as stream:return list(csv.DictReader(stream))


def atomic_csv(path,fields,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    descriptor,temporary=tempfile.mkstemp(prefix=path.name,suffix=".tmp",dir=path.parent)
    try:
        with os.fdopen(descriptor,"w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def harmonic_scales(epsilon,anchor_epsilon):
    epsilon=number(epsilon); anchor_epsilon=number(anchor_epsilon)
    if not 0<epsilon<1 or not 0<anchor_epsilon<1:raise ValueError("epsilon must lie in (0,1)")
    scale1=math.sqrt(2*epsilon*(1+epsilon))/math.sqrt(2*anchor_epsilon*(1+anchor_epsilon))
    return scale1,epsilon/anchor_epsilon


def classify_point(source,anchor):
    reasons=[]; data=number(source.get("Y_Data")); data_error=number(source.get("Y_Data_err"))
    mc=number(source.get("Y_MC")); mc_error=number(source.get("Y_MC_err"))
    ratio=number(source.get("Data_by_MC")); ratio_error=number(source.get("Data_by_MC_err"))
    mc_rel=mc_error/mc if mc>0 and mc_error>=0 else math.inf
    ratio_rel=ratio_error/abs(ratio) if finite(ratio) and ratio!=0 and ratio_error>=0 else math.inf
    reduction_ok=source.get("Reduction_status") in {"OK","WARNING"}
    mc_pass=mc>0 and mc_rel<=MAX_MC_REL_ERR
    ratio_pass=data>0 and finite(data_error) and data_error>=0 and ratio_rel<=MAX_RATIO_REL_ERR
    if not reduction_ok:reasons.append("REDUCTION_POINT_NOT_USABLE")
    if not mc_pass:reasons.append("MC_SUPPORT_FAILED")
    if not ratio_pass:reasons.append("RATIO_PRECISION_OR_POSITIVE_DATA_FAILED")
    value=number(source.get("Virtual_photon_data_xsec")); error=number(source.get("Virtual_photon_data_xsec_err"))
    if not finite(value):reasons.append("NONFINITE_VIRTUAL_PHOTON_VALUE")
    if not finite(error) or error<=0:reasons.append("NONPOSITIVE_VIRTUAL_PHOTON_UNCERTAINTY")
    anchor_ok=bool(anchor and anchor.get("Fit_status") in {"OK","WARNING"})
    if not anchor_ok:reasons.append("PASS5_ANCHOR_UNAVAILABLE")
    warning=finite(ratio) and not RATIO_WARNING_LOW<=ratio<=RATIO_WARNING_HIGH
    if warning:reasons.append("RATIO_OUTSIDE_WARNING_RANGE")
    usable=reduction_ok and mc_pass and ratio_pass and finite(value) and finite(error) and error>0 and anchor_ok
    return {"MC_rel_err":mc_rel,"Ratio_rel_err":ratio_rel,"MC_support_pass":int(mc_pass),
            "Ratio_precision_pass":int(ratio_pass),"Ratio_warning":int(warning),
            "Anchor_fit_status":anchor.get("Fit_status","") if anchor else "MISSING",
            "Anchor_fit_reason":anchor.get("Fit_reason","") if anchor else "PASS5_ANCHOR_MISSING",
            "Used_in_lower_fit":int(usable),"Lower_point_status":"USABLE" if usable else "EXCLUDED",
            "Lower_point_reason":";".join(sorted(set(reasons))),"Lower_fit_policy_version":POLICY_VERSION}


def annotate_points(points,anchors):
    output=[]
    for source in points:
        anchor=anchors.get((source["Epsilon_family_id"],int(source["Pt_bin"])))
        row=dict(source); row.update(classify_point(source,anchor)); output.append(row)
    return output


def fit_lower_curve(points,anchor):
    first=points[0]; result={field:"" for field in FIT_FIELDS}
    for field in (*IDENTITY,"Member_phase","Member_pass",*PROVENANCE,"epsilon","Pt_bin","Pt_low","Pt_high","Pt_high_closed"):
        result[field]=first.get(field,"")
    result.update(Fit_policy_version=POLICY_VERSION,
                  Convention_version=first.get("Convention_version",""),
                  Cross_section_level=first.get("Cross_section_level",""),
                  Virtual_photon_xsec_units=first.get("Virtual_photon_xsec_units",""),
                  Coefficient_units=first.get("Virtual_photon_xsec_units",""),
                  C_0_interpretation="sigma_T+epsilon_p*sigma_L",
                  C_1_interpretation="sqrt(2*epsilon_p*(1+epsilon_p))*sigma_LT",
                  C_2_interpretation="epsilon_p*sigma_TT",
                  N_points_available=sum(finite(r.get("Virtual_photon_data_xsec")) for r in points))
    if not anchor or anchor.get("Fit_status") not in {"OK","WARNING"}:
        result.update(Anchor_fit_status=anchor.get("Fit_status","") if anchor else "MISSING",
                      Anchor_fit_reason=anchor.get("Fit_reason","") if anchor else "PASS5_ANCHOR_MISSING",
                      Fit_status="NOT_FITTED",Fit_reason="PASS5_ANCHOR_UNAVAILABLE");return result
    scale1,scale2=harmonic_scales(first["epsilon"],anchor["epsilon"])
    c1=number(anchor["C_1_pass5"]); c2=number(anchor["C_2_pass5"])
    cov11=number(anchor["Cov_11"]); cov12=number(anchor["Cov_12"]); cov22=number(anchor["Cov_22"])
    fixed_cov=np.array([[scale1*scale1*cov11,scale1*scale2*cov12],
                        [scale1*scale2*cov12,scale2*scale2*cov22]])
    result.update(Anchor_group_id=anchor["Combination_group_id"],Anchor_epsilon=anchor["epsilon"],
                  Anchor_fit_status=anchor["Fit_status"],Anchor_fit_reason=anchor["Fit_reason"],
                  C_1_pass5=c1,C_2_pass5=c2,Scale_C1=scale1,Scale_C2=scale2,
                  C_1_fixed=scale1*c1,C_1_fixed_err=math.sqrt(max(0,fixed_cov[0,0])),
                  C_2_fixed=scale2*c2,C_2_fixed_err=math.sqrt(max(0,fixed_cov[1,1])),
                  Cov_fixed_12=fixed_cov[0,1])
    usable=[row for row in points if int(row["Used_in_lower_fit"])==1]; n=len(usable)
    result["N_points_fitted"]=n
    if n<MIN_POINTS:
        result.update(Fit_status="NOT_FITTED",Fit_reason="INSUFFICIENT_LOWER_PASS_POINTS");return result
    phi=np.asarray([number(row["phipq_center"]) for row in usable]); values=np.asarray([number(row["Virtual_photon_data_xsec"]) for row in usable])
    errors=np.asarray([number(row["Virtual_photon_data_xsec_err"]) for row in usable]); harmonics=np.column_stack((np.cos(phi),np.cos(2*phi)))
    fixed=harmonics@np.array([scale1*c1,scale2*c2]); covariance=np.diag(errors**2)+harmonics@fixed_cov@harmonics.T
    rank=int(np.linalg.matrix_rank(covariance)); condition=float(np.linalg.cond(covariance))
    result.update(Covariance_rank=rank,Covariance_condition=condition,
                  Phi_span_observed_rad=float(phi.max()-phi.min()) if n>1 else 0)
    if rank<n:
        result.update(Fit_status="ERROR",Fit_reason="SINGULAR_TOTAL_COVARIANCE");return result
    inverse=np.linalg.inv(covariance); ones=np.ones(n); residual_base=values-fixed
    denominator=float(ones@inverse@ones);weights=(ones@inverse)/denominator;c0=float(weights@residual_base)
    sensitivity=-(weights@harmonics)@np.diag([scale1,scale2]);anchor_covariance=np.array([[cov11,cov12],[cov12,cov22]])
    lower_variance=float(weights@np.diag(errors**2)@weights);anchor_variance=float(sensitivity@anchor_covariance@sensitivity)
    c0_error=math.sqrt(max(0,lower_variance+anchor_variance))
    residual=residual_base-c0; chi2=float(residual@inverse@residual); ndf=n-1
    probability=float(ROOT.TMath.Prob(chi2,ndf)) if ndf>0 else math.nan; reasons=[]
    if c0<=0:reasons.append("NONPOSITIVE_C0")
    if ndf>0 and probability<.01:reasons.append("POOR_FIT_PROBABILITY")
    result.update(C_0=c0,C_0_err=c0_error,C_0_lower_stat_err=math.sqrt(max(0,lower_variance)),
                  C_0_anchor_err=math.sqrt(max(0,anchor_variance)),dC0_dC1_pass5=sensitivity[0],
                  dC0_dC2_pass5=sensitivity[1],Chi2=chi2,NDF=ndf,Chi2_NDF=chi2/ndf if ndf else "",
                  Fit_probability=probability if ndf else "",Fit_status="WARNING" if reasons else "OK",
                  Fit_reason=";".join(reasons));return result


def graph(rows,usable,color,marker):
    result=ROOT.TGraphErrors();i=0
    for row in rows:
        if bool(int(row["Used_in_lower_fit"]))!=usable or not finite(row.get("Virtual_photon_data_xsec")):continue
        result.SetPoint(i,number(row["phipq_center"]),number(row["Virtual_photon_data_xsec"]));result.SetPointError(i,0,number(row["Virtual_photon_data_xsec_err"],0));i+=1
    result.SetMarkerStyle(marker);result.SetMarkerColor(color);result.SetLineColor(color);return result


def draw_pdf(path,families,anchor_fits):
    path.parent.mkdir(parents=True,exist_ok=True);canvas=ROOT.TCanvas("c_lower","lower",1800,1200);canvas.SetCanvasSize(1800,1200);canvas.Print(str(path)+"[")
    for family,points,fits in families:
        canvas.Clear();canvas.Divide(3,2,.002,.002);owned=[];by=defaultdict(list)
        for row in points:by[(row["Member_pass"],int(row["Pt_bin"]))].append(row)
        fitmap={(row["Member_pass"],int(row["Pt_bin"])):row for row in fits}
        for pt in range(1,6):
            pad=canvas.cd(pt);pad.SetGridx();pad.SetGridy();pad.SetMargin(.15,.04,.15,.11)
            allrows=[row for pas in ("pass3","pass4") for row in by.get((pas,pt),[])]
            usable=[r for r in allrows if int(r["Used_in_lower_fit"])];scale=usable or [r for r in allrows if finite(r.get("Virtual_photon_data_xsec"))]
            lows=[number(r["Virtual_photon_data_xsec"])-number(r["Virtual_photon_data_xsec_err"],0) for r in scale];highs=[number(r["Virtual_photon_data_xsec"])+number(r["Virtual_photon_data_xsec_err"],0) for r in scale]
            low=min(lows+[0]);high=max(highs+[1e-12]);span=max(high-low,1e-12);frame=ROOT.TH1D(f"lf_{abs(hash((family,pt)))}","",60,0,2*math.pi)
            frame.SetMinimum(min(0,low-.1*span));frame.SetMaximum(high+.18*span);frame.SetLineColorAlpha(ROOT.kWhite,0);frame.GetXaxis().SetTitle("#phi_{pq} [rad]");frame.GetYaxis().SetTitle("#sigma [#mub/(GeV sr)]");frame.Draw()
            label=ROOT.TLatex();label.SetNDC();label.SetTextSize(.043);relation="#leq" if pt==5 else "<";label.DrawLatex(.20,.92,f"{PT_EDGES[pt-1]:.1f} #leq P_{{T}} {relation} {PT_EDGES[pt]:.1f} GeV")
            zero=ROOT.TLine(0,0,2*math.pi,0);zero.SetLineColor(ROOT.kViolet+1);zero.SetLineStyle(2);zero.Draw();owned.extend((frame,label,zero))
            not_fitted=[]
            for pas in ("pass3","pass4"):
                rows=by.get((pas,pt),[])
                if not rows:continue
                good=graph(rows,True,PASS_COLORS[pas],20);bad=graph(rows,False,ROOT.kGray+2,24)
                if bad.GetN():bad.Draw("P SAME")
                if good.GetN():good.Draw("P SAME")
                fit=fitmap[(pas,pt)]
                if fit["Fit_status"] in {"OK","WARNING"}:
                    func=ROOT.TF1(f"lfit_{abs(hash((family,pas,pt)))}","[0]+[1]*cos(x)+[2]*cos(2*x)",0,2*math.pi);func.SetParameters(number(fit["C_0"]),number(fit["C_1_fixed"]),number(fit["C_2_fixed"]));func.SetLineColor(PASS_COLORS[pas]);func.SetLineWidth(2);func.Draw("SAME");owned.append(func)
                else:not_fitted.append(pas)
                owned.extend((good,bad))
            if not_fitted:
                notice=ROOT.TLatex();notice.SetNDC();notice.SetTextAlign(22);notice.SetTextSize(.038);notice.SetTextColorAlpha(ROOT.kGray+2,.65)
                notice.DrawLatex(.55,.82,"Not fitted: "+", ".join(not_fitted));owned.append(notice)
        first=points[0];pad=canvas.cd(6);pad.Clear();text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.038)
        lines=["Pass-5 anchored virtual-photon fits",f"{first['Run_type']} / {first['Target']}",f"x={number(first['x']):.3g}, Q^{{2}}={number(first['Q2']):.3g} GeV^{{2}}, z={number(first['z']):.3g}",f"Family: {family}","Only C_{0}^{(p)} is fitted"]
        for i,line in enumerate(lines):text.DrawLatex(.05,.93-.075*i,line)
        text.SetTextSize(.031)
        text.DrawLatex(.05,.52,"C_{1}^{(p)} = #sqrt{#frac{#epsilon_{p}(1+#epsilon_{p})}{#epsilon_{5}(1+#epsilon_{5})}} C_{1}^{(5)}")
        text.DrawLatex(.05,.44,"C_{2}^{(p)} = #frac{#epsilon_{p}}{#epsilon_{5}} C_{2}^{(5)}")
        text.SetTextSize(.038);y=.31
        for pas in ("pass3","pass4"):
            if any(r["Member_pass"]==pas for r in points):
                marker=ROOT.TMarker(.09,y,20);marker.SetNDC(True);marker.SetMarkerColor(PASS_COLORS[pas]);marker.Draw();text.DrawLatex(.14,y-.015,pas);owned.append(marker);y-=.08
        text.SetTextSize(.030);text.DrawLatex(.05,.08,"Gray open points: excluded by point policy");owned.append(text);canvas._lower_owned=owned;canvas.Print(str(path))
    canvas.Print(str(path)+"]");canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced-points",type=Path,default=project/"results/Pass1/Tables/RP_extract_sidis_reduced_xsec/RP_sidis_reduced_xsec_points.csv")
    parser.add_argument("--pass5-fits",type=Path,default=project/"results/Pass1/Tables/RP_fit_sidis_reduced_xsec_pass5/RP_sidis_virtual_photon_pass5_fit_parameters.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_fit_sidis_reduced_xsec_lower_passes")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_fit_sidis_reduced_xsec_lower_passes");parser.add_argument("--no-pdf",action="store_true");args=parser.parse_args()
    anchors={(r["Epsilon_family_id"],int(r["Pt_bin"])):r for r in read_csv(args.pass5_fits)}
    raw=[r for r in read_csv(args.reduced_points) if r["Member_pass"] in {"pass3","pass4"} and r["Family_status"].startswith("LT_READY")]
    points=annotate_points(raw,anchors);groups=defaultdict(list)
    for row in points:groups[(row["Epsilon_family_id"],row["Member_pass"],int(row["Pt_bin"]))].append(row)
    fits=[]
    for (family,pas,pt),rows in sorted(groups.items(),key=lambda item:(item[0][0],PASS_ORDER[item[0][1]],item[0][2])):
        fits.append(fit_lower_curve(sorted(rows,key=lambda r:int(r["phipq_bin"])),anchors.get((family,pt))))
    point_fields=list(points[0]);atomic_csv(args.table_dir/"RP_sidis_virtual_photon_lower_pass_points.csv",point_fields,points);atomic_csv(args.table_dir/"RP_sidis_virtual_photon_lower_pass_points_Problematic.csv",point_fields,[r for r in points if r["Lower_point_status"]!="USABLE"])
    atomic_csv(args.table_dir/"RP_sidis_virtual_photon_lower_pass_fit_parameters.csv",FIT_FIELDS,fits);atomic_csv(args.table_dir/"RP_sidis_virtual_photon_lower_pass_fit_parameters_Problematic.csv",FIT_FIELDS,[r for r in fits if r["Fit_status"]!="OK"])
    if not args.no_pdf:
        byfamily=defaultdict(list);fitfamily=defaultdict(list)
        for row in points:byfamily[row["Epsilon_family_id"]].append(row)
        for row in fits:fitfamily[row["Epsilon_family_id"]].append(row)
        for z in sorted({r["z"] for r in points},key=float):
            families=[(key,byfamily[key],fitfamily[key]) for key in sorted(byfamily) if byfamily[key][0]["z"]==z]
            zdir="z_"+f"{float(z):.2f}".replace(".","p");draw_pdf(args.pdf_dir/zdir/"RP_sidis_virtual_photon_lower_passes_vs_phipq.pdf",families,anchors)
    print("Points:",len(points),Counter(r["Lower_point_status"] for r in points));print("Fits:",len(fits),Counter(r["Fit_status"] for r in fits));return 0


if __name__=="__main__":raise SystemExit(main())
