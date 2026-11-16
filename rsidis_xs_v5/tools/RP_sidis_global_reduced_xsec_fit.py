#!/usr/bin/env python3
"""Globally fit Pass-3/4/5 reduced SIDIS cross sections with four coefficients."""

from __future__ import annotations

import argparse,csv,math,os,tempfile
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import ROOT

from RP_sidis_binning import AZIMUTHAL_PT_EDGES as PT_EDGES

ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray+1);ROOT.gStyle.SetGridStyle(3);ROOT.TGaxis.SetMaxDigits(3)

POLICY_VERSION="RP_GLOBAL_REDUCED_XSEC_FIT_V1";MAX_WEIGHTED_CONDITION=100.;MAX_MC_REL_ERR=.30;MAX_RATIO_REL_ERR=.30;RATIO_WARNING_LOW=.50;RATIO_WARNING_HIGH=1.50;PASS_ORDER={"pass3":3,"pass4":4,"pass5":5}
PASS_COLORS={"pass3":ROOT.kOrange+7,"pass4":ROOT.kGreen+2,"pass5":ROOT.kBlue+1}
IDENTITY=("Epsilon_family_id","Epsilon_family_policy_version","Family_status","Family_reason","Run_type","Target","x","Q2","z","N_members","Pass_pattern","Has_pass3","Has_pass4","Has_pass5","N_distinct_epsilon","Epsilon_min","Epsilon_max","Epsilon_span")
POINT_FIELDS=[*IDENTITY,"Member_phase","Member_pass","Combination_group_id","epsilon","Pt_bin","Pt_low","Pt_high","Pt_high_closed","phipq_bin","phipq_center","Reduced_data_xsec","Reduced_data_xsec_err","Global_point_usable","Global_point_status","Global_point_reason","Global_fit_policy_version"]
FIT_FIELDS=[*IDENTITY,"Pt_bin","Pt_low","Pt_high","Pt_high_closed","N_points_total","N_points_fitted","N_pass3","N_pass4","N_pass5","Weighted_design_rank","Weighted_design_condition","C_0_global","C_0_global_err","C_1_global","C_1_global_err","C_2_global","C_2_global_err","C_3_global","C_3_global_err","F_T","F_T_err","F_L","F_L_err","F_UU_cosphi","F_UU_cosphi_err","F_UU_cos2phi","F_UU_cos2phi_err",*[f"Cov_{i}{j}" for i in range(4) for j in range(i,4)],*[f"Correlation_{i}{j}" for i in range(4) for j in range(i+1,4)],"Chi2","NDF","Chi2_NDF","Fit_probability","Fit_status","Fit_reason","Fit_policy_version"]


def finite(value):
    try:return math.isfinite(float(value))
    except (TypeError,ValueError):return False


def number(value,default=math.nan):return float(value) if finite(value) else default


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as stream:return list(csv.DictReader(stream))


def atomic_csv(path,fields,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);descriptor,temporary=tempfile.mkstemp(prefix=path.name,suffix=".tmp",dir=path.parent)
    try:
        with os.fdopen(descriptor,"w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def design_row(epsilon,phi):
    epsilon=number(epsilon);phi=number(phi)
    return np.array([1.,epsilon,math.sqrt(2*epsilon*(1+epsilon))*math.cos(phi),epsilon*math.cos(2*phi)])


def normalize_point(source):
    reasons=[];data=number(source.get("Y_Data"));data_error=number(source.get("Y_Data_err"));mc=number(source.get("Y_MC"));mc_error=number(source.get("Y_MC_err"));ratio=number(source.get("Data_by_MC"));ratio_error=number(source.get("Data_by_MC_err"))
    mc_rel=mc_error/mc if mc>0 and mc_error>=0 else math.inf;ratio_rel=ratio_error/abs(ratio) if finite(ratio) and ratio!=0 and ratio_error>=0 else math.inf
    reduction_ok=source.get("Reduction_status") in {"OK","WARNING"};mc_pass=mc>0 and mc_rel<=MAX_MC_REL_ERR;ratio_pass=data>0 and finite(data_error) and data_error>=0 and ratio_rel<=MAX_RATIO_REL_ERR
    value=number(source.get("Reduced_data_xsec"));error=number(source.get("Reduced_data_xsec_err"))
    if not reduction_ok:reasons.append("REDUCTION_POINT_NOT_USABLE")
    if not mc_pass:reasons.append("MC_SUPPORT_FAILED")
    if not ratio_pass:reasons.append("RATIO_PRECISION_OR_POSITIVE_DATA_FAILED")
    if not finite(value):reasons.append("NONFINITE_REDUCED_VALUE")
    if not finite(error) or error<=0:reasons.append("NONPOSITIVE_REDUCED_UNCERTAINTY")
    if finite(ratio) and not RATIO_WARNING_LOW<=ratio<=RATIO_WARNING_HIGH:reasons.append("RATIO_OUTSIDE_WARNING_RANGE")
    usable=reduction_ok and mc_pass and ratio_pass and finite(value) and finite(error) and error>0
    row={field:source.get(field,"") for field in POINT_FIELDS};row.update(Global_point_usable=int(usable),Global_point_status="USABLE" if usable else "EXCLUDED",Global_point_reason=";".join(sorted(set(reasons))),Global_fit_policy_version=POLICY_VERSION)
    return row


def fit_global(points):
    first=points[0];result={field:"" for field in FIT_FIELDS}
    for field in (*IDENTITY,"Pt_bin","Pt_low","Pt_high","Pt_high_closed"):result[field]=first.get(field,"")
    usable=[row for row in points if int(row["Global_point_usable"])==1];counts=Counter(row["Member_pass"] for row in usable)
    result.update(N_points_total=len(points),N_points_fitted=len(usable),N_pass3=counts["pass3"],N_pass4=counts["pass4"],N_pass5=counts["pass5"],Fit_policy_version=POLICY_VERSION)
    if len(usable)<4:
        result.update(Fit_status="NOT_FITTED",Fit_reason="INSUFFICIENT_GLOBAL_POINTS");return result
    design=np.vstack([design_row(row["epsilon"],row["phipq_center"]) for row in usable]);values=np.asarray([number(row["Reduced_data_xsec"]) for row in usable]);errors=np.asarray([number(row["Reduced_data_xsec_err"]) for row in usable])
    weighted=design/errors[:,None];rank=int(np.linalg.matrix_rank(weighted));condition=float(np.linalg.cond(weighted));result.update(Weighted_design_rank=rank,Weighted_design_condition=condition)
    if rank<4:
        result.update(Fit_status="NOT_FITTED",Fit_reason="GLOBAL_DESIGN_RANK_BELOW_4");return result
    if condition>MAX_WEIGHTED_CONDITION:
        result.update(Fit_status="NOT_FITTED",Fit_reason="GLOBAL_WEIGHTED_CONDITION_ABOVE_100");return result
    normal=weighted.T@weighted;covariance=np.linalg.inv(normal);parameters=covariance@(design.T@(values/errors**2));prediction=design@parameters
    chi2=float(np.sum(((values-prediction)/errors)**2));ndf=len(usable)-4;probability=float(ROOT.TMath.Prob(chi2,ndf)) if ndf>0 else math.nan
    for i,name in enumerate(("C_0_global","C_1_global","C_2_global","C_3_global")):
        result[name]=parameters[i];result[name+"_err"]=math.sqrt(max(0,covariance[i,i]))
    for name,index in (("F_T",0),("F_L",1),("F_UU_cosphi",2),("F_UU_cos2phi",3)):
        result[name]=parameters[index];result[name+"_err"]=math.sqrt(max(0,covariance[index,index]))
    for i in range(4):
        for j in range(i,4):result[f"Cov_{i}{j}"]=covariance[i,j]
    for i in range(4):
        for j in range(i+1,4):result[f"Correlation_{i}{j}"]=covariance[i,j]/math.sqrt(covariance[i,i]*covariance[j,j])
    reasons=[]
    if ndf>0 and probability<.01:reasons.append("POOR_FIT_PROBABILITY")
    if parameters[0]<0:reasons.append("NEGATIVE_FT")
    if parameters[1]<0:reasons.append("NEGATIVE_FL")
    result.update(Chi2=chi2,NDF=ndf,Chi2_NDF=chi2/ndf if ndf else "",Fit_probability=probability if ndf else "",Fit_status="WARNING" if reasons else "OK",Fit_reason=";".join(reasons));return result


def graph(rows,usable,color,marker):
    result=ROOT.TGraphErrors();i=0
    for row in rows:
        if bool(int(row["Global_point_usable"]))!=usable or not finite(row["Reduced_data_xsec"]):continue
        result.SetPoint(i,number(row["phipq_center"]),number(row["Reduced_data_xsec"]));result.SetPointError(i,0,number(row["Reduced_data_xsec_err"],0));i+=1
    result.SetMarkerStyle(marker);result.SetMarkerColor(color);result.SetLineColor(color);return result


def model_value(fit,epsilon,phi):
    coefficients=np.array([number(fit[f"C_{i}_global"]) for i in range(4)]);return float(design_row(epsilon,phi)@coefficients)


def draw_pdf(path,rows,fits):
    path.parent.mkdir(parents=True,exist_ok=True);canvas=ROOT.TCanvas("c_global","global",1800,1200);canvas.SetCanvasSize(1800,1200);canvas.Print(str(path)+"[")
    byfamily=defaultdict(list);fitfamily=defaultdict(list)
    for row in rows:byfamily[row["Epsilon_family_id"]].append(row)
    for row in fits:fitfamily[row["Epsilon_family_id"]].append(row)
    for family,points in sorted(byfamily.items()):
        canvas.Clear();canvas.Divide(3,2,.002,.002);owned=[];fitmap={int(r["Pt_bin"]):r for r in fitfamily[family]};bypass=defaultdict(list)
        for row in points:bypass[(row["Member_pass"],int(row["Pt_bin"]))].append(row)
        for pt in range(1,6):
            pad=canvas.cd(pt);pad.SetGridx();pad.SetGridy();pad.SetMargin(.15,.04,.15,.11);allrows=[r for pas in ("pass3","pass4","pass5") for r in bypass.get((pas,pt),[])]
            scale=[r for r in allrows if int(r["Global_point_usable"])] or [r for r in allrows if finite(r["Reduced_data_xsec"])]
            lows=[number(r["Reduced_data_xsec"])-number(r["Reduced_data_xsec_err"],0) for r in scale];highs=[number(r["Reduced_data_xsec"])+number(r["Reduced_data_xsec_err"],0) for r in scale];low=min(lows+[0]);high=max(highs+[1e-12]);span=max(high-low,1e-12)
            frame=ROOT.TH1D(f"gf_{abs(hash((family,pt)))}","",60,0,2*math.pi);frame.SetMinimum(min(0,low-.1*span));frame.SetMaximum(high+.18*span);frame.SetLineColorAlpha(ROOT.kWhite,0);frame.GetXaxis().SetTitle("#phi_{pq} [rad]");frame.GetYaxis().SetTitle("#sigma_{red} [GeV^{-2}]");frame.Draw()
            label=ROOT.TLatex();label.SetNDC();label.SetTextSize(.043);relation="#leq" if pt==5 else "<";label.DrawLatex(.20,.92,f"{PT_EDGES[pt-1]:.1f} #leq P_{{T}} {relation} {PT_EDGES[pt]:.1f} GeV")
            zero=ROOT.TLine(0,0,2*math.pi,0);zero.SetLineColor(ROOT.kViolet+1);zero.SetLineStyle(2);zero.Draw();owned.extend((frame,label,zero));fit=fitmap[pt]
            for pas in ("pass3","pass4","pass5"):
                selected=bypass.get((pas,pt),[])
                if not selected:continue
                good=graph(selected,True,PASS_COLORS[pas],20);bad=graph(selected,False,ROOT.kGray+2,24)
                if bad.GetN():bad.Draw("P SAME")
                if good.GetN():good.Draw("P SAME")
                if fit["Fit_status"] in {"OK","WARNING"}:
                    epsilon=number(selected[0]["epsilon"]);func=ROOT.TF1(f"gfit_{abs(hash((family,pas,pt)))}","[0]+[1]*cos(x)+[2]*cos(2*x)",0,2*math.pi)
                    func.SetParameters(number(fit["F_T"])+epsilon*number(fit["F_L"]),math.sqrt(2*epsilon*(1+epsilon))*number(fit["F_UU_cosphi"]),epsilon*number(fit["F_UU_cos2phi"]));func.SetLineColor(PASS_COLORS[pas]);func.SetLineWidth(2);func.Draw("SAME");owned.append(func)
                owned.extend((good,bad))
            if fit["Fit_status"]=="NOT_FITTED":
                notice=ROOT.TLatex();notice.SetNDC();notice.SetTextAlign(22);notice.SetTextSize(.04);notice.SetTextColorAlpha(ROOT.kGray+2,.65);notice.DrawLatex(.55,.82,"Global fit unavailable");owned.append(notice)
        first=points[0];pad=canvas.cd(6);pad.Clear();text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.038)
        lines=["Global reduced-xsec fit",f"{first['Run_type']} / {first['Target']}",f"x={number(first['x']):.3g}, Q^{{2}}={number(first['Q2']):.3g} GeV^{{2}}, z={number(first['z']):.3g}",f"Family: {family}","C_{0}=F_{T}, C_{1}=F_{L}","C_{2}=F_{UU}^{cos#phi}, C_{3}=F_{UU}^{cos2#phi}"]
        for i,line in enumerate(lines):text.DrawLatex(.05,.94-.072*i,line)
        y=.42
        for pas,color in (("pass3",ROOT.kOrange+7),("pass4",ROOT.kGreen+2),("pass5",ROOT.kBlue+1)):
            if any(r["Member_pass"]==pas for r in points):
                marker=ROOT.TMarker(.09,y,20);marker.SetNDC(True);marker.SetMarkerColor(color);marker.Draw();text.DrawLatex(.14,y-.015,pas+" points + global curve");owned.append(marker);y-=.075
        text.SetTextSize(.030);text.DrawLatex(.05,.09,"All passes contribute directly to one fit");text.DrawLatex(.05,.04,POLICY_VERSION);owned.append(text);canvas._global_owned=owned;canvas.Print(str(path))
    canvas.Print(str(path)+"]");canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced-points",type=Path,default=project/"results/Pass1/Tables/RP_extract_sidis_reduced_xsec/RP_sidis_reduced_xsec_points.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_global_reduced_xsec_fit")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_global_reduced_xsec_fit");parser.add_argument("--no-pdf",action="store_true");args=parser.parse_args()
    points=[normalize_point(r) for r in read_csv(args.reduced_points) if r["Family_status"].startswith("LT_READY")];groups=defaultdict(list)
    for row in points:groups[(row["Epsilon_family_id"],int(row["Pt_bin"]))].append(row)
    fits=[fit_global(rows) for _,rows in sorted(groups.items())]
    atomic_csv(args.table_dir/"RP_sidis_global_reduced_xsec_points.csv",POINT_FIELDS,points);atomic_csv(args.table_dir/"RP_sidis_global_reduced_xsec_points_Problematic.csv",POINT_FIELDS,[r for r in points if r["Global_point_status"]!="USABLE"])
    atomic_csv(args.table_dir/"RP_sidis_global_reduced_xsec_fit_parameters.csv",FIT_FIELDS,fits);atomic_csv(args.table_dir/"RP_sidis_global_reduced_xsec_fit_parameters_Problematic.csv",FIT_FIELDS,[r for r in fits if r["Fit_status"]!="OK"])
    if not args.no_pdf:
        for z in sorted({r["z"] for r in points},key=float):
            selected_points=[r for r in points if r["z"]==z];selected_fits=[r for r in fits if r["z"]==z];zdir="z_"+f"{float(z):.2f}".replace(".","p");draw_pdf(args.pdf_dir/zdir/"RP_sidis_global_reduced_xsec_fit.pdf",selected_points,selected_fits)
    print("Points:",len(points),Counter(r["Global_point_status"] for r in points));print("Fits:",len(fits),Counter(r["Fit_status"] for r in fits));return 0


if __name__=="__main__":raise SystemExit(main())
