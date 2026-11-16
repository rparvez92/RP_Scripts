#!/usr/bin/env python3
"""Perform Pass-4/Pass-5 L/T separation with optional Pass-3 validation."""

from __future__ import annotations

import argparse,csv,math,os,tempfile
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import ROOT

from RP_sidis_binning import AZIMUTHAL_PT_EDGES as PT_EDGES

ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray+1);ROOT.gStyle.SetGridStyle(3);ROOT.TGaxis.SetMaxDigits(3)

POLICY_VERSION="RP_PASS5_ANCHORED_LT_V1"
VALID_FITS={"OK","WARNING"}
IDENTITY=("Epsilon_family_id","Epsilon_family_policy_version","Family_status","Family_reason",
          "Run_type","Target","x","Q2","z","N_members","Pass_pattern","Has_pass3",
          "Has_pass4","Has_pass5","N_distinct_epsilon","Epsilon_min","Epsilon_max","Epsilon_span")
FIELDS=[*IDENTITY,"Pt_bin","Pt_low","Pt_high","Pt_high_closed","Pass3_fit_status",
        "Pass4_fit_status","Pass5_fit_status","epsilon_pass3","epsilon_pass4","epsilon_pass5",
        "C_0_pass3","C_0_pass3_err","C_0_pass4","C_0_pass4_err","C_0_pass5","C_0_pass5_err",
        "Cov_C0_pass4_pass5","Cov_C0_pass3_pass4","Cov_C0_pass3_pass5","F_T","F_T_err",
        "F_L","F_L_err","Cov_FT_FL","Correlation_FT_FL","Pass3_predicted_C0",
        "Pass3_prediction_err","Pass3_validation_residual","Pass3_validation_err",
        "Pass3_validation_pull","Pass3_validation_status","Separation_method","Separation_status",
        "Separation_reason","Separation_policy_version"]


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


def harmonic_covariance(anchor):
    return np.array([[number(anchor["Cov_11"]),number(anchor["Cov_12"])],
                     [number(anchor["Cov_12"]),number(anchor["Cov_22"])]])


def lower_gradient(row):
    return np.array([number(row["dC0_dC1_pass5"]),number(row["dC0_dC2_pass5"])])


def lower_anchor_c0_covariance(lower,anchor):
    gradient=lower_gradient(lower)
    return float(gradient@np.array([number(anchor["Cov_01"]),number(anchor["Cov_02"])]))


def lower_lower_covariance(left,right,anchor):
    return float(lower_gradient(left)@harmonic_covariance(anchor)@lower_gradient(right))


def separate(anchor,pass4,pass3=None):
    first=anchor;result={field:"" for field in FIELDS}
    for field in (*IDENTITY,"Pt_bin","Pt_low","Pt_high","Pt_high_closed"):result[field]=first.get(field,"")
    result.update(Pass5_fit_status=anchor.get("Fit_status",""),Pass4_fit_status=pass4.get("Fit_status","") if pass4 else "MISSING",
                  Pass3_fit_status=pass3.get("Fit_status","") if pass3 else "MISSING",Separation_method="PASS4_PASS5_EXACT",
                  Separation_policy_version=POLICY_VERSION)
    if anchor.get("Fit_status") not in VALID_FITS:
        result.update(Separation_status="NOT_SEPARATED",Separation_reason="PASS5_C0_UNAVAILABLE",Pass3_validation_status="UNAVAILABLE");return result
    if not pass4 or pass4.get("Fit_status") not in VALID_FITS:
        result.update(Separation_status="NOT_SEPARATED",Separation_reason="PASS4_C0_UNAVAILABLE",Pass3_validation_status="UNAVAILABLE");return result
    e4,e5=number(pass4["epsilon"]),number(anchor["epsilon"]);c4,c5=number(pass4["C_0"]),number(anchor["C_0_pass5"])
    v4,v5=number(pass4["C_0_err"])**2,number(anchor["C_0_pass5_err"])**2;cov45=lower_anchor_c0_covariance(pass4,anchor)
    result.update(epsilon_pass4=e4,epsilon_pass5=e5,C_0_pass4=c4,C_0_pass4_err=math.sqrt(v4),
                  C_0_pass5=c5,C_0_pass5_err=math.sqrt(v5),Cov_C0_pass4_pass5=cov45)
    design=np.array([[1.,e4],[1.,e5]])
    if abs(np.linalg.det(design))<1e-12:
        result.update(Separation_status="ERROR",Separation_reason="INSUFFICIENT_EPSILON_LEVER_ARM",Pass3_validation_status="UNAVAILABLE");return result
    values=np.array([c4,c5]);covariance=np.array([[v4,cov45],[cov45,v5]])
    if np.linalg.matrix_rank(covariance)<2:
        result.update(Separation_status="ERROR",Separation_reason="SINGULAR_C0_COVARIANCE",Pass3_validation_status="UNAVAILABLE");return result
    inverse_design=np.linalg.inv(design);parameters=inverse_design@values;parameter_covariance=inverse_design@covariance@inverse_design.T
    ft,fl=parameters;ft_error=math.sqrt(max(0,parameter_covariance[0,0]));fl_error=math.sqrt(max(0,parameter_covariance[1,1]))
    correlation=parameter_covariance[0,1]/(ft_error*fl_error) if ft_error>0 and fl_error>0 else math.nan
    reasons=[]
    if ft<0:reasons.append("NEGATIVE_FT")
    if fl<0:reasons.append("NEGATIVE_FL")
    result.update(F_T=ft,F_T_err=ft_error,F_L=fl,F_L_err=fl_error,Cov_FT_FL=parameter_covariance[0,1],
                  Correlation_FT_FL=correlation,Separation_status="WARNING" if reasons else "OK",Separation_reason=";".join(reasons))
    if not pass3 or pass3.get("Fit_status") not in VALID_FITS:
        result["Pass3_validation_status"]="UNAVAILABLE";return result
    e3=number(pass3["epsilon"]);c3=number(pass3["C_0"]);v3=number(pass3["C_0_err"])**2
    cov34=lower_lower_covariance(pass3,pass4,anchor);cov35=lower_anchor_c0_covariance(pass3,anchor)
    cov3y=np.array([cov34,cov35]);cov3parameters=cov3y@inverse_design.T;projection=np.array([1.,e3])
    predicted=float(projection@parameters);prediction_variance=float(projection@parameter_covariance@projection)
    residual=c3-predicted;residual_variance=v3+prediction_variance-2*float(cov3parameters@projection)
    residual_error=math.sqrt(max(0,residual_variance));pull=residual/residual_error if residual_error>0 else math.nan
    validation="TENSION" if finite(pull) and abs(pull)>3 else "CONSISTENT" if finite(pull) else "ERROR"
    result.update(epsilon_pass3=e3,C_0_pass3=c3,C_0_pass3_err=math.sqrt(v3),Cov_C0_pass3_pass4=cov34,
                  Cov_C0_pass3_pass5=cov35,Pass3_predicted_C0=predicted,Pass3_prediction_err=math.sqrt(max(0,prediction_variance)),
                  Pass3_validation_residual=residual,Pass3_validation_err=residual_error,Pass3_validation_pull=pull,
                  Pass3_validation_status=validation)
    if validation=="TENSION":
        found=set(filter(None,result["Separation_reason"].split(";")));found.add("PASS3_VALIDATION_TENSION")
        result["Separation_reason"]=";".join(sorted(found));result["Separation_status"]="WARNING"
    return result


def draw_pdf(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True);canvas=ROOT.TCanvas("c_lt","lt",1800,1200);canvas.SetCanvasSize(1800,1200);canvas.Print(str(path)+"[")
    groups=defaultdict(list)
    for row in rows:groups[row["Epsilon_family_id"]].append(row)
    for family,members in sorted(groups.items()):
        canvas.Clear();canvas.Divide(3,2,.002,.002);owned=[];bypt={int(r["Pt_bin"]):r for r in members}
        for pt in range(1,6):
            row=bypt[pt];pad=canvas.cd(pt);pad.SetGridx();pad.SetGridy();pad.SetMargin(.15,.04,.15,.11)
            points=[]
            for pas,color,marker in (("pass3",ROOT.kOrange+7,24),("pass4",ROOT.kGreen+2,20),("pass5",ROOT.kBlue+1,20)):
                epsilon=number(row.get(f"epsilon_{pas}"));value=number(row.get(f"C_0_{pas}"));error=number(row.get(f"C_0_{pas}_err"))
                if finite(epsilon) and finite(value):points.append((pas,epsilon,value,error,color,marker))
            lows=[v-e for _,_,v,e,_,_ in points];highs=[v+e for _,_,v,e,_,_ in points];low=min(lows+[0]);high=max(highs+[1e-12]);span=max(high-low,1e-12)
            frame=ROOT.TH1D(f"ltf_{abs(hash((family,pt)))}","",50,.15,.9);frame.SetMinimum(min(0,low-.15*span));frame.SetMaximum(high+.2*span);frame.SetLineColorAlpha(ROOT.kWhite,0);frame.GetXaxis().SetTitle("#epsilon");frame.GetYaxis().SetTitle("C_{0}^{(p)} [GeV^{-2}]");frame.Draw()
            label=ROOT.TLatex();label.SetNDC();label.SetTextSize(.043);relation="#leq" if pt==5 else "<";label.DrawLatex(.20,.92,f"{PT_EDGES[pt-1]:.1f} #leq P_{{T}} {relation} {PT_EDGES[pt]:.1f} GeV")
            if row["Separation_status"] in {"OK","WARNING"}:
                line=ROOT.TF1(f"ltline_{abs(hash((family,pt)))}","[0]+[1]*x",.15,.9);line.SetParameters(number(row["F_T"]),number(row["F_L"]));line.SetLineColor(ROOT.kRed+1);line.SetLineWidth(2);line.Draw("SAME");owned.append(line)
            else:
                notice=ROOT.TLatex();notice.SetNDC();notice.SetTextAlign(22);notice.SetTextSize(.04);notice.SetTextColorAlpha(ROOT.kGray+2,.65);notice.DrawLatex(.55,.82,"Not separated");owned.append(notice)
            for pas,epsilon,value,error,color,marker in points:
                graph=ROOT.TGraphErrors();graph.SetPoint(0,epsilon,value);graph.SetPointError(0,0,error);graph.SetMarkerStyle(marker);graph.SetMarkerColor(color);graph.SetLineColor(color);graph.SetMarkerSize(1.15);graph.Draw("P SAME");owned.append(graph)
            owned.extend((frame,label))
        first=members[0];pad=canvas.cd(6);pad.Clear();text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.038)
        lines=["Pass-5 anchored L/T separation",f"{first['Run_type']} / {first['Target']}",f"x={number(first['x']):.3g}, Q^{{2}}={number(first['Q2']):.3g} GeV^{{2}}, z={number(first['z']):.3g}",f"Family: {family}","C_{0}^{(p)} = F_{T} + #epsilon_{p} F_{L}","Pass 4 + Pass 5: separation","Pass 3: validation only"]
        for i,line in enumerate(lines):text.DrawLatex(.05,.94-.068*i,line)
        y=.39
        for pas,color,marker,labeltext in (("pass3",ROOT.kOrange+7,24,"Pass 3 validation"),("pass4",ROOT.kGreen+2,20,"Pass 4 separation"),("pass5",ROOT.kBlue+1,20,"Pass 5 separation")):
            markerobj=ROOT.TMarker(.09,y,marker);markerobj.SetNDC(True);markerobj.SetMarkerColor(color);markerobj.Draw();text.DrawLatex(.14,y-.015,labeltext);owned.append(markerobj);y-=.07
        text.SetTextSize(.030);text.DrawLatex(.05,.08,"Red line: extracted F_{T} + #epsilon F_{L}");text.DrawLatex(.05,.035,POLICY_VERSION);owned.append(text);canvas._lt_owned=owned;canvas.Print(str(path))
    canvas.Print(str(path)+"]");canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pass5-fits",type=Path,default=project/"results/Pass1/Tables/RP_fit_sidis_reduced_xsec_pass5/RP_sidis_reduced_xsec_pass5_fit_parameters.csv")
    parser.add_argument("--lower-fits",type=Path,default=project/"results/Pass1/Tables/RP_fit_sidis_reduced_xsec_lower_passes/RP_sidis_reduced_xsec_lower_pass_fit_parameters.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_pass5_anchored_LT_separation")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_pass5_anchored_LT_separation");parser.add_argument("--no-pdf",action="store_true");args=parser.parse_args()
    anchors={(r["Epsilon_family_id"],int(r["Pt_bin"])):r for r in read_csv(args.pass5_fits)};lower={(r["Epsilon_family_id"],r["Member_pass"],int(r["Pt_bin"])):r for r in read_csv(args.lower_fits)}
    rows=[]
    for (family,pt),anchor in sorted(anchors.items()):rows.append(separate(anchor,lower.get((family,"pass4",pt)),lower.get((family,"pass3",pt))))
    atomic_csv(args.table_dir/"RP_sidis_pass5_anchored_LT_separation.csv",FIELDS,rows);atomic_csv(args.table_dir/"RP_sidis_pass5_anchored_LT_separation_Problematic.csv",FIELDS,[r for r in rows if r["Separation_status"]!="OK"])
    if not args.no_pdf:
        for z in sorted({r["z"] for r in rows},key=float):
            selected=[r for r in rows if r["z"]==z];zdir="z_"+f"{float(z):.2f}".replace(".","p");draw_pdf(args.pdf_dir/zdir/"RP_sidis_pass5_anchored_LT_separation.pdf",selected)
    print("Separations:",len(rows),Counter(r["Separation_status"] for r in rows));print("Pass3 validation:",Counter(r["Pass3_validation_status"] for r in rows));return 0


if __name__=="__main__":raise SystemExit(main())
