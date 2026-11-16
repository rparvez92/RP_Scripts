#!/usr/bin/env python3
"""Plot global reduced-SIDIS coefficients and audit anchored/global agreement."""

from __future__ import annotations

import argparse,csv,math,os,tempfile
from collections import Counter,defaultdict
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray+1);ROOT.gStyle.SetGridStyle(3);ROOT.TGaxis.SetMaxDigits(3)

POLICY_VERSION="RP_REDUCED_STRUCTURE_FUNCTION_PLOTS_V1"
VALID={"OK","WARNING"}
Z_COLORS={.36:ROOT.kRed+1,.50:ROOT.kGreen+2,.52:ROOT.kGreen-7,.67:ROOT.kBlue+1,.90:ROOT.kBlack}
IDENTITY=("Epsilon_family_id","Run_type","Target","x","Q2","z","Pt_bin","Pt_low","Pt_high","Pt_high_closed")
FIELDS=[*IDENTITY,"Pt_center","Global_fit_status","Global_fit_reason","C_0_global","C_0_global_err","C_1_global","C_1_global_err","C_2_global","C_2_global_err","C_3_global","C_3_global_err","Global_F_T","Global_F_T_err","Global_F_L","Global_F_L_err","Anchored_status","Anchored_reason","Anchored_F_T","Anchored_F_T_err","Anchored_F_L","Anchored_F_L_err","Delta_F_T_global_minus_anchored","Delta_F_L_global_minus_anchored","Comparison_status","Comparison_note","Plot_policy_version"]


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


def comparison_rows(global_rows,anchored_rows):
    anchored={(r["Epsilon_family_id"],int(r["Pt_bin"])):r for r in anchored_rows};output=[]
    for global_fit in global_rows:
        key=(global_fit["Epsilon_family_id"],int(global_fit["Pt_bin"]));anchor=anchored.get(key,{})
        global_ok=global_fit.get("Fit_status") in VALID;anchor_ok=anchor.get("Separation_status") in VALID
        row={field:global_fit.get(field,"") for field in IDENTITY};row.update(Pt_center=.5*(number(global_fit["Pt_low"])+number(global_fit["Pt_high"])),
            Global_fit_status=global_fit.get("Fit_status",""),Global_fit_reason=global_fit.get("Fit_reason",""),Anchored_status=anchor.get("Separation_status","MISSING"),Anchored_reason=anchor.get("Separation_reason",""),Plot_policy_version=POLICY_VERSION)
        for i in range(4):row[f"C_{i}_global"]=global_fit.get(f"C_{i}_global","");row[f"C_{i}_global_err"]=global_fit.get(f"C_{i}_global_err","")
        row.update(Global_F_T=global_fit.get("F_T",""),Global_F_T_err=global_fit.get("F_T_err",""),Global_F_L=global_fit.get("F_L",""),Global_F_L_err=global_fit.get("F_L_err",""),Anchored_F_T=anchor.get("F_T",""),Anchored_F_T_err=anchor.get("F_T_err",""),Anchored_F_L=anchor.get("F_L",""),Anchored_F_L_err=anchor.get("F_L_err",""))
        if global_ok and anchor_ok:
            row.update(Delta_F_T_global_minus_anchored=number(global_fit["F_T"])-number(anchor["F_T"]),Delta_F_L_global_minus_anchored=number(global_fit["F_L"])-number(anchor["F_L"]),Comparison_status="AVAILABLE",Comparison_note="SHARED_DATA_CROSS_METHOD_COVARIANCE_NOT_AVAILABLE_NO_PULL")
        else:row.update(Comparison_status="UNAVAILABLE",Comparison_note="ONE_OR_BOTH_METHODS_NOT_AVAILABLE")
        output.append(row)
    return output


def z_color(value):return next((color for z,color in Z_COLORS.items() if abs(number(value)-z)<1e-6),ROOT.kMagenta+1)


def graph(rows,value_field,error_field,color,status):
    result=ROOT.TGraphErrors();i=0
    for row in sorted(rows,key=lambda r:number(r["Pt_center"])):
        if row["Global_fit_status"] != status or not finite(row[value_field]):continue
        result.SetPoint(i,number(row["Pt_center"]),number(row[value_field]));result.SetPointError(i,0,number(row[error_field],0));i+=1
    result.SetMarkerColor(color);result.SetLineColor(color);result.SetMarkerStyle(20 if status == "OK" else 24);return result


def safe_name(key):
    run,target,x,Q2=key;return f"{run}_{target}_x{number(x):.2f}_Q2{number(Q2):.2f}".replace(".","p")


def draw_pdf(path,key,rows):
    path.parent.mkdir(parents=True,exist_ok=True);canvas=ROOT.TCanvas("c_structure","structure",1800,1200);canvas.SetCanvasSize(1800,1200);canvas.Divide(3,2,.002,.002);owned=[]
    panels=(("C_0_global","C_0_global_err","C_{0}=F_{T}"),("C_1_global","C_1_global_err","C_{1}=F_{L}"),("C_2_global","C_2_global_err","C_{2}=F_{UU}^{cos#phi}"),("C_3_global","C_3_global_err","C_{3}=F_{UU}^{cos2#phi}"))
    byz=defaultdict(list)
    for row in rows:byz[row["z"]].append(row)
    for index,(value,error,title) in enumerate(panels,1):
        pad=canvas.cd(index);pad.SetGridx();pad.SetGridy();pad.SetMargin(.15,.04,.15,.11)
        valid=[r for r in rows if r["Global_fit_status"] in VALID and finite(r[value])];lows=[number(r[value])-number(r[error],0) for r in valid];highs=[number(r[value])+number(r[error],0) for r in valid];low=min(lows+[0]);high=max(highs+[1e-12]);span=max(high-low,1e-12)
        frame=ROOT.TH1D(f"sf_{abs(hash((key,index)))}","",50,0,.5);frame.SetMinimum(min(0,low-.15*span));frame.SetMaximum(high+.2*span);frame.SetLineColorAlpha(ROOT.kWhite,0);frame.GetXaxis().SetTitle("P_{T} [GeV]");frame.GetYaxis().SetTitle(title+" [GeV^{-2}]");frame.Draw();zero=ROOT.TLine(0,0,.5,0);zero.SetLineColor(ROOT.kViolet+1);zero.SetLineStyle(2);zero.Draw();owned.extend((frame,zero))
        for z,zrows in sorted(byz.items(),key=lambda item:number(item[0])):
            for status in ("OK","WARNING"):
                g=graph(zrows,value,error,z_color(z),status);g.SetMarkerSize(1.1)
                if g.GetN():g.Draw("P SAME")
                owned.append(g)
    pad=canvas.cd(5);pad.Clear();legend=ROOT.TLegend(.12,.20,.88,.88);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.055)
    for z in sorted(byz,key=number):
        marker=ROOT.TGraph();marker.SetMarkerStyle(20);marker.SetMarkerColor(z_color(z));legend.AddEntry(marker,f"z = {number(z):g}","p");owned.append(marker)
    legend.Draw();owned.append(legend)
    pad=canvas.cd(6);pad.Clear();text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.043);run,target,x,Q2=key
    lines=["Global reduced-SIDIS coefficients",f"{run} / {target}",f"x={number(x):.3g}, Q^{{2}}={number(Q2):.3g} GeV^{{2}}","Filled: OK fit","Open: fit warning","Pass-3/4/5 global extraction"]
    for i,line in enumerate(lines):text.DrawLatex(.08,.90-.10*i,line)
    text.SetTextSize(.032);text.DrawLatex(.08,.12,"Anchored/global comparison is stored in the CSV audit");text.DrawLatex(.08,.06,POLICY_VERSION);owned.append(text);canvas._structure_owned=owned;canvas.Print(str(path));canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--global-fits",type=Path,default=project/"results/Pass1/Tables/RP_sidis_global_reduced_xsec_fit/RP_sidis_global_reduced_xsec_fit_parameters.csv")
    parser.add_argument("--anchored",type=Path,default=project/"results/Pass1/Tables/RP_sidis_pass5_anchored_LT_separation/RP_sidis_pass5_anchored_LT_separation.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_reduced_structure_function_plots")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_reduced_structure_function_plots");parser.add_argument("--no-pdf",action="store_true");args=parser.parse_args()
    rows=comparison_rows(read_csv(args.global_fits),read_csv(args.anchored));atomic_csv(args.table_dir/"RP_sidis_reduced_structure_functions.csv",FIELDS,rows);atomic_csv(args.table_dir/"RP_sidis_reduced_structure_functions_Problematic.csv",FIELDS,[r for r in rows if r["Global_fit_status"]!="OK" or r["Comparison_status"]!="AVAILABLE"])
    if not args.no_pdf:
        groups=defaultdict(list)
        for row in rows:groups[(row["Run_type"],row["Target"],row["x"],row["Q2"])].append(row)
        for key,members in sorted(groups.items()):draw_pdf(args.pdf_dir/f"RP_sidis_reduced_structure_functions_{safe_name(key)}.pdf",key,members)
    print("Rows:",len(rows),"Global:",Counter(r["Global_fit_status"] for r in rows),"Comparison:",Counter(r["Comparison_status"] for r in rows));return 0


if __name__=="__main__":raise SystemExit(main())
