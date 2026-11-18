#!/usr/bin/env python3
"""Plot the separated virtual-photon ratio R=sigma_L/sigma_T."""

from __future__ import annotations

import argparse,csv,math,os,tempfile
from collections import defaultdict
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray+1);ROOT.gStyle.SetGridStyle(3);ROOT.TGaxis.SetMaxDigits(3)
POLICY_VERSION="RP_SIDIS_LT_RATIO_PLOTS_V1";VALID_STATUS={"OK","WARNING"}
X_COLORS={"0.25":ROOT.kBlue+1,"0.31":ROOT.kGreen+2,"0.44":ROOT.kOrange+7}
EXTRA_FIELDS=("Ratio_plot_policy_version","Pt_center","Pt2_low","Pt2_high","Pt2_center","Pt2_half_width")

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

def build_points(rows):
    output=[]
    for source in rows:
        if source.get("Separation_status") not in VALID_STATUS:continue
        if not finite(source.get("R_sigmaL_over_sigmaT")) or not finite(source.get("R_sigmaL_over_sigmaT_err")):continue
        low,high=number(source["Pt_low"]),number(source["Pt_high"]);row=dict(source)
        row.update(Ratio_plot_policy_version=POLICY_VERSION,Pt_center=.5*(low+high),Pt2_low=low*low,
                   Pt2_high=high*high,Pt2_center=.5*(low*low+high*high),Pt2_half_width=.5*(high*high-low*low));output.append(row)
    return output

def axis_range(values,errors=None,minimum_span=.2):
    expanded=[]
    for index,value in enumerate(values):
        error=errors[index] if errors else 0.;expanded.extend((value-error,value+error))
    if not expanded:return -1.,1.
    low,high=min(expanded),max(expanded);span=max(high-low,minimum_span);return low-.12*span,high+.16*span

def frame(name,x_title,x_low,x_high,y_low,y_high):
    hist=ROOT.TH1D(name,"",100,x_low,x_high);hist.SetMinimum(y_low);hist.SetMaximum(y_high);hist.SetLineColorAlpha(ROOT.kWhite,0)
    hist.GetXaxis().SetTitle(x_title);hist.GetYaxis().SetTitle("R = #sigma_{L}/#sigma_{T}");hist.Draw();return hist

def draw_graph(rows,x_field,x_error_field,color,style):
    graph=ROOT.TGraphErrors();expected="OK" if style==20 else "WARNING"
    for row in sorted(rows,key=lambda item:number(item[x_field])):
        if row["Separation_status"]!=expected:continue
        index=graph.GetN();graph.SetPoint(index,number(row[x_field]),number(row["R_sigmaL_over_sigmaT"]));graph.SetPointError(index,number(row.get(x_error_field),0) if x_error_field else 0,number(row["R_sigmaL_over_sigmaT_err"]))
    graph.SetMarkerStyle(style);graph.SetMarkerSize(1.15);graph.SetMarkerColor(color);graph.SetLineColor(color);return graph

def draw_panel(pad,rows,name,x_field,x_error_field,x_title,title,color_by_x=False,x_bounds=None):
    pad.SetGridx();pad.SetGridy();pad.SetMargin(.16,.04,.16,.12);owned=[]
    if not rows:
        xlow,xhigh=x_bounds or (0,1);hist=frame(name,x_title,xlow,xhigh,-1,1);hist.GetXaxis().SetNdivisions(505);notice=ROOT.TLatex();notice.SetNDC();notice.SetTextAlign(22);notice.SetTextColorAlpha(ROOT.kGray+2,.65);notice.DrawLatex(.55,.55,"No separated points");owned.extend((hist,notice));return owned
    xs=[number(r[x_field]) for r in rows];xe=[number(r.get(x_error_field),0) if x_error_field else 0 for r in rows];ys=[number(r["R_sigmaL_over_sigmaT"]) for r in rows];ye=[number(r["R_sigmaL_over_sigmaT_err"]) for r in rows]
    xlow,xhigh=x_bounds if x_bounds else axis_range(xs,xe,.05);ylow,yhigh=axis_range(ys,ye,.5);hist=frame(name,x_title,xlow,xhigh,ylow,yhigh);hist.GetXaxis().SetNdivisions(505);zero=ROOT.TLine(xlow,0,xhigh,0);zero.SetLineColor(ROOT.kViolet+1);zero.SetLineStyle(2);zero.Draw();owned.extend((hist,zero))
    groups=defaultdict(list)
    for row in rows:groups[row["x"] if color_by_x else "all"].append(row)
    for key,members in sorted(groups.items(),key=lambda item:number(item[0],0)):
        color=X_COLORS.get(key,ROOT.kBlue+1) if color_by_x else ROOT.kBlue+1
        for style in (20,24):
            graph=draw_graph(members,x_field,x_error_field,color,style)
            if graph.GetN():graph.Draw("P SAME")
            owned.append(graph)
    label=ROOT.TLatex();label.SetNDC();label.SetTextAlign(22);label.SetTextSize(.043);label.DrawLatex(.55,.92,title);owned.append(label);return owned

def info_pad(pad,heading,lines,show_x_legend=False):
    pad.Clear();owned=[];text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.037);text.DrawLatex(.05,.93,heading)
    for index,line in enumerate(lines):text.DrawLatex(.05,.85-.065*index,line)
    y=.42
    for style,label in ((20,"filled: separation OK"),(24,"open: separation warning")):
        marker=ROOT.TMarker(.09,y,style);marker.SetNDC(True);marker.SetMarkerColor(ROOT.kGray+2);marker.Draw();text.DrawLatex(.14,y-.014,label);owned.append(marker);y-=.07
    if show_x_legend:
        for x,color in sorted(X_COLORS.items(),key=lambda item:float(item[0])):
            marker=ROOT.TMarker(.09,y,20);marker.SetNDC(True);marker.SetMarkerColor(color);marker.Draw();text.DrawLatex(.14,y-.014,f"x = {x}");owned.append(marker);y-=.06
    text.SetTextSize(.029);text.DrawLatex(.05,.06,"Error bars include correlated L/T uncertainty");owned.append(text);return owned

def draw_five_panel_pdf(path,groups,mode):
    path.parent.mkdir(parents=True,exist_ok=True);canvas=ROOT.TCanvas(f"c_{mode}",mode,1800,1200);canvas.SetCanvasSize(1800,1200);canvas.Print(str(path)+"[")
    for identity,members in sorted(groups.items()):
        canvas.Clear();canvas.Divide(3,2,.002,.002);owned=[]
        if mode=="z":
            run,target,x,Q2=identity
            for pt in range(1,6):owned+=draw_panel(canvas.cd(pt),[r for r in members if int(r["Pt_bin"])==pt],f"rz_{abs(hash((identity,pt)))}","z","","z",f"P_{{T}} bin {pt}",x_bounds=(.30,.95))
            owned+=info_pad(canvas.cd(6),"R versus z",[f"{run} / {target}",f"x={number(x):.3g}, Q^{{2}}={number(Q2):.3g} GeV^{{2}}"])
        elif mode=="Q2":
            run,target,z=identity
            for pt in range(1,6):owned+=draw_panel(canvas.cd(pt),[r for r in members if int(r["Pt_bin"])==pt],f"rq_{abs(hash((identity,pt)))}","Q2","","Q^{2} [GeV^{2}]",f"P_{{T}} bin {pt}",True,(2.8,4.7))
            owned+=info_pad(canvas.cd(6),"R versus Q^{2}",[f"{run} / {target}",f"z={number(z):.3g}","x varies with available Q^{2} settings"],True)
        else:
            run,target,x,Q2=identity;zs=sorted({r["z"] for r in members},key=float)
            for index,z in enumerate(zs[:5],1):owned+=draw_panel(canvas.cd(index),[r for r in members if r["z"]==z],f"rp_{abs(hash((identity,z)))}","Pt2_center","Pt2_half_width","P_{T}^{2} [GeV^{2}]",f"z = {number(z):.3g}",x_bounds=(0,.25))
            for index in range(len(zs)+1,6):owned+=draw_panel(canvas.cd(index),[],f"rp_empty_{abs(hash((identity,index)))}","Pt2_center","Pt2_half_width","P_{T}^{2} [GeV^{2}]","",x_bounds=(0,.25))
            owned+=info_pad(canvas.cd(6),"R versus P_{T}^{2}",[f"{run} / {target}",f"x={number(x):.3g}, Q^{{2}}={number(Q2):.3g} GeV^{{2}}"])
        canvas._ratio_owned=owned;canvas.Print(str(path))
    canvas.Print(str(path)+"]");canvas.Close()

def main():
    project=Path(__file__).resolve().parents[1];parser=argparse.ArgumentParser(description=__doc__);source=project/"results/Pass1/Tables/RP_sidis_pass5_anchored_LT_separation/RP_sidis_pass5_anchored_virtual_photon_LT_separation.csv"
    parser.add_argument("--separations",type=Path,default=source);parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_LT_ratio_plots");parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_LT_ratio_plots");parser.add_argument("--no-pdf",action="store_true");args=parser.parse_args()
    source_rows=read_csv(args.separations);points=build_points(source_rows);fields=list(source_rows[0])+list(EXTRA_FIELDS);atomic_csv(args.table_dir/"RP_sidis_LT_ratio_points.csv",fields,points)
    if not args.no_pdf:
        zgroups=defaultdict(list);qgroups=defaultdict(list);ptgroups=defaultdict(list)
        for row in points:
            zgroups[(row["Run_type"],row["Target"],row["x"],row["Q2"])].append(row);qgroups[(row["Run_type"],row["Target"],row["z"])].append(row);ptgroups[(row["Run_type"],row["Target"],row["x"],row["Q2"])].append(row)
        draw_five_panel_pdf(args.pdf_dir/"RP_sidis_LT_ratio_vs_z.pdf",zgroups,"z");draw_five_panel_pdf(args.pdf_dir/"RP_sidis_LT_ratio_vs_Q2.pdf",qgroups,"Q2");draw_five_panel_pdf(args.pdf_dir/"RP_sidis_LT_ratio_vs_Pt2.pdf",ptgroups,"Pt2")
    print(f"ratio points={len(points)}");return 0

if __name__=="__main__":raise SystemExit(main())
