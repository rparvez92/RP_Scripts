#!/usr/bin/env python3
"""Tabulate and plot Pass-5 SIDIS harmonic coefficients versus Pt."""

from __future__ import annotations

import argparse, csv, math, os, tempfile
from collections import defaultdict
from pathlib import Path

import ROOT

import RP_data_mc_compare as dmc
import RP_sidis_azimuthal_analysis as azimuth

ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetGridColor(ROOT.kGray+1); ROOT.gStyle.SetGridStyle(3)

FAMILY_FIELDS=("Phase","Pass","Run_type","Target","x","Q2")
COMMON_FIELDS=[
    "Combination_group_id","Combination_policy_version","Ebeam_GeV","N_contributors",
    "Contributing_setting_ids","Contributing_thpq_deg",*FAMILY_FIELDS,"z",
    "Pt_bin","Pt_low","Pt_high","Pt_high_closed","Pt_center","Pt_half_width",
    "Coverage_binning_version","Coverage_policy_version","Coverage_curve_status",
    "Coverage_curve_reason","Coverage_curve_n_selected","Coverage_curve_min_selected_span_rad",
    "Coverage_design_rank","Coverage_weighted_design_condition",
]
OBS_FIELDS=[]
for prefix,coefficient in (("C","C"),("M","M")):
    OBS_FIELDS += [
        *(item for index in range(3) for item in (f"{coefficient}_{index}",f"{coefficient}_{index}_err")),
        f"{coefficient}_1_over_{coefficient}_0",f"{coefficient}_1_over_{coefficient}_0_err",
        f"{coefficient}_2_over_{coefficient}_0",f"{coefficient}_2_over_{coefficient}_0_err",
        f"{prefix}_Chi2",f"{prefix}_NDF",f"{prefix}_Chi2_NDF",f"{prefix}_Fit_probability",
        f"{prefix}_Fit_status",f"{prefix}_Fit_reason",f"{prefix}_Plotted_as_coefficient",
    ]
FIELDS=COMMON_FIELDS+OBS_FIELDS


def finite(value):
    try: return math.isfinite(float(value))
    except (TypeError,ValueError): return False


def number(value,default=math.nan): return float(value) if finite(value) else default


def ratio_with_covariance(numerator,numerator_error,denominator,denominator_error,covariance):
    values=(numerator,numerator_error,denominator,denominator_error,covariance)
    if not all(finite(value) for value in values) or float(denominator)==0: return math.nan,math.nan
    n,en,d,ed,cov=map(float,values)
    variance=en*en/(d*d)+n*n*ed*ed/(d**4)-2*n*cov/(d**3)
    if variance < -1e-12: return n/d,math.nan
    return n/d,math.sqrt(max(0.0,variance))


def plottable(row):
    try:
        return (row.get("Pass")=="pass5" and row.get("Coverage_curve_status")=="APPROVED"
                and int(row.get("N_points_fitted",0))==15
                and int(row.get("Coverage_curve_n_selected",0))==15
                and math.isclose(float(row.get("Coverage_curve_min_selected_span_rad")),2*math.pi,
                                 rel_tol=0,abs_tol=1e-12)
                and row.get("Fit_status") in {"OK","WARNING"})
    except (TypeError,ValueError): return False


def atomic_csv(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    handle,temporary=tempfile.mkstemp(prefix=path.name,suffix=".tmp",dir=path.parent)
    try:
        with os.fdopen(handle,"w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def build_rows(sigma_rows,multiplicity_rows):
    key=lambda row:(row["Combination_group_id"],int(row["Pt_bin"]))
    sigma={key(row):row for row in sigma_rows}; multiplicity={key(row):row for row in multiplicity_rows}
    if set(sigma)!=set(multiplicity): raise ValueError("sigma and multiplicity fit tables have different group/Pt keys")
    output=[]
    for item in sorted(sigma,key=lambda value:(value[0],value[1])):
        s,m=sigma[item],multiplicity[item]; row={field:"" for field in FIELDS}
        for field in COMMON_FIELDS:
            if field=="Pt_center": row[field]=(number(s["Pt_low"])+number(s["Pt_high"]))/2
            elif field=="Pt_half_width": row[field]=(number(s["Pt_high"])-number(s["Pt_low"]))/2
            else: row[field]=s.get(field,"")
        for source,prefix,coefficient in ((s,"C","C"),(m,"M","M")):
            for index in range(3):
                row[f"{coefficient}_{index}"]=source.get(f"{coefficient}_{index}","")
                row[f"{coefficient}_{index}_err"]=source.get(f"{coefficient}_{index}_err","")
            for index,covariance in ((1,"Cov_01"),(2,"Cov_02")):
                ratio,error=ratio_with_covariance(source.get(f"{coefficient}_{index}"),source.get(f"{coefficient}_{index}_err"),
                                                  source.get(f"{coefficient}_0"),source.get(f"{coefficient}_0_err"),source.get(covariance))
                row[f"{coefficient}_{index}_over_{coefficient}_0"]=ratio
                row[f"{coefficient}_{index}_over_{coefficient}_0_err"]=error
            for field in ("Chi2","NDF","Chi2_NDF","Fit_probability","Fit_status","Fit_reason"):
                row[f"{prefix}_{field}"]=source.get(field,"")
            row[f"{prefix}_Plotted_as_coefficient"]=int(plottable(source))
        output.append(row)
    return output


def safe_name(value): return str(value).replace(".","p").replace("-","m").replace("/","_")


def family_name(family,observable):
    phase,pas,run_type,target,x,q2=family
    return "_".join((phase,pas,run_type,target,"x"+safe_name(x),"Q2"+safe_name(q2),observable))+".pdf"


def displayed_z(row): return not math.isclose(number(row.get("z")),.90,rel_tol=0,abs_tol=1e-8)


def graph(rows,yfield,errorfield,status,color,style):
    result=ROOT.TGraphErrors()
    for row in sorted(rows,key=lambda item:(number(item["z"]),number(item["Pt_center"]))):
        if int(row[status+"_Plotted_as_coefficient"])!=1 or not finite(row.get(yfield)): continue
        expected="OK" if style==20 else "WARNING"
        if row[status+"_Fit_status"]!=expected: continue
        index=result.GetN(); result.SetPoint(index,number(row["Pt_center"]),number(row[yfield]))
        result.SetPointError(index,0,number(row.get(errorfield),0))
    result.SetMarkerStyle(style); result.SetMarkerSize(1.45); result.SetMarkerColor(color); result.SetLineColor(color)
    return result


def y_range(rows,fields,errors=None,include_zero=True):
    values=[]
    for row in rows:
        for index,field in enumerate(fields):
            if not finite(row.get(field)): continue
            error=number(row.get(errors[index]),0) if errors else 0
            values += [number(row[field])-error,number(row[field])+error]
    if include_zero: values.append(0.0)
    if not values: return -1.0,1.0
    low,high=min(values),max(values); span=max(high-low,abs(high)*.1,1e-12)
    return low-.18*span,high+.20*span


def base_frame(name,title,ytitle,low,high):
    frame=ROOT.TH1D(name,f";P_{{t}} [GeV];{ytitle}",5,0,.5)
    frame.SetMinimum(low); frame.SetMaximum(high); frame.GetXaxis().SetNdivisions(505)
    for axis in (frame.GetXaxis(),frame.GetYaxis()):
        axis.SetLabelSize(.046); axis.SetTitleSize(.054)
    frame.GetXaxis().SetTitleOffset(1.15); frame.GetYaxis().SetTitleOffset(1.65)
    frame.SetLineColorAlpha(ROOT.kWhite,0); frame.Draw()
    label=ROOT.TLatex(); label.SetNDC(); label.SetTextFont(62); label.SetTextAlign(22); label.SetTextSize(.064)
    drawn=label.DrawLatex(.52,.875,title)
    return frame,label,drawn


def draw_series(rows,yfield,errorfield,status,low,high,owned,show_not_fitted=False):
    byz=defaultdict(list)
    for row in rows: byz[str(row["z"])].append(row)
    for z,members in sorted(byz.items(),key=lambda item:float(item[0])):
        color=azimuth.z_color(z)
        for style in (20,24):
            item=graph(members,yfield,errorfield,status,color,style)
            if item.GetN(): item.Draw("P SAME")
            owned.append(item)
        if show_not_fitted:
            missing=ROOT.TGraph(); missing.SetMarkerStyle(5); missing.SetMarkerSize(1.55); missing.SetMarkerColor(color)
            y=low+.045*(high-low)
            for row in members:
                if int(row[status+"_Plotted_as_coefficient"])==0: missing.SetPoint(missing.GetN(),number(row["Pt_center"]),y)
            if missing.GetN(): missing.Draw("P SAME")
            owned.append(missing)


def page_header(canvas,family,subtitle):
    phase,pas,run_type,target,x,q2=family
    text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.033)
    first=text.DrawLatex(.045,.955,subtitle)
    second=text.DrawLatex(.045,.918,f"{phase}, {pas}, {run_type}, {target}; x={x}, Q^{{2}}={q2} GeV^{{2}}; full 2#pi coverage (15/15 bins)")
    return text,first,second


def page_legend(canvas,rows,status,include_not_fitted=True):
    canvas.cd(); owned=[]
    text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42); text.SetTextAlign(12); text.SetTextSize(.025)
    z_values=sorted({str(row["z"]) for row in rows},key=float)
    spacing=.19; start=.5-spacing*(len(z_values)-1)/2
    for index,z in enumerate(z_values):
        x=start+index*spacing; marker=ROOT.TMarker(x,.112,20); marker.SetNDC(True)
        marker.SetMarkerSize(1.45); marker.SetMarkerColor(azimuth.z_color(z)); marker.Draw()
        drawn_text=text.DrawLatex(x+.014,.112,f"z = {float(z):g}")
        owned += [marker,drawn_text]
    status_entries=[(.20,20,"filled: fit OK"),(.43,24,"open: fit warning")]
    if include_not_fitted: status_entries.append((.66,5,"colored #times: coverage rejected"))
    for x,style,label in status_entries:
        marker=ROOT.TMarker(x,.067,style); marker.SetNDC(True)
        marker.SetMarkerSize(1.4); marker.SetMarkerColor(ROOT.kGray+2); marker.Draw()
        drawn_text=text.DrawLatex(x+.014,.067,label); owned += [marker,drawn_text]
    owned.append(text); return owned


def draw_absolute(canvas,path,rows,family,observable):
    coefficient,status=("C","C") if observable=="sigma" else ("M","M")
    canvas.Clear(); canvas.Divide(3,1,.002,.002); owned=[]
    for index in range(3):
        field,error=f"{coefficient}_{index}",f"{coefficient}_{index}_err"
        pad=canvas.cd(index+1); pad.SetGrid(); pad.SetMargin(.19,.045,.255,.15)
        eligible=[row for row in rows if int(row[status+"_Plotted_as_coefficient"])==1]
        low,high=y_range(eligible,[field],[error])
        frame,label,drawn=base_frame(f"abs_{observable}_{index}_{abs(hash(family))}",f"{coefficient}_{{{index}}}",f"{coefficient}_{{{index}}}",low,high)
        zero=ROOT.TLine(0,0,.5,0); zero.SetLineStyle(2); zero.SetLineColor(ROOT.kViolet+1); zero.Draw()
        draw_series(rows,field,error,status,low,high,owned,True); owned += [frame,label,drawn,zero]
    canvas.cd(); owned += list(page_header(canvas,family,"Absolute harmonic coefficients"))
    owned += page_legend(canvas,rows,status,True); canvas._owned=owned; canvas.Print(str(path))


def draw_normalized(canvas,path,rows,family,observable):
    coefficient,status=("C","C") if observable=="sigma" else ("M","M")
    canvas.Clear(); canvas.Divide(2,1,.002,.002); owned=[]
    for pad_index,index in enumerate((1,2),1):
        pad=canvas.cd(pad_index); pad.SetGrid(); pad.SetMargin(.17,.045,.255,.15)
        field=f"{coefficient}_{index}_over_{coefficient}_0"; error=field+"_err"
        eligible=[row for row in rows if int(row[status+"_Plotted_as_coefficient"])==1]
        low,high=y_range(eligible,[field],[error]); title=f"{coefficient}_{{{index}}}/{coefficient}_{{0}}"
        frame,label,drawn=base_frame(f"norm_{observable}_{index}_{abs(hash(family))}",title,title,low,high)
        zero=ROOT.TLine(0,0,.5,0); zero.SetLineStyle(2); zero.SetLineColor(ROOT.kViolet+1); zero.Draw()
        draw_series(rows,field,error,status,low,high,owned,True); owned += [frame,label,drawn,zero]
    canvas.cd(); owned += list(page_header(canvas,family,"Normalized harmonic coefficients (covariance propagated)"))
    owned += page_legend(canvas,rows,status,True); canvas._owned=owned; canvas.Print(str(path))


def draw_diagnostics(canvas,path,rows,family,observable):
    status="C" if observable=="sigma" else "M"
    specs=((f"{status}_Chi2_NDF","#chi^{2}/NDF"),(f"{status}_Fit_probability","Fit probability"),("Coverage_weighted_design_condition","Coverage weighted condition"))
    canvas.Clear(); canvas.Divide(3,1,.002,.002); owned=[]
    for index,(field,title) in enumerate(specs,1):
        pad=canvas.cd(index); pad.SetGrid(); pad.SetMargin(.19,.045,.255,.15)
        eligible=[row for row in rows if int(row[status+"_Plotted_as_coefficient"])==1 and finite(row.get(field))]
        low,high=y_range(eligible,[field]); frame,label,drawn=base_frame(f"diag_{observable}_{index}_{abs(hash(family))}",title,title,low,high)
        draw_series(rows,field,"",status,low,high,owned,False); owned += [frame,label,drawn]
    canvas.cd(); owned += list(page_header(canvas,family,"Fit and coverage diagnostics"))
    owned += page_legend(canvas,rows,status,False); canvas._owned=owned; canvas.Print(str(path))


def draw_pdf(path,rows,family,observable):
    path.parent.mkdir(parents=True,exist_ok=True)
    canvas=ROOT.TCanvas("coefficient_canvas","coefficients",2200,1700); canvas.SetCanvasSize(2200,1700)
    canvas.Print(str(path)+"[")
    draw_absolute(canvas,path,rows,family,observable); draw_normalized(canvas,path,rows,family,observable); draw_diagnostics(canvas,path,rows,family,observable)
    canvas.Print(str(path)+"]"); canvas.Close()


def main():
    project=Path(__file__).resolve().parents[1]; parser=argparse.ArgumentParser(description=__doc__)
    source=project/"results/Pass1/Tables/RP_sidis_azimuthal_analysis"
    parser.add_argument("--sigma-fits",type=Path,default=source/"RP_sidis_sigma_vs_phipq_FitParameters.csv")
    parser.add_argument("--multiplicity-fits",type=Path,default=source/"RP_sidis_multiplicity_vs_phipq_FitParameters.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_coefficient_plots")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_coefficient_plots")
    parser.add_argument("--no-pdf",action="store_true"); args=parser.parse_args()
    rows=build_rows(dmc.read_csv(args.sigma_fits),dmc.read_csv(args.multiplicity_fits))
    atomic_csv(args.table_dir/"RP_sidis_coefficients.csv",FIELDS,rows)
    problematic=[row for row in rows if row["C_Fit_status"]!="OK" or row["M_Fit_status"]!="OK"]
    atomic_csv(args.table_dir/"RP_sidis_coefficients_Problematic.csv",FIELDS,problematic)
    if not args.no_pdf:
        families=defaultdict(list)
        for row in rows:
            if displayed_z(row): families[tuple(str(row[field]) for field in FAMILY_FIELDS)].append(row)
        for family,members in sorted(families.items()):
            for observable in ("sigma","multiplicity"):
                output=args.pdf_dir/family_name(family,observable); draw_pdf(output,members,family,observable); print("OK",output.name)
    print(f"rows={len(rows)} problematic={len(problematic)}"); return 0


if __name__=="__main__": raise SystemExit(main())
