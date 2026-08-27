#!/usr/bin/env python3
"""Extract Delta-acceptance one-dimensional SIDIS projections."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import ROOT

import RP_data_mc_compare as dmc
import RP_sidis_differential_common as common

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

FIELDS = [
    *common.IDENTITY_FIELDS, "Acceptance", "Variable", "Bin_index", "Bin_low", "Bin_high", "Bin_center",
    "x_eff", "x_RMS", "Q2_eff", "Q2_RMS", "z_eff", "z_RMS", "pT_eff", "pT_RMS",
    "Y_Data", "Y_Data_err", *(item for r in common.REACTIONS for item in (f"Y_MC_{r}", f"Y_MC_{r}_err")),
    "Y_MC", "Y_MC_err", "C_Y", "C_Y_err", "sigma_sidis_model", "sigma_sidis_data",
    "sigma_sidis_data_err", "sigma_sidis_units", "M_sidis_model", "M_sidis_event_weighted",
    "M_sidis_model_consistency_rel", "M_sidis_data", "M_sidis_data_err", "M_sidis_units",
    "Purity", "Stability", "Recon_minus_original_mean", "Recon_minus_original_RMS",
    "Y_Data_rel_err", "Y_MC_rel_err", "Data_underflow", "Data_overflow", "MC_underflow", "MC_overflow",
    "Y_Data_integral_closure_rel", *(f"Y_MC_{r}_QA_closure_rel" for r in common.REACTIONS),
    "Y_MC_complete", "Projection_status", "Projection_reason",
]


def aggregate(hist, start: int, stop: int) -> tuple[float, float]:
    return common.hist_range(hist, start, stop)


def fine_span(axis, low: float, high: float) -> tuple[int, int]:
    eps = 1e-10 * max(1.0, abs(high - low))
    return axis.FindFixBin(low + eps), axis.FindFixBin(high - eps)


def ratio_error(data, data_err, mc, mc_err):
    if not all(common.finite(item) for item in (data, data_err, mc, mc_err)) or mc == 0:
        return math.nan, math.nan
    return data / mc, math.hypot(data_err / mc, data * mc_err / (mc * mc))


def relative(observed, reference):
    return (observed - reference) / reference if common.finite(observed) and common.finite(reference) and reference != 0 else math.nan


def status_row(identity, status, reason):
    row = {field: "" for field in FIELDS}; row.update(identity)
    row.update({"Acceptance": "Delta-only", "Projection_status": status, "Projection_reason": reason})
    return row


def upstream_integrals(path: Path):
    if not path.is_file(): return {}
    with path.open(newline="") as stream:
        rows=[row for row in csv.DictReader(stream) if row.get("Tier")=="delta" and row.get("Variable")=="hdelta" and row.get("Bin_index")]
    result={}
    for output,input_name in (("Data","Data_final"),*((r,f"MC_{r}") for r in common.REACTIONS)):
        values=[common.number(row.get(input_name)) for row in rows]
        result[output]=sum(values) if values and all(common.finite(v) for v in values) else math.nan
    return result


def graph(rows, xfield, yfield, efield, color, marker):
    valid = [row for row in rows if common.finite(row.get(xfield)) and common.finite(row.get(yfield))]
    result = ROOT.TGraphErrors(len(valid))
    for index, row in enumerate(valid):
        result.SetPoint(index, float(row[xfield]), float(row[yfield]))
        result.SetPointError(index, 0.0, common.number(row.get(efield), 0.0))
    result.SetMarkerStyle(marker); result.SetMarkerColor(color); result.SetLineColor(color)
    return result


def draw_metadata(canvas, pdf, identity, status, reasons, config):
    canvas.Clear(); text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    lines = ["One-dimensional SIDIS extraction", ", ".join(identity.values()),
             "Projection over the other accepted variables; Delta acceptance only.",
             "Y_{MC}=Y_{sidis}+Y_{rho}+Y_{delta}+Y_{exclusive}; C_{Y}=Y_{Data}/Y_{MC}",
             "#sigma_{sidis,data}=C_{Y}#sigma_{sidis,model}; M_{sidis,data}=C_{Y}M_{sidis,model}",
             "Model kinematics are fWeight-weighted reconstructed SIDIS means in each bin.",
             "M_{sidis,model} is calc_semi_xsec sighad; event-level sighad is an independent check.",
             "Migration purity/stability use |fWeight|. No Neff warning is propagated.",
             f"Binning source: {config}", f"Status: {status}", f"Reason: {';'.join(reasons) if reasons else 'none'}"]
    y = .93
    for i, line in enumerate(lines):
        text.SetTextSize(.040 if i == 0 else .028); text.DrawLatex(.055, y, line); y -= .069
    canvas.Print(str(pdf) + "(")


def draw_pages(canvas, pdf, variable, rows, response):
    usable = [row for row in rows if row.get("Variable") == variable and row.get("Bin_index")]
    if not usable: return
    title = common.FINE_BINNING[variable][3]
    for yfield, efield, ylabel, model_field in (
        ("sigma_sidis_data", "sigma_sidis_data_err", "#sigma_{sidis} (#mub/GeV^{2}/sr^{2})", "sigma_sidis_model"),
        ("M_sidis_data", "M_sidis_data_err", "M_{sidis} (GeV^{-2})", "M_sidis_model"),
    ):
        canvas.Clear(); frame = ROOT.TH1D(f"f_{variable}_{yfield}_{id(rows)}", f"{ylabel} vs {title};{title};{ylabel}", 100, usable[0]["Bin_low"], usable[-1]["Bin_high"]); frame.SetStats(False)
        values = [float(r[yfield]) for r in usable if common.finite(r.get(yfield))]
        models = [float(r[model_field]) for r in usable if common.finite(r.get(model_field))]
        ymax = max(values + models + [1e-12]); ymin = min(values + models + [0.0]); frame.SetMinimum(min(0., ymin)*1.15); frame.SetMaximum(ymax*1.25); frame.Draw()
        gd = graph(usable, "Bin_center", yfield, efield, ROOT.kBlack, 20); gm = graph(usable, "Bin_center", model_field, "", ROOT.kBlue+1, 24)
        gd.Draw("P SAME"); gm.Draw("P SAME")
        leg=ROOT.TLegend(.72,.78,.94,.92); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.AddEntry(gd,"Extracted","p"); leg.AddEntry(gm,"Model","p"); leg.Draw(); canvas.SetGrid(); canvas.Print(str(pdf))
    canvas.Clear(); frame=ROOT.TH1D(f"fy_{variable}_{id(rows)}",f"Normalized Yield projection: {title};{title};Normalized Yield",100,usable[0]["Bin_low"],usable[-1]["Bin_high"]); frame.SetStats(False)
    ymax=max([common.number(r.get("Y_Data"),0) for r in usable]+[common.number(r.get("Y_MC"),0) for r in usable]+[1e-12]); frame.SetMinimum(0); frame.SetMaximum(1.35*ymax); frame.Draw()
    colors={"sidis":ROOT.kBlue+1,"rho":ROOT.kGreen+2,"delta":ROOT.kMagenta+1,"exclusive":ROOT.kOrange+7}
    gd=graph(usable,"Bin_center","Y_Data","Y_Data_err",ROOT.kBlack,20); gd.Draw("P SAME"); leg=ROOT.TLegend(.70,.64,.94,.92); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.AddEntry(gd,"Data","p"); owned_graphs=[gd]
    for reaction in common.REACTIONS:
        g=graph(usable,"Bin_center",f"Y_MC_{reaction}",f"Y_MC_{reaction}_err",colors[reaction],24); g.Draw("P SAME"); leg.AddEntry(g,reaction,"p"); owned_graphs.append(g)
    gm=graph(usable,"Bin_center","Y_MC","Y_MC_err",ROOT.kRed+1,21); gm.Draw("P SAME"); leg.AddEntry(gm,"MC","p"); owned_graphs.append(gm); leg.Draw(); canvas.SetGrid(); canvas.Print(str(pdf))
    canvas.Clear(); frame=ROOT.TH1D(f"fc_{variable}_{id(rows)}",f"C_{{Y}} projection: {title};{title};C_{{Y}}",100,usable[0]["Bin_low"],usable[-1]["Bin_high"]); frame.SetStats(False); frame.Draw(); gc=graph(usable,"Bin_center","C_Y","C_Y_err",ROOT.kBlack,20); gc.Draw("P SAME"); canvas.SetGrid(); canvas.Print(str(pdf))
    canvas.Clear(); canvas.SetRightMargin(.14); response.SetTitle(f"SIDIS original vs reconstructed: {title};Original;Reconstructed"); response.Draw("COLZ"); canvas.Print(str(pdf)); canvas.SetRightMargin(.04)
    canvas.Clear(); frame=ROOT.TH1D(f"fq_{variable}_{id(rows)}",f"Migration and statistical precision: {title};{title};Metric",100,usable[0]["Bin_low"],usable[-1]["Bin_high"]); frame.SetStats(False); frame.SetMinimum(0); frame.SetMaximum(1.2); frame.Draw()
    leg=ROOT.TLegend(.70,.70,.94,.92); leg.SetBorderSize(0); leg.SetFillStyle(0); metric_graphs=[]
    for field,color,marker in (("Purity",ROOT.kBlue+1,20),("Stability",ROOT.kGreen+2,21),("Y_Data_rel_err",ROOT.kBlack,22),("Y_MC_rel_err",ROOT.kRed+1,23)):
        g=graph(usable,"Bin_center",field,"",color,marker); g.Draw("PL SAME"); leg.AddEntry(g,field.replace("Y_","").replace("_rel_err"," rel. error"),"lp"); metric_graphs.append(g)
    leg.Draw(); canvas.SetGrid(); canvas.Print(str(pdf))


def draw_closure_page(canvas, pdf, rows, final_status):
    canvas.Clear(); text=ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    lines=["Integrated closure and final status"]
    for variable in common.VARIABLES:
        selected=[row for row in rows if row.get("Variable")==variable]
        if not selected: continue
        first=selected[0]
        lines.append(f"{variable}: Data/input = {common.number(first.get('Y_Data_integral_closure_rel')):.4g}")
        lines.append("  MC/QA: "+", ".join(f"{r}={common.number(first.get(f'Y_MC_{r}_QA_closure_rel')):.4g}" for r in common.REACTIONS))
    reasons=sorted({item for row in rows for item in str(row.get("Projection_reason","")).split(";") if item})
    lines.extend((f"Final status: {final_status}",f"Reasons: {';'.join(reasons) if reasons else 'none'}"))
    y=.92
    for i,line in enumerate(lines): text.SetTextSize(.042 if i==0 else .031); text.DrawLatex(.06,y,line); y-=.065
    canvas.Print(str(pdf))


def draw_integrated_overview(path: Path, integrated_path: Path, models):
    """Plot existing setting-integrated results; central phi is intentionally absent."""
    if not integrated_path.is_file():
        return False
    with integrated_path.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row.get("Extraction_status") in {"OK", "WARNING"}]
    canvas=ROOT.TCanvas("c_sidis_integrated_overview","overview",1800,1100); first=True
    categories=(("PIPLUS","LH2",ROOT.kBlue+1,20),("PIMINUS","LH2",ROOT.kRed+1,21),("PIPLUS","LD2",ROOT.kGreen+2,22),("PIMINUS","LD2",ROOT.kMagenta+1,23))
    for xfield,xtitle in (("z","Central z"),("pT_GeV","Central p_{T} (GeV)")):
        for yfield,ytitle in (("sigma_sidis_data","#sigma_{sidis,data} (#mub/GeV^{2}/sr^{2})"),("M_sidis_data","M_{sidis,data} (GeV^{-2})")):
            canvas.Clear(); frame=ROOT.TH1D(f"ov_{xfield}_{yfield}",f"Setting-integrated overview;{xtitle};{ytitle}",100,0,1); frame.SetStats(False)
            points=[]
            for run_type,target,color,marker in categories:
                selected=[]
                for row in rows:
                    if row.get("Run_type")!=run_type or row.get("Target")!=target: continue
                    key=tuple(row.get(field,"") for field in common.IDENTITY_FIELDS); model=models.get(key,{})
                    x=common.number(model.get(xfield)); y=common.number(row.get(yfield)); e=common.number(row.get(yfield+"_err"),0)
                    if common.finite(x) and common.finite(y): selected.append((x,y,e))
                graph=ROOT.TGraphErrors(len(selected))
                for i,(x,y,e) in enumerate(selected): graph.SetPoint(i,x,y); graph.SetPointError(i,0,e)
                graph.SetMarkerStyle(marker); graph.SetMarkerColor(color); graph.SetLineColor(color); points.append((f"{run_type} {target}",graph))
            all_y=[g.GetPointY(i) for _,g in points for i in range(g.GetN())]
            frame.SetMaximum(1.25*max(all_y+[1e-12])); frame.SetMinimum(min(0.,min(all_y+[0.]))*1.15); frame.Draw()
            legend=ROOT.TLegend(.72,.70,.94,.92); legend.SetBorderSize(0); legend.SetFillStyle(0)
            for label,g in points: g.Draw("P SAME"); legend.AddEntry(g,label,"p")
            legend.Draw(); canvas.SetGrid(); canvas.Print(str(path)+("(" if first else "")); first=False
    if not first: canvas.Print(str(path)+")")
    return not first


def process(key, files, catalog, models, config, calculator, t7_root, data_mc_dir, table_dir, pdf_dir):
    identity=common.identity_from_key(key); stem=common.setting_stem(identity)
    configured=config.get("settings",{}).get(stem)
    if configured is None: rows=[status_row(identity,"PENDING","RECOMMENDED_BINNING_MISSING")]; common.atomic_csv(table_dir/f"{stem}.csv",FIELDS,rows); return rows,"PENDING"
    if configured.get("status") in {"SKIPPED","PENDING","ERROR"}: status=configured["status"]; rows=[status_row(identity,status,configured.get("reason",""))]; common.atomic_csv(table_dir/f"{stem}.csv",FIELDS,rows); return rows,status
    payload,status,reasons=common.prepare_histograms(key,files,catalog,common.FINE_BINNING,t7_root)
    if payload is None: rows=[status_row(identity,status,";".join(reasons))]; common.atomic_csv(table_dir/f"{stem}.csv",FIELDS,rows); return rows,status
    model=models.get(tuple(key));
    if not model or model.get("Model_status")!="OK": rows=[status_row(identity,"PENDING","MODEL_ROW_UNAVAILABLE")]; common.atomic_csv(table_dir/f"{stem}.csv",FIELDS,rows); return rows,"PENDING"
    rows=[]; responses={}; references=upstream_integrals(data_mc_dir/f"{stem}.csv")
    for variable in common.VARIABLES:
        edges=configured["variables"][variable]["edges"]
        fine=payload["data"][variable]["final"]; total=payload["totals"][variable]
        fine_edges=[fine.GetXaxis().GetBinLowEdge(1+i) for i in range(fine.GetNbinsX())]+[fine.GetXaxis().GetBinUpEdge(fine.GetNbinsX())]
        response=common.response_and_moments(payload["simc_rows"]["sidis"],variable,fine_edges,stem+"_extract"); responses[variable]=response
        moments=common.weighted_bin_moments(payload["simc_rows"]["sidis"],variable,edges)
        data_all,_=aggregate(fine,0,fine.GetNbinsX()+1)
        qa={r:common.number(payload["simc_rows"][r].get("SimYield_delta")) for r in common.REACTIONS}
        data_closure=relative(data_all,references.get("Data",math.nan))
        reaction_closures={r:relative(aggregate(payload["components"][r][variable],0,fine.GetNbinsX()+1)[0],qa[r]) for r in common.REACTIONS}
        for ibin,(low,high,moment) in enumerate(zip(edges[:-1],edges[1:],moments),1):
            start,stop=fine_span(fine.GetXaxis(),low,high); dy,de=aggregate(fine,start,stop)
            values={}; errors=[]
            for reaction in common.REACTIONS:
                values[reaction],values[reaction+"_err"]=aggregate(payload["components"][reaction][variable],start,stop); errors.append(values[reaction+"_err"])
            my=sum(values[r] for r in common.REACTIONS); me=math.sqrt(sum(e*e for e in errors)); cy,ce=ratio_error(dy,de,my,me)
            purity,stability=common.coarse_response_metrics(response,start,stop); bias,brms=common.coarse_response_bias_rms(response,start,stop)
            calc=common.calculator_model(calculator,model,moment) if all(common.finite(moment.get(k)) for k in ("x_mean","Q2_mean","z_mean","pt_mean")) else {}
            sigma=common.number(calc.get("sigsemi")); mcalc=common.number(calc.get("sighad")); mevent=moment.get("M_mean",math.nan)
            row={field:"" for field in FIELDS}; row.update(identity); row.update({"Acceptance":"Delta-only","Variable":variable,"Bin_index":ibin,"Bin_low":low,"Bin_high":high,"Bin_center":.5*(low+high),"x_eff":moment["x_mean"],"x_RMS":moment["x_rms"],"Q2_eff":moment["Q2_mean"],"Q2_RMS":moment["Q2_rms"],"z_eff":moment["z_mean"],"z_RMS":moment["z_rms"],"pT_eff":moment["pt_mean"],"pT_RMS":moment["pt_rms"],"Y_Data":dy,"Y_Data_err":de,"Y_MC":my,"Y_MC_err":me,"C_Y":cy,"C_Y_err":ce,"sigma_sidis_model":sigma,"sigma_sidis_data":cy*sigma if common.finite(cy) and common.finite(sigma) else math.nan,"sigma_sidis_data_err":abs(sigma)*ce if common.finite(sigma) and common.finite(ce) else math.nan,"sigma_sidis_units":"ub/GeV^2/sr^2","M_sidis_model":mcalc,"M_sidis_event_weighted":mevent,"M_sidis_model_consistency_rel":relative(mevent,mcalc),"M_sidis_data":cy*mcalc if common.finite(cy) and common.finite(mcalc) else math.nan,"M_sidis_data_err":abs(mcalc)*ce if common.finite(mcalc) and common.finite(ce) else math.nan,"M_sidis_units":"GeV^-2","Purity":purity,"Stability":stability,"Recon_minus_original_mean":bias,"Recon_minus_original_RMS":brms,"Y_Data_rel_err":abs(de/dy) if dy else math.nan,"Y_MC_rel_err":abs(me/my) if my else math.nan,"Data_underflow":fine.GetBinContent(0),"Data_overflow":fine.GetBinContent(fine.GetNbinsX()+1),"MC_underflow":total.GetBinContent(0),"MC_overflow":total.GetBinContent(total.GetNbinsX()+1),"Y_Data_integral_closure_rel":data_closure,"Y_MC_complete":1})
            bin_reasons=[]
            for r in common.REACTIONS: row[f"Y_MC_{r}"]=values[r]; row[f"Y_MC_{r}_err"]=values[r+"_err"]; row[f"Y_MC_{r}_QA_closure_rel"]=reaction_closures[r]
            if my==0 or not common.finite(my): bin_reasons.append("ZERO_OR_NONFINITE_MC")
            if dy<0: bin_reasons.append("NEGATIVE_CORRECTED_DATA")
            if common.finite(data_closure) and abs(data_closure)>0.01: bin_reasons.append("DATA_INTEGRAL_CLOSURE")
            if any(common.finite(v) and abs(v)>1e-5 for v in reaction_closures.values()): bin_reasons.append("MC_QA_INTEGRAL_CLOSURE")
            if not (common.finite(purity) and purity>=common.PURITY_TARGET): bin_reasons.append("LOW_PURITY")
            if not (common.finite(stability) and stability>=common.STABILITY_TARGET): bin_reasons.append("LOW_STABILITY")
            row["Projection_status"]="ERROR" if "ZERO_OR_NONFINITE_MC" in bin_reasons else ("WARNING" if bin_reasons or reasons else "OK"); row["Projection_reason"]=";".join(sorted(set(reasons+bin_reasons))); rows.append(row)
    common.atomic_csv(table_dir/f"{stem}.csv",FIELDS,rows)
    pdf_dir.mkdir(parents=True,exist_ok=True); pdf=pdf_dir/f"{stem}.pdf"; canvas=ROOT.TCanvas(f"c1d_{abs(hash(stem))}","1D",1800,1100); final="ERROR" if any(r["Projection_status"]=="ERROR" for r in rows) else ("WARNING" if any(r["Projection_status"]=="WARNING" for r in rows) else "OK")
    draw_metadata(canvas,pdf,identity,final,reasons,common.BINNING_VERSION)
    for variable in common.VARIABLES: draw_pages(canvas,pdf,variable,rows,responses[variable])
    draw_closure_page(canvas,pdf,rows,final)
    canvas.Print(str(pdf)+")"); return rows,final


def main():
    project=Path(__file__).resolve().parents[1]; parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("input",nargs="?",type=Path); parser.add_argument("--all",action="store_true"); parser.add_argument("--norm-dir",type=Path,default=dmc.DEFAULT_NORM_DIR); parser.add_argument("--data-mc-dir",type=Path,default=project/"results/Tables/RP_data_mc_compare"); parser.add_argument("--simc-catalog",type=Path,default=dmc.DEFAULT_SIMC_CATALOG); parser.add_argument("--model-catalog",type=Path,default=project/"results/Tables/RP_sidis_model/RP_sidis_model.csv"); parser.add_argument("--integrated-summary",type=Path,default=project/"results/Tables/RP_extract_sidis_xsec/RP_extract_sidis_xsec_Summary.csv"); parser.add_argument("--config",type=Path,default=project/"results/Tables/RP_scan_sidis_binning/RP_sidis_binning_v1.json"); parser.add_argument("--calculator",type=Path,default=Path("/Users/radwanparvez/Documents/JLab/simc_gfortran/util/sidisxsec/calc_semi_xsec")); parser.add_argument("--t7-root",type=Path,default=dmc.DEFAULT_T7_ROOT); parser.add_argument("--table-dir",type=Path,default=project/"results/Tables/RP_extract_sidis_1d"); parser.add_argument("--pdf-dir",type=Path,default=project/"results/PDFs/RP_extract_sidis_1d"); args=parser.parse_args()
    if bool(args.input)==bool(args.all): parser.error("provide exactly one input CSV or --all")
    groups,keys=common.discover_keys(args.norm_dir,args.input,args.all); catalog=common.catalog_index(args.simc_catalog); models=common.model_index(args.model_catalog); config=json.loads(args.config.read_text()); counts={}; all_rows=[]
    for i,key in enumerate(keys,1):
        try: rows,status=process(key,groups.get((key[0],key[1],key[2],key[4]),{}),catalog,models,config,args.calculator,args.t7_root,args.data_mc_dir,args.table_dir,args.pdf_dir)
        except Exception as error: identity=common.identity_from_key(key); status="ERROR"; rows=[status_row(identity,status,f"EXTRACTION_FAILED={type(error).__name__}:{error}")]; common.atomic_csv(args.table_dir/f"{common.setting_stem(identity)}.csv",FIELDS,rows)
        all_rows.extend(rows); counts[status]=counts.get(status,0)+1; print(f"[{i}/{len(keys)}] {status}: {common.setting_stem(common.identity_from_key(key))}")
    common.atomic_csv(args.table_dir/"RP_extract_sidis_1d_Summary.csv",FIELDS,all_rows); common.atomic_csv(args.table_dir/"RP_extract_sidis_1d_Problematic.csv",FIELDS,[r for r in all_rows if r.get("Projection_status") in {"WARNING","ERROR"}])
    if args.all: draw_integrated_overview(args.pdf_dir/"RP_extract_sidis_1d_IntegratedOverview.pdf",args.integrated_summary,models)
    print("Status:",counts); return 2 if counts.get("ERROR") else 0


if __name__=="__main__": raise SystemExit(main())
