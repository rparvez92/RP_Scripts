#!/usr/bin/env python3
"""Plot Pass1 setting-integrated SIDIS cross-section and multiplicity overviews."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
import tempfile
from typing import Iterable, Mapping, Sequence

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPaperSize(30.0, 16.6667)
ROOT.gStyle.SetGridColor(ROOT.kGray + 1)
ROOT.gStyle.SetGridStyle(3)
ROOT.gStyle.SetGridWidth(1)
ROOT.TGaxis.SetMaxDigits(4)

IDENTITY = ("Phase", "Pass", "Run_type", "Target", "setting_id")
USABLE = {"OK", "WARNING"}
CANVAS_WIDTH, CANVAS_HEIGHT = 1800, 1000
PAD_LEFT_MARGIN, PAD_RIGHT_MARGIN = 0.23, 0.06
FIELDS = [
    *IDENTITY, "Ebeam_GeV", "x", "Q2", "z", "thpq",
    "nu_GeV", "p_pion_GeV", "Pt_GeV", "Y_Data", "Y_Data_err",
    "Y_MC", "Y_MC_err", "C_Y", "C_Y_err", "sigma_sidis_model",
    "sigma_sidis_data", "sigma_sidis_data_err", "sigma_sidis_units",
    "M_sidis_model", "M_sidis_data", "M_sidis_data_err", "M_sidis_units",
    "Extraction_status", "Extraction_reason", "Overview_status",
    "Overview_reason", "Source_extraction_summary", "Source_model_catalog",
]


def finite(value: object) -> bool:
    try: return math.isfinite(float(value))
    except (TypeError, ValueError): return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def atomic_csv(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def derive_pt(p_pion: float, theta_deg: float) -> float:
    if not all(math.isfinite(value) for value in (p_pion, theta_deg)):
        return math.nan
    return abs(p_pion * math.sin(math.radians(theta_deg)))


def relative_source(path: Path) -> str:
    """Format provenance relative to the RP_Scripts repository root."""
    repository = Path(__file__).resolve().parents[2]
    try:
        return path.resolve().relative_to(repository).as_posix()
    except ValueError:
        return path.as_posix()


def join_catalogs(extraction_rows: Sequence[Mapping[str, str]],
                  model_rows: Sequence[Mapping[str, str]],
                  extraction_path: Path, model_path: Path) -> list[dict[str, object]]:
    models = {tuple(row.get(field, "") for field in IDENTITY): row for row in model_rows}
    output = []
    for extracted in extraction_rows:
        key = tuple(extracted.get(field, "") for field in IDENTITY)
        model = models.get(key)
        row: dict[str, object] = {field: "" for field in FIELDS}
        row.update({field: extracted.get(field, "") for field in IDENTITY})
        for field in (
            "Y_Data", "Y_Data_err", "Y_MC", "Y_MC_err", "C_Y", "C_Y_err",
            "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
            "sigma_sidis_units", "M_sidis_model", "M_sidis_data",
            "M_sidis_data_err", "M_sidis_units", "Extraction_status",
            "Extraction_reason",
        ):
            row[field] = extracted.get(field, "")
        row["Source_extraction_summary"] = relative_source(extraction_path)
        row["Source_model_catalog"] = relative_source(model_path)
        if model is None:
            row.update(Overview_status="ERROR", Overview_reason="MODEL_ROW_MISSING")
            output.append(row); continue
        for output_field, model_field in (
            ("Ebeam_GeV", "Ebeam_GeV"), ("x", "x"), ("Q2", "Q2_GeV2"),
            ("z", "z"), ("thpq", "theta_pq_deg"), ("nu_GeV", "nu_GeV"),
            ("p_pion_GeV", "p_pion_GeV"),
        ):
            row[output_field] = model.get(model_field, "")
        row["Pt_GeV"] = model.get("pT_GeV", "")
        status = extracted.get("Extraction_status", "") or "ERROR"
        reasons = []
        derived = derive_pt(number(model.get("p_pion_GeV")),
                            number(model.get("theta_pq_deg")))
        catalog = number(model.get("pT_GeV"))
        if finite(derived) and finite(catalog) and not math.isclose(
                derived, catalog, rel_tol=1e-10, abs_tol=1e-12):
            status = "ERROR"; reasons.append("PT_CATALOG_MISMATCH")
        if status in USABLE:
            required = (
                "x", "Q2", "z", "thpq", "p_pion_GeV", "Pt_GeV",
                "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
                "M_sidis_model", "M_sidis_data", "M_sidis_data_err",
            )
            if not all(finite(row.get(field)) for field in required):
                status = "ERROR"; reasons.append("NONFINITE_PHYSICS_VALUE")
        row["Overview_status"] = status
        row["Overview_reason"] = ";".join(reasons) if reasons else extracted.get("Extraction_reason", "")
        output.append(row)
    return output


def controlled_groups(rows: Sequence[Mapping[str, object]], variable: str):
    if variable == "z":
        fixed = ("Phase", "Pass", "Run_type", "Target", "x", "Q2", "thpq")
        coordinate = "z"
    elif variable == "Pt":
        fixed = ("Phase", "Pass", "Run_type", "Target", "x", "Q2", "z")
        coordinate = "Pt_GeV"
    else:
        raise ValueError(variable)
    groups = defaultdict(list)
    for row in rows:
        if row.get("Overview_status") in USABLE and finite(row.get(coordinate)):
            groups[tuple(str(row.get(field, "")) for field in fixed)].append(row)
    result = []
    for key, members in groups.items():
        members.sort(key=lambda row: number(row[coordinate]))
        if len({number(row[coordinate]) for row in members}) >= 2:
            result.append((fixed, key, coordinate, members))
    return sorted(result, key=lambda item: item[1])


def padded_range(values: Sequence[float], errors: Sequence[float]) -> tuple[float, float]:
    lows = [value-error for value,error in zip(values,errors) if finite(value) and finite(error)]
    highs = [value+error for value,error in zip(values,errors) if finite(value) and finite(error)]
    if not lows: return 0.0, 1.0
    low, high = min(lows), max(highs)
    span = max(high-low, abs(high)*0.15, abs(low)*0.15, 1e-12)
    return min(0.0, low-0.15*span), high+0.20*span


def make_graph(rows, xfield, yfield, error_field, color, marker):
    graph = ROOT.TGraphErrors(len(rows))
    for index, row in enumerate(rows):
        graph.SetPoint(index, number(row[xfield]), number(row[yfield]))
        graph.SetPointError(index, 0.0, number(row.get(error_field), 0.0) if error_field else 0.0)
    graph.SetMarkerStyle(marker); graph.SetMarkerSize(1.65)
    graph.SetMarkerColor(color); graph.SetLineColor(color); graph.SetLineWidth(2)
    return graph


def fixed_title(fixed: Sequence[str], key: Sequence[str]) -> str:
    labels = {"x":"x=", "Q2":"Q2=", "z":"z=", "thpq":"thpq="}
    values = dict(zip(fixed, key))
    parts = [f"{values['Phase']}, {values['Pass']}",
             f"{values['Run_type']}, {values['Target']}"]
    for field, value in zip(fixed, key):
        if field not in labels: continue
        suffix = " GeV^{2}" if field == "Q2" else (" deg" if field == "thpq" else "")
        parts.append(f"{labels[field]}{value}{suffix}")
    return ", ".join(parts)


def draw_panel(pad, rows, xfield, panel_title, group_title,
               ymodel, ydata, yerror, ylabel):
    pad.cd(); pad.SetLeftMargin(PAD_LEFT_MARGIN); pad.SetRightMargin(PAD_RIGHT_MARGIN)
    pad.SetBottomMargin(0.14); pad.SetTopMargin(0.17)
    pad.SetGridx(True); pad.SetGridy(True)
    xs = [number(row[xfield]) for row in rows]
    xlow, xhigh = (0.0, 1.0) if xfield == "z" else (0.0, 0.6)
    data = [number(row[ydata]) for row in rows]
    errors = [number(row[yerror], 0.0) for row in rows]
    model = [number(row[ymodel]) for row in rows]
    ymin, ymax = padded_range(data+model, errors+[0.0]*len(model))
    frame = ROOT.TH1D(f"frame_{abs(hash((panel_title,ydata,tuple(xs))))}", "", 100,
                      xlow, xhigh)
    frame.SetDirectory(0); frame.SetMinimum(ymin); frame.SetMaximum(ymax)
    frame.GetXaxis().SetTitle("Central z" if xfield == "z" else "Central P_{t} (GeV)")
    # Ten labeled primary divisions give 0.1 steps on both fixed axes. ROOT may
    # draw minor ticks between them, but does not label those intermediate ticks.
    frame.GetXaxis().SetNdivisions(510, True)
    frame.GetYaxis().SetTitle(ylabel); frame.GetYaxis().SetDecimals(True)
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetTitleFont(42); axis.SetLabelFont(42)
        axis.SetTitleSize(0.047); axis.SetLabelSize(0.039)
    frame.GetXaxis().SetTitleOffset(1.25); frame.GetYaxis().SetTitleOffset(1.60)
    frame.Draw()
    extracted = make_graph(rows, xfield, ydata, yerror, ROOT.kBlue+1, 20)
    model_graph = make_graph(rows, xfield, ymodel, "", ROOT.kRed+1, 24)
    extracted.Draw("P SAME"); model_graph.Draw("P SAME")
    legend = ROOT.TLegend(0.69,0.68,0.92,0.81)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.038)
    legend.AddEntry(extracted,"Extracted","lep"); legend.AddEntry(model_graph,"Model","p"); legend.Draw()
    heading = ROOT.TLatex(); heading.SetNDC(); heading.SetTextFont(42)
    heading.SetTextSize(0.046); title_line=heading.DrawLatex(0.17,0.955,panel_title)
    heading.SetTextSize(0.030); group_line=heading.DrawLatex(0.17,0.910,group_title)
    return frame, extracted, model_graph, legend, heading, title_line, group_line


def metadata_page(canvas, pdf: Path, rows, z_groups, pt_groups, sources, quantity):
    canvas.Clear(); text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42)
    counts=Counter(str(row["Overview_status"]) for row in rows)
    is_sigma = quantity == "sigma"
    physics_lines = (
        ["Quantity: extracted and model SIDIS cross section.",
         "Cross-section units: ub/GeV2/sr2."]
        if is_sigma else
        ["Quantity: extracted and model SIDIS multiplicity.",
         "Multiplicity units: GeV^-2."]
    )
    lines = [
        f"Setting-Integrated SIDIS {'Cross-Section' if is_sigma else 'Multiplicity'} Overview",
        "One point represents one complete setting; no event-level physics binning is used.",
        "Selection: Baseline PID + Delta-only acceptance.",
        *physics_lines,
        "P_{t} = |p_{#pi} sin #theta_{pq}|, p_{#pi}=sqrt((z#nu)^{2}-m_{#pi}^{2}), #nu=Q^{2}/(2M_{p}x)",
        "z scans hold x, Q^{2}, and #theta_{pq} fixed.",
        "P_{t} scans hold x, Q^{2}, and z fixed.",
        "Statuses: " + ", ".join(f"{key}={value}" for key,value in sorted(counts.items())),
        f"Controlled scan pages: z={len(z_groups)}, P_{{t}}={len(pt_groups)}",
        f"Extraction source: {relative_source(sources[0])}",
        f"Model source: {relative_source(sources[1])}",
    ]
    y=0.92
    for index,line in enumerate(lines):
        text.SetTextSize(0.043 if index==0 else 0.030)
        text.SetTextFont(62 if index==0 else 42); text.DrawLatex(0.055,y,line); y-=0.072
    canvas.Print(str(pdf))


def draw_unavailable(pad, label: str):
    pad.cd(); pad.SetGrid(1, 1)
    text = ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42); text.SetTextAlign(22)
    text.SetTextSize(0.050); text.DrawLatex(0.5, 0.55, label)
    text.SetTextSize(0.038); text.DrawLatex(0.5, 0.47, "No controlled multi-setting scan")
    return text


def new_page_canvas(name: str) -> ROOT.TCanvas:
    """Create every PDF page from the same locked canvas geometry."""
    canvas = ROOT.TCanvas(name, "overview", CANVAS_WIDTH, CANVAS_HEIGHT)
    canvas.SetCanvasSize(CANVAS_WIDTH, CANVAS_HEIGHT)
    return canvas


def draw_pdf(path: Path, rows, extraction_path: Path, model_path: Path, quantity: str):
    z_groups=controlled_groups(rows,"z"); pt_groups=controlled_groups(rows,"Pt")
    path.parent.mkdir(parents=True,exist_ok=True)
    opener = new_page_canvas(f"c_setting_overview_{quantity}_open")
    opener.Print(str(path)+"["); opener.Close()
    canvas = new_page_canvas(f"c_setting_overview_{quantity}_metadata")
    metadata_page(canvas,path,rows,z_groups,pt_groups,(extraction_path,model_path),quantity)
    canvas.Close()
    if quantity == "sigma":
        ymodel, ydata, yerror = "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err"
        symbol, ylabel = "#sigma_{sidis}", "#sigma_{sidis} (#mub/GeV^{2}/sr^{2})"
    else:
        ymodel, ydata, yerror = "M_sidis_model", "M_sidis_data", "M_sidis_data_err"
        symbol, ylabel = "M_{sidis}", "M_{sidis} (GeV^{-2})"
    for index in range(max(len(z_groups), len(pt_groups))):
        canvas = new_page_canvas(f"c_setting_overview_{quantity}_plot_{index}")
        canvas.Divide(2,1,0.002,0.002); owned=[]
        if index < len(z_groups):
            fixed,key,xfield,members=z_groups[index]
            owned.extend(draw_panel(canvas.cd(1),members,xfield,
                                    f"{symbol} vs central z",fixed_title(fixed,key),
                                    ymodel,ydata,yerror,ylabel))
        else:
            owned.append(draw_unavailable(canvas.cd(1), f"{symbol} vs central z"))
        if index < len(pt_groups):
            fixed,key,xfield,members=pt_groups[index]
            owned.extend(draw_panel(canvas.cd(2),members,xfield,
                                    f"{symbol} vs central P_{{t}}",fixed_title(fixed,key),
                                    ymodel,ydata,yerror,ylabel))
        else:
            owned.append(draw_unavailable(canvas.cd(2), f"{symbol} vs central P_{{t}}"))
        canvas.Print(str(path)); canvas.Close()
    canvas = new_page_canvas(f"c_setting_overview_{quantity}_status")
    text=ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42)
    text.SetTextSize(0.043); text.DrawLatex(0.06,0.91,
        f"Setting-Integrated {'Cross-Section' if quantity == 'sigma' else 'Multiplicity'} Status")
    text.SetTextSize(0.032); text.DrawLatex(0.06,0.82,f"Rows retained: {len(rows)}")
    text.DrawLatex(0.06,0.75,f"Usable physics rows: {sum(row['Overview_status'] in USABLE for row in rows)}")
    text.DrawLatex(0.06,0.68,f"z scan pages: {len(z_groups)}; P_{{t}} scan pages: {len(pt_groups)}")
    text.DrawLatex(0.06,0.61,"Single-point groups remain in the CSV but do not create physics pages.")
    canvas.Print(str(path)); canvas.Print(str(path)+"]"); canvas.Close()


def main() -> int:
    project=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extraction",type=Path,default=project/"results/Pass1/Tables/RP_extract_sidis_xsec/RP_extract_sidis_xsec_Summary.csv")
    parser.add_argument("--model",type=Path,default=project/"results/Pass1/Tables/RP_sidis_model/RP_sidis_model.csv")
    parser.add_argument("--table-dir",type=Path,default=project/"results/Pass1/Tables/RP_sidis_setting_overview")
    parser.add_argument("--pdf-dir",type=Path,default=project/"results/Pass1/PDFs/RP_sidis_setting_overview")
    parser.add_argument("--no-pdf",action="store_true")
    args=parser.parse_args()
    rows=join_catalogs(read_csv(args.extraction),read_csv(args.model),args.extraction,args.model)
    atomic_csv(args.table_dir/"RP_sidis_setting_overview.csv",rows)
    atomic_csv(args.table_dir/"RP_sidis_setting_overview_Problematic.csv",
               [row for row in rows if row["Overview_status"]!="OK"])
    if not args.no_pdf:
        draw_pdf(args.pdf_dir/"RP_sidis_setting_overview_sigma.pdf",
                 rows,args.extraction,args.model,"sigma")
        draw_pdf(args.pdf_dir/"RP_sidis_setting_overview_multiplicity.pdf",
                 rows,args.extraction,args.model,"multiplicity")
    counts=Counter(str(row["Overview_status"]) for row in rows)
    print("Status:",", ".join(f"{key}={value}" for key,value in sorted(counts.items())))
    return 2 if counts.get("ERROR") else 0


if __name__=="__main__": raise SystemExit(main())
