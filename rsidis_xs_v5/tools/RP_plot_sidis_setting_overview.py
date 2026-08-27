#!/usr/bin/env python3
"""Plot setting-integrated SIDIS cross-section and multiplicity overviews."""

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
ROOT.TGaxis.SetMaxDigits(4)

IDENTITY = ("Phase", "Pass", "Run_type", "Target", "Setting")
USABLE = {"OK", "WARNING"}
CANVAS_WIDTH, CANVAS_HEIGHT = 1800, 1000
PAD_LEFT_MARGIN, PAD_RIGHT_MARGIN = 0.23, 0.06
FIELDS = [
    *IDENTITY, "Ebeam_GeV", "x", "Q2_GeV2", "z", "theta_pq_deg",
    "nu_GeV", "p_pion_GeV", "pT_GeV", "pT2_GeV2",
    "Y_Data", "Y_Data_err", "Y_MC", "Y_MC_err", "C_Y", "C_Y_err",
    "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
    "sigma_sidis_units", "M_sidis_model", "M_sidis_data",
    "M_sidis_data_err", "M_sidis_units", "Extraction_status",
    "Extraction_reason", "Overview_status", "Overview_reason",
    "Source_extraction_summary", "Source_model_catalog",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def atomic_csv(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def derive_pt2(p_pion: float, theta_deg: float) -> float:
    if not all(math.isfinite(value) for value in (p_pion, theta_deg)):
        return math.nan
    return (p_pion * math.sin(math.radians(theta_deg))) ** 2


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
        row.update({field: extracted.get(field, "") for field in (
            "Y_Data", "Y_Data_err", "Y_MC", "Y_MC_err", "C_Y", "C_Y_err",
            "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
            "sigma_sidis_units", "M_sidis_model", "M_sidis_data",
            "M_sidis_data_err", "M_sidis_units", "Extraction_status", "Extraction_reason",
        )})
        row["Source_extraction_summary"] = str(extraction_path.resolve())
        row["Source_model_catalog"] = str(model_path.resolve())
        if model is None:
            row.update({"Overview_status": "ERROR", "Overview_reason": "MODEL_ROW_MISSING"})
            output.append(row); continue
        for field in ("Ebeam_GeV", "x", "Q2_GeV2", "z", "theta_pq_deg",
                      "nu_GeV", "p_pion_GeV", "pT_GeV"):
            row[field] = model.get(field, "")
        pt2 = derive_pt2(number(model.get("p_pion_GeV")), number(model.get("theta_pq_deg")))
        row["pT2_GeV2"] = pt2
        status = extracted.get("Extraction_status", "") or "PENDING"
        reasons = []
        catalog_pt = number(model.get("pT_GeV"))
        if finite(pt2) and finite(catalog_pt) and not math.isclose(pt2, catalog_pt**2, rel_tol=1e-10, abs_tol=1e-12):
            status = "ERROR"; reasons.append("PT2_CATALOG_MISMATCH")
        if status in USABLE:
            required = ("x", "Q2_GeV2", "z", "theta_pq_deg", "p_pion_GeV",
                        "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
                        "M_sidis_model", "M_sidis_data", "M_sidis_data_err")
            if not all(finite(row.get(field)) for field in required):
                status = "ERROR"; reasons.append("NONFINITE_PHYSICS_VALUE")
        row["Overview_status"] = status
        row["Overview_reason"] = ";".join(reasons) if reasons else extracted.get("Extraction_reason", "")
        output.append(row)
    return output


def controlled_groups(rows: Sequence[Mapping[str, object]], variable: str):
    if variable == "z":
        fixed = ("Phase", "Pass", "Run_type", "Target", "x", "Q2_GeV2", "theta_pq_deg")
        coordinate = "z"
    elif variable == "pT2":
        fixed = ("Phase", "Pass", "Run_type", "Target", "x", "Q2_GeV2", "z")
        coordinate = "pT2_GeV2"
    else:
        raise ValueError(variable)
    groups = defaultdict(list)
    for row in rows:
        if row.get("Overview_status") not in USABLE or not finite(row.get(coordinate)):
            continue
        groups[tuple(str(row.get(field, "")) for field in fixed)].append(row)
    result = []
    for key, members in groups.items():
        members.sort(key=lambda row: number(row[coordinate]))
        if len({number(row[coordinate]) for row in members}) >= 2:
            result.append((fixed, key, coordinate, members))
    return sorted(result, key=lambda item: item[1])


def padded_range(values: Sequence[float], errors: Sequence[float]) -> tuple[float, float]:
    lows = [value - error for value, error in zip(values, errors) if math.isfinite(value) and math.isfinite(error)]
    highs = [value + error for value, error in zip(values, errors) if math.isfinite(value) and math.isfinite(error)]
    if not lows: return 0.0, 1.0
    low, high = min(lows), max(highs)
    span = max(high - low, abs(high) * 0.15, abs(low) * 0.15, 1e-12)
    return min(0.0, low - 0.15 * span), high + 0.20 * span


def make_graph(rows, xfield, yfield, error_field, color, marker):
    graph = ROOT.TGraphErrors(len(rows))
    for index, row in enumerate(rows):
        graph.SetPoint(index, number(row[xfield]), number(row[yfield]))
        graph.SetPointError(index, 0.0, number(row.get(error_field), 0.0) if error_field else 0.0)
    graph.SetMarkerStyle(marker); graph.SetMarkerSize(1.65)
    graph.SetMarkerColor(color); graph.SetLineColor(color); graph.SetLineWidth(2)
    return graph


def fixed_title(fixed: Sequence[str], key: Sequence[str]) -> str:
    labels = {"Phase":"", "Pass":"", "Run_type":"", "Target":"",
              "x":"x=", "Q2_GeV2":"Q^{2}=", "z":"z=", "theta_pq_deg":"#theta_{pq}="}
    values = dict(zip(fixed, key)); parts = [f"{values['Phase']}, {values['Pass']}",
                                             f"{values['Run_type']}, {values['Target']}"]
    for field, value in zip(fixed, key):
        if field in {"Phase", "Pass", "Run_type", "Target"}: continue
        suffix = " GeV^{2}" if field == "Q2_GeV2" else (" deg" if field == "theta_pq_deg" else "")
        parts.append(f"{labels[field]}{value}{suffix}")
    return ", ".join(parts)


def draw_panel(pad, rows, xfield, panel_title, group_title, ymodel, ydata, yerror, ylabel):
    pad.cd(); pad.SetLeftMargin(PAD_LEFT_MARGIN); pad.SetRightMargin(PAD_RIGHT_MARGIN)
    pad.SetBottomMargin(0.14); pad.SetTopMargin(0.17); pad.SetGrid(1, 1)
    xs = [number(row[xfield]) for row in rows]
    xspan = max(max(xs)-min(xs), abs(max(xs))*0.1, 1e-6)
    data = [number(row[ydata]) for row in rows]; errors = [number(row[yerror], 0.0) for row in rows]
    model = [number(row[ymodel]) for row in rows]
    ymin, ymax = padded_range(data + model, errors + [0.0]*len(model))
    frame = ROOT.TH1D(f"frame_{abs(hash((panel_title,ydata,tuple(xs))))}", "", 100,
                      min(xs)-0.10*xspan, max(xs)+0.10*xspan)
    frame.SetDirectory(0); frame.SetStats(False); frame.SetMinimum(ymin); frame.SetMaximum(ymax)
    frame.GetXaxis().SetTitle("Central z" if xfield == "z" else "Central p_{T}^{2} (GeV^{2})")
    frame.GetYaxis().SetTitle(ylabel); frame.GetYaxis().SetDecimals(True)
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetTitleFont(42); axis.SetLabelFont(42); axis.SetTitleSize(0.047); axis.SetLabelSize(0.039)
    frame.GetXaxis().SetTitleOffset(1.25); frame.GetYaxis().SetTitleOffset(1.60)
    frame.Draw("AXIS")
    extracted = make_graph(rows, xfield, ydata, yerror, ROOT.kBlue+1, 20)
    model_graph = make_graph(rows, xfield, ymodel, "", ROOT.kRed+1, 24)
    extracted.Draw("P SAME"); model_graph.Draw("P SAME")
    legend = ROOT.TLegend(0.69, 0.68, 0.92, 0.81)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.038)
    legend.AddEntry(extracted, "Extracted", "lep"); legend.AddEntry(model_graph, "Model", "p"); legend.Draw()
    heading = ROOT.TLatex(); heading.SetNDC(True); heading.SetTextFont(42)
    heading.SetTextSize(0.046); title_line = heading.DrawLatex(0.17, 0.955, panel_title)
    heading.SetTextSize(0.030); group_line = heading.DrawLatex(0.17, 0.910, group_title)
    return frame, extracted, model_graph, legend, heading, title_line, group_line


def metadata_page(canvas, pdf: Path, rows, z_groups, pt_groups, sources):
    canvas.Clear(); text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    counts = Counter(str(row["Overview_status"]) for row in rows)
    lines = [
        "Setting-Integrated SIDIS Physics Overview",
        "One point represents one complete leaf setting; no event-level physics binning is used.",
        "Selection: Delta acceptance only. Cross-sections and multiplicities are copied unchanged",
        "from RP_extract_sidis_xsec; this tool only joins metadata and derives central pT^{2}.",
        "pT^{2} = (p_{#pi} sin #theta_{pq})^{2}, p_{#pi} = sqrt((z#nu)^{2}-m_{#pi}^{2}), #nu=Q^{2}/(2M_{p}x)",
        "z scans hold x, Q^{2}, and #theta_{pq} fixed; central pT follows the setting kinematics.",
        "pT^{2} scans hold x, Q^{2}, and z fixed.",
        f"Catalog statuses: " + ", ".join(f"{key}={value}" for key,value in sorted(counts.items())),
        f"Controlled scan pages: z={len(z_groups)}, pT^{{2}}={len(pt_groups)}",
        f"Extraction source: {sources[0]}", f"Model source: {sources[1]}",
    ]
    y = 0.92
    for index, line in enumerate(lines):
        text.SetTextSize(0.043 if index == 0 else 0.030)
        text.DrawLatex(0.055, y, line); y -= 0.072
    canvas.Print(str(pdf) + "(")


def draw_pdf(path: Path, rows, extraction_path: Path, model_path: Path):
    z_groups = controlled_groups(rows, "z"); pt_groups = controlled_groups(rows, "pT2")
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas = ROOT.TCanvas("c_setting_overview", "overview", CANVAS_WIDTH, CANVAS_HEIGHT)
    metadata_page(canvas, path, rows, z_groups, pt_groups, (extraction_path, model_path))
    for variable, groups in (("z", z_groups), ("pT2", pt_groups)):
        for fixed, key, xfield, members in groups:
            canvas.Clear(); canvas.Divide(2, 1, 0.002, 0.002)
            group_title = fixed_title(fixed, key)
            owned = []
            owned.extend(draw_panel(canvas.cd(1), members, xfield,
                                    f"#sigma_{{sidis}} vs central {'z' if variable == 'z' else 'p_{T}^{2}'}",
                                    group_title,
                                    "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
                                    "#sigma_{sidis} (#mub/GeV^{2}/sr^{2})"))
            owned.extend(draw_panel(canvas.cd(2), members, xfield,
                                    f"M_{{sidis}} vs central {'z' if variable == 'z' else 'p_{T}^{2}'}",
                                    group_title,
                                    "M_sidis_model", "M_sidis_data", "M_sidis_data_err",
                                    "M_{sidis} (GeV^{-2})"))
            canvas.Print(str(path))
    canvas.Clear(); text=ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    text.SetTextSize(0.043); text.DrawLatex(0.06,0.91,"Setting-Integrated Overview Status")
    text.SetTextSize(0.032); text.DrawLatex(0.06,0.82,f"Rows retained: {len(rows)}")
    text.DrawLatex(0.06,0.75,f"Usable physics rows: {sum(row['Overview_status'] in USABLE for row in rows)}")
    text.DrawLatex(0.06,0.68,f"z scan pages: {len(z_groups)}; pT^{{2}} scan pages: {len(pt_groups)}")
    text.DrawLatex(0.06,0.61,"Single-point and unavailable groups remain in the CSV but do not create physics pages.")
    canvas.Print(str(path) + ")")


def main() -> int:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extraction", type=Path, default=project/"results/Tables/RP_extract_sidis_xsec/RP_extract_sidis_xsec_Summary.csv")
    parser.add_argument("--model", type=Path, default=project/"results/Tables/RP_sidis_model/RP_sidis_model.csv")
    parser.add_argument("--table-dir", type=Path, default=project/"results/Tables/RP_sidis_setting_overview")
    parser.add_argument("--pdf", type=Path, default=project/"results/PDFs/RP_sidis_setting_overview/RP_sidis_setting_overview.pdf")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    rows = join_catalogs(read_csv(args.extraction), read_csv(args.model), args.extraction, args.model)
    atomic_csv(args.table_dir/"RP_sidis_setting_overview.csv", rows)
    atomic_csv(args.table_dir/"RP_sidis_setting_overview_Problematic.csv",
               [row for row in rows if row["Overview_status"] in {"WARNING", "ERROR"}])
    if not args.no_pdf: draw_pdf(args.pdf, rows, args.extraction, args.model)
    counts = Counter(str(row["Overview_status"]) for row in rows)
    print("Status:", ", ".join(f"{key}={value}" for key,value in sorted(counts.items())))
    return 2 if counts.get("ERROR") else 0


if __name__ == "__main__": raise SystemExit(main())
