#!/usr/bin/env python3
"""Scan Delta-acceptance SIDIS occupancy/migration and recommend 1D bins."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import ROOT

import RP_data_mc_compare as dmc
import RP_sidis_differential_common as common

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

FIELDS = [
    *common.IDENTITY_FIELDS, "Variable", "Fine_bin", "Bin_low", "Bin_high",
    "Y_Data", "Y_Data_err", "Y_Data_rel_err",
    "Y_MC_sidis", "Y_MC_sidis_err", "Y_MC_rho", "Y_MC_rho_err",
    "Y_MC_delta", "Y_MC_delta_err", "Y_MC_exclusive", "Y_MC_exclusive_err",
    "Y_MC", "Y_MC_err", "Y_MC_rel_err", "Purity", "Stability",
    "Recon_minus_original_mean", "Recon_minus_original_RMS",
    "Data_underflow", "Data_overflow", "MC_underflow", "MC_overflow",
    "Recommended_bin", "Recommended_low", "Recommended_high",
    "Quality_pass", "Scan_status", "Scan_reason",
]


def bin_value(hist, index):
    return hist.GetBinContent(index), hist.GetBinError(index)


def flow(hist):
    return hist.GetBinContent(0), hist.GetBinContent(hist.GetNbinsX() + 1)


def recommended_bin_for(index: int, details, edges):
    for number, detail in enumerate(details, 1):
        if int(detail["start"]) <= index <= int(detail["stop"]):
            return number, edges[number - 1], edges[number], int(bool(detail["passes"]))
    return "", "", "", 0


def metadata_page(canvas, pdf: Path, identity, status, reasons, config):
    canvas.Clear()
    text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    lines = [
        "SIDIS Differential Occupancy and Migration Scan",
        ", ".join(identity[field] for field in common.IDENTITY_FIELDS),
        "Selection: Delta-only (Data includes PID; MC uses Delta acceptance)",
        "Fine grids: z=(20,0,1), pT^{2}=(20,0,0.6 GeV^{2}), phi_{pq}=(24,0,2pi)",
        "Recommended bins require complete finite MC, MC rel. error <=10%,",
        "Data rel. error <=30% when defined, purity >=50%, and stability >=50%.",
        "Migration matrices use |fWeight|; physics yields use signed fWeight.",
        "No Neff criterion or warning is used.",
        f"Status: {status}",
        f"Reason: {';'.join(reasons) if reasons else 'none'}",
    ]
    for variable in common.VARIABLES:
        item = config.get(variable)
        if item:
            lines.append(f"{variable} recommended edges: " + ", ".join(f"{x:.5g}" for x in item["edges"]))
    y = 0.93
    for index, line in enumerate(lines):
        text.SetTextSize(0.040 if index == 0 else 0.029)
        text.DrawLatex(0.06, y, line); y -= 0.058
    canvas.Print(str(pdf) + "(")


def draw_variable_pages(canvas, pdf: Path, variable: str, payload, response):
    data = payload["data"][variable]["final"]
    components = payload["components"]
    total = payload["totals"][variable]
    canvas.Clear(); canvas.SetLogy(False)
    data.SetTitle(f"Fine-bin occupancy: {variable};{data.GetXaxis().GetTitle()};Normalized Yield")
    data.SetMarkerStyle(20); data.SetMarkerColor(ROOT.kBlack); data.SetLineColor(ROOT.kBlack)
    ymax = max(data.GetMaximum(), total.GetMaximum(), 1e-12) * 1.35
    data.SetMaximum(ymax); data.Draw("E1")
    colors = {"sidis": ROOT.kBlue + 1, "rho": ROOT.kGreen + 2,
              "delta": ROOT.kMagenta + 1, "exclusive": ROOT.kOrange + 7}
    legend = ROOT.TLegend(0.72, 0.68, 0.94, 0.91); legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(data, "Data", "lep")
    for reaction in common.REACTIONS:
        hist = components[reaction][variable]
        hist.SetLineColor(colors[reaction]); hist.SetLineWidth(2); hist.Draw("HIST SAME")
        legend.AddEntry(hist, reaction, "l")
    total.SetLineColor(ROOT.kRed + 1); total.SetLineWidth(3); total.Draw("HIST SAME")
    legend.AddEntry(total, "MC", "l"); legend.Draw(); canvas.SetGrid(True, True)
    canvas.Print(str(pdf))

    canvas.Clear()
    n = data.GetNbinsX()
    frame = ROOT.TH1D(f"scan_metric_{variable}_{abs(hash(str(pdf)))}", f"Statistical precision and migration: {variable};Fine bin;Metric", n, 0.5, n + 0.5)
    frame.SetStats(False); frame.SetMinimum(0); frame.SetMaximum(1.2); frame.Draw()
    graphs = []
    specs = (("Data rel. error", ROOT.kBlack, 20), ("MC rel. error", ROOT.kRed + 1, 21),
             ("Purity", ROOT.kBlue + 1, 22), ("Stability", ROOT.kGreen + 2, 23))
    for label, color, marker in specs:
        graph = ROOT.TGraph(n); graph.SetName(f"g_{variable}_{marker}_{abs(hash(str(pdf)))}")
        for index in range(1, n + 1):
            dv, de = bin_value(data, index); mv, me = bin_value(total, index)
            purity, stability = common.response_metrics(response, index)
            metric = {"Data rel. error": abs(de / dv) if dv else math.nan,
                      "MC rel. error": abs(me / mv) if mv else math.nan,
                      "Purity": purity, "Stability": stability}[label]
            graph.SetPoint(index - 1, index, metric)
        graph.SetMarkerStyle(marker); graph.SetMarkerColor(color); graph.SetLineColor(color)
        graph.Draw("PL SAME"); graphs.append((label, graph))
    legend = ROOT.TLegend(0.72, 0.70, 0.94, 0.91); legend.SetBorderSize(0); legend.SetFillStyle(0)
    for label, graph in graphs: legend.AddEntry(graph, label, "lp")
    legend.Draw(); canvas.SetGrid(True, True); canvas.Print(str(pdf))

    canvas.Clear(); canvas.SetRightMargin(0.14)
    response.SetTitle(f"SIDIS response: {variable};Original;Reconstructed")
    response.Draw("COLZ"); canvas.Print(str(pdf)); canvas.SetRightMargin(0.04)


def process_setting(key, files, catalog, t7_root: Path, table_dir: Path, pdf_dir: Path):
    identity = common.identity_from_key(key); stem = common.setting_stem(identity)
    if "neg999" in identity["Setting"]:
        status, reasons = "SKIPPED", ["SENTINEL_SETTING_IDENTITY"]
        row = {field: "" for field in FIELDS}; row.update(identity)
        row.update({"Scan_status": status, "Scan_reason": reasons[0]})
        common.atomic_csv(table_dir / f"{stem}.csv", FIELDS, [row])
        return [row], {"identity": identity, "status": status, "reason": reasons[0], "variables": {}}, status
    payload, status, reasons = common.prepare_histograms(
        key, files, catalog, common.FINE_BINNING, t7_root
    )
    setting_config = {"identity": identity, "status": status, "reason": ";".join(reasons), "variables": {}}
    if payload is None:
        row = {field: "" for field in FIELDS}; row.update(identity)
        row.update({"Scan_status": status, "Scan_reason": ";".join(reasons)})
        common.atomic_csv(table_dir / f"{stem}.csv", FIELDS, [row])
        return [row], setting_config, status
    rows = []
    responses = {}
    for variable in common.VARIABLES:
        fine = common.FINE_BINNING[variable]
        width = (fine[2] - fine[1]) / fine[0]
        edges = [fine[1] + index * width for index in range(fine[0] + 1)]
        response = common.response_and_moments(
            payload["simc_rows"]["sidis"], variable, edges, stem
        )
        responses[variable] = response
        data = payload["data"][variable]["final"]; total = payload["totals"][variable]
        if variable == "phipq":
            recommended, details = common.recommend_phi_edges(data, total, response)
        else:
            recommended, details = common.recommend_adjacent_edges(data, total, response)
        setting_config["variables"][variable] = {
            "edges": recommended,
            "quality": details,
            "all_bins_pass": all(bool(item["passes"]) for item in details),
        }
        data_under, data_over = flow(data); mc_under, mc_over = flow(total)
        for index in range(1, data.GetNbinsX() + 1):
            data_y, data_e = bin_value(data, index); mc_y, mc_e = bin_value(total, index)
            purity, stability = common.response_metrics(response, index)
            bias, bias_rms = common.response_bias_rms(response, index)
            rec_bin, rec_low, rec_high, quality = recommended_bin_for(index, details, recommended)
            row = {field: "" for field in FIELDS}; row.update(identity)
            row.update({
                "Variable": variable, "Fine_bin": index,
                "Bin_low": data.GetXaxis().GetBinLowEdge(index),
                "Bin_high": data.GetXaxis().GetBinUpEdge(index),
                "Y_Data": data_y, "Y_Data_err": data_e,
                "Y_Data_rel_err": abs(data_e / data_y) if data_y else math.nan,
                "Y_MC": mc_y, "Y_MC_err": mc_e,
                "Y_MC_rel_err": abs(mc_e / mc_y) if mc_y else math.nan,
                "Purity": purity, "Stability": stability,
                "Recon_minus_original_mean": bias, "Recon_minus_original_RMS": bias_rms,
                "Data_underflow": data_under, "Data_overflow": data_over,
                "MC_underflow": mc_under, "MC_overflow": mc_over,
                "Recommended_bin": rec_bin, "Recommended_low": rec_low,
                "Recommended_high": rec_high, "Quality_pass": quality,
                "Scan_status": status, "Scan_reason": ";".join(reasons),
            })
            for reaction in common.REACTIONS:
                value, error = bin_value(payload["components"][reaction][variable], index)
                row[f"Y_MC_{reaction}"] = value; row[f"Y_MC_{reaction}_err"] = error
            rows.append(row)
    table_dir.mkdir(parents=True, exist_ok=True)
    common.atomic_csv(table_dir / f"{stem}.csv", FIELDS, rows)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    pdf = pdf_dir / f"{stem}.pdf"; canvas = ROOT.TCanvas(f"c_scan_{abs(hash(stem))}", "scan", 1800, 1100)
    metadata_page(canvas, pdf, identity, status, reasons, setting_config["variables"])
    for variable in common.VARIABLES:
        draw_variable_pages(canvas, pdf, variable, payload, responses[variable])
    canvas.Print(str(pdf) + ")")
    return rows, setting_config, status


def main() -> int:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--norm-dir", type=Path, default=dmc.DEFAULT_NORM_DIR)
    parser.add_argument("--simc-catalog", type=Path, default=dmc.DEFAULT_SIMC_CATALOG)
    parser.add_argument("--t7-root", type=Path, default=dmc.DEFAULT_T7_ROOT)
    parser.add_argument("--table-dir", type=Path, default=project / "results/Tables/RP_scan_sidis_binning")
    parser.add_argument("--pdf-dir", type=Path, default=project / "results/PDFs/RP_scan_sidis_binning")
    parser.add_argument("--config", type=Path, default=project / "results/Tables/RP_scan_sidis_binning/RP_sidis_binning_v1.json")
    args = parser.parse_args()
    if bool(args.input) == bool(args.all): parser.error("provide exactly one input CSV or --all")
    groups, keys = common.discover_keys(args.norm_dir, args.input, args.all)
    catalog = common.catalog_index(args.simc_catalog)
    config = {
        "version": common.BINNING_VERSION, "selection": "Delta-only",
        "thresholds": {"mc_rel_err_max": common.MC_REL_TARGET,
                       "data_rel_err_max": common.DATA_REL_TARGET,
                       "purity_min": common.PURITY_TARGET,
                       "stability_min": common.STABILITY_TARGET},
        "fine_binning": {key: list(value[:3]) for key, value in common.FINE_BINNING.items()},
        "settings": {},
    }
    counts = {}; all_rows = []
    for index, key in enumerate(keys, 1):
        rows = []
        files = groups.get((key[0], key[1], key[2], key[4]), {})
        try:
            rows, setting_config, status = process_setting(
                key, files, catalog, args.t7_root, args.table_dir, args.pdf_dir
            )
        except Exception as error:
            identity = common.identity_from_key(key); status = "ERROR"
            setting_config = {"identity": identity, "status": status,
                              "reason": f"SCAN_FAILED={type(error).__name__}:{error}", "variables": {}}
            row = {field: "" for field in FIELDS}; row.update(identity)
            row.update({"Scan_status": status, "Scan_reason": setting_config["reason"]})
            rows = [row]
            common.atomic_csv(args.table_dir / f"{common.setting_stem(identity)}.csv", FIELDS, rows)
        stem = common.setting_stem(setting_config["identity"])
        all_rows.extend(rows)
        config["settings"][stem] = setting_config
        counts[status] = counts.get(status, 0) + 1
        print(f"[{index}/{len(keys)}] {status}: {stem}")
    common.atomic_json(args.config, config)
    common.atomic_csv(args.table_dir / "RP_scan_sidis_binning.csv", FIELDS, all_rows)
    common.atomic_csv(
        args.table_dir / "RP_scan_sidis_binning_Problematic.csv", FIELDS,
        [row for row in all_rows if row.get("Scan_status") in {"WARNING", "ERROR"}],
    )
    print("Status:", ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 2 if counts.get("ERROR") else 0


if __name__ == "__main__":
    raise SystemExit(main())
