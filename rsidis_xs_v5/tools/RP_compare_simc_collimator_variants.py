#!/usr/bin/env python3
"""Compare baseline and SHMS-collimator-on data-to-MC results."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import ROOT


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
SETTING_STEM = "phase1_pass4_PIMINUS_LD2_x0p25Q23p3z0p5thpq2p0"
DEFAULT_BASELINE = (
    PROJECT_DIR
    / "results"
    / "Tables"
    / "RP_data_mc_compare"
    / f"{SETTING_STEM}.csv"
)
DEFAULT_VARIANT = (
    PROJECT_DIR
    / "results"
    / "Tables"
    / "RP_data_mc_compare"
    / "collimator_on"
    / f"{SETTING_STEM}.csv"
)
DEFAULT_TABLE = (
    PROJECT_DIR
    / "results"
    / "Tables"
    / "RP_simc_collimator_compare"
    / f"{SETTING_STEM}.csv"
)
DEFAULT_PDF = (
    PROJECT_DIR
    / "results"
    / "PDFs"
    / "RP_simc_collimator_compare"
    / f"{SETTING_STEM}.pdf"
)
VARIABLES = ("hxp", "hyp", "pxp", "pyp")
VARIABLE_LABELS = {
    "hxp": "HMS x'_{tar}",
    "hyp": "HMS y'_{tar}",
    "pxp": "SHMS x'_{tar}",
    "pyp": "SHMS y'_{tar}",
}
EDGE_LIMITS = {"hxp": 0.15, "hyp": 0.10, "pxp": 0.15, "pyp": 0.10}
IDENTITY_COLUMNS = ("Phase", "Pass", "Run_type", "Target", "Setting")
OUTPUT_COLUMNS = [
    *IDENTITY_COLUMNS,
    "Tier",
    "Variable",
    "Bin_index",
    "Bin_low",
    "Bin_high",
    "Is_edge",
    "Data",
    "Data_err",
    "MC_baseline",
    "MC_baseline_err",
    "MC_collimator_on",
    "MC_collimator_on_err",
    "Data_by_MC_baseline",
    "Data_by_MC_baseline_err",
    "Data_by_MC_collimator_on",
    "Data_by_MC_collimator_on_err",
    "MC_collimator_on_by_baseline",
    "MC_collimator_on_by_baseline_err",
    "Chi2_baseline_whole",
    "NDF_baseline_whole",
    "Chi2_collimator_on_whole",
    "NDF_collimator_on_whole",
    "Chi2_baseline_edge",
    "NDF_baseline_edge",
    "Chi2_collimator_on_edge",
    "NDF_collimator_on_edge",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {
            *IDENTITY_COLUMNS,
            "Tier",
            "Variable",
            "Bin_index",
            "Bin_low",
            "Bin_high",
            "Data_final",
            "Data_final_err",
            "MC_total",
            "MC_total_err",
            "MC_complete",
        }
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return list(reader)


def selected_index(
    rows: Sequence[Mapping[str, str]]
) -> Dict[Tuple[str, int], Mapping[str, str]]:
    output = {}
    for row in rows:
        if row.get("Tier") != "delta" or row.get("Variable") not in VARIABLES:
            continue
        if row.get("MC_complete") != "1":
            raise ValueError(
                f"incomplete MC total for {row.get('Variable')} bin {row.get('Bin_index')}"
            )
        key = (row["Variable"], int(row["Bin_index"]))
        if key in output:
            raise ValueError(f"duplicate comparison bin: {key}")
        output[key] = row
    return output


def ratio(
    numerator: float, numerator_error: float,
    denominator: float, denominator_error: float,
) -> Tuple[float | str, float | str]:
    if not all(map(math.isfinite, (
        numerator, numerator_error, denominator, denominator_error
    ))) or denominator == 0.0:
        return "", ""
    value = numerator / denominator
    variance = (
        (numerator_error / denominator) ** 2
        + (numerator * denominator_error / denominator ** 2) ** 2
    )
    return value, math.sqrt(max(variance, 0.0))


def chi2(
    rows: Sequence[Mapping[str, object]], mc_name: str, edge_only: bool
) -> Tuple[float, int]:
    total = 0.0
    count = 0
    for row in rows:
        if edge_only and not row["Is_edge"]:
            continue
        data = float(row["Data"])
        data_error = float(row["Data_err"])
        mc = float(row[mc_name])
        mc_error = float(row[f"{mc_name}_err"])
        variance = data_error * data_error + mc_error * mc_error
        if not all(map(math.isfinite, (data, data_error, mc, mc_error, variance))):
            continue
        if variance <= 0.0:
            continue
        total += (data - mc) ** 2 / variance
        count += 1
    return total, count


def build_comparison(
    baseline_rows: Sequence[Mapping[str, str]],
    variant_rows: Sequence[Mapping[str, str]],
) -> List[Dict[str, object]]:
    baseline = selected_index(baseline_rows)
    variant = selected_index(variant_rows)
    if set(baseline) != set(variant):
        raise ValueError("baseline and collimator-on bin inventories differ")
    output: List[Dict[str, object]] = []
    by_variable: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for key in sorted(baseline, key=lambda item: (VARIABLES.index(item[0]), item[1])):
        variable, bin_index = key
        first, second = baseline[key], variant[key]
        for column in (*IDENTITY_COLUMNS, "Bin_low", "Bin_high"):
            if first[column] != second[column]:
                raise ValueError(f"baseline/variant mismatch for {key}: {column}")
        for column in ("Data_final", "Data_final_err"):
            if not math.isclose(
                float(first[column]), float(second[column]),
                rel_tol=1e-12, abs_tol=1e-15,
            ):
                raise ValueError(f"data changed between variants for {key}: {column}")
        data = float(first["Data_final"])
        data_error = float(first["Data_final_err"])
        baseline_mc = float(first["MC_total"])
        baseline_error = float(first["MC_total_err"])
        variant_mc = float(second["MC_total"])
        variant_error = float(second["MC_total_err"])
        low, high = float(first["Bin_low"]), float(first["Bin_high"])
        center = 0.5 * (low + high)
        data_baseline = ratio(data, data_error, baseline_mc, baseline_error)
        data_variant = ratio(data, data_error, variant_mc, variant_error)
        variant_baseline = ratio(
            variant_mc, variant_error, baseline_mc, baseline_error
        )
        row: Dict[str, object] = {
            **{column: first[column] for column in IDENTITY_COLUMNS},
            "Tier": "delta",
            "Variable": variable,
            "Bin_index": bin_index,
            "Bin_low": low,
            "Bin_high": high,
            "Is_edge": int(abs(center) >= EDGE_LIMITS[variable]),
            "Data": data,
            "Data_err": data_error,
            "MC_baseline": baseline_mc,
            "MC_baseline_err": baseline_error,
            "MC_collimator_on": variant_mc,
            "MC_collimator_on_err": variant_error,
            "Data_by_MC_baseline": data_baseline[0],
            "Data_by_MC_baseline_err": data_baseline[1],
            "Data_by_MC_collimator_on": data_variant[0],
            "Data_by_MC_collimator_on_err": data_variant[1],
            "MC_collimator_on_by_baseline": variant_baseline[0],
            "MC_collimator_on_by_baseline_err": variant_baseline[1],
        }
        output.append(row)
        by_variable[variable].append(row)
    for variable, variable_rows in by_variable.items():
        metrics = {}
        for label, mc_name in (
            ("baseline", "MC_baseline"),
            ("collimator_on", "MC_collimator_on"),
        ):
            for region, edge_only in (("whole", False), ("edge", True)):
                value, ndf = chi2(variable_rows, mc_name, edge_only)
                metrics[f"Chi2_{label}_{region}"] = value
                metrics[f"NDF_{label}_{region}"] = ndf
        for row in variable_rows:
            row.update(metrics)
    return output


def atomic_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(
                {column: row.get(column, "") for column in OUTPUT_COLUMNS}
                for row in rows
            )
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def graph(
    rows: Sequence[Mapping[str, object]],
    value_name: str,
    error_name: str,
    color: int,
    marker: int,
) -> ROOT.TGraphErrors:
    output = ROOT.TGraphErrors(len(rows))
    for index, row in enumerate(rows):
        x = 0.5 * (float(row["Bin_low"]) + float(row["Bin_high"]))
        value = row[value_name]
        error = row[error_name]
        output.SetPoint(index, x, float(value) if finite(value) else math.nan)
        output.SetPointError(index, 0.0, float(error) if finite(error) else 0.0)
    output.SetMarkerStyle(marker)
    output.SetMarkerSize(1.2)
    output.SetMarkerColor(color)
    output.SetLineColor(color)
    return output


def frame_range(
    rows: Sequence[Mapping[str, object]],
    series: Sequence[Tuple[str, str]],
) -> Tuple[float, float]:
    values = []
    for row in rows:
        for value_name, error_name in series:
            if finite(row[value_name]) and finite(row[error_name]):
                value, error = float(row[value_name]), float(row[error_name])
                values.extend((value - error, value + error))
    if not values:
        return 0.0, 1.0
    low, high = min(0.0, min(values)), max(values)
    span = max(high - low, 1e-12)
    return low - 0.08 * span, high + 0.16 * span


def draw_pdf(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    baseline_path: Path,
    variant_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("c_collimator_ab", "SIMC collimator A/B", 1800, 1200)
    canvas.Print(path.as_posix() + "[")
    canvas.Clear()
    text = ROOT.TPaveText(0.05, 0.06, 0.95, 0.94, "NDC")
    text.SetFillStyle(0); text.SetBorderSize(0); text.SetTextFont(42)
    text.SetTextAlign(12); text.SetTextSize(0.028)
    text.AddText("SIMC SHMS collimator-on (A/B) study")
    text.AddText(SETTING_STEM)
    def display_path(source: Path) -> str:
        try:
            return source.resolve().relative_to(PROJECT_DIR.parent).as_posix()
        except ValueError:
            return source.as_posix()

    baseline_display = Path(display_path(baseline_path))
    variant_display = Path(display_path(variant_path))
    text.AddText(f"Baseline CSV directory: {baseline_display.parent.as_posix()}/")
    text.AddText(f"  {baseline_display.name}")
    text.AddText(f"Collimator-on CSV directory: {variant_display.parent.as_posix()}/")
    text.AddText(f"  {variant_display.name}")
    text.AddText("Baseline: default using_SHMScoll=0")
    text.AddText("Variant: using_HMScoll=0, using_SHMScoll=1")
    text.AddText("Tier: Delta-only")
    text.AddText("Edge: |xptar| >= 0.15 or |yptar| >= 0.10")
    #text.AddText("#chi^{2} uses #sigma_{data}^{2} + #sigma_{MC}^{2}; NDF is the number of valid bins.")
    #text.AddText("The ordinary geometric aperture checks are active in both samples.")
    text.Draw()
    canvas.Print(path.as_posix())

    by_variable = defaultdict(list)
    for row in rows:
        by_variable[str(row["Variable"])].append(row)
    for variable in VARIABLES:
        variable_rows = sorted(
            by_variable[variable], key=lambda row: int(row["Bin_index"])
        )
        canvas.Clear()
        upper = ROOT.TPad(f"upper_{variable}", "", 0.0, 0.48, 1.0, 1.0)
        middle = ROOT.TPad(f"middle_{variable}", "", 0.0, 0.24, 1.0, 0.48)
        lower = ROOT.TPad(f"lower_{variable}", "", 0.0, 0.0, 1.0, 0.24)
        for pad in (upper, middle, lower):
            pad.SetLeftMargin(0.09); pad.SetRightMargin(0.03)
            pad.SetGrid(1, 1); pad.Draw()
        low = float(variable_rows[0]["Bin_low"])
        high = float(variable_rows[-1]["Bin_high"])
        upper.cd()
        frame = ROOT.TH1D(f"frame_upper_{variable}", "", 20, low, high)
        ymin, ymax = frame_range(variable_rows, (
            ("Data", "Data_err"),
            ("MC_baseline", "MC_baseline_err"),
            ("MC_collimator_on", "MC_collimator_on_err"),
        ))
        frame.SetMinimum(ymin); frame.SetMaximum(ymax)
        frame.GetYaxis().SetTitle("Normalized Yield")
        frame.GetYaxis().SetTitleSize(0.050)
        frame.GetYaxis().SetLabelSize(0.040)
        frame.GetXaxis().SetLabelSize(0)
        frame.SetTitle(f"Delta-only: {VARIABLE_LABELS[variable]}")
        frame.Draw()
        data_graph = graph(variable_rows, "Data", "Data_err", ROOT.kBlack, 20)
        baseline_graph = graph(
            variable_rows, "MC_baseline", "MC_baseline_err", ROOT.kBlue + 1, 24
        )
        variant_graph = graph(
            variable_rows, "MC_collimator_on", "MC_collimator_on_err",
            ROOT.kRed + 1, 25,
        )
        for item in (data_graph, baseline_graph, variant_graph):
            item.Draw("P SAME")
        legend = ROOT.TLegend(0.72, 0.72, 0.96, 0.91)
        legend.SetFillStyle(0); legend.SetBorderSize(0)
        legend.AddEntry(data_graph, "Corrected data", "lep")
        legend.AddEntry(baseline_graph, "Baseline MC total", "lep")
        legend.AddEntry(variant_graph, "SHMS collimator-on MC total", "lep")
        legend.Draw()

        middle.cd()
        middle_frame = ROOT.TH1D(f"frame_middle_{variable}", "", 20, low, high)
        middle_frame.SetMinimum(0.0); middle_frame.SetMaximum(2.0)
        middle_frame.GetYaxis().SetTitle("Data / MC")
        middle_frame.GetYaxis().SetTitleSize(0.085)
        middle_frame.GetYaxis().SetLabelSize(0.070)
        middle_frame.GetYaxis().SetTitleOffset(0.48)
        middle_frame.GetXaxis().SetLabelSize(0)
        middle_frame.Draw()
        data_baseline = graph(
            variable_rows,
            "Data_by_MC_baseline",
            "Data_by_MC_baseline_err",
            ROOT.kBlue + 1,
            24,
        )
        data_variant = graph(
            variable_rows,
            "Data_by_MC_collimator_on",
            "Data_by_MC_collimator_on_err",
            ROOT.kRed + 1,
            25,
        )
        data_baseline.Draw("P SAME"); data_variant.Draw("P SAME")

        lower.cd()
        lower_frame = ROOT.TH1D(f"frame_lower_{variable}", "", 20, low, high)
        lower_frame.SetMinimum(0.5); lower_frame.SetMaximum(1.5)
        lower_frame.GetYaxis().SetTitle("Collimator-on / baseline")
        lower_frame.GetYaxis().SetTitleSize(0.082)
        lower_frame.GetYaxis().SetLabelSize(0.068)
        lower_frame.GetYaxis().SetTitleOffset(0.50)
        lower_frame.GetXaxis().SetTitle(VARIABLE_LABELS[variable] + " (rad)")
        lower_frame.GetXaxis().SetTitleSize(0.085)
        lower_frame.GetXaxis().SetLabelSize(0.070)
        lower_frame.Draw()
        variant_baseline = graph(
            variable_rows,
            "MC_collimator_on_by_baseline",
            "MC_collimator_on_by_baseline_err",
            ROOT.kBlack,
            20,
        )
        variant_baseline.Draw("P SAME")
        metrics = variable_rows[0]
        note = ROOT.TPaveText(0.11, 0.70, 0.60, 0.92, "NDC")
        note.SetFillStyle(0); note.SetBorderSize(0); note.SetTextFont(42)
        note.SetTextAlign(12); note.SetTextSize(0.075)
        note.AddText(
            "Whole #chi^{2}/NDF: baseline "
            f"{float(metrics['Chi2_baseline_whole']):.2f}/"
            f"{int(metrics['NDF_baseline_whole'])}, on "
            f"{float(metrics['Chi2_collimator_on_whole']):.2f}/"
            f"{int(metrics['NDF_collimator_on_whole'])}"
        )
        note.AddText(
            "Edge #chi^{2}/NDF: baseline "
            f"{float(metrics['Chi2_baseline_edge']):.2f}/"
            f"{int(metrics['NDF_baseline_edge'])}, on "
            f"{float(metrics['Chi2_collimator_on_edge']):.2f}/"
            f"{int(metrics['NDF_collimator_on_edge'])}"
        )
        note.Draw()
        canvas._rp_objects = (
            upper, middle, lower, frame, middle_frame, lower_frame,
            data_graph, baseline_graph, variant_graph, data_baseline,
            data_variant, variant_baseline, legend, note,
        )
        canvas.Print(path.as_posix())
    canvas.Print(path.as_posix() + "]")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--collimator-on", type=Path, default=DEFAULT_VARIANT)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--output-pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    baseline = args.baseline.expanduser().resolve()
    variant = args.collimator_on.expanduser().resolve()
    rows = build_comparison(read_rows(baseline), read_rows(variant))
    output_csv = args.output_csv.expanduser().resolve()
    atomic_csv(output_csv, rows)
    print(f"A/B CSV: {output_csv}")
    if not args.no_pdf:
        output_pdf = args.output_pdf.expanduser().resolve()
        draw_pdf(output_pdf, rows, baseline, variant)
        print(f"A/B PDF: {output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
