#!/usr/bin/env python3
"""Extract setting-integrated SIDIS model-scaled results using Delta acceptance."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile

REACTIONS = ("sidis", "rho", "delta", "exclusive")
IDENTITY = ("Phase", "Pass", "Run_type", "Target", "Setting")
PDF_SIGMA_UNITS = "#mub/GeV^{2}/sr^{2}"
PDF_M_UNITS = "GeV^{-2}"

FIELDS = [
    *IDENTITY, "Acceptance", "Y_Data", "Y_Data_err",
    "Y_MC_sidis", "Y_MC_sidis_err", "Y_MC_rho", "Y_MC_rho_err",
    "Y_MC_delta", "Y_MC_delta_err", "Y_MC_exclusive", "Y_MC_exclusive_err",
    "Y_MC", "Y_MC_err", "Y_MC_sidis_fraction", "Y_MC_rho_fraction",
    "Y_MC_delta_fraction", "Y_MC_exclusive_fraction", "C_Y", "C_Y_err",
    "sigma_sidis_model", "sigma_sidis_data", "sigma_sidis_data_err",
    "sigma_sidis_units", "M_sidis_model", "M_sidis_data", "M_sidis_data_err",
    "M_sidis_units", "Y_Data_input_closure_rel_max",
    "Y_MC_sidis_QA_closure_rel", "Y_MC_rho_QA_closure_rel",
    "Y_MC_delta_QA_closure_rel", "Y_MC_exclusive_QA_closure_rel", "Y_MC_complete",
    "Source_data_mc_csv", "Source_model_catalog", "Extraction_status",
    "Extraction_reason",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def quadrature(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def ratio_with_error(data: float, data_err: float, mc: float,
                     mc_err: float) -> tuple[float, float]:
    if not all(math.isfinite(v) for v in (data, data_err, mc, mc_err)) or mc == 0:
        return math.nan, math.nan
    ratio = data / mc
    error = math.hypot(data_err / mc, data * mc_err / (mc * mc))
    return ratio, error


def relative_difference(observed: float, reference: float) -> float:
    if not all(math.isfinite(v) for v in (observed, reference)) or reference == 0:
        return math.nan
    return (observed - reference) / reference


def integrate_delta_rows(rows: list[dict[str, str]]) -> dict[str, object]:
    selected = [row for row in rows if row.get("Tier") == "delta"
                and row.get("Variable") == "hdelta" and row.get("Bin_index")]
    if not selected:
        metadata = rows[0] if rows else {}
        return {"metadata": metadata, "selected": [], "reason": "DELTA_HDELTA_ROWS_MISSING"}
    values: dict[str, object] = {"metadata": selected[0], "selected": selected, "reason": ""}
    for stem in ("Data_final", *(f"MC_{reaction}" for reaction in REACTIONS), "MC_total"):
        entries = [number(row.get(stem)) for row in selected]
        errors = [number(row.get(f"{stem}_err")) for row in selected]
        values[stem] = sum(entries) if all(math.isfinite(v) for v in entries) else math.nan
        values[f"{stem}_err"] = quadrature(errors) if all(math.isfinite(v) for v in errors) else math.nan
    values["Data_closure_rel"] = max(
        (abs(number(row.get("Data_closure_rel"), 0.0)) for row in selected), default=math.nan
    )
    values["MC_complete"] = all(row.get("MC_complete") == "1" for row in selected)
    return values


def load_model_catalog(path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {tuple(row.get(name, "") for name in IDENTITY): row for row in rows}


def load_qa_catalog(path: Path | None) -> dict[tuple[str, ...], dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        tuple(row.get(name, "") for name in IDENTITY) + (row.get("Reaction", ""),): row
        for row in rows
    }


def extract_one(data_path: Path, model_rows: dict[tuple[str, ...], dict[str, str]],
                qa_rows: dict[tuple[str, ...], dict[str, str]], model_path: Path) -> dict[str, object]:
    with data_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    integrated = integrate_delta_rows(rows)
    metadata = integrated.get("metadata", {})
    identity = tuple(str(metadata.get(name, "")) for name in IDENTITY)
    output: dict[str, object] = {field: "" for field in FIELDS}
    output.update(dict(zip(IDENTITY, identity)))
    output.update({"Acceptance": "Delta-only",
                   "Source_data_mc_csv": str(data_path.resolve()),
                   "Source_model_catalog": str(model_path.resolve())})
    reasons: list[str] = []
    status = "OK"

    model = model_rows.get(identity)
    if model is not None and model.get("Model_status") == "SKIPPED":
        output.update({"Extraction_status": "SKIPPED",
                       "Extraction_reason": f"MODEL_SKIPPED={model.get('Model_reason','')}"})
        return output

    upstream_status = str(metadata.get("Status", ""))
    if not integrated.get("selected"):
        if upstream_status == "SKIPPED":
            status = "SKIPPED"
        elif upstream_status in {"PENDING", ""}:
            status = "PENDING"
        else:
            status = "ERROR"
        reasons.append(str(integrated.get("reason", "DELTA_ROWS_MISSING")))
        output.update({"Extraction_status": status, "Extraction_reason": ";".join(reasons)})
        return output

    if model is None:
        output.update({"Extraction_status": "PENDING", "Extraction_reason": "MODEL_ROW_MISSING"})
        return output
    model_status = model.get("Model_status", "")
    if model_status != "OK":
        output.update({"Extraction_status": model_status or "PENDING",
                       "Extraction_reason": f"MODEL_{model_status or 'MISSING'}={model.get('Model_reason','')}"})
        return output

    data = number(integrated["Data_final"])
    data_err = number(integrated["Data_final_err"])
    component_values = [number(integrated[f"MC_{reaction}"]) for reaction in REACTIONS]
    component_errors = [number(integrated[f"MC_{reaction}_err"]) for reaction in REACTIONS]
    mc_total = sum(component_values) if all(math.isfinite(v) for v in component_values) else math.nan
    mc_total_err = quadrature(component_errors) if all(math.isfinite(v) for v in component_errors) else math.nan

    output.update({"Y_Data": data, "Y_Data_err": data_err,
                   "Y_MC": mc_total, "Y_MC_err": mc_total_err,
                   "Y_Data_input_closure_rel_max": integrated["Data_closure_rel"],
                   "Y_MC_complete": int(bool(integrated["MC_complete"]))})
    for reaction, value, error in zip(REACTIONS, component_values, component_errors):
        output[f"Y_MC_{reaction}"] = value
        output[f"Y_MC_{reaction}_err"] = error
        output[f"Y_MC_{reaction}_fraction"] = value / mc_total if math.isfinite(mc_total) and mc_total != 0 else math.nan
        qa = qa_rows.get(identity + (reaction,))
        reference = number(qa.get("SimYield_delta")) if qa else math.nan
        closure = relative_difference(value, reference)
        output[f"Y_MC_{reaction}_QA_closure_rel"] = closure
        if qa and qa.get("QA_status") == "ERROR":
            status = "ERROR"
            reasons.append(f"Y_MC_{reaction}_QA_ERROR")

    if not integrated["MC_complete"]:
        status = "PENDING"
        reasons.append("MC_TOTAL_INCOMPLETE")
    elif not all(math.isfinite(v) for v in (data, data_err, mc_total, mc_total_err)):
        status = "ERROR"
        reasons.append("NONFINITE_YIELD")
    elif mc_total == 0:
        status = "ERROR"
        reasons.append("ZERO_MC_TOTAL")
    else:
        cy, cy_err = ratio_with_error(data, data_err, mc_total, mc_total_err)
        sigma_model = number(model.get("sigma_sidis_model"))
        multiplicity_model = number(model.get("M_sidis_model"))
        output.update({
            "C_Y": cy, "C_Y_err": cy_err,
            "sigma_sidis_model": sigma_model,
            "sigma_sidis_data": cy * sigma_model,
            "sigma_sidis_data_err": abs(sigma_model) * cy_err,
            "sigma_sidis_units": model.get("sigma_sidis_units", ""),
            "M_sidis_model": multiplicity_model,
            "M_sidis_data": cy * multiplicity_model,
            "M_sidis_data_err": abs(multiplicity_model) * cy_err,
            "M_sidis_units": model.get("M_sidis_units", ""),
        })
        if data < 0 or cy < 0:
            status = "WARNING"
            reasons.append("NEGATIVE_DATA_OR_C_Y")
        data_closure = number(integrated["Data_closure_rel"])
        if math.isfinite(data_closure) and abs(data_closure) > 0.01:
            status = "WARNING" if status == "OK" else status
            reasons.append(f"DATA_CLOSURE={data_closure:.6g}")
        for reaction in REACTIONS:
            closure = number(output[f"Y_MC_{reaction}_QA_closure_rel"])
            if math.isfinite(closure) and abs(closure) > 1e-5:
                status = "WARNING" if status == "OK" else status
                reasons.append(f"Y_MC_{reaction}_QA_CLOSURE={closure:.6g}")

    output.update({"Extraction_status": status, "Extraction_reason": ";".join(reasons)})
    return output


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def draw_pdf(path: Path, row: dict[str, object]) -> None:
    import ROOT  # type: ignore
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    path.parent.mkdir(parents=True, exist_ok=True)
    object_suffix = str(abs(hash(str(path))))
    canvas = ROOT.TCanvas(f"c_sidis_xsec_{object_suffix}", "sidis extraction", 1800, 1000)
    canvas.SetMargin(0.07, 0.04, 0.08, 0.06)
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextFont(42)
    text.SetTextSize(0.031)
    lines = [
        "Delta-Acceptance sidis Cross-Section and Multiplicity",
        f"{row['Phase']}, {row['Pass']}, {row['Run_type']}, {row['Target']}, {row['Setting']}",
        "Data cuts: PID + HMS/SHMS delta acceptance",
        "MC cuts: HMS/SHMS delta acceptance",
        "Y_{MC} = Y_{MC,sidis} + Y_{MC,rho} + Y_{MC,delta} + Y_{MC,exclusive}",
        "C_{Y} = Y_{Data} / Y_{MC}",
        "Cross-section: #sigma_{sidis,data} = C_{Y} #sigma_{sidis,model}",
        "Multiplicity: M_{sidis,data} = C_{Y} M_{sidis,model}",
        "M_{sidis,model} is the calc_semi_xsec sighad output",
        f"Data normalized yield = {number(row['Y_Data']):.6g} +/- {number(row['Y_Data_err']):.3g}",
        f"MC normalized yield = {number(row['Y_MC']):.6g} +/- {number(row['Y_MC_err']):.3g}",
        f"C_Y = {number(row['C_Y']):.6g} +/- {number(row['C_Y_err']):.3g}",
        f"#sigma_{{sidis,model}} = {number(row['sigma_sidis_model']):.6g} {PDF_SIGMA_UNITS}",
        f"#sigma_{{sidis,data}} = {number(row['sigma_sidis_data']):.6g} +/- {number(row['sigma_sidis_data_err']):.3g} {PDF_SIGMA_UNITS}",
        f"M_{{sidis,model}} = {number(row['M_sidis_model']):.6g} {PDF_M_UNITS}",
        f"M_{{sidis,data}} = {number(row['M_sidis_data']):.6g} +/- {number(row['M_sidis_data_err']):.3g} {PDF_M_UNITS}",
        f"Status = {row['Extraction_status']}",
        f"Reason = {row['Extraction_reason'] or 'none'}",
    ]
    y = 0.94
    for index, line in enumerate(lines):
        if index == 0:
            text.SetTextSize(0.042)
        elif index == 1:
            text.SetTextSize(0.033)
        else:
            text.SetTextSize(0.031)
        text.DrawLatex(0.08, y, line)
        y -= 0.051
    canvas.Print(str(path) + "(")

    canvas.Clear()
    frame = ROOT.TH1D(f"h_components_{object_suffix}",
                      "Normalized Yield Composition;Component;Normalized Yield", 6, 0.5, 6.5)
    frame.SetStats(False)
    labels = ("Data", "sidis", "rho", "delta", "exclusive", "MC")
    values = (row["Y_Data"], row["Y_MC_sidis"], row["Y_MC_rho"],
              row["Y_MC_delta"], row["Y_MC_exclusive"], row["Y_MC"])
    errors = (row["Y_Data_err"], row["Y_MC_sidis_err"], row["Y_MC_rho_err"],
              row["Y_MC_delta_err"], row["Y_MC_exclusive_err"], row["Y_MC_err"])
    for index, (label, value, error) in enumerate(zip(labels, values, errors), 1):
        frame.GetXaxis().SetBinLabel(index, label)
        if finite(value):
            frame.SetBinContent(index, float(value))
        if finite(error):
            frame.SetBinError(index, float(error))
    frame.SetMarkerStyle(20)
    frame.SetMarkerSize(1.5)
    frame.SetLineWidth(2)
    ymin = min((number(v) - number(e, 0.0) for v, e in zip(values, errors) if finite(v)), default=0.0)
    ymax = max((number(v) + number(e, 0.0) for v, e in zip(values, errors) if finite(v)), default=1.0)
    span = max(ymax - ymin, abs(ymax) * 0.2, 1e-12)
    frame.SetMinimum(min(0.0, ymin - 0.15 * span))
    frame.SetMaximum(ymax + 0.2 * span)
    frame.Draw("E1")
    canvas.SetGridy(True)
    canvas.Print(str(path))

    canvas.Clear()
    text.SetTextSize(0.038)
    text.DrawLatex(0.08, 0.91, "Integrated Closure and Final Status")
    y = 0.82
    closure_lines = [
        "Data input closure: max_{class} |Y_{event histogram}-Y_{saved}|/|Y_{saved}|",
        f"Data input closure relative max = {number(row['Y_Data_input_closure_rel_max']):.6g}",
        "MC QA closure: (Y_{DMC,reaction}-Y_{QA,reaction})/Y_{QA,reaction}",
    ]
    closure_lines.extend(
        f"{reaction} MC QA closure relative = {number(row[f'Y_MC_{reaction}_QA_closure_rel']):.6g}"
        for reaction in REACTIONS
    )
    closure_lines.extend((f"MC complete = {row['Y_MC_complete']}",
                          f"Final status = {row['Extraction_status']}",
                          f"Reason = {row['Extraction_reason'] or 'none'}"))
    text.SetTextSize(0.033)
    for line in closure_lines:
        text.DrawLatex(0.08, y, line)
        y -= 0.07
    canvas.Print(str(path) + ")")


def main() -> int:
    analysis_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="one RP_data_mc_compare CSV")
    parser.add_argument("--all", action="store_true", help="process all setting CSVs")
    parser.add_argument("--data-dir", type=Path,
                        default=analysis_root / "results/Tables/RP_data_mc_compare")
    parser.add_argument("--model-catalog", type=Path,
                        default=analysis_root / "results/Tables/RP_sidis_model/RP_sidis_model.csv")
    parser.add_argument("--simc-qa", type=Path,
                        default=analysis_root / "results/Tables/RP_extract_simc_info/RP_extract_simc_info.csv")
    parser.add_argument("--table-dir", type=Path,
                        default=analysis_root / "results/Tables/RP_extract_sidis_xsec")
    parser.add_argument("--pdf-dir", type=Path,
                        default=analysis_root / "results/PDFs/RP_extract_sidis_xsec")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    if args.all == bool(args.input):
        parser.error("provide exactly one input CSV or --all")
    paths = sorted(args.data_dir.glob("phase*.csv")) if args.all else [args.input]
    model_rows = load_model_catalog(args.model_catalog)
    qa_rows = load_qa_catalog(args.simc_qa)
    results = []
    for source in paths:
        result = extract_one(source, model_rows, qa_rows, args.model_catalog)
        stem = source.stem
        atomic_csv(args.table_dir / f"{stem}.csv", [result])
        if not args.no_pdf:
            draw_pdf(args.pdf_dir / f"{stem}.pdf", result)
        results.append(result)
        print(f"{stem}: {result['Extraction_status']}")
    return 1 if any(row["Extraction_status"] == "ERROR" for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
