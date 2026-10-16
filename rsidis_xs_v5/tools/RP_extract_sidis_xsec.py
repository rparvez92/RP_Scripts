#!/usr/bin/env python3
"""Extract Pass1 setting-integrated SIDIS model-scaled cross sections."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile

REACTIONS = ("sidis", "rho", "delta", "exclusive")
IDENTITY = ("Phase", "Pass", "Run_type", "Target", "setting_id")
PDF_SIGMA_UNITS = "#mub/GeV^{2}/sr^{2}"
PDF_M_UNITS = "GeV^{-2}"
FIELDS = [
    *IDENTITY, "Acceptance", "Y_Data", "Y_Data_err",
    "Y_MC_sidis", "Y_MC_sidis_err", "Y_MC_rho", "Y_MC_rho_err",
    "Y_MC_delta", "Y_MC_delta_err", "Y_MC_exclusive", "Y_MC_exclusive_err",
    "Y_MC", "Y_MC_err", *(f"Y_MC_{reaction}_fraction" for reaction in REACTIONS),
    "C_Y", "C_Y_err", "sigma_sidis_model", "sigma_sidis_data",
    "sigma_sidis_data_err", "sigma_sidis_units", "M_sidis_model",
    "M_sidis_data", "M_sidis_data_err", "M_sidis_units",
    "Y_Data_input_closure_rel_max",
    *(f"Y_MC_{reaction}_QA_closure_rel" for reaction in REACTIONS),
    "Y_MC_complete", "Source_data_mc_csv", "Source_model_catalog",
    "Extraction_status", "Extraction_reason",
]


def finite(value: object) -> bool:
    try: return math.isfinite(float(value))
    except (TypeError, ValueError): return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def quadrature(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def ratio_with_error(data: float, data_err: float, mc: float,
                     mc_err: float) -> tuple[float, float]:
    if not all(math.isfinite(v) for v in (data, data_err, mc, mc_err)) or mc == 0:
        return math.nan, math.nan
    return data / mc, math.hypot(data_err / mc, data * mc_err / (mc * mc))


def integrate_baseline_rows(rows: list[dict[str, str]]) -> dict[str, object]:
    selected = [row for row in rows if row.get("Variable") == "hdelta"
                and row.get("Bin_index")]
    if not selected:
        return {"metadata": rows[0] if rows else {}, "selected": [],
                "reason": "BASELINE_HDELTA_ROWS_MISSING"}
    result: dict[str, object] = {"metadata": selected[0], "selected": selected,
                                "reason": ""}
    for stem in ("Data_final", *(f"MC_{reaction}" for reaction in REACTIONS)):
        values = [number(row.get(stem)) for row in selected]
        errors = [number(row.get(f"{stem}_err")) for row in selected]
        result[stem] = sum(values) if all(map(math.isfinite, values)) else math.nan
        result[f"{stem}_err"] = quadrature(errors) if all(map(math.isfinite, errors)) else math.nan
    result["Data_closure_rel"] = max(
        (abs(number(row.get("Data_closure_rel"), 0.0)) for row in selected),
        default=math.nan,
    )
    result["MC_complete"] = all(row.get("MC_complete") == "1" for row in selected)
    return result


def load_catalog(path: Path, reaction: bool = False) -> dict[tuple[str, ...], dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return {
        tuple(row.get(field, "") for field in IDENTITY)
        + ((row.get("Reaction", ""),) if reaction else ()): row
        for row in rows
    }


def relative_difference(observed: float, reference: float) -> float:
    if not all(math.isfinite(v) for v in (observed, reference)) or reference == 0:
        return math.nan
    return (observed - reference) / reference


def extract_one(data_path: Path, models: dict[tuple[str, ...], dict[str, str]],
                qa: dict[tuple[str, ...], dict[str, str]], model_path: Path) -> dict[str, object]:
    with data_path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    integrated = integrate_baseline_rows(rows)
    metadata = integrated["metadata"]
    identity = tuple(str(metadata.get(field, "")) for field in IDENTITY)
    output: dict[str, object] = {field: "" for field in FIELDS}
    output.update(dict(zip(IDENTITY, identity)))
    output.update({
        "Acceptance": "Baseline PID + Delta-only acceptance",
        "Source_data_mc_csv": str(data_path.resolve()),
        "Source_model_catalog": str(model_path.resolve()),
    })
    if not integrated["selected"]:
        output.update(Extraction_status="ERROR",
                      Extraction_reason=integrated["reason"])
        return output
    model = models.get(identity)
    if model is None or model.get("Model_status") != "OK":
        model_status = model.get("Model_status", "MISSING") if model else "MISSING"
        output.update(Extraction_status="ERROR",
                      Extraction_reason=f"MODEL_{model_status}")
        return output

    reasons = []
    status = "OK"
    data, data_err = number(integrated["Data_final"]), number(integrated["Data_final_err"])
    values = [number(integrated[f"MC_{reaction}"]) for reaction in REACTIONS]
    errors = [number(integrated[f"MC_{reaction}_err"]) for reaction in REACTIONS]
    mc = sum(values) if all(map(math.isfinite, values)) else math.nan
    mc_err = quadrature(errors) if all(map(math.isfinite, errors)) else math.nan
    output.update(Y_Data=data, Y_Data_err=data_err, Y_MC=mc, Y_MC_err=mc_err,
                  Y_Data_input_closure_rel_max=integrated["Data_closure_rel"],
                  Y_MC_complete=int(bool(integrated["MC_complete"])))
    for reaction, value, error in zip(REACTIONS, values, errors):
        output[f"Y_MC_{reaction}"] = value
        output[f"Y_MC_{reaction}_err"] = error
        output[f"Y_MC_{reaction}_fraction"] = value / mc if finite(mc) and mc != 0 else math.nan
        reference_row = qa.get(identity + (reaction,))
        reference = number(reference_row.get("SimYield")) if reference_row else math.nan
        closure = relative_difference(value, reference)
        output[f"Y_MC_{reaction}_QA_closure_rel"] = closure
        if reference_row and reference_row.get("QA_status") == "ERROR":
            status = "ERROR"; reasons.append(f"Y_MC_{reaction}_QA_ERROR")

    if not integrated["MC_complete"]:
        status = "ERROR"; reasons.append("MC_TOTAL_INCOMPLETE")
    elif not all(map(math.isfinite, (data, data_err, mc, mc_err))):
        status = "ERROR"; reasons.append("NONFINITE_YIELD")
    elif mc == 0:
        status = "ERROR"; reasons.append("ZERO_MC_TOTAL")
    else:
        correction, correction_err = ratio_with_error(data, data_err, mc, mc_err)
        sigma_model = number(model["sigma_sidis_model"])
        multiplicity_model = number(model["M_sidis_model"])
        output.update({
            "C_Y": correction, "C_Y_err": correction_err,
            "sigma_sidis_model": sigma_model,
            "sigma_sidis_data": correction * sigma_model,
            "sigma_sidis_data_err": abs(sigma_model) * correction_err,
            "sigma_sidis_units": model["sigma_sidis_units"],
            "M_sidis_model": multiplicity_model,
            "M_sidis_data": correction * multiplicity_model,
            "M_sidis_data_err": abs(multiplicity_model) * correction_err,
            "M_sidis_units": model["M_sidis_units"],
        })
        if data < 0 or correction < 0:
            status = "WARNING"; reasons.append("NEGATIVE_DATA_OR_C_Y")
        closure = number(integrated["Data_closure_rel"])
        if finite(closure) and abs(closure) > 0.01:
            if status == "OK": status = "WARNING"
            reasons.append(f"DATA_CLOSURE={closure:.6g}")
        for reaction in REACTIONS:
            closure = number(output[f"Y_MC_{reaction}_QA_closure_rel"])
            if finite(closure) and abs(closure) > 1e-5:
                if status == "OK": status = "WARNING"
                reasons.append(f"Y_MC_{reaction}_QA_CLOSURE={closure:.6g}")
    output.update(Extraction_status=status,
                  Extraction_reason=";".join(sorted(set(reasons))))
    return output


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def draw_pdf(path: Path, row: dict[str, object]) -> None:
    import ROOT
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = str(abs(hash(path)))
    canvas = ROOT.TCanvas(f"c_xsec_{suffix}", "SIDIS extraction", 1800, 1200)
    canvas.Print(str(path) + "[")
    text = ROOT.TLatex(); text.SetNDC(); text.SetTextFont(42)
    lines = [
        "RP SIDIS Cross-Section and Multiplicity",
        "/".join(str(row[field]) for field in IDENTITY),
        "Acceptance: Baseline PID + Delta-only acceptance",
        "Y_{MC} = Y_{sidis} + Y_{rho} + Y_{delta} + Y_{exclusive}",
        "C_{Y} = Y_{Data}/Y_{MC}",
        "#sigma_{sidis,data} = C_{Y} #sigma_{sidis,model}",
        "M_{sidis,data} = C_{Y} M_{sidis,model}",
        f"Y_{{Data}} = {number(row['Y_Data']):.6g} #pm {number(row['Y_Data_err']):.3g}",
        f"Y_{{MC}} = {number(row['Y_MC']):.6g} #pm {number(row['Y_MC_err']):.3g}",
        f"C_{{Y}} = {number(row['C_Y']):.6g} #pm {number(row['C_Y_err']):.3g}",
        f"#sigma_{{data}} = {number(row['sigma_sidis_data']):.6g} #pm {number(row['sigma_sidis_data_err']):.3g} {PDF_SIGMA_UNITS}",
        f"M_{{data}} = {number(row['M_sidis_data']):.6g} #pm {number(row['M_sidis_data_err']):.3g} {PDF_M_UNITS}",
        f"Status = {row['Extraction_status']}",
        f"Reason = {row['Extraction_reason'] or 'none'}",
    ]
    y = 0.93
    for index, line in enumerate(lines):
        text.SetTextSize(0.044 if index == 0 else 0.031)
        text.SetTextFont(62 if index == 0 else 42)
        text.DrawLatex(0.07, y, line); y -= 0.061
    canvas.Print(str(path)); canvas.Clear()
    hist = ROOT.TH1D(f"components_{suffix}",
                     "Normalized Yield Composition;Component;Normalized Yield", 6, 0.5, 6.5)
    labels = ("Data", "sidis", "rho", "delta", "exclusive", "MC total")
    values = (row["Y_Data"], *(row[f"Y_MC_{reaction}"] for reaction in REACTIONS), row["Y_MC"])
    errors = (row["Y_Data_err"], *(row[f"Y_MC_{reaction}_err"] for reaction in REACTIONS), row["Y_MC_err"])
    for index, (label, value, error) in enumerate(zip(labels, values, errors), 1):
        hist.GetXaxis().SetBinLabel(index, label)
        if finite(value): hist.SetBinContent(index, float(value))
        if finite(error): hist.SetBinError(index, float(error))
    hist.SetMarkerStyle(20); hist.SetMarkerSize(1.5); hist.SetLineWidth(2)
    hist.Draw("E1"); canvas.SetGridy(); canvas.Print(str(path)); canvas.Clear()
    text.SetTextFont(62); text.SetTextSize(0.044)
    text.DrawLatex(0.07, 0.92, "Integrated Closure and Status")
    text.SetTextFont(42); text.SetTextSize(0.031); y = 0.82
    closure_lines = [f"Data closure maximum = {number(row['Y_Data_input_closure_rel_max']):.6g}"]
    closure_lines += [f"{reaction} MC closure = {number(row[f'Y_MC_{reaction}_QA_closure_rel']):.6g}"
                      for reaction in REACTIONS]
    closure_lines += [f"MC complete = {row['Y_MC_complete']}",
                      f"Status = {row['Extraction_status']}",
                      f"Reason = {row['Extraction_reason'] or 'none'}"]
    for line in closure_lines:
        text.DrawLatex(0.07, y, line); y -= 0.075
    canvas.Print(str(path)); canvas.Print(str(path) + "]")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--data-dir", type=Path, default=root / "results/Pass1/Tables/RP_data_mc_compare")
    parser.add_argument("--model-catalog", type=Path, default=root / "results/Pass1/Tables/RP_sidis_model/RP_sidis_model.csv")
    parser.add_argument("--simc-qa", type=Path, default=root / "results/Pass1/Tables/RP_extract_simc_info/RP_extract_simc_info.csv")
    parser.add_argument("--table-dir", type=Path, default=root / "results/Pass1/Tables/RP_extract_sidis_xsec")
    parser.add_argument("--pdf-dir", type=Path, default=root / "results/Pass1/PDFs/RP_extract_sidis_xsec")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    if args.all == bool(args.input): parser.error("provide exactly one input CSV or --all")
    paths = sorted(args.data_dir.glob("phase*.csv")) if args.all else [args.input]
    models, qa = load_catalog(args.model_catalog), load_catalog(args.simc_qa, True)
    results = []
    for source in paths:
        result = extract_one(source, models, qa, args.model_catalog)
        atomic_csv(args.table_dir / source.name, [result])
        if not args.no_pdf: draw_pdf(args.pdf_dir / f"{source.stem}.pdf", result)
        results.append(result); print(f"{source.stem}: {result['Extraction_status']}")
    return 1 if any(row["Extraction_status"] == "ERROR" for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
