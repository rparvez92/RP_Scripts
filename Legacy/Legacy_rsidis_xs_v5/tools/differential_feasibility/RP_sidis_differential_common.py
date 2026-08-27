#!/usr/bin/env python3
"""Shared Delta-acceptance machinery for SIDIS differential projections."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Dict, Iterable, Mapping, Sequence

import ROOT

import RP_data_mc_compare as dmc
from RP_build_sidis_model import calculator_identity, run_calculator

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

VARIABLES = ("z", "pt2", "phipq")
FINE_BINNING = {
    "z": (20, 0.0, 1.0, "z", ""),
    "pt2": (20, 0.0, 0.6, "p_{T}^{2}", "GeV^{2}"),
    "phipq": (24, 0.0, 2.0 * math.pi, "#phi_{pq}", "rad"),
}
REACTIONS = dmc.REACTIONS
IDENTITY_FIELDS = ("Phase", "Pass", "Run_type", "Target", "Setting")
BINNING_VERSION = "RP_SIDIS_DIFF_V1"
MC_REL_TARGET = 0.10
DATA_REL_TARGET = 0.30
PURITY_TARGET = 0.50
STABILITY_TARGET = 0.50


def locate_bin(value: float, edges: Sequence[float], periodic: bool = False) -> int | None:
    """Return a half-open bin index; the periodic upper edge maps to the first bin."""
    if not finite(value) or len(edges) < 2:
        return None
    low, high = float(edges[0]), float(edges[-1])
    if periodic:
        width = high - low
        if width <= 0: return None
        value = low + ((float(value) - low) % width)
    if value < low or value >= high:
        return None
    for index in range(len(edges) - 1):
        if edges[index] <= value < edges[index + 1]: return index
    return None


def accumulate_nd(events: Iterable[Sequence[float]], weights: Iterable[float],
                  edge_sets: Sequence[Sequence[float]], periodic_axes: Sequence[int] = ()):
    """Sparse event-level N-D accumulator shared by future z×pT²×phi extraction."""
    sums: dict[tuple[int, ...], list[float]] = {}
    flow = 0.0
    for values, weight in zip(events, weights):
        if not finite(weight):
            continue
        key = tuple(locate_bin(float(value), edges, axis in periodic_axes)
                    for axis, (value, edges) in enumerate(zip(values, edge_sets)))
        if any(index is None for index in key):
            flow += float(weight); continue
        cell = sums.setdefault(key, [0.0, 0.0])  # type: ignore[arg-type]
        cell[0] += float(weight); cell[1] += float(weight) ** 2
    return sums, flow


def project_nd(cells: Mapping[tuple[int, ...], Sequence[float]], axis: int):
    result: dict[int, list[float]] = {}
    for key, values in cells.items():
        cell = result.setdefault(key[axis], [0.0, 0.0])
        cell[0] += float(values[0]); cell[1] += float(values[1])
    return result


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def atomic_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(handle, "w") as stream:
            json.dump(json_safe(payload), stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def json_safe(value: object):
    """Convert non-finite diagnostics to JSON null without changing CSV output."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def catalog_index(path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    rows = read_csv(path)
    return {
        tuple(row.get(field, "") for field in IDENTITY_FIELDS) + (row.get("Reaction", ""),): row
        for row in rows
    }


def model_index(path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    return {
        tuple(row.get(field, "") for field in IDENTITY_FIELDS): row
        for row in read_csv(path)
    }


def setting_stem(identity: Mapping[str, str]) -> str:
    return "_".join(identity[field] for field in IDENTITY_FIELDS)


def discover_keys(norm_dir: Path, input_path: Path | None, process_all: bool):
    groups = dmc.discover_norm_files(norm_dir)
    if input_path:
        item = dmc.parse_norm_stem(input_path)
        keys = [(item["Phase"], item["Pass"], item["Run_type"], item["Target"], item["Setting"])]
    elif process_all:
        keys = sorted(
            (*base[:3], target, base[3])
            for base, files in groups.items()
            for target in ("LH2", "LD2")
            if (target, "Elec") in files
        )
    else:
        raise ValueError("provide one normalized-yield CSV or --all")
    return groups, keys


def identity_from_key(key: Sequence[str]) -> dict[str, str]:
    return dict(zip(IDENTITY_FIELDS, key))


def required_data_files(files: Mapping[tuple[str, str], Path], target: str):
    required = {
        "target_e": (target, "Elec"), "target_p": (target, "Pos"),
        "dummy_e": ("DUMMY", "Elec"), "dummy_p": ("DUMMY", "Pos"),
    }
    missing = [name for name, category in required.items() if category not in files]
    return required, missing


def prepare_histograms(
    key: tuple[str, str, str, str, str],
    files: Mapping[tuple[str, str], Path],
    catalog: Mapping[tuple[str, ...], dict[str, str]],
    binning: Mapping[str, tuple[int, float, float, str, str]],
    t7_root: Path,
):
    phase, pass_name, run_type, target, setting = key
    identity = identity_from_key(key)
    required, missing_data = required_data_files(files, target)
    reasons: list[str] = []
    if missing_data:
        return None, "PENDING", ["MISSING_DATA_CLASSES=" + "|".join(missing_data)]
    simc_rows: dict[str, dict[str, str]] = {}
    missing_mc, structural_mc = [], []
    for reaction in REACTIONS:
        row = catalog.get((*key, reaction))
        if row is None or row.get("QA_status") in {"PENDING", "SKIPPED"}:
            missing_mc.append(reaction)
        elif row.get("QA_status") == "ERROR":
            structural_mc.append(reaction)
        else:
            simc_rows[reaction] = row
    if structural_mc:
        return None, "ERROR", ["STRUCTURAL_SIMC_ERROR=" + "|".join(structural_mc)]
    if missing_mc:
        return None, "PENDING", ["INCOMPLETE_MC=" + "|".join(missing_mc)]

    old = dict(dmc.BINNING)
    dmc.BINNING.update(binning)
    try:
        classes = {}
        closures = []
        for name, category in required.items():
            histograms, class_reasons, closure = dmc.class_histograms(
                files[category], t7_root, "delta", tuple(binning),
                f"{setting_stem(identity)}_{name}_diff",
            )
            classes[name] = histograms
            reasons.extend(class_reasons)
            if finite(closure):
                closures.append(float(closure))
        data = {
            variable: dmc.combine_data(classes, variable, dmc.S_DUMMY[(phase, target)])
            for variable in binning
        }
        components = {}
        for reaction, row in simc_rows.items():
            histograms, component_reasons = dmc.simc_histograms(
                row, "delta", tuple(binning), f"{setting_stem(identity)}_{reaction}_diff"
            )
            components[reaction] = histograms
            reasons.extend(component_reasons)
        totals = {
            variable: dmc.add_mc_total(components, variable) for variable in binning
        }
    finally:
        dmc.BINNING.clear()
        dmc.BINNING.update(old)
    status = "WARNING" if reasons else "OK"
    return {
        "identity": identity, "data": data, "components": components,
        "totals": totals, "simc_rows": simc_rows,
        "data_closure_max": max(closures) if closures else math.nan,
    }, status, sorted(set(reasons))


def normalized_phi_expression(branch: str) -> str:
    return f"({branch}<0 ? {branch}+2*TMath::Pi() : ({branch}>=2*TMath::Pi() ? {branch}-2*TMath::Pi() : {branch}))"


def response_and_moments(catalog_row: Mapping[str, str], variable: str,
                         edges: Sequence[float], prefix: str):
    path = Path(catalog_row.get("Recon_root_path", ""))
    if not path.is_file():
        raise FileNotFoundError(path)
    original = {"z": "z", "pt2": "pt2", "phipq": normalized_phi_expression("phipq")}[variable]
    reconstructed = {
        "z": "z_recon", "pt2": "pt2_recon",
        "phipq": normalized_phi_expression("phipq_recon"),
    }[variable]
    edge_array = ROOT.std.vector("double")()
    for edge in edges:
        edge_array.push_back(float(edge))
    nbin = len(edges) - 1
    frame = ROOT.RDataFrame("h10", str(path)).Filter(dmc.SIMC_DELTA)
    frame = frame.Define("_rp_diff_orig", original).Define("_rp_diff_recon", reconstructed)
    frame = frame.Filter(
        "std::isfinite((double)_rp_diff_orig) && _rp_diff_orig>-998 && "
        "std::isfinite((double)_rp_diff_recon) && _rp_diff_recon>-998 && "
        "std::isfinite((double)fWeight)"
    ).Define("_rp_diff_absw", "abs((double)fWeight)")
    response_ptr = frame.Histo2D(
        (f"{prefix}_{variable}_response", "", nbin, edge_array.data(), nbin, edge_array.data()),
        "_rp_diff_orig", "_rp_diff_recon", "_rp_diff_absw",
    )
    response = response_ptr.GetValue().Clone(f"{prefix}_{variable}_response_owned")
    response.SetDirectory(0); response.SetStats(False)
    return response


def response_metrics(response: ROOT.TH2, bin_index: int) -> tuple[float, float]:
    diagonal = response.GetBinContent(bin_index, bin_index)
    reconstructed_total = sum(
        response.GetBinContent(ix, bin_index) for ix in range(1, response.GetNbinsX() + 1)
    )
    generated_total = sum(
        response.GetBinContent(bin_index, iy) for iy in range(1, response.GetNbinsY() + 1)
    )
    purity = diagonal / reconstructed_total if reconstructed_total > 0 else math.nan
    stability = diagonal / generated_total if generated_total > 0 else math.nan
    return purity, stability


def response_bias_rms(response: ROOT.TH2, reconstructed_bin: int) -> tuple[float, float]:
    recon_center = response.GetYaxis().GetBinCenter(reconstructed_bin)
    weight_sum = first = second = 0.0
    for original_bin in range(1, response.GetNbinsX() + 1):
        weight = response.GetBinContent(original_bin, reconstructed_bin)
        residual = recon_center - response.GetXaxis().GetBinCenter(original_bin)
        weight_sum += weight
        first += weight * residual
        second += weight * residual * residual
    if weight_sum <= 0:
        return math.nan, math.nan
    mean = first / weight_sum
    return mean, math.sqrt(max(0.0, second / weight_sum - mean * mean))


def coarse_response_metrics(response: ROOT.TH2, start: int, stop: int) -> tuple[float, float]:
    diagonal = sum(
        response.GetBinContent(ix, iy)
        for ix in range(start, stop + 1) for iy in range(start, stop + 1)
    )
    reconstructed_total = sum(
        response.GetBinContent(ix, iy)
        for ix in range(1, response.GetNbinsX() + 1) for iy in range(start, stop + 1)
    )
    generated_total = sum(
        response.GetBinContent(ix, iy)
        for ix in range(start, stop + 1) for iy in range(1, response.GetNbinsY() + 1)
    )
    purity = diagonal / reconstructed_total if reconstructed_total > 0 else math.nan
    stability = diagonal / generated_total if generated_total > 0 else math.nan
    return purity, stability


def coarse_response_bias_rms(response: ROOT.TH2, start: int, stop: int) -> tuple[float, float]:
    """Weighted reconstructed-minus-original residual for a reconstructed coarse bin."""
    weight_sum = first = second = 0.0
    for ix in range(1, response.GetNbinsX() + 1):
        for iy in range(start, stop + 1):
            weight = response.GetBinContent(ix, iy)
            residual = response.GetYaxis().GetBinCenter(iy) - response.GetXaxis().GetBinCenter(ix)
            weight_sum += weight; first += weight * residual; second += weight * residual * residual
    if weight_sum <= 0:
        return math.nan, math.nan
    mean = first / weight_sum
    return mean, math.sqrt(max(0.0, second / weight_sum - mean * mean))


def hist_range(hist: ROOT.TH1, start: int, stop: int) -> tuple[float, float]:
    value = sum(hist.GetBinContent(index) for index in range(start, stop + 1))
    error = math.sqrt(sum(hist.GetBinError(index) ** 2 for index in range(start, stop + 1)))
    return value, error


def candidate_quality(data: ROOT.TH1, mc: ROOT.TH1, response: ROOT.TH2,
                      start: int, stop: int) -> dict[str, object]:
    data_yield, data_error = hist_range(data, start, stop)
    mc_yield, mc_error = hist_range(mc, start, stop)
    purity, stability = coarse_response_metrics(response, start, stop)
    data_rel = abs(data_error / data_yield) if data_yield != 0 else math.nan
    mc_rel = abs(mc_error / mc_yield) if mc_yield != 0 else math.inf
    passes = (
        math.isfinite(mc_yield) and mc_yield != 0 and math.isfinite(mc_rel)
        and mc_rel <= MC_REL_TARGET
        and (not math.isfinite(data_rel) or data_rel <= DATA_REL_TARGET)
        and math.isfinite(purity) and purity >= PURITY_TARGET
        and math.isfinite(stability) and stability >= STABILITY_TARGET
    )
    return {
        "data_yield": data_yield, "data_error": data_error, "data_rel": data_rel,
        "mc_yield": mc_yield, "mc_error": mc_error, "mc_rel": mc_rel,
        "purity": purity, "stability": stability, "passes": passes,
    }


def recommend_adjacent_edges(data: ROOT.TH1, mc: ROOT.TH1, response: ROOT.TH2) -> tuple[list[float], list[dict[str, object]]]:
    nbin = data.GetNbinsX()
    groups: list[tuple[int, int]] = []
    start = 1
    while start <= nbin:
        chosen = nbin
        for stop in range(start, nbin + 1):
            if candidate_quality(data, mc, response, start, stop)["passes"]:
                chosen = stop
                break
        groups.append((start, chosen))
        start = chosen + 1
    if len(groups) > 1:
        last_quality = candidate_quality(data, mc, response, *groups[-1])
        if not last_quality["passes"]:
            previous = groups[-2]
            groups[-2:] = [(previous[0], groups[-1][1])]
    edges = [data.GetXaxis().GetBinLowEdge(groups[0][0])]
    details = []
    for start, stop in groups:
        edges.append(data.GetXaxis().GetBinUpEdge(stop))
        details.append({"start": start, "stop": stop,
                        **candidate_quality(data, mc, response, start, stop)})
    return edges, details


def recommend_phi_edges(data: ROOT.TH1, mc: ROOT.TH1, response: ROOT.TH2) -> tuple[list[float], list[dict[str, object]]]:
    fine_n = data.GetNbinsX()
    for coarse_n in (24, 12, 8, 6):
        if fine_n % coarse_n:
            continue
        width = fine_n // coarse_n
        groups = [(1 + index * width, (index + 1) * width) for index in range(coarse_n)]
        details = [candidate_quality(data, mc, response, *group) | {"start": group[0], "stop": group[1]}
                   for group in groups]
        if all(detail["passes"] for detail in details) or coarse_n == 6:
            low, high = data.GetXaxis().GetXmin(), data.GetXaxis().GetXmax()
            edges = [low + index * (high - low) / coarse_n for index in range(coarse_n + 1)]
            return edges, details
    raise RuntimeError("no phi binning candidate")


def weighted_bin_moments(catalog_row: Mapping[str, str], variable: str,
                         edges: Sequence[float]) -> list[dict[str, float]]:
    path = Path(catalog_row.get("Recon_root_path", ""))
    if not path.is_file():
        raise FileNotFoundError(path)
    expression = dmc.SIMC_EXPRESSIONS[variable]
    frame = ROOT.RDataFrame("h10", str(path)).Filter(dmc.SIMC_DELTA).Define("_rp_bin_v", expression)
    arrays = frame.AsNumpy(["_rp_bin_v", "fWeight", "xbj_recon", "Q2_recon",
                            "z_recon", "pt2_recon", "sighad"])
    result = []
    for low, high in zip(edges[:-1], edges[1:]):
        sums = {name: 0.0 for name in ("w", "x", "x2", "Q2", "Q22", "z", "z2", "pt", "pt2", "M", "M2")}
        for value, weight, x, q2, z, pt2, sighad in zip(
            arrays["_rp_bin_v"], arrays["fWeight"], arrays["xbj_recon"],
            arrays["Q2_recon"], arrays["z_recon"], arrays["pt2_recon"], arrays["sighad"]
        ):
            if not (math.isfinite(float(value)) and low <= value < high):
                continue
            values = (weight, x, q2, z, pt2, sighad)
            if not all(math.isfinite(float(item)) and float(item) > -998 for item in values):
                continue
            pt = math.sqrt(max(0.0, float(pt2)))
            w = float(weight)
            sums["w"] += w
            for name, item in (("x", x), ("Q2", q2), ("z", z), ("pt", pt), ("M", sighad)):
                item = float(item); sums[name] += w * item; sums[name + "2"] += w * item * item
        row = {"weight_sum": sums["w"]}
        for name in ("x", "Q2", "z", "pt", "M"):
            mean = sums[name] / sums["w"] if sums["w"] != 0 else math.nan
            variance = sums[name + "2"] / sums["w"] - mean * mean if sums["w"] != 0 else math.nan
            row[f"{name}_mean"] = mean
            row[f"{name}_rms"] = math.sqrt(max(0.0, variance)) if math.isfinite(variance) else math.nan
        result.append(row)
    return result


def calculator_model(executable: Path, model_row: Mapping[str, str], moments: Mapping[str, float]):
    target_a, target_z, charge = calculator_identity(model_row["Target"], model_row["Run_type"])
    inputs = (number(model_row["Ebeam_GeV"]), moments["Q2_mean"], moments["x_mean"],
              moments["z_mean"], moments["pt_mean"], target_a, target_z, charge, 1.0)
    phi, sighad, sigsemi = run_calculator(executable, inputs)
    return {"phi": phi, "sighad": sighad, "sigsemi": sigsemi}
