#!/usr/bin/env python3
"""Extract SIMC/recon_hcana normalization metadata and QA diagnostics."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import ROOT


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_SIMC_DIR = PROJECT_DIR.parents[1] / "simc_gfortran"
DEFAULT_MANIFEST = DEFAULT_SIMC_DIR / "infiles" / "RP_Simc" / "simc_input_manifest.csv"
DEFAULT_T7_ROOT = Path("/Volumes/T7/RSIDIS")
DEFAULT_TABLE_DIR = PROJECT_DIR / "results" / "Tables" / "RP_extract_simc_info"
DEFAULT_PDF_DIR = PROJECT_DIR / "results" / "PDFs" / "RP_extract_simc_info"

RUNNABLE_STATUSES = {"GENERATED", "EXISTING_FILE_NOT_OVERWRITTEN"}
REACTIONS = ("sidis", "rho", "delta", "exclusive")
PHASES = ("phase1", "phase2")

DELTA_CUT = (
    "hsdelta > -8.0 && hsdelta < 8.0 && "
    "ssdelta > -10.0 && ssdelta < 22.0"
)
FULL_CUT = (
    DELTA_CUT
    + " && hsxptar > -0.15 && hsxptar < 0.15"
    + " && hsyptar > -0.10 && hsyptar < 0.10"
    + " && ssxptar > -0.15 && ssxptar < 0.15"
    + " && ssyptar > -0.10 && ssyptar < 0.10"
)

REQUIRED_RAW_BRANCHES = {
    "Weight", "hsdelta", "hsxptar", "hsyptar",
    "ssdelta", "ssxptar", "ssyptar",
}
REQUIRED_RECON_BRANCHES = REQUIRED_RAW_BRANCHES | {"Ngen", "normfac", "fWeight"}

KINEMATICS = (
    ("xbj", "xbj", "xbj_recon", False),
    ("Q2", "Q2", "Q2_recon", False),
    ("W", "W", "W_recon", False),
    ("z", "z", "z_recon", False),
    ("thetapq", "thetapq", "thetapq_recon", False),
    ("phipq", "phipq", "phipq_recon", True),
    ("pt2", "pt2", "pt2_recon", False),
    ("missmass", "missmass", "missmass_recon", False),
)

IDENTITY_COLUMNS = [
    "Phase", "Pass", "Reaction", "Run_type", "Target", "Setting",
    "Simulation_variant", "using_HMScoll", "using_SHMScoll",
    "Parent_input", "Random_initialization",
    "Manifest_Ngen", "Generation_status", "Generation_reason",
    "Source_input", "Source_leaf_csv", "Source_template",
]
FILE_COLUMNS = [
    "Hist_path", "Raw_root_path", "Recon_root_path",
    "Hist_exists", "Raw_root_exists", "Recon_root_exists",
    "Hist_size_MB", "Raw_root_size_MB", "Recon_root_size_MB",
    "Raw_tree", "Recon_tree", "Raw_entries", "Recon_entries",
]
NORMALIZATION_COLUMNS = [
    "Hist_Ngen", "Recon_Ngen",
    "Manifest_Ngen_matches_observed",
    "Hist_normfac", "Recon_normfac",
    "Normfac_over_Ngen", "FWeight_identity_max_rel_diff",
]
WEIGHT_COLUMNS = [
    "Weight_sum_nocut", "Weight_sum2_nocut", "Weight_min", "Weight_max",
    "FWeight_zero_count", "FWeight_negative_count", "FWeight_nonfinite_count",
]
TIER_COLUMNS: List[str] = []
for _tier in ("nocut", "delta", "full"):
    TIER_COLUMNS.extend([
        f"N_selected_{_tier}", f"Acceptance_fraction_{_tier}",
        f"Sum_fWeight_{_tier}", f"Sum_fWeight2_{_tier}",
        f"SimYield_{_tier}", f"SimYield_err_{_tier}",
        f"SimYield_rel_err_{_tier}", f"Neff_{_tier}",
    ])
KINEMATIC_COLUMNS: List[str] = [
    "Kinematic_branch_availability", "Kinematic_value_provenance"
]
for _label, _, _, _ in KINEMATICS:
    for _tier in ("delta", "full"):
        KINEMATIC_COLUMNS.extend([
            f"{_label}_{_tier}_simc_weighted_mean",
            f"{_label}_{_tier}_simc_weighted_rms",
            f"{_label}_{_tier}_recon_weighted_mean",
            f"{_label}_{_tier}_recon_weighted_rms",
            f"{_label}_{_tier}_residual_weighted_mean",
            f"{_label}_{_tier}_residual_weighted_rms",
        ])
QA_COLUMNS = [
    "File_flag", "Tree_flag", "Entry_flag", "Ngen_flag", "Normfac_flag",
    "FWeight_flag", "Weight_flag", "MC_precision_flag",
    "QA_status", "QA_reason",
]
OUTPUT_COLUMNS = (
    IDENTITY_COLUMNS + FILE_COLUMNS + NORMALIZATION_COLUMNS + WEIGHT_COLUMNS
    + TIER_COLUMNS + KINEMATIC_COLUMNS + QA_COLUMNS
)

NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def empty_output_row() -> Dict[str, object]:
    return {column: "" for column in OUTPUT_COLUMNS}


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {
            "Phase", "Pass", "Run_type", "Target", "Setting", "Reaction",
            "Ngen", "Output_file", "Generation_status",
        }
        if reader.fieldnames is None:
            raise ValueError(f"manifest has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"manifest is missing columns: {sorted(missing)}")
        return list(reader)


def parse_hist(path: Path) -> Tuple[int | None, float | None]:
    ngen = None
    normfac = None
    ngen_pattern = re.compile(rf"^\s*Ngen(?:\s*\(request\))?\s*=\s*(\d+)", re.I)
    norm_pattern = re.compile(rf"\bnormfac\s*=\s*({NUMBER_RE})", re.I)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if ngen is None:
                match = ngen_pattern.search(line)
                if match:
                    ngen = int(match.group(1))
            if normfac is None:
                match = norm_pattern.search(line)
                if match:
                    normfac = float(match.group(1).replace("D", "E").replace("d", "e"))
    return ngen, normfac


def phase_directory(phase: str) -> str:
    normalized = phase.lower()
    if normalized not in PHASES:
        raise ValueError(f"unsupported phase: {phase}")
    return "Phase" + normalized[-1]


def expected_paths(
    t7_root: Path, manifest_row: Mapping[str, str]
) -> Dict[str, Path]:
    source_input = Path(manifest_row["Output_file"])
    stem = source_input.stem
    simulation = t7_root / phase_directory(manifest_row["Phase"]) / "Simulation"
    run_type = manifest_row["Run_type"]
    output_category = source_input.parent.name

    def choose(category: str, filename: str) -> Path:
        candidate_categories = []
        for value in (output_category, "coin", run_type):
            if value and value not in candidate_categories:
                candidate_categories.append(value)
        candidates = tuple(
            simulation / category / value / filename
            for value in candidate_categories
        )
        return next((candidate for candidate in candidates if candidate.is_file()), candidates[0])

    return {
        "hist": choose("outfiles", f"{stem}.hist"),
        "raw_root": choose("ROOTfiles", f"{stem}.root"),
        "recon_root": choose("ROOTfiles", f"recon_hcana_{stem}.root"),
    }


def manifest_identity(source: Mapping[str, str]) -> Dict[str, object]:
    row = empty_output_row()
    mappings = {
        "Phase": "Phase", "Pass": "Pass", "Reaction": "Reaction",
        "Run_type": "Run_type", "Target": "Target", "Setting": "Setting",
        "Simulation_variant": "Simulation_variant",
        "using_HMScoll": "using_HMScoll",
        "using_SHMScoll": "using_SHMScoll",
        "Parent_input": "Parent_input",
        "Random_initialization": "Random_initialization",
        "Manifest_Ngen": "Ngen", "Generation_status": "Generation_status",
        "Generation_reason": "Diagnostic_reason", "Source_input": "Output_file",
        "Source_leaf_csv": "Source_leaf_csv", "Source_template": "Source_template",
    }
    for output, original in mappings.items():
        row[output] = source.get(original, "")
    return row


def branch_names(tree: ROOT.TTree) -> set[str]:
    return {branch.GetName() for branch in tree.GetListOfBranches()}


def relative_difference(first: float, second: float, floor: float = 1e-30) -> float:
    return abs(first - second) / max(abs(first), abs(second), floor)


def result_value(result: object) -> object:
    return result.GetValue()


def weighted_stats(
    sumw: float, sumwx: float, sumwx2: float
) -> Tuple[float | None, float | None]:
    if not all(math.isfinite(value) for value in (sumw, sumwx, sumwx2)):
        return None, None
    if abs(sumw) <= 0.0:
        return None, None
    mean = sumwx / sumw
    variance = sumwx2 / sumw - mean * mean
    tolerance = 1e-12 * max(1.0, abs(sumwx2 / sumw), mean * mean)
    if variance < 0.0 and abs(variance) <= tolerance:
        variance = 0.0
    if variance < 0.0:
        return mean, None
    return mean, math.sqrt(variance)


def add_tier_metrics(
    row: MutableMapping[str, object],
    tier: str,
    count: int,
    total_entries: int,
    sumw: float,
    sumw2: float,
) -> None:
    row[f"N_selected_{tier}"] = count
    row[f"Acceptance_fraction_{tier}"] = (
        count / total_entries if total_entries > 0 else ""
    )
    row[f"Sum_fWeight_{tier}"] = sumw
    row[f"Sum_fWeight2_{tier}"] = sumw2
    row[f"SimYield_{tier}"] = sumw
    error = math.sqrt(sumw2) if math.isfinite(sumw2) and sumw2 >= 0.0 else math.nan
    row[f"SimYield_err_{tier}"] = error
    row[f"SimYield_rel_err_{tier}"] = (
        error / abs(sumw)
        if math.isfinite(error) and math.isfinite(sumw) and sumw != 0.0
        else ""
    )
    row[f"Neff_{tier}"] = (
        sumw * sumw / sumw2
        if math.isfinite(sumw) and math.isfinite(sumw2) and sumw2 > 0.0
        else ""
    )


def scan_recon_tree(
    tree: ROOT.TTree,
    available: set[str],
    row: MutableMapping[str, object],
    reaction: str,
) -> None:
    dataframe = (
        ROOT.RDataFrame(tree)
        .Define("_fw", "static_cast<double>(fWeight)")
        .Define("_fw2", "_fw * _fw")
        .Define("_raw_weight", "static_cast<double>(Weight)")
        .Define("_raw_weight2", "_raw_weight * _raw_weight")
        .Define(
            "_fw_expected",
            "Ngen > 0 ? static_cast<double>(Weight) * normfac / "
            "static_cast<double>(Ngen) : std::numeric_limits<double>::quiet_NaN()",
        )
        .Define(
            "_fw_identity_rel",
            "std::abs(_fw - _fw_expected) / "
            "std::max({std::abs(_fw), std::abs(_fw_expected), 1e-30})",
        )
    )
    kinematics = list(KINEMATICS)
    provenance = {label: "generator_native" for label, _, _, _ in kinematics}
    if reaction.lower() in {"delta", "exclusive"}:
        derived_dependencies = {
            "xbj": {"Q2", "nu", "Q2_recon", "nu_recon"},
            "z": {"phad", "nu", "nu_recon"},
            "pt2": {"phad", "thetapq", "thetapq_recon"},
        }
        expressions = {
            "xbj": (
                "(std::isfinite(static_cast<double>(Q2)) && Q2>-998 && "
                "std::isfinite(static_cast<double>(nu)) && nu>-998 && "
                "std::abs(static_cast<double>(nu))>1e-12) ? "
                "static_cast<double>(Q2)/(2.0*0.9382720813*static_cast<double>(nu)) : "
                "std::numeric_limits<double>::quiet_NaN()",
                "(std::isfinite(static_cast<double>(Q2_recon)) && Q2_recon>-998 && "
                "std::isfinite(static_cast<double>(nu_recon)) && nu_recon>-998 && "
                "std::abs(static_cast<double>(nu_recon))>1e-12) ? "
                "static_cast<double>(Q2_recon)/(2.0*0.9382720813*static_cast<double>(nu_recon)) : "
                "std::numeric_limits<double>::quiet_NaN()",
            ),
            "z": (
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(nu)) && nu>-998 && "
                "std::abs(static_cast<double>(nu))>1e-12) ? "
                "static_cast<double>(phad)/static_cast<double>(nu) : "
                "std::numeric_limits<double>::quiet_NaN()",
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(nu_recon)) && nu_recon>-998 && "
                "std::abs(static_cast<double>(nu_recon))>1e-12) ? "
                "static_cast<double>(phad)/static_cast<double>(nu_recon) : "
                "std::numeric_limits<double>::quiet_NaN()",
            ),
            "pt2": (
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(thetapq)) && thetapq>-998) ? "
                "std::pow(static_cast<double>(phad)*std::sin(static_cast<double>(thetapq)),2) : "
                "std::numeric_limits<double>::quiet_NaN()",
                "(std::isfinite(static_cast<double>(phad)) && phad>-998 && "
                "std::isfinite(static_cast<double>(thetapq_recon)) && thetapq_recon>-998) ? "
                "std::pow(static_cast<double>(phad)*std::sin(static_cast<double>(thetapq_recon)),2) : "
                "std::numeric_limits<double>::quiet_NaN()",
            ),
        }
        replaced = []
        for label, original, recon, circular in kinematics:
            if label not in expressions:
                replaced.append((label, original, recon, circular))
                continue
            if not derived_dependencies[label].issubset(available):
                replaced.append((label, "__missing_derived__", "__missing_derived__", circular))
                provenance[label] = "derived_missing_dependencies"
                continue
            original_alias = f"_derived_{label}_simc"
            recon_alias = f"_derived_{label}_recon"
            dataframe = dataframe.Define(
                original_alias, expressions[label][0]
            ).Define(recon_alias, expressions[label][1])
            available.update({original_alias, recon_alias})
            replaced.append((label, original_alias, recon_alias, circular))
            provenance[label] = "derived_common_kinematics"
        kinematics = replaced

    row["Kinematic_value_provenance"] = ";".join(
        f"{label}:{provenance[label]}" for label, _, _, _ in kinematics
    )
    delta_frame = dataframe.Filter(DELTA_CUT)
    full_frame = dataframe.Filter(FULL_CUT)

    actions: Dict[str, object] = {
        "count": dataframe.Count(),
        "ngen_min": dataframe.Min("Ngen"),
        "ngen_max": dataframe.Max("Ngen"),
        "normfac_min": dataframe.Min("normfac"),
        "normfac_max": dataframe.Max("normfac"),
        "weight_sum": dataframe.Sum("_raw_weight"),
        "weight_sum2": dataframe.Sum("_raw_weight2"),
        "weight_min": dataframe.Min("_raw_weight"),
        "weight_max": dataframe.Max("_raw_weight"),
        "fw_zero": dataframe.Filter("_fw == 0.0").Count(),
        "fw_negative": dataframe.Filter("_fw < 0.0").Count(),
        "fw_nonfinite": dataframe.Filter("!std::isfinite(_fw)").Count(),
        "identity_max": dataframe.Max("_fw_identity_rel"),
    }
    tier_frames = {"nocut": dataframe, "delta": delta_frame, "full": full_frame}
    for tier, frame in tier_frames.items():
        actions[f"{tier}_count"] = frame.Count()
        actions[f"{tier}_sumw"] = frame.Sum("_fw")
        actions[f"{tier}_sumw2"] = frame.Sum("_fw2")

    kinematic_actions: Dict[Tuple[str, str], Tuple[str, ...]] = {}
    availability_actions: Dict[str, str] = {}
    availability_states: Dict[str, str] = {}
    for index, (label, original, recon, circular) in enumerate(kinematics):
        if original not in available or recon not in available:
            availability_states[label] = "missing"
            continue
        difference = (
            f"std::atan2(std::sin(static_cast<double>({recon})"
            f"-static_cast<double>({original})),"
            f"std::cos(static_cast<double>({recon})"
            f"-static_cast<double>({original})))"
            if circular
            else f"static_cast<double>({recon})-static_cast<double>({original})"
        )
        valid = (
            f"std::isfinite(static_cast<double>({original})) && "
            f"std::isfinite(static_cast<double>({recon})) && "
            f"static_cast<double>({original}) > -998.0 && "
            f"static_cast<double>({recon}) > -998.0 && std::isfinite(_fw)"
        )
        availability_action = f"_kin_available_{index}"
        actions[availability_action] = dataframe.Filter(valid).Count()
        availability_actions[label] = availability_action
        for tier, tier_frame in (("delta", delta_frame), ("full", full_frame)):
            prefix = f"_kin_{tier}_{index}"
            defined = (
                tier_frame
                .Define(f"{prefix}_orig", f"static_cast<double>({original})")
                .Define(f"{prefix}_recon", f"static_cast<double>({recon})")
                .Define(f"{prefix}_diff", difference)
            )
            valid_frame = defined.Filter(valid)
            names = (
                f"{prefix}_count", f"{prefix}_sumw",
                f"{prefix}_orig1", f"{prefix}_orig2",
                f"{prefix}_recon1", f"{prefix}_recon2",
                f"{prefix}_diff1", f"{prefix}_diff2",
            )
            actions[names[0]] = valid_frame.Count()
            actions[names[1]] = valid_frame.Sum("_fw")
            actions[names[2]] = valid_frame.Define(
                f"{prefix}_wo", f"_fw*{prefix}_orig"
            ).Sum(f"{prefix}_wo")
            actions[names[3]] = valid_frame.Define(
                f"{prefix}_wo2", f"_fw*{prefix}_orig*{prefix}_orig"
            ).Sum(f"{prefix}_wo2")
            actions[names[4]] = valid_frame.Define(
                f"{prefix}_wr", f"_fw*{prefix}_recon"
            ).Sum(f"{prefix}_wr")
            actions[names[5]] = valid_frame.Define(
                f"{prefix}_wr2", f"_fw*{prefix}_recon*{prefix}_recon"
            ).Sum(f"{prefix}_wr2")
            actions[names[6]] = valid_frame.Define(
                f"{prefix}_wd", f"_fw*{prefix}_diff"
            ).Sum(f"{prefix}_wd")
            actions[names[7]] = valid_frame.Define(
                f"{prefix}_wd2", f"_fw*{prefix}_diff*{prefix}_diff"
            ).Sum(f"{prefix}_wd2")
            kinematic_actions[(tier, label)] = names

    values = {name: result_value(action) for name, action in actions.items()}
    total_entries = int(values["count"])
    row["_Recon_Ngen_min"] = int(values["ngen_min"])
    row["_Recon_Ngen_max"] = int(values["ngen_max"])
    row["_Recon_normfac_min"] = float(values["normfac_min"])
    row["_Recon_normfac_max"] = float(values["normfac_max"])
    row["Weight_sum_nocut"] = float(values["weight_sum"])
    row["Weight_sum2_nocut"] = float(values["weight_sum2"])
    row["Weight_min"] = float(values["weight_min"])
    row["Weight_max"] = float(values["weight_max"])
    row["FWeight_zero_count"] = int(values["fw_zero"])
    row["FWeight_negative_count"] = int(values["fw_negative"])
    row["FWeight_nonfinite_count"] = int(values["fw_nonfinite"])
    row["FWeight_identity_max_rel_diff"] = float(values["identity_max"])

    for tier in tier_frames:
        add_tier_metrics(
            row, tier, int(values[f"{tier}_count"]), total_entries,
            float(values[f"{tier}_sumw"]), float(values[f"{tier}_sumw2"]),
        )

    for label, action in availability_actions.items():
        availability_states[label] = (
            "available" if int(values[action]) > 0 else "sentinel_only"
        )
    row["Kinematic_branch_availability"] = ";".join(
        f"{label}:{availability_states[label]}"
        for label, _, _, _ in kinematics
    )
    for label, state in availability_states.items():
        row[f"_Kinematic_availability_{label}"] = state
    for (tier, label), names in kinematic_actions.items():
        row[f"_Kinematic_valid_count_{tier}_{label}"] = int(values[names[0]])
        sumw = float(values[names[1]])
        original_stats = weighted_stats(
            sumw, float(values[names[2]]), float(values[names[3]])
        )
        recon_stats = weighted_stats(
            sumw, float(values[names[4]]), float(values[names[5]])
        )
        residual_stats = weighted_stats(
            sumw, float(values[names[6]]), float(values[names[7]])
        )
        for suffix, stats in (
            ("simc", original_stats), ("recon", recon_stats),
            ("residual", residual_stats),
        ):
            row[f"{label}_{tier}_{suffix}_weighted_mean"] = (
                "" if stats[0] is None else stats[0]
            )
            row[f"{label}_{tier}_{suffix}_weighted_rms"] = (
                "" if stats[1] is None else stats[1]
            )


def inspect_manifest_row(
    source: Mapping[str, str], t7_root: Path
) -> Dict[str, object]:
    row = manifest_identity(source)
    generation_status = source.get("Generation_status", "")
    if generation_status not in RUNNABLE_STATUSES:
        row.update({
            "File_flag": "SKIPPED", "QA_status": "SKIPPED",
            "QA_reason": source.get("Diagnostic_reason", generation_status),
        })
        return row

    paths = expected_paths(t7_root, source)
    for key, column in (
        ("hist", "Hist_path"), ("raw_root", "Raw_root_path"),
        ("recon_root", "Recon_root_path"),
    ):
        path = paths[key]
        row[column] = path.as_posix()
        exists_column = {
            "hist": "Hist_exists", "raw_root": "Raw_root_exists",
            "recon_root": "Recon_root_exists",
        }[key]
        size_column = {
            "hist": "Hist_size_MB", "raw_root": "Raw_root_size_MB",
            "recon_root": "Recon_root_size_MB",
        }[key]
        row[exists_column] = int(path.is_file())
        row[size_column] = path.stat().st_size / 1e6 if path.is_file() else ""

    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        row.update({
            "File_flag": "PENDING", "QA_status": "PENDING",
            "QA_reason": "missing=" + ",".join(missing),
        })
        return row
    row["File_flag"] = "OK"

    errors: List[str] = []
    warnings: List[str] = []
    raw_file = ROOT.TFile.Open(paths["raw_root"].as_posix(), "READ")
    recon_file = ROOT.TFile.Open(paths["recon_root"].as_posix(), "READ")
    raw_tree = None if not raw_file or raw_file.IsZombie() else raw_file.Get("h10")
    recon_tree = None if not recon_file or recon_file.IsZombie() else recon_file.Get("h10")
    if not raw_tree:
        errors.append("RAW_H10_UNREADABLE")
    if not recon_tree:
        errors.append("RECON_H10_UNREADABLE")
    if errors:
        row["Tree_flag"] = "ERROR"
        row["QA_status"] = "ERROR"
        row["QA_reason"] = ";".join(errors)
        if raw_file:
            raw_file.Close()
        if recon_file:
            recon_file.Close()
        return row

    row["Tree_flag"] = "OK"
    row["Raw_tree"] = "h10"
    row["Recon_tree"] = "h10"
    raw_entries = int(raw_tree.GetEntries())
    recon_entries = int(recon_tree.GetEntries())
    row["Raw_entries"] = raw_entries
    row["Recon_entries"] = recon_entries
    raw_branches = branch_names(raw_tree)
    recon_branches = branch_names(recon_tree)
    missing_raw = sorted(REQUIRED_RAW_BRANCHES - raw_branches)
    missing_recon = sorted(REQUIRED_RECON_BRANCHES - recon_branches)
    if missing_raw:
        errors.append("RAW_MISSING_BRANCHES=" + "|".join(missing_raw))
    if missing_recon:
        errors.append("RECON_MISSING_BRANCHES=" + "|".join(missing_recon))

    try:
        hist_ngen, hist_normfac = parse_hist(paths["hist"])
    except OSError as error:
        hist_ngen, hist_normfac = None, None
        errors.append(f"HIST_READ_ERROR={error}")
    row["Hist_Ngen"] = "" if hist_ngen is None else hist_ngen
    row["Hist_normfac"] = "" if hist_normfac is None else hist_normfac
    if hist_ngen is None or hist_normfac is None:
        errors.append("HIST_NORMALIZATION_MISSING")
    elif not math.isfinite(hist_normfac) or hist_normfac <= 0.0:
        errors.append("HIST_NORMFAC_INVALID")

    if not errors and recon_entries > 0:
        try:
            scan_recon_tree(
                recon_tree, recon_branches, row, str(source["Reaction"])
            )
        except Exception as error:
            errors.append(f"RECON_SCAN_FAILED={type(error).__name__}:{error}")
    elif recon_entries <= 0:
        errors.append("RECON_TREE_EMPTY")

    if not errors:
        ngen_min = int(row["_Recon_Ngen_min"])
        ngen_max = int(row["_Recon_Ngen_max"])
        if raw_entries != recon_entries:
            errors.append(
                f"RAW_RECON_ENTRY_MISMATCH:raw={raw_entries},recon={recon_entries}"
            )
            row["Entry_flag"] = "ERROR"
        else:
            row["Entry_flag"] = "OK"

        ngen_constant = ngen_min == ngen_max
        if not ngen_constant:
            errors.append(
                f"RECON_NGEN_NONCONSTANT:min={ngen_min},max={ngen_max}"
            )
            row["Ngen_flag"] = "ERROR"
        else:
            row["Recon_Ngen"] = ngen_min
            observed = [raw_entries, recon_entries, hist_ngen, ngen_min]
            if len(set(observed)) != 1:
                errors.append(
                    "ENTRY_NGEN_MISMATCH=" + "|".join(map(str, observed))
                )
                row["Ngen_flag"] = "ERROR"
            else:
                row["Ngen_flag"] = "OK"

        manifest_ngen = int(source["Ngen"]) if source.get("Ngen", "").isdigit() else None
        row["Manifest_Ngen_matches_observed"] = (
            "" if manifest_ngen is None else int(manifest_ngen == hist_ngen)
        )
        if hist_ngen and hist_normfac is not None:
            row["Normfac_over_Ngen"] = hist_normfac / hist_ngen

        norm_min = float(row["_Recon_normfac_min"])
        norm_max = float(row["_Recon_normfac_max"])
        if relative_difference(norm_min, norm_max) > 1e-12:
            errors.append(
                f"RECON_NORMFAC_NONCONSTANT:min={norm_min:.12g},"
                f"max={norm_max:.12g}"
            )
            row["Normfac_flag"] = "ERROR"
        else:
            row["Recon_normfac"] = norm_min
            if relative_difference(norm_min, float(hist_normfac)) > 1e-12:
                errors.append(
                    f"NORMFAC_MISMATCH:hist={float(hist_normfac):.12g},"
                    f"recon={norm_min:.12g}"
                )
                row["Normfac_flag"] = "ERROR"
            else:
                row["Normfac_flag"] = "OK"

        identity = float(row["FWeight_identity_max_rel_diff"])
        nonfinite = int(row["FWeight_nonfinite_count"])
        if not math.isfinite(identity) or identity > 1e-5:
            errors.append(f"FWEIGHT_IDENTITY_MISMATCH={identity:.6g}")
            row["FWeight_flag"] = "ERROR"
        elif nonfinite > 0:
            errors.append(f"FWEIGHT_NONFINITE={nonfinite}")
            row["FWeight_flag"] = "ERROR"
        else:
            row["FWeight_flag"] = "OK"

        if int(row["FWeight_negative_count"]) > 0:
            warnings.append(f"NEGATIVE_FWEIGHT={row['FWeight_negative_count']}")
            row["Weight_flag"] = "WARNING"
        else:
            row["Weight_flag"] = "OK"

        if int(row["N_selected_full"]) <= 0:
            errors.append("FULL_ACCEPTANCE_EMPTY")
            row["MC_precision_flag"] = "ERROR"
        else:
            relative_error = row["SimYield_rel_err_full"]
            effective = row["Neff_full"]
            precision_warnings = []
            if isinstance(relative_error, (float, int)) and relative_error > 0.05:
                precision_warnings.append(f"FULL_RELERR={relative_error:.6g}")
            if isinstance(effective, (float, int)) and effective < 1000.0:
                precision_warnings.append(f"FULL_NEFF={effective:.6g}")
            if precision_warnings:
                warnings.extend(precision_warnings)
                row["MC_precision_flag"] = "WARNING"
            else:
                row["MC_precision_flag"] = "OK"

    if errors:
        for column in ("Entry_flag", "Ngen_flag", "Normfac_flag",
                       "FWeight_flag", "Weight_flag", "MC_precision_flag"):
            if not row[column]:
                row[column] = "NOT_EVALUATED"
        row["QA_status"] = "ERROR"
        row["QA_reason"] = ";".join(errors + warnings)
    elif warnings:
        row["QA_status"] = "WARNING"
        row["QA_reason"] = ";".join(warnings)
    else:
        row["QA_status"] = "OK"
        row["QA_reason"] = ""

    raw_file.Close()
    recon_file.Close()
    return row


def format_csv_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return format(value, ".12g")
    return value


def write_csv_atomic(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    column: format_csv_value(row.get(column, ""))
                    for column in OUTPUT_COLUMNS
                })
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def configure_canvas(canvas: ROOT.TCanvas) -> None:
    canvas.Clear()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.08)
    # Setting identifiers are intentionally verbose; reserve enough space for
    # vertical labels without crowding the x-axis title.
    canvas.SetBottomMargin(0.28)
    canvas.SetTopMargin(0.22)
    canvas.SetGrid(1, 1)


def draw_text_page(
    canvas: ROOT.TCanvas, pdf: Path, title: str, lines: Sequence[str]
) -> None:
    canvas.Clear()
    canvas.SetGrid(0, 0)
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(13)
    latex.SetTextFont(62)
    latex.SetTextSize(0.038)
    latex.DrawLatex(0.07, 0.94, title)
    latex.SetTextFont(42)
    latex.SetTextSize(0.021)
    y = 0.89
    for line in lines:
        latex.DrawLatex(0.08, y, line)
        y -= 0.033
        if y < 0.06:
            canvas.Print(pdf.as_posix())
            canvas.Clear()
            latex.DrawLatex(0.07, 0.94, title + " (continued)")
            y = 0.89
    canvas.Print(pdf.as_posix())


LONG_CATEGORY_LABELS = {
    ("PIPLUS", "LH2"): "PIPLUS LH2",
    ("PIMINUS", "LH2"): "PIMINUS LH2",
    ("PIPLUS", "LD2"): "PIPLUS LD2",
    ("PIMINUS", "LD2"): "PIMINUS LD2",
}
CATEGORY_COLORS = {
    ("PIPLUS", "LH2"): ROOT.kBlue + 1,
    ("PIMINUS", "LH2"): ROOT.kRed + 1,
    ("PIPLUS", "LD2"): ROOT.kGreen + 2,
    ("PIMINUS", "LD2"): ROOT.kMagenta + 1,
}
CATEGORY_ORDER = tuple(CATEGORY_COLORS)
CATEGORY_OFFSETS = {
    category: -0.30 + 0.20 * index
    for index, category in enumerate(CATEGORY_ORDER)
}


def group_rows(
    rows: Sequence[Mapping[str, object]],
    statuses: Sequence[str] = ("OK", "WARNING"),
) -> Dict[Tuple[str, str, str], List[Mapping[str, object]]]:
    groups: Dict[Tuple[str, str, str], List[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        if row["QA_status"] not in statuses:
            continue
        key = (str(row["Phase"]), str(row["Pass"]), str(row["Reaction"]))
        groups[key].append(row)
    for values in groups.values():
        values.sort(key=lambda row: (
            str(row["Setting"]), str(row["Target"]), str(row["Run_type"])
        ))
    return groups


def setting_names(rows: Sequence[Mapping[str, object]]) -> List[str]:
    return sorted({str(row["Setting"]) for row in rows})


def chunk_rows_by_settings(
    rows: Sequence[Mapping[str, object]],
    max_settings: int = 6,
) -> List[List[Mapping[str, object]]]:
    if max_settings <= 0:
        raise ValueError("max_settings must be positive")
    settings = setting_names(rows)
    chunks: List[List[Mapping[str, object]]] = []
    for start in range(0, len(settings), max_settings):
        selected = set(settings[start:start + max_settings])
        chunks.append([
            row for row in rows if str(row["Setting"]) in selected
        ])
    return chunks


def grouped_x(
    row: Mapping[str, object],
    setting_index: Mapping[str, int],
    series_index: int = 0,
    series_count: int = 1,
) -> float:
    category = (str(row["Run_type"]), str(row["Target"]))
    if category not in CATEGORY_OFFSETS:
        raise ValueError(f"unsupported plot category: {category}")
    series_offset = (
        0.0
        if series_count <= 1
        else (series_index - 0.5 * (series_count - 1)) * 0.060
    )
    return (
        float(setting_index[str(row["Setting"])])
        + CATEGORY_OFFSETS[category]
        + series_offset
    )


def categorical_frame(
    rows: Sequence[Mapping[str, object]], title: str, y_title: str,
    low: float, high: float,
) -> ROOT.TH1D:
    settings = setting_names(rows)
    frame = ROOT.TH1D(
        f"frame_{abs(hash((title, y_title, tuple(settings))))}",
        f"{title};Setting;{y_title}", len(settings), 0.5, len(settings) + 0.5,
    )
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(low)
    frame.SetMaximum(high)
    frame.GetXaxis().SetLabelSize(0.020)
    frame.GetXaxis().SetTitleOffset(3.4)
    frame.GetXaxis().CenterTitle(True)
    for index, setting in enumerate(settings, start=1):
        frame.GetXaxis().SetBinLabel(index, setting)
    frame.LabelsOption("v", "X")
    return frame


def draw_setting_separators(
    settings: Sequence[str], low: float, high: float
) -> List[ROOT.TLine]:
    lines = []
    for index in range(1, len(settings)):
        line = ROOT.TLine(index + 0.5, low, index + 0.5, high)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw()
        lines.append(line)
    return lines


def draw_plot_legends(
    series: Sequence[Tuple[str, int]],
) -> List[object]:
    series_legend = ROOT.TLegend(0.17, 0.80, 0.45, 0.93)
    category_legend = ROOT.TLegend(0.48, 0.80, 0.93, 0.93)
    for legend in (series_legend, category_legend):
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.018)
        legend.SetNColumns(2)

    keepalive: List[object] = [series_legend, category_legend]
    for label, marker in series:
        dummy = ROOT.TGraph()
        dummy.SetMarkerStyle(marker)
        dummy.SetMarkerColor(ROOT.kBlack)
        dummy.SetLineColor(ROOT.kBlack)
        series_legend.AddEntry(dummy, label, "p")
        keepalive.append(dummy)
    for category in CATEGORY_ORDER:
        dummy = ROOT.TGraph()
        dummy.SetMarkerStyle(20)
        dummy.SetMarkerColor(CATEGORY_COLORS[category])
        dummy.SetLineColor(CATEGORY_COLORS[category])
        category_legend.AddEntry(dummy, LONG_CATEGORY_LABELS[category], "p")
        keepalive.append(dummy)
    series_legend.Draw()
    category_legend.Draw()
    return keepalive


def padded_range(
    values: Iterable[float], errors: Iterable[float] | None = None,
    nonnegative: bool = False,
) -> Tuple[float, float]:
    value_list = list(values)
    error_list = list(errors) if errors is not None else [0.0] * len(value_list)
    extents = [
        bound
        for value, error in zip(value_list, error_list)
        if math.isfinite(value) and math.isfinite(error)
        for bound in (value - abs(error), value + abs(error))
    ]
    if not extents:
        return 0.0, 1.0
    low, high = min(extents), max(extents)
    span = high - low
    padding = max(0.12 * span, 0.05 * max(abs(low), abs(high)), 1e-12)
    plotted_low = low - padding
    if nonnegative and low >= 0.0:
        plotted_low = max(0.0, plotted_low)
    return plotted_low, high + padding


def draw_graph_page(
    canvas: ROOT.TCanvas,
    pdf: Path,
    rows: Sequence[Mapping[str, object]],
    title: str,
    y_title: str,
    columns: Sequence[Tuple[str, str, str | None]],
    series_markers: Sequence[int] | None = None,
    nonnegative: bool = False,
    logy: bool = False,
) -> None:
    configure_canvas(canvas)
    values: List[float] = []
    errors: List[float] = []
    for _, value_column, error_column in columns:
        for row in rows:
            value = row.get(value_column, "")
            error = row.get(error_column, 0.0) if error_column else 0.0
            if isinstance(value, (float, int)) and isinstance(error, (float, int)):
                values.append(float(value))
                errors.append(float(error))
    low, high = padded_range(values, errors, nonnegative)
    if logy:
        positive = [
            value - abs(error) for value, error in zip(values, errors)
            if value - abs(error) > 0.0
        ]
        if positive:
            low = max(min(positive) * 0.7, 1e-30)
            high = max(high, max(values) * 1.4)
            canvas.SetLogy(True)

    frame = categorical_frame(rows, title, y_title, low, high)
    frame.Draw()
    settings = setting_names(rows)
    setting_index = {setting: index for index, setting in enumerate(settings, start=1)}
    keepalive: List[object] = [frame]
    keepalive.extend(draw_setting_separators(settings, low, high))
    markers = list(series_markers or (20, 24, 22))
    if len(markers) < len(columns):
        raise ValueError("not enough marker styles for plotted series")
    for series_index, (label, value_column, error_column) in enumerate(columns):
        marker = markers[series_index]
        for category in CATEGORY_ORDER:
            color = CATEGORY_COLORS[category]
            x_values, y_values, y_errors = [], [], []
            for row in rows:
                if (row["Run_type"], row["Target"]) != category:
                    continue
                value = row.get(value_column, "")
                error = row.get(error_column, 0.0) if error_column else 0.0
                if not isinstance(value, (float, int)) or not math.isfinite(float(value)):
                    continue
                x_values.append(grouped_x(
                    row, setting_index, series_index, len(columns)
                ))
                y_values.append(float(value))
                y_errors.append(float(error) if isinstance(error, (float, int)) else 0.0)
            if not x_values:
                continue
            graph = ROOT.TGraphErrors(len(x_values))
            for point, (x, y, error) in enumerate(zip(x_values, y_values, y_errors)):
                graph.SetPoint(point, x, y)
                graph.SetPointError(point, 0.0, error)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetMarkerSize(1.15)
            graph.Draw("P SAME")
            keepalive.append(graph)
    keepalive.extend(draw_plot_legends([
        (label, markers[index]) for index, (label, _, _) in enumerate(columns)
    ]))
    canvas.Print(pdf.as_posix())
    canvas.SetLogy(False)


def draw_dual_panel_page(
    canvas: ROOT.TCanvas,
    pdf: Path,
    rows: Sequence[Mapping[str, object]],
    context: str,
) -> None:
    canvas.Clear()
    canvas.Divide(1, 2)
    panels = (
        (
            "Relative MC uncertainty: Delta-only and Full-cut",
            "Relative uncertainty",
            (
                ("Delta-only", "SimYield_rel_err_delta"),
                ("Full-cut", "SimYield_rel_err_full"),
            ),
        ),
        (
            "Effective MC entries: Delta-only and Full-cut",
            "N_{eff}",
            (("Delta-only", "Neff_delta"), ("Full-cut", "Neff_full")),
        ),
    )
    keepalive = []
    for pad_index, (title, y_title, columns) in enumerate(panels, start=1):
        pad = canvas.cd(pad_index)
        pad.SetGrid(1, 1)
        pad.SetBottomMargin(0.24)
        pad.SetTopMargin(0.25)
        values = [
            float(row[column])
            for _, column in columns
            for row in rows
            if isinstance(row.get(column), (float, int))
        ]
        low, high = padded_range(values, nonnegative=True)
        frame = categorical_frame(rows, f"{context}: {title}", y_title, low, high)
        frame.Draw()
        keepalive.append(frame)
        settings = setting_names(rows)
        setting_index = {
            setting: index for index, setting in enumerate(settings, start=1)
        }
        keepalive.extend(draw_setting_separators(settings, low, high))
        markers = (24, 20)
        for series_index, (_, column) in enumerate(columns):
            for category in CATEGORY_ORDER:
                graph = ROOT.TGraph()
                point = 0
                for row in rows:
                    if (row["Run_type"], row["Target"]) != category:
                        continue
                    value = row.get(column, "")
                    if not isinstance(value, (float, int)):
                        continue
                    graph.SetPoint(
                        point,
                        grouped_x(row, setting_index, series_index, len(columns)),
                        float(value),
                    )
                    point += 1
                if point:
                    graph.SetMarkerStyle(markers[series_index])
                    graph.SetMarkerColor(CATEGORY_COLORS[category])
                    graph.SetMarkerSize(1.1)
                    graph.Draw("P SAME")
                    keepalive.append(graph)
        keepalive.extend(draw_plot_legends([
            (label, markers[index]) for index, (label, _) in enumerate(columns)
        ]))
    canvas.Print(pdf.as_posix())


def kinematic_availability(
    row: Mapping[str, object], label: str
) -> str:
    internal = row.get(f"_Kinematic_availability_{label}")
    if isinstance(internal, str) and internal:
        return internal
    encoded = str(row.get("Kinematic_branch_availability", ""))
    for item in encoded.split(";"):
        name, separator, state = item.partition(":")
        if separator and name == label:
            return state
    return "missing"


def unavailable_kinematic_message(
    rows: Sequence[Mapping[str, object]],
    label: str,
    tier: str,
) -> str:
    states = {kinematic_availability(row, label) for row in rows}
    if states == {"missing"}:
        return "Unavailable: branch missing"
    if "available" not in states and "sentinel_only" in states:
        return "Unavailable: all values = -999"
    valid_count = sum(
        int(row.get(f"_Kinematic_valid_count_{tier}_{label}", 0))
        for row in rows
    )
    return (
        "No valid values after selection"
        if valid_count == 0
        else "Weighted statistics unavailable"
    )


def draw_kinematic_page(
    canvas: ROOT.TCanvas,
    pdf: Path,
    rows: Sequence[Mapping[str, object]],
    context: str,
    tier: str,
) -> None:
    tier_label = "Delta-only" if tier == "delta" else "Full-cut"
    canvas.Clear()
    canvas.Divide(3, 3)
    keepalive = []
    for pad_index, (label, _, _, _) in enumerate(KINEMATICS, start=1):
        pad = canvas.cd(pad_index)
        pad.SetGrid(1, 1)
        pad.SetBottomMargin(0.25)
        pad.SetLeftMargin(0.24)
        column = f"{label}_{tier}_residual_weighted_mean"
        values = [
            float(row[column]) for row in rows
            if isinstance(row.get(column), (float, int))
        ]
        if not values:
            message = unavailable_kinematic_message(rows, label, tier)
            pad.SetGrid(0, 0)
            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextAlign(22)
            latex.SetTextFont(62)
            latex.SetTextSize(0.060)
            latex.DrawLatex(
                0.5,
                0.80,
                f"{tier_label}: {label} recon-SIMC",
            )
            latex.SetTextFont(42)
            latex.SetTextSize(0.075)
            latex.DrawLatex(0.5, 0.48, message)
            keepalive.append(latex)
            continue
        low, high = padded_range(values)
        frame = categorical_frame(
            rows,
            f"{context}: {tier_label} kinematic residuals: {label} recon-SIMC",
            "Weighted mean residual",
            low, high,
        )
        frame.GetXaxis().SetLabelSize(0.012)
        frame.GetYaxis().SetTitleOffset(2.2)
        frame.GetYaxis().SetTitleSize(0.04)
        frame.Draw()
        keepalive.append(frame)
        settings = setting_names(rows)
        setting_index = {
            setting: index for index, setting in enumerate(settings, start=1)
        }
        keepalive.extend(draw_setting_separators(settings, low, high))
        for category in CATEGORY_ORDER:
            graph = ROOT.TGraph()
            point = 0
            for row in rows:
                if (row["Run_type"], row["Target"]) != category:
                    continue
                value = row.get(column, "")
                if not isinstance(value, (float, int)):
                    continue
                graph.SetPoint(point, grouped_x(row, setting_index), float(value))
                point += 1
            if point:
                graph.SetMarkerStyle(20)
                graph.SetMarkerColor(CATEGORY_COLORS[category])
                graph.SetMarkerSize(0.9)
                graph.Draw("P SAME")
                keepalive.append(graph)
    legend_pad = canvas.cd(9)
    legend_pad.SetGrid(0, 0)
    legend = ROOT.TLegend(0.08, 0.18, 0.92, 0.82)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.055)
    keepalive.append(legend)
    for category in CATEGORY_ORDER:
        dummy = ROOT.TGraph()
        dummy.SetMarkerStyle(20)
        dummy.SetMarkerColor(CATEGORY_COLORS[category])
        legend.AddEntry(dummy, LONG_CATEGORY_LABELS[category], "p")
        keepalive.append(dummy)
    legend.Draw()

    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(22)
    latex.SetTextFont(62)
    latex.SetTextSize(0.025)
    latex.DrawLatex(0.5, 0.995, f"{tier_label} kinematic residuals")
    latex.SetTextFont(42)
    latex.SetTextSize(0.018)
    latex.DrawLatex(0.5, 0.970, context)
    keepalive.append(latex)
    canvas.Print(pdf.as_posix())


def write_pdf(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    manifest: Path,
    t7_root: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetEndErrorSize(5)
    canvas = ROOT.TCanvas("c_simc_info", "SIMC QA", 1200, 1200)
    canvas.Print(path.as_posix() + "[")
    counts = Counter(str(row["QA_status"]) for row in rows)
    draw_text_page(canvas, path, "SIMC metadata, normalization, and QA", [
        f"Manifest: {manifest}",
        f"Simulation root: {t7_root}",
        "Variants: " + ", ".join(sorted({
            str(row.get("Simulation_variant") or "baseline") for row in rows
        })),
        "Collimator flags: " + ", ".join(sorted({
            (
                f"using_HMScoll={row.get('using_HMScoll') or 'default 0'}, "
                f"using_SHMScoll={row.get('using_SHMScoll') or 'default 0'}"
            )
            for row in rows
        })),
        f"Inventory rows: {len(rows)}",
        ", ".join(f"{status}: {count}" for status, count in sorted(counts.items())),
        "",
        "SIMC yield = sum(fWeight) after Boolean event selection.",
        "MC uncertainty = sqrt(sum(fWeight^{2})).",
        "N_{eff} = (sum fWeight)^{2} / sum(fWeight^{2}).",
        "fWeight = Weight #times normfac / Ngen; don't need to normalize it again.",
        "",
        "Acceptance tiers: nocut, delta-only, and full matched acceptance.",
        "Delta cut: -8 < hsdelta < 8, -10 < ssdelta < 22.",
        "Full-cut angle cuts: -0.15 < HMS/SHMS xptar < 0.15.",
        "Full-cut angle cuts: -0.10 < HMS/SHMS yptar < 0.10.",
        "Kinematic summaries are saved separately for Delta-only and Full-cut.",
        "Delta/exclusive xbj, z, pt2 use common derived original/recon coordinates.",
        "Derived momentum uses physical phad; provenance is retained in the CSV.",
        "Missing mass uses missmass and missmass_recon.",
        "Kinematic availability: available, sentinel_only, or missing.",
        "Sentinel-only branches (for example all -999) retain blank metrics.",
        "WARNING: negative fWeight, full relative error > 5%, or N_{eff} < 1000.",
        "PENDING and SKIPPED inventory rows are not classified as problematic.",
    ])
    status_lines = []
    by_phase_reaction = Counter(
        (str(row["Phase"]), str(row["Reaction"]), str(row["QA_status"]))
        for row in rows
    )
    for (phase, reaction, status), count in sorted(by_phase_reaction.items()):
        status_lines.append(f"{phase}, {reaction}, {status}: {count}")
    draw_text_page(canvas, path, "Inventory status by phase and reaction", status_lines)

    physics_groups = group_rows(rows)
    entry_groups = group_rows(rows, ("OK", "WARNING", "ERROR"))
    for phase, pass_name, reaction in sorted(
        set(physics_groups).union(entry_groups)
    ):
        context = f"{phase}, {pass_name}, {reaction}"
        entry_group = entry_groups.get((phase, pass_name, reaction), [])
        if entry_group:
            entry_chunks = chunk_rows_by_settings(entry_group)
            for chunk_index, entry_chunk in enumerate(entry_chunks, start=1):
                part = (
                    ""
                    if len(entry_chunks) == 1
                    else f" (part {chunk_index}/{len(entry_chunks)})"
                )
                draw_graph_page(
                    canvas,
                    path,
                    entry_chunk,
                    (
                        f"{context}: generated and reconstructed entry counts"
                        f"{part}"
                    ),
                    "Entries",
                    (
                        ("Raw entries", "Raw_entries", None),
                        ("Hist Ngen", "Hist_Ngen", None),
                        ("Recon entries", "Recon_entries", None),
                    ),
                    series_markers=(20, 21, 22),
                    nonnegative=True,
                )
        group = physics_groups.get((phase, pass_name, reaction), [])
        if not group:
            continue
        draw_graph_page(
            canvas, path, group, f"{context}: acceptance fractions", "Fraction",
            (
                ("Delta-only", "Acceptance_fraction_delta", None),
                ("Full-cut", "Acceptance_fraction_full", None),
            ),
            series_markers=(24, 20),
            nonnegative=True,
        )
        draw_graph_page(
            canvas,
            path,
            group,
            f"{context}: normalized yield: Delta-only and Full-cut",
            "sum(fWeight)",
            (
                ("Delta-only", "SimYield_delta", "SimYield_err_delta"),
                ("Full-cut", "SimYield_full", "SimYield_err_full"),
            ),
            series_markers=(24, 20),
            nonnegative=False,
            logy=all(
                float(row[column]) > 0
                for row in group
                for column in ("SimYield_delta", "SimYield_full")
            ),
        )
        draw_dual_panel_page(canvas, path, group, context)
        draw_kinematic_page(canvas, path, group, context, "delta")
        draw_kinematic_page(canvas, path, group, context, "full")
    canvas.Print(path.as_posix() + "]")


def filtered_manifest_rows(
    rows: Iterable[Dict[str, str]], phases: Sequence[str], reactions: Sequence[str]
) -> List[Dict[str, str]]:
    selected = []
    for row in rows:
        if phases and row["Phase"].lower() not in phases:
            continue
        if reactions and row["Reaction"].lower() not in reactions:
            continue
        selected.append(row)
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--t7-root", type=Path, default=DEFAULT_T7_ROOT)
    parser.add_argument("--phase", action="append", choices=PHASES, default=[])
    parser.add_argument("--reaction", action="append", choices=REACTIONS, default=[])
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--pdf-dir", type=Path, default=DEFAULT_PDF_DIR)
    args = parser.parse_args()

    manifest = args.manifest.expanduser().resolve()
    t7_root = args.t7_root.expanduser().resolve()
    table_dir = args.table_dir.expanduser().resolve()
    pdf_dir = args.pdf_dir.expanduser().resolve()
    if not manifest.is_file():
        raise FileNotFoundError(f"manifest not found: {manifest}")
    if not t7_root.is_dir():
        raise FileNotFoundError(f"simulation root unavailable: {t7_root}")

    source_rows = filtered_manifest_rows(
        read_manifest(manifest), args.phase, args.reaction
    )
    rows: List[Dict[str, object]] = []
    for index, source in enumerate(source_rows, start=1):
        row = inspect_manifest_row(source, t7_root)
        rows.append(row)
        print(
            f"[{index}/{len(source_rows)}] {row['QA_status']}: "
            f"{source['Reaction']} {source['Phase']} {source['Setting']}"
        )

    output_csv = table_dir / "RP_extract_simc_info.csv"
    problematic_csv = table_dir / "RP_extract_simc_info_Problematic.csv"
    output_pdf = pdf_dir / "RP_extract_simc_info.pdf"
    write_csv_atomic(output_csv, rows)
    problematic = [
        row for row in rows if row["QA_status"] in {"WARNING", "ERROR"}
    ]
    write_csv_atomic(problematic_csv, problematic)
    if not args.no_pdf:
        write_pdf(output_pdf, rows, manifest, t7_root)

    counts = Counter(str(row["QA_status"]) for row in rows)
    print("Status summary:")
    for status, count in sorted(counts.items()):
        print(f"  {status}: {count}")
    print(f"Master CSV: {output_csv}")
    print(f"Problematic CSV: {problematic_csv}")
    if not args.no_pdf:
        print(f"QA PDF: {output_pdf}")
    return 2 if counts["ERROR"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
