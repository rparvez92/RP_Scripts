#!/usr/bin/env python3
"""Build deterministic groups of SIDIS settings that differ only in nominal thpq."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping

import RP_data_mc_compare as dmc

POLICY_VERSION = "RP_THPQ_COMBINATION_V1"
GROUP_FIELDS = ("Phase", "Pass", "Run_type", "Target", "x", "Q2", "z")
MAP_FIELDS = [
    "Combination_group_id", "Combination_policy_version", *GROUP_FIELDS,
    "Ebeam_GeV", "setting_id", "thpq_nominal_deg", "Contributor_index",
    "N_contributors", "Combination_status", "Combination_reason",
]
SUMMARY_FIELDS = [
    "Combination_group_id", "Combination_policy_version", *GROUP_FIELDS,
    "Ebeam_GeV", "N_contributors", "Contributing_setting_ids",
    "Contributing_thpq_deg", "Combination_status", "Combination_reason",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def token(value: object, digits: int) -> str:
    return f"{float(value):.{digits}f}".replace(".", "p")


def combination_group_id(row: Mapping[str, object]) -> str:
    phase = {"phase1": "1", "phase2": "2"}[str(row["Phase"])]
    pass_number = str(row["Pass"]).removeprefix("pass")
    charge = {"pim": "0", "pip": "1"}[str(row["Run_type"])]
    target = {"LH2": "1", "LD2": "2"}[str(row["Target"])]
    return (
        f"{phase}{pass_number}{charge}{target}"
        f"x{token(row['x'], 2)}Q2{token(row['Q2'], 1)}z{token(row['z'], 2)}"
    )


def build_rows(
    model_rows: Iterable[Mapping[str, str]],
    available: set[tuple[str, str, str, str, str]] | None = None,
) -> tuple[list[dict], list[dict]]:
    candidates = []
    for model in model_rows:
        if model.get("Model_status") != "OK":
            continue
        identity = tuple(model.get(field, "") for field in
                         ("Phase", "Pass", "Run_type", "Target", "setting_id"))
        if available is not None and identity not in available:
            continue
        row = {
            "Phase": model.get("Phase", ""), "Pass": model.get("Pass", ""),
            "Run_type": model.get("Run_type", ""), "Target": model.get("Target", ""),
            "x": model.get("x", ""), "Q2": model.get("Q2_GeV2", ""),
            "z": model.get("z", ""), "Ebeam_GeV": model.get("Ebeam_GeV", ""),
            "setting_id": model.get("setting_id", ""),
            "thpq_nominal_deg": model.get("theta_pq_deg", ""),
        }
        candidates.append(row)
    groups: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    for row in candidates:
        key = tuple(str(row[field]) for field in GROUP_FIELDS)
        groups[key].append(row)
    mapping, summaries = [], []
    seen_ids: dict[str, tuple[str, ...]] = {}
    for key, members in sorted(groups.items()):
        members.sort(key=lambda row: (float(row["thpq_nominal_deg"]), row["setting_id"]))
        group_id = combination_group_id(members[0])
        reasons = []
        energies = {token(row["Ebeam_GeV"], 4) for row in members if finite(row["Ebeam_GeV"])}
        if len(energies) != 1 or len(energies) != len({token(row["Ebeam_GeV"], 4) for row in members}):
            reasons.append("INCONSISTENT_OR_NONFINITE_EBEAM")
        settings = [str(row["setting_id"]) for row in members]
        if len(settings) != len(set(settings)):
            reasons.append("DUPLICATE_SETTING_ID")
        previous = seen_ids.get(group_id)
        if previous is not None and previous != key:
            reasons.append("COMBINATION_GROUP_ID_COLLISION")
        seen_ids[group_id] = key
        status = "ERROR" if reasons else "OK"
        reason = ";".join(reasons)
        for index, member in enumerate(members, 1):
            mapped = {field: "" for field in MAP_FIELDS}
            mapped.update(member, Combination_group_id=group_id,
                          Combination_policy_version=POLICY_VERSION,
                          Contributor_index=index, N_contributors=len(members),
                          Combination_status=status, Combination_reason=reason)
            mapping.append(mapped)
        summary = {field: "" for field in SUMMARY_FIELDS}
        summary.update({field: members[0][field] for field in (*GROUP_FIELDS, "Ebeam_GeV")})
        summary.update(
            Combination_group_id=group_id, Combination_policy_version=POLICY_VERSION,
            N_contributors=len(members), Contributing_setting_ids="|".join(settings),
            Contributing_thpq_deg="|".join(str(row["thpq_nominal_deg"]) for row in members),
            Combination_status=status, Combination_reason=reason,
        )
        summaries.append(summary)
    return mapping, summaries


def atomic_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def main() -> int:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-catalog", type=Path, default=(
        project / "results/Pass1/Tables/RP_sidis_model/RP_sidis_model.csv"))
    parser.add_argument("--norm-dir", type=Path, default=dmc.DEFAULT_NORM_DIR)
    parser.add_argument("--output-dir", type=Path, default=(
        project / "results/Pass1/Tables/RP_thpq_combination"))
    args = parser.parse_args()
    norm_groups = dmc.discover_norm_files(args.norm_dir)
    available = {
        (*base[:3], target, base[3]) for base, files in norm_groups.items()
        for target in ("LH2", "LD2") if (target, "Elec") in files
    }
    mapping, summaries = build_rows(dmc.read_csv(args.model_catalog), available)
    atomic_csv(args.output_dir / "RP_thpq_combination_map.csv", MAP_FIELDS, mapping)
    atomic_csv(args.output_dir / "RP_thpq_combination_groups.csv", SUMMARY_FIELDS, summaries)
    problematic = [row for row in mapping if row["Combination_status"] != "OK"]
    atomic_csv(args.output_dir / "RP_thpq_combination_Problematic.csv", MAP_FIELDS, problematic)
    print(f"Groups: {len(summaries)}; contributors: {len(mapping)}; problematic: {len(problematic)}")
    return 1 if problematic else 0


if __name__ == "__main__":
    raise SystemExit(main())
