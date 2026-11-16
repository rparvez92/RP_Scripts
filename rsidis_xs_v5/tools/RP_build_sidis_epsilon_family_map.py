#!/usr/bin/env python3
"""Build cross-pass epsilon families from the combined-theta SIDIS groups.

This stage maps measurements only.  It does not calculate reduced cross
sections, select phi points, fit harmonics, or perform an L/T separation.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

from RP_sidis_reduced_kinematics import calculate_reduced_kinematics


POLICY_VERSION = "RP_SIDIS_EPSILON_FAMILY_V1"
EXPECTED_ROWS_PER_GROUP = 75
PASS_ORDER = {"pass3": 3, "pass4": 4, "pass5": 5}

FAMILY_KEY_FIELDS = ("Run_type", "Target", "x", "Q2", "z")
SOURCE_GROUP_FIELDS = (
    "Combination_group_id", "Combination_policy_version", "Phase", "Pass",
    "Run_type", "Target", "x", "Q2", "z", "Ebeam_GeV", "N_contributors",
    "Contributing_setting_ids", "Contributing_thpq_deg", "Binning_version",
)
MEMBER_FIELDS = [
    "Epsilon_family_id", "Epsilon_family_policy_version", "Family_status",
    "Family_reason", *FAMILY_KEY_FIELDS, "N_members", "Pass_pattern",
    "Has_pass3", "Has_pass4", "Has_pass5", "N_distinct_epsilon",
    "Epsilon_min", "Epsilon_max", "Epsilon_span", "Member_index",
    "Member_phase", "Member_pass", "Combination_group_id",
    "Combination_policy_version", "Ebeam_GeV", "nu_GeV", "Eprime_GeV",
    "theta_e_rad", "y", "gamma", "epsilon", "N_contributors",
    "Contributing_setting_ids", "Contributing_thpq_deg", "Binning_version",
]
SUMMARY_FIELDS = [
    "Epsilon_family_id", "Epsilon_family_policy_version", "Family_status",
    "Family_reason", *FAMILY_KEY_FIELDS, "N_members", "Pass_pattern",
    "Has_pass3", "Has_pass4", "Has_pass5", "N_distinct_epsilon",
    "Epsilon_min", "Epsilon_max", "Epsilon_span", "Member_group_ids",
    "Member_phases", "Member_passes", "Member_beam_energies_GeV",
    "Member_epsilons",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def compact_number(value: object) -> str:
    text = f"{float(value):.8g}"
    return text.replace("-", "m").replace(".", "p")


def family_id(row: Mapping[str, object]) -> str:
    return (
        f"EF_{row['Run_type']}_{row['Target']}_x{compact_number(row['x'])}"
        f"Q2{compact_number(row['Q2'])}_z{compact_number(row['z'])}"
    )


def atomic_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name, suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def collapse_combination_groups(
    source_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Validate each 75-bin group and return one provenance row per group."""
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in source_rows:
        grouped[str(row["Combination_group_id"])].append(row)
    output = []
    for group_id, rows in sorted(grouped.items()):
        if len(rows) != EXPECTED_ROWS_PER_GROUP:
            raise ValueError(
                f"{group_id} contains {len(rows)} rows; expected "
                f"{EXPECTED_ROWS_PER_GROUP}"
            )
        first = rows[0]
        for field in SOURCE_GROUP_FIELDS:
            if len({str(row.get(field, "")) for row in rows}) != 1:
                raise ValueError(f"{group_id} has inconsistent {field}")
        if str(first["Pass"]) not in PASS_ORDER:
            raise ValueError(f"{group_id} has unsupported pass {first['Pass']}")
        for field in ("x", "Q2", "z", "Ebeam_GeV"):
            if not finite(first[field]):
                raise ValueError(f"{group_id} has invalid {field}")
        output.append({field: first.get(field, "") for field in SOURCE_GROUP_FIELDS})
    return output


def distinct_count(values: Sequence[float], tolerance: float = 1.0e-10) -> int:
    distinct: list[float] = []
    for value in sorted(values):
        if not distinct or not math.isclose(value, distinct[-1], rel_tol=0.0,
                                             abs_tol=tolerance):
            distinct.append(value)
    return len(distinct)


def family_decision(members: Sequence[Mapping[str, object]]) -> dict[str, object]:
    passes = [str(member["Pass"]) for member in members]
    counts = Counter(passes)
    epsilon = [float(member["epsilon"]) for member in members]
    reasons = []
    if any(count > 1 for count in counts.values()):
        reasons.append("DUPLICATE_PASS_MEMBER")
    if "pass4" not in counts:
        reasons.append("MISSING_PASS4")
    if "pass5" not in counts:
        reasons.append("MISSING_PASS5_ANCHOR")
    n_distinct = distinct_count(epsilon)
    if n_distinct < 2:
        reasons.append("INSUFFICIENT_DISTINCT_EPSILON")
    ready = not any(reason in reasons for reason in (
        "DUPLICATE_PASS_MEMBER", "MISSING_PASS4", "MISSING_PASS5_ANCHOR",
        "INSUFFICIENT_DISTINCT_EPSILON",
    ))
    status = (
        "LT_READY_PASS3_VALIDATION" if ready and "pass3" in counts
        else "LT_READY" if ready else "NOT_LT_READY"
    )
    ordered_passes = sorted(counts, key=lambda value: PASS_ORDER[value])
    return {
        "Family_status": status, "Family_reason": ";".join(reasons),
        "N_members": len(members), "Pass_pattern": "|".join(ordered_passes),
        "Has_pass3": int("pass3" in counts), "Has_pass4": int("pass4" in counts),
        "Has_pass5": int("pass5" in counts), "N_distinct_epsilon": n_distinct,
        "Epsilon_min": min(epsilon), "Epsilon_max": max(epsilon),
        "Epsilon_span": max(epsilon) - min(epsilon),
    }


def build_family_map(
    source_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    groups = collapse_combination_groups(source_rows)
    families: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for group in groups:
        kin = calculate_reduced_kinematics(
            beam_energy=float(group["Ebeam_GeV"]), x=float(group["x"]),
            Q2=float(group["Q2"]), z=float(group["z"]), pt=0.0,
        )
        group.update({
            "nu_GeV": kin.nu, "Eprime_GeV": kin.eprime,
            "theta_e_rad": kin.theta_e, "y": kin.y, "gamma": kin.gamma,
            "epsilon": kin.epsilon,
        })
        key = tuple(str(group[field]) for field in FAMILY_KEY_FIELDS)
        families[key].append(group)
    member_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for _, members in sorted(families.items()):
        members.sort(key=lambda row: (PASS_ORDER[str(row["Pass"])],
                                      float(row["Ebeam_GeV"]),
                                      str(row["Combination_group_id"])))
        first = members[0]
        decision = family_decision(members)
        common = {
            "Epsilon_family_id": family_id(first),
            "Epsilon_family_policy_version": POLICY_VERSION,
            **{field: first[field] for field in FAMILY_KEY_FIELDS},
            **decision,
        }
        for index, member in enumerate(members, 1):
            member_rows.append({
                **common, "Member_index": index,
                "Member_phase": member["Phase"], "Member_pass": member["Pass"],
                **{field: member[field] for field in (
                    "Combination_group_id", "Combination_policy_version",
                    "Ebeam_GeV", "nu_GeV", "Eprime_GeV", "theta_e_rad",
                    "y", "gamma", "epsilon", "N_contributors",
                    "Contributing_setting_ids", "Contributing_thpq_deg",
                    "Binning_version",
                )},
            })
        summary_rows.append({
            **common,
            "Member_group_ids": "|".join(str(row["Combination_group_id"]) for row in members),
            "Member_phases": "|".join(str(row["Phase"]) for row in members),
            "Member_passes": "|".join(str(row["Pass"]) for row in members),
            "Member_beam_energies_GeV": "|".join(str(row["Ebeam_GeV"]) for row in members),
            "Member_epsilons": "|".join(f"{float(row['epsilon']):.16g}" for row in members),
        })
    return member_rows, summary_rows


def main() -> int:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--combined-yields", type=Path,
        default=project / "results/Pass1/Tables/RP_combine_thpq_yields/RP_combined_thpq_yields.csv",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=project / "results/Pass1/Tables/RP_sidis_epsilon_family_map",
    )
    args = parser.parse_args()
    members, summaries = build_family_map(read_csv(args.combined_yields))
    atomic_csv(args.output_dir / "RP_sidis_epsilon_family_members.csv",
               MEMBER_FIELDS, members)
    atomic_csv(args.output_dir / "RP_sidis_epsilon_family_summary.csv",
               SUMMARY_FIELDS, summaries)
    atomic_csv(
        args.output_dir / "RP_sidis_epsilon_family_summary_Problematic.csv",
        SUMMARY_FIELDS,
        [row for row in summaries if not str(row["Family_status"]).startswith("LT_READY")],
    )
    print("Members:", len(members), "Families:", len(summaries),
          "Status:", Counter(row["Family_status"] for row in summaries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
