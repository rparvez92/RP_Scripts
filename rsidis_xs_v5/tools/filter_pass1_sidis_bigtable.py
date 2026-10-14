#!/usr/bin/env python3
"""Split the Pass1 SIDIS bigtable into setting-ID CSV files.

Output layout::

    bigtable/Pass1/phaseN/passM/pim|pip/TARGET/Elec|Pos/SETTING_ID/
        phaseN_passM_pim|pip_TARGET_Elec|Pos_SETTING_ID.csv

Only SIDIS rows on LH2, LD2, and Dummy targets are selected.  Leaf files keep
the source columns in their original order; ``run_type`` and ``target`` are
normalized to match their directory labels.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PASS1_DIR = SCRIPT_DIR.parent / "bigtable" / "Pass1"
SOURCE = PASS1_DIR / "rp_rsidis_bigtable_pass1_sidis.csv"
SUMMARY = PASS1_DIR / "filter_pass1_sidis_bigtable_summary.csv"

# The energy values are the metadata distinction between the two run phases.
ENERGY_TO_PHASE_PASS = {
    8.5831: ("phase1", "pass4"),
    10.6716: ("phase1", "pass5"),
    6.4724: ("phase2", "pass3"),
    8.5814: ("phase2", "pass4"),
    10.6759: ("phase2", "pass5"),
}
PHASE_PASSES = {
    "phase1": ("pass4", "pass5"),
    "phase2": ("pass3", "pass4", "pass5"),
}
RUN_TYPE_MAP = {"PI-SIDIS": "pim", "PI+SIDIS": "pip"}
TARGET_MAP = {"LH2": "LH2", "LD2": "LD2", "Dummy": "DUMMY", "DUMMY": "DUMMY"}
RUN_TYPES = ("pim", "pip")
TARGETS = ("DUMMY", "LD2", "LH2")
CHARGES = ("Elec", "Pos")
KINEMATIC_COLUMNS = ("x", "Q2", "z", "thpq")
SUMMARY_COLUMNS = (
    "phase", "pass", "run_type", "target", "charge_category", "setting_id",
    "x", "Q2", "z", "thpq", "row_count", "output_file",
)

BaseCategory = Tuple[str, str, str, str, str]
GroupKey = Tuple[str, str, str, str, str, str]


def phase_pass(raw_energy: str) -> Tuple[str, str] | None:
    try:
        energy = float(raw_energy)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(energy):
        return None
    return ENERGY_TO_PHASE_PASS.get(energy)


def charge_category(raw_hms_p: str) -> str | None:
    try:
        momentum = float(raw_hms_p)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(momentum) or momentum in (0.0, -999.0):
        return None
    return "Elec" if momentum < 0.0 else "Pos"


def expected_categories() -> List[BaseCategory]:
    return [
        (phase, pass_name, run_type, target, charge)
        for phase, passes in PHASE_PASSES.items()
        for pass_name in passes
        for run_type in RUN_TYPES
        for target in TARGETS
        for charge in CHARGES
    ]


def category_dir(category: Sequence[str]) -> Path:
    return PASS1_DIR.joinpath(*category)


def output_path(key: GroupKey) -> Path:
    phase, pass_name, run_type, target, charge, setting_id = key
    filename = "_".join(key) + ".csv"
    return category_dir(key[:5]) / setting_id / filename


def run_sort_key(row: Mapping[str, str]) -> Tuple[int, str]:
    raw = row.get("run", "")
    try:
        return int(raw), raw
    except ValueError:
        return 2**63 - 1, raw


def write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if not SOURCE.is_file():
        raise FileNotFoundError(f"Missing source bigtable: {SOURCE}")

    groups: Dict[GroupKey, List[Dict[str, str]]] = defaultdict(list)
    exclusions: Counter[str] = Counter()
    for category in expected_categories():
        category_dir(category).mkdir(parents=True, exist_ok=True)

    with SOURCE.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"Source bigtable has no header: {SOURCE}")
        fieldnames = list(reader.fieldnames)
        required = {"run", "ebeam", "target", "hms_p", "run_type", "setting_id", *KINEMATIC_COLUMNS}
        missing = required.difference(fieldnames)
        if missing:
            raise ValueError(f"Source bigtable is missing required columns: {sorted(missing)}")

        for row in reader:
            pp = phase_pass(row["ebeam"])
            if pp is None:
                exclusions[f"ebeam:{row['ebeam'] or '<blank>'}"] += 1
                continue
            run_type = RUN_TYPE_MAP.get(row["run_type"].strip())
            if run_type is None:
                exclusions[f"run_type:{row['run_type'] or '<blank>'}"] += 1
                continue
            target = TARGET_MAP.get(row["target"].strip())
            if target is None:
                exclusions[f"target:{row['target'] or '<blank>'}"] += 1
                continue
            charge = charge_category(row["hms_p"])
            if charge is None:
                exclusions[f"hms_p:{row['hms_p'] or '<blank>'}"] += 1
                continue
            setting_id = row["setting_id"].strip()
            if not setting_id or setting_id in {".", ".."} or "/" in setting_id or "\\" in setting_id:
                exclusions["invalid_setting_id"] += 1
                continue
            if any(not row[column].strip() for column in KINEMATIC_COLUMNS):
                exclusions["missing_kinematics"] += 1
                continue

            normalized = dict(row)
            normalized["run_type"] = run_type
            normalized["target"] = target
            groups[(*pp, run_type, target, charge, setting_id)].append(normalized)

    summary_rows = []
    total_rows = 0
    for key in sorted(groups):
        rows = sorted(groups[key], key=run_sort_key)
        run_counts = Counter(row["run"] for row in rows)
        duplicate_runs = sorted(run for run, count in run_counts.items() if count > 1)
        if duplicate_runs:
            raise ValueError(f"Duplicate runs in {'/'.join(key)}: {', '.join(duplicate_runs)}")
        kinematics = {tuple(row[column].strip() for column in KINEMATIC_COLUMNS) for row in rows}
        if len(kinematics) != 1:
            raise ValueError(f"Setting ID has multiple kinematics in {'/'.join(key)}: {sorted(kinematics)}")

        path = output_path(key)
        write_csv(path, fieldnames, rows)
        x, q2, z, thpq = next(iter(kinematics))
        total_rows += len(rows)
        summary_rows.append(dict(zip(SUMMARY_COLUMNS, (*key, x, q2, z, thpq, len(rows), path.relative_to(PASS1_DIR).as_posix()))))

    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {len(groups)} setting-wise CSV files ({total_rows} selected rows)")
    print(f"Wrote summary: {SUMMARY}")
    if exclusions:
        print("Excluded source rows:")
        for reason, count in sorted(exclusions.items()):
            print(f"  {reason}: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
