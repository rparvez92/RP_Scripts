#!/usr/bin/env python3
"""Split phase bigtables into setting-wise analysis CSV files.

Generated hierarchy:

  bigtable/phaseN/passM/RUN_TYPE/TARGET/CHARGE/SETTING/
    phaseN_passM_RUN_TYPE_TARGET_CHARGE_SETTING.csv

Leaf CSVs retain exactly the source bigtable columns and column order. Only
the values of ``run_type`` and ``target`` are normalized for portable naming.
This script creates and overwrites current outputs but never deletes files.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, List, Mapping, MutableMapping, Sequence, Set, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
BIGTABLE_DIR = PROJECT_DIR / "bigtable"

PHASE_CONFIG = {
    "phase1": {
        "source": BIGTABLE_DIR / "rsidis_bigtable_phase1.csv",
        "passes": ("pass4", "pass5"),
    },
    "phase2": {
        "source": BIGTABLE_DIR / "rsidis_bigtable_phase2.csv",
        "passes": ("pass3", "pass4", "pass5"),
    },
}

ENERGY_TO_PASS = {
    6.4490: "pass3",
    8.5831: "pass4",
    10.6716: "pass5",
}

RUN_TYPE_MAP = {
    "PI+SIDIS": "PIPLUS",
    "PI-SIDIS": "PIMINUS",
    "HMSDIS": "HMSDIS",
    "SHMSDIS": "SHMSDIS",
}

RUN_TYPES = ("PIPLUS", "PIMINUS", "HMSDIS", "SHMSDIS")
TARGETS = ("LH2", "LD2", "DUMMY")
CHARGES = ("Elec", "Pos")

TARGET_MAP = {
    "LH2": "LH2",
    "LD2": "LD2",
    "DUMMY": "DUMMY",
    "Dummy": "DUMMY",
}

SETTING_COLUMNS = ("x", "Q2", "z", "thpq")
SUMMARY_COLUMNS = (
    "phase",
    "pass",
    "run_type",
    "target",
    "charge_category",
    "x",
    "Q2",
    "z",
    "thpq",
    "setting",
    "row_count",
    "output_file",
)

BaseCategory = Tuple[str, str, str, str, str]
SettingValues = Tuple[Decimal, Decimal, Decimal, Decimal]
GroupKey = Tuple[BaseCategory, SettingValues]


def parse_pass(raw_energy: str) -> str | None:
    """Return the pass for an exact supported beam-energy value."""
    try:
        energy = float(raw_energy)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(energy):
        return None
    return ENERGY_TO_PASS.get(energy)


def parse_charge(raw_hms_p: str) -> str | None:
    """Classify HMS polarity, rejecting missing/sentinel/zero values."""
    try:
        hms_p = float(raw_hms_p)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(hms_p) or hms_p == 0.0 or hms_p == -999.0:
        return None
    return "Elec" if hms_p < 0.0 else "Pos"


def parse_decimal(raw_value: str) -> Decimal | None:
    try:
        value = Decimal(raw_value.strip())
    except (AttributeError, InvalidOperation):
        return None
    return value if value.is_finite() else None


def parse_setting(row: Mapping[str, str]) -> SettingValues | None:
    values = tuple(parse_decimal(row.get(column, "")) for column in SETTING_COLUMNS)
    if any(value is None for value in values):
        return None
    return values  # type: ignore[return-value]


def decimal_text(value: Decimal) -> str:
    """Canonical display text used in the summary CSV."""
    if value == value.to_integral():
        if value == Decimal("-999"):
            return "-999"
        return f"{value.quantize(Decimal('1'))}.0"
    return format(value.normalize(), "f")


def path_number(value: Decimal) -> str:
    """Convert a setting number into a portable path component."""
    text = decimal_text(value)
    if text.startswith("-"):
        text = "neg" + text[1:]
    return text.replace(".", "p")


def setting_name(values: SettingValues) -> str:
    x_value, q2_value, z_value, thpq_value = values
    return (
        f"x{path_number(x_value)}"
        f"Q2{path_number(q2_value)}"
        f"z{path_number(z_value)}"
        f"thpq{path_number(thpq_value)}"
    )


def run_sort_key(row: Mapping[str, str]) -> Tuple[int, str]:
    raw_run = row.get("run", "")
    try:
        return int(raw_run), raw_run
    except ValueError:
        return (2**63 - 1), raw_run


def expected_base_categories() -> List[BaseCategory]:
    categories: List[BaseCategory] = []
    for phase, config in PHASE_CONFIG.items():
        for pass_name in config["passes"]:
            for run_type in RUN_TYPES:
                for target in TARGETS:
                    for charge in CHARGES:
                        categories.append(
                            (phase, pass_name, run_type, target, charge)
                        )
    return categories


def category_directory(category: BaseCategory) -> Path:
    phase, pass_name, run_type, target, charge = category
    return BIGTABLE_DIR / phase / pass_name / run_type / target / charge


def output_path(category: BaseCategory, values: SettingValues) -> Path:
    phase, pass_name, run_type, target, charge = category
    setting = setting_name(values)
    filename = (
        f"{phase}_{pass_name}_{run_type}_{target}_{charge}_{setting}.csv"
    )
    return category_directory(category) / setting / filename


def read_and_filter_phase(
    phase: str,
    source: Path,
    allowed_passes: Sequence[str],
    grouped_rows: MutableMapping[GroupKey, List[Dict[str, str]]],
    exclusions: Counter,
) -> List[str]:
    if not source.is_file():
        raise FileNotFoundError(f"Missing source bigtable: {source}")

    with source.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"Source bigtable has no header: {source}")

        fieldnames = list(reader.fieldnames)
        required = {
            "run",
            "ebeam",
            "target",
            "hms_p",
            "run_type",
            *SETTING_COLUMNS,
        }
        missing = required.difference(fieldnames)
        if missing:
            raise ValueError(
                f"{source} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            source_run_type = row["run_type"].strip()
            run_type = RUN_TYPE_MAP.get(source_run_type)
            if run_type is None:
                exclusions[(phase, f"run_type:{source_run_type or '<blank>'}")] += 1
                continue

            source_target = row["target"].strip()
            target = TARGET_MAP.get(source_target)
            if target is None:
                exclusions[(phase, f"target:{source_target or '<blank>'}")] += 1
                continue

            pass_name = parse_pass(row["ebeam"])
            if pass_name is None:
                exclusions[(phase, f"ebeam:{row['ebeam'] or '<blank>'}")] += 1
                continue
            if pass_name not in allowed_passes:
                exclusions[(phase, f"pass_not_in_phase:{pass_name}")] += 1
                continue

            charge = parse_charge(row["hms_p"])
            if charge is None:
                exclusions[(phase, f"hms_p:{row['hms_p'] or '<blank>'}")] += 1
                continue

            values = parse_setting(row)
            if values is None:
                exclusions[(phase, "invalid_setting")] += 1
                continue

            normalized = dict(row)
            normalized["run_type"] = run_type
            normalized["target"] = target
            category = (phase, pass_name, run_type, target, charge)
            grouped_rows[(category, values)].append(normalized)

    return fieldnames


def write_leaf_csv(
    path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def check_duplicate_runs(
    category: BaseCategory,
    values: SettingValues,
    rows: Sequence[Mapping[str, str]],
) -> None:
    counts = Counter(row.get("run", "") for row in rows)
    duplicates = sorted(run for run, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(
            f"Duplicate runs in {'/'.join(category)}/{setting_name(values)}: "
            f"{', '.join(duplicates)}"
        )


def existing_leaf_csvs() -> Set[Path]:
    files: Set[Path] = set()
    for phase in PHASE_CONFIG:
        phase_dir = BIGTABLE_DIR / phase
        if phase_dir.is_dir():
            files.update(path.resolve() for path in phase_dir.rglob("*.csv"))
    return files


def main() -> int:
    categories = expected_base_categories()
    grouped_rows: Dict[GroupKey, List[Dict[str, str]]] = defaultdict(list)
    phase_headers: Dict[str, List[str]] = {}
    exclusions: Counter = Counter()

    for category in categories:
        category_directory(category).mkdir(parents=True, exist_ok=True)

    existing_before = existing_leaf_csvs()

    for phase, config in PHASE_CONFIG.items():
        phase_headers[phase] = read_and_filter_phase(
            phase=phase,
            source=config["source"],
            allowed_passes=config["passes"],
            grouped_rows=grouped_rows,
            exclusions=exclusions,
        )

    generated_paths: Set[Path] = set()
    settings_by_category: Dict[BaseCategory, List[SettingValues]] = defaultdict(list)
    for category, values in grouped_rows:
        settings_by_category[category].append(values)

    summary_rows = []
    total_selected = 0
    for category in categories:
        phase, pass_name, run_type, target, charge = category
        values_list = sorted(set(settings_by_category.get(category, [])))

        if not values_list:
            summary_rows.append(
                {
                    "phase": phase,
                    "pass": pass_name,
                    "run_type": run_type,
                    "target": target,
                    "charge_category": charge,
                    "x": "",
                    "Q2": "",
                    "z": "",
                    "thpq": "",
                    "setting": "",
                    "row_count": 0,
                    "output_file": "",
                }
            )
            continue

        for values in values_list:
            rows = sorted(grouped_rows[(category, values)], key=run_sort_key)
            check_duplicate_runs(category, values, rows)

            path = output_path(category, values)
            write_leaf_csv(path, phase_headers[phase], rows)
            generated_paths.add(path.resolve())

            setting = setting_name(values)
            setting_text = [decimal_text(value) for value in values]
            total_selected += len(rows)
            summary_rows.append(
                {
                    "phase": phase,
                    "pass": pass_name,
                    "run_type": run_type,
                    "target": target,
                    "charge_category": charge,
                    "x": setting_text[0],
                    "Q2": setting_text[1],
                    "z": setting_text[2],
                    "thpq": setting_text[3],
                    "setting": setting,
                    "row_count": len(rows),
                    "output_file": path.relative_to(BIGTABLE_DIR).as_posix(),
                }
            )

    summary_path = BIGTABLE_DIR / "filter_bigtable_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=SUMMARY_COLUMNS, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    stale_paths = sorted(existing_before.difference(generated_paths))

    print(f"Wrote {len(generated_paths)} setting-wise CSV files")
    print(f"  selected rows: {total_selected}")
    print(f"  empty charge categories: {sum(not settings_by_category.get(c) for c in categories)}")
    print(f"Wrote summary: {summary_path}")
    if stale_paths:
        print("WARNING: stale CSV files were left untouched:")
        for path in stale_paths:
            print(f"  {path.relative_to(BIGTABLE_DIR)}")
    if exclusions:
        print("Excluded source rows:")
        for (phase, reason), count in sorted(exclusions.items()):
            print(f"  {phase:6s} {reason}: {count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
