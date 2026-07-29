#!/usr/bin/env python3
"""Generate the isolated one-setting SIMC collimator material-model study."""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_SIMC_DIR = PROJECT_DIR.parents[1] / "simc_gfortran"
DEFAULT_BASELINE_MANIFEST = (
    DEFAULT_SIMC_DIR / "infiles" / "RP_Simc" / "simc_input_manifest.csv"
)
DEFAULT_OUTPUT_MANIFEST = (
    DEFAULT_SIMC_DIR
    / "infiles"
    / "RP_Simc"
    / "simc_input_manifest_collimator_on.csv"
)
OUTPUT_DIRECTORY = Path("infiles/RP_Simc/collimator_on")
VARIANT = "collimator_on"
REACTIONS = ("sidis", "rho", "delta", "exclusive")
STUDY_IDENTITY = {
    "Phase": "phase1",
    "Pass": "pass4",
    "Run_type": "PIMINUS",
    "Target": "LD2",
    "Setting": "x0p25Q23p3z0p5thpq2p0",
}
EXPECTED_NGEN = {
    "sidis": 100_000,
    "rho": 10_000,
    "delta": 10_000,
    "exclusive": 10_000,
}
RUNNABLE_STATUSES = {"GENERATED", "EXISTING_FILE_NOT_OVERWRITTEN"}
EXTRA_COLUMNS = (
    "Simulation_variant",
    "Parent_input",
    "using_HMScoll",
    "using_SHMScoll",
    "Random_initialization",
)


def read_manifest(path: Path) -> tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"manifest has no header: {path}")
        required = {
            *STUDY_IDENTITY,
            "Reaction",
            "Ngen",
            "Output_file",
            "Generation_status",
            "Diagnostic_reason",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"manifest is missing columns: {sorted(missing)}")
        return list(reader.fieldnames), list(reader)


def select_study_rows(rows: Iterable[Mapping[str, str]]) -> List[Dict[str, str]]:
    selected = [
        dict(row)
        for row in rows
        if all(row.get(name) == value for name, value in STUDY_IDENTITY.items())
    ]
    by_reaction = {row["Reaction"]: row for row in selected}
    if len(selected) != len(by_reaction):
        raise ValueError("baseline manifest contains duplicate study reactions")
    missing = set(REACTIONS).difference(by_reaction)
    extra = set(by_reaction).difference(REACTIONS)
    if missing or extra:
        raise ValueError(
            f"study inventory mismatch: missing={sorted(missing)}, extra={sorted(extra)}"
        )
    output = []
    for reaction in REACTIONS:
        row = by_reaction[reaction]
        if row["Generation_status"] not in RUNNABLE_STATUSES:
            raise ValueError(
                f"{reaction} baseline is not runnable: {row['Generation_status']}"
            )
        expected = EXPECTED_NGEN[reaction]
        if row.get("Ngen") != str(expected):
            raise ValueError(
                f"{reaction} Ngen mismatch: expected {expected}, found {row.get('Ngen')}"
            )
        output.append(row)
    return output


def assignment_values(text: str, name: str) -> List[str]:
    pattern = re.compile(
        rf"^\s*{re.escape(name)}\s*=\s*(?P<value>[^;\r\n]+)",
        re.MULTILINE,
    )
    return [match.group("value").strip() for match in pattern.finditer(text)]


def set_or_insert_simulate_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<prefix>\s*{re.escape(name)}\s*=\s*)(?P<value>[^;\r\n]+)",
        re.MULTILINE,
    )
    rendered, count = pattern.subn(
        lambda match: match.group("prefix") + value, text
    )
    if count == 1:
        return rendered
    if count > 1:
        raise ValueError(f"multiple {name} assignments in parent input")
    block = re.compile(
        r"(?P<body>^\s*begin\s+parm\s+simulate\b.*?)(?P<end>^\s*end\s+parm\s+simulate\b)",
        re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    matches = list(block.finditer(text))
    if len(matches) != 1:
        raise ValueError(
            f"expected one simulate parameter block, found {len(matches)}"
        )
    match = matches[0]
    indentation = "  "
    insertion = f"{indentation}{name} = {value};  collimator material-model study\n"
    return text[: match.start("end")] + insertion + text[match.start("end") :]


def render_study_input(parent_text: str) -> str:
    rendered = set_or_insert_simulate_assignment(
        parent_text, "using_HMScoll", "0"
    )
    rendered = set_or_insert_simulate_assignment(
        rendered, "using_SHMScoll", "1"
    )
    if assignment_values(rendered, "using_HMScoll") != ["0"]:
        raise ValueError("rendered using_HMScoll is not uniquely zero")
    if assignment_values(rendered, "using_SHMScoll") != ["1"]:
        raise ValueError("rendered using_SHMScoll is not uniquely one")
    random_states = assignment_values(rendered, "random_state_file")
    if any(value.strip("'\" ") for value in random_states):
        raise ValueError("parent input requests a non-empty evolved random state")
    seeds = assignment_values(rendered, "random_seed")
    if seeds and seeds != ["0"]:
        raise ValueError(f"parent input does not use deterministic seed zero: {seeds}")
    return rendered


def study_output_path(parent_output: str) -> Path:
    parent = Path(parent_output)
    if parent.is_absolute() or ".." in parent.parts or parent.suffix != ".inp":
        raise ValueError(f"unsafe parent input path: {parent_output}")
    return OUTPUT_DIRECTORY / f"{parent.stem}_collimatorOn.inp"


def changed_lines(parent_text: str, rendered_text: str) -> List[str]:
    parent_lines = parent_text.splitlines()
    rendered_lines = rendered_text.splitlines()
    additions = [line for line in rendered_lines if line not in parent_lines]
    removals = [line for line in parent_lines if line not in rendered_lines]
    return [f"+{line}" for line in additions] + [f"-{line}" for line in removals]


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def atomic_manifest(
    path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def generate(
    simc_dir: Path,
    baseline_manifest: Path,
    output_manifest: Path,
    write: bool,
    overwrite: bool,
) -> List[Dict[str, str]]:
    baseline_columns, baseline_rows = read_manifest(baseline_manifest)
    selected = select_study_rows(baseline_rows)
    output_rows = []
    for row in selected:
        parent_relative = Path(row["Output_file"])
        parent_path = simc_dir / parent_relative
        if not parent_path.is_file():
            raise FileNotFoundError(f"parent input is missing: {parent_path}")
        parent_text = parent_path.read_text(encoding="utf-8")
        rendered = render_study_input(parent_text)
        changes = changed_lines(parent_text, rendered)
        if len(changes) != 2 or not all(
            name in "\n".join(changes)
            for name in ("using_HMScoll", "using_SHMScoll")
        ):
            raise ValueError(
                f"{row['Reaction']} changed beyond two inserted assignments: {changes}"
            )
        output_relative = study_output_path(row["Output_file"])
        output_path = simc_dir / output_relative
        if output_path == parent_path:
            raise ValueError("experimental output aliases its baseline parent")
        if write:
            if output_path.exists() and not overwrite:
                raise FileExistsError(
                    f"experimental input exists; pass --overwrite: {output_path}"
                )
            atomic_write(output_path, rendered)
        output_row = dict(row)
        output_row.update({
            "Output_file": output_relative.as_posix(),
            "Generation_status": "GENERATED" if write else "WOULD_GENERATE",
            "Diagnostic_reason": (
                "isolated SHMS collimator material-model study; "
                "parent input otherwise unchanged"
            ),
            "Simulation_variant": VARIANT,
            "Parent_input": parent_relative.as_posix(),
            "using_HMScoll": "0",
            "using_SHMScoll": "1",
            "Random_initialization": "deterministic_default_random_seed_0",
        })
        output_rows.append(output_row)
    if write:
        columns = list(baseline_columns)
        for column in EXTRA_COLUMNS:
            if column not in columns:
                columns.append(column)
        atomic_manifest(output_manifest, columns, output_rows)
    return output_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--write", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--simc-dir", type=Path, default=DEFAULT_SIMC_DIR)
    parser.add_argument(
        "--baseline-manifest", type=Path, default=DEFAULT_BASELINE_MANIFEST
    )
    parser.add_argument(
        "--output-manifest", type=Path, default=DEFAULT_OUTPUT_MANIFEST
    )
    args = parser.parse_args()
    if args.overwrite and not args.write:
        parser.error("--overwrite requires --write")
    rows = generate(
        args.simc_dir.expanduser().resolve(),
        args.baseline_manifest.expanduser().resolve(),
        args.output_manifest.expanduser().resolve(),
        args.write,
        args.overwrite,
    )
    for row in rows:
        print(
            f"{row['Generation_status']}: {row['Reaction']} -> {row['Output_file']}"
        )
    print(f"Study rows: {len(rows)}")
    if args.write:
        print(f"Manifest: {args.output_manifest.expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
