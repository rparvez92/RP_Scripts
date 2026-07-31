#!/usr/bin/env python3
"""Generate setting-wise Hall-C SIMC coincidence input files.

The setting-wise bigtable leaves are the source of beam and spectrometer
settings.  Trusted target/reaction templates are read from
``simc_gfortran/infiles/pdbforce`` and are never modified.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_LEAF_ROOT = PROJECT_DIR / "bigtable"
DEFAULT_SIMC_DIR = PROJECT_DIR.parents[1] / "simc_gfortran"

PHASES = ("phase1", "phase2")
PASSES = ("pass3", "pass4", "pass5")
RUN_TYPES = ("PIPLUS", "PIMINUS")
TARGETS = ("LH2", "LD2")
REACTIONS = ("sidis", "rho", "delta", "exclusive")

BEAM_MEV = {
    "pass3": Decimal("6449.0"),
    "pass4": Decimal("8583.1"),
    "pass5": Decimal("10671.6"),
}
BEAM_GEV = {
    "pass3": Decimal("6.4490"),
    "pass4": Decimal("8.5831"),
    "pass5": Decimal("10.6716"),
}
EVENT_COUNTS = {
    "sidis": 100_000,
    "rho": 10_000,
    "delta": 10_000,
    "exclusive": 10_000,
}
TEMPLATE_REACTION_NAMES = {
    "sidis": "sidis",
    "rho": "rho",
    "delta": "delta",
    "exclusive": "excl",
}
TARGET_NUMBERS = {
    "LH2": (Decimal("1"), Decimal("1")),
    "LD2": (Decimal("2"), Decimal("1")),
}

STATUS_GENERATED = "GENERATED"
STATUS_WOULD_GENERATE = "WOULD_GENERATE"
STATUS_INVALID = "SKIPPED_INVALID_KINEMATICS"
STATUS_CONFLICT = "SKIPPED_SETTING_CONFLICT"
STATUS_EXISTS = "EXISTING_FILE_NOT_OVERWRITTEN"
STATUS_TEMPLATE_ERROR = "ERROR_TEMPLATE_VALIDATION"

MANIFEST_COLUMNS = (
    "Phase",
    "Pass",
    "Run_type",
    "Target",
    "Setting",
    "Reaction",
    "Ngen",
    "Source_leaf_csv",
    "Source_template",
    "Output_file",
    "Ebeam_MeV",
    "HMS_P_MeV",
    "HMS_theta_deg",
    "SHMS_P_MeV",
    "SHMS_theta_deg",
    "which_pion",
    "doing_hplus",
    "Generation_status",
    "Diagnostic_reason",
    "Run_count",
    "Runs",
)


@dataclass(frozen=True)
class LeafIdentity:
    phase: str
    pass_name: str
    run_type: str
    target: str
    setting: str


@dataclass
class LeafResult:
    identity: LeafIdentity
    source: Path
    status: str
    reason: str
    ebeam_mev: str = ""
    hms_p_mev: str = ""
    hms_theta_deg: str = ""
    shms_p_mev: str = ""
    shms_theta_deg: str = ""
    runs: Tuple[int, ...] = ()


def parse_decimal(raw: str) -> Decimal | None:
    try:
        value = Decimal(raw.strip())
    except (AttributeError, InvalidOperation):
        return None
    return value if value.is_finite() else None


def format_decimal(value: Decimal, quantum: str) -> str:
    return format(value.quantize(Decimal(quantum), rounding=ROUND_HALF_UP), "f")


def parse_leaf_identity(path: Path, leaf_root: Path) -> LeafIdentity | None:
    try:
        relative = path.relative_to(leaf_root)
    except ValueError:
        return None
    parts = relative.parts
    # phase/pass/run-type/target/Elec/setting/file.csv
    if len(parts) != 7 or parts[4] != "Elec":
        return None
    phase, pass_name, run_type, target, _, setting, filename = parts
    expected = f"{phase}_{pass_name}_{run_type}_{target}_Elec_{setting}.csv"
    if filename != expected:
        return None
    if (
        phase not in PHASES
        or pass_name not in PASSES
        or run_type not in RUN_TYPES
        or target not in TARGETS
    ):
        return None
    return LeafIdentity(phase, pass_name, run_type, target, setting)


def discover_leaves(leaf_root: Path) -> List[Tuple[LeafIdentity, Path]]:
    discovered: List[Tuple[LeafIdentity, Path]] = []
    for path in sorted(leaf_root.glob("phase*/pass*/PI*/L?2/Elec/*/*.csv")):
        identity = parse_leaf_identity(path, leaf_root)
        if identity is not None:
            discovered.append((identity, path))
    return discovered


def validate_leaf(identity: LeafIdentity, path: Path, leaf_root: Path) -> LeafResult:
    required = {
        "run",
        "ebeam",
        "target",
        "run_type",
        "hms_p",
        "hms_th",
        "shms_p",
        "shms_th",
        "x",
        "Q2",
        "z",
        "thpq",
    }
    try:
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ValueError("CSV has no header")
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            rows = list(reader)
    except (OSError, csv.Error, ValueError) as error:
        return LeafResult(identity, path, STATUS_INVALID, str(error))

    if not rows:
        return LeafResult(identity, path, STATUS_INVALID, "leaf contains no rows")

    expected_ebeam = BEAM_GEV[identity.pass_name]
    tuples: Dict[Tuple[str, str, str, str], List[int]] = {}
    runs: List[int] = []
    errors: List[str] = []
    for row_number, row in enumerate(rows, start=2):
        try:
            run = int(row["run"])
        except (TypeError, ValueError):
            errors.append(f"row {row_number}: invalid run")
            continue
        runs.append(run)
        if row["target"].strip() != identity.target:
            errors.append(f"run {run}: target mismatch")
        if row["run_type"].strip() != identity.run_type:
            errors.append(f"run {run}: run_type mismatch")

        ebeam = parse_decimal(row["ebeam"])
        settings = [parse_decimal(row[name]) for name in ("x", "Q2", "z", "thpq")]
        if ebeam is None or ebeam != expected_ebeam:
            errors.append(f"run {run}: ebeam does not match {identity.pass_name}")
        if any(value is None or value == Decimal("-999") for value in settings):
            errors.append(f"run {run}: invalid x/Q2/z/thpq")

        hms_p = parse_decimal(row["hms_p"])
        hms_th = parse_decimal(row["hms_th"])
        shms_p = parse_decimal(row["shms_p"])
        shms_th = parse_decimal(row["shms_th"])
        if any(value is None or value == Decimal("-999") for value in
               (hms_p, hms_th, shms_p, shms_th)):
            errors.append(f"run {run}: invalid spectrometer setting")
            continue
        assert hms_p is not None and hms_th is not None
        assert shms_p is not None and shms_th is not None
        formatted = (
            format_decimal(abs(hms_p) * Decimal("1000"), "0.1"),
            format_decimal(hms_th, "0.001"),
            format_decimal(abs(shms_p) * Decimal("1000"), "0.1"),
            format_decimal(shms_th, "0.001"),
        )
        tuples.setdefault(formatted, []).append(run)

    source_text = path.relative_to(leaf_root).as_posix()
    if errors:
        reason = "; ".join(errors[:12])
        if len(errors) > 12:
            reason += f"; ... {len(errors) - 12} additional error(s)"
        return LeafResult(identity, path, STATUS_INVALID, reason,
                          runs=tuple(sorted(set(runs))))
    if len(tuples) != 1:
        descriptions = [
            f"{values}:runs={'|'.join(map(str, sorted(tuple_runs)))}"
            for values, tuple_runs in sorted(tuples.items())
        ]
        return LeafResult(
            identity,
            path,
            STATUS_CONFLICT,
            "distinct formatted spectrometer tuples: " + "; ".join(descriptions),
            runs=tuple(sorted(set(runs))),
        )

    values = next(iter(tuples))
    return LeafResult(
        identity=identity,
        source=path,
        status="VALID",
        reason=f"validated {source_text}",
        ebeam_mev=format_decimal(BEAM_MEV[identity.pass_name], "0.1"),
        hms_p_mev=values[0],
        hms_theta_deg=values[1],
        shms_p_mev=values[2],
        shms_theta_deg=values[3],
        runs=tuple(sorted(set(runs))),
    )


def reaction_flags(reaction: str, run_type: str) -> Dict[str, int]:
    positive = run_type == "PIPLUS"
    if reaction == "sidis":
        return {
            "doing_pion": 1,
            "which_pion": 0,
            "doing_rho": 0,
            "doing_semi": 1,
            "doing_hplus": int(positive),
        }
    if reaction == "rho":
        return {
            "doing_pion": 0,
            "which_pion": 0,
            "doing_rho": 1,
            "doing_semi": 0,
            "doing_hplus": int(positive),
        }
    if reaction == "exclusive":
        return {
            "doing_pion": 1,
            "which_pion": 0 if positive else 1,
            "doing_rho": 0,
            "doing_semi": 0,
            "doing_hplus": int(positive),
        }
    if reaction == "delta":
        return {
            "doing_pion": 1,
            "which_pion": 2 if positive else 3,
            "doing_rho": 0,
            "doing_semi": 0,
            "doing_hplus": int(positive),
        }
    raise ValueError(f"Unknown reaction: {reaction}")


def excluded_combination(identity: LeafIdentity, reaction: str) -> bool:
    return (
        reaction == "exclusive"
        and identity.run_type == "PIMINUS"
        and identity.target == "LH2"
    )


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<prefix>\s*{re.escape(name)}\s*=\s*)(?P<value>[^;\r\n]+)",
        re.MULTILINE,
    )
    result, count = pattern.subn(lambda match: match.group("prefix") + value, text)
    if count != 1:
        raise ValueError(f"expected exactly one assignment for {name}, found {count}")
    return result


def assignment_value(text: str, name: str) -> str:
    pattern = re.compile(
        rf"^\s*{re.escape(name)}\s*=\s*(?P<value>[^;\r\n]+)",
        re.MULTILINE,
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one assignment for {name}, found {len(matches)}")
    return matches[0].group("value").strip()


def validate_rendered_input(text: str, leaf: LeafResult, reaction: str) -> None:
    """Verify identity, kinematics, target, and reaction flags before writing."""
    flags = reaction_flags(reaction, leaf.identity.run_type)
    expected_text = {
        "ngen": str(EVENT_COUNTS[reaction]),
        "doing_pion": str(flags["doing_pion"]),
        "which_pion": str(flags["which_pion"]),
        "doing_rho": str(flags["doing_rho"]),
        "doing_semi": str(flags["doing_semi"]),
        "doing_hplus": str(flags["doing_hplus"]),
        "electron_arm": "1",
        "hadron_arm": "5",
    }
    for name, expected in expected_text.items():
        actual = assignment_value(text, name)
        if actual != expected:
            raise ValueError(f"{name}: expected {expected}, found {actual}")

    expected_numbers = {
        "Ebeam": Decimal(leaf.ebeam_mev),
        "spec%e%P": Decimal(leaf.hms_p_mev),
        "spec%e%theta": Decimal(leaf.hms_theta_deg),
        "spec%p%P": Decimal(leaf.shms_p_mev),
        "spec%p%theta": Decimal(leaf.shms_theta_deg),
        "Egamma_gen_max": Decimal(leaf.ebeam_mev),
        "targ%A": TARGET_NUMBERS[leaf.identity.target][0],
        "targ%Z": TARGET_NUMBERS[leaf.identity.target][1],
    }
    for name, expected in expected_numbers.items():
        actual = parse_decimal(assignment_value(text, name))
        if actual != expected:
            raise ValueError(f"{name}: expected {expected}, found {actual}")


def render_input(template_text: str, leaf: LeafResult, reaction: str) -> Tuple[str, Dict[str, int]]:
    flags = reaction_flags(reaction, leaf.identity.run_type)
    replacements = {
        "ngen": str(EVENT_COUNTS[reaction]),
        "doing_pion": str(flags["doing_pion"]),
        "which_pion": str(flags["which_pion"]),
        "doing_rho": str(flags["doing_rho"]),
        "doing_semi": str(flags["doing_semi"]),
        "doing_hplus": str(flags["doing_hplus"]),
        "Ebeam": leaf.ebeam_mev,
        "spec%e%P": leaf.hms_p_mev,
        "spec%e%theta": leaf.hms_theta_deg,
        "spec%p%P": leaf.shms_p_mev,
        "spec%p%theta": leaf.shms_theta_deg,
        "Egamma_gen_max": leaf.ebeam_mev,
    }
    rendered = template_text
    for name, value in replacements.items():
        rendered = replace_assignment(rendered, name, value)
    return rendered, flags


def template_path(
    simc_dir: Path, reaction: str, target: str, run_type: str
) -> Path:
    template_reaction = TEMPLATE_REACTION_NAMES[reaction]
    charge = "pip" if run_type == "PIPLUS" else "pim"
    return (
        simc_dir
        / "infiles"
        / "pdbforce"
        / (
            f"mc_{template_reaction}_{target}_{charge}_"
            "e10p7_x0p25_q23p3_z0p36_thpq2p0.inp"
        )
    )


def output_name(identity: LeafIdentity, reaction: str) -> str:
    return (
        f"mc_{reaction}_{identity.phase}_{identity.pass_name}_"
        f"{identity.run_type}_{identity.target}_{identity.setting}.inp"
    )


def manifest_row(
    leaf: LeafResult,
    reaction: str,
    source_template: str,
    output_file: str,
    flags: Mapping[str, int],
    status: str,
    reason: str,
    leaf_root: Path,
) -> Dict[str, str]:
    return {
        "Phase": leaf.identity.phase,
        "Pass": leaf.identity.pass_name,
        "Run_type": leaf.identity.run_type,
        "Target": leaf.identity.target,
        "Setting": leaf.identity.setting,
        "Reaction": reaction,
        "Ngen": str(EVENT_COUNTS[reaction]),
        "Source_leaf_csv": leaf.source.relative_to(leaf_root).as_posix(),
        "Source_template": source_template,
        "Output_file": output_file,
        "Ebeam_MeV": leaf.ebeam_mev,
        "HMS_P_MeV": leaf.hms_p_mev,
        "HMS_theta_deg": leaf.hms_theta_deg,
        "SHMS_P_MeV": leaf.shms_p_mev,
        "SHMS_theta_deg": leaf.shms_theta_deg,
        "which_pion": str(flags.get("which_pion", "")),
        "doing_hplus": str(flags.get("doing_hplus", "")),
        "Generation_status": status,
        "Diagnostic_reason": reason,
        "Run_count": str(len(leaf.runs)),
        "Runs": ";".join(map(str, leaf.runs)),
    }


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


def write_manifest(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=MANIFEST_COLUMNS, lineterminator="\n")
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
    leaf_root: Path,
    simc_dir: Path,
    write: bool,
    overwrite: bool,
) -> Tuple[List[Dict[str, str]], Counter]:
    output_dir = simc_dir / "infiles" / "RP_Simc" / "coin"
    manifest_rows: List[Dict[str, str]] = []
    statuses: Counter = Counter()
    seen_outputs: set[str] = set()
    leaves = discover_leaves(leaf_root)
    if not leaves:
        raise FileNotFoundError(f"No eligible setting leaves found under {leaf_root}")

    for identity, source in leaves:
        leaf = validate_leaf(identity, source, leaf_root)
        for reaction in REACTIONS:
            if excluded_combination(identity, reaction):
                continue
            flags = reaction_flags(reaction, identity.run_type)
            template = template_path(
                simc_dir, reaction, identity.target, identity.run_type
            )
            relative_template = template.relative_to(simc_dir).as_posix()
            output = output_dir / output_name(identity, reaction)
            relative_output = output.relative_to(simc_dir).as_posix()

            if leaf.status != "VALID":
                status = leaf.status
                reason = leaf.reason
            elif not template.is_file():
                status = STATUS_TEMPLATE_ERROR
                reason = f"missing template: {template}"
            elif relative_output in seen_outputs:
                status = STATUS_TEMPLATE_ERROR
                reason = f"duplicate output path: {relative_output}"
            else:
                seen_outputs.add(relative_output)
                try:
                    template_text = template.read_text(encoding="utf-8")
                    rendered, flags = render_input(template_text, leaf, reaction)
                    validate_rendered_input(rendered, leaf, reaction)
                except (OSError, UnicodeError, ValueError) as error:
                    status = STATUS_TEMPLATE_ERROR
                    reason = str(error)
                else:
                    if write:
                        if output.exists() and not overwrite:
                            status = STATUS_EXISTS
                            reason = "use --overwrite to replace this generated input"
                        else:
                            atomic_write(output, rendered)
                            status = STATUS_GENERATED
                            reason = "generated from validated leaf and trusted template"
                    else:
                        status = STATUS_WOULD_GENERATE
                        reason = "dry run: input validated; no file written"

            statuses[status] += 1
            manifest_rows.append(
                manifest_row(
                    leaf,
                    reaction,
                    relative_template,
                    relative_output if status not in {STATUS_INVALID, STATUS_CONFLICT} else "",
                    flags,
                    status,
                    reason,
                    leaf_root,
                )
            )

    manifest_rows.sort(
        key=lambda row: (
            row["Phase"],
            row["Pass"],
            row["Run_type"],
            row["Target"],
            row["Setting"],
            row["Reaction"],
        )
    )
    if write:
        write_manifest(simc_dir / "infiles" / "RP_Simc" / "simc_input_manifest.csv",
                       manifest_rows)
    return manifest_rows, statuses


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate SIMC coincidence inputs from setting-wise bigtable leaves."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="validate without writing")
    mode.add_argument("--write", action="store_true", help="write inputs and manifest")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace existing generated .inp files (requires --write)",
    )
    parser.add_argument(
        "--leaf-root",
        type=Path,
        default=DEFAULT_LEAF_ROOT,
        help=f"bigtable directory containing phase leaves (default: {DEFAULT_LEAF_ROOT})",
    )
    parser.add_argument(
        "--simc-dir",
        type=Path,
        default=DEFAULT_SIMC_DIR,
        help=f"SIMC checkout directory (default: {DEFAULT_SIMC_DIR})",
    )
    args = parser.parse_args()
    if args.overwrite and not args.write:
        parser.error("--overwrite requires --write")

    leaf_root = args.leaf_root.expanduser().resolve()
    simc_dir = args.simc_dir.expanduser().resolve()
    if not leaf_root.is_dir():
        raise FileNotFoundError(f"Leaf root does not exist: {leaf_root}")
    if not (simc_dir / "infiles" / "pdbforce").is_dir():
        raise FileNotFoundError(f"SIMC pdbforce template directory not found under {simc_dir}")

    rows, statuses = generate(leaf_root, simc_dir, args.write, args.overwrite)
    print(f"Eligible setting leaves: {len(rows) // len(REACTIONS)}")
    print(f"Manifest rows: {len(rows)}")
    for status, count in sorted(statuses.items()):
        print(f"  {status}: {count}")
    if args.write:
        print(f"Inputs: {simc_dir / 'infiles' / 'RP_Simc' / 'coin'}")
        print(f"Manifest: {simc_dir / 'infiles' / 'RP_Simc' / 'simc_input_manifest.csv'}")
    else:
        print("Dry run only: no inputs or manifest were written.")
    return 1 if statuses[STATUS_TEMPLATE_ERROR] else 0


if __name__ == "__main__":
    raise SystemExit(main())
