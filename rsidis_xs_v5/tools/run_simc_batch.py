#!/usr/bin/env python3
"""Run generated SIMC inputs serially from simc_input_manifest.csv."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_SIMC_DIR = PROJECT_DIR.parents[1] / "simc_gfortran"
RUNNABLE_GENERATION_STATUSES = {"GENERATED", "EXISTING_FILE_NOT_OVERWRITTEN"}
REACTIONS = ("sidis", "rho", "delta", "exclusive")

STATUS_COLUMNS = (
    "Timestamp_UTC",
    "Phase",
    "Pass",
    "Run_type",
    "Target",
    "Setting",
    "Reaction",
    "Input_file",
    "Requested_Ngen",
    "Observed_Ngen",
    "Batch_status",
    "Diagnostic_reason",
)


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {
            "Phase", "Pass", "Run_type", "Target", "Setting", "Reaction",
            "Ngen", "Output_file", "Generation_status",
        }
        if reader.fieldnames is None:
            raise ValueError("manifest has no header")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"manifest is missing columns: {sorted(missing)}")
        return list(reader)


def select_rows(
    rows: Iterable[Mapping[str, str]],
    phases: Sequence[str],
    reactions: Sequence[str],
    run_types: Sequence[str],
    targets: Sequence[str],
) -> List[Dict[str, str]]:
    selected: List[Dict[str, str]] = []
    seen_inputs: set[str] = set()
    for source_row in rows:
        row = dict(source_row)
        if not row.get("Output_file"):
            continue
        if row.get("Generation_status") not in RUNNABLE_GENERATION_STATUSES:
            continue
        if phases and row["Phase"] not in phases:
            continue
        if reactions and row["Reaction"] not in reactions:
            continue
        if run_types and row["Run_type"] not in run_types:
            continue
        if targets and row["Target"] not in targets:
            continue
        input_file = row["Output_file"]
        if input_file in seen_inputs:
            raise ValueError(f"duplicate runnable input in manifest: {input_file}")
        seen_inputs.add(input_file)
        selected.append(row)
    return selected


def input_relative_to_infiles(row: Mapping[str, str]) -> Path:
    path = Path(row["Output_file"])
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe manifest output path: {path}")
    if not path.parts or path.parts[0] != "infiles":
        raise ValueError(f"manifest output is not below infiles/: {path}")
    relative = Path(*path.parts[1:])
    if len(relative.parts) < 2:
        raise ValueError(f"input lacks a run-type parent directory: {path}")
    return relative


def expected_paths(
    simc_dir: Path,
    t7_root: Path,
    row: Mapping[str, str],
) -> Dict[str, Path]:
    relative_input = input_relative_to_infiles(row)
    stem = relative_input.stem
    run_type = relative_input.parent.name
    phase = row["Phase"].lower()
    if phase not in {"phase1", "phase2"}:
        raise ValueError(f"unsupported phase: {row['Phase']}")
    phase_dir = "Phase" + phase[-1]
    simulation = t7_root / phase_dir / "Simulation"
    return {
        "input": simc_dir / "infiles" / relative_input,
        "hist": simulation / "outfiles" / run_type / f"{stem}.hist",
        "runout": simulation / "runout" / run_type / f"{stem}.out",
        "raw_root": simulation / "ROOTfiles" / run_type / f"{stem}.root",
        "recon_root": simulation / "ROOTfiles" / run_type / f"recon_hcana_{stem}.root",
    }


def read_hist_ngen(path: Path) -> int | None:
    if not path.is_file():
        return None
    pattern = re.compile(
        r"^\s*Ngen(?:\s*\(request\))?\s*=\s*([0-9]+)",
        re.IGNORECASE,
    )
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            for line in stream:
                match = pattern.match(line)
                if match:
                    return int(match.group(1))
    except OSError:
        return None
    return None


def root_has_h10(path: Path, rootls: str = "rootls") -> bool:
    if not path.is_file():
        return False
    result = subprocess.run(
        [rootls, "-1", str(path)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0 and "h10" in result.stdout.splitlines()


def inspect_outputs(
    paths: Mapping[str, Path],
    requested_ngen: int,
    rootls: str = "rootls",
) -> tuple[str, str, int | None]:
    existing = [name for name in ("hist", "runout", "raw_root", "recon_root")
                if paths[name].exists()]
    if not existing:
        return "PENDING", "no output files found", None

    observed_ngen = read_hist_ngen(paths["hist"])
    complete_files = len(existing) == 4
    readable_roots = (
        root_has_h10(paths["raw_root"], rootls)
        and root_has_h10(paths["recon_root"], rootls)
    )
    if complete_files and observed_ngen == requested_ngen and readable_roots:
        return "COMPLETE", "complete output set with readable h10 trees", observed_ngen

    reasons: List[str] = []
    missing = [name for name in ("hist", "runout", "raw_root", "recon_root")
               if name not in existing]
    if missing:
        reasons.append("missing=" + ",".join(missing))
    if observed_ngen != requested_ngen:
        reasons.append(f"Ngen={observed_ngen}, requested={requested_ngen}")
    if complete_files and not readable_roots:
        reasons.append("one or both ROOT files lack readable h10 trees")
    return "PARTIAL_OR_MISMATCH", "; ".join(reasons), observed_ngen


def status_row(
    row: Mapping[str, str],
    input_file: Path,
    requested_ngen: int,
    observed_ngen: int | None,
    status: str,
    reason: str,
) -> Dict[str, str]:
    return {
        "Timestamp_UTC": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "Phase": row["Phase"],
        "Pass": row["Pass"],
        "Run_type": row["Run_type"],
        "Target": row["Target"],
        "Setting": row["Setting"],
        "Reaction": row["Reaction"],
        "Input_file": input_file.as_posix(),
        "Requested_Ngen": str(requested_ngen),
        "Observed_Ngen": "" if observed_ngen is None else str(observed_ngen),
        "Batch_status": status,
        "Diagnostic_reason": reason,
    }


def write_status(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=STATUS_COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="inspect without running SIMC")
    mode.add_argument("--run", action="store_true", help="run selected jobs serially")
    parser.add_argument("--simc-dir", type=Path, default=DEFAULT_SIMC_DIR)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--t7-root", type=Path,
                        default=Path(os.environ.get("SIMC_T7_ROOT", "/Volumes/T7/RSIDIS")))
    parser.add_argument("--phase", action="append", choices=("phase1", "phase2"), default=[])
    parser.add_argument("--reaction", action="append", choices=REACTIONS, default=[])
    parser.add_argument("--run-type", action="append", choices=("PIPLUS", "PIMINUS"), default=[])
    parser.add_argument("--target", action="append", choices=("LH2", "LD2"), default=[])
    parser.add_argument("--ngen", type=int, help="temporary event-count override")
    parser.add_argument("--overwrite", action="store_true",
                        help="permit replacement of partial/mismatched output sets")
    parser.add_argument("--limit", type=int, help="process at most this many selected jobs")
    parser.add_argument("--status-file", type=Path)
    args = parser.parse_args()

    if args.ngen is not None and args.ngen <= 0:
        parser.error("--ngen must be a positive integer")
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")

    simc_dir = args.simc_dir.expanduser().resolve()
    manifest = (args.manifest or
                simc_dir / "infiles" / "RP_Simc" / "simc_input_manifest.csv")
    manifest = manifest.expanduser().resolve()
    t7_root = args.t7_root.expanduser().resolve()
    status_file = (args.status_file or
                   simc_dir / "infiles" / "RP_Simc" / "simc_batch_status.csv")
    status_file = status_file.expanduser().resolve()
    runner = simc_dir / "run_simc_recon.sh"

    if not runner.is_file():
        raise FileNotFoundError(f"SIMC runner not found: {runner}")
    if not t7_root.is_dir():
        raise FileNotFoundError(f"T7 root is unavailable: {t7_root}")

    rows = select_rows(
        read_manifest(manifest), args.phase, args.reaction, args.run_type, args.target
    )
    if args.limit is not None:
        rows = rows[:args.limit]
    if not rows:
        print("No runnable manifest rows matched the requested filters.")
        return 1

    results: List[Dict[str, str]] = []
    counts: Counter = Counter()
    for index, row in enumerate(rows, start=1):
        paths = expected_paths(simc_dir, t7_root, row)
        if not paths["input"].is_file():
            status, reason, observed = "ERROR_MISSING_INPUT", str(paths["input"]), None
        else:
            requested_ngen = args.ngen if args.ngen is not None else int(row["Ngen"])
            status, reason, observed = inspect_outputs(paths, requested_ngen)
            if args.run and status != "COMPLETE":
                if status == "PARTIAL_OR_MISMATCH" and not args.overwrite:
                    status = "BLOCKED_NEEDS_OVERWRITE"
                    reason += "; pass --overwrite to replace this output set"
                else:
                    command = [
                        str(runner), input_relative_to_infiles(row).as_posix(),
                        row["Reaction"], "mpi", "1",
                    ]
                    if args.ngen is not None:
                        command.extend(["--ngen", str(args.ngen)])
                    if args.overwrite:
                        command.append("--overwrite")
                    print(f"[{index}/{len(rows)}] RUN {' '.join(command[1:])}", flush=True)
                    completed = subprocess.run(command, cwd=simc_dir, check=False)
                    if completed.returncode == 0:
                        status, reason, observed = inspect_outputs(paths, requested_ngen)
                        if status != "COMPLETE":
                            status = "ERROR_OUTPUT_VALIDATION"
                    else:
                        status = "ERROR_RUN_FAILED"
                        reason = f"runner exited with status {completed.returncode}"
        requested_ngen = args.ngen if args.ngen is not None else int(row["Ngen"])
        result = status_row(
            row, input_relative_to_infiles(row), requested_ngen, observed, status, reason
        )
        results.append(result)
        counts[status] += 1
        print(f"[{index}/{len(rows)}] {status}: {row['Output_file']} ({reason})")
        if args.run:
            write_status(status_file, results)

    print(f"Selected jobs: {len(rows)}")
    for status, count in sorted(counts.items()):
        print(f"  {status}: {count}")
    if args.run:
        print(f"Status file: {status_file}")
    return 1 if any(status.startswith("ERROR") or status.startswith("BLOCKED")
                    for status in counts) else 0


if __name__ == "__main__":
    raise SystemExit(main())
