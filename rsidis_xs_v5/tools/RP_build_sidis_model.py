#!/usr/bin/env python3
"""Build the central SIDIS model catalog used by the yield-ratio extraction."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Iterable

MP_GEV = 0.9382720813
MPI_GEV = 0.13957039
M_SIDIS_UNITS = "GeV^-2"
SIGMA_SIDIS_UNITS = "ub/GeV^2/sr^2"

FIELDS = [
    "Phase", "Pass", "Run_type", "Target", "Setting", "Ebeam_GeV",
    "x", "Q2_GeV2", "z", "theta_pq_deg", "nu_GeV", "p_pion_GeV",
    "pT_GeV", "Target_A", "Target_Z", "Pion_charge", "Phi_model_rad",
    "M_sidis_model", "M_sidis_units", "sigma_sidis_model",
    "sigma_sidis_units", "Calculator", "Calculator_git_commit",
    "Source_manifest", "Model_status", "Model_reason",
]


def decode_token(token: str) -> float:
    """Decode setting tokens such as 3p3, neg0p8, and neg999."""
    sign = -1.0 if token.startswith("neg") else 1.0
    if token.startswith("neg"):
        token = token[3:]
    return sign * float(token.replace("p", "."))


SETTING_RE = re.compile(
    r"^x(?P<x>(?:neg)?[0-9]+(?:p[0-9]+)?)"
    r"Q2(?P<Q2>(?:neg)?[0-9]+(?:p[0-9]+)?)"
    r"z(?P<z>(?:neg)?[0-9]+(?:p[0-9]+)?)"
    r"thpq(?P<thpq>(?:neg)?[0-9]+(?:p[0-9]+)?)$"
)


def parse_setting(setting: str) -> dict[str, float]:
    match = SETTING_RE.fullmatch(setting)
    if not match:
        raise ValueError(f"invalid setting name: {setting}")
    return {key: decode_token(value) for key, value in match.groupdict().items()}


def derive_kinematics(ebeam: float, x: float, q2: float, z: float,
                      theta_deg: float) -> dict[str, float]:
    values = (ebeam, x, q2, z, theta_deg)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("nonfinite central kinematics")
    if min(ebeam, x, q2, z) <= 0 or any(abs(value + 999.0) < 1e-9 for value in values):
        raise ValueError("sentinel or nonphysical central kinematics")
    nu = q2 / (2.0 * MP_GEV * x)
    eprime = ebeam - nu
    epi = z * nu
    if eprime <= 0:
        raise ValueError("central scattered-electron energy is not positive")
    if epi <= MPI_GEV:
        raise ValueError("central pion energy is below its mass")
    ppion = math.sqrt(epi * epi - MPI_GEV * MPI_GEV)
    # calc_semi_xsec expects the transverse-momentum magnitude.  The sign of
    # theta_pq carries orientation information but must not make pT negative.
    pt = abs(ppion * math.sin(math.radians(theta_deg)))
    if pt < 0 or pt > ppion * (1.0 + 1e-12):
        raise ValueError("derived pT is outside the physical range")
    return {"nu": nu, "ppion": ppion, "pt": pt}


def calculator_identity(target: str, run_type: str) -> tuple[float, float, float]:
    if target == "LH2":
        target_a, target_z = 1.0, 1.0
    elif target == "LD2":
        target_a, target_z = 2.0, 1.0
    else:
        raise ValueError(f"unsupported target {target}")
    if run_type == "PIPLUS":
        charge = 1.0
    elif run_type == "PIMINUS":
        charge = -1.0
    else:
        raise ValueError(f"unsupported run type {run_type}")
    return target_a, target_z, charge


def parse_calculator_output(output: str) -> tuple[float, float, float]:
    """Return phi, sighad, and sigsemi from the final numeric output row."""
    candidates: list[tuple[float, float, float]] = []
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            values = tuple(float(item.replace("D", "E")) for item in parts)
        except ValueError:
            continue
        if all(math.isfinite(value) for value in values):
            candidates.append(values)  # type: ignore[arg-type]
    if not candidates:
        raise ValueError("calculator output contains no phi/sighad/sigsemi row")
    return candidates[-1]


def run_calculator(executable: Path, inputs: Iterable[float]) -> tuple[float, float, float]:
    payload = " ".join(f"{value:.15g}" for value in inputs) + "\n"
    result = subprocess.run(
        [str(executable)], input=payload, text=True, capture_output=True, check=False,
        cwd=executable.parent,
    )
    if result.returncode:
        detail = (result.stderr or result.stdout).strip().replace("\n", " ")
        raise RuntimeError(f"calc_semi_xsec exited {result.returncode}: {detail[:400]}")
    return parse_calculator_output(result.stdout)


def atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def calculator_commit(simc_root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=simc_root, text=True,
        capture_output=True, check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else ""


def build_catalog(manifest: Path, executable: Path) -> list[dict[str, object]]:
    with manifest.open(newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    identities: dict[tuple[str, ...], dict[str, str]] = {}
    for source in source_rows:
        if source.get("Reaction", "").lower() != "sidis":
            continue
        key = tuple(source.get(name, "") for name in
                    ("Phase", "Pass", "Run_type", "Target", "Setting"))
        identities.setdefault(key, source)

    simc_root = executable.resolve().parents[2]
    git_commit = calculator_commit(simc_root)
    rows: list[dict[str, object]] = []
    for key in sorted(identities):
        source = identities[key]
        phase, pass_name, run_type, target, setting = key
        row: dict[str, object] = {name: "" for name in FIELDS}
        row.update({
            "Phase": phase, "Pass": pass_name, "Run_type": run_type,
            "Target": target, "Setting": setting,
            "Calculator": str(executable.resolve()),
            "Calculator_git_commit": git_commit,
            "Source_manifest": str(manifest.resolve()),
            "M_sidis_units": M_SIDIS_UNITS,
            "sigma_sidis_units": SIGMA_SIDIS_UNITS,
        })
        try:
            parsed = parse_setting(setting)
            generation = source.get("Generation_status", "")
            if generation not in {"GENERATED", "OK"}:
                model_status = "SKIPPED" if generation.startswith("SKIPPED") else "PENDING"
                row.update({"Model_status": model_status,
                            "Model_reason": f"sidis_input_status={generation or 'MISSING'}"})
                rows.append(row)
                continue
            ebeam = float(source["Ebeam_MeV"]) / 1000.0
            row.update({"Ebeam_GeV": ebeam, "x": parsed["x"],
                        "Q2_GeV2": parsed["Q2"], "z": parsed["z"],
                        "theta_pq_deg": parsed["thpq"]})
            derived = derive_kinematics(
                ebeam, parsed["x"], parsed["Q2"], parsed["z"], parsed["thpq"]
            )
            row.update({"nu_GeV": derived["nu"], "p_pion_GeV": derived["ppion"],
                        "pT_GeV": derived["pt"]})
            target_a, target_z, charge = calculator_identity(target, run_type)
            row.update({"Target_A": target_a, "Target_Z": target_z,
                        "Pion_charge": charge})
            if not executable.is_file() or not os.access(executable, os.X_OK):
                row.update({"Model_status": "PENDING", "Model_reason": "CALCULATOR_MISSING"})
            else:
                phi, sighad, sigma = run_calculator(
                    executable,
                    (ebeam, parsed["Q2"], parsed["x"], parsed["z"], derived["pt"],
                     target_a, target_z, charge, 1.0),
                )
                if sighad <= 0 or sigma <= 0:
                    raise ValueError("calculator returned nonpositive model values")
                row.update({"Phi_model_rad": phi, "M_sidis_model": sighad,
                            "sigma_sidis_model": sigma, "Model_status": "OK"})
        except ValueError as exc:
            text = str(exc)
            status = "SKIPPED" if "sentinel" in text or "setting name" in text else "ERROR"
            row.update({"Model_status": status, "Model_reason": text})
        except (OSError, RuntimeError) as exc:
            row.update({"Model_status": "ERROR", "Model_reason": str(exc)})
        rows.append(row)
    return rows


def main() -> int:
    here = Path(__file__).resolve()
    analysis_root = here.parents[1]
    jlab_root = here.parents[3]
    default_simc = jlab_root / "simc_gfortran"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=default_simc / "infiles/RP_Simc/simc_input_manifest.csv")
    parser.add_argument("--calculator", type=Path,
                        default=default_simc / "util/sidisxsec/calc_semi_xsec")
    parser.add_argument("--output", type=Path,
                        default=analysis_root / "results/Tables/RP_sidis_model/RP_sidis_model.csv")
    args = parser.parse_args()
    rows = build_catalog(args.manifest, args.calculator)
    atomic_write(args.output, rows)
    counts: dict[str, int] = {}
    for row in rows:
        status = str(row["Model_status"])
        counts[status] = counts.get(status, 0) + 1
    print(f"Wrote {len(rows)} model rows to {args.output}")
    print("Status:", ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 1 if counts.get("ERROR", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
