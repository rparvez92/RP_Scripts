#!/usr/bin/env python3
"""Build the Pass1 setting-ID SIDIS model catalog for yield-ratio extraction."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

MP_GEV = 0.9382720813
MPI_GEV = 0.13957039
IDENTITY = ("Phase", "Pass", "Run_type", "Target", "setting_id")
M_SIDIS_UNITS = "GeV^-2"
SIGMA_SIDIS_UNITS = "ub/GeV^2/sr^2"
FIELDS = [
    *IDENTITY, "Ebeam_GeV", "x", "Q2_GeV2", "z", "theta_pq_deg",
    "nu_GeV", "p_pion_GeV", "pT_GeV", "Target_A", "Target_Z",
    "Pion_charge", "Phi_model_rad", "M_sidis_model", "M_sidis_units",
    "sigma_sidis_model", "sigma_sidis_units", "Calculator",
    "Calculator_git_commit", "Source_manifest", "Model_status", "Model_reason",
]


def derive_kinematics(ebeam: float, x: float, q2: float, z: float,
                      theta_deg: float) -> dict[str, float]:
    values = (ebeam, x, q2, z, theta_deg)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("nonfinite central kinematics")
    if min(ebeam, x, q2, z) <= 0 or any(abs(value + 999.0) < 1e-9 for value in values):
        raise ValueError("sentinel or nonphysical central kinematics")
    nu = q2 / (2.0 * MP_GEV * x)
    pion_energy = z * nu
    if ebeam - nu <= 0 or pion_energy <= MPI_GEV:
        raise ValueError("nonphysical final-state energy")
    pion_momentum = math.sqrt(pion_energy * pion_energy - MPI_GEV * MPI_GEV)
    pt = abs(pion_momentum * math.sin(math.radians(theta_deg)))
    return {"nu": nu, "ppion": pion_momentum, "pt": pt}


def calculator_identity(target: str, run_type: str) -> tuple[float, float, float]:
    target_values = {"LH2": (1.0, 1.0), "LD2": (2.0, 1.0)}
    charges = {"pip": 1.0, "pim": -1.0}
    if target not in target_values or run_type not in charges:
        raise ValueError(f"unsupported identity {target}/{run_type}")
    return (*target_values[target], charges[run_type])


def parse_calculator_rows(output: str) -> list[tuple[float, float, float]]:
    candidates = []
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            values = tuple(float(item.replace("D", "E")) for item in parts)
        except ValueError:
            continue
        if all(math.isfinite(value) for value in values):
            candidates.append(values)
    if not candidates:
        raise ValueError("calculator output contains no phi/sighad/sigsemi row")
    return candidates


def parse_calculator_output(output: str) -> tuple[float, float, float]:
    return parse_calculator_rows(output)[-1]


def run_calculator(executable: Path, values: tuple[float, ...]) -> tuple[float, float, float]:
    return run_calculator_rows(executable, values)[-1]


def run_calculator_rows(
    executable: Path, values: tuple[float, ...]
) -> list[tuple[float, float, float]]:
    result = subprocess.run(
        [str(executable)], input=" ".join(f"{value:.15g}" for value in values) + "\n",
        text=True, capture_output=True, check=False, cwd=executable.parent,
    )
    if result.returncode:
        raise RuntimeError(f"calc_semi_xsec exited {result.returncode}")
    return parse_calculator_rows(result.stdout)


def calculator_commit(executable: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=executable.resolve().parents[2],
        text=True, capture_output=True, check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else ""


def build_catalog(manifest: Path, executable: Path) -> list[dict[str, object]]:
    with manifest.open(newline="", encoding="utf-8-sig") as stream:
        source_rows = list(csv.DictReader(stream))
    sources = {}
    for source in source_rows:
        if source.get("reaction", "").lower() != "sidis":
            continue
        key = (source["phase"], source["pass"], source["charge"],
               source["target"], source["setting_id"])
        sources.setdefault(key, source)
    commit = calculator_commit(executable) if executable.is_file() else ""
    rows = []
    for key in sorted(sources):
        source = sources[key]
        row = {field: "" for field in FIELDS}
        row.update(dict(zip(IDENTITY, key)))
        row.update({
            "Calculator": str(executable.resolve()),
            "Calculator_git_commit": commit,
            "Source_manifest": str(manifest.resolve()),
            "M_sidis_units": M_SIDIS_UNITS,
            "sigma_sidis_units": SIGMA_SIDIS_UNITS,
        })
        try:
            if source.get("validation_status") != "VALID":
                row.update(Model_status="SKIPPED", Model_reason=(
                    f"manifest_status={source.get('validation_status','MISSING')}"
                ))
                rows.append(row)
                continue
            ebeam = float(source["ebeam"])
            x, q2, z, theta = (float(source[name]) for name in
                               ("x", "Q2", "z", "theta_pq"))
            derived = derive_kinematics(ebeam, x, q2, z, theta)
            target_a, target_z, charge = calculator_identity(key[3], key[2])
            row.update({
                "Ebeam_GeV": ebeam, "x": x, "Q2_GeV2": q2, "z": z,
                "theta_pq_deg": theta, "nu_GeV": derived["nu"],
                "p_pion_GeV": derived["ppion"], "pT_GeV": derived["pt"],
                "Target_A": target_a, "Target_Z": target_z,
                "Pion_charge": charge,
            })
            if not executable.is_file() or not os.access(executable, os.X_OK):
                row.update(Model_status="PENDING", Model_reason="CALCULATOR_MISSING")
            else:
                phi, multiplicity, sigma = run_calculator(
                    executable, (ebeam, q2, x, z, derived["pt"], target_a,
                                 target_z, charge, 1.0),
                )
                if multiplicity <= 0 or sigma <= 0:
                    raise ValueError("calculator returned nonpositive model values")
                row.update(Phi_model_rad=phi, M_sidis_model=multiplicity,
                           sigma_sidis_model=sigma, Model_status="OK")
        except (KeyError, ValueError, OSError, RuntimeError) as error:
            row.update(Model_status="ERROR", Model_reason=str(error))
        rows.append(row)
    return rows


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    jlab = root.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=(
        jlab / "JLab_upstream_forks/simc_gfortran/util/batch_job/simc_job_list_pass1.csv"
    ))
    parser.add_argument("--calculator", type=Path,
                        default=jlab / "simc_gfortran/util/sidisxsec/calc_semi_xsec")
    parser.add_argument("--output", type=Path, default=(
        root / "results/Pass1/Tables/RP_sidis_model/RP_sidis_model.csv"
    ))
    args = parser.parse_args()
    rows = build_catalog(args.manifest, args.calculator)
    atomic_csv(args.output, rows)
    counts = {status: sum(row["Model_status"] == status for row in rows)
              for status in ("OK", "PENDING", "SKIPPED", "ERROR")}
    print(f"Wrote {len(rows)} model rows to {args.output}")
    print("Status:", ", ".join(f"{key}={value}" for key, value in counts.items() if value))
    return 1 if counts["ERROR"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
