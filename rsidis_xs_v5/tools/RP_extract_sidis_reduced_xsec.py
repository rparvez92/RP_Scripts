#!/usr/bin/env python3
"""Extract table-only Born SIDIS virtual-photon cross sections."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Mapping, Sequence

from RP_build_sidis_model import calculator_identity, run_calculator_rows
from RP_sidis_binning import (
    AZIMUTHAL_BINNING_VERSION,
    AZIMUTHAL_PHI_EDGES as PHI_EDGES,
    AZIMUTHAL_PT_EDGES as PT_EDGES,
)
from RP_sidis_reduced_kinematics import (
    CONVENTION_VERSION,
    CROSS_SECTION_LEVEL,
    RADIATIVE_TREATMENT,
    TARGET_MASS_CONVENTION,
    calculate_reduced_kinematics,
)

POLICY_VERSION = "RP_EXTRACT_SIDIS_VIRTUAL_PHOTON_XSEC_V2"
MODEL_POLICY = "CALC_SEMI_XSEC_AT_15_PHI_BIN_CENTERS"
FAMILY_FIELDS = (
    "Epsilon_family_id", "Epsilon_family_policy_version", "Family_status",
    "Family_reason", "Run_type", "Target", "x", "Q2", "z", "N_members",
    "Pass_pattern", "Has_pass3", "Has_pass4", "Has_pass5",
    "N_distinct_epsilon", "Epsilon_min", "Epsilon_max", "Epsilon_span",
    "Member_index", "Member_phase", "Member_pass",
)
PROVENANCE_FIELDS = (
    "Combination_group_id", "Combination_policy_version", "Ebeam_GeV",
    "N_contributors", "Contributing_setting_ids", "Contributing_thpq_deg",
)
BIN_FIELDS = (
    "Pt_bin", "Pt_low", "Pt_high", "Pt_high_closed", "phipq_bin",
    "phipq_low", "phipq_high", "phipq_center", "Binning_version",
)
SOURCE_FIELDS = (
    "Y_Data", "Y_Data_err", "N_Data_events", "N_MC_events",
    "Y_MC_sidis", "Y_MC_sidis_err", "Y_MC_rho", "Y_MC_rho_err",
    "Y_MC_delta", "Y_MC_delta_err", "Y_MC_exclusive",
    "Y_MC_exclusive_err", "Y_MC", "Y_MC_err", "Data_by_MC",
    "Data_by_MC_err", "Combination_status", "Combination_reason",
)
KINEMATIC_FIELDS = (
    "Reduced_xsec_policy_version", "Convention_version",
    "Target_mass_convention", "Cross_section_level", "Radiative_treatment",
    "Model_phi_policy", "Model_Pt_GeV", "Target_mass_GeV", "nu_GeV",
    "Eprime_GeV", "theta_e_rad", "y", "gamma", "epsilon", "W2_GeV2",
    "K_Hand_GeV", "Gamma_Hand_per_GeV", "Epi_GeV", "Ppi_GeV",
    "theta_pq_central_rad", "cos_theta_pq", "dEpi_dPpi",
    "SIMC_internal_hadronic_jacobian_GeV",
    "Virtual_photon_reduction_factor_GeV",
)
RESULT_FIELDS = (
    "SIMC_lab_xsec_dEpi", "Extracted_lab_xsec_dEpi",
    "Extracted_lab_xsec_dEpi_err", "Lab_xsec_dEpi_units",
    "SIMC_lab_xsec_dPpi", "Extracted_lab_xsec_dPpi",
    "Extracted_lab_xsec_dPpi_err", "Lab_xsec_dPpi_units",
    "SIMC_multiplicity", "Virtual_photon_model_xsec",
    "Virtual_photon_data_xsec", "Virtual_photon_data_xsec_err",
    "Virtual_photon_xsec_units", "Lab_reconstruction_residual",
    "Reduction_status", "Reduction_reason",
)
POINT_FIELDS = [*FAMILY_FIELDS, *PROVENANCE_FIELDS, *BIN_FIELDS,
                *SOURCE_FIELDS, *KINEMATIC_FIELDS, *RESULT_FIELDS]
SUMMARY_FIELDS = [
    *FAMILY_FIELDS, *PROVENANCE_FIELDS, "epsilon", "N_points", "N_OK",
    "N_WARNING", "N_INVALID", "N_ERROR", "Gamma_Hand_min_per_GeV",
    "Gamma_Hand_max_per_GeV", "dEpi_dPpi_min", "dEpi_dPpi_max",
    "Virtual_photon_reduction_factor_min_GeV",
    "Virtual_photon_reduction_factor_max_GeV",
    "Max_abs_lab_reconstruction_residual", "Summary_status", "Summary_reason",
]


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def number(value: object, default: float = math.nan) -> float:
    return float(value) if finite(value) else default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def atomic_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def calculator_model_values(calculator: Path, member: Mapping[str, object],
                            pt: float) -> list[tuple[float, float, float]]:
    """Return (phi, lab xsec, multiplicity) at all analysis-bin centers."""
    target_a, target_z, charge = calculator_identity(
        str(member["Target"]), str(member["Run_type"])
    )
    rows = run_calculator_rows(
        calculator,
        (number(member["Ebeam_GeV"]), number(member["Q2"]), number(member["x"]),
         number(member["z"]), pt, target_a, target_z, charge,
         float(len(PHI_EDGES) - 1)),
    )
    if len(rows) != len(PHI_EDGES) - 1:
        raise ValueError(f"calculator returned {len(rows)} phi rows; expected {len(PHI_EDGES)-1}")
    output = []
    for index, (phi, multiplicity, sigma) in enumerate(rows):
        expected = 0.5 * (PHI_EDGES[index] + PHI_EDGES[index + 1])
        if not math.isclose(phi, expected, rel_tol=0.0, abs_tol=1.0e-10):
            raise ValueError(f"calculator phi mismatch: {phi} != {expected}")
        if not finite(sigma) or not finite(multiplicity) or sigma <= 0 or multiplicity <= 0:
            raise ValueError("calculator returned a nonpositive or nonfinite model")
        output.append((phi, sigma, multiplicity))
    return output


def validate_sources(rows: Sequence[Mapping[str, object]], group_id: str
                     ) -> list[Mapping[str, object]]:
    expected = (len(PT_EDGES) - 1) * (len(PHI_EDGES) - 1)
    if len(rows) != expected:
        raise ValueError(f"{group_id} contains {len(rows)} rows; expected {expected}")
    ordered = sorted(rows, key=lambda row: (int(row["Pt_bin"]), int(row["phipq_bin"])))
    required = {(pt, phi) for pt in range(1, len(PT_EDGES))
                for phi in range(1, len(PHI_EDGES))}
    observed = {(int(row["Pt_bin"]), int(row["phipq_bin"])) for row in ordered}
    if observed != required:
        raise ValueError(f"{group_id} does not contain the complete Pt/phi grid")
    for field in ("Phase", "Pass", "Run_type", "Target", "x", "Q2", "z", "Ebeam_GeV"):
        if len({str(row[field]) for row in ordered}) != 1:
            raise ValueError(f"{group_id} has inconsistent {field}")
    return ordered


def reduction_decision(source: Mapping[str, object], reduced_data: float,
                       reduced_error: float) -> tuple[str, str]:
    reasons = set(filter(None, str(source.get("Combination_reason", "")).split(";")))
    upstream = str(source.get("Combination_status", ""))
    if upstream == "ERROR":
        reasons.add("UPSTREAM_COMBINATION_ERROR")
        return "ERROR", ";".join(sorted(reasons))
    if not finite(reduced_data) or not finite(reduced_error) or reduced_error < 0:
        reasons.add("NONFINITE_VIRTUAL_PHOTON_RESULT")
        return "INVALID", ";".join(sorted(reasons))
    if reduced_data <= 0:
        reasons.add("NONPOSITIVE_VIRTUAL_PHOTON_DATA")
    status = "WARNING" if upstream == "WARNING" or reasons else "OK"
    return status, ";".join(sorted(reasons))


def extract_member(member: Mapping[str, object],
                   source_rows: Sequence[Mapping[str, object]], calculator: Path,
                   model_loader: Callable = calculator_model_values
                   ) -> list[dict[str, object]]:
    group_id = str(member["Combination_group_id"])
    ordered = validate_sources(source_rows, group_id); first = ordered[0]
    checks = {
        "Member_phase": "Phase", "Member_pass": "Pass", "Run_type": "Run_type",
        "Target": "Target", "x": "x", "Q2": "Q2", "z": "z",
        "Ebeam_GeV": "Ebeam_GeV",
    }
    for member_field, source_field in checks.items():
        if str(member[member_field]) != str(first[source_field]):
            raise ValueError(f"{group_id} family/source mismatch in {source_field}")
    output = []
    for pt_bin, (pt_low, pt_high) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:]), 1):
        pt = 0.5 * (pt_low + pt_high)
        kin = calculate_reduced_kinematics(
            number(member["Ebeam_GeV"]), number(member["x"]),
            number(member["Q2"]), number(member["z"]), pt,
        )
        if not math.isclose(kin.epsilon, number(member["epsilon"]),
                            rel_tol=0.0, abs_tol=2.0e-12):
            raise ValueError(f"{group_id} epsilon-family/kernel mismatch")
        models = model_loader(calculator, member, pt)
        sources = [row for row in ordered if int(row["Pt_bin"]) == pt_bin]
        if len(models) != len(sources):
            raise ValueError(f"{group_id} model/source phi count mismatch")
        for source, (phi, sigma, multiplicity) in zip(sources, models):
            if not math.isclose(phi, number(source["phipq_center"]),
                                rel_tol=0.0, abs_tol=1.0e-10):
                raise ValueError(f"{group_id} model/source phi-center mismatch")
            ratio = number(source.get("Data_by_MC")); ratio_error = number(source.get("Data_by_MC_err"))
            lab_data = ratio * sigma if finite(ratio) else math.nan
            lab_error = ratio_error * sigma if finite(ratio_error) else math.nan
            dEpi_dPpi = kin.dEpi_dPpi
            reduction_factor = kin.virtual_photon_reduction_factor
            model_dPpi = sigma * dEpi_dPpi
            data_dPpi = lab_data * dEpi_dPpi
            error_dPpi = lab_error * dEpi_dPpi
            virtual_model = sigma * reduction_factor
            virtual_data = lab_data * reduction_factor
            virtual_error = lab_error * reduction_factor
            reconstructed = virtual_data / reduction_factor if finite(virtual_data) else math.nan
            residual = reconstructed - lab_data if finite(reconstructed) else math.nan
            status, reason = reduction_decision(source, virtual_data, virtual_error)
            row = {field: member.get(field, "") for field in FAMILY_FIELDS}
            row.update({field: source.get(field, "") for field in (
                *PROVENANCE_FIELDS, *BIN_FIELDS, *SOURCE_FIELDS,
            )})
            row.update({
                "Reduced_xsec_policy_version": POLICY_VERSION,
                "Convention_version": CONVENTION_VERSION,
                "Target_mass_convention": TARGET_MASS_CONVENTION,
                "Cross_section_level": CROSS_SECTION_LEVEL,
                "Radiative_treatment": RADIATIVE_TREATMENT,
                "Model_phi_policy": MODEL_POLICY, "Model_Pt_GeV": pt,
                "Target_mass_GeV": kin.target_mass,
                "nu_GeV": kin.nu, "Eprime_GeV": kin.eprime,
                "theta_e_rad": kin.theta_e, "y": kin.y, "gamma": kin.gamma,
                "epsilon": kin.epsilon, "W2_GeV2": kin.W2,
                "K_Hand_GeV": kin.hand_equivalent_energy,
                "Gamma_Hand_per_GeV": kin.hand_flux,
                "Epi_GeV": kin.pion_energy, "Ppi_GeV": kin.pion_momentum,
                "theta_pq_central_rad": kin.theta_pq,
                "cos_theta_pq": kin.cos_theta_pq,
                "dEpi_dPpi": dEpi_dPpi,
                "SIMC_internal_hadronic_jacobian_GeV": kin.simc_internal_hadronic_jacobian,
                "Virtual_photon_reduction_factor_GeV": reduction_factor,
                "SIMC_lab_xsec_dEpi": sigma,
                "Extracted_lab_xsec_dEpi": lab_data,
                "Extracted_lab_xsec_dEpi_err": lab_error,
                "Lab_xsec_dEpi_units": "ub/GeV2/sr2",
                "SIMC_lab_xsec_dPpi": model_dPpi,
                "Extracted_lab_xsec_dPpi": data_dPpi,
                "Extracted_lab_xsec_dPpi_err": error_dPpi,
                "Lab_xsec_dPpi_units": "ub/GeV2/sr2",
                "SIMC_multiplicity": multiplicity,
                "Virtual_photon_model_xsec": virtual_model,
                "Virtual_photon_data_xsec": virtual_data,
                "Virtual_photon_data_xsec_err": virtual_error,
                "Virtual_photon_xsec_units": "ub/GeV/sr",
                "Lab_reconstruction_residual": residual,
                "Reduction_status": status, "Reduction_reason": reason,
            })
            output.append(row)
    return output


def summarize_member(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    first = rows[0]; counts = Counter(str(row["Reduction_status"]) for row in rows)
    residuals = [abs(number(row["Lab_reconstruction_residual"])) for row in rows]
    max_residual = max(residuals)
    reasons = []
    if counts["ERROR"]: reasons.append("REDUCTION_ERRORS_PRESENT")
    if counts["INVALID"]: reasons.append("INVALID_POINTS_PRESENT")
    if counts["WARNING"]: reasons.append("REDUCTION_WARNINGS_PRESENT")
    if max_residual > 1.0e-12: reasons.append("LAB_RECONSTRUCTION_RESIDUAL")
    row = {field: first.get(field, "") for field in (*FAMILY_FIELDS, *PROVENANCE_FIELDS)}
    row.update({
        "epsilon": first["epsilon"], "N_points": len(rows),
        "N_OK": counts["OK"], "N_WARNING": counts["WARNING"],
        "N_INVALID": counts["INVALID"], "N_ERROR": counts["ERROR"],
        "Gamma_Hand_min_per_GeV": min(number(item["Gamma_Hand_per_GeV"]) for item in rows),
        "Gamma_Hand_max_per_GeV": max(number(item["Gamma_Hand_per_GeV"]) for item in rows),
        "dEpi_dPpi_min": min(number(item["dEpi_dPpi"]) for item in rows),
        "dEpi_dPpi_max": max(number(item["dEpi_dPpi"]) for item in rows),
        "Virtual_photon_reduction_factor_min_GeV": min(number(item["Virtual_photon_reduction_factor_GeV"]) for item in rows),
        "Virtual_photon_reduction_factor_max_GeV": max(number(item["Virtual_photon_reduction_factor_GeV"]) for item in rows),
        "Max_abs_lab_reconstruction_residual": max_residual,
        "Summary_status": "ERROR" if counts["ERROR"] else "WARNING" if reasons else "OK",
        "Summary_reason": ";".join(reasons),
    })
    return row


def main() -> int:
    project = Path(__file__).resolve().parents[1]; jlab = project.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family-members", type=Path, default=project / "results/Pass1/Tables/RP_sidis_epsilon_family_map/RP_sidis_epsilon_family_members.csv")
    parser.add_argument("--combined-yields", type=Path, default=project / "results/Pass1/Tables/RP_combine_thpq_yields/RP_combined_thpq_yields.csv")
    parser.add_argument("--calculator", type=Path, default=jlab / "simc_gfortran/util/sidisxsec/calc_semi_xsec")
    parser.add_argument("--output-dir", type=Path, default=project / "results/Pass1/Tables/RP_extract_sidis_reduced_xsec")
    args = parser.parse_args()
    family_members = read_csv(args.family_members)
    source_groups = defaultdict(list)
    for row in read_csv(args.combined_yields):
        source_groups[row["Combination_group_id"]].append(row)
    mapped = {row["Combination_group_id"] for row in family_members}
    if mapped != set(source_groups):
        raise ValueError("epsilon-family and combined-yield group sets do not match")
    points = []; summaries = []
    for index, member in enumerate(family_members, 1):
        found = extract_member(member, source_groups[str(member["Combination_group_id"])], args.calculator)
        points.extend(found); summaries.append(summarize_member(found))
        print(f"[{index}/{len(family_members)}] {summaries[-1]['Summary_status']}: {member['Combination_group_id']}")
    atomic_csv(args.output_dir / "RP_sidis_reduced_xsec_points.csv", POINT_FIELDS, points)
    atomic_csv(args.output_dir / "RP_sidis_reduced_xsec_points_Problematic.csv", POINT_FIELDS,
               [row for row in points if row["Reduction_status"] != "OK"])
    atomic_csv(args.output_dir / "RP_sidis_reduced_xsec_member_summary.csv", SUMMARY_FIELDS, summaries)
    print("Points:", len(points), Counter(row["Reduction_status"] for row in points))
    print("Members:", len(summaries), Counter(row["Summary_status"] for row in summaries))
    return 1 if any(row["Reduction_status"] == "ERROR" for row in points) else 0


if __name__ == "__main__":
    raise SystemExit(main())
