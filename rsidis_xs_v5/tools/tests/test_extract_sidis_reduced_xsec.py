import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_extract_sidis_reduced_xsec as stage  # noqa: E402
from RP_sidis_reduced_kinematics import calculate_reduced_kinematics  # noqa: E402


def member():
    kin = calculate_reduced_kinematics(10.6759, 0.44, 4.4, 0.52, 0.0)
    return {
        "Epsilon_family_id": "EF_test", "Epsilon_family_policy_version": "test",
        "Family_status": "LT_READY_PASS3_VALIDATION", "Family_reason": "",
        "Run_type": "pip", "Target": "LH2", "x": "0.44", "Q2": "4.4",
        "z": "0.52", "N_members": "3", "Pass_pattern": "pass3|pass4|pass5",
        "Has_pass3": "1", "Has_pass4": "1", "Has_pass5": "1",
        "N_distinct_epsilon": "3", "Epsilon_min": "0.2",
        "Epsilon_max": "0.8", "Epsilon_span": "0.6", "Member_index": "3",
        "Member_phase": "phase2", "Member_pass": "pass5",
        "Combination_group_id": "test_group", "Ebeam_GeV": "10.6759",
        "epsilon": str(kin.epsilon),
    }


def source_rows(status="OK", reason=""):
    rows = []
    for pt_bin, (pt_low, pt_high) in enumerate(zip(stage.PT_EDGES[:-1], stage.PT_EDGES[1:]), 1):
        for phi_bin, (phi_low, phi_high) in enumerate(zip(stage.PHI_EDGES[:-1], stage.PHI_EDGES[1:]), 1):
            rows.append({
                "Combination_group_id": "test_group",
                "Combination_policy_version": "test", "Phase": "phase2",
                "Pass": "pass5", "Run_type": "pip", "Target": "LH2",
                "x": "0.44", "Q2": "4.4", "z": "0.52",
                "Ebeam_GeV": "10.6759", "N_contributors": "1",
                "Contributing_setting_ids": "RPTEST", "Contributing_thpq_deg": "2",
                "Pt_bin": str(pt_bin), "Pt_low": str(pt_low),
                "Pt_high": str(pt_high), "Pt_high_closed": str(int(pt_bin == 5)),
                "phipq_bin": str(phi_bin), "phipq_low": str(phi_low),
                "phipq_high": str(phi_high),
                "phipq_center": str(0.5 * (phi_low + phi_high)),
                "Binning_version": stage.AZIMUTHAL_BINNING_VERSION,
                "Y_Data": "8", "Y_Data_err": "0.8", "Y_MC": "10",
                "Y_MC_err": "1", "Data_by_MC": "0.8", "Data_by_MC_err": "0.1",
                "Combination_status": status, "Combination_reason": reason,
            })
    return rows


def model_loader(_calculator, _member, _pt):
    return [(0.5 * (low + high), 2.0, 0.25)
            for low, high in zip(stage.PHI_EDGES[:-1], stage.PHI_EDGES[1:])]


class ExtractReducedXsecTests(unittest.TestCase):
    def test_complete_member_reduces_all_points(self):
        points = stage.extract_member(member(), source_rows(), Path("unused"), model_loader)
        self.assertEqual(len(points), 75)
        first = points[0]
        self.assertEqual(first["Reduction_status"], "OK")
        self.assertEqual(first["Cross_section_level"], "BORN")
        self.assertAlmostEqual(first["Extracted_lab_xsec_dEpi"], 1.6)
        self.assertAlmostEqual(
            first["Extracted_lab_xsec_dPpi"],
            1.6 * first["dEpi_dPpi"],
        )
        self.assertAlmostEqual(
            first["Virtual_photon_data_xsec"],
            first["Extracted_lab_xsec_dPpi"] / first["Gamma_Hand_per_GeV"],
        )
        self.assertAlmostEqual(
            first["Virtual_photon_data_xsec_err"],
            0.2 * first["Virtual_photon_reduction_factor_GeV"],
        )
        self.assertEqual(first["Virtual_photon_xsec_units"], "ub/GeV/sr")
        self.assertAlmostEqual(first["Lab_reconstruction_residual"], 0.0)

    def test_upstream_warning_remains_warning(self):
        points = stage.extract_member(
            member(), source_rows("WARNING", "ADVISORY"), Path("unused"), model_loader
        )
        self.assertTrue(all(row["Reduction_status"] == "WARNING" for row in points))
        self.assertTrue(all("ADVISORY" in row["Reduction_reason"] for row in points))

    def test_ld2_uses_the_same_central_proton_mass_convention(self):
        ld2_member = member(); ld2_member["Target"] = "LD2"
        ld2_sources = source_rows()
        for row in ld2_sources:
            row["Target"] = "LD2"
        points = stage.extract_member(
            ld2_member, ld2_sources, Path("unused"), model_loader
        )
        self.assertTrue(all(point["Target_mass_GeV"] == points[0]["Target_mass_GeV"]
                            for point in points))
        self.assertEqual(
            points[0]["Target_mass_convention"],
            "CENTRAL_STATIONARY_PROTON_MASS_FOR_LH2_AND_LD2",
        )
        self.assertEqual(
            points[0]["Radiative_treatment"],
            "RADIATED_MC_YIELD_WITH_CENTRAL_BORN_MODEL",
        )

    def test_upstream_structural_error_remains_error(self):
        points = stage.extract_member(
            member(), source_rows("ERROR", "MISSING_INPUT"), Path("unused"), model_loader
        )
        self.assertTrue(all(row["Reduction_status"] == "ERROR" for row in points))
        self.assertTrue(all("UPSTREAM_COMBINATION_ERROR" in row["Reduction_reason"]
                            for row in points))

    def test_nonfinite_ratio_is_invalid_not_computational_error(self):
        rows = source_rows(); rows[0]["Data_by_MC"] = "nan"
        points = stage.extract_member(member(), rows, Path("unused"), model_loader)
        self.assertEqual(points[0]["Reduction_status"], "INVALID")
        self.assertIn("NONFINITE_VIRTUAL_PHOTON_RESULT", points[0]["Reduction_reason"])

    def test_family_source_mismatch_fails_closed(self):
        bad = member(); bad["Member_phase"] = "phase1"
        with self.assertRaisesRegex(ValueError, "mismatch in Phase"):
            stage.extract_member(bad, source_rows(), Path("unused"), model_loader)

    def test_incomplete_grid_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "expected 75"):
            stage.extract_member(member(), source_rows()[:-1], Path("unused"), model_loader)


if __name__ == "__main__":
    unittest.main()
