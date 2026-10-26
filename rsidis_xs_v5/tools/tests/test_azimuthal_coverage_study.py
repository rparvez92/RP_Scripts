import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_azimuthal_coverage_study as coverage  # noqa: E402


class AzimuthalCoverageStudyTests(unittest.TestCase):
    def test_longest_run_is_circular(self):
        mask = [True, True, False, False, True, True, True]
        self.assertEqual(coverage.longest_circular_run(mask), [4, 5, 6, 0, 1])

    def test_tie_uses_first_run_in_bin_order(self):
        self.assertEqual(
            coverage.longest_circular_run([True, True, False, True, True, False]),
            [0, 1],
        )

    def test_full_coverage_has_rank_three(self):
        rank, condition = coverage.design_diagnostics(
            range(len(coverage.PHI_EDGES) - 1)
        )
        self.assertEqual(rank, 3)
        self.assertTrue(math.isfinite(condition))

    def test_near_zero_data_fails_relative_precision(self):
        result = coverage.classify_bin(0.0, 1.0, 10.0, 1.0, [], [], .3, .5, 0, 2)
        self.assertEqual(result["Ratio_precision_pass"], 0)
        self.assertIn("NONPOSITIVE_CORRECTED_DATA", result["Bin_reason"])

    def test_ratio_sanity_is_warning_not_rejection(self):
        result = coverage.classify_bin(25.0, 2.0, 10.0, 1.0, [], [], .3, .5, 0, 2)
        self.assertEqual(result["Independent_pass"], 1)
        self.assertEqual(result["Bin_status"], "WARNING")

    def test_zero_data_uncertainty_is_governed_by_ratio_precision(self):
        result = coverage.classify_bin(10.0, 0.0, 10.0, 1.0, [], [], .3, .5, 0, 2)
        self.assertEqual(result["Independent_pass"], 1)
        self.assertNotIn("DATA_UNCERTAINTY", result["Bin_reason"])

    def test_after_phi_projection_sums_only_surviving_pt_cells(self):
        rows = [
            {"phipq_bin": 1, "Survives_all_requirements": 1,
             "Y_Data": 2, "Y_Data_err": 1, "Y_MC": 4, "Y_MC_err": 2},
            {"phipq_bin": 1, "Survives_all_requirements": 0,
             "Y_Data": 20, "Y_Data_err": 3, "Y_MC": 40, "Y_MC_err": 4},
        ]
        data, mc = coverage.phi_projection(rows, True, "test_after_projection")
        self.assertEqual(data.GetBinContent(1), 2)
        self.assertEqual(data.GetBinError(1), 1)
        self.assertEqual(mc.GetBinContent(1), 4)
        self.assertEqual(mc.GetBinError(1), 2)

    def test_missing_optional_correction_is_warning_not_rejection(self):
        result = coverage.classify_bin(
            10.0, 1.0, 10.0, 1.0, [], ["ZERO_CORRECTION_CLASS=target_p"],
            .3, .5, 0, 2,
        )
        self.assertEqual(result["Independent_pass"], 1)
        self.assertEqual(result["Bin_status"], "WARNING")

    def test_combined_curve_survival_requires_curve_approval(self):
        rows=[]
        for index,(low,high) in enumerate(zip(coverage.PHI_EDGES[:-1],coverage.PHI_EDGES[1:]),1):
            rows.append({"Independent_pass":1,"phipq_center":(low+high)/2,
                         "Data_by_MC_err":.1,"Selected_contiguous":0,
                         "Survives_all_requirements":0})
        common={field:"" for field in coverage.COMBINED_CURVE_FIELDS}
        common.update({"Combination_group_id":"g","Phase":"phase1","Pass":"pass5",
                       "Run_type":"pim","Target":"LH2","x":.25,"Q2":3.3,"z":.36,
                       "Pt_bin":1,"Pt_low":0,"Pt_high":.1,"Pt_high_closed":0,
                       "Binning_version":"v","Policy_version":"p"})
        for row in rows: row.update(common)
        curve=coverage.finalize_combined_curve(rows,math.pi,100)
        self.assertEqual(curve["Curve_status"],"APPROVED")
        self.assertEqual(curve["N_surviving"],15)
        self.assertTrue(all(row["Survives_all_requirements"]==1 for row in rows))


if __name__ == "__main__":
    unittest.main()
