import math
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_sidis_azimuthal_analysis as azimuth  # noqa: E402


def curve(n, normalization=2.0, amplitude=0.25, second=-0.1):
    rows=[]
    for index in range(1, n + 1):
        phi=(index - 0.5) * math.pi / 3.0
        value=normalization * (1 + amplitude * math.cos(phi) + second * math.cos(2 * phi))
        rows.append({"Used_in_fit":1, "phipq_bin":index, "phipq_center":phi,
                     "Extracted_value":value, "Extracted_err":0.1})
    return rows


class SidisAzimuthalAnalysisTests(unittest.TestCase):
    def test_phi_canonicalization_is_half_open(self):
        self.assertEqual(azimuth.canonical_phi(2 * math.pi), 0.0)
        self.assertAlmostEqual(azimuth.canonical_phi(-math.pi / 3), 5 * math.pi / 3)
        self.assertEqual(azimuth.bin_index(0.0, azimuth.PHI_EDGES), 1)
        self.assertIsNone(azimuth.bin_index(2 * math.pi, azimuth.PHI_EDGES))

    def test_only_final_pt_bin_closes_its_upper_edge(self):
        self.assertEqual(azimuth.bin_index(0.1, azimuth.PT_EDGES, True), 2)
        self.assertEqual(azimuth.bin_index(0.5, azimuth.PT_EDGES, True), 5)
        self.assertIsNone(azimuth.bin_index(0.5, azimuth.PT_EDGES, False))

    def test_final_analysis_uses_shared_fifteen_by_five_grid(self):
        self.assertEqual(len(azimuth.PHI_EDGES)-1,15)
        self.assertEqual(len(azimuth.PT_EDGES)-1,5)

    def test_model_uses_all_fifteen_calc_semi_xsec_phi_rows(self):
        centers = [
            (azimuth.PHI_EDGES[index] + azimuth.PHI_EDGES[index + 1]) / 2
            for index in range(15)
        ]
        returned = [
            (phi, 10.0 + index, 20.0 + index)
            for index, phi in enumerate(centers)
        ]
        model = {
            "Target": "LH2", "Run_type": "pip", "Ebeam_GeV": "10.6",
            "Q2_GeV2": "3.1", "x": "0.31", "z": "0.5",
        }
        with patch.object(azimuth, "run_calculator_rows", return_value=returned) as call:
            values = azimuth.model_values(Path("calc_semi_xsec"), model, 0.15)
        self.assertEqual(len(values), 15)
        self.assertEqual(call.call_args.args[1][-1], 15.0)
        self.assertEqual(values[0], (centers[0], 20.0, 10.0))
        self.assertEqual(values[-1], (centers[-1], 34.0, 24.0))

    def test_two_points_are_not_fitted(self):
        result=azimuth.fit_curve(curve(2), "sigma")
        self.assertEqual(result["Fit_status"], "NOT_FITTED")
        self.assertEqual(result["Fit_line_style"], "none")

    def test_three_points_give_dashed_ndf_zero_fit(self):
        result=azimuth.fit_curve(curve(3), "sigma")
        self.assertEqual(result["NDF"], 0)
        self.assertEqual(result["Fit_line_style"], "dashed")
        self.assertAlmostEqual(result["C_0"], 2.0)
        self.assertAlmostEqual(result["C_1"], 0.5)
        self.assertAlmostEqual(result["C_2"], -0.2)

    def test_four_points_give_solid_multiplicity_fit(self):
        result=azimuth.fit_curve(curve(4), "multiplicity")
        self.assertEqual(result["Fit_line_style"], "solid")
        self.assertAlmostEqual(result["M_0"], 2.0)
        self.assertAlmostEqual(result["M_1"], 0.5)
        self.assertAlmostEqual(result["M_2"], -0.2)

    def test_approved_z_palette(self):
        self.assertEqual(azimuth.z_color(.36), azimuth.ROOT.kRed + 1)
        self.assertEqual(azimuth.z_color(.50), azimuth.z_color(.52))
        self.assertEqual(azimuth.z_color(.67), azimuth.ROOT.kBlue + 1)
        self.assertEqual(azimuth.z_color(.90), azimuth.ROOT.kBlack)

    def test_rejection_codes_are_compact_and_ordered(self):
        reason="WEIGHTED_DESIGN_CONDITION_ABOVE_LIMIT;INSUFFICIENT_ANGULAR_SPAN"
        self.assertEqual(azimuth.rejection_codes(reason),[1,3])


if __name__ == "__main__":
    unittest.main()
