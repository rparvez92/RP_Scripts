import math
import sys
import unittest
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(TOOLS))

import RP_sidis_coefficient_plots as coefficients  # noqa: E402


class SidisCoefficientPlotTests(unittest.TestCase):
    def test_ratio_uncertainty_uses_covariance(self):
        ratio,error=coefficients.ratio_with_covariance(2,.2,4,.4,.04)
        self.assertAlmostEqual(ratio,.5)
        expected=.2**2/4**2+2**2*.4**2/4**4-2*2*.04/4**3
        self.assertAlmostEqual(error,math.sqrt(expected))

    def test_full_coverage_fit_is_plottable(self):
        row={"Pass":"pass5","Coverage_curve_status":"APPROVED","N_points_fitted":"15",
             "Coverage_curve_n_selected":"15","Coverage_curve_min_selected_span_rad":str(2*math.pi),"Fit_status":"WARNING"}
        self.assertTrue(coefficients.plottable(row)); row["Coverage_curve_n_selected"]="14"
        self.assertFalse(coefficients.plottable(row))

    def test_pdf_filename_contains_observable_without_subdirectory(self):
        name=coefficients.family_name(("phase1","pass5","pim","LH2","0.25","3.3"),"sigma")
        self.assertEqual(name,"phase1_pass5_pim_LH2_x0p25_Q23p3_sigma.pdf")

    def test_z_point_nine_is_not_displayed(self):
        self.assertFalse(coefficients.displayed_z({"z":"0.90"}))
        self.assertTrue(coefficients.displayed_z({"z":"0.67"}))


if __name__=="__main__": unittest.main()
