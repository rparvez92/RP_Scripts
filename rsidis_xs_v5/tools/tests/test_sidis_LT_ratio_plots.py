import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[1];sys.path.insert(0,str(TOOLS))
import RP_sidis_LT_ratio_plots as stage  # noqa:E402

class LTRatioPlotTests(unittest.TestCase):
    def test_build_points_keeps_only_finite_separations(self):
        base={"Separation_status":"OK","R_sigmaL_over_sigmaT":".5","R_sigmaL_over_sigmaT_err":".1","Pt_low":".1","Pt_high":".2"}
        rows=[dict(base),dict(base,Separation_status="NOT_SEPARATED"),dict(base,R_sigmaL_over_sigmaT="nan")]
        points=stage.build_points(rows);self.assertEqual(len(points),1);self.assertAlmostEqual(points[0]["Pt2_low"],.01);self.assertAlmostEqual(points[0]["Pt2_high"],.04);self.assertAlmostEqual(points[0]["Pt2_center"],.025)

    def test_axis_range_contains_values_and_errors(self):
        low,high=stage.axis_range([1.,2.],[.2,.3]);self.assertLess(low,.8);self.assertGreater(high,2.3)

if __name__=="__main__":unittest.main()
