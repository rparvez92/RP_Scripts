import math,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[1];sys.path.insert(0,str(TOOLS))
import RP_sidis_global_reduced_xsec_fit as stage  # noqa:E402

class GlobalReducedFitTests(unittest.TestCase):
    def test_design_row_matches_physics_prefactors(self):
        row=stage.design_row(.5,.3);self.assertEqual(row[0],1);self.assertEqual(row[1],.5);self.assertAlmostEqual(row[2],math.sqrt(1.5)*math.cos(.3));self.assertAlmostEqual(row[3],.5*math.cos(.6))

    def test_global_fit_recovers_four_coefficients(self):
        coefficients=[2.,.4,.25,-.15];rows=[]
        for pas,epsilon in (("pass3",.3),("pass4",.55),("pass5",.8)):
            for i in range(15):
                phi=(i+.5)*2*math.pi/15;value=sum(a*b for a,b in zip(stage.design_row(epsilon,phi),coefficients));row={field:"test" for field in stage.IDENTITY};row.update(Member_pass=pas,epsilon=str(epsilon),Pt_bin="1",Pt_low="0",Pt_high=".1",Pt_high_closed="0",phipq_center=str(phi),Reduced_data_xsec=str(value),Reduced_data_xsec_err=".05",Global_point_usable=1);rows.append(row)
        fit=stage.fit_global(rows);self.assertEqual(fit["Fit_status"],"OK")
        for i,value in enumerate(coefficients):self.assertAlmostEqual(fit[f"C_{i}_global"],value,places=12)
        self.assertEqual(fit["Weighted_design_rank"],4)

    def test_rank_deficient_single_epsilon_is_not_fitted(self):
        rows=[]
        for i in range(15):
            phi=(i+.5)*2*math.pi/15;row={field:"test" for field in stage.IDENTITY};row.update(Member_pass="pass5",epsilon=".8",Pt_bin="1",Pt_low="0",Pt_high=".1",Pt_high_closed="0",phipq_center=str(phi),Reduced_data_xsec="1",Reduced_data_xsec_err=".1",Global_point_usable=1);rows.append(row)
        fit=stage.fit_global(rows);self.assertEqual(fit["Fit_status"],"NOT_FITTED");self.assertEqual(fit["Fit_reason"],"GLOBAL_DESIGN_RANK_BELOW_4")

    def test_global_point_selection_is_anchor_independent(self):
        source={field:"test" for field in stage.POINT_FIELDS};source.update(Y_Data="1",Y_Data_err=".1",Y_MC="1",Y_MC_err=".1",Data_by_MC="1",Data_by_MC_err=".1",Reduced_data_xsec="1",Reduced_data_xsec_err=".1",Reduction_status="OK")
        self.assertEqual(stage.normalize_point(source)["Global_point_usable"],1)
        source["Data_by_MC_err"]=".4";self.assertEqual(stage.normalize_point(source)["Global_point_usable"],0)

if __name__=="__main__":unittest.main()
