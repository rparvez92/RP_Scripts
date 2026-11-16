import math,sys,unittest
from pathlib import Path
import numpy as np
TOOLS=Path(__file__).resolve().parents[1];sys.path.insert(0,str(TOOLS))
import RP_fit_sidis_reduced_xsec_lower_passes as stage  # noqa:E402

class LowerPassFitTests(unittest.TestCase):
    def test_epsilon_scaling(self):
        s1,s2=stage.harmonic_scales(.4,.8)
        self.assertAlmostEqual(s1,math.sqrt(.4*1.4/(.8*1.8)));self.assertAlmostEqual(s2,.5)

    def test_gls_recovers_c0_with_anchor_covariance(self):
        eps,eps5=.4,.8;s1,s2=stage.harmonic_scales(eps,eps5);c0,c1,c2=2.,.3,-.2
        anchor={"Fit_status":"OK","Fit_reason":"","epsilon":str(eps5),"C_1_pass5":str(c1),"C_2_pass5":str(c2),"Cov_11":".0004","Cov_12":"0","Cov_22":".0004","Combination_group_id":"g5"}
        rows=[]
        for i,(low,high) in enumerate(zip(stage.PHI_EDGES[:-1],stage.PHI_EDGES[1:]),1):
            phi=.5*(low+high);value=c0+s1*c1*math.cos(phi)+s2*c2*math.cos(2*phi)
            row={field:"test" for field in (*stage.IDENTITY,*stage.PROVENANCE)};row.update(Member_phase="phase2",Member_pass="pass4",epsilon=str(eps),Pt_bin="1",Pt_low="0",Pt_high=".1",Pt_high_closed="0",phipq_bin=str(i),phipq_center=str(phi),Reduced_data_xsec=str(value),Reduced_data_xsec_err=".05",Used_in_lower_fit=1)
            rows.append(row)
        fit=stage.fit_lower_curve(rows,anchor);self.assertEqual(fit["Fit_status"],"OK");self.assertAlmostEqual(fit["C_0"],c0,places=12);self.assertGreater(fit["C_0_err"],0)
        self.assertAlmostEqual(fit["C_0_err"]**2,fit["C_0_lower_stat_err"]**2+fit["C_0_anchor_err"]**2,places=14)
        self.assertTrue(math.isfinite(fit["dC0_dC1_pass5"]));self.assertTrue(math.isfinite(fit["dC0_dC2_pass5"]))

    def test_point_policy_and_missing_anchor(self):
        source={"Y_Data":"1","Y_Data_err":".1","Y_MC":"1","Y_MC_err":".1","Data_by_MC":"1","Data_by_MC_err":".1","Reduced_data_xsec":"1","Reduced_data_xsec_err":".1","Reduction_status":"OK"}
        self.assertEqual(stage.classify_point(source,None)["Used_in_lower_fit"],0)
        anchor={"Fit_status":"OK","Fit_reason":""};self.assertEqual(stage.classify_point(source,anchor)["Used_in_lower_fit"],1)
        source["Data_by_MC_err"]=".4";self.assertEqual(stage.classify_point(source,anchor)["Used_in_lower_fit"],0)

if __name__=="__main__":unittest.main()
