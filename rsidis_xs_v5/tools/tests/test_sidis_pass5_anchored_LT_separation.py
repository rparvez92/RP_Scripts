import math,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[1];sys.path.insert(0,str(TOOLS))
import RP_sidis_pass5_anchored_LT_separation as stage  # noqa:E402

def anchor(sigma_t=2.,sigma_l=.5,e=.8):
    row={field:"test" for field in stage.IDENTITY};row.update(Pt_bin="1",Pt_low="0",Pt_high=".1",Pt_high_closed="0",epsilon=str(e),C_0_pass5=str(sigma_t+e*sigma_l),C_0_pass5_err=".1",C_1_pass5=".2",C_2_pass5=".1",Cov_00=".01",Cov_01=".001",Cov_02="-.0005",Cov_11=".0025",Cov_12=".0002",Cov_22=".0016",Fit_status="OK",Fit_reason="",Combination_group_id="g5",Convention_version="RP_VIRTUAL_PHOTON_XSEC_HAND_V1",Cross_section_level="BORN",Virtual_photon_xsec_units="ub/GeV/sr",Coefficient_units="ub/GeV/sr");return row

def lower(pas,e,sigma_t=2.,sigma_l=.5,g1=.1,g2=-.05):
    return {"Member_pass":pas,"epsilon":str(e),"C_0":str(sigma_t+e*sigma_l),"C_0_err":".12","dC0_dC1_pass5":str(g1),"dC0_dC2_pass5":str(g2),"Fit_status":"OK","Fit_reason":""}

class AnchoredLTTests(unittest.TestCase):
    def test_two_point_separation_recovers_sigma_t_sigma_l(self):
        result=stage.separate(anchor(),lower("pass4",.4))
        self.assertEqual(result["Separation_status"],"OK");self.assertAlmostEqual(result["sigma_T"],2.);self.assertAlmostEqual(result["sigma_L"],.5);self.assertGreater(result["sigma_T_err"],0);self.assertGreater(result["sigma_L_err"],0)
        self.assertAlmostEqual(result["R_sigmaL_over_sigmaT"],.25);self.assertGreater(result["R_sigmaL_over_sigmaT_err"],0)
        self.assertEqual(result["Response_units"],"ub/GeV/sr");self.assertEqual(result["Cross_section_level"],"BORN")

    def test_pass3_is_validation_not_separation_input(self):
        result=stage.separate(anchor(),lower("pass4",.4),lower("pass3",.2))
        self.assertEqual(result["Pass3_validation_status"],"CONSISTENT");self.assertAlmostEqual(result["Pass3_validation_residual"],0,places=14)

    def test_pass3_tension_is_warning(self):
        p3=lower("pass3",.2);p3["C_0"]="4"
        result=stage.separate(anchor(),lower("pass4",.4),p3)
        self.assertEqual(result["Pass3_validation_status"],"TENSION");self.assertEqual(result["Separation_status"],"WARNING")

    def test_missing_pass4_is_not_separated(self):
        result=stage.separate(anchor(),None);self.assertEqual(result["Separation_status"],"NOT_SEPARATED");self.assertEqual(result["Separation_reason"],"PASS4_C0_UNAVAILABLE")

if __name__=="__main__":unittest.main()
