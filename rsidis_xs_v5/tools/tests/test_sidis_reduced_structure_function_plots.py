import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[1];sys.path.insert(0,str(TOOLS))
import RP_sidis_reduced_structure_function_plots as stage  # noqa:E402

class StructureFunctionPlotTests(unittest.TestCase):
    def test_comparison_keeps_shared_data_pull_undefined(self):
        global_row={"Epsilon_family_id":"f","Run_type":"pip","Target":"LH2","x":".25","Q2":"3.3","z":".5","Pt_bin":"1","Pt_low":"0","Pt_high":".1","Pt_high_closed":"0","Fit_status":"OK","Fit_reason":"","F_T":"2","F_T_err":".2","F_L":".5","F_L_err":".1",**{f"C_{i}_global":str(i) for i in range(4)},**{f"C_{i}_global_err":".1" for i in range(4)}}
        anchor={"Epsilon_family_id":"f","Pt_bin":"1","Separation_status":"OK","Separation_reason":"","F_T":"1.8","F_T_err":".3","F_L":".6","F_L_err":".2"}
        row=stage.comparison_rows([global_row],[anchor])[0];self.assertEqual(row["Comparison_status"],"AVAILABLE");self.assertAlmostEqual(row["Delta_F_T_global_minus_anchored"],.2);self.assertIn("NO_PULL",row["Comparison_note"]);self.assertNotIn("Pull",row)

    def test_unavailable_method_disables_comparison(self):
        global_row={"Epsilon_family_id":"f","Run_type":"pip","Target":"LH2","x":".25","Q2":"3.3","z":".5","Pt_bin":"1","Pt_low":"0","Pt_high":".1","Pt_high_closed":"0","Fit_status":"NOT_FITTED","Fit_reason":"rank"}
        row=stage.comparison_rows([global_row],[])[0];self.assertEqual(row["Comparison_status"],"UNAVAILABLE")

if __name__=="__main__":unittest.main()
