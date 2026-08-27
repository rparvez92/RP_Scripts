import csv
import importlib
import math
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
overview = importlib.import_module("RP_plot_sidis_setting_overview")


def model(setting, z, theta, pt, status="OK"):
    return {"Phase":"phase1","Pass":"pass4","Run_type":"PIMINUS","Target":"LD2",
            "Setting":setting,"Ebeam_GeV":"8.5831","x":"0.25","Q2_GeV2":"3.3",
            "z":str(z),"theta_pq_deg":str(theta),"nu_GeV":"7.0",
            "p_pion_GeV":str(pt/math.sin(math.radians(theta))),"pT_GeV":str(pt),
            "Model_status":status}


def extraction(setting, sigma, multiplicity, status="OK"):
    return {"Phase":"phase1","Pass":"pass4","Run_type":"PIMINUS","Target":"LD2",
            "Setting":setting,"sigma_sidis_model":str(sigma/2),"sigma_sidis_data":str(sigma),
            "sigma_sidis_data_err":"0.1","sigma_sidis_units":"ub/GeV^2/sr^2",
            "M_sidis_model":str(multiplicity/2),"M_sidis_data":str(multiplicity),
            "M_sidis_data_err":"0.01","M_sidis_units":"GeV^-2","Y_Data":"1",
            "Y_Data_err":"0.1","Y_MC":"2","Y_MC_err":"0.2","C_Y":"0.5",
            "C_Y_err":"0.1","Extraction_status":status,"Extraction_reason":""}


class SettingOverviewTests(unittest.TestCase):
    def test_pt2_uses_degree_to_radian_conversion(self):
        self.assertAlmostEqual(overview.derive_pt2(3.0, 30.0), 2.25)

    def test_exact_join_and_upstream_values_are_unchanged(self):
        models=[model("a",.5,2,.12)]
        extracted=[extraction("a",4.5,.3)]
        rows=overview.join_catalogs(extracted,models,Path("e.csv"),Path("m.csv"))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]["sigma_sidis_data"],"4.5")
        self.assertEqual(rows[0]["M_sidis_data"],"0.3")
        self.assertAlmostEqual(float(rows[0]["pT2_GeV2"]),.12**2)

    def test_missing_exact_model_identity_is_error(self):
        rows=overview.join_catalogs([extraction("missing",1,.1)],[],Path("e"),Path("m"))
        self.assertEqual(rows[0]["Overview_status"],"ERROR")
        self.assertEqual(rows[0]["Overview_reason"],"MODEL_ROW_MISSING")

    def test_controlled_z_and_pt2_groups(self):
        models=[model("z1",.36,2,.09),model("z2",.50,2,.12),
                model("angle",.50,5.2,.31),model("single",.67,8.5,.4)]
        rows=overview.join_catalogs([extraction(m["Setting"],2+i,.2+i/10) for i,m in enumerate(models)],models,Path("e"),Path("m"))
        zgroups=overview.controlled_groups(rows,"z")
        ptgroups=overview.controlled_groups(rows,"pT2")
        self.assertEqual(len(zgroups),1)
        self.assertEqual([float(r["z"]) for r in zgroups[0][3]],[.36,.50])
        self.assertEqual(len(ptgroups),1)
        self.assertEqual(len(ptgroups[0][3]),2)

    def test_pending_and_single_points_remain_catalog_only(self):
        models=[model("ok",.5,2,.12),model("pending",.67,2,.16)]
        rows=overview.join_catalogs([extraction("ok",2,.2),extraction("pending",2,.2,"PENDING")],models,Path("e"),Path("m"))
        self.assertEqual(len(rows),2)
        self.assertEqual(overview.controlled_groups(rows,"z"),[])

    def test_plot_geometry_guards_axis_titles(self):
        self.assertEqual((overview.CANVAS_WIDTH,overview.CANVAS_HEIGHT),(1800,1000))
        self.assertGreaterEqual(overview.PAD_LEFT_MARGIN,.22)
        self.assertGreaterEqual(overview.PAD_RIGHT_MARGIN,.06)

    def test_no_event_or_full_cut_interface(self):
        names=" ".join(overview.FIELDS)
        self.assertNotIn("Bin_",names); self.assertNotIn("Full",names)
        self.assertNotIn("Purity",names); self.assertNotIn("Stability",names)


if __name__ == "__main__": unittest.main()
