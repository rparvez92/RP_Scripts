import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_plot_sidis_setting_overview as overview  # noqa: E402


def model(setting, z, theta, pt):
    return {
        "Phase":"phase1", "Pass":"pass4", "Run_type":"pim", "Target":"LD2",
        "setting_id":setting, "Ebeam_GeV":"8.5831", "x":"0.25",
        "Q2_GeV2":"3.3", "z":str(z), "theta_pq_deg":str(theta),
        "nu_GeV":"7.0", "p_pion_GeV":str(pt/math.sin(math.radians(theta))),
        "pT_GeV":str(pt), "Model_status":"OK",
    }


def extraction(setting, sigma, multiplicity, status="OK"):
    return {
        "Phase":"phase1", "Pass":"pass4", "Run_type":"pim", "Target":"LD2",
        "setting_id":setting, "sigma_sidis_model":str(sigma/2),
        "sigma_sidis_data":str(sigma), "sigma_sidis_data_err":"0.1",
        "sigma_sidis_units":"ub/GeV^2/sr^2", "M_sidis_model":str(multiplicity/2),
        "M_sidis_data":str(multiplicity), "M_sidis_data_err":"0.01",
        "M_sidis_units":"GeV^-2", "Y_Data":"1", "Y_Data_err":"0.1",
        "Y_MC":"2", "Y_MC_err":"0.2", "C_Y":"0.5", "C_Y_err":"0.1",
        "Extraction_status":status, "Extraction_reason":"",
    }


class SettingOverviewTests(unittest.TestCase):
    def test_pt_uses_degree_to_radian_conversion_and_magnitude(self):
        self.assertAlmostEqual(overview.derive_pt(3.0, 30.0), 1.5)
        self.assertAlmostEqual(overview.derive_pt(3.0, -30.0), 1.5)

    def test_join_uses_setting_id_and_preserves_extraction(self):
        rows = overview.join_catalogs(
            [extraction("GT001", 4.5, .3)], [model("GT001", .5, 2, .12)],
            Path("e.csv"), Path("m.csv"),
        )
        self.assertEqual(rows[0]["setting_id"], "GT001")
        self.assertEqual(rows[0]["sigma_sidis_data"], "4.5")
        self.assertAlmostEqual(float(rows[0]["Pt_GeV"]), .12)

    def test_controlled_z_and_pt_groups(self):
        models = [model("GT001",.36,2,.09), model("GT002",.50,2,.12),
                  model("GT003",.50,5.2,.31)]
        rows = overview.join_catalogs(
            [extraction(row["setting_id"],2+i,.2+i/10) for i,row in enumerate(models)],
            models, Path("e"), Path("m"),
        )
        self.assertEqual(len(overview.controlled_groups(rows,"z")), 1)
        self.assertEqual(len(overview.controlled_groups(rows,"Pt")), 1)

    def test_schema_has_no_pt2_full_or_tier(self):
        names = " ".join(overview.FIELDS)
        self.assertNotIn("pT2", names)
        self.assertNotIn("Full", names)
        self.assertNotIn("Tier", names)
        self.assertIn("Pt_GeV", overview.FIELDS)
        self.assertIn("Q2", overview.FIELDS)
        self.assertIn("thpq", overview.FIELDS)
        self.assertNotIn("Q2_GeV2", overview.FIELDS)
        self.assertNotIn("theta_pq_deg", overview.FIELDS)

    def test_source_path_is_relative_to_rp_scripts(self):
        source = TOOLS.parent / "results/Pass1/example.csv"
        self.assertEqual(
            overview.relative_source(source),
            "rsidis_xs_v5/results/Pass1/example.csv",
        )

    def test_fixed_plot_ranges(self):
        source = Path(overview.__file__).read_text()
        self.assertIn('(0.0, 1.0) if xfield == "z" else (0.0, 0.6)', source)
        self.assertIn("SetNdivisions(510, True)", source)

    def test_each_pdf_page_uses_locked_fresh_canvas(self):
        canvas = overview.new_page_canvas("page_geometry_test")
        self.assertEqual(canvas.GetWw(), overview.CANVAS_WIDTH)
        self.assertEqual(canvas.GetWh(), overview.CANVAS_HEIGHT)
        canvas.Close()


if __name__ == "__main__":
    unittest.main()
