import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_build_thpq_combination as combination  # noqa: E402


class ThpqCombinationTests(unittest.TestCase):
    def model(self, setting, theta, target="LH2", z="0.36"):
        return {
            "Phase":"phase1", "Pass":"pass5", "Run_type":"pim",
            "Target":target, "setting_id":setting, "Ebeam_GeV":"10.6716",
            "x":"0.25", "Q2_GeV2":"3.3", "z":z,
            "theta_pq_deg":str(theta), "Model_status":"OK",
        }

    def test_requested_human_readable_group_id(self):
        row={"Phase":"phase1","Pass":"pass5","Run_type":"pim","Target":"LH2",
             "x":.25,"Q2":3.3,"z":.36}
        self.assertEqual(combination.combination_group_id(row),
                         "1501x0p25Q23p3z0p36")

    def test_pip_ld2_mapping_and_fixed_numeric_format(self):
        row={"Phase":"phase2","Pass":"pass5","Run_type":"pip","Target":"LD2",
             "x":.31,"Q2":3.1,"z":.5}
        self.assertEqual(combination.combination_group_id(row),
                         "2512x0p31Q23p1z0p50")

    def test_theta_settings_combine_but_z_and_target_do_not(self):
        rows=[self.model("A",2),self.model("B",5.8),self.model("C",8),
              self.model("D",2,z="0.50"),self.model("E",2,target="LD2")]
        mapping,groups=combination.build_rows(rows)
        self.assertEqual(len(groups),3)
        combined=next(row for row in groups if row["z"]=="0.36" and row["Target"]=="LH2")
        self.assertEqual(combined["N_contributors"],3)
        self.assertEqual(combined["Contributing_setting_ids"],"A|B|C")
        self.assertEqual(combined["Contributing_thpq_deg"],"2|5.8|8")
        self.assertTrue(all(row["Combination_status"]=="OK" for row in mapping))


if __name__ == "__main__": unittest.main()
