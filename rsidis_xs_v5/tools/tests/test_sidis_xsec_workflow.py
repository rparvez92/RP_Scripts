import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_build_sidis_model as model  # noqa: E402
import RP_extract_sidis_xsec as extract  # noqa: E402


class SidisXsecWorkflowTests(unittest.TestCase):
    def test_current_identity_schema(self):
        self.assertEqual(
            extract.IDENTITY,
            ("Phase", "Pass", "Run_type", "Target", "setting_id"),
        )
        self.assertNotIn("Tier", extract.FIELDS)
        self.assertFalse(any("full" in field.lower() for field in extract.FIELDS))

    def test_model_pt_is_a_magnitude(self):
        positive = model.derive_kinematics(10.6716, 0.25, 3.3, 0.5, 2.0)
        negative = model.derive_kinematics(10.6716, 0.25, 3.3, 0.5, -2.0)
        self.assertAlmostEqual(positive["pt"], negative["pt"])
        self.assertGreater(positive["pt"], 0)

    def test_integrates_untiered_baseline_rows(self):
        rows = []
        for index in range(2):
            row = {
                "Variable": "hdelta", "Bin_index": str(index + 1),
                "Data_final": str(2 + index), "Data_final_err": "1",
                "MC_complete": "1", "Data_closure_rel": "0.001",
            }
            for reaction, value in zip(extract.REACTIONS, (1, 2, 3, 4)):
                row[f"MC_{reaction}"] = str(value + index)
                row[f"MC_{reaction}_err"] = "0.5"
            rows.append(row)
        result = extract.integrate_baseline_rows(rows)
        self.assertEqual(result["Data_final"], 5.0)
        self.assertAlmostEqual(result["Data_final_err"], math.sqrt(2))
        self.assertEqual(
            sum(result[f"MC_{reaction}"] for reaction in extract.REACTIONS),
            24.0,
        )

    def test_ratio_error_uses_data_and_mc_uncertainty(self):
        value, error = extract.ratio_with_error(12.0, 3.0, 4.0, 0.5)
        self.assertEqual(value, 3.0)
        self.assertAlmostEqual(error, math.hypot(3.0 / 4.0, 12.0 * 0.5 / 16.0))


if __name__ == "__main__":
    unittest.main()
