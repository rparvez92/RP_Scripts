import csv
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


model = load("RP_build_sidis_model")
extract = load("RP_extract_sidis_xsec")


class ModelTests(unittest.TestCase):
    def test_setting_parser(self):
        parsed = model.parse_setting("x0p25Q23p3z0p67thpqneg0p8")
        self.assertEqual(parsed, {"x": 0.25, "Q2": 3.3, "z": 0.67, "thpq": -0.8})

    def test_pt_is_derived_magnitude_not_theta(self):
        result = model.derive_kinematics(10.6716, 0.25, 3.3, 0.5, 2.0)
        expected = result["ppion"] * math.sin(math.radians(2.0))
        self.assertAlmostEqual(result["pt"], expected)
        self.assertNotAlmostEqual(result["pt"], 2.0)
        negative = model.derive_kinematics(10.6716, 0.25, 3.3, 0.5, -2.0)
        self.assertAlmostEqual(negative["pt"], expected)

    def test_sentinel_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "sentinel"):
            model.derive_kinematics(10.6716, -999, 3.3, 0.5, 2.0)

    def test_calculator_parser(self):
        phi, sighad, sigma = model.parse_calculator_output(
            "header\n 3.141592653589793 0.2955564 7.63286E-003\n"
        )
        self.assertAlmostEqual(phi, math.pi)
        self.assertAlmostEqual(sighad, 0.2955564)
        self.assertAlmostEqual(sigma, 0.00763286)

    def test_target_and_charge_mapping(self):
        self.assertEqual(model.calculator_identity("LH2", "PIPLUS"), (1.0, 1.0, 1.0))
        self.assertEqual(model.calculator_identity("LD2", "PIMINUS"), (2.0, 1.0, -1.0))


class ExtractionTests(unittest.TestCase):
    def test_ratio_and_error(self):
        value, error = extract.ratio_with_error(12.0, 3.0, 4.0, 0.5)
        self.assertAlmostEqual(value, 3.0)
        self.assertAlmostEqual(error, math.hypot(3.0 / 4.0, 12.0 * 0.5 / 16.0))

    def test_integrates_delta_only_and_total_is_component_sum(self):
        rows = []
        for index in range(2):
            row = {"Tier": "delta", "Variable": "hdelta", "Bin_index": str(index + 1),
                   "Data_final": str(2 + index), "Data_final_err": "1",
                   "MC_complete": "1", "Data_closure_rel": "0"}
            for reaction, value in zip(extract.REACTIONS, (1, 2, 3, 4)):
                row[f"MC_{reaction}"] = str(value + index)
                row[f"MC_{reaction}_err"] = "0.5"
            row["MC_total"] = "10"
            row["MC_total_err"] = "1"
            rows.append(row)
        rows.append({"Tier": "full", "Variable": "hdelta", "Bin_index": "1",
                     "Data_final": "999"})
        result = extract.integrate_delta_rows(rows)
        self.assertEqual(result["Data_final"], 5.0)
        self.assertAlmostEqual(result["Data_final_err"], math.sqrt(2))
        components = sum(result[f"MC_{reaction}"] for reaction in extract.REACTIONS)
        self.assertEqual(components, 24.0)

    def test_total_mc_not_sidis_is_denominator(self):
        value, _ = extract.ratio_with_error(20.0, 1.0, 10.0 + 2.0 + 3.0 + 5.0, 1.0)
        self.assertEqual(value, 1.0)
        self.assertNotEqual(value, 20.0 / 10.0)

    def test_zero_denominator(self):
        value, error = extract.ratio_with_error(1.0, 1.0, 0.0, 1.0)
        self.assertTrue(math.isnan(value))
        self.assertTrue(math.isnan(error))

    def test_public_fields_reserve_R(self):
        self.assertIn("C_Y", extract.FIELDS)
        self.assertNotIn("R", extract.FIELDS)
        self.assertFalse(any("full" in field.lower() for field in extract.FIELDS))

    def test_model_skipped_precedes_missing_delta_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "phase2_pass3_PIPLUS_LH2_xneg999Q2neg999zneg999thpqneg999.csv"
            with source.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=[*extract.IDENTITY, "Status"])
                writer.writeheader()
                writer.writerow(dict(zip(extract.IDENTITY,
                    ("phase2", "pass3", "PIPLUS", "LH2", "xneg999Q2neg999zneg999thpqneg999")),
                    Status="PENDING"))
            identity = ("phase2", "pass3", "PIPLUS", "LH2", "xneg999Q2neg999zneg999thpqneg999")
            models = {identity: {"Model_status": "SKIPPED", "Model_reason": "sentinel"}}
            row = extract.extract_one(source, models, {}, Path("model.csv"))
            self.assertEqual(row["Extraction_status"], "SKIPPED")


if __name__ == "__main__":
    unittest.main()
