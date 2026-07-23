#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import generate_simc_inputs as generator  # noqa: E402
import run_simc_batch as batch  # noqa: E402


class GeneratorTests(unittest.TestCase):
    def leaf(self, run_type: str = "PIMINUS", target: str = "LD2"):
        return generator.LeafResult(
            identity=generator.LeafIdentity("phase1", "pass4", run_type, target,
                                            "x0p25Q23p3z0p5thpq2p0"),
            source=Path("leaf.csv"),
            status="VALID",
            reason="test",
            ebeam_mev="8583.1",
            hms_p_mev="1531.0",
            hms_theta_deg="29.045",
            shms_p_mev="3632.0",
            shms_theta_deg="7.865",
        )

    def test_charge_and_reaction_flags(self):
        self.assertEqual(generator.reaction_flags("sidis", "PIMINUS")["doing_hplus"], 0)
        self.assertEqual(generator.reaction_flags("rho", "PIPLUS")["doing_hplus"], 1)
        self.assertEqual(generator.reaction_flags("exclusive", "PIMINUS")["which_pion"], 1)
        self.assertEqual(generator.reaction_flags("delta", "PIPLUS")["which_pion"], 2)
        self.assertEqual(generator.reaction_flags("delta", "PIMINUS")["which_pion"], 3)

    def test_rendered_input_validation(self):
        leaf = self.leaf()
        assignments = {
            "ngen": "1", "doing_pion": "0", "which_pion": "0", "doing_rho": "0",
            "doing_semi": "0", "doing_hplus": "1", "electron_arm": "1",
            "hadron_arm": "5", "Ebeam": "1", "spec%e%P": "1",
            "spec%e%theta": "1", "spec%p%P": "1", "spec%p%theta": "1",
            "Egamma_gen_max": "1", "targ%A": "2", "targ%Z": "1",
        }
        template = "\n".join(f"{name} = {value};" for name, value in assignments.items())
        rendered, _ = generator.render_input(template, leaf, "sidis")
        generator.validate_rendered_input(rendered, leaf, "sidis")
        broken = generator.replace_assignment(rendered, "electron_arm", "5")
        with self.assertRaisesRegex(ValueError, "electron_arm"):
            generator.validate_rendered_input(broken, leaf, "sidis")


class BatchTests(unittest.TestCase):
    def row(self):
        return {
            "Phase": "phase1", "Pass": "pass4", "Run_type": "PIMINUS",
            "Target": "LD2", "Setting": "x0p25Q23p3z0p5thpq2p0",
            "Reaction": "sidis", "Ngen": "100000",
            "Output_file": "infiles/RP_Simc/coin/example_phase1.inp",
            "Generation_status": "GENERATED",
        }

    def test_selection_and_safe_relative_input(self):
        row = self.row()
        selected = batch.select_rows([row], ["phase1"], ["sidis"], ["PIMINUS"], ["LD2"])
        self.assertEqual(len(selected), 1)
        self.assertEqual(batch.input_relative_to_infiles(row),
                         Path("RP_Simc/coin/example_phase1.inp"))
        row["Output_file"] = "../outside.inp"
        with self.assertRaisesRegex(ValueError, "unsafe"):
            batch.input_relative_to_infiles(row)

    def test_hist_ngen_and_output_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = {
                "hist": root / "test.hist", "runout": root / "test.out",
                "raw_root": root / "test.root", "recon_root": root / "recon.root",
            }
            self.assertEqual(batch.inspect_outputs(paths, 100)[0], "PENDING")
            paths["hist"].write_text("Ngen (request) = 100\n", encoding="utf-8")
            self.assertEqual(batch.read_hist_ngen(paths["hist"]), 100)
            status, reason, observed = batch.inspect_outputs(paths, 100)
            self.assertEqual(status, "PARTIAL_OR_MISMATCH")
            self.assertEqual(observed, 100)
            self.assertIn("missing=", reason)

    def test_manifest_schema(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "manifest.csv"
            with manifest.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=self.row().keys())
                writer.writeheader()
                writer.writerow(self.row())
            self.assertEqual(len(batch.read_manifest(manifest)), 1)


if __name__ == "__main__":
    unittest.main()
