from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import RP_compare_simc_collimator_variants as comparison  # noqa: E402
import RP_extract_simc_info as extractor  # noqa: E402
import RP_generate_simc_collimator_study as generator  # noqa: E402


class CollimatorStudyTests(unittest.TestCase):
    def test_inserted_flags_are_unique_and_other_text_is_preserved(self):
        parent = """begin parm experiment
  ngen = 100000;
end parm experiment
begin parm simulate
  using_rad = 1;
end parm simulate
"""
        rendered = generator.render_study_input(parent)
        self.assertEqual(
            generator.assignment_values(rendered, "using_HMScoll"), ["0"]
        )
        self.assertEqual(
            generator.assignment_values(rendered, "using_SHMScoll"), ["1"]
        )
        changes = generator.changed_lines(parent, rendered)
        self.assertEqual(len(changes), 2)
        self.assertIn("using_rad = 1;", rendered)

    def test_existing_flags_are_replaced_without_duplicate_assignments(self):
        parent = """begin parm simulate
  using_HMScoll = 1;
  using_SHMScoll = 0;
end parm simulate
"""
        rendered = generator.render_study_input(parent)
        self.assertEqual(
            generator.assignment_values(rendered, "using_HMScoll"), ["0"]
        )
        self.assertEqual(
            generator.assignment_values(rendered, "using_SHMScoll"), ["1"]
        )

    def test_study_inventory_requires_one_row_per_reaction(self):
        rows = []
        for reaction, ngen in generator.EXPECTED_NGEN.items():
            rows.append({
                **generator.STUDY_IDENTITY,
                "Reaction": reaction,
                "Ngen": str(ngen),
                "Generation_status": "GENERATED",
            })
        selected = generator.select_study_rows(rows)
        self.assertEqual(
            [row["Reaction"] for row in selected], list(generator.REACTIONS)
        )
        with self.assertRaises(ValueError):
            generator.select_study_rows(rows[:-1])

    def test_extractor_uses_input_parent_as_output_category(self):
        with tempfile.TemporaryDirectory() as temporary:
            paths = extractor.expected_paths(
                Path(temporary),
                {
                    "Phase": "phase1",
                    "Run_type": "PIMINUS",
                    "Output_file": (
                        "infiles/RP_Simc/collimator_on/"
                        "mc_sidis_phase1_pass4_PIMINUS_LD2_"
                        "x0p25Q23p3z0p5thpq2p0_collimatorOn.inp"
                    ),
                },
            )
            self.assertEqual(paths["hist"].parent.name, "collimator_on")
            self.assertEqual(paths["raw_root"].parent.name, "collimator_on")

    @staticmethod
    def dmc_rows(mc_scale: float):
        rows = []
        for variable in comparison.VARIABLES:
            for index in range(1, 21):
                low = -0.20 + (index - 1) * 0.02
                high = low + 0.02
                if variable in {"hyp", "pyp"}:
                    low = -0.15 + (index - 1) * 0.015
                    high = low + 0.015
                rows.append({
                    "Phase": "phase1",
                    "Pass": "pass4",
                    "Run_type": "PIMINUS",
                    "Target": "LD2",
                    "Setting": "x0p25Q23p3z0p5thpq2p0",
                    "Tier": "delta",
                    "Variable": variable,
                    "Bin_index": str(index),
                    "Bin_low": str(low),
                    "Bin_high": str(high),
                    "Data_final": "10",
                    "Data_final_err": "2",
                    "MC_total": str(5 * mc_scale),
                    "MC_total_err": "1",
                    "MC_complete": "1",
                })
        return rows

    def test_ab_arithmetic_edge_flags_and_chi2(self):
        rows = comparison.build_comparison(
            self.dmc_rows(1.0), self.dmc_rows(2.0)
        )
        self.assertEqual(len(rows), 80)
        first = rows[0]
        self.assertAlmostEqual(float(first["Data_by_MC_baseline"]), 2.0)
        self.assertAlmostEqual(float(first["Data_by_MC_collimator_on"]), 1.0)
        self.assertAlmostEqual(
            float(first["MC_collimator_on_by_baseline"]), 2.0
        )
        self.assertGreater(sum(int(row["Is_edge"]) for row in rows), 0)
        self.assertGreater(float(first["Chi2_baseline_whole"]), 0.0)
        self.assertEqual(float(first["Chi2_collimator_on_whole"]), 0.0)
        self.assertGreater(int(first["NDF_baseline_edge"]), 0)

    def test_data_change_between_variants_is_rejected(self):
        baseline = self.dmc_rows(1.0)
        variant = self.dmc_rows(1.0)
        variant[0]["Data_final"] = "11"
        with self.assertRaises(ValueError):
            comparison.build_comparison(baseline, variant)


if __name__ == "__main__":
    unittest.main()
