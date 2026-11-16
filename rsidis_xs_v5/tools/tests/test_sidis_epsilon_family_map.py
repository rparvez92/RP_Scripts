import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_build_sidis_epsilon_family_map as family  # noqa: E402


def group_rows(group_id, phase, pas, energy, run_type="pip", target="LH2",
               x="0.44", Q2="4.4", z="0.52"):
    common = {
        "Combination_group_id": group_id,
        "Combination_policy_version": "test-combination", "Phase": phase,
        "Pass": pas, "Run_type": run_type, "Target": target, "x": x,
        "Q2": Q2, "z": z, "Ebeam_GeV": str(energy), "N_contributors": "1",
        "Contributing_setting_ids": "RPTEST", "Contributing_thpq_deg": "2.0",
        "Binning_version": "test-binning",
    }
    return [{**common, "Pt_bin": str(index // 15 + 1),
             "phipq_bin": str(index % 15 + 1)} for index in range(75)]


class SidisEpsilonFamilyMapTests(unittest.TestCase):
    def test_pass4_pass5_family_is_lt_ready(self):
        rows = group_rows("g4", "phase1", "pass4", 8.5831)
        rows += group_rows("g5", "phase1", "pass5", 10.6716)
        members, summaries = family.build_family_map(rows)
        self.assertEqual(len(members), 2)
        self.assertEqual(len(summaries), 1)
        summary = summaries[0]
        self.assertEqual(summary["Family_status"], "LT_READY")
        self.assertEqual(summary["Pass_pattern"], "pass4|pass5")
        self.assertGreater(summary["Epsilon_span"], 0)

    def test_pass3_is_validation_member_for_anchored_family(self):
        rows = group_rows("g3", "phase2", "pass3", 6.4724)
        rows += group_rows("g4", "phase2", "pass4", 8.5814)
        rows += group_rows("g5", "phase2", "pass5", 10.6759)
        _, summaries = family.build_family_map(rows)
        self.assertEqual(summaries[0]["Family_status"],
                         "LT_READY_PASS3_VALIDATION")
        self.assertEqual(summaries[0]["N_distinct_epsilon"], 3)

    def test_phase_is_member_metadata_not_family_key(self):
        rows = group_rows("g4", "phase1", "pass4", 8.5831)
        rows += group_rows("g5", "phase2", "pass5", 10.6759)
        members, summaries = family.build_family_map(rows)
        self.assertEqual(len(summaries), 1)
        self.assertEqual([row["Member_phase"] for row in members],
                         ["phase1", "phase2"])

    def test_family_without_pass5_anchor_is_not_ready(self):
        rows = group_rows("g3", "phase2", "pass3", 6.4724)
        rows += group_rows("g4", "phase2", "pass4", 8.5814)
        _, summaries = family.build_family_map(rows)
        self.assertEqual(summaries[0]["Family_status"], "NOT_LT_READY")
        self.assertIn("MISSING_PASS5_ANCHOR", summaries[0]["Family_reason"])

    def test_duplicate_pass_is_rejected(self):
        rows = group_rows("g4a", "phase1", "pass4", 8.58)
        rows += group_rows("g4b", "phase2", "pass4", 8.60)
        rows += group_rows("g5", "phase2", "pass5", 10.67)
        _, summaries = family.build_family_map(rows)
        self.assertEqual(summaries[0]["Family_status"], "NOT_LT_READY")
        self.assertIn("DUPLICATE_PASS_MEMBER", summaries[0]["Family_reason"])

    def test_incomplete_combination_grid_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "expected 75"):
            family.build_family_map(group_rows("bad", "phase1", "pass5", 10.6)[:-1])


if __name__ == "__main__":
    unittest.main()
