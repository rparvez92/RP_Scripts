#!/usr/bin/env python3

import math
import sys
import tempfile
import unittest
from array import array
from pathlib import Path

import ROOT


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import RP_extract_simc_info as extractor  # noqa: E402


class ExtractorTests(unittest.TestCase):
    def manifest_row(self, status="GENERATED"):
        return {
            "Phase": "phase1",
            "Pass": "pass4",
            "Run_type": "PIPLUS",
            "Target": "LH2",
            "Setting": "x0p25Q23p3z0p5thpq2p0",
            "Reaction": "sidis",
            "Ngen": "1001",
            "Output_file": (
                "infiles/RP_Simc/coin/"
                "mc_sidis_phase1_pass4_PIPLUS_LH2_x0p25Q23p3z0p5thpq2p0.inp"
            ),
            "Generation_status": status,
            "Diagnostic_reason": "test status",
            "Source_leaf_csv": "leaf.csv",
            "Source_template": "template.inp",
        }

    def output_paths(self, root):
        stem = Path(self.manifest_row()["Output_file"]).stem
        simulation = root / "Phase1" / "Simulation"
        hist_dir = simulation / "outfiles" / "coin"
        root_dir = simulation / "ROOTfiles" / "coin"
        hist_dir.mkdir(parents=True)
        root_dir.mkdir(parents=True)
        return {
            "hist": hist_dir / f"{stem}.hist",
            "raw": root_dir / f"{stem}.root",
            "recon": root_dir / f"recon_hcana_{stem}.root",
        }

    def write_tree_pair(
        self, paths, entries=1001, hist_ngen=None, recon_ngen=None,
        normfac=100.0, fweight_scale=1.0, negative=False,
        nonconstant_ngen=False, nonconstant_normfac=False,
        sentinel_variables=(), omit_kinematics=(), outside_selection=False,
    ):
        hist_ngen = entries if hist_ngen is None else hist_ngen
        recon_ngen = entries if recon_ngen is None else recon_ngen
        paths["hist"].write_text(
            f"Ngen (request) = {hist_ngen}\n"
            f"normfac = {normfac:.12e}\n",
            encoding="utf-8",
        )

        required_floats = (
            "Weight", "hsdelta", "hsxptar", "hsyptar",
            "ssdelta", "ssxptar", "ssyptar",
        )
        kinematic_pairs = tuple(
            pair for pair in (
            ("xbj", "xbj_recon"), ("Q2", "Q2_recon"), ("W", "W_recon"),
            ("z", "z_recon"), ("thetapq", "thetapq_recon"),
            ("phipq", "phipq_recon"), ("pt2", "pt2_recon"),
            )
            if pair[0] not in omit_kinematics
        )

        raw_file = ROOT.TFile(paths["raw"].as_posix(), "RECREATE")
        raw_tree = ROOT.TTree("h10", "raw")
        raw_values = {name: array("f", [0.0]) for name in required_floats}
        for name, value in raw_values.items():
            raw_tree.Branch(name, value, f"{name}/F")
        for original, _ in kinematic_pairs:
            value = array("f", [0.0])
            raw_values[original] = value
            raw_tree.Branch(original, value, f"{original}/F")

        recon_file = ROOT.TFile(paths["recon"].as_posix(), "RECREATE")
        recon_tree = ROOT.TTree("h10", "recon")
        recon_values = {name: array("f", [0.0]) for name in required_floats}
        for name, value in recon_values.items():
            recon_tree.Branch(name, value, f"{name}/F")
        ngen_value = array("i", [recon_ngen])
        normfac_value = array("d", [normfac])
        fweight_value = array("f", [0.0])
        recon_tree.Branch("Ngen", ngen_value, "Ngen/I")
        recon_tree.Branch("normfac", normfac_value, "normfac/D")
        recon_tree.Branch("fWeight", fweight_value, "fWeight/F")
        for original, reconstructed in kinematic_pairs:
            original_value = array("f", [0.0])
            recon_value = array("f", [0.0])
            recon_values[original] = original_value
            recon_values[reconstructed] = recon_value
            recon_tree.Branch(original, original_value, f"{original}/F")
            recon_tree.Branch(reconstructed, recon_value, f"{reconstructed}/F")

        for index in range(entries):
            ngen_value[0] = (
                recon_ngen + 1
                if nonconstant_ngen and index == entries - 1
                else recon_ngen
            )
            normfac_value[0] = (
                normfac * 2.0
                if nonconstant_normfac and index == entries - 1
                else normfac
            )
            weight = -1.0 if negative and index == 0 else 1.0
            for values in (raw_values, recon_values):
                values["Weight"][0] = weight
                values["hsdelta"][0] = 9.0 if outside_selection else 0.0
                values["hsxptar"][0] = 0.0
                values["hsyptar"][0] = 0.0
                values["ssdelta"][0] = 0.0
                values["ssxptar"][0] = 0.0
                values["ssyptar"][0] = 0.0
            for offset, (original, reconstructed) in enumerate(kinematic_pairs):
                value = (
                    -999.0
                    if original in sentinel_variables
                    else 0.1 * (offset + 1) + 1e-5 * index
                )
                raw_values[original][0] = value
                recon_values[original][0] = value
                recon_values[reconstructed][0] = (
                    -999.0 if original in sentinel_variables else value + 0.01
                )
            expected = weight * normfac_value[0] / ngen_value[0]
            fweight_value[0] = expected * fweight_scale
            raw_tree.Fill()
            recon_tree.Fill()

        raw_file.cd()
        raw_tree.Write()
        raw_file.Close()
        recon_file.cd()
        recon_tree.Write()
        recon_file.Close()

    def test_hist_parser_and_exact_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            paths["hist"].write_text(
                "Ngen (request) = 123\nnormfac = 0.452308D+11\n",
                encoding="utf-8",
            )
            self.assertEqual(extractor.parse_hist(paths["hist"]), (123, 4.52308e10))
            resolved = extractor.expected_paths(root, self.manifest_row())
            self.assertEqual(resolved["hist"], paths["hist"])
            self.assertFalse(resolved["hist"].name.startswith("._"))

    def test_skipped_and_pending_are_not_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            skipped = extractor.inspect_manifest_row(
                self.manifest_row("SKIPPED_INVALID_KINEMATICS"), root
            )
            pending = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(skipped["QA_status"], "SKIPPED")
            self.assertEqual(pending["QA_status"], "PENDING")

    def test_healthy_normalization_and_weight_sums(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "OK", result["QA_reason"])
            self.assertEqual(result["Raw_entries"], 1001)
            self.assertEqual(result["Recon_entries"], 1001)
            self.assertEqual(result["Recon_Ngen"], 1001)
            self.assertEqual(result["Recon_normfac"], 100.0)
            self.assertEqual(result["N_selected_nocut"], 1001)
            self.assertEqual(result["N_selected_full"], 1001)
            self.assertAlmostEqual(result["SimYield_full"], 100.0, places=5)
            self.assertAlmostEqual(result["SimYield_err_full"], 100 / math.sqrt(1001), places=5)
            self.assertAlmostEqual(result["Neff_full"], 1001.0, places=3)
            self.assertLess(result["FWeight_identity_max_rel_diff"], 1e-5)
            for tier in ("delta", "full"):
                self.assertAlmostEqual(
                    result[f"Q2_{tier}_simc_weighted_mean"], 0.205, places=5
                )
                self.assertAlmostEqual(
                    result[f"Q2_{tier}_residual_weighted_mean"],
                    0.01,
                    places=5,
                )
            self.assertNotIn("Q2_residual_weighted_mean", result)
            self.assertAlmostEqual(
                result["Raw_root_size_MB"],
                paths["raw"].stat().st_size / 1e6,
            )
            self.assertNotIn("Raw_root_size_bytes", result)
            self.assertNotIn("N_selected_all", result)
            self.assertNotIn("Weight_sum_all", result)
            self.assertIn("Weight_sum_nocut", result)

    def test_grouped_setting_positions(self):
        rows = [
            {"Setting": "setting_b", "Run_type": "PIPLUS", "Target": "LH2"},
            {"Setting": "setting_a", "Run_type": "PIMINUS", "Target": "LH2"},
            {"Setting": "setting_a", "Run_type": "PIPLUS", "Target": "LH2"},
            {"Setting": "setting_a", "Run_type": "PIMINUS", "Target": "LD2"},
        ]
        settings = extractor.setting_names(rows)
        self.assertEqual(settings, ["setting_a", "setting_b"])
        setting_index = {
            setting: index for index, setting in enumerate(settings, start=1)
        }
        setting_a = [row for row in rows if row["Setting"] == "setting_a"]
        category_positions = [
            extractor.grouped_x(row, setting_index) for row in setting_a
        ]
        self.assertTrue(all(0.5 < value < 1.5 for value in category_positions))
        self.assertEqual(len(set(category_positions)), len(category_positions))
        raw_x = extractor.grouped_x(setting_a[0], setting_index, 0, 3)
        hist_x = extractor.grouped_x(setting_a[0], setting_index, 1, 3)
        recon_x = extractor.grouped_x(setting_a[0], setting_index, 2, 3)
        self.assertAlmostEqual(hist_x - raw_x, 0.060)
        self.assertAlmostEqual(recon_x - hist_x, 0.060)

    def test_setting_chunks_never_split_categories(self):
        for count, expected_chunks in ((1, 1), (6, 1), (7, 2), (10, 2)):
            rows = [
                {
                    "Setting": f"setting_{setting:02d}",
                    "Run_type": run_type,
                    "Target": target,
                }
                for setting in range(count)
                for run_type, target in extractor.CATEGORY_ORDER
            ]
            chunks = extractor.chunk_rows_by_settings(rows)
            self.assertEqual(len(chunks), expected_chunks)
            seen = []
            for chunk in chunks:
                names = extractor.setting_names(chunk)
                self.assertLessEqual(len(names), 6)
                for name in names:
                    self.assertEqual(
                        sum(row["Setting"] == name for row in chunk), 4
                    )
                seen.extend(names)
            self.assertEqual(seen, extractor.setting_names(rows))

    def test_kinematic_availability_states(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(
                paths,
                sentinel_variables=("xbj",),
                omit_kinematics=("z",),
            )
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(
                extractor.kinematic_availability(result, "xbj"),
                "sentinel_only",
            )
            self.assertEqual(
                extractor.kinematic_availability(result, "z"), "missing"
            )
            self.assertEqual(
                extractor.kinematic_availability(result, "Q2"), "available"
            )
            self.assertEqual(result["xbj_full_simc_weighted_mean"], "")
            self.assertEqual(
                extractor.unavailable_kinematic_message(
                    [result], "xbj", "full"
                ),
                "Unavailable: all values = -999",
            )
            self.assertEqual(
                extractor.unavailable_kinematic_message(
                    [result], "z", "full"
                ),
                "Unavailable: branch missing",
            )

    def test_no_valid_kinematics_after_selection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, outside_selection=True)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(
                extractor.kinematic_availability(result, "Q2"), "available"
            )
            self.assertEqual(result["Q2_delta_simc_weighted_mean"], "")
            self.assertEqual(
                extractor.unavailable_kinematic_message(
                    [result], "Q2", "delta"
                ),
                "No valid values after selection",
            )

    def test_identity_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, fweight_scale=1.2)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "ERROR")
            self.assertIn("FWEIGHT_IDENTITY_MISMATCH", result["QA_reason"])

    def test_negative_weight_is_warning(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, negative=True)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "WARNING")
            self.assertIn("NEGATIVE_FWEIGHT", result["QA_reason"])

    def test_entry_ngen_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, hist_ngen=1000)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "ERROR")
            self.assertIn("ENTRY_NGEN_MISMATCH", result["QA_reason"])

    def test_nonconstant_recon_ngen_is_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, nonconstant_ngen=True)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "ERROR")
            self.assertEqual(result["Recon_Ngen"], "")
            self.assertIn(
                "RECON_NGEN_NONCONSTANT:min=1001,max=1002",
                result["QA_reason"],
            )

    def test_nonconstant_recon_normfac_is_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.output_paths(root)
            self.write_tree_pair(paths, nonconstant_normfac=True)
            result = extractor.inspect_manifest_row(self.manifest_row(), root)
            self.assertEqual(result["QA_status"], "ERROR")
            self.assertEqual(result["Recon_normfac"], "")
            self.assertIn(
                "RECON_NORMFAC_NONCONSTANT:min=100,max=200",
                result["QA_reason"],
            )


if __name__ == "__main__":
    unittest.main()
