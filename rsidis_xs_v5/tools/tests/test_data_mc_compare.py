#!/usr/bin/env python3

import csv
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path

import ROOT

TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import RP_data_mc_compare as compare  # noqa: E402


class DataMcComparisonTests(unittest.TestCase):
    def hist(self, name, value, error):
        result = compare.make_hist(name, "xbj")
        result.SetBinContent(1, value)
        result.SetBinError(1, error)
        return result

    def test_dummy_scales(self):
        self.assertAlmostEqual(compare.S_DUMMY[("phase1", "LH2")], 1 / 3.5274)
        self.assertAlmostEqual(compare.S_DUMMY[("phase1", "LD2")], 1 / 3.7825)
        self.assertAlmostEqual(compare.S_DUMMY[("phase2", "LH2")], 1 / 3.9031)
        self.assertAlmostEqual(compare.S_DUMMY[("phase2", "LD2")], 1 / 4.1854)

    def test_nminus1_removes_only_named_acceptance_cut(self):
        hms = compare.nminus1_cut("hxp", False)
        self.assertNotIn("H_gtr_th>", hms)
        self.assertIn("H_gtr_ph>", hms)
        self.assertIn("H_gtr_dp>-8", hms)
        hdelta = compare.nminus1_cut("hdelta", False)
        self.assertNotIn("H_gtr_dp>-8", hdelta)
        self.assertIn("P_gtr_dp>-10", hdelta)
        simc = compare.nminus1_cut("pyp", True)
        self.assertNotIn("ssyptar>", simc)
        self.assertIn("hsxptar>", simc)
        delta_tier = compare.nminus1_cut("hdelta", False, "delta")
        self.assertNotIn("H_gtr_dp>-8", delta_tier)
        self.assertIn("P_gtr_dp>-10", delta_tier)
        self.assertNotIn("H_gtr_th>", delta_tier)

    def test_all_variables_use_twenty_bins_and_requested_ranges(self):
        self.assertTrue(compare.BINNING)
        self.assertTrue(all(values[0] == 20 for values in compare.BINNING.values()))
        self.assertEqual(compare.BINNING["xbj"][:3], (20, 0.0, 0.6))
        self.assertEqual(compare.BINNING["Q2"][:3], (20, 0.0, 6.0))
        self.assertEqual(compare.BINNING["W"][:3], (20, 0.0, 4.0))
        self.assertEqual(compare.BINNING["thetapq"][:3], (20, 0.0, 0.2))
        self.assertEqual(compare.BINNING["pt2"][:3], (20, 0.0, 0.6))
        self.assertEqual(compare.BINNING["hxp"][:3], (20, -0.20, 0.20))
        self.assertEqual(compare.BINNING["hyp"][:3], (20, -0.15, 0.15))
        self.assertEqual(compare.BINNING["pxp"][:3], (20, -0.20, 0.20))
        self.assertEqual(compare.BINNING["pyp"][:3], (20, -0.15, 0.15))
        self.assertEqual((compare.CANVAS_WIDTH, compare.CANVAS_HEIGHT), (1800, 1200))

    def test_data_and_simc_use_matching_wide_angle_acceptance(self):
        expected = {
            "hxp": (-0.15, 0.15), "hyp": (-0.10, 0.10),
            "pxp": (-0.15, 0.15), "pyp": (-0.10, 0.10),
        }
        for variable, limits in expected.items():
            self.assertEqual(compare.DATA_ANGLE_CUTS[variable][1:], limits)
            self.assertEqual(compare.SIMC_ANGLE_CUTS[variable][1:], limits)

    def test_metadata_run_lists_wrap_without_losing_runs(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "signal.csv"
            expected = list(range(24000, 24040))
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["Run"])
                writer.writeheader()
                for run in expected:
                    writer.writerow({"Run": run})
            lines = compare.run_list_lines({("LH2", "Elec"): path}, "LH2")
            self.assertTrue(all(
                len(line) <= compare.METADATA_LINE_CHARS for line in lines
            ))
            observed = [
                int(value) for value in re.findall(
                    r"\b\d{5}\b", "\n".join(lines)
                )
            ]
            self.assertEqual(observed, expected)

    def test_neff_warning_filter_keeps_other_diagnostics(self):
        reason = "FULL_NEFF=42;FULL_REL_MC_ERR=0.1;DELTA_NEFF=11"
        self.assertEqual(
            compare.filtered_simc_warning(reason), "FULL_REL_MC_ERR=0.1"
        )

    def test_signal_stability_joins_leaf_current_by_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stem = "phase1_pass4_PIPLUS_LH2_Elec_x0p25Q23p3z0p5thpq2p0"
            norm = root / f"{stem}.csv"
            with norm.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=[
                    "Run", "Norm_status",
                    "RP_Normyield_delta", "RP_Normyield_err_delta",
                    "RP_Normyield_full", "RP_Normyield_err_full",
                ])
                writer.writeheader()
                writer.writerow({
                    "Run": 100, "Norm_status": "OK",
                    "RP_Normyield_delta": 2.5, "RP_Normyield_err_delta": 0.2,
                    "RP_Normyield_full": 2.0, "RP_Normyield_err_full": 0.15,
                })
                writer.writerow({
                    "Run": 101, "Norm_status": "OK",
                    "RP_Normyield_delta": 2.4, "RP_Normyield_err_delta": 0.2,
                    "RP_Normyield_full": 1.9, "RP_Normyield_err_full": 0.15,
                })
            bigtable = root / "bigtable"
            leaf = (
                bigtable / "phase1" / "pass4" / "PIPLUS" / "LH2"
                / "Elec" / "x0p25Q23p3z0p5thpq2p0" / f"{stem}.csv"
            )
            leaf.parent.mkdir(parents=True)
            with leaf.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["run", "BCM2_I"])
                writer.writeheader()
                writer.writerow({"run": 100, "BCM2_I": 28.5})
                writer.writerow({"run": 101, "BCM2_I": -999})
            points, reasons = compare.signal_stability_points(norm, bigtable)
            self.assertEqual(points["delta"], [(100.0, 28.5, 2.5, 0.2)])
            self.assertEqual(points["full"], [(100.0, 28.5, 2.0, 0.15)])
            self.assertEqual(reasons, ["BCM2_I_UNAVAILABLE_RUN_101"])

    def test_positron_and_dummy_subtraction_error_propagation(self):
        classes = {
            "target_e": {"xbj": self.hist("te", 10, 2)},
            "target_p": {"xbj": self.hist("tp", 3, 1)},
            "dummy_e": {"xbj": self.hist("de", 5, 1.5)},
            "dummy_p": {"xbj": self.hist("dp", 1, 0.5)},
        }
        result = compare.combine_data(classes, "xbj", 0.25)
        self.assertAlmostEqual(result["final"].GetBinContent(1), 6.0)
        expected_variance = 2**2 + 1**2 + 0.25**2 * (1.5**2 + 0.5**2)
        self.assertAlmostEqual(
            result["final"].GetBinError(1), math.sqrt(expected_variance)
        )

    def test_reaction_total_adds_values_and_variances(self):
        components = {
            reaction: {"xbj": self.hist(reaction, index + 1, 0.1 * (index + 1))}
            for index, reaction in enumerate(compare.REACTIONS)
        }
        total = compare.add_mc_total(components, "xbj")
        self.assertAlmostEqual(total.GetBinContent(1), 10.0)
        self.assertAlmostEqual(
            total.GetBinError(1),
            math.sqrt(sum((0.1 * index) ** 2 for index in range(1, 5))),
        )

    def test_derived_simc_expressions(self):
        native = compare.simc_expressions("sidis")
        derived = compare.simc_expressions("delta")
        self.assertEqual(native["xbj"], "xbj_recon")
        self.assertIn("Q2_recon", derived["xbj"])
        self.assertIn("phad/nu_recon", derived["z"])
        self.assertIn("sin(thetapq_recon)", derived["pt2"])
        self.assertEqual(derived["missmass"], "missmass_recon")

    def test_incomplete_mc_never_publishes_total_ratio_or_pull(self):
        data = {
            key: self.hist("data_" + key, value, 0.2)
            for key, value in (
                ("target_e", 10), ("target_p", 1), ("target_sub", 9),
                ("dummy_e", 2), ("dummy_p", 0.5), ("dummy_sub", 1.5),
                ("dummy_scaled", 0.4), ("final", 8.6),
            )
        }
        components = {
            reaction: {"xbj": self.hist("mc_" + reaction, 1, 0.1)}
            for reaction in compare.REACTIONS
        }
        total = compare.add_mc_total(components, "xbj")
        rows = compare.rows_for_variable(
            {
                "Phase": "phase1", "Pass": "pass4", "Run_type": "PIPLUS",
                "Target": "LH2", "Setting": "setting",
            },
            "delta", "xbj", data, components, total, 0.25,
            "PENDING", ["INCOMPLETE_MC=exclusive"], 0.0, False,
            {"sidis", "rho", "delta"},
        )
        first = rows[0]
        self.assertEqual(first["MC_complete"], 0)
        self.assertEqual(first["MC_total"], "")
        self.assertTrue(math.isnan(first["Data_by_MC"]))
        self.assertTrue(math.isnan(first["Pull"]))
        self.assertEqual(first["MC_exclusive"], "")
        self.assertEqual(first["Provisional_flag"], 1)


if __name__ == "__main__":
    unittest.main()
