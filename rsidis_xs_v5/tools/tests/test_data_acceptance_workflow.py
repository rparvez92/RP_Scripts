#!/usr/bin/env python3
"""Regression tests for the data-acceptance result schemas."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import RP_get_coin_normyield_Summary as norm_summary  # noqa: E402
import RP_get_good_coin_ev_Summary as good_summary  # noqa: E402


class DataAcceptanceSummaryTests(unittest.TestCase):
    def test_goodcoin_uses_only_baseline_delta_acceptance(self):
        project = TOOLS_DIR.parent
        goodcoin = (project / "macros" / "RP_get_good_coin_ev.C").read_text()
        xsec = (project / "macros" / "TableCoinXsec.C").read_text()
        for branch in ("H_gtr_th", "P_gtr_th"):
            self.assertNotIn(f"{branch}>-0.15", goodcoin)
            self.assertNotIn(f"{branch}<0.15", goodcoin)
            self.assertIn(f"({branch}>-0.15)", xsec)
            self.assertIn(f"({branch}<0.15)", xsec)
        for branch in ("H_gtr_ph", "P_gtr_ph"):
            self.assertNotIn(f"{branch}>-0.10", goodcoin)
            self.assertNotIn(f"{branch}<0.10", goodcoin)
            self.assertIn(f"({branch}>-0.10)", xsec)
            self.assertIn(f"({branch}<0.10)", xsec)
        for branch in ("hsxptar", "ssxptar"):
            self.assertIn(f"({branch}>-0.15)", xsec)
            self.assertIn(f"({branch}<0.15)", xsec)
        for branch in ("hsyptar", "ssyptar"):
            self.assertIn(f"({branch}>-0.10)", xsec)
            self.assertIn(f"({branch}<0.10)", xsec)

    def test_goodcoin_pass1_configuration(self):
        project = TOOLS_DIR.parent
        goodcoin = (project / "macros" / "RP_get_good_coin_ev.C").read_text()
        self.assertIn("{-6, -5, -4, -3, -2, 2, 3, 4, 5, 6}", goodcoin)
        self.assertIn('row.rfPeriod = input.rfSignal;', goodcoin)
        self.assertIn('results/Pass1/PDFs/RP_get_good_coin_ev/', goodcoin)
        self.assertIn('results/Pass1/Tables/RP_get_good_coin_ev/', goodcoin)
        self.assertIn('/Volumes/T7/RSIDIS/Pass1/Data/ROOTfiles/SkimmedData', goodcoin)
        self.assertIn('"RP_Goodcoin,RP_Goodcoin_err,CoinCuts_entries,ransubcoin,"', goodcoin)
        self.assertNotIn("RP_Goodcoin_full", goodcoin)

    def test_goodcoin_baseline_row_is_valid(self):
        row = {
            "CT_method": "TWO_STAGE_GAUSSIAN_FIT",
            "CT_POS_ref_status": "NOT_APPLICABLE",
            "CT_POS_ref_count": "0",
            "CTmean": "51.0",
            "CTmean_err": "0.1",
            "CTsigma": "0.4",
            "CTsigma_err": "0.05",
            "CTmean_residual": "0.2",
            "CTsigma_residual": "0.1",
            "N_rndm_peak": "10",
            "RP_Goodcoin": "1000",
            "RP_Goodcoin_err": "35",
            "ransubcoin_by_RP_Goodcoin": "1.05",
            "Fit_status": "OK",
        }
        flagged = good_summary.flag_row(row)
        self.assertEqual(flagged["Problematic"], "FALSE")

    def test_missing_goodcoin_is_problematic(self):
        row = {
            "CT_method": "TWO_STAGE_GAUSSIAN_FIT",
            "CT_POS_ref_status": "NOT_APPLICABLE",
            "CT_POS_ref_count": "0",
            "CTmean": "51.0",
            "CTmean_err": "0.1",
            "CTsigma": "0.4",
            "CTsigma_err": "0.05",
            "CTmean_residual": "0.2",
            "CTsigma_residual": "0.1",
            "N_rndm_peak": "10",
            "RP_Goodcoin": "nan",
            "ransubcoin_by_RP_Goodcoin": "1.05",
            "Fit_status": "OK",
        }
        flagged = good_summary.flag_row(row)
        self.assertEqual(flagged["Yield_flag"], "TRUE")

    def test_normyield_dual_tier_row_is_valid(self):
        row = {
            "Norm_status": "OK",
            "Norm_factor": "0.01",
            "RP_Normyield_delta": "10.0",
            "RP_Normyield_err_delta": "0.4",
            "RP_Normyield_full": "6.0",
            "RP_Normyield_err_full": "0.3",
            "RP_Normyield_full_by_delta": "0.6",
            "normyield": "10.5",
            "normyield_by_RP_Normyield_delta": "1.05",
        }
        flagged = norm_summary.flag_row(row)
        self.assertEqual(flagged["Problematic"], "FALSE")


if __name__ == "__main__":
    unittest.main()
