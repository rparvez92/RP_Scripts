#!/usr/bin/env python3
"""Regression tests for explicit Delta-only and Full-cut data schemas."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

import RP_get_coin_normyield_Summary as norm_summary  # noqa: E402
import RP_get_good_coin_ev_Summary as good_summary  # noqa: E402


class DataAcceptanceSummaryTests(unittest.TestCase):
    def test_active_cpp_stages_use_matching_wide_angle_acceptance(self):
        project = TOOLS_DIR.parent
        goodcoin = (project / "macros" / "RP_get_good_coin_ev.C").read_text()
        xsec = (project / "macros" / "TableCoinXsec.C").read_text()
        for branch in ("H_gtr_th", "P_gtr_th"):
            self.assertIn(f"{branch}>-0.15", goodcoin)
            self.assertIn(f"{branch}<0.15", goodcoin)
            self.assertIn(f"({branch}>-0.15)", xsec)
            self.assertIn(f"({branch}<0.15)", xsec)
        for branch in ("H_gtr_ph", "P_gtr_ph"):
            self.assertIn(f"{branch}>-0.10", goodcoin)
            self.assertIn(f"{branch}<0.10", goodcoin)
            self.assertIn(f"({branch}>-0.10)", xsec)
            self.assertIn(f"({branch}<0.10)", xsec)
        for branch in ("hsxptar", "ssxptar"):
            self.assertIn(f"({branch}>-0.15)", xsec)
            self.assertIn(f"({branch}<0.15)", xsec)
        for branch in ("hsyptar", "ssyptar"):
            self.assertIn(f"({branch}>-0.10)", xsec)
            self.assertIn(f"({branch}<0.10)", xsec)

    def test_goodcoin_dual_tier_row_is_valid(self):
        row = {
            "CT_method": "TWO_STAGE_GAUSSIAN_FIT",
            "CT_POS_ref_status": "NOT_APPLICABLE",
            "CT_POS_ref_count": "0",
            "CTmean": "51.0",
            "CTmean_err": "0.1",
            "CTsigma": "0.4",
            "CTsigma_err": "0.05",
            "CTmean_residual": "0.2",
            "N_rndm_peak": "6",
            "RP_Goodcoin_delta": "1000",
            "RP_Goodcoin_err_delta": "35",
            "RP_Goodcoin_full": "600",
            "RP_Goodcoin_err_full": "27",
            "RP_Goodcoin_full_by_delta": "0.6",
            "ransubcoin_by_RP_Goodcoin_delta": "1.05",
            "Fit_status": "OK",
        }
        flagged = good_summary.flag_row(row)
        self.assertEqual(flagged["Problematic"], "FALSE")

    def test_missing_full_goodcoin_is_problematic(self):
        row = {
            "CT_method": "TWO_STAGE_GAUSSIAN_FIT",
            "CT_POS_ref_status": "NOT_APPLICABLE",
            "CT_POS_ref_count": "0",
            "CTmean": "51.0",
            "CTmean_err": "0.1",
            "CTsigma": "0.4",
            "CTsigma_err": "0.05",
            "CTmean_residual": "0.2",
            "N_rndm_peak": "6",
            "RP_Goodcoin_delta": "1000",
            "RP_Goodcoin_full": "nan",
            "ransubcoin_by_RP_Goodcoin_delta": "1.05",
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
