import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_sidis_binning as binning  # noqa: E402


class SidisBinningTests(unittest.TestCase):
    def test_shared_analysis_binning_is_fifteen_by_five(self):
        self.assertEqual(len(binning.AZIMUTHAL_PHI_EDGES) - 1, 15)
        self.assertEqual(len(binning.AZIMUTHAL_PT_EDGES) - 1, 5)
        self.assertAlmostEqual(binning.AZIMUTHAL_PHI_EDGES[-1], 2 * math.pi)
        self.assertEqual(binning.AZIMUTHAL_PT_EDGES, (0, .1, .2, .3, .4, .5))

    def test_phi_is_half_open_and_final_pt_edge_can_close(self):
        self.assertEqual(binning.canonical_phi(2 * math.pi), 0.0)
        self.assertAlmostEqual(binning.canonical_phi(-math.pi / 2), 3 * math.pi / 2)
        self.assertEqual(binning.bin_index(0, binning.AZIMUTHAL_PHI_EDGES), 1)
        self.assertIsNone(binning.bin_index(2 * math.pi, binning.AZIMUTHAL_PHI_EDGES))
        self.assertEqual(binning.bin_index(.5, binning.AZIMUTHAL_PT_EDGES, True), 5)


if __name__ == "__main__":
    unittest.main()
