import importlib
import json
import math
import sys
import tempfile
import unittest
from array import array
import csv
from pathlib import Path

import ROOT

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
common = importlib.import_module("RP_sidis_differential_common")
extract = importlib.import_module("RP_extract_sidis_1d")


def hist(name, values, errors=None):
    result = ROOT.TH1D(name, "", len(values), 0.0, float(len(values)))
    for i, value in enumerate(values, 1):
        result.SetBinContent(i, value)
        result.SetBinError(i, (errors or [0.0] * len(values))[i - 1])
    return result


def response(name, matrix):
    result = ROOT.TH2D(name, "", len(matrix), 0, len(matrix), len(matrix), 0, len(matrix))
    for ix, row in enumerate(matrix, 1):
        for iy, value in enumerate(row, 1):
            result.SetBinContent(ix, iy, value)
    return result


class DifferentialWorkflowTests(unittest.TestCase):
    def test_response_purity_stability_and_abs_weight_semantics(self):
        matrix = response("r_metrics", [[8, 2], [1, 9]])
        purity, stability = common.response_metrics(matrix, 1)
        self.assertTrue(math.isclose(purity, 8 / 9))
        self.assertTrue(math.isclose(stability, 8 / 10))

    def test_adjacent_merge_and_failed_bin_is_reported(self):
        data = hist("d_merge", [10, 10], [8, 8])
        mc = hist("m_merge", [10, 10], [2, 2])
        matrix = response("r_merge", [[9, 1], [1, 9]])
        edges, details = common.recommend_adjacent_edges(data, mc, matrix)
        self.assertEqual(edges, [0.0, 2.0])
        self.assertEqual(len(details), 1)
        self.assertFalse(details[0]["passes"])

    def test_periodic_phi_reduces_globally_to_allowed_equal_bins(self):
        data = hist("d_phi", [100] * 24, [1] * 24)
        mc = hist("m_phi", [100] * 24, [1] * 24)
        matrix = response("r_phi", [[100 if i == j else 0 for j in range(24)] for i in range(24)])
        edges, details = common.recommend_phi_edges(data, mc, matrix)
        self.assertEqual(len(details), 24)
        self.assertTrue(math.isclose(edges[0], 0.0) and math.isclose(edges[-1], 24.0))
        self.assertEqual(len({round(edges[i + 1] - edges[i], 12) for i in range(24)}), 1)

    def test_ratio_and_extraction_arithmetic(self):
        ratio, error = extract.ratio_error(12.0, 3.0, 4.0, 1.0)
        self.assertEqual(ratio, 3.0)
        self.assertTrue(math.isclose(error, math.hypot(3 / 4, 12 / 16)))
        self.assertTrue(math.isclose(ratio * 2.5, 7.5))

    def test_json_writer_replaces_nonfinite_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            common.atomic_json(path, {"bad": math.nan, "good": 1.0})
            self.assertEqual(json.loads(path.read_text()), {"bad": None, "good": 1.0})

    def test_schema_has_no_full_cut_or_neff(self):
        joined = " ".join(extract.FIELDS)
        self.assertNotIn("Full", joined); self.assertNotIn("full", joined)
        self.assertNotIn("Neff", joined); self.assertNotIn("NEFF", joined)
        self.assertIn("C_Y", extract.FIELDS); self.assertIn("Y_MC", extract.FIELDS)

    def test_toy_root_response_matrix_and_periodic_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "toy.root"
            output = ROOT.TFile(str(path), "RECREATE"); tree = ROOT.TTree("h10", "toy")
            branches = {name: array("f", [0.0]) for name in ("hsdelta", "ssdelta", "z", "z_recon", "phipq", "phipq_recon", "fWeight")}
            for name, value in branches.items(): tree.Branch(name, value, f"{name}/F")
            events = ((0.10, 0.12, -0.01, 2*math.pi-0.01, 1.0),
                      (0.40, 0.42, 2*math.pi-0.01, 0.01, -2.0))
            for original, reconstructed, phi, phi_recon, weight in events:
                branches["hsdelta"][0]=0; branches["ssdelta"][0]=0
                branches["z"][0]=original; branches["z_recon"][0]=reconstructed
                branches["phipq"][0]=phi; branches["phipq_recon"][0]=phi_recon; branches["fWeight"][0]=weight
                tree.Fill()
            tree.Write(); output.Close()
            row={"Recon_root_path":str(path)}
            z_response=common.response_and_moments(row,"z",[0,.25,.5],"toy")
            self.assertAlmostEqual(z_response.Integral(),3.0)  # abs(1)+abs(-2)
            phi_response=common.response_and_moments(row,"phipq",[0,math.pi,2*math.pi],"toyphi")
            self.assertAlmostEqual(phi_response.Integral(),3.0)

    def test_upstream_integrated_closure_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"dmc.csv"
            fields=["Tier","Variable","Bin_index","Data_final",*(f"MC_{r}" for r in common.REACTIONS)]
            with path.open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
                for index in (1,2):
                    writer.writerow({"Tier":"delta","Variable":"hdelta","Bin_index":index,"Data_final":index,**{f"MC_{r}":2*index for r in common.REACTIONS}})
            values=extract.upstream_integrals(path)
            self.assertEqual(values["Data"],3.0)
            self.assertTrue(all(values[r]==6.0 for r in common.REACTIONS))

    def test_three_dimensional_accumulator_and_projection(self):
        edges=([0,.5,1],[0,.3,.6],[0,math.pi,2*math.pi])
        events=[(.2,.1,0.0),(.7,.4,2*math.pi),(.2,.1,-.01),(1.2,.1,0)]
        cells,flow=common.accumulate_nd(events,[1,2,-.5,4],edges,periodic_axes=(2,))
        self.assertEqual(cells[(0,0,0)],[1.0,1.0])
        self.assertEqual(cells[(1,1,0)],[2.0,4.0])
        self.assertEqual(cells[(0,0,1)],[-.5,.25])
        self.assertEqual(flow,4.0)
        projection=common.project_nd(cells,0)
        self.assertEqual(projection[0],[.5,1.25])
