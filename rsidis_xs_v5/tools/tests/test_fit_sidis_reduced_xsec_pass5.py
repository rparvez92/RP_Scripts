import math
import sys
import unittest
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(TOOLS))
import RP_fit_sidis_reduced_xsec_pass5 as stage  # noqa: E402


def points():
    rows=[]
    for index,(low,high) in enumerate(zip(stage.PHI_EDGES[:-1],stage.PHI_EDGES[1:]),1):
        phi=.5*(low+high); error=.05
        rows.append({**{field:"test" for field in (*stage.FAMILY_FIELDS,*stage.PROVENANCE_FIELDS)},
                     "Combination_group_id":"g5","Pt_bin":"1","phipq_bin":str(index),
                     "phipq_center":str(phi),"Reduced_data_xsec":str(2+.3*math.cos(phi)-.1*math.cos(2*phi)),
                     "Reduced_data_xsec_err":str(error),"Reduction_status":"OK",
                     "Member_pass":"pass5","Family_status":"LT_READY","epsilon":"0.7",
                     "Pt_low":"0","Pt_high":"0.1","Pt_high_closed":"0"})
    return rows


def coverage():
    bins={}; curves={}
    for row in points():
        key=("g5",1,int(row["phipq_bin"])); bins[key]={"Policy_version":"coverage-v1","Bin_status":"PASS","Bin_reason":"","Survives_all_requirements":"1"}
    curves[("g5",1)]={"Curve_status":"APPROVED","Curve_reason":"","N_selected":"15","Surviving_mask":"1"*15,"Selected_mask":"1"*15,"Min_selected_span_rad":str(2*math.pi)}
    return bins,curves


class ReducedPass5FitTests(unittest.TestCase):
    def test_approved_full_ring_recovers_coefficients(self):
        bins,curves=coverage(); annotated=stage.annotate_points(points(),bins,curves)
        fit=stage.fit_curve(annotated)
        self.assertEqual(fit["Fit_status"],"OK")
        self.assertAlmostEqual(fit["C_0_pass5"],2.0,places=12)
        self.assertAlmostEqual(fit["C_1_pass5"],.3,places=12)
        self.assertAlmostEqual(fit["C_2_pass5"],-.1,places=12)
        self.assertEqual(fit["Weighted_design_rank"],3)

    def test_rejected_coverage_is_not_fitted(self):
        bins,curves=coverage(); curves[("g5",1)]["Curve_status"]="REJECTED"
        annotated=stage.annotate_points(points(),bins,curves); fit=stage.fit_curve(annotated)
        self.assertTrue(all(row["Used_in_fit"]==0 for row in annotated))
        self.assertEqual(fit["Fit_status"],"NOT_FITTED")

    def test_invalid_reduction_point_is_excluded(self):
        rows=points(); rows[0]["Reduction_status"]="INVALID"
        bins,curves=coverage(); annotated=stage.annotate_points(rows,bins,curves)
        self.assertEqual(annotated[0]["Used_in_fit"],0)
        self.assertIn("REDUCTION_POINT_NOT_USABLE",annotated[0]["Fit_point_reason"])


if __name__=="__main__": unittest.main()
