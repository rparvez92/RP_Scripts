import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

import RP_sidis_reduced_kinematics as reduced  # noqa: E402


class SidisReducedKinematicsTests(unittest.TestCase):
    def setUp(self):
        self.kin = reduced.calculate_reduced_kinematics(
            beam_energy=10.6716, x=0.25, Q2=3.3, z=0.5, pt=0.25
        )

    def test_simc_constants_are_reproduced(self):
        self.assertAlmostEqual(reduced.ALPHA, 1 / 137.0359895)
        self.assertEqual(reduced.HBARC_GEV_FM, 0.197327053)
        self.assertEqual(reduced.PROTON_MASS_GEV, 0.93827231)
        self.assertEqual(reduced.CHARGED_PION_MASS_GEV, 0.13957018)

    def test_epsilon_forms_are_equivalent(self):
        from_y = reduced.epsilon_from_y(self.kin.y, self.kin.gamma)
        from_angle = reduced.epsilon_from_angle(
            self.kin.theta_e, self.kin.nu, self.kin.Q2
        )
        self.assertAlmostEqual(from_y, from_angle, places=13)
        self.assertAlmostEqual(self.kin.epsilon, 0.590101, places=6)

    def test_hadronic_jacobian_matches_finite_difference(self):
        delta = 1.0e-7
        cosine = self.kin.cos_theta_pq
        momentum2 = self.kin.hadron_momentum ** 2
        pt2_plus = momentum2 * (1.0 - (cosine + delta) ** 2)
        pt2_minus = momentum2 * (1.0 - (cosine - delta) ** 2)
        derivative = abs((pt2_plus - pt2_minus) / (2.0 * delta))
        numerical = derivative / self.kin.nu
        self.assertAlmostEqual(numerical, self.kin.hadronic_jacobian, places=8)

    def test_w1_w2_to_ft_fl_identity(self):
        f1, f2 = 0.43, 0.31
        transverse, longitudinal = reduced.ft_fl_from_f1_f2(
            f1, f2, self.kin.x, self.kin.gamma
        )
        left = reduced.w1_w2_bracket(
            f1, f2, reduced.PROTON_MASS_GEV, self.kin.nu, self.kin.theta_e
        )
        right = self.kin.lab_structure_factor * (
            transverse + self.kin.epsilon * longitudinal
        )
        self.assertAlmostEqual(left, right, places=14)

    def test_reduced_normalization_factorization(self):
        direct = reduced.reduced_normalization(
            self.kin.eprime, self.kin.Q2, self.kin.hadronic_jacobian,
            self.kin.theta_e, self.kin.nu, self.kin.epsilon, self.kin.gamma,
        )
        self.assertAlmostEqual(direct, self.kin.reduced_normalization, places=14)
        self.assertAlmostEqual(
            self.kin.reduced_normalization,
            self.kin.base_normalization * self.kin.lab_structure_factor,
            places=14,
        )
        self.assertGreater(self.kin.reduced_normalization, 0.0)

    def test_cross_section_reduction_round_trip(self):
        value, error = 1.7e-3, 2.1e-4
        reduced_value, reduced_error = reduced.reduce_cross_section(
            value, error, self.kin.reduced_normalization
        )
        reconstructed = reduced.reconstruct_lab_cross_section(
            reduced_value, self.kin.reduced_normalization
        )
        self.assertAlmostEqual(reconstructed, value, places=15)
        self.assertAlmostEqual(reduced_error / reduced_value, error / value, places=15)

    def test_record_is_auditable_and_root_independent(self):
        record = self.kin.as_dict()
        self.assertEqual(record["convention_version"], reduced.CONVENTION_VERSION)
        self.assertEqual(record["beam_energy"], 10.6716)
        self.assertIn("reduced_normalization", record)

    def test_unphysical_kinematics_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "scattered energy"):
            reduced.calculate_reduced_kinematics(3.0, 0.25, 3.3, 0.5, 0.2)
        with self.assertRaisesRegex(ValueError, "pt cannot exceed"):
            reduced.calculate_reduced_kinematics(10.6716, 0.25, 3.3, 0.1, 2.0)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            reduced.reduce_cross_section(1.0, -0.1, 2.0)


if __name__ == "__main__":
    unittest.main()
