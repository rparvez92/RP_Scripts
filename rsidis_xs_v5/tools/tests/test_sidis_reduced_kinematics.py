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

    def test_simc_constants_and_conventions(self):
        self.assertAlmostEqual(reduced.ALPHA, 1 / 137.0359895)
        self.assertEqual(reduced.PROTON_MASS_GEV, 0.93827231)
        self.assertEqual(reduced.CHARGED_PION_MASS_GEV, 0.13957018)
        self.assertEqual(self.kin.target_mass, reduced.PROTON_MASS_GEV)
        self.assertEqual(self.kin.cross_section_level, "BORN")

    def test_epsilon_forms_are_equivalent(self):
        from_y = reduced.epsilon_from_y(self.kin.y, self.kin.gamma)
        from_angle = reduced.epsilon_from_angle(
            self.kin.theta_e, self.kin.nu, self.kin.Q2
        )
        self.assertAlmostEqual(from_y, from_angle, places=13)
        self.assertAlmostEqual(self.kin.epsilon, 0.590101, places=6)

    def test_invariant_mass_and_hand_energy(self):
        expected_w2 = (reduced.PROTON_MASS_GEV ** 2
                       + 2 * reduced.PROTON_MASS_GEV * self.kin.nu
                       - self.kin.Q2)
        self.assertAlmostEqual(self.kin.W2, expected_w2, places=14)
        expected_k = ((self.kin.W2 - reduced.PROTON_MASS_GEV ** 2)
                      / (2 * reduced.PROTON_MASS_GEV))
        self.assertAlmostEqual(self.kin.hand_equivalent_energy, expected_k, places=14)
        self.assertAlmostEqual(
            self.kin.hand_equivalent_energy, self.kin.nu * (1 - self.kin.x), places=14
        )

    def test_hand_flux_matches_direct_formula(self):
        expected = (reduced.ALPHA / (2 * math.pi ** 2)
                    * self.kin.eprime / self.kin.beam_energy / self.kin.Q2
                    / (1 - self.kin.epsilon)
                    * self.kin.hand_equivalent_energy)
        self.assertAlmostEqual(self.kin.hand_flux, expected, places=16)
        self.assertGreater(self.kin.hand_flux, 0)

    def test_energy_momentum_relation_and_jacobian(self):
        self.assertAlmostEqual(
            self.kin.pion_energy ** 2 - self.kin.pion_momentum ** 2,
            reduced.CHARGED_PION_MASS_GEV ** 2,
            places=14,
        )
        self.assertAlmostEqual(
            self.kin.dEpi_dPpi,
            self.kin.pion_momentum / self.kin.pion_energy,
            places=15,
        )
        delta = 1.0e-6
        energy_plus = math.sqrt((self.kin.pion_momentum + delta) ** 2
                                + reduced.CHARGED_PION_MASS_GEV ** 2)
        energy_minus = math.sqrt((self.kin.pion_momentum - delta) ** 2
                                 + reduced.CHARGED_PION_MASS_GEV ** 2)
        numerical = (energy_plus - energy_minus) / (2 * delta)
        self.assertAlmostEqual(self.kin.dEpi_dPpi, numerical, places=9)

    def test_simc_internal_jacobian_is_audit_only(self):
        expected = (2 * self.kin.pion_momentum ** 2 * self.kin.cos_theta_pq
                    / self.kin.nu)
        self.assertAlmostEqual(
            self.kin.simc_internal_hadronic_jacobian, expected, places=14
        )
        self.assertAlmostEqual(
            self.kin.virtual_photon_reduction_factor,
            self.kin.dEpi_dPpi / self.kin.hand_flux,
            places=14,
        )

    def test_cross_section_conversion_round_trip(self):
        value, error = 1.7e-3, 2.1e-4
        virtual, virtual_error = reduced.convert_lab_to_virtual_photon(
            value, error, self.kin.virtual_photon_reduction_factor
        )
        reconstructed = reduced.reconstruct_lab_cross_section(
            virtual, self.kin.virtual_photon_reduction_factor
        )
        self.assertAlmostEqual(reconstructed, value, places=15)
        self.assertAlmostEqual(virtual_error / virtual, error / value, places=15)

    def test_record_is_auditable_and_root_independent(self):
        record = self.kin.as_dict()
        self.assertEqual(record["convention_version"], reduced.CONVENTION_VERSION)
        self.assertEqual(record["target_mass_convention"], reduced.TARGET_MASS_CONVENTION)
        self.assertIn("hand_flux", record)
        self.assertIn("dEpi_dPpi", record)

    def test_unphysical_kinematics_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "scattered energy"):
            reduced.calculate_reduced_kinematics(3.0, 0.25, 3.3, 0.5, 0.2)
        with self.assertRaisesRegex(ValueError, "pt cannot exceed"):
            reduced.calculate_reduced_kinematics(10.6716, 0.25, 3.3, 0.1, 2.0)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            reduced.convert_lab_to_virtual_photon(1.0, -0.1, 2.0)
        with self.assertRaisesRegex(ValueError, r"W\^2"):
            reduced.invariant_mass_squared(1.0, 3.0)


if __name__ == "__main__":
    unittest.main()
