#!/usr/bin/env python3
"""Central Hand-flux kinematics for the Born SIDIS virtual-photon cross section.

All energies and momenta are in GeV. ``Q2`` is positive Q^2 in GeV^2. The
target mass defaults to the proton mass, matching central ``calc_semi_xsec``
SIDIS kinematics for both LH2 and LD2.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


CONVENTION_VERSION = "RP_VIRTUAL_PHOTON_XSEC_HAND_V1"
TARGET_MASS_CONVENTION = "CENTRAL_STATIONARY_PROTON_MASS_FOR_LH2_AND_LD2"
CROSS_SECTION_LEVEL = "BORN"
RADIATIVE_TREATMENT = "RADIATED_MC_YIELD_WITH_CENTRAL_BORN_MODEL"

ALPHA = 1.0 / 137.0359895
PROTON_MASS_GEV = 0.93827231
CHARGED_PION_MASS_GEV = 0.13957018


def _finite(name: str, value: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _positive(name: str, value: float) -> float:
    result = _finite(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def energy_transfer(x: float, Q2: float,
                    target_mass: float = PROTON_MASS_GEV) -> float:
    """Return nu = Q^2/(2 M x)."""
    return _positive("Q2", Q2) / (
        2.0 * _positive("target_mass", target_mass) * _positive("x", x)
    )


def scattered_energy(beam_energy: float, nu: float) -> float:
    """Return E' = E - nu, requiring a physical scattered electron."""
    result = _positive("beam_energy", beam_energy) - _positive("nu", nu)
    if result <= 0.0:
        raise ValueError("scattered energy E' must be positive")
    return result


def electron_scattering_angle(beam_energy: float, eprime: float,
                              Q2: float) -> float:
    """Return theta_e from Q^2 = 4 E E' sin^2(theta_e/2)."""
    sine_squared = _positive("Q2", Q2) / (
        4.0 * _positive("beam_energy", beam_energy) * _positive("eprime", eprime)
    )
    if sine_squared > 1.0 + 1.0e-12:
        raise ValueError("Q^2 is incompatible with E and E'")
    return 2.0 * math.asin(math.sqrt(min(1.0, sine_squared)))


def inelasticity(beam_energy: float, nu: float) -> float:
    """Return y = nu/E and require 0 < y < 1."""
    result = _positive("nu", nu) / _positive("beam_energy", beam_energy)
    if not 0.0 < result < 1.0:
        raise ValueError("inelasticity y must satisfy 0 < y < 1")
    return result


def gamma_factor(x: float, Q2: float,
                 target_mass: float = PROTON_MASS_GEV) -> float:
    """Return gamma = 2 M x/Q used in the epsilon identity."""
    return (2.0 * _positive("target_mass", target_mass) * _positive("x", x)
            / math.sqrt(_positive("Q2", Q2)))


def epsilon_from_y(y: float, gamma: float) -> float:
    """Return virtual-photon polarization epsilon from y and gamma."""
    y = _finite("y", y)
    gamma = _finite("gamma", gamma)
    if not 0.0 < y < 1.0:
        raise ValueError("inelasticity y must satisfy 0 < y < 1")
    if gamma < 0.0:
        raise ValueError("gamma must be nonnegative")
    gamma2 = gamma * gamma
    numerator = 1.0 - y - 0.25 * gamma2 * y * y
    denominator = 1.0 - y + 0.5 * y * y + 0.25 * gamma2 * y * y
    result = numerator / denominator
    if not 0.0 < result < 1.0:
        raise ValueError("virtual-photon epsilon must satisfy 0 < epsilon < 1")
    return result


def epsilon_from_angle(theta_e: float, nu: float, Q2: float) -> float:
    """Return epsilon from the electron scattering angle."""
    theta_e = _finite("theta_e", theta_e)
    if not 0.0 <= theta_e < math.pi:
        raise ValueError("theta_e must satisfy 0 <= theta_e < pi")
    result = 1.0 / (
        1.0 + 2.0 * (1.0 + _positive("nu", nu) ** 2 / _positive("Q2", Q2))
        * math.tan(theta_e / 2.0) ** 2
    )
    if not 0.0 < result < 1.0:
        raise ValueError("virtual-photon epsilon must satisfy 0 < epsilon < 1")
    return result


def invariant_mass_squared(nu: float, Q2: float,
                           target_mass: float = PROTON_MASS_GEV) -> float:
    """Return W^2 = M^2 + 2 M nu - Q^2."""
    mass = _positive("target_mass", target_mass)
    result = mass * mass + 2.0 * mass * _positive("nu", nu) - _positive("Q2", Q2)
    if result <= mass * mass:
        raise ValueError("W^2 must exceed M^2 for a positive Hand equivalent energy")
    return result


def hand_equivalent_energy(W2: float,
                           target_mass: float = PROTON_MASS_GEV) -> float:
    """Return K_Hand = (W^2-M^2)/(2M)."""
    mass = _positive("target_mass", target_mass)
    result = (_positive("W2", W2) - mass * mass) / (2.0 * mass)
    if result <= 0.0:
        raise ValueError("Hand equivalent photon energy must be positive")
    return result


def hand_flux(beam_energy: float, eprime: float, Q2: float, epsilon: float,
              equivalent_energy: float) -> float:
    """Return Gamma_Hand in GeV^-1."""
    beam = _positive("beam_energy", beam_energy)
    eprime = _positive("eprime", eprime)
    epsilon = _finite("epsilon", epsilon)
    if not 0.0 < epsilon < 1.0:
        raise ValueError("epsilon must satisfy 0 < epsilon < 1")
    return (ALPHA / (2.0 * math.pi ** 2) * eprime / beam
            / _positive("Q2", Q2) / (1.0 - epsilon)
            * _positive("equivalent_energy", equivalent_energy))


def hadron_kinematics(z: float, nu: float, pt: float,
                      hadron_mass: float = CHARGED_PION_MASS_GEV
                      ) -> tuple[float, float, float, float]:
    """Return (E_pi,p_pi,theta_pq,cos(theta_pq)) at central kinematics."""
    energy = _positive("z", z) * _positive("nu", nu)
    mass = _positive("hadron_mass", hadron_mass)
    pt = _finite("pt", pt)
    if pt < 0.0:
        raise ValueError("pt must be nonnegative")
    if energy <= mass:
        raise ValueError("pion energy z*nu must exceed the pion mass")
    momentum = math.sqrt(energy * energy - mass * mass)
    if pt > momentum + 1.0e-12:
        raise ValueError("pt cannot exceed the pion momentum")
    theta = math.asin(min(1.0, pt / momentum))
    return energy, momentum, theta, math.cos(theta)


def energy_to_momentum_jacobian(pion_momentum: float,
                                pion_energy: float) -> float:
    """Return dE_pi/dP_pi = P_pi/E_pi (dimensionless)."""
    momentum = _positive("pion_momentum", pion_momentum)
    energy = _positive("pion_energy", pion_energy)
    if momentum >= energy:
        raise ValueError("massive-pion momentum must be smaller than its energy")
    return momentum / energy


def simc_internal_hadronic_jacobian(pion_momentum: float,
                                    cos_theta_pq: float, nu: float) -> float:
    """Return SIMC's internal d(z,Pt^2)/d(E_pi,cos theta) Jacobian in GeV."""
    momentum = _positive("pion_momentum", pion_momentum)
    cosine = _finite("cos_theta_pq", cos_theta_pq)
    if not 0.0 < cosine <= 1.0:
        raise ValueError("central convention requires 0 < cos(theta_pq) <= 1")
    return 2.0 * momentum * momentum * cosine / _positive("nu", nu)


def virtual_photon_reduction_factor(dEpi_dPpi: float,
                                    gamma_hand: float) -> float:
    """Return (dE_pi/dP_pi)/Gamma_Hand in GeV."""
    return _positive("dEpi_dPpi", dEpi_dPpi) / _positive("gamma_hand", gamma_hand)


def convert_lab_to_virtual_photon(value: float, error: float,
                                  reduction_factor: float) -> tuple[float, float]:
    """Convert a Born lab dE_pi cross section and error to a gamma* cross section."""
    value = _finite("value", value)
    error = _finite("error", error)
    if error < 0.0:
        raise ValueError("cross-section error must be nonnegative")
    factor = _positive("reduction_factor", reduction_factor)
    return value * factor, error * factor


def reconstruct_lab_cross_section(virtual_photon_value: float,
                                  reduction_factor: float) -> float:
    """Reconstruct the Born lab dE_pi cross section."""
    return _finite("virtual_photon_value", virtual_photon_value) / _positive(
        "reduction_factor", reduction_factor
    )


@dataclass(frozen=True)
class ReducedKinematics:
    """Audited central kinematics for one Hand-flux reduction point."""

    convention_version: str
    target_mass_convention: str
    cross_section_level: str
    radiative_treatment: str
    target_mass: float
    beam_energy: float
    x: float
    Q2: float
    z: float
    pt: float
    nu: float
    eprime: float
    theta_e: float
    y: float
    gamma: float
    epsilon: float
    W2: float
    hand_equivalent_energy: float
    hand_flux: float
    pion_energy: float
    pion_momentum: float
    theta_pq: float
    cos_theta_pq: float
    dEpi_dPpi: float
    simc_internal_hadronic_jacobian: float
    virtual_photon_reduction_factor: float

    def as_dict(self) -> dict[str, float | str]:
        return asdict(self)


def calculate_reduced_kinematics(
    beam_energy: float,
    x: float,
    Q2: float,
    z: float,
    pt: float,
    target_mass: float = PROTON_MASS_GEV,
    hadron_mass: float = CHARGED_PION_MASS_GEV,
) -> ReducedKinematics:
    """Calculate the complete central Hand-flux conversion."""
    beam_energy = _positive("beam_energy", beam_energy)
    x = _positive("x", x)
    Q2 = _positive("Q2", Q2)
    z = _positive("z", z)
    pt = _finite("pt", pt)
    if pt < 0.0:
        raise ValueError("pt must be nonnegative")
    target_mass = _positive("target_mass", target_mass)
    nu = energy_transfer(x, Q2, target_mass)
    eprime = scattered_energy(beam_energy, nu)
    theta_e = electron_scattering_angle(beam_energy, eprime, Q2)
    y = inelasticity(beam_energy, nu)
    gamma = gamma_factor(x, Q2, target_mass)
    epsilon = epsilon_from_y(y, gamma)
    epsilon_angle = epsilon_from_angle(theta_e, nu, Q2)
    if not math.isclose(epsilon, epsilon_angle, rel_tol=2.0e-12, abs_tol=2.0e-12):
        raise ArithmeticError("inconsistent epsilon calculations")
    W2 = invariant_mass_squared(nu, Q2, target_mass)
    equivalent_energy = hand_equivalent_energy(W2, target_mass)
    flux = hand_flux(beam_energy, eprime, Q2, epsilon, equivalent_energy)
    pion_energy, pion_momentum, theta_pq, cosine = hadron_kinematics(
        z, nu, pt, hadron_mass
    )
    dEpi_dPpi = energy_to_momentum_jacobian(pion_momentum, pion_energy)
    internal_jacobian = simc_internal_hadronic_jacobian(
        pion_momentum, cosine, nu
    )
    reduction_factor = virtual_photon_reduction_factor(dEpi_dPpi, flux)
    return ReducedKinematics(
        convention_version=CONVENTION_VERSION,
        target_mass_convention=TARGET_MASS_CONVENTION,
        cross_section_level=CROSS_SECTION_LEVEL,
        radiative_treatment=RADIATIVE_TREATMENT,
        target_mass=target_mass, beam_energy=beam_energy, x=x, Q2=Q2,
        z=z, pt=pt, nu=nu, eprime=eprime, theta_e=theta_e, y=y,
        gamma=gamma, epsilon=epsilon, W2=W2,
        hand_equivalent_energy=equivalent_energy, hand_flux=flux,
        pion_energy=pion_energy, pion_momentum=pion_momentum,
        theta_pq=theta_pq, cos_theta_pq=cosine,
        dEpi_dPpi=dEpi_dPpi,
        simc_internal_hadronic_jacobian=internal_jacobian,
        virtual_photon_reduction_factor=reduction_factor,
    )
