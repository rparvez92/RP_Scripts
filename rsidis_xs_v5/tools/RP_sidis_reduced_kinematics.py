#!/usr/bin/env python3
"""Pure kinematics and normalization for the SIMC-native reduced SIDIS cross section.

All energies and momenta accepted by this module are in GeV.  ``Q2`` denotes
the positive quantity Q^2 in GeV^2.  The implementation mirrors the central,
non-Fermi kinematics used by ``calc_semi_xsec`` and does no file or ROOT I/O.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


CONVENTION_VERSION = "RP_REDUCED_XSEC_SIMC_LAB_V1"

# Exact values used by simc_gfortran/constants.inc, converted from MeV where
# necessary.  Keeping these values here makes numerical comparisons with SIMC
# reproducible rather than dependent on a separately updated constants package.
ALPHA = 1.0 / 137.0359895
HBARC_GEV_FM = 0.197327053
PROTON_MASS_GEV = 0.93827231
CHARGED_PION_MASS_GEV = 0.13957018
MICROBARN_PER_FM2 = 1.0e4


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
    """Return nu = Q^2/(2 M x) in GeV."""
    return _positive("Q2", Q2) / (
        2.0 * _positive("target_mass", target_mass) * _positive("x", x)
    )


def scattered_energy(beam_energy: float, nu: float) -> float:
    """Return E' = E - nu in GeV and require physical scattered energy."""
    result = _positive("beam_energy", beam_energy) - _positive("nu", nu)
    if result <= 0.0:
        raise ValueError("scattered energy E' must be positive")
    return result


def electron_scattering_angle(beam_energy: float, eprime: float,
                              Q2: float) -> float:
    """Return theta_e from Q^2 = 4 E E' sin^2(theta_e/2), in radians."""
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
    """Return gamma = 2 M x / Q."""
    return (2.0 * _positive("target_mass", target_mass) * _positive("x", x)
            / math.sqrt(_positive("Q2", Q2)))


def epsilon_from_y(y: float, gamma: float) -> float:
    """Return virtual-photon epsilon from y and gamma."""
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
    """Return epsilon = [1+2(1+nu^2/Q^2)tan^2(theta_e/2)]^-1."""
    theta_e = _finite("theta_e", theta_e)
    if not 0.0 <= theta_e < math.pi:
        raise ValueError("theta_e must satisfy 0 <= theta_e < pi")
    result = 1.0 / (
        1.0 + 2.0 * (1.0 + _positive("nu", nu) ** 2 / _positive("Q2", Q2))
        * math.tan(theta_e / 2.0) ** 2
    )
    if not 0.0 < result <= 1.0:
        raise ValueError("virtual-photon epsilon must satisfy 0 < epsilon <= 1")
    return result


def hadron_kinematics(z: float, nu: float, pt: float,
                      hadron_mass: float = CHARGED_PION_MASS_GEV
                      ) -> tuple[float, float, float, float]:
    """Return (E_h, p_h, theta_pq, cos(theta_pq)) for central kinematics."""
    z = _positive("z", z)
    nu = _positive("nu", nu)
    pt = _finite("pt", pt)
    mass = _positive("hadron_mass", hadron_mass)
    if pt < 0.0:
        raise ValueError("pt must be nonnegative")
    energy = z * nu
    if energy <= mass:
        raise ValueError("hadron energy z*nu must exceed the hadron mass")
    momentum = math.sqrt(energy * energy - mass * mass)
    if pt > momentum + 1.0e-12:
        raise ValueError("pt cannot exceed the hadron momentum")
    sine = min(1.0, pt / momentum)
    theta = math.asin(sine)
    return energy, momentum, theta, math.cos(theta)


def hadronic_jacobian(hadron_momentum: float, cos_theta_pq: float,
                      nu: float) -> float:
    """Return |d(z,Pt^2)/d(E_h,cos(theta_pq))| in GeV."""
    momentum = _positive("hadron_momentum", hadron_momentum)
    cosine = _finite("cos_theta_pq", cos_theta_pq)
    if not 0.0 < cosine <= 1.0:
        raise ValueError("central SIMC convention requires 0 < cos(theta_pq) <= 1")
    return 2.0 * momentum * momentum * cosine / _positive("nu", nu)


def base_normalization(eprime: float, Q2: float, jacobian: float) -> float:
    """Return (hbar*c)^2 1e4 [4 alpha^2 E'^2/Q^4] J_h, in microbarn*GeV."""
    return (
        HBARC_GEV_FM ** 2 * MICROBARN_PER_FM2
        * 4.0 * ALPHA ** 2 * _positive("eprime", eprime) ** 2
        / _positive("Q2", Q2) ** 2
        * _positive("jacobian", jacobian)
    )


def lab_structure_factor(theta_e: float, nu: float, epsilon: float,
                         gamma: float) -> float:
    """Return cos^2(theta_e/2)/[nu epsilon (1+gamma^2)], in GeV^-1."""
    theta_e = _finite("theta_e", theta_e)
    if not 0.0 <= theta_e < math.pi:
        raise ValueError("theta_e must satisfy 0 <= theta_e < pi")
    epsilon = _finite("epsilon", epsilon)
    gamma = _finite("gamma", gamma)
    if not 0.0 < epsilon <= 1.0:
        raise ValueError("epsilon must satisfy 0 < epsilon <= 1")
    if gamma < 0.0:
        raise ValueError("gamma must be nonnegative")
    return (math.cos(theta_e / 2.0) ** 2
            / (_positive("nu", nu) * epsilon * (1.0 + gamma * gamma)))


def reduced_normalization(eprime: float, Q2: float, jacobian: float,
                          theta_e: float, nu: float, epsilon: float,
                          gamma: float) -> float:
    """Return K_red = K_base*A_lab in microbarn."""
    return base_normalization(eprime, Q2, jacobian) * lab_structure_factor(
        theta_e, nu, epsilon, gamma
    )


def reduce_cross_section(value: float, error: float,
                         normalization: float) -> tuple[float, float]:
    """Divide a laboratory cross section and statistical error by K_red."""
    value = _finite("value", value)
    error = _finite("error", error)
    if error < 0.0:
        raise ValueError("cross-section error must be nonnegative")
    normalization = _positive("normalization", normalization)
    return value / normalization, error / normalization


def reconstruct_lab_cross_section(reduced_value: float,
                                  normalization: float) -> float:
    """Reconstruct the laboratory cross section from a reduced value."""
    return _finite("reduced_value", reduced_value) * _positive(
        "normalization", normalization
    )


def ft_fl_from_f1_f2(f1: float, f2: float, x: float,
                     gamma: float) -> tuple[float, float]:
    """Return F_T=2xF1 and F_L=(1+gamma^2)F2-2xF1."""
    f1 = _finite("f1", f1)
    f2 = _finite("f2", f2)
    x = _positive("x", x)
    gamma = _finite("gamma", gamma)
    if gamma < 0.0:
        raise ValueError("gamma must be nonnegative")
    transverse = 2.0 * x * f1
    longitudinal = (1.0 + gamma * gamma) * f2 - transverse
    return transverse, longitudinal


def w1_w2_bracket(f1: float, f2: float, target_mass: float, nu: float,
                  theta_e: float) -> float:
    """Return W2 cos^2(theta/2)+2 W1 sin^2(theta/2)."""
    f1 = _finite("f1", f1)
    f2 = _finite("f2", f2)
    mass = _positive("target_mass", target_mass)
    nu = _positive("nu", nu)
    theta_e = _finite("theta_e", theta_e)
    if not 0.0 <= theta_e < math.pi:
        raise ValueError("theta_e must satisfy 0 <= theta_e < pi")
    return (f2 / nu * math.cos(theta_e / 2.0) ** 2
            + 2.0 * f1 / mass * math.sin(theta_e / 2.0) ** 2)


@dataclass(frozen=True)
class ReducedKinematics:
    """Central kinematics and normalization for one (x,Q2,z,Pt,E) point."""

    convention_version: str
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
    hadron_energy: float
    hadron_momentum: float
    theta_pq: float
    cos_theta_pq: float
    hadronic_jacobian: float
    base_normalization: float
    lab_structure_factor: float
    reduced_normalization: float

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
    """Calculate every audited central quantity entering K_red."""
    beam_energy = _positive("beam_energy", beam_energy)
    x = _positive("x", x)
    Q2 = _positive("Q2", Q2)
    z = _positive("z", z)
    pt = _finite("pt", pt)
    if pt < 0.0:
        raise ValueError("pt must be nonnegative")
    nu = energy_transfer(x, Q2, target_mass)
    eprime = scattered_energy(beam_energy, nu)
    theta_e = electron_scattering_angle(beam_energy, eprime, Q2)
    y = inelasticity(beam_energy, nu)
    gamma = gamma_factor(x, Q2, target_mass)
    epsilon = epsilon_from_y(y, gamma)
    epsilon_angle = epsilon_from_angle(theta_e, nu, Q2)
    if not math.isclose(epsilon, epsilon_angle, rel_tol=2.0e-12, abs_tol=2.0e-12):
        raise ArithmeticError("inconsistent epsilon calculations")
    hadron_energy, hadron_momentum, theta_pq, cosine = hadron_kinematics(
        z, nu, pt, hadron_mass
    )
    jacobian = hadronic_jacobian(hadron_momentum, cosine, nu)
    base = base_normalization(eprime, Q2, jacobian)
    structure = lab_structure_factor(theta_e, nu, epsilon, gamma)
    return ReducedKinematics(
        convention_version=CONVENTION_VERSION,
        beam_energy=beam_energy, x=x, Q2=Q2, z=z, pt=pt,
        nu=nu, eprime=eprime, theta_e=theta_e, y=y, gamma=gamma,
        epsilon=epsilon, hadron_energy=hadron_energy,
        hadron_momentum=hadron_momentum, theta_pq=theta_pq,
        cos_theta_pq=cosine, hadronic_jacobian=jacobian,
        base_normalization=base, lab_structure_factor=structure,
        reduced_normalization=base * structure,
    )
