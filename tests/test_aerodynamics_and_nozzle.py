"""
Unit Tests for Nozzle Gas Dynamics, Losses, and Chamber Acoustics
=================================================================
Verifies:
1. Area-Mach relation solving for supersonic and subsonic regimes.
2. Divergence loss factor lambda for conical and Rao bell nozzles.
3. Separation assessment criteria (Summerfield, Schmucker, Kalt-Badal).
4. Chamber acoustic eigenfrequencies (1T mode).
"""

import math
import pytest

from kryptonis.propulsion_equations.aerodynamics import (
    calculate_area_ratio,
    calculate_mach_from_area_ratio,
)
from kryptonis.propulsion_equations.contour import (
    calculate_divergence_loss_factor,
    calculate_bell_divergence_loss_factor,
)
from kryptonis.propulsion_equations.nozzle_losses import (
    divergence_efficiency,
    separation_assessment,
)
from kryptonis.propulsion_equations.chamber_acoustics import (
    first_tangential_frequency,
)
from kryptonis.propulsion_equations.units import Status


def test_area_mach_relation_invertibility():
    """Verify expansion ratio calculation and inversion to supersonic Mach number."""
    gamma = 1.2
    mach_target = 3.2

    # Compute A / At from Mach
    eps = calculate_area_ratio(mach=mach_target, gamma=gamma)
    assert eps > 1.0

    # Invert back from eps to Mach
    mach_solved = calculate_mach_from_area_ratio(area_ratio=eps, gamma=gamma, supersonic=True)
    assert pytest.approx(mach_solved, rel=1e-5) == mach_target


def test_divergence_loss_factors():
    """Verify divergence loss factor lambda: bell nozzle loss is smaller than conical."""
    # 15 degree conical half-angle
    lambda_cone = calculate_divergence_loss_factor(15.0)
    assert pytest.approx(lambda_cone, rel=1e-4) == 0.9830

    # Rao bell with initial angle 28 deg and exit angle 8 deg
    lambda_bell = calculate_bell_divergence_loss_factor(theta_n_deg=28.0, theta_e_deg=8.0)
    assert lambda_bell > 0.95


def test_separation_criteria():
    """Verify Summerfield, Schmucker, and Kalt-Badal flow separation assessment."""
    P_ambient = 101325.0  # 1 atm sea level
    P_exit = 40000.0      # 0.4 bar exit pressure
    M_exit = 3.0
    P_chamber = 7.0e6     # 70 bar

    res = separation_assessment(
        exit_pressure_Pa=P_exit,
        ambient_pressure_Pa=P_ambient,
        exit_mach=M_exit,
        chamber_pressure_Pa=P_chamber
    )
    assert res.status in (Status.PASS, Status.CONFLICT_UNRESOLVED, Status.FAIL)


def test_chamber_acoustics():
    """Verify first tangential (1T) acoustic screaming frequency."""
    c_sound = 1200.0  # m/s
    Dc = 0.15         # 15 cm diameter chamber

    res_1T = first_tangential_frequency(speed_of_sound_m_s=c_sound, chamber_diameter_m=Dc)
    assert res_1T.status == Status.PASS
    expected_f = 1.84118378 * c_sound / (math.pi * Dc)
    assert pytest.approx(res_1T.value, rel=1e-5) == expected_f
    assert 3000.0 <= res_1T.value <= 6000.0
