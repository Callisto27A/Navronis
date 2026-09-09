"""
Unit Tests for Thrust Chamber Sizing & Combustion Equations
============================================================
Verifies:
1. Throat area and diameter sizing from mass flow and chamber pressure.
2. Contraction ratio and chamber diameter geometry.
3. Vandenkerckhove isentropic choked flow function Gamma(gamma).
4. Ideal characteristic exhaust velocity c*.
5. Chamber volume and cylindrical length calculations.
"""

import math
import pytest

from kryptonis.propulsion_equations.chamber import (
    throat_area,
    throat_diameter,
    chamber_diameter,
    chamber_volume,
    convergent_volume,
    cylinder_length,
    vandenkerckhove,
    c_star_ideal,
)
from kryptonis.propulsion_equations.combustion import (
    assumed_c_star_efficiency,
    chamber_bulk_residence_time,
)
from kryptonis.propulsion_equations.units import Status


def test_vandenkerckhove_gamma_function():
    """Verify Vandenkerckhove function Gamma(gamma) across known analytical values."""
    gamma = 1.2
    res = vandenkerckhove(gamma)
    assert res.status == Status.PASS
    expected = math.sqrt(gamma) * math.pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)))
    assert pytest.approx(res.value, rel=1e-6) == expected
    assert pytest.approx(res.value, rel=1e-4) == 0.6485


def test_throat_area_and_diameter():
    """Verify throat sizing for a 30 kN class LOX/CH4 engine."""
    m_dot = 10.0
    c_star = 1780.0
    P_c = 7.0e6

    res_At = throat_area(mass_flow_kg_s=m_dot, c_star_m_s=c_star, chamber_pressure_Pa=P_c)
    assert res_At.status == Status.PASS
    expected_At = (m_dot * c_star) / P_c
    assert pytest.approx(res_At.value, rel=1e-7) == expected_At

    res_Dt = throat_diameter(throat_area_m2=res_At.value)
    assert res_Dt.status == Status.PASS
    expected_Dt = math.sqrt(4.0 * expected_At / math.pi)
    assert pytest.approx(res_Dt.value, rel=1e-7) == expected_Dt


def test_chamber_geometry_and_volume():
    """Verify chamber volume Vc = L* * At and cylindrical length Lc."""
    At = 0.0055  # m^2 (~83.6 mm diameter)
    Dt = math.sqrt(4 * At / math.pi)
    L_star = 1.0  # 1.0 m characteristic length
    CR = 3.0     # Contraction ratio Ac / At
    Ac = At * CR

    res_Vc = chamber_volume(l_star_m=L_star, throat_area_m2=At)
    assert res_Vc.status == Status.PASS
    assert pytest.approx(res_Vc.value, rel=1e-7) == At * L_star

    res_Dc = chamber_diameter(throat_diameter_m=Dt, contraction_ratio_=CR)
    assert res_Dc.status == Status.PASS
    assert pytest.approx(res_Dc.value, rel=1e-7) == Dt * math.sqrt(CR)

    # Convergent cone volume with 20 degree half-angle
    half_angle_deg = 20.0
    res_Vconv = convergent_volume(
        throat_diameter_m=Dt,
        chamber_diameter_m=res_Dc.value,
        half_angle_deg=half_angle_deg
    )
    assert res_Vconv.status == Status.PASS
    assert res_Vconv.value > 0.0
    assert res_Vconv.value < res_Vc.value

    # Cylindrical length
    res_Lc = cylinder_length(
        chamber_volume_m3=res_Vc.value,
        convergent_volume_m3=res_Vconv.value,
        chamber_area_m2=Ac
    )
    assert res_Lc.status == Status.PASS
    assert res_Lc.value > 0.0


def test_ideal_characteristic_velocity():
    """Verify c* calculation from molar mass and chamber temperature."""
    gamma = 1.2
    M_w = 0.022  # 22 g/mol = 0.022 kg/mol
    T_c = 3400.0 # K

    res = c_star_ideal(
        gamma=gamma,
        molar_mass_kg_per_mol=M_w,
        chamber_temperature_K=T_c
    )
    assert res.status == Status.PASS
    assert 1500.0 <= res.value <= 2000.0


def test_residence_time():
    """Verify chamber stay time / residence time."""
    Vc = 0.0055
    rho_gas = 5.5  # kg/m^3
    mdot = 10.0    # kg/s

    tau = chamber_bulk_residence_time(
        chamber_volume_m3=Vc,
        density_kg_m3=rho_gas,
        mdot_kg_s=mdot
    )
    assert tau.status == Status.PASS
    expected_tau = (Vc * rho_gas) / mdot
    assert pytest.approx(tau.value, rel=1e-6) == expected_tau
