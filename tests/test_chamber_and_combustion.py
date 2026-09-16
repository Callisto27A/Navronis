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


def test_storable_propellant_l_star_mappings():
    """Verify L* resolution for storable/hypergolic fuels according to SP-125 Table 4-1."""
    from kryptonis.propulsion_equations.chamber import characteristic_length

    # MMH and N2H4 map onto N2O4/hydrazine-base (30-35 inches -> midpoint 32.5 in = 0.8255 m)
    res_mmh = characteristic_length(fuel="MMH")
    assert res_mmh.status == Status.UNVALIDATED_ASSUMPTION
    assert pytest.approx(res_mmh.value, rel=1e-4) == 32.5 * 0.0254

    res_n2h4 = characteristic_length(fuel="N2H4")
    assert res_n2h4.status == Status.UNVALIDATED_ASSUMPTION
    assert pytest.approx(res_n2h4.value, rel=1e-4) == 32.5 * 0.0254

    # Ethanol has no SP-125 Table 4-1 row, must return INSUFFICIENT_EVIDENCE
    res_eth = characteristic_length(fuel="ETHANOL")
    assert res_eth.status == Status.INSUFFICIENT_EVIDENCE
    assert math.isnan(res_eth.value)
    assert "ETHANOL" in res_eth.notes


def test_storable_propellants_combustor_design():
    """Verify closed-form analytical sizing with N2O4/MMH, Hydrazine, and N2O/Ethanol."""
    from kryptonis.propulsion_equations.combustor import CombustorDesign

    # 1. N2O4 / MMH Hypergolic in-space / RCS engine
    des_mmh = CombustorDesign(
        thrust=10000.0,
        chamber_pressure=15.0e5,
        mixture_ratio=1.65,
        propellant="N2O4/MMH",
    )
    res_mmh = des_mmh.solve()
    assert res_mmh.throat_diameter > 0.0
    assert res_mmh.chamber_diameter > res_mmh.throat_diameter
    assert 1600.0 <= res_mmh.c_star_ideal <= 1850.0
    assert pytest.approx(res_mmh.design.mixture_ratio, rel=1e-4) == 1.65

    # 2. Monopropellant Hydrazine catalytic thruster (Tc ~ 1200 K)
    des_n2h4 = CombustorDesign(
        thrust=500.0,
        chamber_pressure=10.0e5,
        mixture_ratio=0.0,
        propellant="HYDRAZINE",
    )
    res_n2h4 = des_n2h4.solve()
    assert res_n2h4.throat_diameter > 0.0
    assert 1200.0 <= res_n2h4.c_star_ideal <= 1450.0
    assert res_n2h4.oxidizer_mass_flow == 0.0
    assert res_n2h4.fuel_mass_flow > 0.0

    # 3. Green storable N2O / Ethanol thruster
    des_green = CombustorDesign(
        thrust=5000.0,
        chamber_pressure=25.0e5,
        mixture_ratio=4.5,
        propellant="N2O/ETHANOL",
    )
    res_green = des_green.solve()
    assert res_green.throat_diameter > 0.0
    assert 1500.0 <= res_green.c_star_ideal <= 1750.0


def test_storable_propellants_cli_injector_sizing():
    """Verify injector sizing CLI handler accepts storable propellants."""
    from kryptonis.propulsion_equations.cli import run_injector_sizing

    ret_mmh = run_injector_sizing(thrust_n=4000.0, pc_bar=15.0, propellants="N2O4/MMH", injector_type="swirl")
    assert ret_mmh == 0

    ret_green = run_injector_sizing(thrust_n=3000.0, pc_bar=20.0, propellants="N2O/ETHANOL", injector_type="pintle")
    assert ret_green == 0

