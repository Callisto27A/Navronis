"""
Unit tests for Component 2: Injector Head & Atomization Elements.
================================================================
Verifies orifice continuity, chugging stiffness criteria, shear coaxial momentum flux ratios,
pintle TMR spray angles, and unlike impinging doublet Rupe balance.
"""

import math
import pytest

from kryptonis.propulsion_equations.injector import (
    orifice_area,
    orifice_diameter,
    orifice_velocity,
    injector_pressure_drop_stiffness,
    size_shear_coaxial,
    size_swirl_coaxial,
    size_pintle_injector,
    size_impinging_doublet,
    InjectorDesign,
    CITATIONS,
)


def test_orifice_hydraulics():
    """Verify Torricelli-Bernoulli orifice continuity and jet velocity."""
    m_dot = 2.5     # kg/s
    rho = 1000.0    # kg/m³
    delta_p = 5.0e5 # 5 bar
    cd = 0.80

    # Theoretical jet velocity V = Cd * sqrt(2 * dP / rho) = 0.8 * sqrt(1000) = 25.298 m/s
    v_jet = orifice_velocity(cd, delta_p, rho)
    assert v_jet == pytest.approx(0.80 * math.sqrt(2.0 * 5.0e5 / 1000.0), rel=1e-5)

    # Area A = m_dot / (Cd * sqrt(2 * rho * dP))
    a = orifice_area(m_dot, cd, rho, delta_p)
    assert a > 0.0
    d = orifice_diameter(a)
    assert d > 0.0
    assert d == pytest.approx(math.sqrt(4.0 * a / math.pi), rel=1e-5)

    # Negative inputs fail closed
    with pytest.raises(ValueError):
        orifice_area(-1.0, cd, rho, delta_p)
    with pytest.raises(ValueError):
        orifice_velocity(cd, -100.0, rho)


def test_chugging_stiffness_margin():
    """Verify NASA SP-194 chugging decoupling criterion (dP / Pc >= 0.15)."""
    pc = 30.0e5  # 30 bar
    dp_safe = 6.0e5  # 6 bar (20%)
    dp_risky = 1.0e5 # 1 bar (3.3%)

    ratio_safe = injector_pressure_drop_stiffness(dp_safe, pc)
    assert ratio_safe == pytest.approx(0.20)
    assert ratio_safe >= 0.15  # Meets stiffness requirement

    ratio_risky = injector_pressure_drop_stiffness(dp_risky, pc)
    assert ratio_risky < 0.15


def test_shear_coaxial_sizing():
    """Size 19-element LOX/CH4 shear coaxial injector head."""
    # 30 kN engine at Pc = 30 bar: m_dot ~ 10.5 kg/s, O/F = 3.5
    m_dot_ox = 8.167    # kg/s LOX
    m_dot_fuel = 2.333  # kg/s Methane
    rho_ox = 1141.0     # kg/m³ (LOX)
    rho_fuel = 422.0    # kg/m³ (liquid CH4)
    dp_ox = 6.0e5       # 6 bar
    dp_fuel = 6.0e5     # 6 bar

    res = size_shear_coaxial(
        m_dot_ox=m_dot_ox,
        m_dot_fuel=m_dot_fuel,
        rho_ox=rho_ox,
        rho_fuel=rho_fuel,
        delta_p_ox=dp_ox,
        delta_p_fuel=dp_fuel,
        n_elements=19,
        cd_ox=0.75,
        cd_fuel=0.80,
    )

    assert res["n_elements"] == 19
    assert res["post_id_mm"] > 1.0
    assert res["post_od_mm"] > res["post_id_mm"]
    assert res["annulus_id_mm"] > res["post_od_mm"]
    assert res["annular_gap_mm"] > 0.1

    # Velocity checks
    assert 10.0 < res["v_ox_m_s"] < 60.0
    assert 15.0 < res["v_fuel_m_s"] < 100.0

    # Momentum flux ratio J = (rho_g * V_g^2) / (rho_l * V_l^2)
    assert res["momentum_flux_ratio_J"] > 0.0
    assert res["velocity_ratio_VR"] > 1.0

    # Droplet SMD
    assert 10.0 < res["smd_um"] < 250.0

    # Literature provenance check
    assert "Yang" in res["provenance"]["sizing"]
    assert "Lefebvre" in res["provenance"]["atomization"]


def test_swirl_coaxial_sizing():
    """Verify RD-170/NK-33 centrifugal swirl coaxial injector sizing."""
    m_dot_ox = 8.167
    m_dot_fuel = 2.333
    rho_ox = 1141.0
    rho_fuel = 422.0
    dp = 6.0e5

    res = size_swirl_coaxial(
        m_dot_liquid=m_dot_ox,
        m_dot_gas=m_dot_fuel,
        rho_liquid=rho_ox,
        rho_gas=rho_fuel,
        delta_p_liquid=dp,
        delta_p_gas=dp,
        n_elements=19,
        geometric_swirl_k=3.0,
    )

    assert res["n_elements"] == 19
    assert res["orifice_diameter_mm"] > 1.0
    assert res["gas_core_diameter_mm"] > 0.5
    assert res["liquid_film_thickness_mm"] > 0.05
    assert 20.0 < res["spray_half_angle_deg"] < 80.0
    assert res["geometric_swirl_K"] == 3.0
    assert 5.0 < res["smd_um"] < 150.0
    assert "Bazarov" in res["provenance"]["swirl_mechanics"]
    assert "Lefebvre" in res["provenance"]["atomization"]


def test_pintle_injector_sizing():
    """Verify Apollo/Merlin style central pintle injector sizing."""
    m_dot_fuel_ann = 2.333  # kg/s annular fuel
    m_dot_ox_rad = 8.167    # kg/s radial ox
    rho_fuel = 422.0
    rho_ox = 1141.0
    dp = 6.0e5
    d_pintle = 0.022  # 22 mm

    res = size_pintle_injector(
        m_dot_annular=m_dot_fuel_ann,
        m_dot_radial=m_dot_ox_rad,
        rho_annular=rho_fuel,
        rho_radial=rho_ox,
        delta_p_annular=dp,
        delta_p_radial=dp,
        pintle_diameter=d_pintle,
        n_radial_slots=16,
    )

    assert res["pintle_diameter_mm"] == 22.0
    assert res["annular_gap_thickness_mm"] > 0.05
    assert res["radial_slot_height_mm"] > 0.1
    assert res["total_momentum_ratio_TMR"] > 0.0

    # Spray angle beta is physically bounded between 20 deg and 75 deg
    assert 20.0 < res["spray_half_angle_deg"] < 75.0
    assert "Dressler" in res["provenance"]["tmr"]


def test_impinging_doublet_sizing():
    """Verify unlike doublet impingement and Rupe momentum balance."""
    m_dot_1 = 8.167
    m_dot_2 = 2.333
    rho_1 = 1141.0
    rho_2 = 422.0
    dp = 6.0e5

    res = size_impinging_doublet(
        m_dot_1=m_dot_1,
        m_dot_2=m_dot_2,
        rho_1=rho_1,
        rho_2=rho_2,
        delta_p_1=dp,
        delta_p_2=dp,
        n_elements=16,
        impingement_half_angle_deg=30.0,
    )

    assert res["n_elements"] == 16
    assert res["orifice_diameter_1_mm"] > 0.5
    assert res["orifice_diameter_2_mm"] > 0.5
    assert res["free_jet_length_mm"] > 0.0
    assert res["impingement_total_angle_deg"] == 60.0
    assert res["smd_um"] > 5.0
    assert "Rupe" in res["provenance"]["rupe_criterion"]


def test_injector_design_facade():
    """Verify high-level InjectorDesign class across all modes."""
    # Coaxial
    coax = InjectorDesign(
        injector_type="coaxial",
        chamber_pressure=30.0e5,
        mass_flow_ox=8.0,
        mass_flow_fuel=2.5,
        rho_ox=1141.0,
        rho_fuel=422.0,
        delta_p_ratio=0.20,
        n_elements=19,
    ).solve()
    assert coax["injector_type"] == "coaxial"
    assert coax["chugging_margin_adequate"] is True
    assert coax["delta_p_bar"] == pytest.approx(6.0)

    # Swirl
    swirl = InjectorDesign(
        injector_type="swirl",
        chamber_pressure=30.0e5,
        mass_flow_ox=8.0,
        mass_flow_fuel=2.5,
        rho_ox=1141.0,
        rho_fuel=422.0,
        delta_p_ratio=0.20,
        n_elements=19,
    ).solve()
    assert swirl["injector_type"] == "swirl"
    assert swirl["geometric_swirl_K"] == 3.0
    assert swirl["spray_half_angle_deg"] > 30.0

    # Pintle
    pintle = InjectorDesign(
        injector_type="pintle",
        chamber_pressure=30.0e5,
        mass_flow_ox=8.0,
        mass_flow_fuel=2.5,
        rho_ox=1141.0,
        rho_fuel=422.0,
        delta_p_ratio=0.20,
        pintle_diameter=0.024,
    ).solve()
    assert pintle["injector_type"] == "pintle"
    assert pintle["total_momentum_ratio_TMR"] > 0.0

    # Impinging
    imp = InjectorDesign(
        injector_type="impinging",
        chamber_pressure=30.0e5,
        mass_flow_ox=8.0,
        mass_flow_fuel=2.5,
        rho_ox=1141.0,
        rho_fuel=422.0,
        delta_p_ratio=0.20,
        n_elements=16,
    ).solve()
    assert imp["injector_type"] == "impinging"
    assert imp["rupe_momentum_parameter"] > 0.0

    # Low delta_p_ratio fails closed
    with pytest.raises(ValueError, match="chugging"):
        InjectorDesign(
            injector_type="coaxial",
            chamber_pressure=30.0e5,
            mass_flow_ox=8.0,
            mass_flow_fuel=2.5,
            rho_ox=1141.0,
            rho_fuel=422.0,
            delta_p_ratio=0.02, # dangerously low
        ).solve()
