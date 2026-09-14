"""
Unit Tests for Day 3: Regenerative Cooling Channels & 1D Conjugate Thermal Analysis
=====================================================================================
Verifies:
1. Milled channel circumferential geometry & perimeter packing limits.
2. Coolant hydraulics, Reynolds number, and Haaland friction factor.
3. Fin efficiency and enhanced convective heat transfer area.
4. 1D three-resistance conjugate thermal closure (T_bulk < T_wc < T_wg < T_aw).
5. Darcy-Weisbach pressure drop and coolant bulk temperature rise.
6. Thermal compressive stress and yield safety margins.
7. Input bounds validation and self-explaining literature provenance.
"""

import math
import pytest

from kryptonis.propulsion_equations.regen_channel import (
    RegenCoolingJacket,
    RegenChannelResult,
)


@pytest.fixture
def nominal_jacket():
    return RegenCoolingJacket(
        throat_diameter_m=0.04265,
        chamber_diameter_m=0.09148,
        chamber_length_m=0.2359,
        mass_flow_coolant_kg_s=2.87,
        chamber_pressure_pa=12.0e6,
        gas_recovery_temp_k=3400.0,
        gas_throat_htc_w_m2k=12500.0,
        n_channels=80,
        channel_height_m=0.0018,
        fin_thickness_m=0.0008,
        wall_thickness_m=0.0015,
        coolant_type="CH4",
        liner_material="CuCrZr",
    )


def test_regen_channel_geometry_and_perimeter_fit(nominal_jacket):
    """Verify channel widths, aspect ratios, and hydraulic diameters."""
    res = nominal_jacket.solve()

    assert res.n_channels == 80
    assert res.channel_width_mm > 0.5  # Realistic microchannel width (~0.87 mm)
    assert res.channel_width_mm < 2.0
    assert res.aspect_ratio > 1.0       # Rectangular finned aspect ratio (~2.0)
    assert res.hydraulic_diameter_mm > 0.5
    assert res.hydraulic_diameter_mm < res.channel_width_mm + res.channel_height_mm


def test_regen_coolant_hydraulics(nominal_jacket):
    """Verify coolant flow velocity, Reynolds number, and Haaland friction factor."""
    res = nominal_jacket.solve()

    assert res.coolant_velocity_m_s > 20.0  # Typical rocket coolant speed: 30-80 m/s
    assert res.coolant_velocity_m_s < 120.0
    assert res.reynolds_number > 50000.0   # Fully turbulent
    assert 0.015 < res.friction_factor < 0.08 # Realistic friction factor with roughness
    assert res.prandtl_number > 0.5


def test_fin_efficiency_and_area_enhancement(nominal_jacket):
    """Verify fin efficiency is strictly between 0 and 1, enhancing effective HTC."""
    res = nominal_jacket.solve()

    assert 0.1 < res.fin_efficiency < 1.0
    assert res.enhanced_coolant_htc_W_m2K > res.coolant_htc_W_m2K


def test_conjugate_heat_transfer_closure(nominal_jacket):
    """Verify 3-resistance thermal equilibrium: T_bulk < T_wc < T_wg < T_aw."""
    res = nominal_jacket.solve()

    assert res.hot_gas_wall_temp_K < res.gas_recovery_temp_K
    assert res.coolant_wall_temp_K < res.hot_gas_wall_temp_K
    assert res.coolant_wall_temp_K > 120.0  # Above coolant inlet bulk temperature
    assert res.hot_gas_wall_temp_K < 850.0  # Below CuCrZr maximum service temperature
    assert res.peak_heat_flux_MW_m2 > 10.0  # Typical high-pressure throat flux: 20-50 MW/m2


def test_pressure_drop_and_temperature_rise(nominal_jacket):
    """Verify reasonable coolant delta P and sensible bulk temperature rise."""
    res = nominal_jacket.solve()

    assert 10.0 < res.coolant_pressure_drop_bar < 100.0  # High pressure cooling jacket
    assert 15.0 < res.coolant_temp_rise_K < 150.0        # Meaningful pre-heat without cracking
    assert res.boiling_margin_adequate is True


def test_thermal_stress_and_safety_margin(nominal_jacket):
    """Verify thermal compressive hoop/axial stress and safety margin."""
    res = nominal_jacket.solve()

    assert res.thermal_stress_MPa > 50.0
    assert res.yield_safety_margin > -0.5  # Structurally valid design domain


def test_regen_invalid_inputs_rejected():
    """Verify non-physical channel counts and negative dimensions fail explicitly."""
    # Channel count exceeds throat circumference
    with pytest.raises(ValueError, match=r"exceeds pi \* Dt"):
        RegenCoolingJacket(
            throat_diameter_m=0.04265,
            chamber_diameter_m=0.09148,
            chamber_length_m=0.2359,
            mass_flow_coolant_kg_s=2.87,
            chamber_pressure_pa=12.0e6,
            gas_recovery_temp_k=3400.0,
            gas_throat_htc_w_m2k=12500.0,
            n_channels=300,  # 300 * 0.8 mm = 240 mm > 134 mm perimeter!
            fin_thickness_m=0.0008,
        ).solve()

    # Negative wall thickness
    with pytest.raises(ValueError, match="strictly positive"):
        RegenCoolingJacket(
            throat_diameter_m=0.04265,
            chamber_diameter_m=0.09148,
            chamber_length_m=0.2359,
            mass_flow_coolant_kg_s=2.87,
            chamber_pressure_pa=12.0e6,
            gas_recovery_temp_k=3400.0,
            gas_throat_htc_w_m2k=12500.0,
            wall_thickness_m=-0.001,
        )


def test_regen_explain_provenance(nominal_jacket):
    """Verify .explain() returns exact governing equations and literature citations."""
    res = nominal_jacket.solve()

    exp_wg = res.explain("hot_gas_wall_temp_K")
    assert "T_wg = T_aw - q / h_g" in exp_wg
    assert "Fourier 1D Conduction" in exp_wg

    exp_dp = res.explain("coolant_pressure_drop_bar")
    assert "Darcy-Weisbach" in exp_dp
    assert "NASA SP-8087" in exp_dp

    exp_f = res.explain("friction_factor")
    assert "Haaland, S. E. (1983)" in exp_f
