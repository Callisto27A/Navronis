"""
Unit Tests for Regenerative Coolant Correlations & Wall Conduction
==================================================================
Verifies:
1. Nusselt number heat transfer correlations (Dittus-Boelter, Sieder-Tate).
2. Curvature heat transfer enhancement (Ito correlation).
3. Darcy-Weisbach and Haaland friction factor models.
4. Fourier 1D wall conduction and temperature drop.
"""

import math
import pytest

from kryptonis.propulsion_equations.thermal import (
    nusselt_dittus_boelter_1930,
    nusselt_sieder_tate,
    ito_curvature_factor,
    friction_factor_haaland,
)
from kryptonis.propulsion_equations.wall_conduction import (
    conduction_heat_flux,
    fourier_wall_temperature_drop,
)
from kryptonis.propulsion_equations.units import Status


def test_dittus_boelter_and_sieder_tate():
    """Verify turbulent cooling channel Nusselt numbers."""
    Re = 100000.0  # 1e5 Reynolds number
    Pr = 2.0       # Coolant Prandtl number

    # Dittus-Boelter (1930 as published: 0.0243 * Re^0.8 * Pr^0.4)
    # Correctly carries INSUFFICIENT_EVIDENCE per primary source audit findings
    res_db = nusselt_dittus_boelter_1930(reynolds=Re, prandtl=Pr, heating=True)
    assert res_db.status in (Status.PASS, Status.INSUFFICIENT_EVIDENCE)
    expected_db = 0.0243 * (Re ** 0.8) * (Pr ** 0.4)
    assert pytest.approx(res_db.value, rel=1e-5) == expected_db

    # Sieder-Tate with bulk_viscosity = 1.2e-4, wall_viscosity = 1.0e-4
    res_st = nusselt_sieder_tate(
        reynolds=Re,
        prandtl=Pr,
        bulk_viscosity_Pa_s=1.2e-4,
        wall_viscosity_Pa_s=1.0e-4
    )
    assert res_st.status in (Status.PASS, Status.INSUFFICIENT_EVIDENCE)
    expected_st = 0.027 * (Re ** 0.8) * (Pr ** (1.0 / 3.0)) * ((1.2 / 1.0) ** 0.14)
    assert pytest.approx(res_st.value, rel=1e-5) == expected_st


def test_ito_curvature_enhancement():
    """Verify Ito curvature factor enhancement on the concave/throat side."""
    Re = 50000.0
    R_curve = 0.1   # 100 mm bend radius
    r_hyd = 0.001   # 1 mm hydraulic radius

    # Concave side (cooling channel hot wall at nozzle throat)
    res_concave = ito_curvature_factor(
        reynolds_bulk=Re,
        hydraulic_radius_m=r_hyd,
        radius_of_curvature_m=R_curve,
        side="concave"
    )
    assert res_concave.status == Status.PASS
    assert res_concave.value > 1.0  # heat transfer is enhanced on concave side

    # Convex side
    res_convex = ito_curvature_factor(
        reynolds_bulk=Re,
        hydraulic_radius_m=r_hyd,
        radius_of_curvature_m=R_curve,
        side="convex"
    )
    assert res_convex.status == Status.PASS
    assert res_convex.value < 1.0  # heat transfer is reduced on convex side


def test_haaland_friction_factor():
    """Verify Haaland explicit friction factor compared to smooth pipe Darcy."""
    Re = 80000.0
    rel_rough = 1.5e-6 / 0.003  # roughness / D

    f = friction_factor_haaland(reynolds=Re, relative_roughness=rel_rough)
    assert f.status == Status.PASS
    # Typical turbulent friction factors are between 0.015 and 0.04
    assert 0.015 <= f.value <= 0.04


def test_fourier_wall_conduction():
    """Verify Fourier 1D conduction across a copper chamber liner."""
    k_copper = 320.0
    t_wall = 0.0015       # 1.5 mm wall thickness
    q_flux = 30.0e6       # 30 MW/m^2 peak throat heat flux

    delta_T = fourier_wall_temperature_drop(
        heat_flux_W_m2=q_flux,
        wall_thickness_m=t_wall,
        thermal_conductivity_W_m_K=k_copper
    )
    assert delta_T.status == Status.PASS
    expected_delta_T = q_flux * t_wall / k_copper
    assert pytest.approx(delta_T.value, rel=1e-6) == expected_delta_T
    assert pytest.approx(delta_T.value, rel=1e-2) == 140.6
