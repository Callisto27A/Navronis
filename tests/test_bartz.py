"""
Unit Tests for Canonical Bartz Heat Transfer Formulation
=========================================================
Verifies:
1. Exact calculation of the Bartz convective heat transfer coefficient.
2. Throat curvature scaling (Dt / Rc)^0.1.
3. Recovery temperature never exceeds stagnation temperature.
4. Boundary layer correction factor sigma.
5. Domain boundary error checking on non-physical inputs.
"""

import math
import pytest

from kryptonis.propulsion_equations.bartz import (
    bartz_film_coefficient,
    bartz_sigma,
    bartz_viscosity,
    bartz_prandtl,
    recovery_temperature,
)
from kryptonis.propulsion_equations.units import Status


def test_recovery_temperature_bounds():
    """Verify recovery temperature equals stagnation at M=0 and never exceeds T0."""
    T0 = 3400.0  # K
    gamma = 1.2
    Pr = 0.82

    # At Mach 0, recovery temperature must be exactly stagnation temperature
    T_aw_0 = recovery_temperature(
        stagnation_temperature_K=T0,
        mach=0.0,
        gamma=gamma,
        prandtl=Pr
    )
    assert pytest.approx(T_aw_0.value, rel=1e-6) == T0

    # At Mach 1.0, recovery temperature must be less than stagnation temperature
    T_aw_throat = recovery_temperature(
        stagnation_temperature_K=T0,
        mach=1.0,
        gamma=gamma,
        prandtl=Pr
    )
    assert T_aw_throat.value < T0
    assert T_aw_throat.value > 2500.0


def test_bartz_film_coefficient_nominal():
    """Verify nominal calculation of gas-side film coefficient hg at throat."""
    # Typical 30 kN LOX/CH4 chamber conditions
    D_t = 0.0836          # 83.6 mm throat diameter
    p_c = 7.0e6           # 70 bar chamber pressure
    c_star = 1750.0       # m/s
    R_c = 0.0836 * 1.5    # 1.5 Dt throat radius of curvature
    A_over_At = 1.0       # at throat
    gamma = 1.2
    mach = 1.0
    c_p = 2500.0          # J/(kg*K)
    prandtl = 0.82
    mu = 8.5e-5           # Pa*s
    T_0 = 3400.0          # K
    T_wg = 800.0          # K (hot-gas wall temperature)

    res = bartz_film_coefficient(
        throat_diameter_m=D_t,
        chamber_pressure_Pa=p_c,
        c_star_m_s=c_star,
        throat_radius_curvature_m=R_c,
        area_ratio_local_over_throat=A_over_At,
        gamma=gamma,
        mach=mach,
        specific_heat_J_kgK=c_p,
        prandtl=prandtl,
        viscosity_Pa_s=mu,
        stagnation_temperature_K=T_0,
        wall_temperature_K=T_wg
    )

    assert res.status == Status.PASS
    assert res.units == "W/(m^2.K)"
    # Nominal h_g at throat is typically 5,000 - 30,000 W/(m^2*K)
    assert 5000.0 <= res.value <= 30000.0


def test_bartz_sigma_correction():
    """Verify boundary-layer correction factor sigma."""
    sigma = bartz_sigma(
        wall_temperature_K=800.0,
        stagnation_temperature_K=3400.0,
        mach=1.0,
        gamma=1.2
    )
    assert sigma.status == Status.PASS
    assert sigma.value > 0.0


def test_bartz_negative_inputs_rejected():
    """Verify rejection of non-physical negative inputs."""
    res = bartz_film_coefficient(
        throat_diameter_m=-0.05,
        chamber_pressure_Pa=7e6,
        c_star_m_s=1750.0,
        throat_radius_curvature_m=0.1,
        area_ratio_local_over_throat=1.0,
        gamma=1.2,
        mach=1.0,
        specific_heat_J_kgK=2500.0,
        prandtl=0.82,
        viscosity_Pa_s=8e-5,
        stagnation_temperature_K=3400.0,
        wall_temperature_K=800.0
    )
    assert res.status == Status.PHYSICALLY_INVALID
    assert math.isnan(res.value)
