"""Unit and benchmark tests for bi-centrifugal liquid-liquid swirl coaxial injectors.

Incorporates:
- Abramovich (1963) inviscid swirl core mechanics (phi, mu, characteristic A).
- Soltani et al. (2005) Phase Doppler Anemometry (PDA) experimental benchmark.
- Interfacial liquid-liquid shear atomization Sauter Mean Diameter (SMD).
"""

import math
import pytest

from kryptonis.propulsion_equations.injector import (
    abramovich_phi_from_A,
    abramovich_mu_from_phi,
    calculate_coaxial_swirl_injector,
)


def test_abramovich_inviscid_phi_mu():
    """Verify Abramovich root-solver and discharge coefficient against theoretical limits."""
    # A = 3.0 gives phi ~ 0.4208
    phi_3 = abramovich_phi_from_A(3.0)
    assert 0.40 < phi_3 < 0.45
    
    # Discharge coefficient mu should be ~ 0.217
    mu_3 = abramovich_mu_from_phi(phi_3)
    assert 0.20 < mu_3 < 0.24
    
    # Higher swirl A -> smaller fullness phi (thinner liquid sheet, larger gas core)
    phi_6 = abramovich_phi_from_A(6.0)
    assert phi_6 < phi_3
    
    # Check invalid inputs
    with pytest.raises(ValueError):
        abramovich_phi_from_A(-1.0)
    with pytest.raises(ValueError):
        abramovich_mu_from_phi(2.5)


def test_soltani_bicentrifugal_smd_benchmark():
    """
    Benchmarks the liquid-liquid swirl injector model against PDA experimental data 
    from Soltani et al. (2005) for water cold-flow.
    """
    # Soltani operating conditions for max Re case:
    # Inner: m_dot = 130 kg/hr, delta_p = 0.8 bar
    # Outer: m_dot = 300 kg/hr, delta_p = 3.0 bar
    inputs = {
        "m_dot_inner": 130.0 / 3600.0,  # Convert kg/hr to kg/s (~0.0361 kg/s)
        "delta_p_inner": 0.8e5,         # 0.8 bar to Pa
        "m_dot_outer": 300.0 / 3600.0,  # Convert kg/hr to kg/s (~0.0833 kg/s)
        "delta_p_outer": 3.0e5,         # 3.0 bar to Pa
        "rho_inner": 998.0,             # Water density
        "rho_outer": 998.0,             # Water density
        "A_inner": 3.0,                 # Nominal swirl characteristic
        "A_outer": 3.0,
        "t_wall": 0.001075,             # Based on Soltani (d_o - d_in) / 2
    }
    
    result = calculate_coaxial_swirl_injector(**inputs)
    
    # 1. Test Geometry Sanity Checks
    # Inner orifice radius R_n1 should be ~2.05 mm (diameter ~4.1 mm)
    assert result["inner_stage"]["R_n1_mm"] == pytest.approx(2.05, rel=0.05)
    # Outer orifice radius R_n2 should be ~3.84 mm
    assert result["outer_stage"]["R_n2_mm"] == pytest.approx(3.84, rel=0.05)
    
    # 2. Test SMD Benchmark against Soltani et al. (2005) PDA experimental measurements
    # Target SMD from Soltani et al. (2005) at Re_i = 57,218 and Re_o = 30,292: 138.8 um
    expected_smd_um = 138.8 
    
    assert result["smd_um"] is not None
    # Verifies within 5% of experimental PDA data (paper allows 15% error margin)
    assert result["smd_um"] == pytest.approx(expected_smd_um, rel=0.05)
    
    # 3. Verify provenance
    assert "Soltani et al. (2005)" in result["provenance"]["benchmark_dataset"]
    assert "Abramovich" in result["provenance"]["swirl_theory"]


def test_calculate_coaxial_swirl_cryogenic_hydrocarbon():
    """Verify calculate_coaxial_swirl_injector under representative LOX/RP-1 conditions."""
    res = calculate_coaxial_swirl_injector(
        m_dot_inner=0.5,
        delta_p_inner=0.5e6,
        m_dot_outer=0.2,
        delta_p_outer=0.4e6,
        rho_inner=1141.0,
        rho_outer=810.0,
        A_inner=3.0,
        A_outer=2.5,
        t_wall=0.001,
    )
    assert res["inner_stage"]["R_n1_mm"] > 0
    assert res["outer_stage"]["R_n2_mm"] > res["inner_stage"]["R_n1_mm"]
    assert 50.0 < res["smd_um"] < 300.0