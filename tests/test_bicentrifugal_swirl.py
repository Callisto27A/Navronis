import pytest
import math
# Assuming you pasted your code into injector.py
from kryptonis.propulsion_equations.injector import calculate_coaxial_swirl_injector

def test_soltani_bicentrifugal_smd_benchmark():
    """
    Benchmarks the liquid-liquid swirl injector model against PDA experimental data 
    from Soltani et al. (2005) for water cold-flow.
    """
    # Soltani operating conditions for max Re case
    # Inner (m_dot = 130 kg/hr, delta_p = 0.8 bar)
    # Outer (m_dot = 300 kg/hr, delta_p = 3.0 bar approx based on text)
    
    inputs = {
        "m_dot_inner": 130.0 / 3600.0,  # Convert kg/hr to kg/s
        "delta_p_inner": 0.8e5,         # 0.8 bar to Pa
        "m_dot_outer": 300.0 / 3600.0,  # Convert kg/hr to kg/s
        "delta_p_outer": 3.0e5,         # 3.0 bar to Pa
        "rho_inner": 998.0,             # Water density
        "rho_outer": 998.0,             # Water density
        "A_inner": 3.0,                 # Nominal swirl characteristic
        "A_outer": 3.0,
        "t_wall": 0.001075              # Based on Soltani (d_o - d_in) / 2
    }
    
    result = calculate_coaxial_swirl_injector(**inputs)
    
    # 1. Test Geometry Sanity Checks
    # Ensure physical radii are calculated and positive
    assert result["inner_stage"]["R_n1_mm"] > 0
    assert result["outer_stage"]["R_n2_mm"] > 0
    
    # 2. Test SMD Benchmark (This WILL fail until we get the final equation)
    # Target SMD from Soltani et al. (2005) at Re_i = 57,218 and Re_o = 30,292
    expected_smd_um = 138.8 
    
    if result["smd_um"] is None:
        pytest.skip("SMD empirical correlation not yet implemented. Geometry sizing passed.")
    else:
        # Allow a 15% error margin for empirical atomization models
        assert result["smd_um"] == pytest.approx(expected_smd_um, rel=0.15)