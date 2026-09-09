"""
Nozzle Aerodynamics Module
Handles 1D isentropic gas dynamics, Area-Mach relations, and thrust coefficients.
"""

import math

try:
    import sympy as sp
except ImportError:
    sp = None


def derive_area_mach_relation():
    """
    Symbolically derives the Area-Mach relationship from 1D mass continuity.
    
    Returns: sympy Equality object
    """
    if sp is None:
        return "SymPy not installed."
        
    A, A_t, M, gamma = sp.symbols('A A_t M gamma')
    
    # Area ratio A / A_t
    term1 = 1 / M
    term2 = (2 + (gamma - 1) * M**2) / (gamma + 1)
    power = (gamma + 1) / (2 * (gamma - 1))
    
    eq_am = sp.Eq(A / A_t, term1 * (term2 ** power))
    
    return eq_am


def calculate_area_ratio(mach, gamma):
    """
    Calculates the expansion ratio (A/A_t) required to achieve a specific Mach number.
    """
    if mach <= 0:
        raise ValueError("Mach number must be > 0.")
        
    term1 = 1.0 / mach
    term2 = (2.0 + (gamma - 1.0) * mach**2) / (gamma + 1.0)
    power = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    
    area_ratio = term1 * math.pow(term2, power)
    return area_ratio


def calculate_mach_from_area_ratio(area_ratio, gamma, supersonic=True):
    """
    Numerically solves the Area-Mach relation to find the local Mach number 
    given a geometric expansion ratio. 
    
    supersonic: True for diverging section (M>1), False for converging section (M<1).
    """
    if area_ratio < 1.0:
        raise ValueError("Area ratio A/A_t cannot be less than 1.0")
        
    # Simple bisection root finder
    def f(M):
        return calculate_area_ratio(M, gamma) - area_ratio
        
    if supersonic:
        low = 1.0
        high = 10.0 # Upper bound for realistic rockets
    else:
        low = 0.001
        high = 1.0
        
    # Bisection
    for _ in range(100):
        mid = (low + high) / 2.0
        if f(mid) * f(low) > 0:
            low = mid
        else:
            high = mid
            
        if (high - low) < 1e-6:
            break
            
    return (low + high) / 2.0


def calculate_thrust_coefficient(gamma, p_c, p_e, p_a, area_ratio):
    """
    Calculates the ideal Thrust Coefficient (Cf).
    
    p_c: Chamber pressure
    p_e: Nozzle exit static pressure
    p_a: Ambient atmospheric pressure
    """
    term1 = (2.0 * gamma**2) / (gamma - 1.0)
    term2 = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0))
    term3 = 1.0 - math.pow(p_e / p_c, (gamma - 1.0) / gamma)
    
    momentum_thrust = math.sqrt(term1 * term2 * term3)
    pressure_thrust = (p_e / p_c - p_a / p_c) * area_ratio
    
    cf = momentum_thrust + pressure_thrust
    return cf
