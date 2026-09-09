"""
Nozzle Contour Module
Handles physical (X,Y) coordinates of the Rao Bell nozzle and 2D divergence losses.
"""

import math

try:
    import sympy as sp
except ImportError:
    sp = None


def derive_divergence_loss():
    """
    Symbolically derives the 2D divergence loss factor (lambda) for a conical nozzle.
    Integrates the axial momentum vector across the spherical exit cap.
    
    Returns: sympy Equality object
    """
    if sp is None:
        return "SymPy not installed."
        
    Lambda, alpha = sp.symbols('Lambda alpha')
    
    # Lambda = (1 + cos(alpha)) / 2
    eq_lam = sp.Eq(Lambda, (1 + sp.cos(alpha)) / 2)
    
    return eq_lam


def calculate_divergence_loss_factor(exit_angle_deg):
    """
    Conical-nozzle divergence loss factor, lambda = (1 + cos(alpha)) / 2.

    VALID ONLY FOR CONICAL NOZZLES: it assumes uniform spherical source flow
    at the exit plane. For contoured (Rao bell) nozzles use
    calculate_bell_divergence_loss_factor instead — applying this factor to a
    bell contour underpredicts delivered Isp and biases optimization toward
    conical shapes.
    """
    alpha_rad = math.radians(exit_angle_deg)
    return (1.0 + math.cos(alpha_rad)) / 2.0


def calculate_bell_divergence_loss_factor(theta_n_deg, theta_e_deg):
    """
    Divergence loss factor for a parabolic (Rao-style) bell nozzle, using the
    standard engineering approximation of evaluating the conical formula at the
    MEAN of the initial parabola angle theta_n and the exit angle theta_e:

        lambda_bell = (1 + cos((theta_n + theta_e)/2)) / 2

    Because bell exit angles are shallow (theta_e ~ 2-12 deg), this yields a
    markedly smaller loss than a conical nozzle of equal length. The exact
    value requires a method-of-characteristics exit-plane integration (not yet
    implemented — see knowledge_base.json entry).

    Reference: Sutton & Biblarz, Rocket Propulsion Elements (bell-nozzle
    divergence-loss approximation); Rao, G.V.R., Jet Propulsion 28, 1958.
    """
    if theta_e_deg > theta_n_deg:
        raise ValueError("Bell exit angle must not exceed the initial wall angle.")
    mean_angle_rad = math.radians(0.5 * (theta_n_deg + theta_e_deg))
    return (1.0 + math.cos(mean_angle_rad)) / 2.0


def calculate_displacement_thickness(x_m, velocity_m_s, density_kg_m3, mu_Pa_s):
    """
    Turbulent flat-plate boundary-layer DISPLACEMENT thickness estimate at
    running length x (roadmap 1.2 'Nozzle Boundary Layer Displacement'):

        delta   = 0.37 x / Re_x^0.2      (Schlichting 1/7th-power turbulent BL)
        delta*  = delta / 8              (1/7th-power-profile displacement ratio)

    First-order estimate for nozzle-wall boundary layers; a compressible
    integral method (or CFD) supersedes it at high fidelity.
    Reference: Schlichting, Boundary-Layer Theory, 7th ed.
    """
    if x_m <= 0:
        return 0.0
    Re_x = density_kg_m3 * velocity_m_s * x_m / mu_Pa_s
    if Re_x <= 0:
        return 0.0
    delta = 0.37 * x_m / Re_x ** 0.2
    return delta / 8.0


def calculate_effective_area(A_geometric_m2, wall_radius_m, delta_star_m):
    """
    Boundary-layer-corrected effective flow area of an axisymmetric station:

        A_eff = A_geo - 2*pi*R*delta*

    Ignoring delta* overpredicts the achieved expansion ratio and exit
    pressure (flow-separation misprediction near sea level). Guarded against
    non-physical (delta* so large the annulus closes).
    """
    A_eff = A_geometric_m2 - 2.0 * math.pi * wall_radius_m * delta_star_m
    if A_eff <= 0:
        raise ValueError("Displacement thickness exceeds geometric area — "
                         "boundary-layer model out of range.")
    return A_eff


def calculate_conical_nozzle_length(R_t, R_e, alpha_deg=15.0):
    """
    Calculates the axial length of a pure 15-degree conical nozzle.
    This serves as the benchmark for the fractional bell length.
    """
    alpha_rad = math.radians(alpha_deg)
    L_cone = (R_e - R_t) / math.tan(alpha_rad)
    return L_cone


def generate_rao_bell_coordinates(R_t, R_e, fractional_length=0.8, theta_n_deg=30.0, theta_e_deg=8.0, num_points=50):
    """
    Generates the (X, Y) coordinates for a parabolic Rao bell nozzle using a quadratic Bezier curve.
    
    R_t: Throat radius
    R_e: Exit radius (derived from expansion ratio)
    fractional_length: Fraction of a 15-deg cone (usually 80%)
    theta_n_deg: Initial expansion angle right after the throat
    theta_e_deg: Exit lip angle
    
    Returns: List of (X, Y) tuples from throat to exit.
    """
    L_cone = calculate_conical_nozzle_length(R_t, R_e)
    L_bell = fractional_length * L_cone
    
    # Point 0: Throat
    P0_x = 0.0
    P0_y = R_t
    
    # Point 2: Exit
    P2_x = L_bell
    P2_y = R_e
    
    # Point 1: Control Point (intersection of initial angle tangent and exit angle tangent)
    m1 = math.tan(math.radians(theta_n_deg))
    m2 = math.tan(math.radians(theta_e_deg))
    
    # Line 1: y - P0_y = m1 * (x - P0_x) => y = m1*x + R_t
    # Line 2: y - P2_y = m2 * (x - P2_x) => y = m2*(x - L_bell) + R_e
    
    # Intersection:
    # m1*x + R_t = m2*x - m2*L_bell + R_e
    # x(m1 - m2) = R_e - R_t - m2*L_bell
    P1_x = (R_e - R_t - m2 * L_bell) / (m1 - m2)
    P1_y = m1 * P1_x + R_t
    
    # Generate quadratic bezier points
    coordinates = []
    for i in range(num_points + 1):
        t = i / float(num_points)
        
        # B(t) = (1-t)^2 P0 + 2(1-t)t P1 + t^2 P2
        x = ((1 - t)**2 * P0_x) + (2 * (1 - t) * t * P1_x) + (t**2 * P2_x)
        y = ((1 - t)**2 * P0_y) + (2 * (1 - t) * t * P1_y) + (t**2 * P2_y)
        
        coordinates.append((x, y))
        
    return coordinates
