"""
Example 01: Liquid Rocket Thrust Chamber Analytical Sizing
===========================================================
Sizes a 30 kN class LOX/CH4 thrust chamber from first-principles analytical
equations (NASA SP-125 / Huzel & Huang):
- Throat area (At) and throat diameter (Dt)
- Contraction ratio (eps_c), chamber area (Ac), and chamber diameter (Dc)
- Characteristic length (L*) and total chamber volume (Vc)
- Convergent cone volume and cylindrical chamber length (Lc)
- Bulk residence / stay time (tau_s)
"""

import math
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
from kryptonis.propulsion_equations.combustion import chamber_bulk_residence_time

def main():
    print("=" * 65)
    print("Kryptonis: Liquid Rocket Thrust Chamber Analytical Sizing")
    print("=" * 65)

    # 1. Mission Design Inputs
    thrust_N = 30000.0          # 30 kN vacuum thrust
    P_c = 7.0e6                 # 7.0 MPa (70 bar) chamber pressure
    gamma = 1.2                 # Specific heat ratio of LOX/CH4 combustion products
    M_w = 0.022                 # Molar mass: 22 g/mol = 0.022 kg/mol
    T_c = 3450.0                # Chamber combustion stagnation temperature [K]
    CR = 3.0                    # Contraction ratio Ac / At
    L_star = 1.05               # Characteristic length [m] (typical LOX/CH4: 0.9 - 1.2 m)

    # 2. Compute Vandenkerckhove Function Gamma(gamma)
    gamma_fn = vandenkerckhove(gamma)
    print(f"Combustion Gas Gamma:     {gamma:.2f}")
    print(f"Vandenkerckhove Gamma:    {gamma_fn.value:.4f}")

    # 3. Compute Ideal Characteristic Exhaust Velocity c*
    res_cstar = c_star_ideal(
        gamma=gamma,
        molar_mass_kg_per_mol=M_w,
        chamber_temperature_K=T_c
    )
    c_star = res_cstar.value
    print(f"Ideal c* Velocity:        {c_star:.1f} m/s")

    # 4. Total Propellant Mass Flow Rate
    # Assuming ideal vacuum Cf ~ 1.75
    Cf_est = 1.75
    At_est = thrust_N / (P_c * Cf_est)
    m_dot = (P_c * At_est) / c_star
    print(f"Estimated Mass Flow Rate: {m_dot:.2f} kg/s")

    # 5. Throat Sizing (At, Dt)
    res_At = throat_area(mass_flow_kg_s=m_dot, c_star_m_s=c_star, chamber_pressure_Pa=P_c)
    res_Dt = throat_diameter(throat_area_m2=res_At.value)
    At = res_At.value
    Dt = res_Dt.value
    print(f"\n--- Throat Geometry ---")
    print(f"Throat Area (At):         {At * 1e4:.2f} cm^2 ({At:.6e} m^2)")
    print(f"Throat Diameter (Dt):     {Dt * 1e3:.2f} mm")

    # 6. Chamber Sizing (Ac, Dc)
    res_Dc = chamber_diameter(throat_diameter_m=Dt, contraction_ratio_=CR)
    Dc = res_Dc.value
    Ac = At * CR
    print(f"\n--- Chamber Geometry ---")
    print(f"Contraction Ratio (eps_c):{CR:.1f}")
    print(f"Chamber Area (Ac):        {Ac * 1e4:.2f} cm^2")
    print(f"Chamber Diameter (Dc):    {Dc * 1e3:.2f} mm")

    # 7. Chamber Volume (Vc) and Convergent Volume (Vconv)
    res_Vc = chamber_volume(l_star_m=L_star, throat_area_m2=At)
    Vc = res_Vc.value
    print(f"Characteristic Length L*: {L_star:.2f} m")
    print(f"Total Chamber Volume Vc:  {Vc * 1e6:.1f} cm^3 ({Vc * 1e3:.3f} liters)")

    conv_half_angle = 25.0  # degrees
    res_Vconv = convergent_volume(
        throat_diameter_m=Dt,
        chamber_diameter_m=Dc,
        half_angle_deg=conv_half_angle
    )
    Vconv = res_Vconv.value
    print(f"Convergent Cone Volume:   {Vconv * 1e6:.1f} cm^3 (cone half-angle: {conv_half_angle} deg)")

    # 8. Cylindrical Section Length (Lc)
    res_Lc = cylinder_length(
        chamber_volume_m3=Vc,
        convergent_volume_m3=Vconv,
        chamber_area_m2=Ac
    )
    Lc = res_Lc.value
    print(f"Cylindrical Length (Lc):  {Lc * 1e3:.2f} mm ({Lc * 100 / Dt:.1f}% of Dt)")

    # 9. Chamber Residence / Stay Time
    # Mean chamber gas density: rho ~ Pc / (R * Tc)
    R_spec = 8314.46 / (M_w * 1000.0)
    rho_c = P_c / (R_spec * T_c)
    res_tau = chamber_bulk_residence_time(chamber_volume_m3=Vc, density_kg_m3=rho_c, mdot_kg_s=m_dot)
    print(f"\nMean Chamber Density:     {rho_c:.2f} kg/m^3")
    print(f"Bulk Residence Time:      {res_tau.value * 1e3:.2f} milliseconds")
    print("=" * 65)

if __name__ == "__main__":
    main()
