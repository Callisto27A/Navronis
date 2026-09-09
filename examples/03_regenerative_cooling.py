"""
Example 03: Regenerative Cooling Channel Fluid Dynamics
========================================================
Calculates the coolant-side convective heat transfer coefficient (h_c),
Nusselt number, Ito curvature enhancement factor, and hydraulic pressure drop
for supercritical methane flowing through milled regenerative channels.
"""

from kryptonis.propulsion_equations.thermal import (
    nusselt_sieder_tate,
    nusselt_gnielinski,
    ito_curvature_factor,
    friction_factor_haaland,
)

def main():
    print("=" * 65)
    print("Kryptonis: Regenerative Cooling Channel Fluid Dynamics")
    print("=" * 65)

    # 1. Channel Geometry at Nozzle Throat
    width_m = 0.0015       # 1.5 mm rib-to-rib channel width
    height_m = 0.0035      # 3.5 mm channel height (aspect ratio ~ 2.33)
    A_flow = width_m * height_m
    perimeter = 2.0 * (width_m + height_m)
    D_h = 4.0 * A_flow / perimeter   # Hydraulic diameter ~ 2.1 mm
    r_hyd = D_h / 2.0
    R_curve = 0.083        # Throat concave contour radius of curvature (~83 mm)

    print(f"Channel Width:            {width_m * 1e3:.2f} mm")
    print(f"Channel Height:           {height_m * 1e3:.2f} mm")
    print(f"Hydraulic Diameter (Dh):  {D_h * 1e3:.2f} mm")

    # 2. Supercritical Methane Coolant Properties at 12 MPa (120 bar)
    # Bulk coolant temperature: 150 K
    rho_coolant = 380.0    # kg/m^3 (liquid-dense supercritical CH4)
    v_coolant = 25.0       # 25 m/s flow velocity
    mu_bulk = 8.5e-5       # Pa*s
    mu_wall = 6.0e-5       # Pa*s (hotter near wall -> lower viscosity)
    k_coolant = 0.160      # W/(m*K)
    Cp_coolant = 3600.0    # J/(kg*K)

    # Dimensionless numbers
    Re = (rho_coolant * v_coolant * D_h) / mu_bulk
    Pr = (Cp_coolant * mu_bulk) / k_coolant

    print(f"\nCoolant Flow Velocity:    {v_coolant:.1f} m/s")
    print(f"Reynolds Number (Re):     {Re:.2e} (Fully Turbulent)")
    print(f"Prandtl Number (Pr):      {Pr:.2f}")

    # 3. Nusselt Number (Sieder-Tate with viscosity ratio)
    res_st = nusselt_sieder_tate(
        reynolds=Re,
        prandtl=Pr,
        bulk_viscosity_Pa_s=mu_bulk,
        wall_viscosity_Pa_s=mu_wall
    )
    Nu_st = res_st.value
    print(f"\nSieder-Tate Nusselt (Nu): {Nu_st:.1f}")

    # 4. Ito Throat Curvature Enhancement Factor (concave hot gas side)
    res_ito = ito_curvature_factor(
        reynolds_bulk=Re,
        hydraulic_radius_m=r_hyd,
        radius_of_curvature_m=R_curve,
        side="concave"
    )
    ito_factor = res_ito.value
    Nu_enhanced = Nu_st * ito_factor

    print(f"Ito Curvature Factor:     {ito_factor:.3f} (+{(ito_factor - 1.0) * 100:.1f}% enhancement)")
    print(f"Enhanced Throat Nusselt:  {Nu_enhanced:.1f}")

    # 5. Coolant Convective Film Coefficient: h_c = Nu * k / Dh
    h_c = Nu_enhanced * k_coolant / D_h
    print(f"Coolant Film Coeff (hc):  {h_c:.1f} W/(m^2*K)")

    # 6. Darcy Friction Factor & Pressure Drop per 10 cm Channel Length
    rel_rough = 1.5e-6 / D_h  # 1.5 um typical milled surface roughness
    res_f = friction_factor_haaland(reynolds=Re, relative_roughness=rel_rough)
    f_darcy = res_f.value
    print(f"\nHaaland Friction Factor:  {f_darcy:.4f}")

    L_channel = 0.10  # 10 cm section around throat
    delta_P = f_darcy * (L_channel / D_h) * (0.5 * rho_coolant * (v_coolant ** 2))
    print(f"Frictional Delta P (10cm):{delta_P / 1e5:.2f} bar ({delta_P / 1e3:.1f} kPa)")
    print("=" * 65)

if __name__ == "__main__":
    main()
