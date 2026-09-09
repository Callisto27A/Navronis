"""
Example 02: Bartz Hot-Gas Heat Transfer & Wall Conduction
==========================================================
Demonstrates the canonical Bartz equation (1957) and Fourier 1D thermal
conduction for a high-pressure rocket combustion chamber throat:
- Evaluates throat convective film coefficient (h_g)
- Determines boundary layer recovery temperature (T_aw)
- Solves for peak throat heat flux (q_throat)
- Calculates temperature drop (Delta T) across a 1.2 mm GRCop-42 copper liner
"""

from kryptonis.propulsion_equations.bartz import (
    bartz_film_coefficient,
    recovery_temperature,
)
from kryptonis.propulsion_equations.wall_conduction import fourier_wall_temperature_drop

def main():
    print("=" * 65)
    print("Kryptonis: Bartz Hot-Gas Heat Transfer & Wall Conduction")
    print("=" * 65)

    # 1. Throat Operating Conditions
    Dt = 0.05584          # Throat diameter: 55.84 mm [m]
    Rc = Dt * 1.5         # Upstream curvature radius (1.5 * Dt)
    Pc = 7.0e6            # 70 bar chamber pressure [Pa]
    c_star = 1760.7       # Characteristic velocity [m/s]
    gamma = 1.2
    T0 = 3450.0           # Combustion stagnation temperature [K]
    mach_throat = 1.0     # Choked sonic throat
    Pr = 0.82             # Prandtl number
    mu = 8.8e-5           # Dynamic viscosity [Pa*s]
    Cp = 2600.0           # Specific heat [J/(kg*K)]
    T_wg = 850.0          # Target hot-gas wall temperature [K] (typical for GRCop-42)

    print(f"Throat Diameter (Dt):     {Dt * 1e3:.2f} mm")
    print(f"Chamber Pressure (Pc):    {Pc / 1e6:.1f} MPa ({Pc / 1e5:.0f} bar)")
    print(f"Combustion Temp (T0):     {T0:.1f} K")
    print(f"Hot-Gas Wall Temp (Twg):  {T_wg:.1f} K")

    # 2. Compute Recovery Temperature T_aw
    res_Taw = recovery_temperature(
        stagnation_temperature_K=T0,
        mach=mach_throat,
        gamma=gamma,
        prandtl=Pr
    )
    T_aw = res_Taw.value
    print(f"\nAdiabatic Wall Temp (Taw): {T_aw:.1f} K (Recovery factor Pr^0.33: {Pr**0.33:.3f})")

    # 3. Compute Bartz Gas-Side Film Coefficient h_g
    res_hg = bartz_film_coefficient(
        throat_diameter_m=Dt,
        chamber_pressure_Pa=Pc,
        c_star_m_s=c_star,
        throat_radius_curvature_m=Rc,
        area_ratio_local_over_throat=1.0,
        gamma=gamma,
        mach=mach_throat,
        specific_heat_J_kgK=Cp,
        prandtl=Pr,
        viscosity_Pa_s=mu,
        stagnation_temperature_K=T0,
        wall_temperature_K=T_wg
    )
    h_g = res_hg.value
    print(f"Bartz Coefficient (hg):   {h_g:.1f} W/(m^2*K)")

    # 4. Peak Throat Heat Flux: q = h_g * (T_aw - T_wg)
    q_throat = h_g * (T_aw - T_wg)
    print(f"Peak Throat Heat Flux:    {q_throat / 1e6:.2f} MW/m^2")

    # 5. Wall Conduction across GRCop-42 Liner
    t_wall = 0.0012       # 1.2 mm liner thickness
    k_grcop42 = 330.0     # High-conductivity Cu-Cr-Nb alloy [W/(m*K)]

    res_dt = fourier_wall_temperature_drop(
        heat_flux_W_m2=q_throat,
        wall_thickness_m=t_wall,
        thermal_conductivity_W_m_K=k_grcop42
    )
    delta_T_wall = res_dt.value
    T_wc = T_wg - delta_T_wall

    print(f"\n--- Liner Conduction (GRCop-42) ---")
    print(f"Liner Thickness:          {t_wall * 1e3:.2f} mm")
    print(f"Thermal Conductivity (k): {k_grcop42:.1f} W/(m*K)")
    print(f"Wall Temperature Drop:    {delta_T_wall:.1f} K")
    print(f"Coolant-Side Wall Temp:   {T_wc:.1f} K")
    print("=" * 65)

if __name__ == "__main__":
    main()
