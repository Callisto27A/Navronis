"""
Example 06: Regenerative Cooling Channel Sizing & 1D Conjugate Thermal Marching (Day 3)
========================================================================================
Demonstrates the full end-to-end design of a 30 kN Methalox rocket thrust chamber cooling jacket:
1. Sizing milled rectangular cooling channels at the sonic throat.
2. Evaluating coolant turbulent convection via Gnielinski (1976) and Haaland (1983).
3. Computing fin efficiency enhancement across rib walls.
4. Solving the 3-resistance conjugate wall energy balance for Twg and Twc without guessing.
5. Evaluating Darcy-Weisbach coolant pressure drop, temperature rise, and thermal stress.
"""

from kryptonis.propulsion_equations import (
    CombustorDesign,
    RegenCoolingJacket,
)


def main():
    print("=" * 80)
    print("NAVRONIS PROPULSION: DAY 3 REGENERATIVE COOLING CHANNEL DESIGN BENCHMARK")
    print("=" * 80)

    # 1. Size the combustor first (Day 1)
    engine = CombustorDesign(
        thrust=30000.0,            # 30 kN
        chamber_pressure=12.0e6,    # 120 bar
        mixture_ratio=2.6,         # LOX / LCH4
        propellant="LOX/CH4",
        convergent_half_angle_deg=30.0,
    )
    c_res = engine.solve()

    mdot_fuel = c_res.total_mass_flow / (1.0 + engine.mixture_ratio)
    print(f"Throat Diameter (Dt):        {c_res.throat_diameter * 1000:.2f} mm")
    print(f"Combustor Diameter (Dc):     {c_res.chamber_diameter * 1000:.2f} mm")
    print(f"Total Chamber Length (Lc):   {c_res.chamber_length * 1000:.2f} mm")
    print(f"Fuel Mass Flow (Coolant):    {mdot_fuel:.3f} kg/s (LCH4)")
    print("-" * 80)

    # 2. Size the regenerative cooling jacket (Day 3)
    jacket = RegenCoolingJacket(
        throat_diameter_m=c_res.throat_diameter,
        chamber_diameter_m=c_res.chamber_diameter,
        chamber_length_m=c_res.chamber_length,
        mass_flow_coolant_kg_s=mdot_fuel,
        chamber_pressure_pa=12.0e6,
        gas_recovery_temp_k=3400.0,
        gas_throat_htc_w_m2k=c_res.heat_flux_screen / 2500.0, # Approximate throat gas HTC
        n_channels=80,
        channel_height_m=0.0018,    # 1.8 mm channel depth
        fin_thickness_m=0.0008,     # 0.8 mm rib fin width
        wall_thickness_m=0.0015,    # 1.5 mm liner thickness
        coolant_type="CH4",
        liner_material="CuCrZr",
    )
    j_res = jacket.solve()

    # 3. Report Results
    print("REGENERATIVE COOLING JACKET SPECIFICATIONS:")
    print(f"  Milled Channel Count (N):    {j_res.n_channels} passages")
    print(f"  Throat Channel Width (wc):   {j_res.channel_width_mm:.3f} mm")
    print(f"  Channel Height (hc):         {j_res.channel_height_mm:.3f} mm")
    print(f"  Channel Aspect Ratio (AR):   {j_res.aspect_ratio:.2f}")
    print(f"  Hydraulic Diameter (Dh):     {j_res.hydraulic_diameter_mm:.3f} mm")
    print("-" * 80)
    print("COOLANT HYDRAULICS & HEAT TRANSFER:")
    print(f"  Coolant Velocity (vc):       {j_res.coolant_velocity_m_s:.1f} m/s")
    print(f"  Reynolds Number (Re):        {j_res.reynolds_number:.0f} (Fully Turbulent)")
    print(f"  Haaland Friction Factor (f): {j_res.friction_factor:.4f}")
    print(f"  Coolant Base HTC (hc):       {j_res.coolant_htc_W_m2K:.1f} W/m²-K (Gnielinski 1976)")
    print(f"  Fin Efficiency (eta_fin):    {j_res.fin_efficiency * 100:.1f}%")
    print(f"  Enhanced Effective HTC:      {j_res.enhanced_coolant_htc_W_m2K:.1f} W/m²-K")
    print("-" * 80)
    print("CONJUGATE THERMAL EQUILIBRIUM (NO GUESSES):")
    print(f"  Hot-Gas Wall Temp (Twg):     {j_res.hot_gas_wall_temp_K:.1f} K ({j_res.hot_gas_wall_temp_K - 273.15:.1f} °C)")
    print(f"  Coolant Wall Temp (Twc):     {j_res.coolant_wall_temp_K:.1f} K ({j_res.coolant_wall_temp_K - 273.15:.1f} °C)")
    print(f"  Peak Throat Heat Flux (q):   {j_res.peak_heat_flux_MW_m2:.2f} MW/m²")
    print(f"  Coolant Pressure Drop (Delta P):  {j_res.coolant_pressure_drop_bar:.2f} bar")
    print(f"  Coolant Bulk Temp Rise (Delta T): {j_res.coolant_temp_rise_K:.1f} K")
    print(f"  Thermal Compressive Stress:  {j_res.thermal_stress_MPa:.1f} MPa")
    print(f"  Yield Safety Margin (MS):    {j_res.yield_safety_margin:.2f} (PASS >= 0)")
    print(f"  Boiling/Coking Margin:       {'PASS' if j_res.boiling_margin_adequate else 'FAIL'}")
    print("=" * 80)

    # 4. Provenance query
    print("\nQUERYING ENGINEERING CITATION:")
    print(j_res.explain("hot_gas_wall_temp_K"))


if __name__ == "__main__":
    main()
