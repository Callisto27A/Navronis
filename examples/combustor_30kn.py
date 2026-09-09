"""
Kryptonis Propulsion: Day 1 Combustor Release Example
=====================================================
Sizes a 30 kN LOX/LCH4 combustion chamber at Pc = 120 bar,
explains the provenance of every major dimension, and exports
JSON, CSV, and a cross-section plot.
"""

from kryptonis.propulsion_equations import CombustorDesign

def main():
    print("=" * 75)
    print("  KRYPTONIS PROPULSION: DAY 1 -- COMBUSTOR PRELIMINARY DESIGN")
    print("=" * 75)

    # 1. Instantiate the user-facing design specification
    design = CombustorDesign(
        thrust=30000.0,            # 30 kN
        chamber_pressure=12.0e6,    # 120 bar (12.0 MPa)
        mixture_ratio=2.6,         # LOX / LCH4 stoichiometric-ish ratio
        propellant="LOX/CH4",
        convergent_half_angle_deg=30.0,
    )

    # 2. Solve the closed-form analytical equations
    result = design.solve()

    # 3. Inspect the calculated dimensions
    print("\n[1] MASS FLOW & PROPELLANT ALLOCATION")
    print(f"  Total Mass Flow:          {result.total_mass_flow:.2f} kg/s")
    print(f"  Fuel Mass Flow (CH4):     {result.fuel_mass_flow:.2f} kg/s")
    print(f"  Oxidizer Mass Flow (LOX): {result.oxidizer_mass_flow:.2f} kg/s")

    print("\n[2] GEOMETRY (NASA SP-125 / HUZEL & HUANG)")
    print(f"  Throat Diameter:          {result.throat_diameter * 1e3:.2f} mm")
    print(f"  Throat Area:              {result.throat_area * 1e4:.2f} cm^2")
    print(f"  Contraction Ratio (eps_c):{result.contraction_ratio:.2f}")
    print(f"  Chamber Diameter:         {result.chamber_diameter * 1e3:.2f} mm")
    print(f"  Chamber Volume:           {result.chamber_volume * 1e3:.3f} liters")
    print(f"  Cylindrical Length:       {result.cylindrical_length * 1e3:.2f} mm")
    print(f"  Total Chamber Length:     {result.chamber_length * 1e3:.2f} mm")
    print(f"  Bulk Stay Time:           {result.stay_time * 1e3:.2f} ms")

    print("\n[3] PRELIMINARY SCREENING (THERMAL / ACOUSTIC / STRUCTURAL)")
    print(f"  Peak Throat Heat Flux:    {result.heat_flux_screen / 1e6:.2f} MW/m^2 (Bartz 1957)")
    print(f"  1st Tangential Mode (1T): {result.acoustic_modes['1T_Hz']:.1f} Hz (NASA SP-194)")
    print(f"  1st Radial Mode (1R):     {result.acoustic_modes['1R_Hz']:.1f} Hz")
    print(f"  1st Longitudinal (1L):    {result.acoustic_modes['1L_Hz']:.1f} Hz")
    print(f"  Min Wall Thickness:       {result.wall_thickness_screen * 1e3:.2f} mm (ASME Sec VIII)")
    print(f"  Recommended Wall:         {result.wall_thickness_recommended * 1e3:.2f} mm (+25% margin)")

    # 4. Inspect Source Provenance
    print("\n[4] PROVENANCE DEMONSTRATION")
    print("-" * 55)
    print(result.explain("throat_diameter"))
    print("-" * 55)

    # 5. Export results
    json_path = "combustor_30kn.json"
    csv_path = "combustor_30kn.csv"
    plot_path = "combustor_30kn.png"

    result.to_json(json_path)
    result.to_csv(csv_path)
    print(f"\n[SAVED] Structured JSON: {json_path}")
    print(f"[SAVED] Parameter CSV:   {csv_path}")

    try:
        result.plot(save_path=plot_path, show=False)
        print(f"[SAVED] 2D Contour Plot: {plot_path}")
    except Exception as e:
        print(f"[NOTE] Matplotlib plot skipped: {e}")

    print("\n" + "=" * 75)
    print("  DAY 1 COMBUSTOR DESIGN COMPLETE. READY FOR GITHUB.")
    print("=" * 75)

if __name__ == "__main__":
    main()
