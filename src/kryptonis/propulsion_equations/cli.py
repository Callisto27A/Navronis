"""
Command-line interface (CLI) for Kryptonis Propulsion Equations.
================================================================
Provides immediate, zero-code rocket engine analytical sizing directly from the terminal.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

from kryptonis.propulsion_equations.chamber import (
    throat_area,
    throat_diameter,
    chamber_diameter,
    chamber_volume,
    convergent_volume,
    cylinder_length,
    vandenkerckhove,
    c_star_ideal,
    contraction_ratio,
)
from kryptonis.propulsion_equations.combustion import chamber_bulk_residence_time
from kryptonis.propulsion_equations.chamber_acoustics import (
    first_tangential_frequency,
    first_radial_frequency,
    first_longitudinal_frequency,
)


def run_chamber_sizing(
    thrust_n: float,
    pc_bar: float,
    propellants: str = "LOX/RP-1",
    cr_supplied: float | None = None,
    l_star: float | None = None,
    c_star_supplied: float | None = None,
    cf_est: float = 1.75,
    conv_half_angle_deg: float = 30.0,
    yield_strength_mpa: float = 280.0,
    safety_factor: float = 1.5,
    expansion_ratio: float = 20.0,
    nozzle_type: str = "bell",
    do_plot: bool = False,
    plot_save: str | None = None,
    export_json_path: str | None = None,
    export_csv_path: str | None = None,
    export_cadquery_path: str | None = None,
) -> int:
    pc_pa = pc_bar * 1.0e5

    prop_defaults = {
        "LOX/RP-1": {"gamma": 1.22, "mw": 0.0235, "tc": 3600.0, "l_star": 1.05, "c_star": 1780.0},
        "LOX/CH4": {"gamma": 1.20, "mw": 0.0220, "tc": 3450.0, "l_star": 1.00, "c_star": 1820.0},
        "LOX/LH2": {"gamma": 1.23, "mw": 0.0150, "tc": 3250.0, "l_star": 0.85, "c_star": 2350.0},
    }

    norm_prop = propellants.upper().replace("METHANE", "CH4").replace("KEROSENE", "RP-1")
    defaults = prop_defaults.get(norm_prop, prop_defaults["LOX/RP-1"])

    gamma = defaults["gamma"]
    mw = defaults["mw"]
    tc = defaults["tc"]
    effective_l_star = l_star if l_star is not None else defaults["l_star"]

    gamma_res = vandenkerckhove(gamma)
    if c_star_supplied is not None:
        c_star = c_star_supplied
    else:
        c_star_res = c_star_ideal(gamma=gamma, molar_mass_kg_per_mol=mw, chamber_temperature_K=tc)
        c_star = c_star_res.value

    at_est = thrust_n / (pc_pa * cf_est)
    mdot = (pc_pa * at_est) / c_star

    at_res = throat_area(mass_flow_kg_s=mdot, c_star_m_s=c_star, chamber_pressure_Pa=pc_pa)
    at_m2 = at_res.value
    dt_res = throat_diameter(throat_area_m2=at_m2)
    dt_m = dt_res.value

    if cr_supplied is not None:
        cr = cr_supplied
    else:
        cr_res = contraction_ratio(throat_diameter_m=dt_m)
        cr = cr_res.value

    dc_res = chamber_diameter(throat_diameter_m=dt_m, contraction_ratio_=cr)
    dc_m = dc_res.value
    ac_m2 = at_m2 * cr

    vc_res = chamber_volume(l_star_m=effective_l_star, throat_area_m2=at_m2)
    vc_m3 = vc_res.value

    vconv_res = convergent_volume(
        throat_diameter_m=dt_m,
        chamber_diameter_m=dc_m,
        half_angle_deg=conv_half_angle_deg,
    )
    vconv_m3 = vconv_res.value

    lconv_m = (dc_m - dt_m) / (2.0 * math.tan(math.radians(conv_half_angle_deg)))

    lcyl_res = cylinder_length(
        chamber_volume_m3=vc_m3,
        convergent_volume_m3=vconv_m3,
        chamber_area_m2=ac_m2,
    )
    lcyl_m = lcyl_res.value if not math.isnan(lcyl_res.value) else 0.0
    ltot_m = lcyl_m + lconv_m

    r_spec = 8314.4626 / (mw * 1000.0)
    rho_c = pc_pa / (r_spec * tc)
    tau_res = chamber_bulk_residence_time(chamber_volume_m3=vc_m3, density_kg_m3=rho_c, mdot_kg_s=mdot)

    a_sound = math.sqrt(gamma * r_spec * tc)
    f_1t = first_tangential_frequency(speed_of_sound_m_s=a_sound, chamber_diameter_m=dc_m).value
    f_1r = first_radial_frequency(speed_of_sound_m_s=a_sound, chamber_diameter_m=dc_m).value
    f_1l = first_longitudinal_frequency(speed_of_sound_m_s=a_sound, chamber_length_m=ltot_m).value

    allowable_stress_pa = (yield_strength_mpa * 1.0e6) / safety_factor
    chamber_radius_m = dc_m / 2.0
    t_wall_min_m = (pc_pa * chamber_radius_m) / allowable_stress_pa
    t_wall_rec_m = t_wall_min_m * 1.25

    print("=" * 78)
    print("           KRYPTONIS PROPULSION ENGINE SIZER -- COMPONENT 1: COMBUSTOR")
    print("=" * 78)
    print("Inputs:")
    print(f"  Thrust:                {thrust_n / 1e3:.2f} kN ({thrust_n:.0f} N)")
    print(f"  Chamber Pressure (Pc): {pc_bar:.2f} bar ({pc_pa / 1e6:.2f} MPa)")
    print(f"  Propellants:           {norm_prop} (gamma={gamma:.2f}, Tc={tc:.0f} K, Mw={mw*1e3:.1f} g/mol)")
    print(f"  Characteristic L*:     {effective_l_star:.2f} m")
    print(f"  Estimated Mass Flow:   {mdot:.2f} kg/s (assuming Cf={cf_est:.2f}, c*={c_star:.1f} m/s)")
    print("-" * 78)
    print("1. THROAT & CHAMBER SIZING (NASA SP-125 / Huzel & Huang)")
    print("-" * 78)
    print(f"  Throat Area (At):         {at_m2 * 1e4:8.2f} cm^2 ({at_m2:.6e} m^2)")
    print(f"  Throat Diameter (Dt):     {dt_m * 1e3:8.2f} mm")
    print(f"  Contraction Ratio (eps_c):{cr:8.2f}")
    print(f"  Chamber Area (Ac):        {ac_m2 * 1e4:8.2f} cm^2")
    print(f"  Chamber Diameter (Dc):    {dc_m * 1e3:8.2f} mm")
    print(f"  Chamber Volume (Vc):      {vc_m3 * 1e3:8.3f} liters ({vc_m3 * 1e6:.1f} cm^3)")
    print("-" * 78)
    print("2. AXIAL PROFILE & STAY TIME")
    print("-" * 78)
    print(f"  Convergent Half-Angle:    {conv_half_angle_deg:8.1f} deg")
    print(f"  Convergent Length (Lconv):{lconv_m * 1e3:8.2f} mm")
    print(f"  Cylindrical Barrel (Lcyl):{lcyl_m * 1e3:8.2f} mm")
    print(f"  Total Chamber Length (Lc):{ltot_m * 1e3:8.2f} mm")
    print(f"  Mean Stay / Res. Time:    {tau_res.value * 1e3:8.2f} ms")
    print("-" * 78)
    print("3. ACOUSTIC STABILITY FREQUENCIES (NASA SP-194)")
    print("-" * 78)
    print(f"  Chamber Speed of Sound:   {a_sound:8.1f} m/s")
    print(f"  1st Tangential Mode (1T): {f_1t:8.1f} Hz  [Primary transverse buzz]")
    print(f"  1st Radial Mode (1R):     {f_1r:8.1f} Hz")
    print(f"  1st Longitudinal (1L):    {f_1l:8.1f} Hz  [Chugging / organ-pipe coupling]")
    print("-" * 78)
    print("4. MECHANICAL INTEGRITY (ASME Sec VIII / Thin-Shell Hoop)")
    print("-" * 78)
    print(f"  Material Yield Strength:  {yield_strength_mpa:8.1f} MPa (SF: {safety_factor:.1f})")
    print(f"  Allowable Stress:         {allowable_stress_pa / 1e6:8.2f} MPa")
    print(f"  Min Wall Thickness:       {t_wall_min_m * 1e3:8.2f} mm")
    print(f"  Recommended Thickness:    {t_wall_rec_m * 1e3:8.2f} mm (+25% machining margin)")
    print("=" * 78)

    # --- Profile generation, plotting, and export ---
    needs_profile = do_plot or export_json_path or export_csv_path or export_cadquery_path

    if needs_profile:
        from kryptonis.propulsion_equations.profile import generate_chamber_profile
        profile = generate_chamber_profile(
            throat_diameter_m=dt_m,
            chamber_diameter_m=dc_m,
            cylindrical_length_m=lcyl_m,
            convergent_half_angle_deg=conv_half_angle_deg,
            expansion_ratio=expansion_ratio,
            nozzle_type=nozzle_type,
        )
        print(f"  Profile generated: {len(profile.x)} points, "
              f"total length {profile.total_length_m * 1e3:.1f} mm")

    if do_plot:
        from kryptonis.propulsion_equations.plotting import plot_chamber_profile
        title = (f"{norm_prop}  {thrust_n/1e3:.0f} kN  Pc={pc_bar:.0f} bar  "
                 f"eps={expansion_ratio:.0f}")
        plot_chamber_profile(
            profile,
            title=title,
            wall_thickness_mm=t_wall_rec_m * 1e3,
            save_path=plot_save,
            show=(plot_save is None),
        )
        if plot_save:
            print(f"  Plot saved: {os.path.abspath(plot_save)}")

    sizing_results = {
        "thrust_kN": thrust_n / 1e3,
        "chamber_pressure_bar": pc_bar,
        "propellants": norm_prop,
        "c_star_m_s": c_star,
        "mass_flow_kg_s": mdot,
        "throat_diameter_mm": dt_m * 1e3,
        "chamber_diameter_mm": dc_m * 1e3,
        "chamber_volume_liters": vc_m3 * 1e3,
        "stay_time_ms": tau_res.value * 1e3,
        "convergent_length_mm": lconv_m * 1e3,
        "cylindrical_length_mm": lcyl_m * 1e3,
        "acoustic_1T_Hz": f_1t,
        "acoustic_1R_Hz": f_1r,
        "acoustic_1L_Hz": f_1l,
        "min_wall_thickness_mm": t_wall_min_m * 1e3,
        "recommended_wall_thickness_mm": t_wall_rec_m * 1e3,
    }

    if export_json_path:
        from kryptonis.propulsion_equations.export import export_json
        p = export_json(profile, export_json_path, sizing_results=sizing_results)
        print(f"  JSON exported: {p}")

    if export_csv_path:
        from kryptonis.propulsion_equations.export import export_csv
        p = export_csv(profile, export_csv_path)
        print(f"  CSV exported:  {p}")

    if export_cadquery_path:
        from kryptonis.propulsion_equations.export import export_cadquery_script
        p = export_cadquery_script(profile, export_cadquery_path)
        print(f"  CadQuery script exported: {p}")

    print("=" * 78)
    print("Execution complete. 100% first-principles closed-form analytical solutions.")
    print("=" * 78)

    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="kryptonis-chamber",
        description="Authority-controlled analytical sizing for liquid rocket engine thrust chambers.",
    )
    parser.add_argument(
        "--thrust", "-t", type=float, default=30000.0,
        help="Thrust in Newtons (default: 30000 N = 30 kN)",
    )
    parser.add_argument(
        "--pc", "-p", type=float, default=70.0,
        help="Chamber pressure in bar (default: 70.0 bar = 7.0 MPa)",
    )
    parser.add_argument(
        "--propellants", default="LOX/RP-1",
        choices=["LOX/RP-1", "LOX/CH4", "LOX/LH2"],
        help="Propellant combination (default: LOX/RP-1)",
    )
    parser.add_argument(
        "--cr", type=float, default=None,
        help="Contraction ratio (Ac/At). If omitted, evaluates Humble correlation.",
    )
    parser.add_argument(
        "--lstar", type=float, default=None,
        help="Characteristic chamber length L* in meters (default: propellant standard)",
    )
    parser.add_argument(
        "--expansion-ratio", type=float, default=20.0,
        help="Nozzle expansion ratio A_e/A_t (default: 20)",
    )
    parser.add_argument(
        "--nozzle", default="bell", choices=["bell", "conical"],
        help="Nozzle contour type (default: bell)",
    )
    parser.add_argument(
        "--yield-strength", type=float, default=280.0,
        help="Liner material yield strength in MPa (default: 280.0 MPa)",
    )
    parser.add_argument(
        "--safety-factor", type=float, default=1.5,
        help="Structural safety factor (default: 1.5)",
    )

    # Visualization & export
    parser.add_argument(
        "--plot", action="store_true",
        help="Display a matplotlib cross-section plot of the engine",
    )
    parser.add_argument(
        "--plot-save", type=str, default=None, metavar="FILE",
        help="Save the plot to a file (PNG, PDF, SVG) instead of displaying",
    )
    parser.add_argument(
        "--export-json", type=str, default=None, metavar="FILE",
        help="Export full geometry + sizing to JSON (for CadQuery / FreeCAD)",
    )
    parser.add_argument(
        "--export-csv", type=str, default=None, metavar="FILE",
        help="Export (x, r) profile points to CSV",
    )
    parser.add_argument(
        "--export-cadquery", type=str, default=None, metavar="FILE",
        help="Export a ready-to-run CadQuery .py script that generates STEP",
    )

    args = parser.parse_args()

    sys.exit(
        run_chamber_sizing(
            thrust_n=args.thrust,
            pc_bar=args.pc,
            propellants=args.propellants,
            cr_supplied=args.cr,
            l_star=args.lstar,
            yield_strength_mpa=args.yield_strength,
            safety_factor=args.safety_factor,
            expansion_ratio=args.expansion_ratio,
            nozzle_type=args.nozzle,
            do_plot=args.plot or (args.plot_save is not None),
            plot_save=args.plot_save,
            export_json_path=args.export_json,
            export_csv_path=args.export_csv,
            export_cadquery_path=args.export_cadquery,
        )
    )


if __name__ == "__main__":
    main()
