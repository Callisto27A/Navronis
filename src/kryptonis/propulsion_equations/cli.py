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


def run_injector_sizing(
    thrust_n: float,
    pc_bar: float,
    propellants: str = "LOX/CH4",
    injector_type: str = "coaxial",
    n_elements: int = 19,
    delta_p_ratio: float = 0.20,
) -> int:
    from kryptonis.propulsion_equations.injector import InjectorDesign

    pc_pa = pc_bar * 1e5
    # Standard propellant properties (density kg/m³, typical Isp sea-level s, nominal O/F)
    prop_table = {
        "LOX/CH4": {"rho_ox": 1141.0, "rho_f": 422.0, "isp": 295.0, "of": 3.5},
        "LOX/RP-1": {"rho_ox": 1141.0, "rho_f": 810.0, "isp": 285.0, "of": 2.6},
        "LOX/LH2": {"rho_ox": 1141.0, "rho_f": 71.0, "isp": 390.0, "of": 6.0},
    }
    norm = propellants.upper().replace("METHANE", "CH4").replace("KEROSENE", "RP-1")
    pinfo = prop_table.get(norm, prop_table["LOX/CH4"])

    # Mass flow estimate: m_dot = Thrust / (Isp * 9.80665)
    m_dot_total = thrust_n / (pinfo["isp"] * 9.80665)
    of = pinfo["of"]
    m_dot_fuel = m_dot_total / (1.0 + of)
    m_dot_ox = m_dot_total - m_dot_fuel

    des = InjectorDesign(
        injector_type=injector_type,
        chamber_pressure=pc_pa,
        mass_flow_ox=m_dot_ox,
        mass_flow_fuel=m_dot_fuel,
        rho_ox=pinfo["rho_ox"],
        rho_fuel=pinfo["rho_f"],
        delta_p_ratio=delta_p_ratio,
        n_elements=n_elements,
    )
    res = des.solve()

    print("=" * 78)
    print(f"NAVRONIS PROPULSION -- INJECTOR SIZING REPORT: {injector_type.upper()}")
    print(f"Propellant: {propellants} | Thrust: {thrust_n/1e3:.1f} kN | Pc: {pc_bar:.1f} bar")
    print("=" * 78)
    print(f"Mass Flow: Total = {m_dot_total:.3f} kg/s (LOX: {m_dot_ox:.3f} kg/s, Fuel: {m_dot_fuel:.3f} kg/s)")
    print(f"Injector Delta P: {res['delta_p_bar']:.2f} bar ({res['delta_p_ratio']*100:.1f}% Pc)")
    print(f"Chugging Decoupling Margin: {'PASS (>=15%)' if res['chugging_margin_adequate'] else 'FAIL (<15%)'}")
    print("-" * 78)

    if injector_type in {"coaxial", "shear_coaxial"}:
        print(f"Elements:                  {res['n_elements']}")
        print(f"Liquid Post ID:            {res['post_id_mm']:.2f} mm")
        print(f"Liquid Post OD:            {res['post_od_mm']:.2f} mm")
        print(f"Gas Sleeve ID:             {res['annulus_id_mm']:.2f} mm")
        print(f"Annular Gap:               {res['annular_gap_mm']:.2f} mm")
        print(f"Liquid Ox Velocity:        {res['v_ox_m_s']:.2f} m/s")
        print(f"Gas/Fuel Velocity:         {res['v_fuel_m_s']:.2f} m/s")
        print(f"Momentum Flux Ratio J:     {res['momentum_flux_ratio_J']:.2f}  [Target: 2.0 - 20.0]")
        print(f"Velocity Ratio VR:         {res['velocity_ratio_VR']:.2f}")
        print(f"Recess Length:             {res['recess_length_mm']:.2f} mm")
        print(f"Droplet SMD (D32):         {res['smd_um']:.1f} um")
    elif injector_type in {"swirl", "swirl_coaxial"}:
        print(f"Elements:                  {res['n_elements']}")
        print(f"Centrifugal Orifice Diam:  {res['orifice_diameter_mm']:.2f} mm")
        print(f"Central Gas Core Diam:     {res['gas_core_diameter_mm']:.2f} mm")
        print(f"Liquid Film Thickness:     {res['liquid_film_thickness_mm']:.3f} mm")
        print(f"Spray Cone Half-Angle:     {res['spray_half_angle_deg']:.1f} deg")
        print(f"Tangential Inlet Diam:     {res['tangential_inlet_diameter_mm']:.2f} mm")
        print(f"Coaxial Gas Annulus Gap:   {res['annular_gas_gap_mm']:.2f} mm")
        print(f"Gas Velocity:              {res['gas_velocity_m_s']:.2f} m/s")
        print(f"Geometric Swirl K:         {res['geometric_swirl_K']:.2f}")
        print(f"Droplet SMD (D32):         {res['smd_um']:.1f} um")
    elif injector_type in {"pintle", "pintle_injector"}:
        print(f"Pintle Diameter:           {res['pintle_diameter_mm']:.1f} mm")
        print(f"Annular Gap Thickness:     {res['annular_gap_thickness_mm']:.3f} mm")
        print(f"Annular Fuel Velocity:     {res['annular_velocity_m_s']:.2f} m/s")
        print(f"Radial Slot Height:        {res['radial_slot_height_mm']:.3f} mm")
        print(f"Radial Ox Velocity:        {res['radial_velocity_m_s']:.2f} m/s")
        print(f"Total Momentum Ratio TMR:  {res['total_momentum_ratio_TMR']:.3f}")
        print(f"Spray Cone Half-Angle:     {res['spray_half_angle_deg']:.1f} deg")
    else:
        print(f"Elements:                  {res['n_elements']}")
        print(f"Oxidizer Orifice Diam:     {res['orifice_diameter_1_mm']:.2f} mm")
        print(f"Fuel Orifice Diam:         {res['orifice_diameter_2_mm']:.2f} mm")
        print(f"Oxidizer Jet Velocity:     {res['jet_velocity_1_m_s']:.2f} m/s")
        print(f"Fuel Jet Velocity:         {res['jet_velocity_2_m_s']:.2f} m/s")
        print(f"Rupe Momentum Parameter:   {res['rupe_momentum_parameter']:.2f}")
        print(f"Free Jet Length:           {res['free_jet_length_mm']:.2f} mm")
        print(f"Droplet SMD (D32):         {res['smd_um']:.1f} um")

    print("=" * 78)
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="navronis",
        description="Authority-controlled analytical sizing for liquid rocket engine thrust chambers and injectors.",
    )
    parser.add_argument(
        "--subsystem", default="chamber", choices=["chamber", "injector"],
        help="Subsystem to size: 'chamber' or 'injector'",
    )
    parser.add_argument(
        "--injector-type", default="coaxial", choices=["coaxial", "swirl", "pintle", "impinging"],
        help="4 Canonical injector families: 'coaxial', 'swirl', 'pintle', or 'impinging' (for --subsystem injector)",
    )
    parser.add_argument(
        "--elements", type=int, default=19,
        help="Number of injector elements (default: 19)",
    )
    parser.add_argument(
        "--delta-p-ratio", type=float, default=0.20,
        help="Injector pressure drop ratio delta_p / Pc (default: 0.20 = 20%)",
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

    if args.subsystem == "injector":
        sys.exit(
            run_injector_sizing(
                thrust_n=args.thrust,
                pc_bar=args.pc,
                propellants=args.propellants,
                injector_type=args.injector_type,
                n_elements=args.elements,
                delta_p_ratio=args.delta_p_ratio,
            )
        )

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
