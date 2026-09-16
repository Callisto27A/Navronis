"""
Kryptonis Propulsion: High-Level Nozzle Aerodynamics & Altitude Engine (Day 4 Release)
======================================================================================
Provides an authority-controlled, primary-literature-traceable preliminary supersonic
nozzle design and altitude performance analysis engine.

Features:
  - 1D isentropic Area-Mach supersonic inversion (NASA SP-125)
  - Rao (1958) 80% parabolic bell contour coordinates via quadratic Bézier spline
  - 2D divergence efficiency factor (lambda) for conical and Rao bell geometries
  - Turbulent boundary layer displacement thickness (delta*) and effective flow area
  - Stark (2005) tri-criteria flow separation assessment (Summerfield 1954, Schmucker 1973, Kalt-Badal 1965)
  - U.S. Standard Atmosphere (1976) barometric altitude trajectory solver (0 to 86 km)
  - Delivered thrust and specific impulse across altitude from sea level to vacuum
  - Self-explaining engineering provenance and JSON/CSV export capabilities
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from kryptonis.propulsion_equations.aerodynamics import (
    calculate_area_ratio,
    calculate_mach_from_area_ratio,
    calculate_thrust_coefficient,
)
from kryptonis.propulsion_equations.chamber import (
    get_thermochemical_preset,
    c_star_ideal,
)
from kryptonis.propulsion_equations.contour import (
    calculate_divergence_loss_factor,
    calculate_bell_divergence_loss_factor,
    calculate_displacement_thickness,
    calculate_effective_area,
    calculate_conical_nozzle_length,
    generate_rao_bell_coordinates,
)
from kryptonis.propulsion_equations.nozzle_losses import (
    divergence_efficiency,
    separation_assessment,
    boundary_layer_displacement,
    delivered_thrust_coefficient,
    SEA_LEVEL_PA,
    STARK_2005,
)
from kryptonis.propulsion_equations.units import Result, Status


def standard_atmosphere(altitude_m: float) -> Tuple[float, float, float]:
    """
    Computes ambient atmospheric static pressure, temperature, and density
    using the 1976 U.S. Standard Atmosphere model from sea level to 86 km.

    Parameters
    ----------
    altitude_m : float
        Geometric altitude above mean sea level in meters (0 <= z <= 86000 m).

    Returns
    -------
    (pressure_Pa, temperature_K, density_kg_m3) : Tuple[float, float, float]
        Atmospheric state at the specified altitude.
    """
    if altitude_m < 0.0:
        altitude_m = 0.0

    # Atmospheric constants (1976 Standard Atmosphere)
    P0 = 101325.0       # Sea level standard pressure [Pa]
    T0 = 288.15         # Sea level standard temperature [K]
    R_air = 287.05287   # Gas constant for dry air [J/(kg*K)]
    g0 = 9.80665        # Standard gravitational acceleration [m/s^2]
    R_earth = 6356766.0 # Earth radius for geopotential conversion [m]

    # Convert geometric altitude to geopotential altitude [m]
    h = (R_earth * altitude_m) / (R_earth + altitude_m)

    # Base heights, lapse rates, base temperatures, and base pressures for layers
    layers = [
        # (h_base [m], lapse_rate L [K/m], T_base [K], P_base [Pa])
        (0.0,     -0.0065, 288.15, 101325.0),
        (11000.0,  0.0,    216.65, 22632.1),
        (20000.0,  0.0010, 216.65, 5474.89),
        (32000.0,  0.0028, 228.65, 868.019),
        (47000.0,  0.0,    270.65, 110.906),
        (51000.0, -0.0028, 270.65, 66.9389),
        (71000.0, -0.0020, 214.65, 3.95642),
        (84852.0,  0.0,    186.95, 0.37338),
    ]

    current_layer = layers[0]
    for layer in layers:
        if h >= layer[0]:
            current_layer = layer
        else:
            break

    h_base, L, T_base, P_base = current_layer
    dh = min(h - h_base, 84852.0 - h_base)

    if abs(L) < 1.0e-9:
        # Isothermal layer
        T = T_base
        P = P_base * math.exp(-g0 * dh / (R_air * T_base))
    else:
        # Linear gradient layer
        T = T_base + L * dh
        P = P_base * math.pow(T / T_base, -g0 / (L * R_air))

    rho = P / (R_air * T) if T > 0 else 0.0
    return (P, T, rho)


@dataclass
class NozzleDesign:
    """
    Preliminary supersonic liquid rocket nozzle aerodynamic sizing and altitude performance specification.
    """

    expansion_ratio: float = 20.0
    chamber_pressure: float = 7.0e6
    ambient_pressure: Optional[float] = None
    altitude: float = 0.0
    throat_diameter: Optional[float] = None
    throat_area: Optional[float] = None
    thrust: Optional[float] = None
    c_star: Optional[float] = None
    propellant: str = "LOX/RP-1"
    gamma: Optional[float] = None
    chamber_temperature: Optional[float] = None
    molar_mass: Optional[float] = None
    nozzle_type: str = "bell"
    divergent_half_angle_deg: float = 15.0
    initial_wall_angle_deg: float = 30.0
    exit_wall_angle_deg: float = 8.0
    bell_fractional_length: float = 0.80

    def solve(self) -> "NozzleResult":
        """
        Executes the closed-form supersonic aerodynamics, contour sizing, loss assessment,
        and altitude performance marching sequence.
        """
        if self.expansion_ratio <= 1.0:
            raise ValueError(f"Expansion ratio must be > 1.0, got {self.expansion_ratio}")
        if self.chamber_pressure <= 0.0:
            raise ValueError(f"Chamber pressure must be positive, got {self.chamber_pressure}")

        # Thermochemical defaults
        preset = get_thermochemical_preset(self.propellant)
        gamma = self.gamma if self.gamma is not None else preset["gamma"]
        tc = self.chamber_temperature if self.chamber_temperature is not None else preset["tc"]
        mw = self.molar_mass if self.molar_mass is not None else preset["mw"]

        if self.c_star is not None:
            cstar_val = self.c_star
        else:
            cstar_res = c_star_ideal(gamma=gamma, molar_mass_kg_per_mol=mw, chamber_temperature_K=tc)
            cstar_val = cstar_res.value

        # Ambient pressure from altitude if not supplied
        if self.ambient_pressure is not None:
            pa = self.ambient_pressure
        else:
            pa, _, _ = standard_atmosphere(self.altitude)

        # 1. Throat Sizing
        if self.throat_area is not None:
            at = self.throat_area
            dt = math.sqrt(4.0 * at / math.pi)
        elif self.throat_diameter is not None:
            dt = self.throat_diameter
            at = 0.25 * math.pi * dt ** 2
        elif self.thrust is not None:
            cf_est = 1.65
            at = self.thrust / (self.chamber_pressure * cf_est)
            dt = math.sqrt(4.0 * at / math.pi)
        else:
            dt = 0.050
            at = 0.25 * math.pi * dt ** 2

        rt = dt / 2.0
        ae = at * self.expansion_ratio
        de = dt * math.sqrt(self.expansion_ratio)
        re = de / 2.0

        # 2. 1D Supersonic Area-Mach Inversion (NASA SP-125)
        exit_mach = calculate_mach_from_area_ratio(self.expansion_ratio, gamma, supersonic=True)

        # 3. Isentropic Static State at Exit Plane
        t_total_ratio = 1.0 + 0.5 * (gamma - 1.0) * (exit_mach ** 2)
        p_ratio = math.pow(t_total_ratio, -gamma / (gamma - 1.0))
        pe = self.chamber_pressure * p_ratio
        te = tc / t_total_ratio

        if mw < 0.100:  # already in kg/mol
            R_spec = 8.314462618 / mw
        else:
            R_spec = 8314.462618 / mw

        speed_of_sound_exit = math.sqrt(gamma * R_spec * te)
        ve = exit_mach * speed_of_sound_exit

        # 4. Ideal Thrust Coefficient (NASA SP-125 Eq. 1-33)
        cf_ideal = calculate_thrust_coefficient(gamma, self.chamber_pressure, pe, pa, self.expansion_ratio)
        cf_vac = calculate_thrust_coefficient(gamma, self.chamber_pressure, pe, 0.0, self.expansion_ratio)
        cf_sl = calculate_thrust_coefficient(gamma, self.chamber_pressure, pe, SEA_LEVEL_PA, self.expansion_ratio)

        # 5. Contour Length & Divergence Efficiency
        l_cone = calculate_conical_nozzle_length(rt, re, alpha_deg=15.0)

        if self.nozzle_type.lower() == "bell":
            l_nozzle = self.bell_fractional_length * l_cone
            div_res = divergence_efficiency(
                contour="bell",
                exit_half_angle_deg=self.exit_wall_angle_deg,
                initial_wall_angle_deg=self.initial_wall_angle_deg,
            )
            lambda_div = float(div_res.value)
            contour_coords = generate_rao_bell_coordinates(
                rt, re,
                fractional_length=self.bell_fractional_length,
                theta_n_deg=self.initial_wall_angle_deg,
                theta_e_deg=self.exit_wall_angle_deg,
                num_points=60,
            )
        else:
            l_nozzle = calculate_conical_nozzle_length(rt, re, alpha_deg=self.divergent_half_angle_deg)
            div_res = divergence_efficiency(
                contour="conical",
                exit_half_angle_deg=self.divergent_half_angle_deg,
            )
            lambda_div = float(div_res.value)
            contour_coords = [
                (i * l_nozzle / 60.0, rt + i * (re - rt) / 60.0)
                for i in range(61)
            ]

        # 6. Boundary Layer Displacement Thickness (Schlichting 1/7th-power turbulent BL)
        rho_e = pe / (R_spec * te) if (R_spec * te) > 0 else 0.1
        mu_e = 4.0e-5
        bl_res = boundary_layer_displacement(
            running_length_m=l_nozzle,
            velocity_m_s=ve,
            density_kg_m3=rho_e,
            viscosity_Pa_s=mu_e,
            geometric_area_m2=ae,
            wall_radius_m=re,
        )
        delta_star = bl_res.value["displacement_thickness_m"] if bl_res.value else 0.0
        area_deficit = bl_res.value["area_deficit_fraction"] if bl_res.value else 0.0

        # 7. Stark (2005) Tri-Criteria Flow Separation Assessment
        sep_res = separation_assessment(
            exit_pressure_Pa=pe,
            ambient_pressure_Pa=pa,
            exit_mach=exit_mach,
            chamber_pressure_Pa=self.chamber_pressure,
        )
        sep_verdicts = sep_res.value if isinstance(sep_res.value, dict) else {}

        # 8. Delivered Performance Assembly
        deliv_res = delivered_thrust_coefficient(
            ideal_thrust_coefficient=cf_ideal,
            divergence=div_res,
            separation=sep_res,
            boundary_layer=bl_res,
        )
        cf_deliv = float(deliv_res.value)

        g0 = 9.80665
        mdot_throat = (self.chamber_pressure * at) / cstar_val

        thrust_deliv = cf_deliv * self.chamber_pressure * at
        isp_deliv = (cstar_val * cf_deliv) / g0

        thrust_ideal = cf_ideal * self.chamber_pressure * at
        isp_ideal = (cstar_val * cf_ideal) / g0

        deliv_vac_res = delivered_thrust_coefficient(
            ideal_thrust_coefficient=cf_vac,
            divergence=div_res,
            boundary_layer=bl_res,
        )
        cf_vac_deliv = float(deliv_vac_res.value)
        thrust_vac = cf_vac_deliv * self.chamber_pressure * at
        isp_vac = (cstar_val * cf_vac_deliv) / g0

        sep_sl_res = separation_assessment(
            exit_pressure_Pa=pe,
            ambient_pressure_Pa=SEA_LEVEL_PA,
            exit_mach=exit_mach,
            chamber_pressure_Pa=self.chamber_pressure,
        )
        deliv_sl_res = delivered_thrust_coefficient(
            ideal_thrust_coefficient=cf_sl,
            divergence=div_res,
            separation=sep_sl_res,
            boundary_layer=bl_res,
        )
        cf_sl_deliv = float(deliv_sl_res.value)
        thrust_sl = cf_sl_deliv * self.chamber_pressure * at
        isp_sl = (cstar_val * cf_sl_deliv) / g0

        if sep_res.status == Status.PASS:
            sep_summary = "ATTACHED (All criteria pass: Summerfield, Schmucker, Kalt-Badal)"
        elif sep_res.status == Status.FAIL:
            sep_summary = "SEPARATED (All criteria predict flow separation)"
        elif sep_res.status == Status.CONFLICT_UNRESOLVED:
            sep_summary = "CONFLICT_UNRESOLVED (Summerfield/Kalt-Badal vs Schmucker disagree)"
        else:
            sep_summary = str(sep_res.status)

        return NozzleResult(
            design=self,
            expansion_ratio=self.expansion_ratio,
            throat_diameter_m=dt,
            throat_area_m2=at,
            exit_diameter_m=de,
            exit_area_m2=ae,
            nozzle_length_m=l_nozzle,
            exit_mach=exit_mach,
            exit_pressure_Pa=pe,
            exit_temperature_K=te,
            exit_velocity_m_s=ve,
            ambient_pressure_Pa=pa,
            altitude_m=self.altitude,
            mass_flow_kg_s=mdot_throat,
            c_star_m_s=cstar_val,
            gamma=gamma,
            ideal_thrust_coefficient=cf_ideal,
            delivered_thrust_coefficient=cf_deliv,
            divergence_efficiency=lambda_div,
            boundary_layer_displacement_m=delta_star,
            boundary_layer_area_deficit=area_deficit,
            thrust_ideal_N=thrust_ideal,
            thrust_delivered_N=thrust_deliv,
            isp_ideal_s=isp_ideal,
            isp_delivered_s=isp_deliv,
            vacuum_thrust_coefficient=cf_vac_deliv,
            vacuum_thrust_N=thrust_vac,
            vacuum_isp_s=isp_vac,
            sea_level_thrust_coefficient=cf_sl_deliv,
            sea_level_thrust_N=thrust_sl,
            sea_level_isp_s=isp_sl,
            separation_status=sep_summary,
            separation_verdicts=sep_verdicts,
            contour_coordinates=contour_coords,
            provenance_results={
                "divergence": div_res,
                "separation": sep_res,
                "boundary_layer": bl_res,
                "delivered_cf": deliv_res,
            },
        )


@dataclass
class NozzleResult:
    """
    Complete aerodynamic, geometric, and altitude performance results of a supersonic nozzle.
    """

    design: NozzleDesign
    expansion_ratio: float
    throat_diameter_m: float
    throat_area_m2: float
    exit_diameter_m: float
    exit_area_m2: float
    nozzle_length_m: float
    exit_mach: float
    exit_pressure_Pa: float
    exit_temperature_K: float
    exit_velocity_m_s: float
    ambient_pressure_Pa: float
    altitude_m: float
    mass_flow_kg_s: float
    c_star_m_s: float
    gamma: float
    ideal_thrust_coefficient: float
    delivered_thrust_coefficient: float
    divergence_efficiency: float
    boundary_layer_displacement_m: float
    boundary_layer_area_deficit: float
    thrust_ideal_N: float
    thrust_delivered_N: float
    isp_ideal_s: float
    isp_delivered_s: float
    vacuum_thrust_coefficient: float
    vacuum_thrust_N: float
    vacuum_isp_s: float
    sea_level_thrust_coefficient: float
    sea_level_thrust_N: float
    sea_level_isp_s: float
    separation_status: str
    separation_verdicts: Dict[str, Any]
    contour_coordinates: List[Tuple[float, float]]
    provenance_results: Dict[str, Result] = field(default_factory=dict)

    @property
    def throat_diameter_mm(self) -> float:
        return self.throat_diameter_m * 1000.0

    @property
    def exit_diameter_mm(self) -> float:
        return self.exit_diameter_m * 1000.0

    @property
    def nozzle_length_mm(self) -> float:
        return self.nozzle_length_m * 1000.0

    @property
    def exit_pressure_bar(self) -> float:
        return self.exit_pressure_Pa / 1.0e5

    @property
    def exit_pressure_kPa(self) -> float:
        return self.exit_pressure_Pa / 1000.0

    @property
    def ambient_pressure_bar(self) -> float:
        return self.ambient_pressure_Pa / 1.0e5

    @property
    def thrust_delivered_kN(self) -> float:
        return self.thrust_delivered_N / 1000.0

    @property
    def thrust_vac_kN(self) -> float:
        return self.vacuum_thrust_N / 1000.0

    @property
    def thrust_sl_kN(self) -> float:
        return self.sea_level_thrust_N / 1000.0

    def summary(self) -> str:
        """Human-readable engineering report."""
        lines = [
            "=" * 78,
            "NAVRONIS PROPULSION: SUPERSONIC NOZZLE & ALTITUDE PERFORMANCE (DAY 4)",
            "=" * 78,
            f"Contour Type:             {self.design.nozzle_type.upper()} "
            f"({'80% Rao Parabolic Bell' if self.design.nozzle_type == 'bell' else f'{self.design.divergent_half_angle_deg} deg Conical'})",
            f"Propellant:               {self.design.propellant} (gamma = {self.gamma:.3f}, c* = {self.c_star_m_s:.1f} m/s)",
            f"Chamber Pressure (Pc):    {self.design.chamber_pressure / 1e5:.1f} bar ({self.design.chamber_pressure / 1e6:.2f} MPa)",
            f"Ambient Pressure (Pa):    {self.ambient_pressure_bar:.3f} bar ({self.ambient_pressure_Pa/1000:.1f} kPa) at z = {self.altitude_m:.0f} m",
            f"Expansion Ratio (eps):    {self.expansion_ratio:.1f}",
            "-" * 78,
            "GEOMETRIC DIMENSIONS:",
            f"  Throat Diameter (Dt):   {self.throat_diameter_mm:.2f} mm",
            f"  Exit Diameter (De):     {self.exit_diameter_mm:.2f} mm",
            f"  Nozzle Axial Length:    {self.nozzle_length_mm:.1f} mm",
            "-" * 78,
            "1D SUPERSONIC EXIT STATE:",
            f"  Exit Mach Number (Me):  {self.exit_mach:.3f} (M > 1)",
            f"  Exit Static Pressure:   {self.exit_pressure_kPa:.2f} kPa ({self.exit_pressure_bar:.4f} bar)",
            f"  Exit Static Temp (Te):  {self.exit_temperature_K:.1f} K ({self.exit_temperature_K - 273.15:.1f} deg C)",
            f"  Exit Velocity (ve):     {self.exit_velocity_m_s:.1f} m/s",
            f"  Pressure Ratio (Pe/Pa): {self.exit_pressure_Pa / self.ambient_pressure_Pa:.3f}",
            "-" * 78,
            "EFFICIENCY & LOSS FACTORS:",
            f"  Divergence Factor (lambda): {self.divergence_efficiency:.5f} (Thrust loss: {100*(1-self.divergence_efficiency):.2f}%)",
            f"  Boundary Layer delta*:      {self.boundary_layer_displacement_m * 1e6:.1f} um (Area deficit: {100*self.boundary_layer_area_deficit:.2f}%)",
            f"  Delivered Cf Factor:        {self.delivered_thrust_coefficient:.4f} (vs Ideal Cf: {self.ideal_thrust_coefficient:.4f})",
            "-" * 78,
            "FLOW SEPARATION ASSESSMENT (Stark 2005):",
            f"  Overall Status:         {self.separation_status}",
        ]
        for crit, v in self.separation_verdicts.items():
            lines.append(f"    - {crit:<23}: {v['verdict']:<9} (p_sep = {v['p_sep_Pa']/1000:6.1f} kPa, margin = {v['margin_ratio']:.2f})")
        lines.extend([
            "-" * 78,
            "DELIVERED PERFORMANCE ACROSS FLIGHT REGIMES:",
            f"  At Operating Altitude:  F = {self.thrust_delivered_kN:.2f} kN | Isp = {self.isp_delivered_s:.1f} s",
            f"  At Sea Level (0 km):    F = {self.thrust_sl_kN:.2f} kN | Isp = {self.sea_level_isp_s:.1f} s",
            f"  In Vacuum (Space):      F = {self.thrust_vac_kN:.2f} kN | Isp = {self.vacuum_isp_s:.1f} s",
            "=" * 78,
        ])
        return "\n".join(lines)

    def explain(self, parameter: str) -> str:
        """
        Retrieves mathematical formula, primary source citation, assumptions,
        and validity boundaries for a given computed parameter.
        """
        key = parameter.lower().strip()
        if "divergence" in key or "lambda" in key:
            res = self.provenance_results.get("divergence")
            if res:
                return (
                    f"PARAMETER: Divergence Efficiency (lambda)\n"
                    f"EQUATION ID: {res.equation_id}\n"
                    f"FORMULA: lambda_bell = (1 + cos((theta_n + theta_e)/2)) / 2  (Rao 1958)\n"
                    f"PRIMARY SOURCE: {res.source}\n"
                    f"STATUS: {res.status.name}\n"
                    f"NOTES: {res.notes}"
                )
        elif "separation" in key:
            res = self.provenance_results.get("separation")
            if res:
                return (
                    f"PARAMETER: Flow Separation Assessment\n"
                    f"EQUATION ID: {res.equation_id}\n"
                    f"PRIMARY CITATION: {STARK_2005}\n"
                    f"CRITERIA EVALUATED: Summerfield (1954), Schmucker (1973), Kalt-Badal (1965)\n"
                    f"STATUS: {res.status.name}\n"
                    f"NOTES: {res.notes}"
                )
        elif "boundary" in key or "displacement" in key or "bl" in key:
            res = self.provenance_results.get("boundary_layer")
            if res:
                return (
                    f"PARAMETER: Boundary Layer Displacement Thickness\n"
                    f"FORMULA: delta = 0.37*x / Re^0.2, delta* = delta / 8\n"
                    f"PRIMARY SOURCE: {res.source}\n"
                    f"STATUS: {res.status.name}\n"
                    f"NOTES: {res.notes}"
                )
        elif "cf" in key or "thrust_coefficient" in key:
            res = self.provenance_results.get("delivered_cf")
            if res:
                return (
                    f"PARAMETER: Delivered Thrust Coefficient\n"
                    f"FORMULA: Cf_deliv = lambda * (1 - delta_A/A) * Cf_ideal\n"
                    f"PRIMARY CITATIONS: NASA SP-125 (Ideal Cf), Rao (1958), Schlichting (1979)\n"
                    f"STATUS: {res.status.name}\n"
                    f"NOTES: {res.notes}"
                )
        return f"Parameter '{parameter}' not found in nozzle provenance registry."

    def to_dict(self) -> Dict[str, Any]:
        """Convert results to dictionary."""
        return {
            "expansion_ratio": self.expansion_ratio,
            "throat_diameter_mm": self.throat_diameter_mm,
            "exit_diameter_mm": self.exit_diameter_mm,
            "nozzle_length_mm": self.nozzle_length_mm,
            "exit_mach": self.exit_mach,
            "exit_pressure_kPa": self.exit_pressure_kPa,
            "exit_velocity_m_s": self.exit_velocity_m_s,
            "ambient_pressure_kPa": self.ambient_pressure_Pa / 1000.0,
            "altitude_m": self.altitude_m,
            "divergence_efficiency": self.divergence_efficiency,
            "boundary_layer_displacement_um": self.boundary_layer_displacement_m * 1e6,
            "delivered_thrust_coefficient": self.delivered_thrust_coefficient,
            "thrust_delivered_kN": self.thrust_delivered_kN,
            "isp_delivered_s": self.isp_delivered_s,
            "vacuum_thrust_kN": self.thrust_vac_kN,
            "vacuum_isp_s": self.vacuum_isp_s,
            "sea_level_thrust_kN": self.thrust_sl_kN,
            "sea_level_isp_s": self.sea_level_isp_s,
            "separation_status": self.separation_status,
            "separation_verdicts": self.separation_verdicts,
        }

    def export_json(self, filepath: str) -> None:
        """Export nozzle results to JSON file."""
        data = self.to_dict()
        data["contour_coordinates"] = self.contour_coordinates
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    def export_csv(self, filepath: str) -> None:
        """Export (X, Y) contour coordinates to CSV."""
        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["x_m", "y_m", "x_mm", "y_mm"])
            for x, y in self.contour_coordinates:
                writer.writerow([x, y, x * 1000.0, y * 1000.0])
