"""
Kryptonis Propulsion: High-Level Combustor Design API (Day 1 Release)
====================================================================
Provides an intuitive, source-traceable preliminary liquid rocket combustion chamber
sizing API with explicit provenance, units, and export capabilities.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

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
    convergent_length,
)
from kryptonis.propulsion_equations.combustion import chamber_bulk_residence_time
from kryptonis.propulsion_equations.chamber_acoustics import (
    first_tangential_frequency,
    first_radial_frequency,
    first_longitudinal_frequency,
)
from kryptonis.propulsion_equations.bartz import bartz_film_coefficient, recovery_temperature
from kryptonis.propulsion_equations.profile import generate_chamber_profile, ChamberProfile
from kryptonis.propulsion_equations.export import export_json, export_csv
from kryptonis.propulsion_equations.plotting import plot_chamber_profile


@dataclass
class CombustorDesign:
    """Preliminary liquid-rocket combustion chamber sizing specification.

    Parameters
    ----------
    thrust : float
        Design thrust in Newtons (e.g. 30,000 for 30 kN).
    chamber_pressure : float
        Total combustion chamber pressure in Pa (e.g. 7.0e6 for 70 bar, 12.0e6 for 120 bar).
    mixture_ratio : float
        Oxidizer-to-fuel mass ratio (O/F), e.g. 2.6 for LOX/CH4, 2.5 for LOX/RP-1, 6.0 for LOX/LH2.
    propellant : str
        Propellant key, e.g. "LOX/LCH4", "LOX/CH4", "LOX/RP-1", "LOX/LH2".
    c_star : float, optional
        Characteristic exhaust velocity in m/s. If omitted, computed from thermochemical defaults.
    contraction_ratio : float, optional
        Contraction area ratio (Ac / At). If omitted, evaluated from the Humble correlation.
    characteristic_length : float, optional
        Chamber characteristic length L* in meters. If omitted, uses standard propellant defaults.
    convergent_half_angle_deg : float
        Convergent cone half-angle in degrees (default: 30.0 deg).
    wall_yield_strength : float
        Liner material yield strength in Pa (default: 280.0 MPa for OFHC Copper / Cu-alloy).
    safety_factor : float
        Structural safety factor for thin-shell hoop stress (default: 1.5).
    """

    thrust: float
    chamber_pressure: float
    mixture_ratio: float
    propellant: str = "LOX/CH4"
    c_star: Optional[float] = None
    contraction_ratio: Optional[float] = None
    characteristic_length: Optional[float] = None
    convergent_half_angle_deg: float = 30.0
    wall_yield_strength: float = 280.0e6
    safety_factor: float = 1.5

    def solve(self) -> "CombustorResult":
        """Execute the closed-form analytical sizing sequence."""
        norm_prop = self.propellant.upper().replace("LCH4", "CH4").replace("KEROSENE", "RP-1")
        
        prop_defaults = {
            "LOX/RP-1": {"gamma": 1.22, "mw": 0.0235, "tc": 3600.0, "l_star": 1.05, "c_star": 1780.0},
            "LOX/CH4": {"gamma": 1.20, "mw": 0.0220, "tc": 3450.0, "l_star": 1.00, "c_star": 1820.0},
            "LOX/LH2": {"gamma": 1.23, "mw": 0.0150, "tc": 3250.0, "l_star": 0.85, "c_star": 2350.0},
        }
        defaults = prop_defaults.get(norm_prop, prop_defaults["LOX/CH4"])

        gamma = defaults["gamma"]
        mw = defaults["mw"]
        tc = defaults["tc"]
        effective_l_star = self.characteristic_length if self.characteristic_length is not None else defaults["l_star"]

        # 1. c* velocity
        if self.c_star is not None:
            cstar_val = self.c_star
        else:
            cstar_res = c_star_ideal(gamma=gamma, molar_mass_kg_per_mol=mw, chamber_temperature_K=tc)
            cstar_val = cstar_res.value

        # 2. Total and split mass flow rates (estimated via standard expansion Cf ~ 1.75)
        cf_est = 1.75
        at_est = self.thrust / (self.chamber_pressure * cf_est)
        total_mdot = (self.chamber_pressure * at_est) / cstar_val
        fuel_mdot = total_mdot / (1.0 + self.mixture_ratio)
        ox_mdot = total_mdot - fuel_mdot

        # 3. Choked Throat Sizing
        at_res = throat_area(mass_flow_kg_s=total_mdot, c_star_m_s=cstar_val, chamber_pressure_Pa=self.chamber_pressure)
        at_val = at_res.value
        dt_res = throat_diameter(throat_area_m2=at_val)
        dt_val = dt_res.value

        # 4. Contraction Ratio & Chamber Cross-Section
        if self.contraction_ratio is not None:
            cr_val = self.contraction_ratio
        else:
            cr_res = contraction_ratio(throat_diameter_m=dt_val)
            cr_val = cr_res.value

        dc_res = chamber_diameter(throat_diameter_m=dt_val, contraction_ratio_=cr_val)
        dc_val = dc_res.value
        ac_val = at_val * cr_val

        # 5. Chamber Volume & Convergent Geometry
        vc_res = chamber_volume(l_star_m=effective_l_star, throat_area_m2=at_val)
        vc_val = vc_res.value

        vconv_res = convergent_volume(
            throat_diameter_m=dt_val,
            chamber_diameter_m=dc_val,
            half_angle_deg=self.convergent_half_angle_deg,
        )
        vconv_val = vconv_res.value

        lconv_val = (dc_val - dt_val) / (2.0 * math.tan(math.radians(self.convergent_half_angle_deg)))

        # 6. Cylindrical Length & Residence Time
        lcyl_res = cylinder_length(
            chamber_volume_m3=vc_val,
            convergent_volume_m3=vconv_val,
            chamber_area_m2=ac_val,
        )
        lcyl_val = lcyl_res.value if not math.isnan(lcyl_res.value) else 0.0
        ltot_val = lcyl_val + lconv_val

        r_spec = 8314.4626 / (mw * 1000.0)
        rho_c = self.chamber_pressure / (r_spec * tc)
        tau_res = chamber_bulk_residence_time(chamber_volume_m3=vc_val, density_kg_m3=rho_c, mdot_kg_s=total_mdot)

        # 7. Acoustic Stability Frequencies (Bessel Mode Roots)
        a_sound = math.sqrt(gamma * r_spec * tc)
        f_1t = first_tangential_frequency(speed_of_sound_m_s=a_sound, chamber_diameter_m=dc_val).value
        f_1r = first_radial_frequency(speed_of_sound_m_s=a_sound, chamber_diameter_m=dc_val).value
        f_1l = first_longitudinal_frequency(speed_of_sound_m_s=a_sound, chamber_length_m=ltot_val).value

        # 8. Preliminary Thermal Screening (Bartz Peak Throat Heat Flux)
        hg_res = bartz_film_coefficient(
            throat_diameter_m=dt_val,
            chamber_pressure_Pa=self.chamber_pressure,
            c_star_m_s=cstar_val,
            throat_radius_curvature_m=dt_val * 1.5,
            area_ratio_local_over_throat=1.0,
            gamma=gamma,
            mach=1.0,
            specific_heat_J_kgK=2600.0,
            prandtl=0.82,
            viscosity_Pa_s=8.8e-5,
            stagnation_temperature_K=tc,
            wall_temperature_K=850.0,
        )
        taw_res = recovery_temperature(stagnation_temperature_K=tc, mach=1.0, gamma=gamma, prandtl=0.82)
        q_throat = hg_res.value * (taw_res.value - 850.0)

        # 9. Mechanical Wall Thickness Screening (Thin-Shell Hoop)
        allowable_stress = self.wall_yield_strength / self.safety_factor
        t_wall_min = (self.chamber_pressure * (dc_val / 2.0)) / allowable_stress
        t_wall_rec = t_wall_min * 1.25

        provenance_records = {
            "total_mass_flow": {
                "equation": "mdot = (P_c * A_t_est) / c*",
                "source": "NASA SP-125 Eq. (1-9)",
                "status": "PRELIMINARY_DESIGN",
                "assumptions": "Isentropic choked flow with estimated thrust coefficient Cf ~ 1.75",
            },
            "throat_diameter": {
                "equation": "D_t = sqrt(4 * A_t / pi)",
                "source": "NASA SP-125 Chapter 4 / Huzel & Huang",
                "status": "PRELIMINARY_DESIGN",
                "assumptions": "Axisymmetric circular sonic throat",
            },
            "chamber_diameter": {
                "equation": "D_c = D_t * sqrt(eps_c)",
                "source": "Geometric identity from contraction ratio",
                "status": "PRELIMINARY_DESIGN",
                "assumptions": "Cylindrical combustion chamber",
            },
            "chamber_length": {
                "equation": "L_c = L_cyl + L_conv",
                "source": "NASA SP-125 Eq. (4-5)",
                "status": "PRELIMINARY_DESIGN",
                "assumptions": "Straight cylindrical barrel + truncated conical convergent section",
            },
            "heat_flux_screen": {
                "equation": "q_throat = h_g * (T_aw - T_wg)",
                "source": "Bartz, D.R. (1957) Jet Propulsion 27(1)",
                "status": "SCREENING_ONLY",
                "assumptions": "Isentropic core flow, sonic throat, assumed gas wall temp 850 K",
            },
            "acoustic_modes": {
                "equation": "f_1T = 1.8412 * a / (pi * D_c)",
                "source": "NASA SP-194 Liquid Rocket Combustion Instability",
                "status": "SCREENING_ONLY",
                "assumptions": "Rigid cylinder cavity, uniform gas, zero mean flow",
            },
            "wall_thickness_screen": {
                "equation": "t_w = P_c * R_c / (sigma_yield / SF)",
                "source": "ASME Section VIII Div 1 / Thin-wall hoop stress",
                "status": "SCREENING_ONLY",
                "assumptions": "Thin cylinder, uniform internal pressure, zero thermal gradient stress",
            },
        }

        return CombustorResult(
            design=self,
            total_mass_flow=total_mdot,
            fuel_mass_flow=fuel_mdot,
            oxidizer_mass_flow=ox_mdot,
            throat_area=at_val,
            throat_diameter=dt_val,
            chamber_area=ac_val,
            chamber_diameter=dc_val,
            contraction_ratio=cr_val,
            characteristic_length=effective_l_star,
            chamber_volume=vc_val,
            convergent_length=lconv_val,
            cylindrical_length=lcyl_val,
            chamber_length=ltot_val,
            stay_time=tau_res.value,
            c_star_ideal=cstar_val,
            heat_flux_screen=q_throat,
            acoustic_modes={"1T_Hz": f_1t, "1R_Hz": f_1r, "1L_Hz": f_1l},
            wall_thickness_screen=t_wall_min,
            wall_thickness_recommended=t_wall_rec,
            provenance=provenance_records,
        )


@dataclass
class CombustorResult:
    """Sizing results for the preliminary combustion chamber."""

    design: CombustorDesign
    total_mass_flow: float
    fuel_mass_flow: float
    oxidizer_mass_flow: float
    throat_area: float
    throat_diameter: float
    chamber_area: float
    chamber_diameter: float
    contraction_ratio: float
    characteristic_length: float
    chamber_volume: float
    convergent_length: float
    cylindrical_length: float
    chamber_length: float
    stay_time: float
    c_star_ideal: float
    heat_flux_screen: float
    acoustic_modes: Dict[str, float]
    wall_thickness_screen: float
    wall_thickness_recommended: float
    provenance: Dict[str, Dict[str, str]] = field(default_factory=dict)

    def explain(self, quantity_name: str) -> str:
        """Return the mathematical and literature provenance of an output."""
        rec = self.provenance.get(quantity_name)
        if not rec:
            return f"No provenance recorded for '{quantity_name}'."
        
        val = getattr(self, quantity_name, "N/A")
        if isinstance(val, float):
            val_str = f"{val:.4e}" if (abs(val) < 1e-3 or abs(val) > 1e4) else f"{val:.4f}"
        else:
            val_str = str(val)

        return (
            f"Quantity:     {quantity_name}\n"
            f"Value:        {val_str}\n"
            f"Equation:     {rec['equation']}\n"
            f"Source:       {rec['source']}\n"
            f"Status:       {rec['status']}\n"
            f"Assumptions:  {rec['assumptions']}"
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert sizing results to a serializable dictionary."""
        return {
            "inputs": {
                "thrust_N": self.design.thrust,
                "chamber_pressure_Pa": self.design.chamber_pressure,
                "mixture_ratio": self.design.mixture_ratio,
                "propellant": self.design.propellant,
                "convergent_half_angle_deg": self.design.convergent_half_angle_deg,
            },
            "combustor": {
                "throat_diameter_mm": self.throat_diameter * 1e3,
                "throat_area_cm2": self.throat_area * 1e4,
                "chamber_diameter_mm": self.chamber_diameter * 1e3,
                "chamber_area_cm2": self.chamber_area * 1e4,
                "contraction_ratio": self.contraction_ratio,
                "chamber_length_mm": self.chamber_length * 1e3,
                "cylindrical_length_mm": self.cylindrical_length * 1e3,
                "convergent_length_mm": self.convergent_length * 1e3,
                "chamber_volume_liters": self.chamber_volume * 1e3,
                "characteristic_length_m": self.characteristic_length,
                "stay_time_ms": self.stay_time * 1e3,
                "c_star_m_s": self.c_star_ideal,
            },
            "mass_flow": {
                "total_kg_s": self.total_mass_flow,
                "fuel_kg_s": self.fuel_mass_flow,
                "oxidizer_kg_s": self.oxidizer_mass_flow,
            },
            "thermal_screen": {
                "peak_throat_heat_flux_MW_m2": self.heat_flux_screen / 1e6,
            },
            "acoustic_screen": self.acoustic_modes,
            "structural_screen": {
                "min_wall_thickness_mm": self.wall_thickness_screen * 1e3,
                "recommended_wall_thickness_mm": self.wall_thickness_recommended * 1e3,
            },
            "provenance": self.provenance,
        }

    def to_json(self, path: str) -> None:
        """Export results to a structured JSON file."""
        with open(path, "w", encoding="utf-8") as fp:
            json.dump(self.to_dict(), fp, indent=2)

    def to_csv(self, path: str) -> None:
        """Export summary parameters to a CSV file."""
        d = self.to_dict()
        rows = []
        for section, values in d.items():
            if isinstance(values, dict) and section != "provenance":
                for k, v in values.items():
                    rows.append((section, k, str(v)))
        
        with open(path, "w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(["Section", "Parameter", "Value"])
            writer.writerows(rows)

    def plot(self, save_path: Optional[str] = None, show: bool = True) -> None:
        """Generate and display the 2D cross-section plot of the chamber."""
        prof = generate_chamber_profile(
            throat_diameter_m=self.throat_diameter,
            chamber_diameter_m=self.chamber_diameter,
            convergent_half_angle_deg=self.design.convergent_half_angle_deg,
            cylindrical_length_m=self.cylindrical_length,
            expansion_ratio=1.0,
            nozzle_type="conical",
        )
        plot_chamber_profile(
            prof,
            title=f"Kryptonis Combustor Profile ({self.design.thrust/1e3:.0f} kN {self.design.propellant})",
            show_dimensions=True,
            show_wall_thickness=True,
            wall_thickness_mm=self.wall_thickness_recommended * 1e3,
            save_path=save_path,
            show=show,
        )
