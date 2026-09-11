"""
Navronis Propulsion: Injector Head & Atomization Analytical Sizing Engine
========================================================================
Focuses strictly on the 4 canonical liquid rocket engine injector families:
1. Shear Coaxial Injectors (SSME, RL10, Vulcain, Raptor - NASA SP-125, Yang 2004)
2. Swirl Coaxial Injectors (RD-170, RD-180, NK-33 - Bazarov & Yang 1998, Lefebvre 1989)
3. Pintle Injectors (Apollo LMDE, Merlin 1D, TRW TR-106 - Dressler 2000, Heister 2019)
4. Unlike Impinging Doublets (Apollo SPS, Titan, Viking - Rupe 1956, Ingebo 1958)

Includes base orifice hydraulics and NASA SP-194 chugging decoupling stability criterion.
All formulations feature strict SI units, sanity validation, and explicit literature citations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# --------------------------------------------------------------------------
# Literature References
# --------------------------------------------------------------------------
CITATIONS = {
    "orifice_hydraulics": "NASA SP-125 (1971), Huzel & Huang, Eq. (4-17) p. 118",
    "chugging_stiffness": "NASA SP-194 (1972), Liquid Propellant Rocket Combustion Instability, p. 55",
    "shear_coaxial": "Yang, V. et al. (2004), Liquid Rocket Thrust Chambers, AIAA Progress in Astro, Ch. 3",
    "coaxial_atomization": "Lorenzetto, P. R. & Lefebvre, A. H. (1977), Measurements of Drop Size on Coaxial Airblast Atomizers, AIAA J.",
    "pintle_tmr": "Dressler, G. A. & Bauer, J. M. (2000), TRW Pintle Engine Heritage and Performance, AIAA-2000-3871",
    "pintle_spray_angle": "Heister, S. D. et al. (2019), Rocket Propulsion, Cambridge University Press, Ch. 6",
    "swirl_coaxial": "Bazarov, V. G. & Yang, V. (1998), Liquid-Propellant Rocket Engine Injectors, AIAA J. Prop. & Power, 14(5)",
    "swirl_atomization": "Lefebvre, A. H. (1989), Atomization and Sprays, Hemisphere Publishing, Eq. (6.33) p. 215",
    "rupe_mixing": "Rupe, J. H. (1956), The Liquid Phase Mixing of a Pair of Impinging Streams, JPL Report No. 20-195",
    "ingebo_smd": "Ingebo, R. D. (1958), Drop-Size Distributions for Impinging Jet Atomizers, NACA TN-4222",
}


# --------------------------------------------------------------------------
# 1. Base Orifice Hydraulics & Stiffness
# --------------------------------------------------------------------------

def orifice_area(mass_flow: float, cd: float, rho: float, delta_p: float) -> float:
    """Calculate single or total orifice area from Torricelli-Bernoulli continuity.
    
    A_o = m_dot / (Cd * sqrt(2 * rho * delta_p))
    """
    if mass_flow <= 0 or cd <= 0 or rho <= 0 or delta_p <= 0:
        raise ValueError("mass_flow, cd, rho, and delta_p must be positive and finite")
    return mass_flow / (cd * math.sqrt(2.0 * rho * delta_p))


def orifice_diameter(area: float) -> float:
    """Circular orifice diameter from cross-sectional area."""
    if area <= 0:
        raise ValueError("area must be positive")
    return math.sqrt(4.0 * area / math.pi)


def orifice_velocity(cd: float, delta_p: float, rho: float) -> float:
    """Liquid or gas jet discharge velocity: V = Cd * sqrt(2 * delta_p / rho)."""
    if cd <= 0 or delta_p <= 0 or rho <= 0:
        raise ValueError("cd, delta_p, and rho must be positive")
    return cd * math.sqrt(2.0 * delta_p / rho)


def injector_pressure_drop_stiffness(delta_p: float, pc: float) -> float:
    """Calculate injector pressure drop ratio: delta_p / Pc."""
    if pc <= 0:
        raise ValueError("Chamber pressure Pc must be positive")
    return delta_p / pc


# --------------------------------------------------------------------------
# 2. Shear Coaxial Injector Sizing (LOX/CH4, LOX/LH2)
# --------------------------------------------------------------------------

@dataclass
class ShearCoaxialElement:
    """Single shear coaxial injector element geometry and performance."""
    element_id: int
    m_dot_ox_per_elem: float       # kg/s
    m_dot_fuel_per_elem: float     # kg/s
    d_ox_inner: float              # Central liquid post ID (m)
    d_ox_outer: float              # Central liquid post OD (m)
    t_post: float                  # Post wall thickness (m)
    d_ann_inner: float             # Annulus sleeve ID (m)
    a_ox: float                    # Liquid post area (m²)
    a_ann: float                   # Gas annular area (m²)
    v_ox: float                    # Liquid velocity (m/s)
    v_fuel: float                  # Gas velocity (m/s)
    momentum_flux_ratio_J: float   # (rho_g * v_g^2) / (rho_l * v_l^2)
    velocity_ratio_VR: float       # v_g / v_l
    recess_length: float           # Recess length (m)
    recess_ratio_RL: float         # L_recess / d_ox_inner
    smd_um: float                  # Sauter Mean Diameter D32 (micrometers)


def size_shear_coaxial(
    m_dot_ox: float,
    m_dot_fuel: float,
    rho_ox: float,
    rho_fuel: float,
    delta_p_ox: float,
    delta_p_fuel: float,
    cd_ox: float = 0.75,
    cd_fuel: float = 0.80,
    n_elements: int = 19,
    post_wall_thickness: float = 0.0006,  # 0.6 mm
    post_recess_ratio: float = 1.5,       # L_recess / D_post_id
    surface_tension_ox: float = 0.013,    # N/m (liquid oxygen at ~90K)
    viscosity_ox: float = 1.9e-4,         # Pa*s
) -> Dict[str, Any]:
    """Size an array of shear coaxial injector elements with atomization SMD.
    
    Literature Source:
    - Yang (2004), Liquid Rocket Thrust Chambers, Ch. 3
    - Lorenzetto & Lefebvre (1977), Drop size in coaxial airblast atomizers
    """
    if n_elements <= 0:
        raise ValueError("Number of elements must be >= 1")

    m_ox_elem = m_dot_ox / n_elements
    m_f_elem = m_dot_fuel / n_elements

    # Central liquid post area & diameter
    a_ox = orifice_area(m_ox_elem, cd_ox, rho_ox, delta_p_ox)
    d_ox_id = orifice_diameter(a_ox)
    d_ox_od = d_ox_id + 2.0 * post_wall_thickness

    # Gas annular area
    a_fuel = orifice_area(m_f_elem, cd_fuel, rho_fuel, delta_p_fuel)
    # a_fuel = pi/4 * (d_ann_id^2 - d_ox_od^2) => d_ann_id = sqrt(d_ox_od^2 + 4 * a_fuel / pi)
    d_ann_id = math.sqrt(d_ox_od**2 + (4.0 * a_fuel / math.pi))
    annular_gap = (d_ann_id - d_ox_od) / 2.0

    v_ox = m_ox_elem / (rho_ox * a_ox)
    v_fuel = m_f_elem / (rho_fuel * a_fuel)

    # Momentum flux ratio J = (rho_g * V_g^2) / (rho_l * V_l^2)
    dyn_pres_fuel = rho_fuel * (v_fuel**2)
    dyn_pres_ox = rho_ox * (v_ox**2)
    J = dyn_pres_fuel / dyn_pres_ox
    velocity_ratio = v_fuel / v_ox

    # Recess length
    recess_length = post_recess_ratio * d_ox_id

    # Sauter Mean Diameter D32 via Lorenzetto-Lefebvre correlation
    # We_g = rho_g * (V_g - V_l)^2 * D_l / sigma
    rel_v = max(v_fuel - v_ox, 10.0)
    we_g = (rho_fuel * (rel_v**2) * d_ox_id) / surface_tension_ox
    oh = viscosity_ox / math.sqrt(rho_ox * surface_tension_ox * d_ox_id)
    
    # SMD formula (um)
    smd_m = d_ox_id * 0.48 * (1.0 / max(we_g, 1.0))**0.4 * (1.0 + 1.0 / max(velocity_ratio, 1.0))**0.4 \
            + 0.15 * d_ox_id * (oh**0.5)
    smd_um = smd_m * 1.0e6

    return {
        "n_elements": n_elements,
        "m_dot_ox_total_kg_s": m_dot_ox,
        "m_dot_fuel_total_kg_s": m_dot_fuel,
        "post_id_mm": d_ox_id * 1e3,
        "post_od_mm": d_ox_od * 1e3,
        "annulus_id_mm": d_ann_id * 1e3,
        "annular_gap_mm": annular_gap * 1e3,
        "v_ox_m_s": v_ox,
        "v_fuel_m_s": v_fuel,
        "momentum_flux_ratio_J": J,
        "velocity_ratio_VR": velocity_ratio,
        "recess_length_mm": recess_length * 1e3,
        "smd_um": smd_um,
        "provenance": {
            "sizing": CITATIONS["shear_coaxial"],
            "atomization": CITATIONS["coaxial_atomization"],
        }
    }


# --------------------------------------------------------------------------
# 3. Swirl Coaxial Injector Sizing (RD-170, RD-180, NK-33 Style)
# --------------------------------------------------------------------------

def size_swirl_coaxial(
    m_dot_liquid: float,
    m_dot_gas: float,
    rho_liquid: float,
    rho_gas: float,
    delta_p_liquid: float,
    delta_p_gas: float,
    n_elements: int = 19,
    geometric_swirl_k: float = 3.0,
    n_tangential_inlets: int = 3,
    post_wall_thickness: float = 0.0008,
    surface_tension_liquid: float = 0.013,
    viscosity_liquid: float = 1.9e-4,
) -> Dict[str, Any]:
    """Size centrifugal pressure-swirl coaxial elements with hollow conical liquid sheet.
    
    Literature Source:
    - Bazarov & Yang (1998), Liquid-Propellant Rocket Engine Injectors, AIAA
    - Lefebvre (1989), Atomization and Sprays, Eq. (6.33)
    """
    if n_elements <= 0:
        raise ValueError("Number of elements must be >= 1")
    if geometric_swirl_k <= 0.5:
        raise ValueError("Geometric swirl characteristic K must be > 0.5")

    m_liq_elem = m_dot_liquid / n_elements
    m_gas_elem = m_dot_gas / n_elements

    # Discharge coefficient Cd via Abramovich/Bazarov swirl theory
    cd_swirl = 0.35 / (geometric_swirl_k ** 0.35)

    # Orifice throat area and diameter
    a_orifice = orifice_area(m_liq_elem, cd_swirl, rho_liquid, delta_p_liquid)
    d_orifice = orifice_diameter(a_orifice)

    # Coefficient of nozzle filling phi (gas core ratio)
    phi = 1.0 / (1.0 + 0.35 * geometric_swirl_k)
    d_core = d_orifice * math.sqrt(max(1.0 - phi, 0.01))
    film_thickness = 0.5 * (d_orifice - d_core)

    # Spray half-cone angle: tan(theta) = 0.45 * K
    tan_theta = 0.45 * geometric_swirl_k
    spray_half_angle_deg = math.degrees(math.atan(tan_theta))

    # Tangential inlet port sizing
    a_inlet_total = a_orifice / (cd_swirl * math.sqrt(2.0))
    d_inlet_port = math.sqrt(4.0 * (a_inlet_total / n_tangential_inlets) / math.pi)

    # Coaxial gas annulus sizing
    cd_gas = 0.82
    a_gas = orifice_area(m_gas_elem, cd_gas, rho_gas, delta_p_gas)
    d_post_od = d_orifice + 2.0 * post_wall_thickness
    d_gas_sleeve_id = math.sqrt(d_post_od**2 + (4.0 * a_gas / math.pi))
    annular_gas_gap = 0.5 * (d_gas_sleeve_id - d_post_od)

    v_gas = m_gas_elem / (rho_gas * a_gas)
    v_liq_axial = m_liq_elem / (rho_liquid * a_orifice * phi)

    # Sauter Mean Diameter (SMD D32) via Lefebvre (1989) Eq. (6.33)
    sigma = surface_tension_liquid
    mu_l = viscosity_liquid
    t_f = max(film_thickness, 1e-6)
    term1 = 4.52 * (((sigma * (mu_l**2)) / (rho_gas * (delta_p_liquid**2)))**0.25) * (t_f**0.25)
    term2 = 0.39 * (((sigma * rho_liquid) / (rho_gas * delta_p_liquid))**0.25) * (t_f**0.75)
    smd_m = term1 + term2
    smd_um = smd_m * 1.0e6

    return {
        "n_elements": n_elements,
        "orifice_diameter_mm": d_orifice * 1e3,
        "gas_core_diameter_mm": d_core * 1e3,
        "liquid_film_thickness_mm": film_thickness * 1e3,
        "spray_half_angle_deg": spray_half_angle_deg,
        "tangential_inlet_diameter_mm": d_inlet_port * 1e3,
        "gas_sleeve_id_mm": d_gas_sleeve_id * 1e3,
        "annular_gas_gap_mm": annular_gas_gap * 1e3,
        "gas_velocity_m_s": v_gas,
        "liquid_axial_velocity_m_s": v_liq_axial,
        "geometric_swirl_K": geometric_swirl_k,
        "smd_um": smd_um,
        "provenance": {
            "swirl_mechanics": CITATIONS["swirl_coaxial"],
            "atomization": CITATIONS["swirl_atomization"],
        }
    }


# --------------------------------------------------------------------------
# 4. Pintle Injector Sizing (Merlin, Starship, Apollo LMDE Style)
# --------------------------------------------------------------------------

def size_pintle_injector(
    m_dot_annular: float,
    m_dot_radial: float,
    rho_annular: float,
    rho_radial: float,
    delta_p_annular: float,
    delta_p_radial: float,
    pintle_diameter: float,
    cd_annular: float = 0.85,
    cd_radial: float = 0.75,
    n_radial_slots: Optional[int] = None,
    slot_width_to_height_ratio: float = 1.0,
) -> Dict[str, Any]:
    """Size a central pintle injector with annular and radial sheet geometry.
    
    Parameters
    ----------
    m_dot_annular : float
        Mass flow rate through outer annular gap (kg/s). E.g. Fuel for wall cooling.
    m_dot_radial : float
        Mass flow rate through central pintle radial slots/gap (kg/s). E.g. Oxidizer.
    pintle_diameter : float
        Outer tip diameter of the pintle shaft D_p in meters (e.g. 0.020 for 20 mm).
    """
    if pintle_diameter <= 0:
        raise ValueError("pintle_diameter must be positive")

    # Annular gap sizing: A_ann = pi * D_p * t_ann
    a_ann = orifice_area(m_dot_annular, cd_annular, rho_annular, delta_p_annular)
    t_ann = a_ann / (math.pi * pintle_diameter)
    v_ann = m_dot_annular / (rho_annular * a_ann)

    # Radial orifice sizing: A_rad = n_slots * w_s * h_s
    a_rad = orifice_area(m_dot_radial, cd_radial, rho_radial, delta_p_radial)
    v_rad = m_dot_radial / (rho_radial * a_rad)

    if n_radial_slots is not None and n_radial_slots > 0:
        # Discrete slots: a_rad / n_slots = w_s * h_s = (slot_ratio * h_s) * h_s
        a_slot = a_rad / n_radial_slots
        h_slot = math.sqrt(a_slot / slot_width_to_height_ratio)
        w_slot = slot_width_to_height_ratio * h_slot
    else:
        # Continuous annular radial slit around circumference
        h_slot = a_rad / (math.pi * pintle_diameter)
        w_slot = math.pi * pintle_diameter
        n_radial_slots = 1

    # Total Momentum Ratio (TMR)
    # TMR = (m_dot_annular * V_annular) / (m_dot_radial * V_radial)
    mom_ann = m_dot_annular * v_ann
    mom_rad = m_dot_radial * v_rad
    tmr = mom_ann / mom_rad

    # Resultant spray half-cone angle beta: beta = arccos(1 / (1 + TMR))
    # Or Heister Eq: beta_deg = degrees(arccos(1 / (1 + TMR)))
    cos_beta = 1.0 / (1.0 + tmr)
    beta_rad = math.acos(min(max(cos_beta, -1.0), 1.0))
    beta_deg = math.degrees(beta_rad)

    return {
        "pintle_diameter_mm": pintle_diameter * 1e3,
        "annular_gap_thickness_mm": t_ann * 1e3,
        "annular_velocity_m_s": v_ann,
        "radial_slot_height_mm": h_slot * 1e3,
        "radial_slot_width_mm": w_slot * 1e3,
        "radial_velocity_m_s": v_rad,
        "total_momentum_ratio_TMR": tmr,
        "spray_half_angle_deg": beta_deg,
        "n_radial_slots": n_radial_slots,
        "provenance": {
            "tmr": CITATIONS["pintle_tmr"],
            "spray_angle": CITATIONS["pintle_spray_angle"],
        }
    }


# --------------------------------------------------------------------------
# 4. Impinging Jet Doublets & Triplets Sizing
# --------------------------------------------------------------------------

def size_impinging_doublet(
    m_dot_1: float,
    m_dot_2: float,
    rho_1: float,
    rho_2: float,
    delta_p_1: float,
    delta_p_2: float,
    n_elements: int = 12,
    impingement_half_angle_deg: float = 30.0,
    cd_1: float = 0.72,
    cd_2: float = 0.72,
    surface_tension_1: float = 0.025,
) -> Dict[str, Any]:
    """Size unlike impinging doublet elements (1-on-1) with Rupe momentum balance.
    
    Literature Source:
    - Rupe (1956), Liquid Phase Mixing of Impinging Streams
    - Ingebo (1958), Drop-Size Distributions for Impinging Jets
    """
    if n_elements <= 0:
        raise ValueError("n_elements must be >= 1")

    m1_elem = m_dot_1 / n_elements
    m2_elem = m_dot_2 / n_elements

    a1 = orifice_area(m1_elem, cd_1, rho_1, delta_p_1)
    a2 = orifice_area(m2_elem, cd_2, rho_2, delta_p_2)
    d1 = orifice_diameter(a1)
    d2 = orifice_diameter(a2)

    v1 = m1_elem / (rho_1 * a1)
    v2 = m2_elem / (rho_2 * a2)

    # Rupe momentum balance parameter: (rho1 * v1^2 * d1^2) / (rho2 * v2^2 * d2^2)
    # An ideal unlike doublet has Rupe parameter = 1.0 (centered spray fan)
    rupe_param = (rho_1 * (v1**2) * (d1**2)) / (rho_2 * (v2**2) * (d2**2))

    # Free jet pre-impingement length: L_free = 6 * d1
    l_free = 6.0 * d1
    
    # Element orifice separation on face: S = 2 * L_free * sin(theta)
    theta_rad = math.radians(impingement_half_angle_deg)
    orifice_separation = 2.0 * l_free * math.sin(theta_rad)

    # Sauter Mean Diameter D32 via Ingebo relation
    # Relative velocity at collision: V_rel = 2 * V_mean * sin(theta)
    v_mean = 0.5 * (v1 + v2)
    v_coll = 2.0 * v_mean * math.sin(theta_rad)
    we = (rho_1 * (v_coll**2) * d1) / surface_tension_1
    re = (rho_1 * v_coll * d1) / 1.0e-3
    smd_m = d1 * 0.027 * (max(we / max(re, 1.0), 1e-4)**-0.25)
    smd_um = smd_m * 1.0e6

    return {
        "n_elements": n_elements,
        "orifice_diameter_1_mm": d1 * 1e3,
        "orifice_diameter_2_mm": d2 * 1e3,
        "jet_velocity_1_m_s": v1,
        "jet_velocity_2_m_s": v2,
        "rupe_momentum_parameter": rupe_param,
        "free_jet_length_mm": l_free * 1e3,
        "orifice_face_separation_mm": orifice_separation * 1e3,
        "impingement_total_angle_deg": impingement_half_angle_deg * 2.0,
        "smd_um": smd_um,
        "provenance": {
            "rupe_criterion": CITATIONS["rupe_mixing"],
            "atomization": CITATIONS["ingebo_smd"],
        }
    }


# --------------------------------------------------------------------------
# 5. High-Level Injector Design Container (Day 2 Facade)
# --------------------------------------------------------------------------

@dataclass
class InjectorDesign:
    """Preliminary liquid-rocket injector head sizing specification.
    
    Parameters
    ----------
    injector_type : str
        Type of injector: "coaxial", "pintle", or "impinging".
    chamber_pressure : float
        Total combustion chamber pressure in Pa (e.g. 3.0e6 for 30 bar).
    mass_flow_ox : float
        Oxidizer mass flow rate in kg/s.
    mass_flow_fuel : float
        Fuel mass flow rate in kg/s.
    rho_ox : float
        Liquid oxidizer density in kg/m³ (e.g. 1141 for LOX).
    rho_fuel : float
        Fuel density in kg/m³ (e.g. 422 for liquid methane, 810 for RP-1, 16.0 for gaseous CH4).
    delta_p_ratio : float
        Pressure drop ratio delta_p / Pc (default: 0.20 = 20%).
    n_elements : int
        Number of coaxial or impinging elements (default: 19).
    pintle_diameter : float
        Diameter of central pintle in meters (for pintle injectors, default: 0.022 m).
    """

    injector_type: str
    chamber_pressure: float
    mass_flow_ox: float
    mass_flow_fuel: float
    rho_ox: float
    rho_fuel: float
    delta_p_ratio: float = 0.20
    n_elements: int = 19
    pintle_diameter: float = 0.022
    cd_ox: float = 0.75
    cd_fuel: float = 0.80

    def solve(self) -> Dict[str, Any]:
        """Execute closed-form injector sizing and atomization evaluation."""
        if self.chamber_pressure <= 0:
            raise ValueError("chamber_pressure must be positive")
        if self.mass_flow_ox <= 0 or self.mass_flow_fuel <= 0:
            raise ValueError("mass flow rates must be positive")
        if self.delta_p_ratio < 0.05:
            raise ValueError("delta_p_ratio below 0.05 risks catastrophic chugging instability")

        delta_p_ox = self.delta_p_ratio * self.chamber_pressure
        delta_p_fuel = self.delta_p_ratio * self.chamber_pressure

        t_type = self.injector_type.lower().strip()
        if t_type in {"coaxial", "shear_coaxial"}:
            res = size_shear_coaxial(
                m_dot_ox=self.mass_flow_ox,
                m_dot_fuel=self.mass_flow_fuel,
                rho_ox=self.rho_ox,
                rho_fuel=self.rho_fuel,
                delta_p_ox=delta_p_ox,
                delta_p_fuel=delta_p_fuel,
                cd_ox=self.cd_ox,
                cd_fuel=self.cd_fuel,
                n_elements=self.n_elements,
            )
        elif t_type in {"swirl", "swirl_coaxial", "centrifugal"}:
            res = size_swirl_coaxial(
                m_dot_liquid=self.mass_flow_ox,
                m_dot_gas=self.mass_flow_fuel,
                rho_liquid=self.rho_ox,
                rho_gas=self.rho_fuel,
                delta_p_liquid=delta_p_ox,
                delta_p_gas=delta_p_fuel,
                n_elements=self.n_elements,
            )
        elif t_type in {"pintle", "pintle_injector"}:
            # Standard arrangement: Fuel annular (outer), Oxidizer radial (inner)
            res = size_pintle_injector(
                m_dot_annular=self.mass_flow_fuel,
                m_dot_radial=self.mass_flow_ox,
                rho_annular=self.rho_fuel,
                rho_radial=self.rho_ox,
                delta_p_annular=delta_p_fuel,
                delta_p_radial=delta_p_ox,
                pintle_diameter=self.pintle_diameter,
                cd_annular=self.cd_fuel,
                cd_radial=self.cd_ox,
                n_radial_slots=20,
            )
        elif t_type in {"impinging", "doublet", "impinging_doublet"}:
            res = size_impinging_doublet(
                m_dot_1=self.mass_flow_ox,
                m_dot_2=self.mass_flow_fuel,
                rho_1=self.rho_ox,
                rho_2=self.rho_fuel,
                delta_p_1=delta_p_ox,
                delta_p_2=delta_p_fuel,
                n_elements=self.n_elements,
                cd_1=self.cd_ox,
                cd_2=self.cd_fuel,
            )
        else:
            raise ValueError(
                f"Unsupported injector_type '{self.injector_type}'. "
                "Choose from the 4 canonical types: 'coaxial', 'swirl', 'pintle', or 'impinging'."
            )

        res["injector_type"] = t_type
        res["chamber_pressure_bar"] = self.chamber_pressure / 1e5
        res["delta_p_ratio"] = self.delta_p_ratio
        res["delta_p_bar"] = delta_p_ox / 1e5
        res["chugging_margin_adequate"] = self.delta_p_ratio >= 0.15
        return res
