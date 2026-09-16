"""
Regenerative Cooling Channels & 1D Conjugate Thermal Analysis (Component 3)
=============================================================================
Provides first-principles analytical sizing and thermal-hydraulic evaluation
for liquid rocket engine milled cooling jackets.

Governing Physics:
1. Milled Channel Geometry:
   - Channel width: w_c = (pi * D - N * t_fin) / N
   - Aspect ratio: AR = h_c / w_c
   - Hydraulic diameter: D_h = 2 * w_c * h_c / (w_c + h_c)
   - Flow area per channel: A_flow = w_c * h_c

2. Coolant Hydraulics & Convection:
   - Coolant velocity: v_c = (m_dot_coolant / N) / (rho_c * A_flow)
   - Reynolds number: Re = rho_c * v_c * D_h / mu_c
   - Prandtl number: Pr = Cp_c * mu_c / k_c
   - Friction factor: Haaland (1983) correlation (explicit Colebrook-White)
   - Heat transfer: Gnielinski (1976) correlation with Filonenko friction factor
   - Fin effectiveness: eta_fin = tanh(m * h_c) / (m * h_c), m = sqrt(2 * h_c / (k_w * t_fin))
   - Enhanced HTC: h_c,eff = h_c * (w_c + 2 * eta_fin * h_c) / (w_c + t_fin)

3. 1D Conjugate Wall Thermal Closure:
   - Gas-side recovery temperature (Bartz 1957): T_aw
   - Cylindrical/Planar wall conduction resistance: R_wall = t_w / k_w
   - Three-resistance equilibrium:
     q = (T_aw - T_bulk) / (1/h_g + R_wall + 1/h_c,eff)
     T_wg = T_aw - q / h_g
     T_wc = T_bulk + q / h_c,eff

4. Thermal-Structural Margins:
   - Darcy-Weisbach channel pressure drop: delta_P = f * (L / D_h) * (rho * v^2 / 2)
   - Coolant bulk temperature rise: delta_T_bulk = Q_total / (m_dot * Cp)
   - 1D thermal compressive stress: sigma_th = E * alpha * (T_wg - T_wc) / (2 * (1 - nu))
   - Yield margin of safety: MS = (S_y / sigma_th) - 1.0

Primary Literature Provenance:
- NASA SP-8087 (1974): Liquid Rocket Engine Fluid-Cooled Combustion Chambers
- Huzel & Huang (NASA SP-125, 1992): Modern Engineering for Design of LRE
- Gnielinski, V. (1976): New Equations for Heat and Mass Transfer in Turbulent Pipe and Channel Flow
- Haaland, S. E. (1983): Simple and Explicit Formulas for the Friction Factor in Turbulent Flow
- Bartz, D. R. (1957): A Simple Equation for Rapid Estimation of Rocket Nozzle HTC
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional


@dataclass
class RegenChannelResult:
    """Immutable container for regenerative cooling jacket sizing outputs."""

    # Channel geometry
    n_channels: int
    channel_width_mm: float
    channel_height_mm: float
    fin_thickness_mm: float
    aspect_ratio: float
    hydraulic_diameter_mm: float
    wall_thickness_mm: float

    # Coolant flow & hydraulics
    coolant_velocity_m_s: float
    reynolds_number: float
    prandtl_number: float
    friction_factor: float
    coolant_htc_W_m2K: float
    fin_efficiency: float
    enhanced_coolant_htc_W_m2K: float

    # Conjugate thermal solution
    gas_recovery_temp_K: float
    gas_htc_W_m2K: float
    hot_gas_wall_temp_K: float
    coolant_wall_temp_K: float
    peak_heat_flux_MW_m2: float

    # Hydraulic & thermal margins
    coolant_pressure_drop_bar: float
    coolant_temp_rise_K: float
    thermal_stress_MPa: float
    yield_safety_margin: float
    boiling_margin_adequate: bool

    # Meta
    liner_material: str
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        raise KeyError(f"Unknown result attribute: {key}")

    def explain(self, quantity_name: str) -> str:
        """Return the mathematical equation and primary literature source for any quantity."""
        explanations = {
            "channel_width_mm": (
                "w_c = (pi * D - N * t_fin) / N",
                "NASA SP-8087 / Huzel & Huang (NASA SP-125 Chapter 4)",
                "Circumferential milled channel spacing allowing for structural rib fins.",
            ),
            "hydraulic_diameter_mm": (
                "D_h = 4 * A_flow / P_wetted = 2 * w_c * h_c / (w_c + h_c)",
                "Standard internal duct hydraulics",
                "Equivalent hydraulic diameter for rectangular cooling passages.",
            ),
            "friction_factor": (
                "1/sqrt(f) = -1.8 * log10(((eps/D_h)/3.7)^1.11 + 6.9/Re)",
                "Haaland, S. E. (1983), J. Fluids Eng. 105(1)",
                "Explicit approximation to the implicit Colebrook-White equation with surface roughness.",
            ),
            "coolant_htc_W_m2K": (
                "Nu = ((f/8)*(Re - 1000)*Pr) / (1 + 12.7*sqrt(f/8)*(Pr^(2/3) - 1)), h_c = Nu*k/D_h",
                "Gnielinski, V. (1976), Int. Chem. Eng. 16(2) / Filonenko (1954)",
                "Canonical turbulent heat transfer correlation for forced convection in ducts.",
            ),
            "fin_efficiency": (
                "eta_fin = tanh(m * h_c) / (m * h_c), where m = sqrt(2 * h_c / (k_w * t_fin))",
                "Incropera & DeWitt, Fundamentals of Heat and Mass Transfer",
                "1D longitudinal fin heat transfer enhancement accounting for temperature drop along rib.",
            ),
            "hot_gas_wall_temp_K": (
                "T_wg = T_aw - q / h_g, where q = (T_aw - T_bulk) / (1/h_g + t_w/k_w + 1/h_c,eff)",
                "Fourier 1D Conduction + Bartz (1957) + Gnielinski (1976)",
                "Closed-form 3-resistance conjugate energy balance. Zero arbitrary assertions.",
            ),
            "coolant_pressure_drop_bar": (
                "delta_P = f * (L_channel / D_h) * (rho * v_c^2 / 2)",
                "Darcy-Weisbach equation (NASA SP-8087 Section 3)",
                "Frictional hydraulic head loss across the regenerative jacket.",
            ),
            "thermal_stress_MPa": (
                "sigma_th = E * alpha * (T_wg - T_wc) / (2 * (1 - nu))",
                "Dieter, G. E., Mechanical Metallurgy / Manson-Coffin LCF",
                "Restrained thermal expansion compressive stress across liner thickness.",
            ),
        }

        if quantity_name not in explanations:
            return f"Quantity '{quantity_name}' has no registered provenance entry."

        eq, src, note = explanations[quantity_name]
        val = getattr(self, quantity_name, "N/A")
        return (
            f"Quantity:     {quantity_name}\n"
            f"Value:        {val}\n"
            f"Equation:     {eq}\n"
            f"Source:       {src}\n"
            f"Description:  {note}\n"
        )


class RegenCoolingJacket:
    """Analytical sizer and conjugate heat transfer evaluator for regenerative cooling jackets.

    Attributes:
        throat_diameter_m: Sonic throat inner diameter (m)
        chamber_diameter_m: Combustor barrel inner diameter (m)
        chamber_length_m: Total combustor length (m)
        mass_flow_coolant_kg_s: Total coolant mass flow rate (kg/s)
        chamber_pressure_pa: Combustion chamber pressure (Pa)
        gas_recovery_temp_k: Adiabatic flame recovery temperature (K)
        gas_throat_htc_w_m2k: Hot-gas Bartz heat transfer coefficient at throat (W/m²-K)
        n_channels: Number of circumferential channels (default: 80)
        channel_height_m: Channel milled depth (m, default: 1.8 mm)
        fin_thickness_m: Channel rib fin width (m, default: 0.8 mm)
        wall_thickness_m: Hot-gas liner wall thickness (m, default: 1.5 mm)
        coolant_type: Coolant fluid name ('CH4', 'RP-1', 'LH2')
        coolant_inlet_temp_k: Bulk coolant inlet temperature (K)
        liner_material: Liner alloy name ('CuCrZr', 'GRCop-42', 'OFHC', 'Inconel-718')
        surface_roughness_m: Equivalent sand-grain surface roughness (default: 15 um)
    """

    COOLANT_PROPERTIES = {
        "CH4": {
            "rho": 422.0,       # kg/m³
            "cp": 3500.0,       # J/kg-K
            "k": 0.18,          # W/m-K
            "mu": 1.1e-4,       # Pa-s
            "t_inlet_default": 120.0,
            "t_limit": 650.0,   # Coking / pyrolysis threshold
        },
        "RP-1": {
            "rho": 810.0,
            "cp": 2050.0,
            "k": 0.12,
            "mu": 1.5e-3,
            "t_inlet_default": 290.0,
            "t_limit": 560.0,   # Coking limit
        },
        "LH2": {
            "rho": 71.0,
            "cp": 14200.0,
            "k": 0.19,
            "mu": 1.4e-5,
            "t_inlet_default": 40.0,
            "t_limit": 300.0,
        },
        "MMH": {
            "rho": 880.0,       # kg/m³
            "cp": 2900.0,       # J/kg-K
            "k": 0.20,          # W/m-K
            "mu": 8.5e-4,       # Pa-s
            "t_inlet_default": 298.0,
            "t_limit": 450.0,   # Thermal decomposition threshold
        },
        "HYDRAZINE": {
            "rho": 1004.0,      # kg/m³
            "cp": 3080.0,       # J/kg-K
            "k": 0.25,          # W/m-K
            "mu": 9.7e-4,       # Pa-s
            "t_inlet_default": 298.0,
            "t_limit": 420.0,   # Exothermic decomposition threshold
        },
        "ETHANOL": {
            "rho": 789.0,       # kg/m³
            "cp": 2440.0,       # J/kg-K
            "k": 0.17,          # W/m-K
            "mu": 1.2e-3,       # Pa-s
            "t_inlet_default": 293.0,
            "t_limit": 400.0,   # Vaporization / boiling limit
        },
    }

    MATERIAL_PROPERTIES = {
        "CuCrZr": {
            "k_w": 320.0,       # W/m-K
            "e_mod": 120.0e9,   # Pa
            "alpha": 17.0e-6,   # 1/K
            "nu": 0.33,
            "yield_strength": 280.0e6, # Pa at elevated temp
            "t_max_service": 850.0,    # K
        },
        "GRCop-42": {
            "k_w": 340.0,
            "e_mod": 130.0e9,
            "alpha": 17.5e-6,
            "nu": 0.33,
            "yield_strength": 310.0e6,
            "t_max_service": 950.0,
        },
        "OFHC": {
            "k_w": 380.0,
            "e_mod": 115.0e9,
            "alpha": 16.5e-6,
            "nu": 0.34,
            "yield_strength": 180.0e6,
            "t_max_service": 650.0,
        },
        "Inconel-718": {
            "k_w": 15.0,
            "e_mod": 200.0e9,
            "alpha": 13.0e-6,
            "nu": 0.29,
            "yield_strength": 900.0e6,
            "t_max_service": 1200.0,
        },
    }

    def __init__(
        self,
        *,
        throat_diameter_m: float,
        chamber_diameter_m: float,
        chamber_length_m: float,
        mass_flow_coolant_kg_s: float,
        chamber_pressure_pa: float,
        gas_recovery_temp_k: float,
        gas_throat_htc_w_m2k: float,
        n_channels: int = 80,
        channel_height_m: float = 0.0018,
        fin_thickness_m: float = 0.0008,
        wall_thickness_m: float = 0.0015,
        coolant_type: str = "CH4",
        coolant_inlet_temp_k: Optional[float] = None,
        liner_material: str = "CuCrZr",
        surface_roughness_m: float = 15.0e-6,
    ) -> None:
        if throat_diameter_m <= 0:
            raise ValueError(f"Throat diameter must be positive, got {throat_diameter_m}")
        if chamber_diameter_m <= throat_diameter_m:
            raise ValueError(f"Chamber diameter ({chamber_diameter_m}) must exceed throat ({throat_diameter_m})")
        if mass_flow_coolant_kg_s <= 0:
            raise ValueError(f"Coolant mass flow must be positive, got {mass_flow_coolant_kg_s}")
        if n_channels < 10:
            raise ValueError(f"Channel count must be at least 10, got {n_channels}")
        if fin_thickness_m <= 0 or channel_height_m <= 0 or wall_thickness_m <= 0:
            raise ValueError("Channel dimensions must be strictly positive")

        self.dt = float(throat_diameter_m)
        self.dc = float(chamber_diameter_m)
        self.lc = float(chamber_length_m)
        self.mdot_c = float(mass_flow_coolant_kg_s)
        self.pc = float(chamber_pressure_pa)
        self.t_aw = float(gas_recovery_temp_k)
        self.h_g = float(gas_throat_htc_w_m2k)
        self.n_ch = int(n_channels)
        self.h_c = float(channel_height_m)
        self.t_fin = float(fin_thickness_m)
        self.t_w = float(wall_thickness_m)
        self.coolant_type = (
            coolant_type.upper()
            .replace("LCH4", "CH4")
            .replace("METHANE", "CH4")
            .replace("KEROSENE", "RP-1")
            .replace("N2H4", "HYDRAZINE")
        )
        self.liner_material = liner_material
        self.eps = float(surface_roughness_m)

        c_props = self.COOLANT_PROPERTIES.get(self.coolant_type, self.COOLANT_PROPERTIES["CH4"])
        self.t_bulk = coolant_inlet_temp_k if coolant_inlet_temp_k is not None else c_props["t_inlet_default"]
        self.rho_c = c_props["rho"]
        self.cp_c = c_props["cp"]
        self.k_c = c_props["k"]
        self.mu_c = c_props["mu"]
        self.t_c_limit = c_props["t_limit"]

        m_props = self.MATERIAL_PROPERTIES.get(self.liner_material, self.MATERIAL_PROPERTIES["CuCrZr"])
        self.k_w = m_props["k_w"]
        self.e_mod = m_props["e_mod"]
        self.alpha_w = m_props["alpha"]
        self.nu_w = m_props["nu"]
        self.s_y = m_props["yield_strength"]
        self.t_w_max = m_props["t_max_service"]

    def solve(self) -> RegenChannelResult:
        """Solve regenerative cooling jacket hydraulics and conjugate wall thermal closure."""
        # 1. Throat station channel geometry
        circ_throat = math.pi * self.dt
        total_fin_width = self.n_ch * self.t_fin
        if total_fin_width >= circ_throat:
            raise ValueError(
                f"Channels cannot physically fit on throat perimeter! "
                f"N * t_fin = {total_fin_width*1000:.1f} mm exceeds pi * Dt = {circ_throat*1000:.1f} mm"
            )

        w_c = (circ_throat - total_fin_width) / self.n_ch
        aspect_ratio = self.h_c / w_c
        a_flow_per_ch = w_c * self.h_c
        p_wetted = 2.0 * (w_c + self.h_c)
        d_h = 4.0 * a_flow_per_ch / p_wetted

        # 2. Coolant hydraulics
        m_dot_per_ch = self.mdot_c / self.n_ch
        v_c = m_dot_per_ch / (self.rho_c * a_flow_per_ch)
        re = self.rho_c * v_c * d_h / self.mu_c
        pr = self.cp_c * self.mu_c / self.k_c

        # Friction factor via Haaland (1983)
        inv_sqrt_f = -1.8 * math.log10(((self.eps / d_h) / 3.7) ** 1.11 + 6.9 / max(re, 1000.0))
        f = (1.0 / inv_sqrt_f) ** 2

        # Gnielinski (1976) forced convection
        f_filon = (0.790 * math.log(max(re, 2300.0)) - 1.64) ** -2
        nu = ((f_filon / 8.0) * (re - 1000.0) * pr) / (
            1.0 + 12.7 * math.sqrt(f_filon / 8.0) * (pr ** (2.0 / 3.0) - 1.0)
        )
        h_c_base = max(nu * self.k_c / d_h, 1000.0)

        # 3. Fin efficiency & enhanced HTC
        m_fin = math.sqrt(2.0 * h_c_base / (self.k_w * self.t_fin))
        eta_fin = math.tanh(m_fin * self.h_c) / (m_fin * self.h_c)
        a_eff = w_c + 2.0 * eta_fin * self.h_c
        area_enhancement = a_eff / (w_c + self.t_fin)
        h_c_eff = h_c_base * area_enhancement

        # 4. 1D Conjugate Wall Thermal Energy Balance
        # Three thermal resistances in series:
        # R_gas = 1/h_g, R_wall = t_w/k_w, R_coolant = 1/h_c_eff
        r_wall = self.t_w / self.k_w
        r_total = (1.0 / self.h_g) + r_wall + (1.0 / h_c_eff)

        q = (self.t_aw - self.t_bulk) / r_total
        t_wg = self.t_aw - (q / self.h_g)
        t_wc = self.t_bulk + (q / h_c_eff)

        # 5. Pressure drop across chamber and throat jacket
        l_channel = self.lc * 1.15  # 15% curvature/manifold path elongation
        delta_p_pa = f * (l_channel / d_h) * 0.5 * self.rho_c * (v_c ** 2)
        delta_p_bar = delta_p_pa / 1.0e5

        # 6. Coolant bulk temperature rise
        wetted_area_chamber = math.pi * self.dt * self.lc * 0.65  # Integral profile area factor
        total_q_dot = q * wetted_area_chamber
        delta_t_bulk = total_q_dot / (self.mdot_c * self.cp_c)

        # 7. Thermal stress and safety margin to yield
        delta_t_wall = t_wg - t_wc
        sigma_th_pa = (self.e_mod * self.alpha_w * delta_t_wall) / (2.0 * (1.0 - self.nu_w))
        sigma_th_mpa = sigma_th_pa / 1.0e6
        margin_sy = (self.s_y / sigma_th_pa) - 1.0

        boiling_margin = t_wc < self.t_c_limit and t_wg < self.t_w_max

        return RegenChannelResult(
            n_channels=self.n_ch,
            channel_width_mm=round(w_c * 1000.0, 3),
            channel_height_mm=round(self.h_c * 1000.0, 3),
            fin_thickness_mm=round(self.t_fin * 1000.0, 3),
            aspect_ratio=round(aspect_ratio, 2),
            hydraulic_diameter_mm=round(d_h * 1000.0, 3),
            wall_thickness_mm=round(self.t_w * 1000.0, 3),
            coolant_velocity_m_s=round(v_c, 1),
            reynolds_number=round(re, 0),
            prandtl_number=round(pr, 2),
            friction_factor=round(f, 4),
            coolant_htc_W_m2K=round(h_c_base, 1),
            fin_efficiency=round(eta_fin, 3),
            enhanced_coolant_htc_W_m2K=round(h_c_eff, 1),
            gas_recovery_temp_K=round(self.t_aw, 1),
            gas_htc_W_m2K=round(self.h_g, 1),
            hot_gas_wall_temp_K=round(t_wg, 1),
            coolant_wall_temp_K=round(t_wc, 1),
            peak_heat_flux_MW_m2=round(q / 1.0e6, 2),
            coolant_pressure_drop_bar=round(delta_p_bar, 2),
            coolant_temp_rise_K=round(delta_t_bulk, 1),
            thermal_stress_MPa=round(sigma_th_mpa, 1),
            yield_safety_margin=round(margin_sy, 2),
            boiling_margin_adequate=boiling_margin,
            liner_material=self.liner_material,
            metadata={
                "coolant_type": self.coolant_type,
                "reynolds": re,
                "aspect_ratio": aspect_ratio,
                "roughness_um": self.eps * 1e6,
            },
        )
