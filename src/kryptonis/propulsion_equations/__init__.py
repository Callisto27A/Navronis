"""
Kryptonis Propulsion Equations: Authority-Controlled Liquid Rocket Analytical Engine
===================================================================================

A production-grade, mathematically verified library implementing the canonical
first-principles equations for liquid rocket propulsion, thrust chamber sizing,
gas dynamics, hot-gas convective heat transfer (Bartz), regenerative cooling,
wall conduction, and nozzle expansion losses.
"""

from kryptonis.propulsion_equations.units import (
    Result,
    Status,
    EvidenceLevel,
    Verification,
    Validation,
    Assumption,
    to_inch,
    to_cm,
    to_mm,
    to_rankine,
    kelvin_from_rankine,
)

from kryptonis.propulsion_equations.combustion import (
    assumed_c_star_efficiency,
    l_star_requirement,
    chamber_bulk_residence_time,
    sp125_stay_time,
)

from kryptonis.propulsion_equations.chamber import (
    throat_area,
    throat_diameter,
    chamber_diameter,
    chamber_volume,
    convergent_volume,
    cylinder_length,
    cylindrical_length,
    vandenkerckhove,
    c_star_ideal,
)

from kryptonis.propulsion_equations.aerodynamics import (
    calculate_area_ratio,
    calculate_mach_from_area_ratio,
)

from kryptonis.propulsion_equations.contour import (
    calculate_divergence_loss_factor,
    calculate_bell_divergence_loss_factor,
    calculate_displacement_thickness,
)

from kryptonis.propulsion_equations.bartz import (
    bartz_film_coefficient,
    bartz_sigma,
    bartz_viscosity,
    bartz_prandtl,
    recovery_temperature,
)

from kryptonis.propulsion_equations.thermal import (
    CoolantCorrelation,
    nusselt_mcadams,
    nusselt_dittus_boelter_1930,
    nusselt_sieder_tate,
    nusselt_gnielinski,
    nusselt_taylor_tn_d4332,
    ito_curvature_factor,
    friction_factor_haaland,
    friction_factor_colebrook,
    solve_wall_temperature,
    axial_thermal_march,
    hot_spot,
)

from kryptonis.propulsion_equations.wall_conduction import (
    conduction_heat_flux,
    cylindrical_wall_correction,
    fourier_wall_temperature_drop,
)

from kryptonis.propulsion_equations.nozzle_losses import (
    divergence_efficiency,
    separation_assessment,
    delivered_thrust_coefficient,
)

from kryptonis.propulsion_equations.chamber_acoustics import (
    first_tangential_frequency,
    first_radial_frequency,
    first_longitudinal_frequency,
)

from kryptonis.propulsion_equations.combustor import (
    CombustorDesign,
    CombustorResult,
)

__version__ = "0.1.0"
__all__ = [
    # Units & result containers
    "Result", "Status", "EvidenceLevel", "Verification", "Validation", "Assumption",
    "to_inch", "to_cm", "to_mm", "to_rankine", "kelvin_from_rankine",
    # Combustion
    "assumed_c_star_efficiency", "l_star_requirement", "chamber_bulk_residence_time", "sp125_stay_time",
    # Chamber sizing
    "throat_area", "throat_diameter", "chamber_diameter", "chamber_volume",
    "convergent_volume", "cylinder_length", "cylindrical_length", "vandenkerckhove", "c_star_ideal",
    # Gas dynamics & Aerodynamics
    "calculate_area_ratio", "calculate_mach_from_area_ratio",
    # Contour & losses
    "calculate_divergence_loss_factor", "calculate_bell_divergence_loss_factor", "calculate_displacement_thickness",
    # Bartz heat transfer
    "bartz_film_coefficient", "bartz_sigma", "bartz_viscosity", "bartz_prandtl", "recovery_temperature",
    # Thermal & Coolant correlations
    "CoolantCorrelation", "nusselt_mcadams", "nusselt_dittus_boelter_1930",
    "nusselt_sieder_tate", "nusselt_gnielinski", "nusselt_taylor_tn_d4332",
    "ito_curvature_factor", "friction_factor_haaland", "friction_factor_colebrook",
    "solve_wall_temperature", "axial_thermal_march", "hot_spot",
    # Wall conduction
    "conduction_heat_flux", "cylindrical_wall_correction", "fourier_wall_temperature_drop",
    # Nozzle losses
    "divergence_efficiency", "separation_assessment", "delivered_thrust_coefficient",
    # Acoustics
    "first_tangential_frequency", "first_radial_frequency", "first_longitudinal_frequency",
    # High-level Combustor API (Day 1 Release)
    "CombustorDesign", "CombustorResult",
]
