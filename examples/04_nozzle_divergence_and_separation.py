"""
Example 04: Nozzle Supersonic Expansion, Divergence Loss & Separation
=====================================================================
Calculates:
- Supersonic Mach number and static pressure along the nozzle expansion contour
- Divergence loss efficiency factor (lambda) for conical vs. Rao bell nozzles
- Flow separation margin at sea-level ambient conditions using Summerfield and
  Schmucker published criteria
"""

import math
from kryptonis.propulsion_equations.aerodynamics import (
    calculate_area_ratio,
    calculate_mach_from_area_ratio,
)
from kryptonis.propulsion_equations.contour import (
    calculate_divergence_loss_factor,
    calculate_bell_divergence_loss_factor,
)
from kryptonis.propulsion_equations.nozzle_losses import separation_assessment

def main():
    print("=" * 65)
    print("Kryptonis: Nozzle Expansion, Divergence Loss & Separation")
    print("=" * 65)

    # 1. Nozzle Inputs
    gamma = 1.2
    P_c = 7.0e6           # 7.0 MPa chamber pressure
    eps_exit = 35.0       # Exit expansion area ratio Ae / At
    P_ambient_sl = 101325.0  # 1 atm sea-level ambient pressure

    print(f"Combustion Gamma:         {gamma:.2f}")
    print(f"Chamber Pressure (Pc):    {P_c / 1e6:.1f} MPa ({P_c / 1e5:.0f} bar)")
    print(f"Expansion Ratio (eps):    {eps_exit:.1f}")

    # 2. Solve Exit Mach from Area-Mach Relation
    mach_exit = calculate_mach_from_area_ratio(
        area_ratio=eps_exit,
        gamma=gamma,
        supersonic=True
    )
    print(f"Solved Exit Mach Number:  {mach_exit:.3f}")

    # 3. Static Exit Pressure P_e = P_c / (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))
    pressure_ratio = (1.0 + 0.5 * (gamma - 1.0) * (mach_exit ** 2)) ** (gamma / (gamma - 1.0))
    P_exit = P_c / pressure_ratio
    print(f"Static Exit Pressure (Pe):{P_exit / 1e3:.2f} kPa ({P_exit / 1e5:.3f} bar)")
    print(f"Exit Pressure Ratio Pe/Pa:{P_exit / P_ambient_sl:.3f}")

    # 4. Divergence Loss Comparison
    # Conical nozzle with 15 deg half-angle
    lambda_cone = calculate_divergence_loss_factor(15.0)
    # Rao bell with 28 deg initial angle and 8 deg exit angle
    lambda_bell = calculate_bell_divergence_loss_factor(theta_n_deg=28.0, theta_e_deg=8.0)

    print(f"\n--- Divergence Efficiency ---")
    print(f"15-deg Conical Lambda:    {lambda_cone:.4f} (Divergence loss: {(1 - lambda_cone) * 100:.2f}%)")
    print(f"Rao Bell Contour Lambda:  {lambda_bell:.4f} (Divergence loss: {(1 - lambda_bell) * 100:.2f}%)")
    print(f"Efficiency Gain with Bell:+{(lambda_bell - lambda_cone) * 100:.2f}% delivered Isp")

    # 5. Sea-Level Flow Separation Assessment
    print(f"\n--- Sea-Level Flow Separation Assessment ---")
    res_sep = separation_assessment(
        exit_pressure_Pa=P_exit,
        ambient_pressure_Pa=P_ambient_sl,
        exit_mach=mach_exit,
        chamber_pressure_Pa=P_c
    )
    print(f"Separation Status:        {res_sep.status.value}")
    if "Summerfield_1954" in res_sep.inputs:
        p_sep_summerfield = res_sep.inputs["Summerfield_1954"]["p_sep_Pa"]
        print(f"Summerfield Limit P_sep:  {p_sep_summerfield / 1e3:.1f} kPa")
        if P_exit < p_sep_summerfield:
            print("Verdict: Pe < P_sep -> Flow will SEPARATE at sea level! (Requires altitude start or smaller eps)")
        else:
            print("Verdict: Pe > P_sep -> Flow remains safely attached at sea level.")
    print("=" * 65)

if __name__ == "__main__":
    main()
