"""
THE canonical Bartz gas-side film coefficient. One implementation.
==================================================================

This module replaces six independent implementations that disagreed by 58.9 %
on one common input. The audit established that the disagreement was almost
entirely CALL-SITE DEFAULTS, not rival equations: fed identical mu, Pr, sigma
and curvature, the algebraic forms collapse to a 7.8 % spread, and four of the
six agree exactly.

So the fix is not "choose the right equation". It is: one equation, and NO
SILENT DEFAULTS. Every property that used to be hardcoded at a call site
(sigma=0.7, Pr=0.5, mu=9e-5) is now either supplied by the caller or supplied
by a NAMED, RETURNED `Assumption` carrying `UNVALIDATED_ASSUMPTION`.

EVIDENCE STATE, STATED ONCE AND NOT PROMOTED
--------------------------------------------
Kryptonis level: **E1**.
  verification : DIMENSIONAL (PASS) + ALGEBRAIC (PASS)
  validation   : NONE -- never compared against a measured chamber here.

Bartz is E3/E4 in the literature. That is the literature's evidence. It does
not transfer. See docs/physics/bartz/BARTZ_CANONICAL_SPEC.md.

THE CURVATURE EXPONENT IS RESOLVED FROM PRIMARY EVIDENCE
---------------------------------------------------------
NASA SP-125 (Huzel & Huang), eq. 4-13 -- the document this repository cites
81 times -- legibly shows the curvature term (D_t/R)^0.1. The OCR fragment
"(Dt)O.!]" with "R = Radius of curvature of nozzle contour at throat" on the
following line reconstructs unambiguously as (Dt/R)^0.1. The single 0.2 reading
came from EUCASS-6291, which the Bartz spec flags as "transcription unreliable".
The page image was independently reviewed on 2026-09-06: printed p.100,
PDF page 109 (one-based), Eq.4-13. Printed p.101/PDF page110 also gives the
closed sigma expression explicitly as Eq.4-14. Numerical exponents are
unchanged; source-image hashes and mappings are recorded with the authority.

Primary-source resolution record: BARTZ_PRIMARY_SOURCE_RESOLUTION.json
"""

from __future__ import annotations

import math

from kryptonis.propulsion_equations.units import (
    Assumption, EvidenceLevel, NOT_REPORTED, Result, Status, Validation,
    Verification, molar_mass_g_per_mol, to_rankine, worst,
)

__all__ = [
    "BARTZ_COEFFICIENT", "CURVATURE_EXPONENT", "VISCOSITY_COEFF_SI",
    "bartz_prandtl", "bartz_viscosity", "bartz_sigma", "recovery_temperature",
    "bartz_film_coefficient", "SOURCE", "LOCATOR", "VALIDITY",
]

SOURCE = ("Bartz, D. R., 'A Simple Equation for Rapid Estimation of Rocket "
          "Nozzle Convective Heat Transfer Coefficients', Jet Propulsion "
          "27(1), Jan 1957, pp. 49-51")
LOCATOR = "Eq. [7], p.50 (main); Eq. [4] p.49 (sigma); Eq. [8] (Pr)"
SP125_SOURCE = "Huzel and Huang, NASA SP-125, Design of Liquid Propellant Rocket Engines, second edition"
SP125_SIGMA_LOCATOR = "printed p.101, PDF page110 (one-based), Eq.4-14; visually reviewed 2026-09-06"
VALIDITY = ("turbulent attached boundary layer AT THE THROAT; D_t/r_c <= 3; "
            "contraction and expansion angles varying by no more than 50%")

#: Dimensionless. PROVED dimensionless by tools/bartz_dimensions.py, so it is
#: unit-system independent and MUST NOT be converted. Bartz's original
#: (p_c*g/c*) carries g_c, the US-customary lbf<->lbm bridge, identically 1
#: in SI.
BARTZ_COEFFICIENT = 0.026

#: Confirmed by primary source (NASA SP-125 eq. 4-13, legible): 0.1.
CURVATURE_EXPONENT = 0.1

#: Bartz's own viscosity relation is mu = 46.6e-10 * M^0.5 * T[degR]^0.6 in
#: lb/(in.s). That coefficient IS dimensional. Converting:
#:     1 lb/(in.s) = 0.45359237/0.0254 = 17.857967 Pa.s
#:     (degR/K)^0.6 = 1.8^0.6         =  1.422864
#:     46.6e-10 * 17.857967 * 1.422864 = 1.184081e-7
#: The pre-existing orchestrator constant 1.184e-7 matches to 0.0068%. It was
#: never a magic number; the provenance had simply never been written down.
VISCOSITY_COEFF_US = 46.6e-10          # lb/(in.s), M in g/mol, T in degR
VISCOSITY_COEFF_SI = 1.1840806e-7      # Pa.s,      M in g/mol, T in K

_LB_PER_IN_S_TO_PA_S = 0.45359237 / 0.0254


# ---------------------------------------------------------------------------
# Companion relations from Bartz's own paper
# ---------------------------------------------------------------------------

def bartz_prandtl(gamma: float) -> Result:
    r"""Prandtl number from gamma. Bartz Eq. [8].

    .. math::  \mathrm{Pr} = \frac{4\gamma}{9\gamma - 5}

    At gamma = 1.2 this gives 0.8276. Two prior implementations hardcoded
    Pr = 0.5, which is +35.3 % on h_g and is not a combustion-gas Prandtl
    number; a third hardcoded 0.7, which is +10.6 %.
    """
    if not (1.0 < gamma <= 1.67):
        return Result(float("nan"), "-", "CAN-BARTZ-PR",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      source_locator="Eq. [8]",
                      notes=f"gamma={gamma} is outside (1, 1.67]; below 1 no "
                            f"gas exists and above 5/3 exceeds monatomic")
    denom = 9.0 * gamma - 5.0
    return Result(4.0 * gamma / denom, "-", "CAN-BARTZ-PR",
                  status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE, source=SOURCE,
                  source_locator="Eq. [8]",
                  validity_domain="combustion gas; Bartz's own relation",
                  uncertainty=NOT_REPORTED, inputs={"gamma": gamma})


def bartz_viscosity(molar_mass_kg_per_mol: float,
                    stagnation_temperature_K: float) -> Result:
    r"""Gas viscosity from Bartz's own relation, in SI.

    .. math::
        \mu = 46.6\times10^{-10}\,M^{0.5}\,T^{0.6}\ [\mathrm{lb/(in\cdot s)}],
        \quad T\ \mathrm{in}\ ^\circ R

    Converted once, explicitly, to SI. **M must be kg/mol and T must be K** --
    the conversion is done inside, so the caller never has to know that the
    original is in g/mol and rankine. That is the entire point of doing it
    here rather than at five call sites.
    """
    if molar_mass_kg_per_mol <= 0 or stagnation_temperature_K <= 0:
        return Result(float("nan"), "Pa.s", "CAN-BARTZ-VISC",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      source_locator="unnumbered viscosity relation",
                      notes="molar mass and temperature must be positive")
    # Guard the unit error that the SI constant silently depends on: a molar
    # mass in kg/mol passed as if it were g/mol is a 31.6x error in mu.
    m_g = molar_mass_g_per_mol(molar_mass_kg_per_mol)
    status = Status.PASS
    note = ""
    if not (1.0 <= m_g <= 200.0):
        status = Status.PHYSICALLY_INVALID
        note = (f"molar mass {m_g:.4g} g/mol is outside 1-200; the caller has "
                f"almost certainly passed g/mol where kg/mol was required")
    mu = VISCOSITY_COEFF_SI * m_g ** 0.5 * stagnation_temperature_K ** 0.6
    return Result(mu, "Pa.s", "CAN-BARTZ-VISC", status=status,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE, source=SOURCE,
                  source_locator="unnumbered viscosity relation, p.50",
                  validity_domain="rocket combustion gas at stagnation",
                  uncertainty=NOT_REPORTED, notes=note,
                  inputs={"M_kg_per_mol": molar_mass_kg_per_mol,
                          "T0_K": stagnation_temperature_K})


def bartz_viscosity_us(molar_mass_g_per_mol_: float,
                       temperature_rankine: float) -> float:
    """The ORIGINAL US-customary form, kept so the SI form can be tested
    against it rather than against itself. Returns lb/(in.s)."""
    return (VISCOSITY_COEFF_US * molar_mass_g_per_mol_ ** 0.5
            * temperature_rankine ** 0.6)


def bartz_sigma(*, wall_temperature_K: float, stagnation_temperature_K: float,
                gamma: float, mach: float) -> Result:
    r"""Property-variation correction across the boundary layer.

    .. math::
        \sigma = \left[\tfrac12\frac{T_{wg}}{T_0}
                 \left(1+\tfrac{\gamma-1}{2}M^2\right)+\tfrac12\right]^{-0.68}
                 \left[1+\tfrac{\gamma-1}{2}M^2\right]^{-0.12}

    **sigma MAY EXCEED 1 AND THAT IS PHYSICAL.** A cold wall raises the
    near-wall density and RAISES the film coefficient: at T_wg/T_0 = 0.2,
    sigma = 1.42. A prior implementation asserted ``0 < sigma <= 1`` and so
    raised ValueError at every realistic regeneratively cooled throat, which
    made the wall-temperature solver unable to run at all. Not clamped, not
    guarded -- returned.

    The closed T_wg/T_0 form is explicitly printed in NASA SP-125 Eq.4-14,
    p.101. This source mapping does not claim experimental validation or
    extend applicability to films, reacting boundary layers or arbitrary fluids.
    """
    if wall_temperature_K <= 0 or stagnation_temperature_K <= 0:
        return Result(float("nan"), "-", "CAN-BARTZ-SIGMA",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      source_locator="Eq. [4] (raw form)",
                      notes="temperatures must be positive")
    if mach < 0:
        return Result(float("nan"), "-", "CAN-BARTZ-SIGMA",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes="Mach number cannot be negative")
    if not (1.0 < gamma <= 1.67):
        return Result(float("nan"), "-", "CAN-BARTZ-SIGMA",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes=f"gamma={gamma} outside (1, 1.67]")
    if wall_temperature_K > stagnation_temperature_K:
        return Result(float("nan"), "-", "CAN-BARTZ-SIGMA",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes="the wall cannot be hotter than stagnation")

    rec = 1.0 + 0.5 * (gamma - 1.0) * mach * mach
    ratio = wall_temperature_K / stagnation_temperature_K
    sigma = ((0.5 * ratio * rec + 0.5) ** -0.68) * (rec ** -0.12)
    return Result(sigma, "-", "CAN-BARTZ-SIGMA", status=Status.PASS,
                  evidence_level=EvidenceLevel.E1,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE, source=SP125_SOURCE,
                  source_locator=SP125_SIGMA_LOCATOR,
                  validity_domain="property-variation correction used with SP-125 Eq.4-13; "
                                  "turbulent gas boundary layer, positive temperatures, "
                                  "0 < Tw/T0 <= 1; no film/recombination/radiation correction",
                  uncertainty=NOT_REPORTED,
                  inputs={"T_wg_K": wall_temperature_K,
                          "T0_K": stagnation_temperature_K,
                          "gamma": gamma, "mach": mach})


def recovery_temperature(*, stagnation_temperature_K: float, prandtl: float,
                         gamma: float, mach: float) -> Result:
    r"""Adiabatic wall (recovery) temperature for a turbulent layer.

    .. math::
        T_{aw} = T_0\,\frac{1 + r\frac{\gamma-1}{2}M^2}
                           {1 + \frac{\gamma-1}{2}M^2},
        \qquad r = \mathrm{Pr}^{1/3}

    The prior implementation used ``T_aw = T_0 * Pr^(1/3)``, which is not the
    recovery relation -- it is the recovery FACTOR applied directly to the
    stagnation temperature. At the throat (M=1, gamma=1.2, Pr=0.5) that gives
    2699 K where the correct value is 3336 K, understating the driving
    temperature difference by 25 % and therefore the heat flux by 25 %.
    """
    if stagnation_temperature_K <= 0:
        return Result(float("nan"), "K", "CAN-RECOVERY-T",
                      status=Status.PHYSICALLY_INVALID,
                      notes="stagnation temperature must be positive")
    if prandtl <= 0:
        return Result(float("nan"), "K", "CAN-RECOVERY-T",
                      status=Status.PHYSICALLY_INVALID,
                      notes="Prandtl number must be positive")
    if not (1.0 < gamma <= 1.67) or mach < 0:
        return Result(float("nan"), "K", "CAN-RECOVERY-T",
                      status=Status.PHYSICALLY_INVALID,
                      notes="gamma must be in (1, 1.67] and M >= 0")
    r = prandtl ** (1.0 / 3.0)
    half = 0.5 * (gamma - 1.0) * mach * mach
    t_aw = stagnation_temperature_K * (1.0 + r * half) / (1.0 + half)
    return Result(t_aw, "K", "CAN-RECOVERY-T", status=Status.PASS,
                  evidence_level=EvidenceLevel.E1,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="standard compressible boundary-layer theory; "
                         "turbulent recovery factor r = Pr^(1/3)",
                  source_locator="e.g. Schlichting, Boundary-Layer Theory; "
                                 "NOT_OPENED in this repository",
                  validity_domain="turbulent boundary layer",
                  uncertainty=NOT_REPORTED,
                  inputs={"T0_K": stagnation_temperature_K, "Pr": prandtl,
                          "gamma": gamma, "mach": mach})


# ---------------------------------------------------------------------------
# THE canonical film coefficient
# ---------------------------------------------------------------------------

def bartz_film_coefficient(
    *,
    throat_diameter_m: float,
    chamber_pressure_Pa: float,
    c_star_m_s: float,
    specific_heat_J_kgK: float,
    # --- properties: supply them, or accept a NAMED assumption -------------
    prandtl: float | None = None,
    viscosity_Pa_s: float | None = None,
    sigma: float | None = None,
    # --- geometry ----------------------------------------------------------
    throat_radius_curvature_m: float | None = None,
    area_ratio_local_over_throat: float = 1.0,
    # --- inputs needed to DERIVE the properties if not supplied ------------
    gamma: float | None = None,
    molar_mass_kg_per_mol: float | None = None,
    stagnation_temperature_K: float | None = None,
    wall_temperature_K: float | None = None,
    mach: float | None = None,
) -> Result:
    r"""Gas-side convective film coefficient. THE canonical implementation.

    .. math::
        h_g = \frac{0.026}{D_t^{0.2}}
              \left(\frac{\mu^{0.2} c_p}{\mathrm{Pr}^{0.6}}\right)
              \left(\frac{p_c}{c^*}\right)^{0.8}
              \left(\frac{D_t}{R_c}\right)^{0.1}
              \left(\frac{A_t}{A}\right)^{0.9}\sigma

    All SI. The coefficient 0.026 is dimensionless (proved), so no unit
    conversion is applied to it and none may be.

    NO SILENT DEFAULTS. If `prandtl`, `viscosity_Pa_s` or `sigma` is omitted,
    it is derived from Bartz's own relations where the inputs allow, and the
    derivation is recorded as a named `Assumption` and the result carries
    `UNVALIDATED_ASSUMPTION`. If it cannot be derived, the result is
    `INSUFFICIENT_EVIDENCE` -- never a hardcoded 0.7 or 0.5.
    """
    inputs = dict(D_t_m=throat_diameter_m, p_c_Pa=chamber_pressure_Pa,
                  c_star_m_s=c_star_m_s, cp_J_kgK=specific_heat_J_kgK,
                  A_over_At=area_ratio_local_over_throat)
    bad = []
    if throat_diameter_m <= 0:
        bad.append("throat diameter must be > 0")
    if chamber_pressure_Pa <= 0:
        bad.append("chamber pressure must be > 0")
    if c_star_m_s <= 0:
        bad.append("c* must be > 0")
    if specific_heat_J_kgK <= 0:
        bad.append("cp must be > 0")
    if area_ratio_local_over_throat < 1.0:
        bad.append("A/A_t must be >= 1; no area is smaller than the throat")
    if throat_radius_curvature_m is not None and throat_radius_curvature_m <= 0:
        bad.append("throat radius of curvature must be > 0")
    if bad:
        return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      source_locator=LOCATOR, notes="; ".join(bad),
                      inputs=inputs)

    assumptions: list[Assumption] = []
    statuses: list[Status] = [Status.PASS]

    # --- Prandtl -----------------------------------------------------------
    if prandtl is None:
        if gamma is None:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.INSUFFICIENT_EVIDENCE, source=SOURCE,
                          source_locator=LOCATOR, inputs=inputs,
                          notes="Prandtl number not supplied and gamma not "
                                "given, so Bartz Eq.[8] cannot supply it "
                                "either. A hardcoded 0.5 or 0.7 is not an "
                                "answer -- 0.5 is +35% on h_g.")
        pr_res = bartz_prandtl(gamma)
        if pr_res.status is Status.PHYSICALLY_INVALID:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.PHYSICALLY_INVALID, source=SOURCE,
                          notes=pr_res.notes, inputs=inputs)
        prandtl = float(pr_res)
        assumptions.append(Assumption(
            "prandtl", round(prandtl, 6), "-",
            "derived from Bartz Eq.[8] Pr = 4g/(9g-5); not supplied by caller",
            Status.UNVALIDATED_ASSUMPTION))
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
    if prandtl <= 0:
        return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes="Prandtl number must be > 0", inputs=inputs)
    if prandtl > 5.0:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)

    # --- viscosity ---------------------------------------------------------
    if viscosity_Pa_s is None:
        if molar_mass_kg_per_mol is None or stagnation_temperature_K is None:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.INSUFFICIENT_EVIDENCE, source=SOURCE,
                          source_locator=LOCATOR, inputs=inputs,
                          notes="viscosity not supplied and (M, T0) not given, "
                                "so Bartz's own viscosity relation cannot "
                                "supply it either.")
        mu_res = bartz_viscosity(molar_mass_kg_per_mol,
                                 stagnation_temperature_K)
        if mu_res.status is Status.PHYSICALLY_INVALID:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.PHYSICALLY_INVALID, source=SOURCE,
                          notes=mu_res.notes, inputs=inputs)
        viscosity_Pa_s = float(mu_res)
        assumptions.append(Assumption(
            "viscosity_Pa_s", f"{viscosity_Pa_s:.4e}", "Pa.s",
            "derived from Bartz's own viscosity relation; not supplied",
            Status.UNVALIDATED_ASSUMPTION))
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
    if viscosity_Pa_s <= 0:
        return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes="viscosity must be > 0", inputs=inputs)

    # --- sigma -------------------------------------------------------------
    if sigma is None:
        have = (wall_temperature_K is not None
                and stagnation_temperature_K is not None
                and gamma is not None and mach is not None)
        if not have:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.INSUFFICIENT_EVIDENCE, source=SOURCE,
                          source_locator=LOCATOR, inputs=inputs,
                          notes="sigma not supplied and (T_wg, T0, gamma, M) "
                                "not given. Two prior implementations "
                                "hardcoded sigma=0.7 where the canonical "
                                "value at a cooled throat is ~1.35 -- a "
                                "-48% error, the single largest term in the "
                                "58.9% cross-implementation spread.")
        s_res = bartz_sigma(wall_temperature_K=wall_temperature_K,
                            stagnation_temperature_K=stagnation_temperature_K,
                            gamma=gamma, mach=mach)
        if s_res.status is Status.PHYSICALLY_INVALID:
            return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                          status=Status.PHYSICALLY_INVALID, source=SOURCE,
                          notes=s_res.notes, inputs=inputs)
        sigma = float(s_res)
        assumptions.append(Assumption(
            "sigma", round(sigma, 6), "-",
            "computed from the Bartz property-variation correction",
            Status.PASS))
    if sigma <= 0:
        return Result(float("nan"), "W/(m^2.K)", "CAN-BARTZ-HG",
                      status=Status.PHYSICALLY_INVALID, source=SOURCE,
                      notes="sigma must be > 0. NOTE: sigma > 1 is PHYSICAL "
                            "at a cold wall and is NOT rejected.", inputs=inputs)

    # --- curvature ---------------------------------------------------------
    if throat_radius_curvature_m is None:
        curvature_term = 1.0
        dt_over_rc = 1.0
        assumptions.append(Assumption(
            "throat_radius_curvature_m", "omitted", "m",
            "curvature term set to 1 because R_c was not supplied; this is "
            "R_c = D_t, not 'no curvature effect'",
            Status.UNVALIDATED_ASSUMPTION))
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
    else:
        dt_over_rc = throat_diameter_m / throat_radius_curvature_m
        curvature_term = dt_over_rc ** CURVATURE_EXPONENT
        if dt_over_rc > 3.0:
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)

    # --- the correlation ---------------------------------------------------
    h_g = (
        (BARTZ_COEFFICIENT / throat_diameter_m ** 0.2)
        * (viscosity_Pa_s ** 0.2 * specific_heat_J_kgK / prandtl ** 0.6)
        * (chamber_pressure_Pa / c_star_m_s) ** 0.8
        * curvature_term
        * (1.0 / area_ratio_local_over_throat) ** 0.9
        * sigma
    )

    # Applied away from the throat, the coefficient's anchor no longer holds.
    if area_ratio_local_over_throat > 1.0 + 1e-12:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)

    notes = []
    if Status.OUT_OF_CORRELATION_RANGE in statuses:
        notes.append("outside the stated Bartz domain (throat, D_t/r_c <= 3)")

    return Result(
        h_g, "W/(m^2.K)", "CAN-BARTZ-HG", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source=SOURCE, source_locator=LOCATOR, validity_domain=VALIDITY,
        uncertainty="NOT_REPORTED by the primary. The commonly quoted "
                    "+/-20-30% is NOT traced to any source and must not be "
                    "cited.",
        assumptions=tuple(assumptions), notes="; ".join(notes),
        inputs={**inputs, "Pr": prandtl, "mu_Pa_s": viscosity_Pa_s,
                "sigma": sigma, "D_t_over_R_c": dt_over_rc},
    )
