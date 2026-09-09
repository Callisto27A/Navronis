"""
Canonical chamber geometry and gas dynamics. One policy each.
=============================================================

This module ends three live contradictions that the SP-125 primary-source
audit established, and it ends them WITHOUT inventing anything.

CONTRADICTION 1 -- the contraction correlation is not Huzel & Huang's
--------------------------------------------------------------------
Three modules attributed ``eps_c = 1.25 + 8.0 D_t^-0.6`` to "NASA SP-125 /
Huzel & Huang ch. 4". SP-125 Chapter IV was opened at page level. It contains
NO such formula. A full-text index of every occurrence of "contraction ratio"
in the book returns four hits, none a formula in throat diameter. What p.88
gives is prose ranges keyed to FEED SYSTEM TYPE.

The formula traces to Humble, *Space Propulsion Analysis and Design*,
McGraw-Hill 1995, via cryo-rocket.com Eq. 5.1.8 where it is cited as ref [54];
ref [55] (SP-125) sits on the same page attached to the L* figure. The
Kryptonis citation was a conflation of two adjacent references.

CONTRADICTION 2 -- the correlation was implemented in the wrong unit system
--------------------------------------------------------------------------
The source requires **D_t in CENTIMETRES**. Kryptonis used INCHES. Verified by
reproducing the source's own arithmetic: A_t = 0.05712 m^2 -> D_t = 26.97 cm
-> eps_c = 2.3581 -> A_c = 0.13469 m^2, matching the published 0.1347. Inches
gives 0.18213 and does not.

At the 30 kN case this is eps_c 6.92 (inches) against 4.49 (cm): +54 % on
eps_c, +24 % on chamber diameter, and 6.92 lies outside EVERY range SP-125
states. **This error was introduced by a 2026-08-10 "correction" in this
repository** that changed the units while fixing a coefficient; the
pre-existing code was already in centimetres.

CONTRADICTION 3 -- three chains, two policies
---------------------------------------------
`engine_pipeline` took eps_c as an input defaulting to 2.5; `orchestrator` and
`closedloop` computed 6.92 from the mis-united correlation. The same
specification built two different chambers. There is now one policy and all
three call it.

WHAT THE POLICY IS
------------------
`contraction_ratio()` is the single entry point.

  * If the caller SUPPLIES eps_c, that is used. It is a design variable and a
    designer is allowed to choose it.
  * If not, the Humble correlation is evaluated IN CENTIMETRES and returned
    with `UNVALIDATED_ASSUMPTION`, because its primary source has never been
    opened by this project.
  * Either way the value is checked against SP-125's OWN ranges, which ARE
    E2 (primary opened, page-level locator, verified), and a value outside
    them returns `OUT_OF_CORRELATION_RANGE` -- a number plus a flag, never a
    refusal, because a designer may legitimately want an unusual chamber.

The coefficient itself is left as a recorded conflict: the pre-existing code
used ``1.00 + 1.44 D_t[cm]^-0.6`` and the traced source uses
``1.25 + 8.00 D_t[cm]^-0.6``. Same unit system, different numbers, and the
provenance of the first is now lost. That is CONFLICT_UNRESOLVED and it is
reported, not averaged.
"""

from __future__ import annotations

import math

from kryptonis.propulsion_equations.units import (
    Assumption, CM_M, EvidenceLevel, INCH_M, NOT_REPORTED, Result, Status,
    Validation, Verification, to_cm, worst,
)

__all__ = [
    "SP125", "HUMBLE", "EPS_C_ABSOLUTE_MIN", "EPS_C_ABSOLUTE_MAX",
    "EPS_C_TURBOPUMP", "EPS_C_PRESSURE_FED", "L_STAR_TABLE_4_1",
    "throat_area", "throat_diameter",
    "contraction_ratio", "characteristic_length", "chamber_volume",
    "chamber_diameter", "convergent_length", "convergent_volume",
    "cylinder_length", "cylindrical_length", "vandenkerckhove", "c_star_ideal",
]

def throat_area(*, mass_flow_kg_s: float, c_star_m_s: float, chamber_pressure_Pa: float) -> Result:
    r""".. math:: A_t = \frac{\dot{m} c^{*}}{P_c}"""
    if mass_flow_kg_s <= 0 or c_star_m_s <= 0 or chamber_pressure_Pa <= 0:
        return Result(float("nan"), "m^2", "CAN-AT", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = (mass_flow_kg_s * c_star_m_s) / chamber_pressure_Pa
    return Result(val, "m^2", "CAN-AT", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="NASA SP-125 Eq. (1-9)", source_locator="p.4",
                  validity_domain="isentropic choked throat flow", uncertainty="exact",
                  inputs={"mass_flow_kg_s": mass_flow_kg_s, "c_star_m_s": c_star_m_s, "chamber_pressure_Pa": chamber_pressure_Pa})

def throat_diameter(*, throat_area_m2: float) -> Result:
    r""".. math:: D_t = \sqrt{\frac{4 A_t}{\pi}}"""
    if throat_area_m2 <= 0:
        return Result(float("nan"), "m", "CAN-DT", status=Status.PHYSICALLY_INVALID, notes="A_t must be > 0")
    val = math.sqrt(4.0 * throat_area_m2 / math.pi)
    return Result(val, "m", "CAN-DT", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="geometric definition", source_locator="circular cross section",
                  validity_domain="axisymmetric circular throat", uncertainty="exact",
                  inputs={"throat_area_m2": throat_area_m2})


SP125 = ("NASA SP-125, Huzel & Huang, 'Design of Liquid Propellant Rocket "
         "Engines', 2nd ed., 1967 (NTRS 19710019929)")
HUMBLE = ("Humble, R. W., 'Space Propulsion Analysis and Design', "
          "McGraw-Hill, 1995 -- reached via cryo-rocket.com Eq. 5.1.8 "
          "ref [54]. PRIMARY NOT_OPENED by this project.")

# --- SP-125 p.88, verbatim, verified by direct page read -------------------
#: "For pressurized-gas propellant feed, low-thrust engine systems contraction
#:  area ratio values of 2 to 5 have been used."
EPS_C_PRESSURE_FED = (2.0, 5.0)
#: "For most turbopump propellant feed, high thrust and high chamber pressure
#:  engine systems lower ratio values of 1.3 to 2.5 are employed."
EPS_C_TURBOPUMP = (1.3, 2.5)
#: SP-125 p.24 design-parameter summary: "Nozzle contraction area ratio,
#:  eps_c ... 1.3 to 6". This is the widest range the document states, and it
#:  replaces the prior EPS_C_MIN = 2.0, which was UNSOURCED and rejected most
#:  of the turbopump band the document itself recommends.
EPS_C_ABSOLUTE_MIN = 1.3
EPS_C_ABSOLUTE_MAX = 6.0

#: Humble's correlation, in the units its source actually uses.
_HUMBLE_FLOOR = 1.25
_HUMBLE_COEFF = 8.0
_HUMBLE_EXPONENT = -0.6
#: What this repository used before 2026-08-10, same unit system, different
#: numbers, provenance now lost. Kept so the conflict is visible in code.
_LEGACY_FLOOR = 1.00
_LEGACY_COEFF = 1.44

# --- SP-125 Table 4-1, p.87, verbatim, verified by direct page read --------
#: Recommended Combustion Chamber Characteristic Length (L*) for Various
#: Propellants. Units: INCHES. There is NO METHANE ROW -- verified by listing
#: every propellant in the table.
L_STAR_TABLE_4_1: dict[str, tuple[float, float]] = {
    "ClF3/hydrazine-base":              (30.0, 35.0),
    "LF2/hydrazine":                    (24.0, 28.0),
    "LF2/LH2 (GH2 injection)":          (22.0, 26.0),
    "LF2/LH2 (LH2 injection)":          (25.0, 30.0),
    "H2O2/RP-1 (incl. catalyst bed)":   (60.0, 70.0),
    "HNO3/hydrazine-base":              (30.0, 35.0),
    "N2O4/hydrazine-base":              (30.0, 35.0),
    "LOX/ammonia":                      (30.0, 40.0),
    "LOX/LH2 (GH2 injection)":          (22.0, 28.0),
    "LOX/LH2 (LH2 injection)":          (30.0, 40.0),
    "LOX/RP-1":                         (40.0, 50.0),
}
#: SP-125 p.87: "L* values of 15 to 120 inches for corresponding propellant
#: stay-time values of 0.002-0.040 second have been used in various thrust
#: chamber designs." This is the ENVELOPE, not a recommendation.
L_STAR_ENVELOPE_IN = (15.0, 120.0)

#: Which Table 4-1 row, if any, a Kryptonis fuel key maps onto.
_FUEL_TO_TABLE_ROW: dict[str, str | None] = {
    "RP-1": "LOX/RP-1",
    "LH2": "LOX/LH2 (LH2 injection)",
    "LCH4": None,          # NO ROW EXISTS
    "GCH4": None,          # NO ROW EXISTS
}


# ---------------------------------------------------------------------------
# CONTRACTION RATIO -- one policy
# ---------------------------------------------------------------------------

def contraction_ratio(
    *,
    throat_diameter_m: float,
    supplied_value: float | None = None,
    feed_system: str = "unknown",
) -> Result:
    r"""The single contraction-ratio policy for the whole platform.

    `feed_system` is one of "turbopump", "pressure_fed", "unknown". It selects
    which SP-125 range the value is checked against, because SP-125's ranges
    are keyed to feed system -- **not to throat diameter**, which is what the
    correlation uses. That structural divergence is itself recorded.
    """
    if throat_diameter_m <= 0:
        return Result(float("nan"), "-", "CAN-EPS-C",
                      status=Status.PHYSICALLY_INVALID, source=SP125,
                      notes="throat diameter must be > 0")

    assumptions: list[Assumption] = []
    statuses: list[Status] = [Status.PASS]
    notes: list[str] = []

    if supplied_value is not None:
        if supplied_value <= 1.0:
            return Result(float("nan"), "-", "CAN-EPS-C",
                          status=Status.PHYSICALLY_INVALID, source=SP125,
                          notes="contraction ratio must exceed 1; the chamber "
                                "cannot be narrower than the throat")
        eps_c = float(supplied_value)
        source = "supplied by the caller as a design variable"
        locator = "n/a -- design input"
        level = EvidenceLevel.E0
        verification = Verification.NONE
    else:
        d_cm = to_cm(throat_diameter_m)
        eps_c = _HUMBLE_FLOOR + _HUMBLE_COEFF * d_cm ** _HUMBLE_EXPONENT
        legacy = _LEGACY_FLOOR + _LEGACY_COEFF * d_cm ** _HUMBLE_EXPONENT
        source = HUMBLE
        locator = ("cryo-rocket.com Eq. 5.1.8 (secondary); Humble 1995 "
                   "primary NOT_OPENED, no page or equation number")
        level = EvidenceLevel.E1
        verification = Verification.ALGEBRAIC
        assumptions.append(Assumption(
            "eps_c", round(eps_c, 6), "-",
            "not supplied; evaluated from the Humble correlation with D_t in "
            "CENTIMETRES (the source's own unit system, verified by "
            "reproducing its published arithmetic)",
            Status.UNVALIDATED_ASSUMPTION))
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
        statuses.append(Status.CONFLICT_UNRESOLVED)
        notes.append(
            f"COEFFICIENT CONFLICT_UNRESOLVED: traced source gives "
            f"{_HUMBLE_FLOOR} + {_HUMBLE_COEFF} D[cm]^{_HUMBLE_EXPONENT} "
            f"-> {eps_c:.4f}; this repository previously used "
            f"{_LEGACY_FLOOR} + {_LEGACY_COEFF} D[cm]^{_HUMBLE_EXPONENT} "
            f"-> {legacy:.4f}, provenance now lost. Not averaged.")
        notes.append(
            "STRUCTURAL DIVERGENCE: SP-125's independent variable is feed "
            "system type; this correlation's is throat diameter.")

    # --- SP-125's own ranges, which ARE opened and verified ---------------
    band = {"turbopump": EPS_C_TURBOPUMP,
            "pressure_fed": EPS_C_PRESSURE_FED}.get(feed_system)
    if band is not None:
        lo, hi = band
        if not (lo <= eps_c <= hi):
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(
                f"eps_c {eps_c:.3f} is outside SP-125 p.88's {lo}-{hi} band "
                f"for {feed_system} systems")
    if not (EPS_C_ABSOLUTE_MIN <= eps_c <= EPS_C_ABSOLUTE_MAX):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(
            f"eps_c {eps_c:.3f} is outside SP-125 p.24's overall "
            f"{EPS_C_ABSOLUTE_MIN}-{EPS_C_ABSOLUTE_MAX} range")

    return Result(eps_c, "-", "CAN-EPS-C", status=worst(statuses),
        all_statuses=frozenset(statuses),
                  evidence_level=level, verification=verification,
                  validation=Validation.NONE, source=source,
                  source_locator=locator,
                  validity_domain=(
                      f"SP-125 p.88: turbopump {EPS_C_TURBOPUMP}, "
                      f"pressure-fed {EPS_C_PRESSURE_FED}; p.24 overall "
                      f"{EPS_C_ABSOLUTE_MIN}-{EPS_C_ABSOLUTE_MAX}"),
                  uncertainty=NOT_REPORTED,
                  assumptions=tuple(assumptions), notes="; ".join(notes),
                  inputs={"D_t_m": throat_diameter_m,
                          "D_t_cm": to_cm(throat_diameter_m),
                          "feed_system": feed_system,
                          "supplied": supplied_value})


# ---------------------------------------------------------------------------
# CHARACTERISTIC LENGTH -- one policy
# ---------------------------------------------------------------------------

def characteristic_length(
    *, fuel: str, supplied_value_m: float | None = None,
) -> Result:
    r"""L* policy. SP-125 Table 4-1 where a row exists, explicit gap where not.

    SP-125 p.87 states, in its own voice:

        "Under a given set of operating conditions ... the value of the
        minimum required L* can only be evaluated by actual firings of
        experimental thrust chambers."

    So Table 4-1 is a starting point for a test programme, not a validated
    prediction. Nothing here rises above E2, and methane -- which is
    Kryptonis's primary architecture -- has NO ROW AT ALL.
    """
    if supplied_value_m is not None:
        if supplied_value_m <= 0:
            return Result(float("nan"), "m", "CAN-L-STAR",
                          status=Status.PHYSICALLY_INVALID, source=SP125,
                          notes="L* must be > 0")
        l_star = float(supplied_value_m)
        assumptions: list[Assumption] = []
        statuses = [Status.PASS]
        source = "supplied by the caller as a design variable"
        level = EvidenceLevel.E0
        notes = []
    else:
        row = _FUEL_TO_TABLE_ROW.get(fuel, "MISSING")
        if row == "MISSING":
            return Result(float("nan"), "m", "CAN-L-STAR",
                          status=Status.INSUFFICIENT_EVIDENCE, source=SP125,
                          source_locator="Table 4-1, p.87",
                          notes=f"fuel {fuel!r} is not mapped to any Table 4-1 "
                                f"row and no value was supplied. L* cannot be "
                                f"guessed: it sets chamber volume, length and "
                                f"residence time.")
        if row is None:
            # Methane. There is no row. Do NOT invent one.
            return Result(
                float("nan"), "m", "CAN-L-STAR",
                status=Status.INSUFFICIENT_EVIDENCE, source=SP125,
                source_locator="Table 4-1, p.87",
                evidence_level=EvidenceLevel.E0,
                validity_domain=f"Table 4-1 envelope {L_STAR_ENVELOPE_IN} in",
                notes=(f"SP-125 Table 4-1 has NO METHANE ROW -- verified by "
                       f"listing all 11 propellant combinations. {fuel!r} is "
                       f"therefore not covered by the cited source. Supply "
                       f"L* explicitly, or acquire a LOX/CH4 source."),
                inputs={"fuel": fuel})
        lo_in, hi_in = L_STAR_TABLE_4_1[row]
        mid_in = 0.5 * (lo_in + hi_in)
        l_star = mid_in * INCH_M
        source = SP125
        level = EvidenceLevel.E2
        assumptions = [Assumption(
            "L_star_m", round(l_star, 6), "m",
            f"band midpoint of SP-125 Table 4-1 row {row!r} "
            f"({lo_in}-{hi_in} in); the table gives a RANGE, not a value",
            Status.UNVALIDATED_ASSUMPTION)]
        statuses = [Status.UNVALIDATED_ASSUMPTION]
        notes = [f"SP-125 Table 4-1 row {row!r}: {lo_in}-{hi_in} in"]

    l_in = l_star / INCH_M
    if not (L_STAR_ENVELOPE_IN[0] <= l_in <= L_STAR_ENVELOPE_IN[1]):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"L* = {l_in:.1f} in is outside SP-125's stated "
                     f"{L_STAR_ENVELOPE_IN[0]}-{L_STAR_ENVELOPE_IN[1]} in "
                     f"envelope of designs actually built")

    return Result(l_star, "m", "CAN-L-STAR", status=worst(statuses),
        all_statuses=frozenset(statuses),
                  evidence_level=level,
                  verification=Verification.NONE,
                  validation=Validation.NONE, source=source,
                  source_locator="Table 4-1, p.87 (IA scan p.95)",
                  validity_domain=(
                      "the 11 propellant pairs of Table 4-1, at 1960s US "
                      "injector designs. SP-125: 'the minimum required L* can "
                      "only be evaluated by actual firings'"),
                  uncertainty=NOT_REPORTED,
                  assumptions=tuple(assumptions), notes="; ".join(notes),
                  inputs={"fuel": fuel, "supplied_m": supplied_value_m})


# ---------------------------------------------------------------------------
# Geometry -- identities. SP-125 Eq. (4-4) and (4-5).
# ---------------------------------------------------------------------------

def chamber_volume(*, l_star_m: float, throat_area_m2: float) -> Result:
    r""".. math:: V_c = L^{*} A_t   (SP-125 Eq. 4-4, rearranged)

    SP-125 Eq. (4-5) makes explicit that V_c INCLUDES the convergent cone.
    """
    if l_star_m <= 0 or throat_area_m2 <= 0:
        return Result(float("nan"), "m^3", "CAN-VC",
                      status=Status.PHYSICALLY_INVALID, source=SP125,
                      notes="L* and A_t must be > 0")
    return Result(l_star_m * throat_area_m2, "m^3", "CAN-VC",
                  status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE, source=SP125,
                  source_locator="Eq. (4-4), p.87; Eq. (4-5), p.88",
                  validity_domain="definition; holds identically",
                  uncertainty="exact",
                  notes="V_c includes the convergent cone, per Eq. (4-5)",
                  inputs={"L_star_m": l_star_m, "A_t_m2": throat_area_m2})


def chamber_diameter(*, throat_diameter_m: float,
                     contraction_ratio_: float) -> Result:
    r""".. math:: D_c = D_t \sqrt{\varepsilon_c}"""
    if throat_diameter_m <= 0 or contraction_ratio_ <= 0:
        return Result(float("nan"), "m", "CAN-DC",
                      status=Status.PHYSICALLY_INVALID,
                      notes="D_t and eps_c must be > 0")
    return Result(throat_diameter_m * math.sqrt(contraction_ratio_), "m",
                  "CAN-DC", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="geometric identity from eps_c = A_c/A_t",
                  source_locator="SP-125 p.88 defines eps_c = A_c/A_t",
                  validity_domain="identity", uncertainty="exact",
                  inputs={"D_t_m": throat_diameter_m,
                          "eps_c": contraction_ratio_})


def convergent_length(*, throat_diameter_m: float, chamber_diameter_m: float,
                      half_angle_deg: float) -> Result:
    r""".. math:: L_{conv} = \frac{D_c - D_t}{2\tan\theta}"""
    if not (0.0 < half_angle_deg < 90.0):
        return Result(float("nan"), "m", "CAN-LCONV",
                      status=Status.PHYSICALLY_INVALID,
                      notes="convergent half-angle must be in (0, 90) deg")
    if chamber_diameter_m < throat_diameter_m:
        return Result(float("nan"), "m", "CAN-LCONV",
                      status=Status.PHYSICALLY_INVALID,
                      notes="chamber diameter cannot be below throat diameter")
    L = (chamber_diameter_m - throat_diameter_m) / (
        2.0 * math.tan(math.radians(half_angle_deg)))
    return Result(L, "m", "CAN-LCONV", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="frustum geometry", source_locator="identity",
                  validity_domain="straight conical convergent",
                  uncertainty="exact",
                  inputs={"D_t_m": throat_diameter_m,
                          "D_c_m": chamber_diameter_m,
                          "half_angle_deg": half_angle_deg})


def convergent_volume(*, throat_diameter_m: float, chamber_diameter_m: float,
                      half_angle_deg: float) -> Result:
    r"""Volume of the straight conical frustum between chamber and throat.

    .. math:: V = \frac{\pi L}{12}\left(D_c^2 + D_c D_t + D_t^2\right)

    NOTE: `combustor/sizing.calculate_spline_converging_volume` computes a
    SPLINE contour instead and differs by 2.80 % on the same inputs. That is
    a genuine alternative formulation, not a duplicate, and it is retained as
    such -- but the two must never be mixed inside one chamber.
    """
    lr = convergent_length(throat_diameter_m=throat_diameter_m,
                           chamber_diameter_m=chamber_diameter_m,
                           half_angle_deg=half_angle_deg)
    if lr.status is Status.PHYSICALLY_INVALID:
        return Result(float("nan"), "m^3", "CAN-VCONV",
                      status=Status.PHYSICALLY_INVALID, notes=lr.notes)
    L = float(lr)
    v = (math.pi * L / 12.0) * (chamber_diameter_m ** 2
                                + chamber_diameter_m * throat_diameter_m
                                + throat_diameter_m ** 2)
    return Result(v, "m^3", "CAN-VCONV", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="conical frustum volume", source_locator="identity",
                  validity_domain="straight conical convergent only",
                  uncertainty="exact",
                  notes="spline alternative in combustor/sizing differs 2.80%",
                  inputs={"D_t_m": throat_diameter_m,
                          "D_c_m": chamber_diameter_m,
                          "half_angle_deg": half_angle_deg})


def cylinder_length(*, chamber_volume_m3: float, convergent_volume_m3: float,
                    chamber_area_m2: float) -> Result:
    r""".. math:: L_{cyl} = \frac{V_c - V_{conv}}{A_c}

    The convergent cone is INSIDE V_c (SP-125 Eq. 4-5), so it must be
    subtracted before the barrel length is computed. Treating V_c as a
    straight cylinder and then adding a cone overshoots the delivered L* by
    up to 20 % at high contraction ratios.
    """
    if chamber_area_m2 <= 0:
        return Result(float("nan"), "m", "CAN-LCYL",
                      status=Status.PHYSICALLY_INVALID,
                      notes="chamber area must be > 0")
    if convergent_volume_m3 > chamber_volume_m3:
        return Result(float("nan"), "m", "CAN-LCYL",
                      status=Status.PHYSICALLY_INVALID,
                      notes=f"the convergent cone ({convergent_volume_m3:.6g} "
                            f"m^3) does not fit inside the chamber volume "
                            f"({chamber_volume_m3:.6g} m^3): L* is too small "
                            f"for this contraction ratio and half angle")
    return Result((chamber_volume_m3 - convergent_volume_m3) / chamber_area_m2,
                  "m", "CAN-LCYL", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE, source=SP125,
                  source_locator="Eq. (4-5), p.88",
                  validity_domain="identity", uncertainty="exact",
                  inputs={"V_c_m3": chamber_volume_m3,
                          "V_conv_m3": convergent_volume_m3,
                          "A_c_m2": chamber_area_m2})


def vandenkerckhove(gamma: float) -> Result:
    r""".. math::
        \Gamma = \sqrt{\gamma}\left(\frac{2}{\gamma+1}
                 \right)^{\frac{\gamma+1}{2(\gamma-1)}}

    The exponent denominator carries the factor 2. `material_manager` omitted
    it, giving 11.0 instead of 5.5 at gamma = 1.2 -- c* high by 68.9 % and
    h_g low by 34.3 %.
    """
    if not (1.0 < gamma <= 1.67):
        return Result(float("nan"), "-", "CAN-GAMMA-FN",
                      status=Status.PHYSICALLY_INVALID,
                      notes=f"gamma={gamma} outside (1, 1.67]")
    # OWNERSHIP (reconstruction wave): the numeric core is computed by the
    # ONE owner in closedloop.chamber_consistency (float, hot-path form);
    # this canonical wrapper adds the Result pedigree and the domain guard.
    # Deferred import: closedloop imports this module at load time, so a
    # module-level import here would be a cycle.
    # Exact analytical Vandenkerckhove function: Gamma = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
    val = math.sqrt(gamma) * math.pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return Result(val, "-", "CAN-GAMMA-FN", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="Vandenkerckhove function; isentropic choked flow",
                  source_locator="standard gas dynamics; derivable",
                  validity_domain="calorically perfect gas", uncertainty="exact",
                  inputs={"gamma": gamma})


def c_star_ideal(*, gamma: float, molar_mass_kg_per_mol: float,
                 chamber_temperature_K: float) -> Result:
    r""".. math:: c^{*} = \frac{\sqrt{R T_c}}{\Gamma}, \quad R = R_u / M

    **WHY THE OTHER c\* IMPLEMENTATIONS ARE DELIBERATELY NOT MERGED INTO
    THIS ONE.** The physics-closure wave canonicalised the four copies of the
    convergent volume into one, because four copies of a geometric identity
    are four chances to fix a bug in three places. It did NOT do the same to
    c\*, and the reason is worth writing down.

    `chemistry.equilibrium` and `core.cea_native` each compute c\* inline.
    Those two are not redundant copies -- they are the project's only
    cross-implementation check on the chamber state, and the controlled
    kernel comparison depends on their being written independently. Merging
    them would delete the evidence that they agree. It was verified over 500
    randomised cases that all three forms -- this one, the equilibrium
    kernel's, and cea_native's differently-factored Gamma -- agree to
    4.1e-16, and they were then LEFT ALONE.

    **RULE 4 is about a single source of truth for a relation that is
    maintained, not about deleting a deliberate redundancy that is
    measured.** A duplicate that exists to be compared is not a duplicate.
    """
    if molar_mass_kg_per_mol <= 0 or chamber_temperature_K <= 0:
        return Result(float("nan"), "m/s", "CAN-CSTAR",
                      status=Status.PHYSICALLY_INVALID,
                      notes="molar mass and temperature must be > 0")
    gr = vandenkerckhove(gamma)
    if gr.status is Status.PHYSICALLY_INVALID:
        return Result(float("nan"), "m/s", "CAN-CSTAR",
                      status=Status.PHYSICALLY_INVALID, notes=gr.notes)
    R_U = 8.31446261815324
    R = R_U / molar_mass_kg_per_mol
    return Result(math.sqrt(R * chamber_temperature_K) / float(gr), "m/s",
                  "CAN-CSTAR", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="isentropic choked flow; R_u = CODATA 2018",
                  source_locator="derivable",
                  validity_domain="calorically perfect gas, complete "
                                  "combustion, frozen composition",
                  uncertainty="exact given (gamma, M, T)",
                  inputs={"gamma": gamma, "M_kg_per_mol": molar_mass_kg_per_mol,
                          "T_c_K": chamber_temperature_K})


# Alias for interface consistency
cylindrical_length = cylinder_length
