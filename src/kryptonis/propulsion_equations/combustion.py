"""
Combustion completion. The honest answer is UNKNOWN, and this module says so.
=============================================================================

WHAT THIS MODULE IS FOR
-----------------------
`solve_chamber_consistency` has been returning `eta_c* = 0.98` from a hardcoded
branch that is ALWAYS taken (every modern chamber is supercritical on LOX, so
the d^2-law fallback beneath it is dead code). The number is not wrong so much
as unaccountable: nothing downstream can tell that it was assumed rather than
computed.

This module does not replace it with a better number. **It replaces it with the
same number plus the truth about where it came from.**

WHY NOT eta_c* = f(L*)
----------------------
It is tempting, because SP-125 Figure 4-7 plots c* against L*. Four things
rule it out, and they were checked rather than assumed:

  1. **Figure 4-7 is ONE chamber.** Caption, verbatim: *"Effect of L* on c*
     value of experimental thrust chamber."* Singular. SP-125 names neither its
     propellant, nor its chamber pressure, nor its injector, nor its thrust --
     every one of those was checked on pp.87-88 and is NOT_STATED.

  2. **SP-125 says in its own voice that L* is not predictive.** p.87:
     *"Under a given set of operating conditions, such as type of propellants,
     mixture ratio, chamber pressure, injector design, and chamber geometry,
     the value of the minimum required L* can only be evaluated by actual
     firings of experimental thrust chambers."*

  3. **A later NASA source refutes the mapping directly.** Hoggatt & Leahy,
     NASA MSFC, NTRS 20090034473, from fired subscale test data:
     *"It is impossible to precisely correlate C*, C* eta, or Isp to L' or L*
     alone, since test data clearly shows that Injector Type, Injector Density,
     Momentum Ratio, Fuel Injection Temperature, Chamber Pressure, and Mixture
     Ratio also affect these performance values."*

  4. **The measured LOX/CH4 data that exist contradict any single curve.** The
     one fired LOX/GCH4 article with a stated L* (Kim et al., EUCASS 2019-574,
     L* = 1.35 m) reports eta_c* of **72.2-82.7 %**. Kryptonis's 30 kN case
     runs L* = 1.0 m and assumes 0.98. Those two facts cannot both sit on one
     eta(L*) curve.

So: **ETA_LSTAR_MAPPING_NOT_JUSTIFIED**, on primary-source evidence, and the
correct engineering answer is that eta_c* is UNKNOWN for this platform.

WHAT CHANGES IN BEHAVIOUR
-------------------------
Numerically, nothing. `assumed_c_star_efficiency()` returns the same 0.98 the
closed loop has always used, so no design moves and no test result changes.
What changes is that it now returns a `Result` carrying
`INSUFFICIENT_EVIDENCE`, the four reasons above, and the exact experiment that
would close the gap. A caller that wants a number still gets one; a caller that
wants to know whether to trust it can now find out.
"""

from __future__ import annotations

import math

from kryptonis.propulsion_equations.units import (
    Assumption, EvidenceLevel, NOT_REPORTED, Result, Status, Validation,
    Verification, worst,
)

__all__ = [
    "LEGACY_ASSUMED_ETA", "SP125_LSTAR_ENVELOPE_M", "LOX_CH4_LSTAR_EVIDENCE",
    "assumed_c_star_efficiency", "l_star_requirement",
    "chamber_bulk_residence_time", "sp125_stay_time",
    "evaluate_combustion_completion",
]

#: The value the closed loop has always used. Preserved EXACTLY so that
#: recording the truth about it does not silently move any design.
LEGACY_ASSUMED_ETA = 0.98

#: SP-125 p.87, VERBATIM: "L* values of 15 to 120 inches ... have been used in
#: various thrust chamber designs." An envelope of practice, not a range of
#: validity, and not a recommendation.
SP125_LSTAR_ENVELOPE_M = (15.0 * 0.0254, 120.0 * 0.0254)

#: Every LOX/CH4 L* this project has actually traced to a primary source that
#: STATES it or from which it is directly computable. Recorded as data because
#: the SPREAD is the finding: 0.72-3.55 m across sub-kN hardware.
LOX_CH4_LSTAR_EVIDENCE: tuple[dict, ...] = (
    dict(key="KIM-2019", l_star_m=1.35, stated=True,
         oxidiser="liquid", fuel="gaseous", thrust_N=(88.0, 121.0),
         chamber_pressure_MPa=(1.04, 1.37), mixture_ratio=(3.68, 5.38),
         c_star_efficiency_pct=(72.2, 82.7),
         citation="Kim J.S. et al., EUCASS 2019, DOI 10.13009/EUCASS2019-574",
         locator="Table 1 (L*), Table 2 and Figure 5 (c*)",
         note="the only fired LOX/CH4 article found that STATES an L*"),
    dict(key="HULKA-2010-PSU", l_star_m=1.98, stated=False,
         oxidiser="liquid", fuel="gaseous", thrust_N=None,
         chamber_pressure_MPa=None, mixture_ratio=(2.6, 3.0),
         c_star_efficiency_pct=None,
         citation="Hulka & Jones, AIAA 2010-6801 (NTRS 20100034924)",
         locator="Table I, p.3: Dc = 1.0 in, Lc = 14.49 in, CR = 5.38",
         note="COMPUTED from stated geometry, cylindrical approximation, so "
              "an UPPER BOUND -- the convergent taper is not subtracted. The "
              "authors do not perform this calculation."),
    dict(key="KANG-2022", l_star_m=3.55, stated=False,
         oxidiser="liquid", fuel="gaseous", thrust_N=(545.0, 545.0),
         chamber_pressure_MPa=(2.5, 3.5), mixture_ratio=(3.0, 3.0),
         c_star_efficiency_pct=None,
         citation="Kang C. et al., EUCASS 2022, DOI 10.13009/EUCASS2022-7247",
         locator="p.4: Dc = 45.75 mm, faceplate-to-throat 177.38 mm, "
                 "Dt = 10.23 mm",
         note="COMPUTED, same cylindrical-approximation upper-bound caveat"),
)


# ---------------------------------------------------------------------------
# eta_c*
# ---------------------------------------------------------------------------

def assumed_c_star_efficiency(
    *, supplied_value: float | None = None,
    l_star_m: float | None = None,
    propellant: str = "unknown",
) -> Result:
    r"""c\* efficiency. **UNKNOWN unless the caller supplies it.**

    Returns `LEGACY_ASSUMED_ETA` with `INSUFFICIENT_EVIDENCE` when nothing is
    supplied -- the same number the loop has always used, now carrying its own
    provenance instead of hiding inside an `if` branch.

    `l_star_m` is accepted ONLY so the returned notes can state how far the
    request sits from the one measured LOX/CH4 point. **It does not influence
    the value.** Making it do so would be exactly the unjustified eta(L*)
    mapping this module exists to refuse.
    """
    if supplied_value is not None:
        if not (0.0 < supplied_value <= 1.0):
            return Result(
                float("nan"), "-", "CAN-ETA-CSTAR",
                status=Status.PHYSICALLY_INVALID,
                notes=f"c* efficiency {supplied_value} is outside (0, 1]. "
                      f"Above 1 the chamber would out-perform equilibrium; at "
                      f"or below 0 no combustion occurred.",
                inputs={"supplied": supplied_value})
        return Result(
            float(supplied_value), "-", "CAN-ETA-CSTAR", status=Status.PASS,
            evidence_level=EvidenceLevel.E0,
            verification=Verification.NONE, validation=Validation.NONE,
            source="supplied by the caller",
            source_locator="n/a -- caller-supplied",
            validity_domain="whatever the caller's own evidence supports",
            uncertainty=NOT_REPORTED,
            notes="Kryptonis does not check a supplied efficiency against any "
                  "measurement, because it holds none for this propellant at "
                  "this scale.",
            inputs={"supplied": supplied_value, "propellant": propellant})

    statuses = [Status.INSUFFICIENT_EVIDENCE, Status.UNVALIDATED_ASSUMPTION]
    notes = [
        "eta_c* is UNKNOWN for this platform. The value returned is the "
        f"pre-existing hardcoded {LEGACY_ASSUMED_ETA}, preserved so that "
        "recording the truth about it does not silently move any design.",
        "ETA_LSTAR_MAPPING_NOT_JUSTIFIED: SP-125 Figure 4-7 is ONE unnamed "
        "chamber (propellant, Pc, injector and thrust all NOT_STATED on "
        "pp.87-88), and Hoggatt & Leahy (NTRS 20090034473) state from fired "
        "test data that c* eta cannot be correlated to L* alone.",
    ]
    if l_star_m is not None:
        kim = LOX_CH4_LSTAR_EVIDENCE[0]
        notes.append(
            f"the only fired LOX/CH4 article with a stated L* "
            f"({kim['citation']}) ran L* = {kim['l_star_m']} m at "
            f"{kim['thrust_N'][0]:.0f}-{kim['thrust_N'][1]:.0f} N and measured "
            f"eta_c* = {kim['c_star_efficiency_pct'][0]}-"
            f"{kim['c_star_efficiency_pct'][1]} %, which is far below the "
            f"{LEGACY_ASSUMED_ETA} assumed here. The requested L* is "
            f"{l_star_m:.3f} m.")

    return Result(
        LEGACY_ASSUMED_ETA, "-", "CAN-ETA-CSTAR",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE, validation=Validation.NONE,
        source="NO SOURCE. This is a carried-forward assumption.",
        source_locator="none",
        validity_domain="VALIDITY_DOMAIN_UNKNOWN",
        uncertainty=NOT_REPORTED,
        assumptions=(Assumption(
            "c_star_efficiency", LEGACY_ASSUMED_ETA, "-",
            "hardcoded in solve_chamber_consistency since before this audit; "
            "no source, no measurement, no derivation",
            Status.INSUFFICIENT_EVIDENCE),),
        notes=" | ".join(notes),
        inputs={"l_star_m": l_star_m, "propellant": propellant})


# ---------------------------------------------------------------------------
# L*
# ---------------------------------------------------------------------------

def l_star_requirement(*, propellant: str,
                       supplied_value_m: float | None = None) -> Result:
    r"""Required L\*. **Not computable.** SP-125 says so in its own voice.

    L\* is a DESIGN VARIABLE selected from practice, not a quantity derived
    from a combustion model. SP-125 p.87: *"the value of the minimum required
    L\* can only be evaluated by actual firings of experimental thrust
    chambers."*

    For LOX/CH4 there is additionally no Table 4-1 row, and the traced fired
    evidence spans **0.72-3.55 m across sub-kN hardware** -- a 4.9x spread,
    none of it at engine scale.
    """
    # Normalise every separator a caller might plausibly use. Missing one is
    # how "LOX/CH4" silently fell through to the non-methane branch and got
    # told to consult a Table 4-1 row that does not exist.
    _norm = propellant.upper()
    for _sep in ("-", "_", "/", " ", ".", "+"):
        _norm = _norm.replace(_sep, "")
    is_methane = ("CH4" in _norm or "METHANE" in _norm or "LNG" in _norm
                  or "NATURALGAS" in _norm)
    statuses = [Status.INSUFFICIENT_EVIDENCE]
    notes = []
    if is_methane:
        notes.append(
            "NO_DIRECT_SOURCE_FOR_LOX_CH4_LSTAR: SP-125 Table 4-1 has no "
            "methane row (all 11 rows enumerated and checked).")
        spread = [e["l_star_m"] for e in LOX_CH4_LSTAR_EVIDENCE]
        notes.append(
            f"traced fired LOX/CH4 evidence spans {min(spread):.2f}-"
            f"{max(spread):.2f} m ({max(spread)/min(spread):.1f}x), all at "
            f"sub-kN thrust; 2 of the 3 are COMPUTED upper bounds from stated "
            f"geometry, not values the authors report.")
    else:
        notes.append("see SP-125 Table 4-1 via "
                     "kryptonis.canonical.chamber.characteristic_length")

    value = float(supplied_value_m) if supplied_value_m is not None else float("nan")
    if supplied_value_m is not None:
        if supplied_value_m <= 0:
            return Result(float("nan"), "m", "CAN-LSTAR-REQ",
                          status=Status.PHYSICALLY_INVALID,
                          notes="L* must be > 0")
        lo, hi = SP125_LSTAR_ENVELOPE_M
        if not (lo <= supplied_value_m <= hi):
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(
                f"L* = {supplied_value_m:.3f} m is outside SP-125's stated "
                f"envelope of built designs, {lo:.3f}-{hi:.3f} m")

    return Result(
        value, "m", "CAN-LSTAR-REQ",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE, validation=Validation.NONE,
        source="NASA SP-125 p.87 states L* is not predictable; "
               "LOX/CH4 evidence traced separately",
        source_locator="SP-125 p.87 (IA scan p.95)",
        validity_domain=f"SP-125 envelope of practice {SP125_LSTAR_ENVELOPE_M}",
        uncertainty=NOT_REPORTED, notes=" | ".join(notes),
        inputs={"propellant": propellant, "supplied_m": supplied_value_m})


# ---------------------------------------------------------------------------
# Residence time -- TWO quantities, not one
# ---------------------------------------------------------------------------

def chamber_bulk_residence_time(*, chamber_volume_m3: float,
                                density_kg_m3: float,
                                mdot_kg_s: float) -> Result:
    r"""Bulk residence time at a SINGLE STATED STATION.

    .. math:: \tau = \frac{V_c\,\rho}{\dot m}

    **CORRECTED BY THE PHYSICS-CLOSURE WAVE.** This docstring previously said
    that `solve_chamber_consistency` computes this using the *station-2 static
    density -- the end of the barrel*. **It does not, and it never did.** The
    live solver uses

        rho = P_inj / (R * T_chamber_actual)

    which is the **INJECTOR-FACE STAGNATION** density: no Rayleigh
    correction, no Mach correction. At the 30 kN reference point the two
    stations differ by **1.91 %** in density and therefore in residence time
    -- of which the Rayleigh pressure drop is -2.08 % and the
    stagnation-temperature factor +0.17 %.

    The drift was invisible to the residual suite BY CONSTRUCTION rather than
    by accident: residual R7 compares `V rho_static/mdot` against `V/(v A_c)`
    and **both sides use the static station**, so it closes to 1e-9 as an
    identity between two static-station forms and never touches the
    stagnation-station number the report prints.

    Nothing about the live computation is changed here. What is corrected is
    the STATEMENT, because a registry entry that describes code which does not
    exist is worse than no entry at all.

    The station is therefore whatever the CALLER's `density_kg_m3` is, and the
    caller must know which one it passed. This is a well-defined quantity and
    it is NOT the same as SP-125 Eq. (4-3); see
    `docs/research/RESIDENCE_TIME_DEFINITION_AUDIT.md` and
    `docs/audit/CHAMBER_DENSITY_RESIDENCE_CLOSURE.md`.
    """
    if chamber_volume_m3 <= 0:
        return Result(float("nan"), "s", "CAN-TAU-BULK",
                      status=Status.PHYSICALLY_INVALID,
                      notes="chamber volume must be > 0")
    if density_kg_m3 <= 0 or mdot_kg_s <= 0:
        return Result(float("nan"), "s", "CAN-TAU-BULK",
                      status=Status.PHYSICALLY_INVALID,
                      notes="density and mass flow must be > 0")
    return Result(chamber_volume_m3 * density_kg_m3 / mdot_kg_s, "s",
                  "CAN-TAU-BULK", status=Status.PASS,
                  evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC,
                  validation=Validation.NONE,
                  source="mass-residence identity at one stated station",
                  source_locator="identity",
                  validity_domain="the station whose density was supplied",
                  uncertainty="exact given the inputs",
                  notes="the STATION MUST BE STATED. This is not SP-125 "
                        "Eq. (4-3), whose 'average specific volume' is a "
                        "different and undefined average.",
                  inputs={"V_c_m3": chamber_volume_m3,
                          "rho_kg_m3": density_kg_m3,
                          "mdot_kg_s": mdot_kg_s})


def sp125_stay_time(*, chamber_volume_m3: float,
                    average_specific_volume_m3_kg: float | None,
                    mdot_kg_s: float) -> Result:
    r"""SP-125 Eq. (4-3) stay time. **The source does not define its average.**

    .. math:: V_c = \dot W_{tc}\,\bar{V}\,t_s

    SP-125 p.86 defines the symbols and no more:
        Vc  = chamber volume, ft^3
        Wtc = propellant mass flow rate, lb/sec
        V   = average specific volume, ft^3/lb
        ts  = propellant stay time, sec

    **Nothing on that page says what the average is taken OVER** -- not over
    which stations, not whether it is mass- or volume-weighted, not whether the
    propellant is counted as liquid at the injector face or as product gas.
    That was checked directly and the answer is NOT_STATED.

    So this function REFUSES to guess. Supply the average yourself, or get
    `INSUFFICIENT_EVIDENCE`.
    """
    if chamber_volume_m3 <= 0 or mdot_kg_s <= 0:
        return Result(float("nan"), "s", "CAN-TAU-SP125",
                      status=Status.PHYSICALLY_INVALID,
                      notes="chamber volume and mass flow must be > 0")
    if average_specific_volume_m3_kg is None:
        return Result(
            float("nan"), "s", "CAN-TAU-SP125",
            status=Status.INSUFFICIENT_EVIDENCE,
            evidence_level=EvidenceLevel.E0,
            source="NASA SP-125 Eq. (4-3), p.86",
            source_locator="p.86 (IA scan p.94)",
            validity_domain="VALIDITY_DOMAIN_UNKNOWN",
            notes="SP-125 defines V-bar only as 'average specific volume, "
                  "ft^3/lb' and NEVER states what it averages over. Choosing "
                  "an averaging convention here would be inventing the "
                  "definition, which is the one thing this audit forbids. "
                  "RESIDENCE_TIME_DEFINITION_UNRESOLVED.",
            inputs={"V_c_m3": chamber_volume_m3, "mdot_kg_s": mdot_kg_s})
    if average_specific_volume_m3_kg <= 0:
        return Result(float("nan"), "s", "CAN-TAU-SP125",
                      status=Status.PHYSICALLY_INVALID,
                      notes="specific volume must be > 0")
    ts = chamber_volume_m3 / (mdot_kg_s * average_specific_volume_m3_kg)
    return Result(
        ts, "s", "CAN-TAU-SP125",
        status=Status.UNVALIDATED_ASSUMPTION,
        all_statuses=frozenset({Status.UNVALIDATED_ASSUMPTION}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="NASA SP-125 Eq. (4-3), p.86",
        source_locator="p.86 (IA scan p.94)",
        validity_domain="VALIDITY_DOMAIN_UNKNOWN -- the averaging convention "
                        "is the caller's, not the source's",
        uncertainty=NOT_REPORTED,
        assumptions=(Assumption(
            "average_specific_volume_m3_kg", average_specific_volume_m3_kg,
            "m^3/kg",
            "supplied by the caller; SP-125 does not define what this is an "
            "average of", Status.UNVALIDATED_ASSUMPTION),),
        notes="SP-125 also states 'the stay time ts is independent of the "
              "combustion chamber geometry', which is a property no bulk "
              "V*rho/mdot form has.",
        inputs={"V_c_m3": chamber_volume_m3, "mdot_kg_s": mdot_kg_s,
                "V_bar_m3_kg": average_specific_volume_m3_kg})


# ---------------------------------------------------------------------------
# The contract
# ---------------------------------------------------------------------------

def evaluate_combustion_completion(
    *,
    propellant: str = "unknown",
    l_star_m: float | None = None,
    chamber_volume_m3: float | None = None,
    density_kg_m3: float | None = None,
    mdot_kg_s: float | None = None,
    supplied_c_star_efficiency: float | None = None,
) -> dict:
    r"""The combustion-completion contract.

    **NOT CONNECTED TO PRODUCTION.** Phase 13 of the brief: connect only if the
    evidence supports a model. It does not, so this returns
    `INSUFFICIENT_EVIDENCE` and the closed loop keeps its existing behaviour.

    Returns a dict of `Result` objects so that each quantity carries its own
    status rather than the caller receiving one blended verdict.
    """
    eta = assumed_c_star_efficiency(
        supplied_value=supplied_c_star_efficiency, l_star_m=l_star_m,
        propellant=propellant)
    lstar = l_star_requirement(propellant=propellant,
                               supplied_value_m=l_star_m)
    if (chamber_volume_m3 is not None and density_kg_m3 is not None
            and mdot_kg_s is not None):
        tau = chamber_bulk_residence_time(
            chamber_volume_m3=chamber_volume_m3,
            density_kg_m3=density_kg_m3, mdot_kg_s=mdot_kg_s)
    else:
        tau = Result(float("nan"), "s", "CAN-TAU-BULK",
                     status=Status.INSUFFICIENT_EVIDENCE,
                     notes="chamber volume, density and mass flow are all "
                           "required, and the STATION of the density must be "
                           "stated by the caller")
    overall = worst([eta.status, lstar.status, tau.status])
    return {
        "c_star_efficiency": eta,
        "l_star_requirement": lstar,
        "residence_time": tau,
        "status": overall,
        "model_available": False,
        "reason": (
            "ETA_LSTAR_MAPPING_NOT_JUSTIFIED and "
            "RESIDENCE_TIME_DEFINITION_UNRESOLVED. No reduced-order "
            "combustion-completion closure is supported by the evidence this "
            "project holds. See docs/analytical/"
            "COMBUSTION_COMPLETION_MODEL_SPEC.md."),
        "missing_experiment": (
            "a fired LOX/CH4 chamber at 10-100 kN reporting: absolute measured "
            "c* AND the exact theoretical baseline, L*, injector type, Pc, MR, "
            "throat area, instrument uncertainty, and heat-loss correction "
            "status. No source examined reports all of these together."),
    }
