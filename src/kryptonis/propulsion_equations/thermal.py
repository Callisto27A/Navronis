r"""
THE canonical combustor thermal chain.
======================================

One chain, one implementation, one place where every constant is either traced
to a primary source or labelled as having none:

    T_0, M -> h_g (Bartz)  ->  q = h_g (T_aw - T_wg)
                                     |
                             wall conduction (CYLINDRICAL, or planar with a
                             stated and BOUNDED approximation error)
                                     |
                             q = h_c (T_wc - T_bulk)   <- coolant-side closure
                                     |
                             energy balance on the coolant (ENTHALPY, not cp dT)

Every previous version of this chain in this repository was broken in the same
way: the wall temperature was never solved, it was *asserted*. Three separate
production sites asserted three different values --

    orchestrator          T_wg = min(900, 0.22 T_c + 100)
    material_manager      T_wg = 0.40 T_c            (and sigma frozen there)
    screens               Q_total = q_throat * A * 0.35

-- none of which has a source, and the three disagree with each other by more
than the effect any correlation choice in this module has. This module exists
so that a wall temperature is either SOLVED or declared UNKNOWN, and never
guessed silently.


WHAT THE RESEARCH ACTUALLY ESTABLISHED (2026-08-20)
---------------------------------------------------

1. ``Nu = 0.023 Re^0.8 Pr^0.4`` IS NOT THE DITTUS-BOELTER EQUATION.
   Winterton (1998) establishes that Dittus & Boelter (1930) published
   ``0.0243 Re^0.8 Pr^0.4`` for heating and ``0.0265 Re^0.8 Pr^0.3`` for
   cooling -- two different coefficients. The universally cited 0.023 form is
   McAdams, *Heat Transmission*, 2nd ed. (1942), for which McAdams "did not
   give a clear reference". This repository cited "Dittus & Boelter, 1930"
   next to the 0.023 form. That is a documented misattribution and is
   corrected here.

2. THE ONE PRIMARY SOURCE THAT WAS ACTUALLY OPENED is Taylor, NASA TN D-4332
   (1968), and it is the strongest coolant-side evidence this project has: a
   fitted correlation with STATED ranges, a STATED accuracy, and -- unusually --
   its author's OWN explicit disclaimer of the near-critical region. It is
   fitted on HYDROGEN. Applying it to methane is out of its fluid domain and
   says so.

3. TAYLOR, NASA TM X-52437 (1968) gives the curvature correction actually used
   in rocket cooling literature: ``[Re_b (r/R)^2]^{+/-0.05}``, positive on the
   concave wall, negative on the convex, attributed to Ito (1959). JAXA's own
   regenerative-cooling model (Negishi et al.) uses this exact factor. So the
   curvature term is NOT invented -- but its validity range is NOT_REPORTED.

4. THERE IS NO CRITICAL HEAT FLUX ABOVE THE CRITICAL PRESSURE. Pioro & Mokry
   (2011), verbatim: *"At supercritical pressures there is no liquid-vapour
   phase transition; therefore, there is no such phenomenon as Critical Heat
   Flux (CHF) or dryout."* Methane at any realistic regenerative channel
   pressure is above its 4.599 MPa critical pressure. A Zuber pool-boiling CHF
   margin computed there is not conservative, not optimistic, and not wrong by
   some percentage -- it is a category error. The analogous failure mode is
   named HEAT TRANSFER DETERIORATION and it is a different phenomenon.

5. NASA SP-8087, the design-criteria monograph for fluid-cooled combustion
   chambers, CONTAINS NO HEAT-TRANSFER EQUATIONS. It was opened and checked.
   It must not be cited as the source of a correlation.

6. THE BEST AVAILABLE UNCERTAINTY for supercritical regenerative cooling is
   NOT +/-10%. Locke & Landrum (2008) evaluated seven correlations against
   2992 measured points; the best performer needed a +/-56% band to cover 95%
   of the data. Anything tighter is not supported by the open literature.


WHAT COULD NOT BE OPENED, AND THEREFORE IS NOT ASSERTED
-------------------------------------------------------
Gnielinski (1976), Sieder & Tate (1936), Colebrook (1938-39), Haaland (1983),
Petukhov (1970) and every dedicated Pizzarelli supercritical-methane paper are
all behind 403s or robots.txt. Their validity ranges here are recorded as
SECONDHAND with the disagreement between secondary sources shown, not averaged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Sequence

from kryptonis.propulsion_equations.units import (
    Assumption,
    EvidenceLevel,
    NOT_REPORTED,
    Result,
    Status,
    UNKNOWN_DOMAIN,
    Validation,
    Verification,
    worst,
)

__all__ = [
    # sources
    "TAYLOR_TN_D_4332", "TAYLOR_TM_X_52437", "MCADAMS_1942",
    "DITTUS_BOELTER_1930", "WINTERTON_1998", "GNIELINSKI_1976",
    "SIEDER_TATE_1936", "HAALAND_1983", "COLEBROOK_1939",
    "LOCKE_LANDRUM_2008", "PIORO_MOKRY_2011", "LAMANNA_2025",
    # gas side
    "gas_side_heat_flux", "assumed_hot_wall_temperature",
    "LEGACY_HOT_WALL_FRACTION", "LEGACY_HOT_WALL_OFFSET_K",
    "LEGACY_HOT_WALL_CAP_K",
    # wall
    "WallGeometry", "planar_wall_resistance", "cylindrical_wall_resistance",
    "wall_resistance", "mean_wall_conductivity",
    # friction
    "friction_factor_haaland", "friction_factor_colebrook",
    # coolant side
    "CoolantCorrelation", "nusselt_mcadams", "nusselt_dittus_boelter_1930",
    "nusselt_sieder_tate", "nusselt_gnielinski", "nusselt_taylor_tn_d4332",
    "ito_curvature_factor", "coolant_side_htc",
    # regime
    "CRITICAL_POINT", "supercritical_regime",
    "critical_heat_flux_applicability",
    # closure
    "energy_balance", "WallSolution", "solve_wall_temperature",
    "ThermalStation", "axial_thermal_march", "hot_spot",
]


# ===========================================================================
# SOURCES.  `access` is what I could actually open, not what exists.
# ===========================================================================

@dataclass(frozen=True)
class _Src:
    citation: str
    access: str          # OPENED | NOT_OPENED | OPENED_MIRROR
    note: str = ""

    def __str__(self) -> str:
        return self.citation


TAYLOR_TN_D_4332 = _Src(
    "Taylor, M. F., 'Correlation of Local Heat-Transfer Coefficients for "
    "Single-Phase Turbulent Flow of Hydrogen in Tubes with Temperature Ratios "
    "to 23', NASA TN D-4332, NASA Lewis Research Center, January 1968",
    "OPENED",
    "full PDF read from NTRS 19680004771. THE primary coolant-side source of "
    "this project.")

TAYLOR_TM_X_52437 = _Src(
    "Taylor, M. F., 'A Method of Predicting Heat Transfer Coefficients in the "
    "Cooling Passages of NERVA and Phoebus-2 Rocket Nozzles', NASA TM X-52437, "
    "1968, Eqs. (4) and (5)",
    "OPENED",
    "full PDF read from NTRS 19680014393. Curvature factor attributed by "
    "Taylor to Ito, H., J. Basic Eng. 81(2), June 1959, pp.123-134 (Ito "
    "itself NOT_OPENED).")

MCADAMS_1942 = _Src(
    "McAdams, W. H., 'Heat Transmission', 2nd ed., McGraw-Hill, 1942 -- the "
    "ACTUAL source of Nu = 0.023 Re^0.8 Pr^0.4",
    "NOT_OPENED",
    "known only through Winterton (1998). McAdams gave no clear reference for "
    "the modification.")

DITTUS_BOELTER_1930 = _Src(
    "Dittus, F. W. & Boelter, L. M. K., 'Heat Transfer in Automobile Radiators "
    "of the Tubular Type', University of California Publications in "
    "Engineering, Vol. 2, No. 13, pp.443-461, 1930",
    "NOT_OPENED",
    "what they actually published, per Winterton: 0.0243 Re^0.8 Pr^0.4 "
    "(heating) and 0.0265 Re^0.8 Pr^0.3 (cooling). NOT 0.023.")

WINTERTON_1998 = _Src(
    "Winterton, R. H. S., 'Where did the Dittus and Boelter equation come "
    "from?', Int. J. Heat Mass Transfer 41(4-5), 1998, pp.809-810",
    "OPENED_MIRROR",
    "publisher PDF is robots-disallowed on ScienceDirect; two independent web "
    "mirrors of the full note were read and agree on every number. They are "
    "copies of ONE document and are counted as one source, not two.")

GNIELINSKI_1976 = _Src(
    "Gnielinski, V., 'New equations for heat and mass transfer in turbulent "
    "pipe and channel flow', International Chemical Engineering 16(2), 1976, "
    "pp.359-368",
    "NOT_OPENED",
    "no accessible full text. His 2013 IJHMT restatement also NOT_OPENED. "
    "Every range below is secondhand and the secondary sources DISAGREE.")

SIEDER_TATE_1936 = _Src(
    "Sieder, E. N. & Tate, G. E., 'Heat Transfer and Pressure Drop of Liquids "
    "in Tubes', Ind. Eng. Chem. 28(12), 1936, pp.1429-1435, "
    "DOI 10.1021/ie50324a027",
    "NOT_OPENED",
    "ACS Publications returns 403. Whether it was fitted on heating or "
    "cooling, and on which fluids, is NOT_ANSWERED.")

HAALAND_1983 = _Src(
    "Haaland, S. E., 'Simple and Explicit Formulas for the Friction Factor in "
    "Turbulent Pipe Flow', Trans. ASME J. Fluids Eng. 105(1), 1983, pp.89-90",
    "NOT_OPENED",
    "ASME returns 403; OSTI returns 502. The EQUATION FORM is settled by three "
    "independent agreeing renderings. Haaland's OWN stated accuracy is "
    "NOT_REPORTED -- the '~1.5%' this repository previously carried could not "
    "be traced to him and has been withdrawn.")

COLEBROOK_1939 = _Src(
    "Colebrook, C. F., 'Turbulent Flow in Pipes, with Particular Reference to "
    "the Transition Region Between the Smooth and Rough Pipe Laws', J. Inst. "
    "Civil Engineers 11(4), 1938-39, pp.133-156",
    "NOT_OPENED",
    "ICE Virtual Library returns 403. Exact printed form and stated range "
    "NOT_REPORTED.")

LOCKE_LANDRUM_2008 = _Src(
    "Locke, J. M. & Landrum, D. B., 'Study of Heat Transfer Correlations for "
    "Supercritical Hydrogen in Regenerative Cooling Channels', J. Propulsion "
    "and Power 24(1), 2008, pp.94-103, DOI 10.2514/1.22496",
    "OPENED",
    "seven correlations vs 2992 measured points. Best performer (Schacht & "
    "Quentmeyer, NASA TN D-7207) needs +/-56% to cover 95% of the data.")

PIORO_MOKRY_2011 = _Src(
    "Pioro, I. & Mokry, S., 'Thermophysical Properties at Critical and "
    "Supercritical Conditions', in Heat Transfer -- Theoretical Analysis, "
    "Experimental Investigations and Industrial Systems, InTech, 2011, Sec. 2",
    "OPENED",
    "verbatim: 'At supercritical pressures there is no liquid-vapour phase "
    "transition; therefore, there is no such phenomenon as Critical Heat Flux "
    "(CHF) or dryout. Only within a certain range of parameters a "
    "deteriorated heat transfer may occur.'")

LAMANNA_2025 = _Src(
    "Lamanna, G. & Steinhausen, C., 'Modelling of heat addition to "
    "near-critical and supercritical fluids', Scientific Reports, 2025, "
    "DOI 10.1038/s41598-025-32459-z",
    "OPENED",
    "verbatim: 'HTD denotes a diminished capacity of SCFs to absorb heat and "
    "it is typically encountered in regenerative cooling applications'.")

NO_SOURCE = "NO SOURCE -- this value has no primary-source basis"


# ===========================================================================
# 1. GAS-SIDE HEAT FLUX
# ===========================================================================
#
# The legacy hot-wall assumption, moved here from three anonymous module
# constants in `pipeline/orchestrator.py`.  It is NOT endorsed and NOT a
# default: nothing in this module calls it.  It exists so that the one
# production site that still depends on it names what it is depending on.
LEGACY_HOT_WALL_FRACTION = 0.22
LEGACY_HOT_WALL_OFFSET_K = 100.0
LEGACY_HOT_WALL_CAP_K = 900.0


def assumed_hot_wall_temperature(stagnation_temperature_K: float) -> Result:
    """The legacy hot-wall guess, QUARANTINED. Not a model. Not a default.

    ``T_wg = min(900, 0.22 T_0 + 100)`` had no source, no derivation and no
    uncertainty. It is preserved here NUMERICALLY UNCHANGED so that removing
    it from production does not silently move a number, but it now returns
    `INSUFFICIENT_EVIDENCE` at `E0` and says in its own result that it is not
    a wall-temperature solution.

    The honest replacement is `solve_wall_temperature`, which needs a
    coolant-side closure. Where the caller has one, use that. Where the caller
    does not, the correct answer is that the wall temperature is UNKNOWN --
    and this function returns a status that says exactly that, so a gate can
    refuse it.
    """
    if stagnation_temperature_K <= 0:
        return Result(float("nan"), "K", "CAN-THERM-TWG-ASSUMED",
                      status=Status.PHYSICALLY_INVALID,
                      notes="stagnation temperature must be > 0")
    value = min(LEGACY_HOT_WALL_CAP_K,
                LEGACY_HOT_WALL_FRACTION * stagnation_temperature_K
                + LEGACY_HOT_WALL_OFFSET_K)
    return Result(
        value, "K", "CAN-THERM-TWG-ASSUMED",
        status=Status.INSUFFICIENT_EVIDENCE,
        all_statuses=frozenset({Status.INSUFFICIENT_EVIDENCE,
                                Status.UNVALIDATED_ASSUMPTION}),
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE, validation=Validation.NONE,
        source=NO_SOURCE, source_locator="none",
        validity_domain=UNKNOWN_DOMAIN, uncertainty=NOT_REPORTED,
        assumptions=(Assumption(
            "hot_wall_temperature_K", round(value, 1), "K",
            f"legacy guess min({LEGACY_HOT_WALL_CAP_K}, "
            f"{LEGACY_HOT_WALL_FRACTION} T_0 + {LEGACY_HOT_WALL_OFFSET_K}); "
            f"NO SOURCE; not a solution of the wall energy balance",
            Status.INSUFFICIENT_EVIDENCE),),
        notes="THIS IS NOT A WALL-TEMPERATURE SOLUTION. It is the pre-2026-08-20 "
              "guess, preserved so its removal is traceable. Two other "
              "production sites used 0.40 T_0 instead, a disagreement of "
              "~450 K on the same quantity, which is larger than the effect "
              "of any correlation choice in this module.",
        inputs={"T0_K": stagnation_temperature_K})


def gas_side_heat_flux(*, film_coefficient_W_m2K: float,
                       recovery_temperature_K: float,
                       wall_temperature_K: float,
                       wall_temperature_source: str) -> Result:
    r"""Gas-side convective heat flux.

    .. math:: q = h_g\,(T_{aw} - T_{wg})

    `wall_temperature_source` is REQUIRED and is not decoration. If it begins
    with ``ASSUMED`` the result carries `INSUFFICIENT_EVIDENCE`, because a flux
    computed at a guessed wall temperature is a guess whatever the precision of
    `h_g`. There is no default: omitting the wall temperature is a TypeError,
    not a 0.22 T_0.
    """
    inputs = {"h_g_W_m2K": film_coefficient_W_m2K,
              "T_aw_K": recovery_temperature_K,
              "T_wg_K": wall_temperature_K,
              "T_wg_source": wall_temperature_source}
    if not wall_temperature_source:
        raise ValueError(
            "wall_temperature_source must say where T_wg came from. This "
            "argument exists because three production sites used three "
            "different unsourced wall temperatures.")
    if film_coefficient_W_m2K <= 0:
        return Result(float("nan"), "W/m^2", "CAN-THERM-QGAS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="film coefficient must be > 0")
    if recovery_temperature_K <= 0 or wall_temperature_K <= 0:
        return Result(float("nan"), "W/m^2", "CAN-THERM-QGAS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="temperatures must be > 0 K")

    statuses = [Status.PASS]
    notes: list[str] = []
    assumptions: list[Assumption] = []

    if wall_temperature_source.upper().startswith("ASSUMED"):
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        assumptions.append(Assumption(
            "wall_temperature_K", round(wall_temperature_K, 1), "K",
            f"NOT SOLVED: {wall_temperature_source}",
            Status.INSUFFICIENT_EVIDENCE))
        notes.append("T_wg was assumed, not solved, so this flux is an "
                     "assumption dressed as a calculation")

    q = film_coefficient_W_m2K * (recovery_temperature_K - wall_temperature_K)

    if wall_temperature_K >= recovery_temperature_K:
        statuses.append(Status.INVALID_MODEL_REGIME)
        notes.append(
            "T_wg >= T_aw, so heat flows INTO the gas. That state is real "
            "(a film-cooled or heated wall) but this chain models a cooled "
            "chamber wall and its correlations are not validated for it.")

    return Result(
        q, "W/m^2", "CAN-THERM-QGAS", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="Newton's law of cooling with the adiabatic-wall (recovery) "
               "temperature as the driving potential",
        source_locator="definitional",
        validity_domain="T_wg < T_aw",
        uncertainty="inherits the uncertainty of h_g and of T_aw; the Bartz "
                    "primary reports NONE",
        assumptions=tuple(assumptions), notes="; ".join(notes), inputs=inputs)


# ===========================================================================
# 2. WALL CONDUCTION
# ===========================================================================

class WallGeometry(str, Enum):
    CYLINDRICAL = "CYLINDRICAL"
    PLANAR = "PLANAR"


def planar_wall_resistance(*, thickness_m: float,
                           conductivity_W_mK: float) -> Result:
    r"""Planar conduction resistance per unit area, :math:`R = t/k`."""
    if thickness_m <= 0 or conductivity_W_mK <= 0:
        return Result(float("nan"), "m^2.K/W", "CAN-THERM-RWALL-PLANAR",
                      status=Status.PHYSICALLY_INVALID,
                      notes="thickness and conductivity must be > 0",
                      inputs={"t_m": thickness_m, "k_W_mK": conductivity_W_mK})
    return Result(thickness_m / conductivity_W_mK, "m^2.K/W",
                  "CAN-THERM-RWALL-PLANAR", status=Status.PASS,
                  evidence_level=EvidenceLevel.E1,
                  verification=Verification.ALGEBRAIC,
                  source="one-dimensional steady Fourier conduction through a "
                         "slab",
                  source_locator="definitional",
                  validity_domain="plane wall, or a curved wall where t << r",
                  inputs={"t_m": thickness_m, "k_W_mK": conductivity_W_mK})


def cylindrical_wall_resistance(*, inner_radius_m: float, thickness_m: float,
                                conductivity_W_mK: float) -> Result:
    r"""Cylindrical conduction resistance REFERRED TO THE INNER (HOT) SURFACE.

    .. math:: R_i = \frac{r_i \ln(r_o/r_i)}{k}

    Referring it to the inner surface is what makes it directly addable to
    :math:`1/h_g`, which is also an inner-surface resistance. Mixing an
    inner-surface :math:`1/h_g` with an outer-surface coolant resistance
    without an area ratio is a silent error of :math:`r_o/r_i`; that is why
    the reference surface is in the name and in the returned notes.
    """
    inputs = {"r_i_m": inner_radius_m, "t_m": thickness_m,
              "k_W_mK": conductivity_W_mK}
    if inner_radius_m <= 0 or thickness_m <= 0 or conductivity_W_mK <= 0:
        return Result(float("nan"), "m^2.K/W", "CAN-THERM-RWALL-CYL",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="radius, thickness and conductivity must be > 0")
    r_o = inner_radius_m + thickness_m
    r = inner_radius_m * math.log(r_o / inner_radius_m) / conductivity_W_mK
    return Result(r, "m^2.K/W", "CAN-THERM-RWALL-CYL", status=Status.PASS,
                  evidence_level=EvidenceLevel.E1,
                  verification=Verification.ALGEBRAIC,
                  source="one-dimensional steady Fourier conduction through a "
                         "cylindrical shell",
                  source_locator="definitional",
                  validity_domain="axisymmetric shell, radial conduction only",
                  notes="referred to the INNER surface, so it adds directly to "
                        "1/h_g",
                  inputs={**inputs, "r_o_m": r_o})


def wall_resistance(*, thickness_m: float, conductivity_W_mK: float,
                    inner_radius_m: float | None = None,
                    geometry: WallGeometry = WallGeometry.CYLINDRICAL
                    ) -> Result:
    """Wall resistance with the planar approximation explicitly BOUNDED.

    A combustor liner is a cylindrical shell. Using the planar form on it is
    permitted -- it is often the right engineering call -- but it is never
    silent here. The exact ratio is

        R_planar / R_cyl = x / ln(1 + x),   x = t / r_i

    which is > 1 for all x > 0, so the planar form OVERSTATES the resistance
    and therefore UNDERSTATES the heat flux and the wall temperature drop.
    The leading term is x/2, i.e. roughly ``t / (2 r_i)`` in relative error.
    That number is computed and returned, and if it exceeds 5% the result
    carries `UNVALIDATED_ASSUMPTION`.
    """
    if geometry is WallGeometry.CYLINDRICAL:
        if inner_radius_m is None:
            return Result(float("nan"), "m^2.K/W", "CAN-THERM-RWALL",
                          status=Status.INSUFFICIENT_EVIDENCE,
                          notes="cylindrical resistance needs the inner "
                                "radius. Silently falling back to the planar "
                                "form is the error this function exists to "
                                "prevent.",
                          inputs={"t_m": thickness_m})
        return cylindrical_wall_resistance(inner_radius_m=inner_radius_m,
                                           thickness_m=thickness_m,
                                           conductivity_W_mK=conductivity_W_mK)

    planar = planar_wall_resistance(thickness_m=thickness_m,
                                    conductivity_W_mK=conductivity_W_mK)
    if planar.status is Status.PHYSICALLY_INVALID or inner_radius_m is None:
        return planar

    exact = cylindrical_wall_resistance(inner_radius_m=inner_radius_m,
                                        thickness_m=thickness_m,
                                        conductivity_W_mK=conductivity_W_mK)
    if exact.status is Status.PHYSICALLY_INVALID:
        return planar
    err = float(planar) / float(exact) - 1.0
    statuses = [Status.PASS]
    if abs(err) > 0.05:
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
    return Result(
        float(planar), "m^2.K/W", "CAN-THERM-RWALL-PLANAR-ON-CYLINDER",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        source="planar Fourier conduction applied to a cylindrical shell",
        source_locator="definitional; the error bound is exact, not estimated",
        validity_domain=f"t/r_i = {thickness_m / inner_radius_m:.4f}",
        uncertainty=f"EXACT: overstates the wall resistance by "
                    f"{err * 100:.2f}% relative to the cylindrical form, and "
                    f"therefore understates q and the wall temperature drop by "
                    f"the same factor",
        assumptions=(Assumption(
            "wall_geometry", "PLANAR", "-",
            f"planar form used on a cylindrical liner; exact error "
            f"{err * 100:.2f}%",
            Status.UNVALIDATED_ASSUMPTION if abs(err) > 0.05 else Status.PASS),),
        notes=f"cylindrical value is {float(exact):.6e} m^2.K/W",
        inputs={"t_m": thickness_m, "k_W_mK": conductivity_W_mK,
                "r_i_m": inner_radius_m, "relative_error": err})


def mean_wall_conductivity(*, material_property: str, hot_side_K: float,
                           cold_side_K: float, samples: int = 9,
                           restrict_to_valid_range: bool = False) -> Result:
    r"""Integral-mean conductivity across the wall.

    .. math:: \bar k = \frac{1}{T_h - T_c}\int_{T_c}^{T_h} k(T)\,dT

    This is the conductivity that makes :math:`q = \bar k (T_h - T_c)/t`
    EXACT for a temperature-dependent :math:`k`, not an approximation --
    Fourier's law integrates directly. Evaluating :math:`k` at a single
    temperature is the approximation, and for GRCop-84 across a 400 K wall
    drop it is worth several percent.

    If ANY sample temperature falls outside the source's stated range, the
    result is `OUT_OF_CORRELATION_RANGE` and the value is NaN. No
    extrapolation. That is the caller's cue to reduce the span or acquire
    data, not to fill in a number.

    `restrict_to_valid_range` narrows the integration to the part of the wall
    the source actually covers and returns the mean over THAT interval,
    permanently flagged `OUT_OF_CORRELATION_RANGE`. It is not a softening of
    the rule: nothing is extrapolated, and the returned notes state exactly
    which fraction of the wall had data. The iterative wall solver needs it
    because its intermediate iterates wander outside the range even when the
    converged answer does not -- refusing there would make a solvable problem
    unsolvable for a purely numerical reason.
    """
    from kryptonis.materials.temperature_dependent import (
        TEMPERATURE_DEPENDENT, evaluate_property)

    inputs = {"material_property": material_property,
              "T_hot_K": hot_side_K, "T_cold_K": cold_side_K}
    if hot_side_K <= 0 or cold_side_K <= 0:
        return Result(float("nan"), "W/(m.K)", "CAN-THERM-KBAR",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="temperatures must be > 0 K")
    if samples < 3 or samples % 2 == 0:
        samples = max(3, samples + 1)

    lo, hi = sorted((cold_side_K, hot_side_K))
    restricted_note = ""
    forced: list[Status] = []
    if restrict_to_valid_range:
        curve = TEMPERATURE_DEPENDENT.get(material_property)
        if curve is not None:
            vlo, vhi = curve.valid_range_K
            span = hi - lo
            nlo, nhi = max(lo, vlo), min(hi, vhi)
            if nhi <= nlo:
                return Result(
                    float("nan"), curve.units, "CAN-THERM-KBAR",
                    status=Status.OUT_OF_CORRELATION_RANGE,
                    all_statuses=frozenset({Status.OUT_OF_CORRELATION_RANGE}),
                    evidence_level=EvidenceLevel.E2,
                    source=curve.source_document,
                    source_locator=curve.source_locator,
                    validity_domain=f"{vlo:.0f}-{vhi:.0f} K",
                    notes=f"the whole wall [{lo:.1f}, {hi:.1f}] K lies outside "
                          f"the source's {vlo:.0f}-{vhi:.0f} K range. NO "
                          f"EXTRAPOLATION IS PERFORMED.",
                    inputs=inputs)
            if nlo > lo + 1e-9 or nhi < hi - 1e-9:
                covered = (nhi - nlo) / span if span > 0 else 1.0
                restricted_note = (
                    f"integration RESTRICTED to [{nlo:.1f}, {nhi:.1f}] K, the "
                    f"part of the wall the source covers ({covered * 100:.0f}% "
                    f"of the {lo:.1f}-{hi:.1f} K span). The remainder is "
                    f"NOT extrapolated and NOT included.")
                forced.append(Status.OUT_OF_CORRELATION_RANGE)
                lo, hi = nlo, nhi

    if hi - lo < 1e-9:
        return evaluate_property(material_property, lo)

    h = (hi - lo) / (samples - 1)
    vals: list[float] = []
    statuses: list[Status] = [Status.PASS, *forced]
    src = ""
    loc = ""
    unc = NOT_REPORTED
    dom = UNKNOWN_DOMAIN
    for i in range(samples):
        r = evaluate_property(material_property, lo + i * h)
        statuses.extend(r.all_statuses)
        src, loc = r.source, r.source_locator
        unc, dom = r.uncertainty, r.validity_domain
        if r.status in (Status.OUT_OF_CORRELATION_RANGE,
                        Status.INSUFFICIENT_EVIDENCE,
                        Status.PHYSICALLY_INVALID):
            return Result(
                float("nan"), "W/(m.K)", "CAN-THERM-KBAR",
                status=r.status, all_statuses=frozenset(statuses),
                evidence_level=r.evidence_level, source=src,
                source_locator=loc, validity_domain=dom, uncertainty=unc,
                notes=f"at T = {lo + i * h:.1f} K: {r.notes}", inputs=inputs)
        vals.append(float(r))

    # composite Simpson
    total = vals[0] + vals[-1]
    for i in range(1, samples - 1):
        total += (4.0 if i % 2 else 2.0) * vals[i]
    k_bar = (h / 3.0) * total / (hi - lo)

    return Result(
        k_bar, "W/(m.K)", "CAN-THERM-KBAR", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E2,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=src, source_locator=loc, validity_domain=dom, uncertainty=unc,
        notes="; ".join(x for x in (
            f"composite Simpson over {samples} points on "
            f"[{lo:.1f}, {hi:.1f}] K; k at the mean temperature would be a "
            f"different and less exact quantity", restricted_note) if x),
        inputs={**inputs, "samples": samples,
                "integrated_over_K": (lo, hi)})


# ===========================================================================
# 3. FRICTION FACTORS
# ===========================================================================

def friction_factor_haaland(*, reynolds: float,
                            relative_roughness: float) -> Result:
    r"""Haaland's explicit Darcy friction factor.

    .. math::
        \frac{1}{\sqrt f} = -1.8\log_{10}\!\left[
            \left(\frac{\varepsilon}{3.7 D}\right)^{1.11} + \frac{6.9}{Re}
        \right]

    THE FORM IS SETTLED, THE ACCURACY IS NOT. Haaland's paper is behind a 403.
    Three independent renderings (a Purdue course PDF, a Cengel & Cimbala
    textbook problem, and a third calculator source) all print the roughness
    term as :math:`(\varepsilon/3.7D)^{1.11}` -- the exponent applies to the
    WHOLE ratio, not to :math:`\varepsilon/D` alone. Two further renderings
    disagree and are demonstrably corrupted (one drops the 1.11 entirely).

    Haaland's OWN stated maximum deviation from Colebrook is `NOT_REPORTED`.
    This repository previously carried "accuracy within ~1.5%" attributed to
    Haaland; that figure could not be traced to him from any accessible source
    and has been WITHDRAWN. Third-party evaluations exist (Jaric et al.,
    FME Transactions: max +1.420% / -1.314% over Re 4e3-1e8) and are cited
    AS third-party, not as Haaland's claim.
    """
    inputs = {"Re": reynolds, "eps_over_D": relative_roughness}
    if reynolds <= 0:
        return Result(float("nan"), "-", "CAN-THERM-F-HAALAND",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Reynolds number must be > 0")
    if relative_roughness < 0:
        return Result(float("nan"), "-", "CAN-THERM-F-HAALAND",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="relative roughness cannot be negative")
    statuses = [Status.PASS]
    notes: list[str] = []
    if reynolds < 4000.0:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append("Re < 4000: this is a TURBULENT correlation and no "
                     "laminar fallback is substituted here, because silently "
                     "returning 64/Re hides the regime change from the caller")
    inner = (relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds
    if inner <= 0:
        return Result(float("nan"), "-", "CAN-THERM-F-HAALAND",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="argument of the logarithm is non-positive")
    f = (-1.8 * math.log10(inner)) ** -2
    return Result(
        f, "-", "CAN-THERM-F-HAALAND", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(HAALAND_1983), source_locator="p.89; PRIMARY NOT_OPENED "
                                                 "(ASME 403)",
        validity_domain="Haaland states NONE that could be retrieved; a "
                        "third-party review evaluated it over Re 4e3-1e8, "
                        "eps/D 1e-6 to 5e-2",
        uncertainty="Haaland's own claim: NOT_REPORTED. Third-party (Jaric "
                    "et al., FME Trans.): max +1.420% / -1.314% vs Colebrook. "
                    "Do not attribute a figure to Haaland.",
        notes="; ".join(notes), inputs=inputs)


def friction_factor_colebrook(*, reynolds: float, relative_roughness: float,
                              tol: float = 1e-10,
                              max_iter: int = 200) -> Result:
    r"""Colebrook-White implicit friction factor, solved by fixed-point.

    .. math::
        \frac{1}{\sqrt f} = -2\log_{10}\!\left(
            \frac{\varepsilon}{3.7D} + \frac{2.51}{Re\sqrt f}\right)

    The primary (J. Inst. Civil Engineers, 1938-39) is behind a 403 and could
    not be opened, so the exact printed form and Colebrook's own stated range
    are `NOT_REPORTED` and this implementation is `E1` on algebra only. It is
    seeded from Haaland, which is what makes it converge in a handful of
    iterations.
    """
    inputs = {"Re": reynolds, "eps_over_D": relative_roughness}
    if reynolds <= 0 or relative_roughness < 0:
        return Result(float("nan"), "-", "CAN-THERM-F-COLEBROOK",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re must be > 0 and roughness >= 0")
    statuses = [Status.PASS]
    if reynolds < 4000.0:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)

    seed = friction_factor_haaland(reynolds=reynolds,
                                   relative_roughness=relative_roughness)
    f = float(seed) if seed.usable and float(seed) > 0 else 0.02
    converged = False
    for _ in range(max_iter):
        rhs = -2.0 * math.log10(relative_roughness / 3.7
                                + 2.51 / (reynolds * math.sqrt(f)))
        if rhs <= 0:
            break
        f_new = rhs ** -2
        if abs(f_new - f) < tol:
            f = f_new
            converged = True
            break
        f = f_new
    if not converged:
        statuses.append(Status.FAIL)
    return Result(
        f, "-", "CAN-THERM-F-COLEBROOK", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(COLEBROOK_1939),
        source_locator="PRIMARY NOT_OPENED (ICE Virtual Library 403)",
        validity_domain=UNKNOWN_DOMAIN + " -- Colebrook's own stated range "
                                         "could not be retrieved",
        uncertainty=NOT_REPORTED,
        notes="" if converged else "fixed-point iteration did not converge",
        inputs={**inputs, "converged": converged})


# ===========================================================================
# 4. COOLANT-SIDE NUSSELT CORRELATIONS
# ===========================================================================
#
# Every one of these is a SEPARATE named function. There is no
# `nusselt(kind=...)` dispatcher that quietly picks one, because the choice
# among them is the single largest modelling decision on the coolant side and
# it must appear at the call site.

class CoolantCorrelation(str, Enum):
    """Which coolant-side correlation. There is deliberately NO default."""

    #: Nu = 0.023 Re^0.8 Pr^0.4. Correctly attributed to McAdams (1942).
    MCADAMS_1942 = "MCADAMS_1942"
    #: What Dittus & Boelter actually published in 1930.
    DITTUS_BOELTER_1930 = "DITTUS_BOELTER_1930"
    SIEDER_TATE_1936 = "SIEDER_TATE_1936"
    GNIELINSKI_1976 = "GNIELINSKI_1976"
    #: The only one with a primary source this project has actually opened.
    TAYLOR_TN_D_4332 = "TAYLOR_TN_D_4332"


def nusselt_mcadams(*, reynolds: float, prandtl: float,
                    heating: bool = True) -> Result:
    r"""``Nu = 0.023 Re^0.8 Pr^n``, n = 0.4 heating / 0.3 cooling.

    **ATTRIBUTION CORRECTED.** This is McAdams (1942), not Dittus & Boelter
    (1930). Winterton (1998) establishes that the 1930 paper published
    ``0.0243 Re^0.8 Pr^0.4`` (heating) and ``0.0265 Re^0.8 Pr^0.3`` (cooling),
    and that McAdams simplified the coefficient without giving a clear
    reference. Winterton, verbatim: *"It would be better to say the Dittus and
    Boelter equation, as introduced by McAdams."*

    The commonly quoted range ``2500 < Re < 1e5, 0.7 < Pr < 120`` is
    TEXTBOOK-LINEAGE and is NOT traceable to either original. It is applied
    here as a range check and labelled as what it is.
    """
    inputs = {"Re": reynolds, "Pr": prandtl, "heating": heating}
    if reynolds <= 0 or prandtl <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-MCADAMS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and Pr must be > 0")
    n = 0.4 if heating else 0.3
    nu = 0.023 * reynolds ** 0.8 * prandtl ** n
    statuses = [Status.PASS,
                # The attribution is settled; the RANGE is not sourced at all.
                Status.UNVALIDATED_ASSUMPTION]
    notes = ["validity range is textbook-lineage, NOT traceable to McAdams "
             "1942 or to Dittus & Boelter 1930"]
    if not (10_000.0 <= reynolds <= 1.0e5):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside the quoted 1e4-1e5")
    if not (0.7 <= prandtl <= 120.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Pr = {prandtl:.3f} outside the quoted 0.7-120")
    return Result(
        nu, "-", "CAN-THERM-NU-MCADAMS", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(MCADAMS_1942),
        source_locator="PRIMARY NOT_OPENED; attribution established by "
                       + str(WINTERTON_1998),
        validity_domain="quoted 1e4 < Re < 1e5, 0.7 < Pr < 120 -- SOURCE "
                        "UNTRACEABLE",
        uncertainty=NOT_REPORTED,
        assumptions=(Assumption("validity_range", "textbook-lineage", "-",
                                "no primary source states it",
                                Status.UNVALIDATED_ASSUMPTION),),
        notes="; ".join(notes), inputs=inputs)


def nusselt_dittus_boelter_1930(*, reynolds: float, prandtl: float,
                                heating: bool = True) -> Result:
    """What Dittus & Boelter ACTUALLY published: 0.0243/0.4 or 0.0265/0.3.

    Provided so the misattribution is falsifiable in code rather than only in
    a comment. Note the coefficients differ between heating and cooling; the
    0.023 form collapses them into one, which is precisely the simplification
    McAdams made and did not reference.
    """
    inputs = {"Re": reynolds, "Pr": prandtl, "heating": heating}
    if reynolds <= 0 or prandtl <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-DB1930",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and Pr must be > 0")
    coeff, n = (0.0243, 0.4) if heating else (0.0265, 0.3)
    nu = coeff * reynolds ** 0.8 * prandtl ** n
    return Result(
        nu, "-", "CAN-THERM-NU-DB1930", status=Status.INSUFFICIENT_EVIDENCE,
        all_statuses=frozenset({Status.INSUFFICIENT_EVIDENCE}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(DITTUS_BOELTER_1930),
        source_locator="PRIMARY NOT_OPENED; coefficients reported by "
                       + str(WINTERTON_1998),
        validity_domain=UNKNOWN_DOMAIN + " -- the 1930 paper's own stated "
                                         "range is NOT_REPORTED",
        uncertainty=NOT_REPORTED,
        notes="INSUFFICIENT_EVIDENCE because the 1930 paper itself was not "
              "opened: these coefficients are reported by Winterton, whose "
              "publisher version was also not opened (two agreeing mirrors "
              "of one document count as one source).",
        inputs={**inputs, "coefficient": coeff, "pr_exponent": n})


def nusselt_sieder_tate(*, reynolds: float, prandtl: float,
                        bulk_viscosity_Pa_s: float,
                        wall_viscosity_Pa_s: float) -> Result:
    r"""``Nu = 0.027 Re^0.8 Pr^(1/3) (mu_b/mu_w)^0.14``.

    The primary is behind a 403. Whether it was fitted on HEATING or COOLING,
    and on which fluids, is `NOT_ANSWERED` -- and that matters, because the
    viscosity-ratio exponent is the whole point of the correlation. Secondhand
    range: Re >= 1e4, 0.7 <= Pr <= 16700.
    """
    inputs = {"Re": reynolds, "Pr": prandtl,
              "mu_b_Pa_s": bulk_viscosity_Pa_s,
              "mu_w_Pa_s": wall_viscosity_Pa_s}
    if (reynolds <= 0 or prandtl <= 0 or bulk_viscosity_Pa_s <= 0
            or wall_viscosity_Pa_s <= 0):
        return Result(float("nan"), "-", "CAN-THERM-NU-SIEDERTATE",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re, Pr and both viscosities must be > 0")
    nu = (0.027 * reynolds ** 0.8 * prandtl ** (1.0 / 3.0)
          * (bulk_viscosity_Pa_s / wall_viscosity_Pa_s) ** 0.14)
    statuses = [Status.PASS]
    notes: list[str] = []
    if reynolds < 1.0e4:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} below the secondhand 1e4 floor")
    if not (0.7 <= prandtl <= 16700.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Pr = {prandtl:.3f} outside the secondhand 0.7-16700")
    return Result(
        nu, "-", "CAN-THERM-NU-SIEDERTATE", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(SIEDER_TATE_1936),
        source_locator="PRIMARY NOT_OPENED (ACS 403)",
        validity_domain="SECONDHAND: Re >= 1e4, 0.7 <= Pr <= 16700",
        uncertainty="SECONDHAND: '10% deviation' quoted for water; not "
                    "verified against the original",
        notes="; ".join(notes), inputs=inputs)


def nusselt_gnielinski(*, reynolds: float, prandtl: float,
                       friction_factor: float | None = None,
                       diameter_over_length: float | None = None) -> Result:
    r"""``Nu = (f/8)(Re-1000)Pr / [1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)]``.

    **THE FRICTION FACTOR IS AN UNRESOLVED CONFLICT AND IS NOT DEFAULTED.**
    Gnielinski's own paper could not be opened. Three secondary sources
    attribute the ``f`` he specified three different ways:

        (1.82 log10 Re - 1.64)^-2   attributed to Filonenko by two sources
        (1.82 log10 Re - 1.64)^-2   attributed to Petukhov by a third
        (1.80 log10 Re - 1.50)^-2   attributed to Konakov by the VDI Heat
                                    Atlas -- Gnielinski's own venue

    The first two are the same algebra with different names; the third is a
    DIFFERENT NUMBER. Averaging them or picking the popular one would be
    exactly the majority-vote this project forbids. So `friction_factor` is a
    REQUIRED argument: supply the one you can defend. Omitting it returns
    `CONFLICT_UNRESOLVED` with both candidates and the spread between them.

    The quoted Re and Pr ranges also conflict across secondary sources
    (Re upper bound 1e6 / 1.5e6 / 5e9; Pr upper bound 200 / 1000 / 2000). The
    NARROWEST of each is applied, because a range check that is wider than the
    evidence is not a check.
    """
    inputs = {"Re": reynolds, "Pr": prandtl, "f": friction_factor,
              "d_over_L": diameter_over_length}
    if reynolds <= 0 or prandtl <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-GNIELINSKI",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and Pr must be > 0")

    if friction_factor is None:
        f_filonenko = (1.82 * math.log10(reynolds) - 1.64) ** -2
        f_konakov = (1.80 * math.log10(reynolds) - 1.50) ** -2
        spread = abs(f_konakov / f_filonenko - 1.0) * 100.0
        return Result(
            float("nan"), "-", "CAN-THERM-NU-GNIELINSKI",
            status=Status.CONFLICT_UNRESOLVED,
            all_statuses=frozenset({Status.CONFLICT_UNRESOLVED}),
            evidence_level=EvidenceLevel.E0,
            source=str(GNIELINSKI_1976),
            source_locator="PRIMARY NOT_OPENED; f attribution conflicts "
                           "across secondary sources",
            validity_domain=UNKNOWN_DOMAIN,
            uncertainty=NOT_REPORTED,
            notes=f"f not supplied. Filonenko/Petukhov form gives "
                  f"{f_filonenko:.6f}; Konakov form (VDI Heat Atlas, "
                  f"Gnielinski's own venue) gives {f_konakov:.6f}. Spread "
                  f"{spread:.2f}% on f, ~{spread:.2f}% on Nu. Choose one "
                  f"explicitly; this function will not choose for you.",
            inputs={**inputs, "f_filonenko": f_filonenko,
                    "f_konakov": f_konakov, "spread_percent": spread})

    if friction_factor <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-GNIELINSKI",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="friction factor must be > 0")

    statuses = [Status.PASS]
    notes: list[str] = []
    f8 = friction_factor / 8.0
    denom = 1.0 + 12.7 * math.sqrt(f8) * (prandtl ** (2.0 / 3.0) - 1.0)
    if denom <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-GNIELINSKI",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="denominator non-positive; Pr is too small for "
                            "this form at this friction factor")
    nu = f8 * (reynolds - 1000.0) * prandtl / denom

    if diameter_over_length is not None:
        if diameter_over_length < 0:
            return Result(float("nan"), "-", "CAN-THERM-NU-GNIELINSKI",
                          status=Status.PHYSICALLY_INVALID, inputs=inputs,
                          notes="d/L cannot be negative")
        nu *= 1.0 + diameter_over_length ** (2.0 / 3.0)
        statuses.append(Status.CONFLICT_UNRESOLVED)
        notes.append("the entrance factor [1 + (d/L)^(2/3)] is applied as a "
                     "MULTIPLIER, which is how two of three secondary "
                     "renderings print it; a third places it inside the "
                     "denominator. Unresolved without the primary.")

    if reynolds < 3000.0 or reynolds > 1.0e6:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside the NARROWEST secondhand "
                     f"range 3e3-1e6 (sources also quote 2.3e3-1e6, "
                     f"3e3-1.5e6 and 3e3-5e9)")
    if prandtl < 0.5 or prandtl > 200.0:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Pr = {prandtl:.3f} outside the NARROWEST secondhand "
                     f"range 0.5-200 (sources also quote 0.5-1000 and "
                     f"0.5-2000)")

    return Result(
        nu, "-", "CAN-THERM-NU-GNIELINSKI", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=str(GNIELINSKI_1976),
        source_locator="PRIMARY NOT_OPENED; form and ranges are SECONDHAND "
                       "and the secondary sources disagree",
        validity_domain="NARROWEST secondhand: 3e3 <= Re <= 1e6, "
                        "0.5 <= Pr <= 200",
        uncertainty="SECONDHAND: '~90% of ~800 experimental values within "
                    "+/-20%'; not verified against the original",
        notes="; ".join(notes), inputs=inputs)


def nusselt_taylor_tn_d4332(*, reynolds_bulk: float, prandtl_bulk: float,
                            wall_temperature_K: float,
                            bulk_temperature_K: float,
                            axial_over_diameter: float,
                            fluid: str = "H2") -> Result:
    r"""Taylor, NASA TN D-4332 (1968), THE ONE PRIMARY SOURCE ACTUALLY OPENED.

    .. math::
        Nu_b = 0.023\,Re_b^{0.8}\,Pr_b^{0.4}\,
               \left(\frac{T_s}{T_b}\right)^{-\left(0.57 - \frac{1.59}{x/D}\right)}

    Properties at BULK temperature. Fitted on **hydrogen**.

    Taylor's own stated envelope, from the report:

        T_s/T_b   1.1 to 23
        x/D       2 to 252
        Re_b      7 500 to 13 800 000
        T_s       114 to 5 630 degR  (63 to 3 130 K)
        q         0.036 to 27.6 Btu/(s.in^2)  =  0.059 to 45.7 MW/m^2

    Stated accuracy: *"87 percent of the 3674 calculated heat-transfer
    coefficients deviated less than +/-25 percent from the measured values"*
    -- in region 1.

    AND ITS AUTHOR'S OWN NEAR-CRITICAL DISCLAIMER, which is the most valuable
    sentence in the coolant-side literature this project has: the equation
    *"does not predict heat-transfer coefficients with acceptable accuracy in
    the near-critical pressure and temperature region"*, where only **40
    percent** of predictions fell within +/-25 percent.

    Applying it to methane is out of its fluid domain. That is flagged, not
    forbidden -- but it is never silent.
    """
    inputs = {"Re_b": reynolds_bulk, "Pr_b": prandtl_bulk,
              "T_s_K": wall_temperature_K, "T_b_K": bulk_temperature_K,
              "x_over_D": axial_over_diameter, "fluid": fluid}
    if reynolds_bulk <= 0 or prandtl_bulk <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-TAYLOR",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and Pr must be > 0")
    if wall_temperature_K <= 0 or bulk_temperature_K <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-TAYLOR",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="temperatures must be > 0 K")
    if axial_over_diameter <= 0:
        return Result(float("nan"), "-", "CAN-THERM-NU-TAYLOR",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="x/D must be > 0; the exponent contains 1.59/(x/D)")

    statuses = [Status.PASS]
    notes: list[str] = []

    ratio = wall_temperature_K / bulk_temperature_K
    exponent = -(0.57 - 1.59 / axial_over_diameter)
    nu = (0.023 * reynolds_bulk ** 0.8 * prandtl_bulk ** 0.4
          * ratio ** exponent)

    if fluid.upper().replace("-", "").replace("_", "") not in ("H2", "GH2",
                                                              "LH2",
                                                              "HYDROGEN"):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(
            f"fluid = {fluid!r}. TN D-4332 is fitted on HYDROGEN. Applying it "
            f"to another coolant is outside its fluid domain; the report "
            f"contains no methane data.")
    if not (1.1 <= ratio <= 23.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"T_s/T_b = {ratio:.2f} outside Taylor's stated 1.1-23")
    if not (2.0 <= axial_over_diameter <= 252.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"x/D = {axial_over_diameter:.2f} outside Taylor's "
                     f"stated 2-252")
    if not (7500.0 <= reynolds_bulk <= 1.38e7):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re_b = {reynolds_bulk:.3e} outside Taylor's stated "
                     f"7.5e3-1.38e7")
    if not (63.0 <= wall_temperature_K <= 3130.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"T_s = {wall_temperature_K:.0f} K outside Taylor's "
                     f"stated 63-3130 K")

    return Result(
        nu, "-", "CAN-THERM-NU-TAYLOR", status=worst(statuses),
        all_statuses=frozenset(statuses),
        # E1, NOT E3. Taylor validated TAYLOR'S correlation against 3674 of
        # HIS OWN measurements. This implementation has been compared to no
        # measurement at all -- the data are in figures that were never
        # digitised. A source's experimental pedigree is not this code's
        # evidence level. It is recorded in `uncertainty` and `source_locator`
        # instead, which is where it belongs.
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source=str(TAYLOR_TN_D_4332),
        source_locator="NTRS 19680004771, full PDF OPENED; equation number "
                       "and page NOT_RECORDED from the fetch",
        validity_domain="HYDROGEN; T_s/T_b 1.1-23; x/D 2-252; Re_b "
                        "7.5e3-1.38e7; T_s 63-3130 K; q 0.059-45.7 MW/m^2. "
                        "EXCLUDED by the author near the critical point.",
        uncertainty="87% of 3674 points within +/-25% (region 1). Only 40% "
                    "within +/-25% in the near-critical region, which Taylor "
                    "excludes explicitly.",
        notes="; ".join(notes),
        inputs={**inputs, "T_s_over_T_b": ratio, "temperature_exponent": exponent})


def ito_curvature_factor(*, reynolds_bulk: float, hydraulic_radius_m: float,
                         radius_of_curvature_m: float,
                         side: str) -> Result:
    r"""Curvature enhancement/suppression factor for a curved cooling passage.

    .. math:: \phi_{curv} = \left[Re_b\left(\frac{r}{R}\right)^2\right]^{\pm 0.05}

    ``+`` on the CONCAVE (outer) wall, ``-`` on the CONVEX (inner) wall.

    This is Taylor, NASA TM X-52437 (1968), Eqs. (4) and (5) -- opened -- who
    attributes the factor to Ito, H., *"Friction Factors for Turbulent Flow in
    Curved Pipes"*, J. Basic Eng. 81(2), June 1959, pp.123-134 (Ito itself
    NOT_OPENED). It is not a Kryptonis invention and it is not a generic Dean
    correction: JAXA's own regenerative-cooling model (Negishi et al., EUCASS)
    uses this exact factor as its Eq. (14), citing the same Taylor report.

    Taylor states NO numerical bounds on Re, r/R or the factor itself --
    `validity_domain` is therefore UNKNOWN, and that is a real gap, not an
    oversight in this docstring. He states only that ~90% of predicted
    coefficients fell within +/-20% of experiment, validated against curved-tube
    hydrogen data (Aerojet-General Rep. 2551, 1963, NOT_OPENED).
    """
    inputs = {"Re_b": reynolds_bulk, "r_h_m": hydraulic_radius_m,
              "R_c_m": radius_of_curvature_m, "side": side}
    key = side.strip().upper()
    if key not in ("CONCAVE", "CONVEX"):
        return Result(float("nan"), "-", "CAN-THERM-CURV-ITO",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="side must be 'CONCAVE' or 'CONVEX'. There is no "
                            "sign-free version of this factor: Taylor's Eqs. "
                            "(4) and (5) differ ONLY in that sign.")
    if reynolds_bulk <= 0 or hydraulic_radius_m <= 0 or radius_of_curvature_m <= 0:
        return Result(float("nan"), "-", "CAN-THERM-CURV-ITO",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re, r and R must be > 0")
    sign = 1.0 if key == "CONCAVE" else -1.0
    arg = reynolds_bulk * (hydraulic_radius_m / radius_of_curvature_m) ** 2
    phi = arg ** (sign * 0.05)
    return Result(
        phi, "-", "CAN-THERM-CURV-ITO", status=Status.PASS,
        # E1 for the same reason as CAN-THERM-NU-TAYLOR: the SOURCE is
        # validated, this transcription of it is not.
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source=str(TAYLOR_TM_X_52437),
        source_locator="Eq. (4) concave / Eq. (5) convex; curvature factor "
                       "attributed by Taylor to Ito (1959)",
        validity_domain=UNKNOWN_DOMAIN + " -- Taylor states no numerical "
                                         "bounds on the curvature term",
        uncertainty="~90% of predicted coefficients within +/-20% of "
                    "experiment (whole method, not the factor alone)",
        notes=f"{key} wall, exponent {sign * 0.05:+.2f}; independently used "
              f"as Eq. (14) of JAXA's regenerative-cooling model",
        inputs={**inputs, "Re_(r/R)^2": arg})


def coolant_side_htc(*, correlation: CoolantCorrelation,
                     conductivity_W_mK: float, hydraulic_diameter_m: float,
                     reynolds: float, prandtl: float,
                     heating: bool = True,
                     bulk_viscosity_Pa_s: float | None = None,
                     wall_viscosity_Pa_s: float | None = None,
                     friction_factor: float | None = None,
                     diameter_over_length: float | None = None,
                     wall_temperature_K: float | None = None,
                     bulk_temperature_K: float | None = None,
                     axial_over_diameter: float | None = None,
                     fluid: str = NOT_REPORTED,
                     curvature_radius_m: float | None = None,
                     curvature_side: str | None = None) -> Result:
    r"""``h_c = Nu k / D_h``, with the correlation chosen EXPLICITLY.

    `correlation` has no default. Every previous coolant-side call in this
    repository silently used the 0.023 form; choosing it is now a visible act
    with a visible attribution.

    The Ito curvature factor is applied only when BOTH `curvature_radius_m`
    and `curvature_side` are given -- there is no "assume straight" and no
    "assume concave".
    """
    inputs = {"correlation": correlation.value, "k_W_mK": conductivity_W_mK,
              "D_h_m": hydraulic_diameter_m, "Re": reynolds, "Pr": prandtl}
    if conductivity_W_mK <= 0 or hydraulic_diameter_m <= 0:
        return Result(float("nan"), "W/(m^2.K)", "CAN-THERM-HC",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="conductivity and hydraulic diameter must be > 0")

    if correlation is CoolantCorrelation.MCADAMS_1942:
        nu_res = nusselt_mcadams(reynolds=reynolds, prandtl=prandtl,
                                 heating=heating)
    elif correlation is CoolantCorrelation.DITTUS_BOELTER_1930:
        nu_res = nusselt_dittus_boelter_1930(reynolds=reynolds,
                                             prandtl=prandtl, heating=heating)
    elif correlation is CoolantCorrelation.SIEDER_TATE_1936:
        if bulk_viscosity_Pa_s is None or wall_viscosity_Pa_s is None:
            return Result(float("nan"), "W/(m^2.K)", "CAN-THERM-HC",
                          status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                          notes="Sieder-Tate needs both bulk and wall "
                                "viscosity. Its viscosity ratio IS the "
                                "correlation; defaulting it to 1 turns it "
                                "into a differently-coefficiented McAdams.")
        nu_res = nusselt_sieder_tate(reynolds=reynolds, prandtl=prandtl,
                                     bulk_viscosity_Pa_s=bulk_viscosity_Pa_s,
                                     wall_viscosity_Pa_s=wall_viscosity_Pa_s)
    elif correlation is CoolantCorrelation.GNIELINSKI_1976:
        nu_res = nusselt_gnielinski(reynolds=reynolds, prandtl=prandtl,
                                    friction_factor=friction_factor,
                                    diameter_over_length=diameter_over_length)
    elif correlation is CoolantCorrelation.TAYLOR_TN_D_4332:
        if (wall_temperature_K is None or bulk_temperature_K is None
                or axial_over_diameter is None):
            return Result(float("nan"), "W/(m^2.K)", "CAN-THERM-HC",
                          status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                          notes="Taylor needs T_s, T_b and x/D. The "
                                "temperature-ratio term with its x/D-dependent "
                                "exponent is what distinguishes it from the "
                                "0.023 form.")
        nu_res = nusselt_taylor_tn_d4332(
            reynolds_bulk=reynolds, prandtl_bulk=prandtl,
            wall_temperature_K=wall_temperature_K,
            bulk_temperature_K=bulk_temperature_K,
            axial_over_diameter=axial_over_diameter, fluid=fluid)
    else:  # pragma: no cover - enum is exhaustive
        raise ValueError(f"unknown correlation {correlation!r}")

    if not nu_res.usable or math.isnan(float(nu_res)):
        return Result(float("nan"), "W/(m^2.K)", "CAN-THERM-HC",
                      status=nu_res.status,
                      all_statuses=nu_res.all_statuses,
                      evidence_level=nu_res.evidence_level,
                      source=nu_res.source,
                      source_locator=nu_res.source_locator,
                      validity_domain=nu_res.validity_domain,
                      uncertainty=nu_res.uncertainty,
                      notes=nu_res.notes, inputs={**inputs, **nu_res.inputs})

    statuses = list(nu_res.all_statuses)
    assumptions = list(nu_res.assumptions)
    notes = [nu_res.notes] if nu_res.notes else []
    nu = float(nu_res)

    if curvature_radius_m is not None and curvature_side is not None:
        curv = ito_curvature_factor(reynolds_bulk=reynolds,
                                    hydraulic_radius_m=hydraulic_diameter_m / 2.0,
                                    radius_of_curvature_m=curvature_radius_m,
                                    side=curvature_side)
        if not curv.usable:
            return Result(float("nan"), "W/(m^2.K)", "CAN-THERM-HC",
                          status=Status.PHYSICALLY_INVALID, inputs=inputs,
                          notes=curv.notes)
        nu *= float(curv)
        statuses.extend(curv.all_statuses)
        notes.append(f"Ito curvature factor {float(curv):.4f} applied "
                     f"({curvature_side})")
    elif curvature_radius_m is not None or curvature_side is not None:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append("curvature radius and side must BOTH be given; the sign "
                     "of the exponent is not inferable from the radius alone")

    h_c = nu * conductivity_W_mK / hydraulic_diameter_m
    return Result(
        h_c, "W/(m^2.K)", "CAN-THERM-HC", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=nu_res.evidence_level,
        verification=nu_res.verification, validation=nu_res.validation,
        source=nu_res.source, source_locator=nu_res.source_locator,
        validity_domain=nu_res.validity_domain,
        uncertainty=nu_res.uncertainty,
        assumptions=tuple(assumptions), notes="; ".join(n for n in notes if n),
        inputs={**inputs, "Nu": nu})


# ===========================================================================
# 5. REGIME: SUPERCRITICAL COOLANT, AND WHY CHF IS THE WRONG QUESTION
# ===========================================================================

#: Critical points. Widely tabulated; NIST Chemistry WebBook SRD 69 was NOT
#: opened in this session, so these carry LOCATOR_UNVERIFIED rather than a
#: page. They are used only for REGIME CLASSIFICATION -- no correlation
#: coefficient depends on them.
CRITICAL_POINT: dict[str, dict[str, float]] = {
    "CH4":  {"p_crit_Pa": 4.5992e6, "T_crit_K": 190.564},
    "H2":   {"p_crit_Pa": 1.2964e6, "T_crit_K": 33.145},
    "O2":   {"p_crit_Pa": 5.043e6,  "T_crit_K": 154.581},
}
_FLUID_ALIAS = {
    "CH4": "CH4", "LCH4": "CH4", "GCH4": "CH4", "METHANE": "CH4", "LNG": "CH4",
    "H2": "H2", "LH2": "H2", "GH2": "H2", "HYDROGEN": "H2",
    "O2": "O2", "LOX": "O2", "GOX": "O2", "OXYGEN": "O2",
}


def _canon_fluid(fluid: str) -> str | None:
    key = fluid.strip().upper().replace("-", "").replace("_", "").replace("/", "")
    return _FLUID_ALIAS.get(key)


def supercritical_regime(*, fluid: str, pressure_Pa: float,
                         bulk_temperature_K: float,
                         wall_temperature_K: float | None = None) -> Result:
    """Classify the coolant regime, and flag heat-transfer deterioration risk.

    Returns one of ``SUBCRITICAL`` / ``SUPERCRITICAL`` / ``NEAR_CRITICAL`` as
    the value. The status is what matters:

    * Above p_crit with the bulk below and the wall above T_crit, the coolant
      crosses the pseudo-critical (Widom) line inside the boundary layer, and
      that is the documented condition for **heat transfer deterioration**
      (Lamanna & Steinhausen 2025). Every bulk-property correlation in this
      module -- McAdams, Sieder-Tate, Gnielinski, Taylor -- is reported to
      MISPREDICT there. Locke & Landrum (2008) found bulk-reference
      correlations *"overpredict the heat transfer in the region close to the
      critical point"* and *"along the pseudocritical temperature line"*.
    * Taylor excludes this region in his own report.

    So the status is `INVALID_MODEL_REGIME`, not a correction factor. There is
    no traceable deterioration correlation in this project: the Pizzarelli
    papers that contain one are all behind 403s.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "T_wall_K": wall_temperature_K}
    key = _canon_fluid(fluid)
    if key is None:
        return Result("UNKNOWN", "-", "CAN-THERM-REGIME",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      evidence_level=EvidenceLevel.E0,
                      notes=f"no critical point recorded for {fluid!r}; the "
                            f"regime cannot be classified and must not be "
                            f"assumed subcritical")
    if pressure_Pa <= 0 or bulk_temperature_K <= 0:
        return Result("INVALID", "-", "CAN-THERM-REGIME",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="pressure and temperature must be > 0")

    p_c = CRITICAL_POINT[key]["p_crit_Pa"]
    t_c = CRITICAL_POINT[key]["T_crit_K"]
    p_r = pressure_Pa / p_c
    statuses = [Status.PASS]
    notes = [f"p/p_crit = {p_r:.3f}, T_crit = {t_c:.1f} K"]

    if pressure_Pa < p_c:
        label = "SUBCRITICAL"
        notes.append("below the critical pressure: a liquid-vapour phase "
                     "boundary exists, so boiling and CHF are meaningful "
                     "questions here")
    else:
        label = "SUPERCRITICAL"
        if 0.98 <= p_r <= 1.5:
            label = "NEAR_CRITICAL"
            statuses.append(Status.INVALID_MODEL_REGIME)
            notes.append("0.98 <= p/p_crit <= 1.5: property variation is at "
                         "its most violent and Taylor excludes this region in "
                         "his own report")
        crosses = (wall_temperature_K is not None
                   and bulk_temperature_K < t_c < wall_temperature_K)
        if crosses:
            statuses.append(Status.INVALID_MODEL_REGIME)
            notes.append(
                "the boundary layer spans the pseudo-critical temperature "
                "(T_bulk < T_crit < T_wall). This is the documented condition "
                "for HEAT TRANSFER DETERIORATION. Bulk-property correlations "
                "are reported to overpredict h here; no traceable "
                "deterioration correlation is available in this project.")

    return Result(
        label, "-", "CAN-THERM-REGIME", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.NONE, validation=Validation.NONE,
        source=f"{LOCKE_LANDRUM_2008}; {LAMANNA_2025}; {PIORO_MOKRY_2011}",
        source_locator="critical-point values are widely tabulated; NIST SRD "
                       "69 NOT_OPENED in this session (LOCATOR_UNVERIFIED). "
                       "They classify the regime only; no correlation "
                       "coefficient depends on them.",
        validity_domain="pure-component critical point; mixtures NOT_HANDLED",
        uncertainty="Locke & Landrum: the best of seven correlations needed "
                    "+/-56% to cover 95% of 2992 supercritical hydrogen "
                    "points. That is the realistic band, not +/-20%.",
        notes="; ".join(notes), inputs={**inputs, "p_over_p_crit": p_r})


def critical_heat_flux_applicability(*, fluid: str,
                                     pressure_Pa: float) -> Result:
    """Is a boiling-crisis (CHF/DNB) model even defined at this pressure?

    Above the critical pressure the answer is NO, and this is not a matter of
    accuracy. Pioro & Mokry (2011), verbatim:

        *"At supercritical pressures there is no liquid-vapour phase
        transition; therefore, there is no such phenomenon as Critical Heat
        Flux (CHF) or dryout. Only within a certain range of parameters a
        deteriorated heat transfer may occur."*

    A regenerative methane channel runs far above methane's 4.599 MPa critical
    pressure. A Zuber pool-boiling CHF margin computed there is a CATEGORY
    ERROR, not a conservative estimate -- Zuber's correlation is built from
    the latent heat, the surface tension and the density difference, all three
    of which cease to exist above the critical point.

    The analogous failure mode is HEAT TRANSFER DETERIORATION, and it is a
    different phenomenon with a temperature-based rather than a flux-based
    onset. Use `supercritical_regime` for that.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa}
    key = _canon_fluid(fluid)
    if key is None:
        return Result(False, "-", "CAN-THERM-CHF-APPLICABLE",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      evidence_level=EvidenceLevel.E0,
                      notes=f"no critical pressure recorded for {fluid!r}; "
                            f"applicability cannot be decided")
    if pressure_Pa <= 0:
        return Result(False, "-", "CAN-THERM-CHF-APPLICABLE",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="pressure must be > 0")
    p_c = CRITICAL_POINT[key]["p_crit_Pa"]
    if pressure_Pa >= p_c:
        return Result(
            False, "-", "CAN-THERM-CHF-APPLICABLE",
            status=Status.INVALID_MODEL_REGIME,
            all_statuses=frozenset({Status.INVALID_MODEL_REGIME}),
            evidence_level=EvidenceLevel.E1,
            verification=Verification.NONE,
            validation=Validation.NONE,
            source=str(PIORO_MOKRY_2011),
            source_locator="Section 2, 'Historical note'",
            validity_domain=f"p >= p_crit = {p_c / 1e6:.3f} MPa for {key}",
            uncertainty="not applicable -- this is a categorical result, not "
                        "a numerical one",
            notes=f"p = {pressure_Pa / 1e6:.3f} MPa >= p_crit "
                  f"{p_c / 1e6:.3f} MPa. THERE IS NO CRITICAL HEAT FLUX. A "
                  f"Zuber margin here is a category error: h_fg, sigma and "
                  f"(rho_l - rho_v) all vanish. The relevant failure mode is "
                  f"HEAT TRANSFER DETERIORATION.",
            inputs={**inputs, "p_crit_Pa": p_c})
    return Result(
        True, "-", "CAN-THERM-CHF-APPLICABLE", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        validation=Validation.NONE,
        source=str(PIORO_MOKRY_2011), source_locator="Section 2",
        validity_domain=f"p < p_crit = {p_c / 1e6:.3f} MPa for {key}",
        notes="a phase boundary exists, so a boiling-crisis model is defined. "
              "Whether the SPECIFIC Zuber pool-boiling correlation applies to "
              "a forced-convection channel is a separate question this "
              "function does not answer.",
        inputs={**inputs, "p_crit_Pa": p_c})


# ===========================================================================
# 6. COOLANT ENERGY BALANCE
# ===========================================================================

def energy_balance(*, heat_rate_W: float, mass_flow_kg_s: float,
                   inlet_enthalpy_J_kg: float | None = None,
                   inlet_temperature_K: float | None = None,
                   specific_heat_J_kgK: float | None = None,
                   temperature_of_enthalpy: Callable[[float], float] | None = None
                   ) -> Result:
    r"""Coolant state after absorbing ``Q``. ENTHALPY FIRST.

    .. math:: h_{out} = h_{in} + \dot Q / \dot m

    Enthalpy is the conserved quantity; temperature is not. For a coolant that
    crosses the pseudo-critical line, :math:`c_p` varies by an order of
    magnitude across the channel, so ``dT = Q/(m cp)`` with a single ``cp`` is
    not a small error -- it is the wrong equation. That form is still
    available, because sometimes it is all the caller has, but it returns
    `UNVALIDATED_ASSUMPTION` and says so.

    Returns the outlet enthalpy when enthalpy is supplied, or the outlet
    temperature when only ``cp`` is.
    """
    inputs = {"Q_W": heat_rate_W, "mdot_kg_s": mass_flow_kg_s}
    if mass_flow_kg_s <= 0:
        return Result(float("nan"), "J/kg", "CAN-THERM-ENERGY",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="coolant mass flow must be > 0")

    if inlet_enthalpy_J_kg is not None:
        h_out = inlet_enthalpy_J_kg + heat_rate_W / mass_flow_kg_s
        statuses = [Status.PASS]
        notes = []
        t_out = None
        if temperature_of_enthalpy is not None:
            t_out = temperature_of_enthalpy(h_out)
        else:
            statuses.append(Status.INSUFFICIENT_EVIDENCE)
            notes.append("no h -> T inversion supplied, so the outlet "
                         "TEMPERATURE is UNKNOWN. The enthalpy is exact; the "
                         "temperature is not derivable without a real-fluid "
                         "equation of state.")
        return Result(
            h_out, "J/kg", "CAN-THERM-ENERGY", status=worst(statuses),
            all_statuses=frozenset(statuses),
            evidence_level=EvidenceLevel.E1,
            verification=Verification.ALGEBRAIC,
            source="steady-flow energy equation, no work, no kinetic-energy "
                   "change",
            source_locator="definitional",
            validity_domain="steady state; all of Q enters the coolant",
            notes="; ".join(notes),
            inputs={**inputs, "h_in_J_kg": inlet_enthalpy_J_kg,
                    "T_out_K": t_out})

    if specific_heat_J_kgK is None or inlet_temperature_K is None:
        return Result(float("nan"), "J/kg", "CAN-THERM-ENERGY",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes="supply either (inlet enthalpy) or (inlet "
                            "temperature and cp). Neither was given, and "
                            "there is no default coolant.")
    if specific_heat_J_kgK <= 0 or inlet_temperature_K <= 0:
        return Result(float("nan"), "K", "CAN-THERM-ENERGY",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="cp and inlet temperature must be > 0")
    dt = heat_rate_W / (mass_flow_kg_s * specific_heat_J_kgK)
    return Result(
        inlet_temperature_K + dt, "K", "CAN-THERM-ENERGY",
        status=Status.UNVALIDATED_ASSUMPTION,
        all_statuses=frozenset({Status.UNVALIDATED_ASSUMPTION}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        source="steady-flow energy equation with a CONSTANT specific heat",
        source_locator="definitional",
        validity_domain="cp effectively constant over the temperature rise -- "
                        "FALSE for a coolant crossing its pseudo-critical line",
        uncertainty="unbounded near the pseudo-critical line, where cp for "
                    "methane varies by an order of magnitude",
        assumptions=(Assumption("specific_heat_J_kgK", specific_heat_J_kgK,
                                "J/(kg.K)",
                                "held constant over the whole temperature "
                                "rise; enthalpy was not available",
                                Status.UNVALIDATED_ASSUMPTION),),
        notes=f"dT = {dt:.1f} K. Prefer the enthalpy form.",
        inputs={**inputs, "T_in_K": inlet_temperature_K})


# ===========================================================================
# 7. THE WALL-TEMPERATURE SOLVE -- what replaces the hardcoded fractions
# ===========================================================================

@dataclass
class WallSolution:
    """A solved wall state, or an honest statement that it was not solved."""

    hot_wall_K: float
    cold_wall_K: float
    heat_flux_W_m2: float
    gas_film_coefficient_W_m2K: float
    coolant_film_coefficient_W_m2K: float
    wall_resistance_m2K_W: float
    mean_conductivity_W_mK: float
    iterations: int
    converged: bool
    result: Result
    trace: list[dict[str, Any]] = field(default_factory=list)


def solve_wall_temperature(
    *,
    recovery_temperature_K: float,
    coolant_bulk_temperature_K: float,
    gas_film_coefficient: Callable[[float], Result],
    coolant_film_coefficient_W_m2K: float,
    coolant_evidence: Result | None = None,
    wall_thickness_m: float,
    inner_radius_m: float | None = None,
    geometry: WallGeometry = WallGeometry.CYLINDRICAL,
    conductivity_property: str | None = None,
    constant_conductivity_W_mK: float | None = None,
    relaxation: float = 0.5,
    tolerance_K: float = 1e-6,
    max_iterations: int = 200,
) -> WallSolution:
    r"""Solve the three-resistance wall energy balance. THE replacement for
    ``T_wg = 0.22 T_0 + 100`` and ``T_wg = 0.40 T_0``.

    .. math::
        q = \frac{T_{aw} - T_{bulk}}{1/h_g + R_{wall} + 1/h_c},\qquad
        T_{wg} = T_{aw} - q/h_g,\qquad T_{wc} = T_{bulk} + q/h_c

    Two couplings make this a fixed point rather than one line of algebra, and
    both were missing from the previous implementation in
    ``material_manager/engine.py``:

    * :math:`h_g` depends on :math:`T_{wg}` through the Bartz property-variation
      correction :math:`\sigma`. That implementation computed :math:`\sigma`
      ONCE at the arbitrary starting guess ``0.40 T_c`` and never updated it
      inside the loop, so the converged answer still carried the guess.
    * :math:`k` depends on the wall temperatures. Here the INTEGRAL MEAN
      across the wall is used, which is exact for Fourier conduction rather
      than a mid-temperature approximation.

    `gas_film_coefficient` is a callable of the hot-wall temperature so that
    the caller keeps ownership of every Bartz argument; this function never
    constructs a Bartz call itself and therefore cannot smuggle in a default.

    Non-convergence is `FAIL`, not a silently returned last iterate.

    `coolant_evidence` is the Phase-14 gate added by the supercritical
    methane wave. Pass the `Result` from
    `canonical.coolant.coolant_side_h_c` alongside its float, and if that
    result carries `INSUFFICIENT_EVIDENCE`, `CONFLICT_UNRESOLVED` or
    `INVALID_MODEL_REGIME`, **this solver refuses to run**. The reason is
    specific: for supercritical methane in a rough rocket channel there is
    currently no correlation in the open literature simultaneously inside its
    fluid domain, its roughness domain and its Reynolds domain, and the
    previous behaviour -- silently applying the smooth-tube 0.023 form -- is
    wrong there by a factor of ~2.9 against the fully-rough branch. A wall
    temperature computed on an unsupported h_c is a precise number with no
    physics behind it.

    Omitting `coolant_evidence` preserves the old behaviour for callers that
    supply an h_c from elsewhere; they simply get no gate.
    """
    trace: list[dict[str, Any]] = []
    statuses: list[Status] = [Status.PASS]
    notes: list[str] = []

    if coolant_evidence is not None:
        _blocking = (Status.INSUFFICIENT_EVIDENCE, Status.CONFLICT_UNRESOLVED,
                     Status.INVALID_MODEL_REGIME, Status.PHYSICALLY_INVALID)
        _hit = [s for s in _blocking if coolant_evidence.has(s)]
        if _hit:
            res = Result(
                float("nan"), "K", "CAN-THERM-WALLSOLVE",
                status=Status.INSUFFICIENT_EVIDENCE,
                all_statuses=frozenset({Status.INSUFFICIENT_EVIDENCE,
                                        *coolant_evidence.all_statuses}),
                evidence_level=EvidenceLevel.E0,
                source=coolant_evidence.source,
                source_locator=coolant_evidence.source_locator,
                validity_domain=coolant_evidence.validity_domain,
                notes="THE WALL TEMPERATURE IS NOT SOLVED because the "
                      "coolant-side model is not supported at this state ("
                      + ", ".join(s.value for s in _hit) + "). No correlation "
                      "is substituted. Coolant-side reason: "
                      + (coolant_evidence.notes or "(none given)"),
                inputs={"T_aw_K": recovery_temperature_K,
                        "T_bulk_K": coolant_bulk_temperature_K,
                        "coolant_equation_id": coolant_evidence.equation_id,
                        "coolant_status": coolant_evidence.status.value})
            return WallSolution(float("nan"), float("nan"), float("nan"),
                                float("nan"), coolant_film_coefficient_W_m2K,
                                float("nan"), float("nan"), 0, False, res,
                                trace)
        statuses.extend(coolant_evidence.all_statuses)

    if recovery_temperature_K <= coolant_bulk_temperature_K:
        res = Result(float("nan"), "K", "CAN-THERM-WALLSOLVE",
                     status=Status.PHYSICALLY_INVALID,
                     notes="T_aw must exceed the coolant bulk temperature for "
                           "heat to flow into the coolant")
        return WallSolution(float("nan"), float("nan"), float("nan"),
                            float("nan"), coolant_film_coefficient_W_m2K,
                            float("nan"), float("nan"), 0, False, res, trace)
    if coolant_film_coefficient_W_m2K <= 0 or wall_thickness_m <= 0:
        res = Result(float("nan"), "K", "CAN-THERM-WALLSOLVE",
                     status=Status.PHYSICALLY_INVALID,
                     notes="h_c and wall thickness must be > 0")
        return WallSolution(float("nan"), float("nan"), float("nan"),
                            float("nan"), coolant_film_coefficient_W_m2K,
                            float("nan"), float("nan"), 0, False, res, trace)
    if conductivity_property is None and constant_conductivity_W_mK is None:
        res = Result(float("nan"), "K", "CAN-THERM-WALLSOLVE",
                     status=Status.INSUFFICIENT_EVIDENCE,
                     notes="supply either a temperature-dependent "
                           "conductivity_property (preferred) or an explicit "
                           "constant_conductivity_W_mK. There is no default "
                           "wall conductivity: the previous code used a bare "
                           "330 W/(m.K) for GRCop with no source and no "
                           "temperature dependence.")
        return WallSolution(float("nan"), float("nan"), float("nan"),
                            float("nan"), coolant_film_coefficient_W_m2K,
                            float("nan"), float("nan"), 0, False, res, trace)

    # Start from the midpoint of the physically admissible interval. This is a
    # numerical seed, not a model assumption: the fixed point is unique on
    # (T_bulk, T_aw) because q is monotone decreasing in T_wg, and the answer
    # is independent of the seed (asserted by test).
    t_wg = 0.5 * (recovery_temperature_K + coolant_bulk_temperature_K)
    t_wc = coolant_bulk_temperature_K
    q = float("nan")
    h_g = float("nan")
    r_wall = float("nan")
    k_bar = float("nan")
    converged = False
    it = 0

    for it in range(1, max_iterations + 1):
        hg_res = gas_film_coefficient(t_wg)
        if not hg_res.usable or float(hg_res) <= 0:
            res = Result(float("nan"), "K", "CAN-THERM-WALLSOLVE",
                         status=Status.PHYSICALLY_INVALID,
                         notes=f"gas film coefficient unusable at T_wg = "
                               f"{t_wg:.1f} K: {hg_res.notes}")
            return WallSolution(float("nan"), float("nan"), float("nan"),
                                float("nan"), coolant_film_coefficient_W_m2K,
                                float("nan"), float("nan"), it, False, res,
                                trace)
        h_g = float(hg_res)
        statuses.extend(hg_res.all_statuses)

        if conductivity_property is not None:
            # RESTRICT, not extrapolate: intermediate iterates wander outside
            # the property's measured range even when the converged answer
            # does not. The restriction is permanently flagged, so a wall that
            # genuinely sits outside the data still reports
            # OUT_OF_CORRELATION_RANGE in the final result.
            k_res = mean_wall_conductivity(
                material_property=conductivity_property,
                hot_side_K=t_wg, cold_side_K=t_wc,
                restrict_to_valid_range=True)
            if not k_res.usable or math.isnan(float(k_res)):
                res = Result(
                    float("nan"), "K", "CAN-THERM-WALLSOLVE",
                    status=k_res.status, all_statuses=k_res.all_statuses,
                    source=k_res.source, source_locator=k_res.source_locator,
                    notes=f"wall conductivity unavailable over "
                          f"[{t_wc:.1f}, {t_wg:.1f}] K: {k_res.notes}. NO "
                          f"EXTRAPOLATED VALUE IS SUBSTITUTED. The iterate "
                          f"reached at abort is NOT returned as a wall "
                          f"temperature.",
                    inputs={"aborted_iterate_T_wg_K": t_wg,
                            "aborted_iterate_T_wc_K": t_wc,
                            "iterations": it})
                return WallSolution(float("nan"), float("nan"), float("nan"),
                                    h_g, coolant_film_coefficient_W_m2K,
                                    float("nan"), float("nan"), it, False,
                                    res, trace)
            k_bar = float(k_res)
            statuses.extend(k_res.all_statuses)
        else:
            k_bar = float(constant_conductivity_W_mK)
            statuses.append(Status.UNVALIDATED_ASSUMPTION)

        rw = wall_resistance(thickness_m=wall_thickness_m,
                             conductivity_W_mK=k_bar,
                             inner_radius_m=inner_radius_m,
                             geometry=geometry)
        if not rw.usable or math.isnan(float(rw)):
            res = Result(float("nan"), "K", "CAN-THERM-WALLSOLVE",
                         status=rw.status, all_statuses=rw.all_statuses,
                         notes=rw.notes)
            return WallSolution(t_wg, t_wc, float("nan"), h_g,
                                coolant_film_coefficient_W_m2K, float("nan"),
                                k_bar, it, False, res, trace)
        r_wall = float(rw)
        statuses.extend(rw.all_statuses)

        r_total = 1.0 / h_g + r_wall + 1.0 / coolant_film_coefficient_W_m2K
        q = (recovery_temperature_K - coolant_bulk_temperature_K) / r_total
        t_wg_new = recovery_temperature_K - q / h_g
        t_wc_new = coolant_bulk_temperature_K + q / coolant_film_coefficient_W_m2K

        d = max(abs(t_wg_new - t_wg), abs(t_wc_new - t_wc))
        t_wg += relaxation * (t_wg_new - t_wg)
        t_wc += relaxation * (t_wc_new - t_wc)
        if d < tolerance_K:
            converged = True
            break

    if not converged:
        statuses.append(Status.FAIL)
        notes.append(f"fixed point did not converge in {max_iterations} "
                     f"iterations; the last iterate is NOT returned as an "
                     f"answer")

    if t_wc > t_wg:
        statuses.append(Status.PHYSICALLY_INVALID)
        notes.append("cold-side wall temperature exceeds the hot-side value")

    trace.append({"h_g_W_m2K": h_g,
                  "h_c_W_m2K": coolant_film_coefficient_W_m2K,
                  "R_wall_m2K_W": r_wall, "k_bar_W_mK": k_bar,
                  "iterations": it, "converged": converged})

    res = Result(
        t_wg if converged else float("nan"), "K", "CAN-THERM-WALLSOLVE",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.INTERNAL_CONSISTENCY,
        validation=Validation.NONE,
        source="steady one-dimensional energy balance across three "
               "resistances in series",
        source_locator="definitional; the EVIDENCE lives in the three "
                       "constituent results (h_g, k(T), h_c), not here",
        validity_domain="steady, one-dimensional, no axial conduction, no "
                        "channel-land fin effect, no radiation",
        uncertainty="dominated by h_g (Bartz reports NONE) and h_c "
                    "(+/-56% for the best supercritical correlation in the "
                    "open literature)",
        notes="; ".join(notes),
        inputs={"T_aw_K": recovery_temperature_K,
                "T_bulk_K": coolant_bulk_temperature_K,
                "h_g_W_m2K": h_g,
                "h_c_W_m2K": coolant_film_coefficient_W_m2K,
                "R_wall_m2K_W": r_wall, "k_bar_W_mK": k_bar,
                "q_W_m2": q, "T_wg_K": t_wg, "T_wc_K": t_wc,
                "iterations": it, "converged": converged,
                "geometry": geometry.value})

    if not converged:
        # A non-converged iterate is not a wall temperature. Returning it
        # would let a caller read a plausible number out of a failed solve,
        # which is how the guessed wall temperatures got into production in
        # the first place. The iterate survives in `res.inputs` for debugging.
        return WallSolution(float("nan"), float("nan"), float("nan"), h_g,
                            coolant_film_coefficient_W_m2K, r_wall, k_bar,
                            it, False, res, trace)
    return WallSolution(t_wg, t_wc, q, h_g, coolant_film_coefficient_W_m2K,
                        r_wall, k_bar, it, converged, res, trace)


# ===========================================================================
# 8. AXIAL MARCH AND HOT-SPOT DETECTION
# ===========================================================================

@dataclass
class ThermalStation:
    """One axial station. `area_ratio` is A_local/A_throat, >= 1 everywhere."""

    name: str
    axial_position_m: float
    area_ratio: float
    mach: float
    #: True only where the boundary layer is the one Bartz correlated.
    is_throat: bool = False


def axial_thermal_march(
    stations: Sequence[ThermalStation], *,
    station_flux: Callable[[ThermalStation], Result],
) -> list[dict[str, Any]]:
    """March the chain along the chamber, flagging EVERY station honestly.

    Bartz's correlation is anchored at the THROAT. Its ``(A_t/A)^0.9`` term is
    an extension away from the throat, not a validated one -- the primary
    states the domain as the throat boundary layer. So every station with
    ``area_ratio > 1`` is marked `OUT_OF_CORRELATION_RANGE` here, and that
    marking survives into the trace rather than being averaged away by a
    chamber-mean heat load.

    This function deliberately does NOT integrate the flux into a total heat
    rate. Integrating an out-of-domain flux profile produces a single
    confident-looking number whose provenance is gone -- which is exactly what
    ``Q_total = q_throat * A_wetted * 0.35`` did.
    """
    out: list[dict[str, Any]] = []
    for st in stations:
        r = station_flux(st)
        statuses = set(r.all_statuses)
        notes = [r.notes] if r.notes else []
        if st.area_ratio > 1.0 + 1e-12 and not st.is_throat:
            statuses.add(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(
                f"station {st.name!r} is at A/A_t = {st.area_ratio:.3f}; the "
                f"Bartz primary states the throat boundary layer as its "
                f"domain, so this is an extension of the correlation, not an "
                f"application of it")
        out.append({
            "station": st.name,
            "x_m": st.axial_position_m,
            "area_ratio": st.area_ratio,
            "mach": st.mach,
            "is_throat": st.is_throat,
            "value": r.value,
            "units": r.units,
            "status": worst(statuses).value,
            "all_statuses": sorted(s.value for s in statuses),
            "equation_id": r.equation_id,
            "source": r.source,
            "notes": "; ".join(n for n in notes if n),
        })
    return out


def hot_spot(march: Sequence[dict[str, Any]]) -> Result:
    """Locate the maximum heat flux. DOES NOT ASSUME IT IS THE THROAT.

    The throat maximises ``h_g`` but not necessarily ``q = h_g (T_aw - T_wg)``:
    the recovery temperature falls with Mach number, the wall temperature
    varies with the local cooling, and the coolant heats up along its own path.
    NASA SP-8087 -- which was opened -- records the opposite experience
    directly: *"The weakest link in the analysis is the analytical description
    of the gas-side thermal conditions, especially in the region just below the
    injector"*, with *"continuing problems of tube burnouts just below the
    injector"*. The burnouts were at the injector end, not the throat.

    So the hot spot is found by search, and if it is not at the throat the
    result says so.
    """
    usable = [m for m in march
              if isinstance(m.get("value"), (int, float))
              and not math.isnan(float(m["value"]))]
    if not usable:
        return Result(float("nan"), "W/m^2", "CAN-THERM-HOTSPOT",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="no station produced a usable flux")
    peak = max(usable, key=lambda m: float(m["value"]))
    statuses = {Status[s] for s in peak["all_statuses"]}
    notes = [f"peak at station {peak['station']!r}, x = {peak['x_m']:.4f} m, "
             f"A/A_t = {peak['area_ratio']:.3f}"]
    if not peak["is_throat"]:
        statuses.add(Status.OUT_OF_CORRELATION_RANGE)
        notes.append("THE HOT SPOT IS NOT AT THE THROAT. The correlation is "
                     "anchored at the throat, so the peak is being read from "
                     "an extrapolated region and the location is less "
                     "trustworthy than the fact that it is not where a "
                     "throat-only screen would have looked.")
    throats = [m for m in usable if m["is_throat"]]
    if throats:
        q_t = float(throats[0]["value"])
        if q_t > 0:
            notes.append(f"peak/throat = {float(peak['value']) / q_t:.3f}")
    return Result(
        float(peak["value"]), peak.get("units", "W/m^2"),
        "CAN-THERM-HOTSPOT", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="search over the marched stations; no assumption that the "
               "throat is the hot spot",
        source_locator="see NASA SP-8087 p.15 on injector-end tube burnouts",
        validity_domain="only as good as the station set supplied",
        notes="; ".join(notes),
        inputs={"n_stations": len(march), "n_usable": len(usable),
                "peak_station": peak["station"]})
