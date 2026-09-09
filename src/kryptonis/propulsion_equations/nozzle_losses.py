"""
Nozzle losses and separation. The three published criteria DISAGREE, and this
module reports the disagreement instead of hiding it behind a choice.
==============================================================================

WHAT THIS MODULE IS FOR
-----------------------
The physics-closure audit established four things about the nozzle, by
measurement rather than by reading:

  1. The live thrust coefficient is **IDEAL and VACUUM**. `_thrust_coefficient`
     takes `p_ambient_over_pc` with a default of 0.0 and the sole production
     call site passed two arguments, so the ambient term never fired.

  2. **No divergence loss is applied anywhere on a live path**, although the
     arithmetic has existed in `kryptonis/nozzle/contour.py` since before the
     audit with no caller at all. For a 15-degree cone that is **1.70 % of
     C_F** -- larger than the entire equilibrium-kernel disagreement seen
     through gamma (~0.13 % on throat area).

  3. **No separation criterion runs.** At the 30 kN reference case (eps = 40,
     Pc = 10 MPa) the exit pressure is 20.99 kPa against a sea-level ambient
     of 101.3 kPa: `P_e/P_a = 0.207`, roughly HALF the Summerfield attachment
     limit.

  4. **The three published separation criteria contradict each other on this
     nozzle.** Summerfield and Kalt-Badal say separated; Schmucker says
     attached. Their predicted separation pressures span a factor of ~21.

WHY NO CRITERION IS SELECTED
----------------------------
Point 4 is the whole reason this module returns a SPREAD rather than a
boolean. Implementing one criterion and calling the question closed would be
exactly the false closure the audit programme exists to prevent. When the
criteria disagree the status is `CONFLICT_UNRESOLVED` and the caller is handed
every verdict, every source and the spread between them.

**`separation_assessment` will never tell you a nozzle is safe.** It tells you
what three published authors would each have said, and whether they agree.

SOURCE ACCESS, STATED HONESTLY
------------------------------
All three criteria were obtained from ONE secondary source that this project
DID open:

    Stark, R., "Flow Separation in Rocket Nozzles, a Simple Criteria",
    AIAA 2005-3940 (DLR elib 49253).

The three originals -- Summerfield 1954, Schmucker 1973, Kalt & Badal 1965 --
were **NOT opened**. Every criterion here therefore carries `E1` and
`access_status = NOT_OPENED` on its primary, and no amount of agreement
between them may promote that. Stark also states that criteria "like
Schmucker's under predict the separation location", which is consistent with
Schmucker being the outlier in the reference case -- and that consistency is
an observation, not a validation.

ONE ASSUMPTION, MADE EXPLICIT
-----------------------------
Every criterion is written on the WALL STATIC PRESSURE AT SEPARATION,
``p_sep``. Applying it to the EXIT pressure ``P_e`` tests only whether
separation would occur AT THE EXIT PLANE. That is the standard screening use
and it is an **assumption**, not an identity: a nozzle can separate upstream of
its exit while its (hypothetical, attached) exit pressure still passes. Every
`Result` from this module carries that assumption explicitly.

WHAT IS STILL ABSENT
--------------------
**Kinetic (finite-rate recombination) loss is NOT implemented and is not
estimated here.** It requires a finite-rate expansion solver that does not
exist in this repository. It is recorded as absent so that a caller summing
these losses knows the sum is incomplete.
"""

from __future__ import annotations

import math

from kryptonis.propulsion_equations.units import (
    Assumption, EvidenceLevel, NOT_REPORTED, Result, Status, Validation,
    Verification, worst,
)
from kryptonis.propulsion_equations.contour import (
    calculate_bell_divergence_loss_factor, calculate_displacement_thickness,
    calculate_divergence_loss_factor, calculate_effective_area,
)

__all__ = [
    "STARK_2005", "SEPARATION_CRITERIA", "SEA_LEVEL_PA",
    "divergence_efficiency", "separation_assessment",
    "boundary_layer_displacement", "delivered_thrust_coefficient",
    "KINETIC_LOSS_ABSENT",
]

#: Standard sea-level ambient pressure [Pa]. Used only where a caller asks for
#: it BY NAME; nothing in this module defaults to it.
SEA_LEVEL_PA = 101_325.0

#: The one document this project opened for any of the separation criteria.
STARK_2005 = (
    "Stark, R., 'Flow Separation in Rocket Nozzles, a Simple Criteria', "
    "AIAA 2005-3940 (DLR elib 49253). SECONDARY SOURCE, opened. The three "
    "criteria below are quoted from it; their originals were NOT opened.")

#: A note carried on every result that sums losses, so that nobody reads the
#: sum as complete.
KINETIC_LOSS_ABSENT = (
    "KINETIC_LOSS_NOT_IMPLEMENTED: finite-rate recombination loss during "
    "expansion is absent from this repository and is NOT included in any "
    "efficiency returned here. The delivered coefficient is therefore an "
    "UPPER BOUND with respect to that mechanism.")


def _summerfield(_mach_exit: float, _pc_over_pa: float, bound: float) -> float:
    """p_sep/p_a, constant. `bound` selects the 0.35 or 0.40 end."""
    return bound


def _schmucker(mach_exit: float, _pc_over_pa: float, _b: float) -> float:
    """p_sep/p_a = 0.64 Ma_sep^-1.88."""
    return 0.64 * mach_exit ** -1.88


def _kalt_badal(_mach_exit: float, pc_over_pa: float, _b: float) -> float:
    """p_sep/p_a = (1/5) (p_c/p_a)^(1/3)."""
    return 0.2 * pc_over_pa ** (1.0 / 3.0)


#: key -> (callable, bound-arg, author, locator, note)
SEPARATION_CRITERIA: dict[str, tuple] = {
    "summerfield_1954_upper": (
        _summerfield, 0.40,
        "Summerfield M., Foster C., Swan W., 'Flow Separation in "
        "Overexpanded Supersonic Exhaust Nozzles'",
        "Jet Propulsion 24(9), 319-321, 1954",
        "the upper end of the 0.35-0.40 band Stark quotes"),
    "summerfield_1954_lower": (
        _summerfield, 0.35,
        "Summerfield M., Foster C., Swan W., 'Flow Separation in "
        "Overexpanded Supersonic Exhaust Nozzles'",
        "Jet Propulsion 24(9), 319-321, 1954",
        "the lower end of the same band"),
    "schmucker_1973": (
        _schmucker, 0.0,
        "Schmucker R., 'Stroemungsvorgaenge beim Betrieb ueberexpandierter "
        "Duesen chemischer Raketentriebwerke, Teil 1: Stroemungsabloesung'",
        "Bericht TB-7, TU Munich, 1973",
        "Mach-dependent. Stark states that criteria like this one UNDER "
        "PREDICT the separation location."),
    "kalt_badal_1965": (
        _kalt_badal, 0.0,
        "Kalt S., Badal D., 'Conical Rocket Nozzle Performance under "
        "Flow-Separated Conditions'",
        "J. Spacecraft and Rockets 2(3), 447-449, 1965",
        "pressure-ratio dependent; the most conservative of the three here"),
}

_EXIT_PLANE_ASSUMPTION = Assumption(
    "separation_tested_at_the_exit_plane", 1.0, "-",
    "every published criterion is written on the WALL STATIC PRESSURE AT "
    "SEPARATION p_sep. This module applies it to the EXIT pressure P_e, "
    "which tests only whether separation occurs AT the exit plane. A nozzle "
    "can separate upstream while its attached-flow exit pressure still "
    "passes. Standard screening practice, and an assumption.",
    Status.UNVALIDATED_ASSUMPTION)


# ==========================================================================
# DIVERGENCE
# ==========================================================================

def divergence_efficiency(
    *,
    contour: str,
    exit_half_angle_deg: float | None = None,
    initial_wall_angle_deg: float | None = None,
) -> Result:
    r"""Divergence loss factor :math:`\lambda`, conical or bell.

    .. math::
        \lambda_{cone} = \frac{1 + \cos\alpha}{2}
        \qquad
        \lambda_{bell} \approx
            \frac{1 + \cos\left(\frac{\theta_n+\theta_e}{2}\right)}{2}

    The arithmetic lives in :mod:`kryptonis.nozzle.contour` and is CALLED from
    here rather than repeated -- there is one implementation of each formula in
    the repository and this is the canonical entry point to it.

    `contour` is required and has no default. Applying the conical factor to a
    bell contour underpredicts delivered Isp and biases an optimiser toward
    cones, which is a documented failure mode of the underlying module, so the
    caller must say which shape it has.

    Parameters
    ----------
    contour
        ``"conical"`` or ``"bell"``.
    exit_half_angle_deg
        Conical: the half-angle alpha. Bell: the exit wall angle theta_e.
    initial_wall_angle_deg
        Bell only: the initial parabola angle theta_n.
    """
    shape = (contour or "").strip().lower()
    if shape not in ("conical", "bell"):
        return Result(
            value=None, units="-", equation_id="nozzle_divergence_efficiency",
            status=Status.INSUFFICIENT_EVIDENCE,
            source="n/a", validity_domain="VALIDITY_DOMAIN_UNKNOWN",
            notes=(f"contour must be 'conical' or 'bell', not {contour!r}. "
                   f"There is no shape-agnostic divergence factor: the "
                   f"conical form assumes uniform spherical source flow at "
                   f"the exit plane, which a bell contour does not have."))

    if exit_half_angle_deg is None:
        return Result(
            value=None, units="-", equation_id="nozzle_divergence_efficiency",
            status=Status.INSUFFICIENT_EVIDENCE, source="n/a",
            notes="exit_half_angle_deg is required; there is no default "
                  "cone angle that is safe to assume.")

    if shape == "conical":
        lam = calculate_divergence_loss_factor(exit_half_angle_deg)
        return Result(
            value=lam, units="-", equation_id="nozzle_divergence_conical",
            status=Status.PASS,
            evidence_level=EvidenceLevel.E2,
            verification=Verification.ALGEBRAIC,
            validation=Validation.NONE,
            source="exact integral of uniform spherical source flow over a "
                   "conical exit plane; first principles",
            source_locator="derived; see kryptonis/nozzle/contour.py",
            validity_domain="CONICAL contours only. Assumes uniform "
                            "spherical source flow at the exit plane.",
            uncertainty="EXACT for the source-flow model; the MODEL FORM "
                        "itself is an idealisation and its error is "
                        "NOT_REPORTED",
            assumptions=(
                Assumption("uniform_spherical_source_flow", 1.0, "-",
                           "the exit velocity vector is radial from a point "
                           "source; real conical nozzles depart from this "
                           "near the wall", Status.UNVALIDATED_ASSUMPTION),),
            inputs={"exit_half_angle_deg": exit_half_angle_deg},
            notes=f"lambda = (1 + cos {exit_half_angle_deg:.3f} deg)/2 = "
                  f"{lam:.6f}; a thrust loss of {100 * (1 - lam):.3f} %.")

    if initial_wall_angle_deg is None:
        return Result(
            value=None, units="-", equation_id="nozzle_divergence_bell",
            status=Status.INSUFFICIENT_EVIDENCE, source="n/a",
            notes="a bell contour needs BOTH theta_n (initial wall angle) "
                  "and theta_e (exit angle). Only theta_e was supplied.")
    try:
        lam = calculate_bell_divergence_loss_factor(
            initial_wall_angle_deg, exit_half_angle_deg)
    except ValueError as exc:
        return Result(
            value=None, units="-", equation_id="nozzle_divergence_bell",
            status=Status.PHYSICALLY_INVALID, source="n/a", notes=str(exc))
    return Result(
        value=lam, units="-", equation_id="nozzle_divergence_bell",
        status=Status.UNVALIDATED_ASSUMPTION,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source="Rao, G.V.R., 'Exhaust Nozzle Contour for Optimum Thrust', "
               "Jet Propulsion 28, 1958; the MEAN-ANGLE approximation is the "
               "standard engineering shortcut, not Rao's own result",
        source_locator="NOT_OPENED",
        validity_domain="parabolic (Rao-style) bell contours. The mean-angle "
                        "form is an APPROXIMATION; the exact value needs a "
                        "method-of-characteristics exit-plane integration, "
                        "which this repository does not have.",
        uncertainty="NOT_REPORTED. The mean-angle shortcut has no published "
                    "error bound in any source this project has opened.",
        assumptions=(
            Assumption("mean_angle_approximation", 1.0, "-",
                       "the conical formula is evaluated at (theta_n+theta_e)"
                       "/2. This is a shortcut, not a derivation, and it is "
                       "why the status is UNVALIDATED_ASSUMPTION rather than "
                       "PASS.", Status.UNVALIDATED_ASSUMPTION),),
        inputs={"initial_wall_angle_deg": initial_wall_angle_deg,
                "exit_half_angle_deg": exit_half_angle_deg},
        notes=f"lambda_bell = {lam:.6f}; a thrust loss of "
              f"{100 * (1 - lam):.3f} %. Method-of-characteristics NOT "
              f"implemented.")


# ==========================================================================
# SEPARATION
# ==========================================================================

def separation_assessment(
    *,
    exit_pressure_Pa: float,
    ambient_pressure_Pa: float,
    exit_mach: float,
    chamber_pressure_Pa: float,
) -> Result:
    """Every published criterion, applied to the same nozzle. NO winner.

    Returns a `Result` whose value is a dict of per-criterion verdicts. The
    STATUS is the finding:

    * ``PASS`` -- every criterion says the flow stays attached.
    * ``FAIL`` -- every criterion says it separates.
    * ``CONFLICT_UNRESOLVED`` -- **they disagree, and this module will not
      break the tie.** Choosing among criteria that contradict each other
      requires evidence this repository does not have.

    `ambient_pressure_Pa` is required and has no default. A separation
    assessment without an ambient pressure is not a conservative assessment,
    it is a meaningless one.
    """
    if ambient_pressure_Pa is None or ambient_pressure_Pa <= 0.0:
        return Result(
            value=None, units="-", equation_id="nozzle_separation",
            status=Status.INSUFFICIENT_EVIDENCE,
            source=STARK_2005,
            notes="ambient pressure is required and must be positive. In "
                  "vacuum there is no separation criterion to apply, and the "
                  "caller should say so explicitly rather than pass 0.")
    if exit_mach <= 1.0:
        return Result(
            value=None, units="-", equation_id="nozzle_separation",
            status=Status.INVALID_MODEL_REGIME, source=STARK_2005,
            notes=f"exit Mach {exit_mach:.4f} is not supersonic; these "
                  f"criteria describe shock-induced separation in an "
                  f"OVEREXPANDED SUPERSONIC nozzle and have no referent here.")

    pc_over_pa = chamber_pressure_Pa / ambient_pressure_Pa
    verdicts: dict[str, dict] = {}
    attached: list[bool] = []
    ratios: list[float] = []
    for key, (fn, bound, author, locator, note) in SEPARATION_CRITERIA.items():
        ratio = fn(exit_mach, pc_over_pa, bound)
        p_sep = ratio * ambient_pressure_Pa
        is_attached = exit_pressure_Pa >= p_sep
        attached.append(is_attached)
        ratios.append(ratio)
        verdicts[key] = {
            "p_sep_over_p_a": ratio,
            "p_sep_Pa": p_sep,
            "verdict": "ATTACHED" if is_attached else "SEPARATED",
            "margin_ratio": exit_pressure_Pa / p_sep if p_sep > 0 else None,
            "author": author,
            "source_locator": locator,
            "access_status": "NOT_OPENED -- quoted from Stark 2005",
            "note": note,
        }

    spread = max(ratios) / min(ratios) if min(ratios) > 0 else float("inf")
    if all(attached):
        status = Status.PASS
        headline = ("every published criterion says the flow stays attached "
                    "at the exit plane")
    elif not any(attached):
        status = Status.FAIL
        headline = ("every published criterion says this nozzle separates at "
                    "the stated ambient pressure")
    else:
        status = Status.CONFLICT_UNRESOLVED
        agree_sep = [k for k, v in verdicts.items()
                     if v["verdict"] == "SEPARATED"]
        agree_att = [k for k, v in verdicts.items()
                     if v["verdict"] == "ATTACHED"]
        headline = (
            f"THE PUBLISHED CRITERIA DISAGREE. Separated: "
            f"{', '.join(sorted(agree_sep))}. Attached: "
            f"{', '.join(sorted(agree_att))}. Predicted separation pressures "
            f"span a factor of {spread:.1f}. This module does NOT break the "
            f"tie: choosing among contradictory criteria needs evidence that "
            f"is not in this repository.")

    return Result(
        value=verdicts, units="-", equation_id="nozzle_separation",
        status=status,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.NONE,
        validation=Validation.NONE,
        source=STARK_2005,
        source_locator="AIAA 2005-3940; the three originals NOT_OPENED",
        validity_domain="overexpanded supersonic nozzles. Each criterion's "
                        "own domain is contested -- Stark states that "
                        "criteria like Schmucker's under predict the "
                        "separation location.",
        uncertainty=f"the criteria span a factor of {spread:.1f} in "
                    f"p_sep/p_a on this nozzle. That spread IS the "
                    f"uncertainty and no narrower figure is available.",
        assumptions=(_EXIT_PLANE_ASSUMPTION,),
        inputs={"exit_pressure_Pa": exit_pressure_Pa,
                "ambient_pressure_Pa": ambient_pressure_Pa,
                "exit_mach": exit_mach,
                "chamber_pressure_Pa": chamber_pressure_Pa},
        notes=headline)


# ==========================================================================
# BOUNDARY LAYER
# ==========================================================================

def boundary_layer_displacement(
    *,
    running_length_m: float,
    velocity_m_s: float,
    density_kg_m3: float,
    viscosity_Pa_s: float,
    geometric_area_m2: float,
    wall_radius_m: float,
) -> Result:
    r"""Displacement thickness and the effective area it leaves.

    .. math::
        \delta = \frac{0.37\,x}{Re_x^{0.2}}, \qquad
        \delta^* = \delta/8, \qquad
        A_{eff} = A_{geo} - 2\pi R \delta^*

    Connects `kryptonis/nozzle/contour.py`, which has carried this arithmetic
    with **no caller anywhere in the repository** since before the audit.

    The 1/7th-power turbulent flat-plate form is INCOMPRESSIBLE. A nozzle
    boundary layer is strongly compressible and strongly cooled, both of which
    change the profile. The status is therefore `UNVALIDATED_ASSUMPTION`, not
    `PASS`, however cleanly the arithmetic runs.
    """
    try:
        d_star = calculate_displacement_thickness(
            running_length_m, velocity_m_s, density_kg_m3, viscosity_Pa_s)
        a_eff = calculate_effective_area(
            geometric_area_m2, wall_radius_m, d_star)
    except ValueError as exc:
        return Result(
            value=None, units="m2",
            equation_id="nozzle_boundary_layer_effective_area",
            status=Status.OUT_OF_CORRELATION_RANGE,
            source="Schlichting, 'Boundary-Layer Theory', 7th ed.",
            source_locator="NOT_OPENED",
            notes=f"{exc} The displacement thickness has closed the annulus, "
                  f"which means the flat-plate form is being applied far "
                  f"outside any regime it describes.")

    re_x = density_kg_m3 * velocity_m_s * running_length_m / viscosity_Pa_s
    return Result(
        value={"displacement_thickness_m": d_star,
               "effective_area_m2": a_eff,
               "area_deficit_fraction": 1.0 - a_eff / geometric_area_m2,
               "reynolds_x": re_x},
        units="m2", equation_id="nozzle_boundary_layer_effective_area",
        status=Status.UNVALIDATED_ASSUMPTION,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.DIMENSIONAL,
        validation=Validation.NONE,
        source="Schlichting, 'Boundary-Layer Theory', 7th ed.; 1/7th-power "
               "turbulent flat-plate boundary layer",
        source_locator="NOT_OPENED",
        validity_domain="INCOMPRESSIBLE turbulent flat plate. A nozzle "
                        "boundary layer is compressible and strongly cooled; "
                        "neither effect is in this form. Re_x should exceed "
                        "~5e5 for the turbulent form to apply at all "
                        f"(here Re_x = {re_x:.3e}).",
        uncertainty="NOT_REPORTED. A compressible integral method or CFD "
                    "supersedes this; no error bound against either exists "
                    "in this repository.",
        assumptions=(
            Assumption("incompressible_flat_plate", 1.0, "-",
                       "the 0.37 x Re^-0.2 form and the delta/8 displacement "
                       "ratio are both incompressible flat-plate results "
                       "applied to a compressible, cooled, curved wall",
                       Status.UNVALIDATED_ASSUMPTION),),
        inputs={"running_length_m": running_length_m,
                "velocity_m_s": velocity_m_s,
                "density_kg_m3": density_kg_m3,
                "viscosity_Pa_s": viscosity_Pa_s,
                "geometric_area_m2": geometric_area_m2,
                "wall_radius_m": wall_radius_m},
        notes=f"delta* = {d_star * 1e6:.2f} um; the effective area is "
              f"{100 * (1 - a_eff / geometric_area_m2):.4f} % below the "
              f"geometric area.")


# ==========================================================================
# THE ASSEMBLED COEFFICIENT
# ==========================================================================

def delivered_thrust_coefficient(
    *,
    ideal_thrust_coefficient: float,
    divergence: Result | None = None,
    separation: Result | None = None,
    boundary_layer: Result | None = None,
) -> Result:
    """Assemble the ideal coefficient with whichever losses were supplied.

    **This is deliberately additive-by-omission.** A loss that is not supplied
    is NOT silently assumed to be unity in the notes -- it is named as absent,
    so that a delivered coefficient can never be mistaken for a complete one.

    Kinetic loss is absent from the repository entirely and is named as such
    on every result.
    """
    cf = float(ideal_thrust_coefficient)
    applied: list[str] = []
    absent: list[str] = []
    statuses = [Status.PASS]

    if divergence is not None and divergence.value is not None:
        cf *= float(divergence.value)
        applied.append(f"divergence lambda = {float(divergence.value):.6f}")
        statuses.append(divergence.status)
    else:
        absent.append("DIVERGENCE (no contour supplied)")

    if boundary_layer is not None and boundary_layer.value is not None:
        deficit = float(boundary_layer.value["area_deficit_fraction"])
        cf *= (1.0 - deficit)
        applied.append(f"boundary-layer area deficit = {100 * deficit:.4f} %")
        statuses.append(boundary_layer.status)
    else:
        absent.append("BOUNDARY LAYER (not evaluated)")

    if separation is not None:
        statuses.append(separation.status)
        if separation.status is Status.FAIL:
            absent.append(
                "SEPARATED FLOW -- every criterion says this nozzle "
                "separates, and NO separated-flow thrust model exists here. "
                "The coefficient below describes attached flow that will not "
                "occur.")
        elif separation.status is Status.CONFLICT_UNRESOLVED:
            absent.append(
                "SEPARATION UNRESOLVED -- the published criteria disagree "
                "about this nozzle, so whether the coefficient below applies "
                "at all is undecided.")
    else:
        absent.append("SEPARATION (not evaluated -- no ambient pressure)")

    absent.append("KINETIC/finite-rate recombination loss (ABSENT FROM THE "
                  "REPOSITORY)")

    return Result(
        value=cf, units="-", equation_id="nozzle_delivered_thrust_coefficient",
        status=worst(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source="assembled from the component losses; the assembly itself is "
               "a product of independent efficiency factors, which is the "
               "standard engineering treatment and NOT a derived result",
        source_locator="derived",
        validity_domain="only as wide as the NARROWEST component supplied; "
                        "see each component Result.",
        uncertainty="NOT_REPORTED for the assembly. The components' "
                    "uncertainties are not independent and were not combined.",
        assumptions=(
            Assumption("losses_multiply_independently", 1.0, "-",
                       "efficiency factors are applied as a product. Real "
                       "loss mechanisms interact -- a thickened boundary "
                       "layer changes the effective divergence angle -- and "
                       "that coupling is not modelled.",
                       Status.UNVALIDATED_ASSUMPTION),),
        inputs={"ideal_thrust_coefficient": ideal_thrust_coefficient},
        notes=("APPLIED: " + ("; ".join(applied) if applied else "nothing")
               + ". NOT APPLIED: " + "; ".join(absent) + ". "
               + KINETIC_LOSS_ABSENT))
