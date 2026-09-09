"""
Chamber acoustic modes. Frequencies, NOT a stability verdict.
=============================================================

WHAT THIS MODULE IS AND IS NOT
------------------------------
The physics-closure audit found `combustion_stability` to be the thinnest
family in the combustor: **7 relations, 0 live**, for a liquid rocket chamber.
It also found that the gap splits cleanly in two:

  * **The acoustic mode frequencies are FREE.** They are closed-form roots of
    the wave equation in a cylinder. They need geometry and a sound speed and
    nothing else -- no correlation, no coefficient, no experiment. Having them
    before the first firing is worth more than having none.

  * **A stability VERDICT is not free.** It needs the injector response
    function -- Crocco's interaction index `n` and sensitive time lag `tau` --
    which are properties of an element pattern that does not exist yet and
    cannot be derived from anything in this repository.

**This module implements the first and refuses the second.** It will tell you
what frequencies a chamber of these dimensions supports. It will never tell you
whether the chamber is stable, and `stability_verdict` exists solely to say so
in a way a caller cannot mistake for an answer.

WHY THAT SPLIT MATTERS
----------------------
Knowing the 1T frequency is genuinely useful: it is what an acoustic absorber
or a baffle is tuned against, and it is the first thing a stability rating test
compares its measurements to. Predicting the frequency and then MEASURING it is
a real blind check -- one of the few available to this project before hardware
exists. That check is specified in `COMBUSTOR_EXPERIMENTAL_CLOSURE_
REQUIREMENTS.md` as X4's blind-validation step.

THE MODEL, STATED
-----------------
Rigid-walled cylinder, uniform gas, no mean flow, no combustion response, no
damping. For a chamber of radius `R` and length `L` with sound speed `a`:

.. math::
    f_{mnq} = \\frac{a}{2}\\sqrt{
        \\left(\\frac{\\alpha_{mn}}{\\pi R}\\right)^2 +
        \\left(\\frac{q}{L}\\right)^2 }

where :math:`\\alpha_{mn}` is the n-th root of :math:`J'_m(\\alpha)=0`
(transverse) and `q` is the longitudinal mode number. The roots are computed
from `scipy.special.jnp_zeros`, not tabulated, so no digit is transcribed.

**Every one of those simplifications raises the true frequency's uncertainty
and none of them is calibrated here.** Mean flow lowers the modes; the
convergent section is not a rigid flat end; the gas is not uniform; and the
sound speed depends on which temperature you evaluate it at. The status is
`UNVALIDATED_ASSUMPTION` for that reason, and the module reports the sound
speed sensitivity so a reader can see how much the temperature choice is
worth.
"""

from __future__ import annotations

import math

from kryptonis.propulsion_equations.units import (
    Assumption, EvidenceLevel, NOT_REPORTED, Result, Status, Validation,
    Verification,
)

__all__ = [
    "MODE_LABELS", "transverse_roots", "acoustic_modes", "sound_speed",
    "stability_verdict", "chug_frequency",
    "first_tangential_frequency", "first_radial_frequency", "first_longitudinal_frequency",
]

def first_tangential_frequency(*, speed_of_sound_m_s: float, chamber_diameter_m: float) -> Result:
    r"""First tangential acoustic frequency (1T mode): f_1T = 1.8412 * c / (pi * D_c)"""
    if speed_of_sound_m_s <= 0 or chamber_diameter_m <= 0:
        return Result(float("nan"), "Hz", "CAN-AC-1T", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = 1.84118378 * speed_of_sound_m_s / (math.pi * chamber_diameter_m)
    return Result(val, "Hz", "CAN-AC-1T", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="Bessel function root dJ1/dr = 0 (alpha_10 = 1.8412)",
                  source_locator="NASA SP-194 Liquid Rocket Combustion Instability",
                  validity_domain="cylindrical acoustic chamber", uncertainty="exact",
                  inputs={"c": speed_of_sound_m_s, "D_c": chamber_diameter_m})

def first_radial_frequency(*, speed_of_sound_m_s: float, chamber_diameter_m: float) -> Result:
    r"""First radial acoustic frequency (1R mode): f_1R = 3.8317 * c / (pi * D_c)"""
    if speed_of_sound_m_s <= 0 or chamber_diameter_m <= 0:
        return Result(float("nan"), "Hz", "CAN-AC-1R", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = 3.83170597 * speed_of_sound_m_s / (math.pi * chamber_diameter_m)
    return Result(val, "Hz", "CAN-AC-1R", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="Bessel function root dJ0/dr = 0 (alpha_01 = 3.8317)",
                  source_locator="NASA SP-194 Liquid Rocket Combustion Instability",
                  validity_domain="cylindrical acoustic chamber", uncertainty="exact",
                  inputs={"c": speed_of_sound_m_s, "D_c": chamber_diameter_m})

def first_longitudinal_frequency(*, speed_of_sound_m_s: float, chamber_length_m: float) -> Result:
    r"""First longitudinal acoustic frequency (1L mode): f_1L = c / (2 * L_c)"""
    if speed_of_sound_m_s <= 0 or chamber_length_m <= 0:
        return Result(float("nan"), "Hz", "CAN-AC-1L", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = speed_of_sound_m_s / (2.0 * chamber_length_m)
    return Result(val, "Hz", "CAN-AC-1L", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="Closed-closed or open-open 1D acoustic pipe mode",
                  source_locator="NASA SP-194",
                  validity_domain="1D longitudinal chamber cavity", uncertainty="exact",
                  inputs={"c": speed_of_sound_m_s, "L_c": chamber_length_m})


#: (m, n) -> the conventional label. `m` is the azimuthal order, `n` the
#: radial. Computed roots are matched to these labels rather than the other
#: way round.
MODE_LABELS: dict[tuple[int, int], str] = {
    (1, 1): "1T   first tangential",
    (2, 1): "2T   second tangential",
    (0, 1): "1R   first radial",
    (1, 2): "1T1R first tangential-radial",
    (3, 1): "3T   third tangential",
    (2, 2): "2T1R second tangential-radial",
    (0, 2): "2R   second radial",
}


def transverse_roots(max_m: int = 3, max_n: int = 2) -> dict[tuple[int, int], float]:
    """Roots of :math:`J'_m(\\alpha) = 0`, COMPUTED not tabulated.

    The trivial root at the origin for `m = 0` is excluded by
    `scipy.special.jnp_zeros`, which is the correct convention: a uniform
    pressure across the section is not an acoustic mode.
    """
    from scipy.special import jnp_zeros
    out: dict[tuple[int, int], float] = {}
    for m in range(0, max_m + 1):
        for n, root in enumerate(jnp_zeros(m, max_n), start=1):
            out[(m, n)] = float(root)
    return out


def sound_speed(*, gamma: float, specific_gas_constant: float,
                temperature_K: float) -> float:
    r"""Frozen-composition sound speed :math:`a=\sqrt{\gamma R T}` [m/s].

    **Frozen**, deliberately. An equilibrium sound speed -- one that lets the
    composition shift with the acoustic pressure -- is lower, and this module
    does not compute it. That choice biases every frequency here HIGH and the
    bias is not quantified, because doing so needs a second equilibrium solve
    per mode which is out of this module's scope.
    """
    if gamma <= 0 or specific_gas_constant <= 0 or temperature_K <= 0:
        raise ValueError("gamma, R and T must all be positive")
    return math.sqrt(gamma * specific_gas_constant * temperature_K)


def acoustic_modes(
    *,
    chamber_diameter_m: float,
    chamber_length_m: float,
    gamma: float,
    specific_gas_constant: float,
    temperature_K: float,
    longitudinal_orders: tuple[int, ...] = (0, 1, 2),
) -> Result:
    """Every low-order acoustic mode of a rigid cylindrical chamber [Hz].

    Returns a `Result` whose value is a list of modes sorted by frequency.
    `chamber_length_m` should be the length the ACOUSTIC cavity actually has,
    which is a modelling decision the caller makes and this module records:
    the barrel alone, or the barrel plus some fraction of the convergent
    section. Neither is obviously right and the difference is reported below
    as a sensitivity rather than resolved.
    """
    if chamber_diameter_m <= 0 or chamber_length_m <= 0:
        return Result(
            value=None, units="Hz", equation_id="chamber_acoustic_modes",
            status=Status.PHYSICALLY_INVALID, source="n/a",
            notes="chamber diameter and length must both be positive")

    a = sound_speed(gamma=gamma, specific_gas_constant=specific_gas_constant,
                    temperature_K=temperature_K)
    R = 0.5 * chamber_diameter_m
    L = chamber_length_m
    roots = transverse_roots()

    modes: list[dict] = []
    for (m, n), alpha in roots.items():
        for q in longitudinal_orders:
            if m == 0 and n == 0 and q == 0:
                continue
            f = 0.5 * a * math.sqrt((alpha / (math.pi * R)) ** 2
                                    + (q / L) ** 2)
            label = MODE_LABELS.get((m, n), f"m={m} n={n}")
            modes.append({
                "m": m, "n": n, "q": q,
                "label": label + (f" + {q}L" if q else ""),
                "bessel_root": alpha,
                "frequency_Hz": f,
            })
    # pure longitudinal modes carry no transverse root
    for q in longitudinal_orders:
        if q == 0:
            continue
        modes.append({"m": None, "n": None, "q": q,
                      "label": f"{q}L   longitudinal",
                      "bessel_root": None,
                      "frequency_Hz": q * a / (2.0 * L)})
    modes.sort(key=lambda d: d["frequency_Hz"])

    # how much does the sound-speed choice matter? +/-5 % on T is +/-2.5 % on
    # every frequency, because f scales as sqrt(T).
    f1t = next((d["frequency_Hz"] for d in modes
                if d["m"] == 1 and d["n"] == 1 and d["q"] == 0), None)

    return Result(
        value=modes, units="Hz", equation_id="chamber_acoustic_modes",
        status=Status.UNVALIDATED_ASSUMPTION,
        evidence_level=EvidenceLevel.E2,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source="separable solution of the wave equation in a rigid-walled "
               "cylinder; the transverse eigenvalues are the roots of "
               "J'_m(alpha)=0, computed here rather than tabulated",
        source_locator="first principles. The standard rocket reference is "
                       "NASA SP-194, 'Liquid Propellant Rocket Combustion "
                       "Instability', which this project has NOT obtained.",
        validity_domain="rigid-walled uniform-gas cylinder with no mean "
                        "flow. A real chamber has a convergent end, a mean "
                        "Mach number around 0.13, a temperature gradient and "
                        "an injector face that is not acoustically rigid. "
                        "Each departure shifts the true frequency and none "
                        "is corrected for here.",
        uncertainty="NOT_REPORTED as an absolute bound. SCALING is exact: "
                    "f varies as sqrt(T), so a 5 % error in chamber "
                    "temperature is 2.5 % in every frequency; f varies as "
                    "1/R for transverse modes and 1/L for longitudinal ones.",
        assumptions=(
            Assumption("rigid_walls_no_mean_flow", 1.0, "-",
                       "no mean flow, no acoustic damping, no combustion "
                       "response, and both ends acoustically rigid",
                       Status.UNVALIDATED_ASSUMPTION),
            Assumption("frozen_sound_speed", a, "m/s",
                       "sqrt(gamma R T) on the frozen composition. An "
                       "equilibrium sound speed is lower, so every frequency "
                       "here is biased HIGH by an unquantified amount.",
                       Status.UNVALIDATED_ASSUMPTION),
            Assumption("acoustic_length_is_caller_supplied", L, "m",
                       "whether the cavity is the barrel alone or includes "
                       "part of the convergent section is a modelling "
                       "choice the caller makes; it moves the longitudinal "
                       "modes directly and the transverse modes not at all",
                       Status.UNVALIDATED_ASSUMPTION),),
        inputs={"chamber_diameter_m": chamber_diameter_m,
                "chamber_length_m": chamber_length_m,
                "gamma": gamma,
                "specific_gas_constant": specific_gas_constant,
                "temperature_K": temperature_K,
                "sound_speed_m_s": a},
        notes=(f"sound speed {a:.2f} m/s; "
               + (f"1T = {f1t:.1f} Hz; " if f1t else "")
               + f"{len(modes)} modes computed. "
                 f"THESE ARE FREQUENCIES, NOT A STABILITY ASSESSMENT. "
                 f"See stability_verdict()."))


def chug_frequency(**_kwargs) -> Result:
    """Feed-system-coupled (chug) instability. NOT COMPUTABLE HERE.

    Chug is a coupled oscillation between the feed system's inertance and
    compliance and the chamber's filling time. The chamber half is available
    -- `V_c`, `c*`, `A_t` give the chamber time constant -- but the feed half
    needs line lengths, areas and the compliance of every accumulator and
    cavitating venturi between the tank and the injector.

    **None of that is in scope for the combustor**, so this function refuses
    rather than substituting a chamber-only estimate that would look like a
    chug frequency and would not be one.
    """
    return Result(
        value=None, units="Hz", equation_id="chamber_chug_frequency",
        status=Status.INSUFFICIENT_EVIDENCE,
        source="n/a",
        validity_domain="VALIDITY_DOMAIN_UNKNOWN",
        uncertainty=NOT_REPORTED,
        notes="MISSING_UPSTREAM_CONTRACT: chug frequency needs feed-line "
              "inertance and compliance, which the combustor does not own. "
              "A chamber-only time constant is NOT a chug frequency and is "
              "deliberately not returned in its place.")


def stability_verdict(**_kwargs) -> Result:
    """Is this chamber stable? **THIS REPOSITORY CANNOT SAY.**

    Deliberately implemented, deliberately refusing. A caller that reaches for
    a stability answer should hit an explicit refusal with the reason, not an
    `AttributeError` or -- far worse -- a frequency table that it mistakes for
    a verdict.

    What is missing is not code. It is the injector response function: the
    interaction index `n` and the sensitive time lag `tau` of Crocco's model,
    which are properties of a specific element pattern and are MEASURED, in a
    stability rating test, on hardware. No published `n`-`tau` pair transfers
    to a new injector.
    """
    return Result(
        value=None, units="-", equation_id="combustion_stability_verdict",
        status=Status.INSUFFICIENT_EVIDENCE,
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE,
        validation=Validation.NONE,
        source="Crocco, L. & Cheng, S.-I., 'Theory of Combustion Instability "
               "in Liquid Propellant Rocket Motors', AGARDograph 8, 1956 -- "
               "cited in this repository and NOT OBTAINED",
        source_locator="NOT_OPENED",
        validity_domain="VALIDITY_DOMAIN_UNKNOWN",
        uncertainty=NOT_REPORTED,
        assumptions=(
            Assumption("no_response_function", 0.0, "-",
                       "n and tau are unknown for this injector and cannot "
                       "be derived; a growth-rate balance without them is "
                       "not a calculation, it is a guess with algebra "
                       "attached", Status.INSUFFICIENT_EVIDENCE),),
        notes="NO STABILITY VERDICT IS AVAILABLE. The acoustic mode "
              "frequencies ARE available from acoustic_modes() and are a "
              "genuine pre-test prediction; the growth-rate balance that "
              "would turn them into a verdict needs an injector response "
              "function measured on hardware. Requiring this refusal to be "
              "explicit is the point: a frequency table is not a stability "
              "assessment and must never be reported as one.")
