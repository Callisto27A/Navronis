r"""
THE canonical coolant-side model for supercritical methane.
===========================================================

This module exists because the thermal wave measured a **4.7x spread** across
four coolant-side correlations on one channel state, and none of them was
defensible end to end. The research that followed did not close that spread by
picking a winner. It found something more useful and more uncomfortable:

    THE SPREAD WAS BEING MEASURED IN THE WRONG PLACE.

All four candidates -- McAdams, Dittus-Boelter, Gnielinski, Taylor -- are
SMOOTH-TUBE correlations. Additively manufactured rocket cooling channels are
not smooth and are not close to smooth:

    Hartsfield et al. (2021), LPBF Inconel 718, MEASURED:
        Ra    = 7.48 - 14.28 um
        eps   = 15.87 - 21.24 um   (equivalent sand-grain, from pressure drop)
        eps/D = 0.013 - 0.105

That is one to two orders of magnitude above the eps/D ~ 1e-4 where a pipe
stays hydraulically smooth at Re ~ 1e5. Sapienza's LRE cooling group (Latini,
Fiore & Nasuti) accordingly models these channels as FULLY ROUGH -- Colebrook
for friction, **Dipprey & Sabersky** for heat transfer -- and not with
Gnielinski at all.

So the ordering of the open questions is:

    smooth vs rough            TENS OF PERCENT
    Ra -> eps conversion       0.4x to 16x, sources structurally incompatible
    Filonenko vs Konakov       ~1.1% at Re = 1e5

The wave began by trying to settle the third. It is settled (see below), and it
turned out not to matter much. The first two now dominate, and only one of them
has an answer.


WHAT THE RESEARCH SETTLED
-------------------------

1.  **THE GNIELINSKI FRICTION-FACTOR "CONFLICT" WAS NOT A CONFLICT.**
    Gnielinski used *both*, in different publications, and the secondary
    sources were each quoting him correctly:

        Forsch. Ing.-Wes. 41(1) 1975 (the German original of the 1976
        Int. Chem. Eng. paper) -- Schrifttum [34] cites FILONENKO, and
        Konakov does not appear in that reference list at all.

        VDI Heat Atlas 2nd Eng. ed. 2010, Ch. G1, written by Gnielinski --
        reference [17] cites KONAKOV, and Filonenko does not appear in
        that reference list at all.

    Kryptonis labels its implementation "Gnielinski (1976)". The correct
    companion for that label is therefore **Filonenko**. The mismatch to fix
    was the LABEL, not the number.

2.  **THERE IS A SOURCE-BACKED HTD ONSET CRITERION**, and it is a threshold on
    the heat-flux-to-mass-flux ratio, not a dimensionless group:

        (q_w / G)_tr = 43.2e-6 * p_in + 31.4      [J/kg], p_in in [Pa]

    Urbano & Nasuti (2013). **NOTE THE ATTRIBUTION**: DOI 10.2514/1.T4001 is
    Urbano, A. & Nasuti, F., J. Thermophysics and Heat Transfer 27(2),
    pp.298-308, 2013 -- it is NOT a Pizzarelli paper, and this project's own
    task brief had it wrong. Confirmed against Crossref, OpenAlex and the
    reference list of Haemisch et al. (2019).

3.  **AND THERE IS EXPERIMENTAL METHANE HTD DATA.** Haemisch, Suslov &
    Oschwald, J. Propulsion & Power 35(4), 819-826 (2019), open at DLR elib.
    Subscale LOX/CH4 chamber, 126 channels, inverse heat-conduction from 80
    thermocouples. Their measured result, verbatim: HTD occurs *"for pressures
    below P_in ~ 70 bar, that is P_in/P_cr ~ 1.5"*, and *"for the higher
    pressure test cases, P_in >= 100 bar, no peak is observable, and thus they
    represent a normal heat transfer."*

4.  **HTD IS NOT A BUOYANCY PHENOMENON AT ROCKET MASS FLUXES.** Haemisch's
    hardware runs at ~1e4 kg/(m^2.s) and his paper never mentions buoyancy or
    Grashof number. The mechanism he describes is a low-density gas-like
    near-wall layer acting as a thermal barrier, plus the flow acceleration it
    induces (an M-shaped velocity profile). Pizzarelli's own 2016 abstract
    states his correlation is for deteriorated flow *"with negligible buoyancy
    effects."* Buoyancy-driven HTD appears in the literature only at
    low-mass-flux vertical tubes, two orders of magnitude below rocket
    conditions.

    **NO opened source gave an exact definition of a buoyancy parameter or an
    acceleration parameter with a threshold for supercritical methane.** So
    this module implements the q/G threshold and NOT a dimensionless-group
    criterion, because that is what the evidence supports.

5.  **Schacht & Quentmeyer's coefficient is 0.023, not 0.025.** NASA
    TN D-7207 was opened. Its own new correlation is
    ``St* Pr*^0.6 = 0.023 Re*^-0.2``. The string "0.025" appears in that
    report's INTRODUCTION, describing the prior practice it set out to
    improve on. Locke & Landrum (2008) are widely cited as giving 0.025;
    that paper could not be opened to determine whether it is a
    mis-transcription or a deliberate re-fit. **8.7% is at stake and it is
    recorded as CONFLICT_UNRESOLVED.**


WHAT THE RESEARCH DID *NOT* SETTLE
----------------------------------

*   **Ra -> equivalent sand-grain roughness.** Three open primary sources give
    structurally incompatible answers spanning 0.4x to 16x. Converting a
    metrology Ra into a Colebrook eps introduces an uncertainty larger than
    every other uncertainty in this chain, and not one of the three sources
    quantifies it. `sand_grain_from_Ra` returns CONFLICT_UNRESOLVED with all
    three candidates and no choice.

*   **The deteriorated heat-transfer coefficient.** The q/G criterion is a
    BINARY ONSET FLAG. It does not tell you what h_c becomes once HTD starts.
    The correlation that would (Pizzarelli, Numer. Heat Transfer A 69(3),
    2016) is behind a 403 and **zero coefficients from it are available**.
    So `htd_onset` can say "you are in trouble here" and cannot say how much.

*   **Any rocket-relevant validation of Dipprey & Sabersky.** It was fitted on
    DISTILLED WATER at Pr 1.20-5.94, eps/D <= 0.0488. Sapienza apply it at
    eps/D up to 0.21 -- a 4x extrapolation -- and state no heat-transfer error
    figure for it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

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
    "COOLPROP_AVAILABLE", "PropertySource", "PROPERTY_PROVENANCE",
    "CoolantState", "coolant_state", "pseudo_critical_temperature",
    "enthalpy_to_temperature", "coolant_regime_map",
    "URBANO_NASUTI_2013", "HAEMISCH_2019", "htd_onset",
    "RoughnessRegime", "roughness_regime", "sand_grain_from_Ra",
    "friction_factor_filonenko", "friction_factor_konakov",
    "GnielinskiVintage", "gnielinski_friction_factor",
    "nusselt_petukhov_1970", "nusselt_dipprey_sabersky",
    "nusselt_stimpson_am", "nusselt_schacht_quentmeyer",
    "CoolantCorrelationChoice", "select_coolant_correlation",
    "coolant_side_h_c", "friction_factor_from_pressure_drop",
    "hydraulic_roughness_from_pressure_drop",
    "normal_supercritical_ht", "deteriorated_supercritical_ht",
    "CoolantHeatTransfer", "coolant_heat_transfer",
    "uncertainty_decomposition",
]

try:                                          # pragma: no cover - env dependent
    from CoolProp.CoolProp import PropsSI as _PropsSI
    COOLPROP_AVAILABLE = True
except Exception:                             # pragma: no cover
    _PropsSI = None
    COOLPROP_AVAILABLE = False


# ===========================================================================
# SOURCES
# ===========================================================================

SETZMANN_WAGNER_1991 = (
    "Setzmann, U. & Wagner, W., 'A New Equation of State and Tables of "
    "Thermodynamic Properties for Methane Covering the Range from the Melting "
    "Line to 625 K at Pressures up to 1000 MPa', J. Phys. Chem. Ref. Data "
    "20(6), 1991, p.1061 ff., DOI 10.1063/1.555898")

FRIEND_ELY_INGHAM_1989 = (
    "Friend, D. G., Ely, J. F. & Ingham, H., 'Thermophysical Properties of "
    "Methane', J. Phys. Chem. Ref. Data 18(2), 1989, pp.583-638, "
    "DOI 10.1063/1.555828")

QUINONES_DEITERS_2006 = (
    "Quinones-Cisneros, S. E. & Deiters, U. K., 'Generalization of the "
    "Friction Theory for Viscosity Modeling', J. Phys. Chem. B 110, 2006, "
    "pp.12820-12834 -- a GENERALIZED model, not a methane-specific "
    "reference correlation")

URBANO_NASUTI_2013 = (
    "Urbano, A. & Nasuti, F., 'Onset of Heat Transfer Deterioration in "
    "Supercritical Methane Flow Channels', J. Thermophysics and Heat Transfer "
    "27(2), 2013, pp.298-308, DOI 10.2514/1.T4001")

HAEMISCH_2019 = (
    "Haemisch, J., Suslov, D. & Oschwald, M., 'Experimental Study of Methane "
    "Heat Transfer Deterioration in a Subscale Combustion Chamber', "
    "J. Propulsion and Power 35(4), 2019, pp.819-826, DOI 10.2514/1.B37394")

DIPPREY_SABERSKY_1962 = (
    "Dipprey, D. F. & Sabersky, R. H., 'Heat and Momentum Transfer in Smooth "
    "and Rough Tubes at Various Prandtl Numbers', JPL Technical Report "
    "No. 32-269, Caltech, 6 June 1962 (NTRS 19630000458); journal version "
    "Int. J. Heat Mass Transfer 6, 1963, pp.329-353")

NIKURADSE_1933 = (
    "Nikuradse, J., 'Laws of Flow in Rough Pipes', NACA TM 1292 (translation "
    "of VDI-Forschungsheft 361, 1933), regime boundaries at Eq. group (3.21)")

STIMPSON_2017 = (
    "Stimpson, C. K., Snyder, J. C., Thole, K. A. & Mongillo, D., 'Roughness "
    "Effects on Flow and Heat Transfer for Additively Manufactured Channels', "
    "J. Turbomachinery 139(2), 2017, 021003 (ASME GT2016-58093)")

HARTSFIELD_2021 = (
    "Hartsfield, C. R., Shelton, T. E., Cobb, R. G., Kemnitz, R. A. & "
    "Weber, R. W., J. Aerospace Engineering, 2021, "
    "DOI 10.1061/(ASCE)AS.1943-5525.0001325 -- measured LPBF Inconel 718 "
    "internal-channel roughness")

FLACK_SCHULTZ_2010 = (
    "Flack, K. A. & Schultz, M. P., 'Review of Hydraulic Roughness Scales in "
    "the Fully Rough Regime', J. Fluids Engineering 132(4), 2010, 041203")

ADAMS_2012 = (
    "Adams, T., Grant, C. & Watson, H., 'A Simple Algorithm to Relate "
    "Measured Surface Roughness to Equivalent Sand-grain Roughness', "
    "Int. J. Mechanical Engineering and Mechatronics 1(1), 2012, pp.66-71")

SCHACHT_QUENTMEYER_1973 = (
    "Schacht, R. L. & Quentmeyer, R. J., 'Coolant-Side Heat-Transfer Rates "
    "for a Hydrogen-Oxygen Rocket and a New Technique for Data Correlation', "
    "NASA TN D-7207, NASA Lewis Research Center, March 1973 "
    "(NTRS 19730010241)")

PETUKHOV_1970 = (
    "Petukhov, B. S., 'Heat Transfer and Friction in Turbulent Pipe Flow with "
    "Variable Physical Properties', Advances in Heat Transfer 6, 1970, p.503 "
    "ff., Eqs. (48) and (49)")

FILONENKO_1954 = (
    "Filonenko, G. K., 'Hydraulic resistance in pipes' (orig. Russian), "
    "Teploenergetika 1(4), 1954, pp.40-44")

KONAKOV_1946 = (
    "Konakov, P. K., 'A new correlation for the friction coefficient in "
    "smooth tubes', Berichte der Akademie der Wissenschaften der UdSSR "
    "51(7), 1946, pp.503-506")

COLEBROOK_1939 = (
    "Colebrook, C. F., 'Turbulent Flow in Pipes, with Particular Reference "
    "to the Transition Region Between the Smooth and Rough Pipe Laws', "
    "J. Inst. Civil Engineers 11(4), 1938-39, pp.133-156 -- PRIMARY "
    "NOT_OPENED (ICE 403; Emerald CAPTCHA wall)")

LATINI_2022 = (
    "Latini, B., Fiore, M. & Nasuti, F., 'Analysis of coolant flow and heat "
    "transfer in highly rough channels for LRE', EUCASS 2022, paper 6138, "
    "DOI 10.13009/EUCASS2022-6138")


# ===========================================================================
# 1. REAL-FLUID PROPERTIES -- provenance PER PROPERTY, not per library
# ===========================================================================

class PropertySource(str, Enum):
    """How a property was obtained. `IDEAL_GAS` is a failure mode here."""

    REFERENCE_EOS = "REFERENCE_EOS"
    REFERENCE_TRANSPORT = "REFERENCE_TRANSPORT"
    GENERALIZED_MODEL = "GENERALIZED_MODEL"
    DERIVED = "DERIVED"
    IDEAL_GAS = "IDEAL_GAS"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class _Prov:
    source: str
    locator: str
    kind: PropertySource
    uncertainty: str
    domain: str
    caveat: str = ""


#: What CoolProp ACTUALLY uses for methane, property by property. This matters:
#: a single "properties from CoolProp" label would hide that the viscosity
#: comes from a generalized friction-theory model rather than a
#: methane-specific reference correlation, and that the thermal conductivity's
#: critical enhancement is only valid in a 11 K window that a regenerative
#: channel runs straight through.
PROPERTY_PROVENANCE: dict[str, _Prov] = {
    "rho": _Prov(
        SETZMANN_WAGNER_1991, "reference EOS; CoolProp HEOS backend",
        PropertySource.REFERENCE_EOS,
        "+/-0.03% for p < 12 MPa and T < 350 K; +/-0.03% to +/-0.15% for "
        "higher p and T (author's own statement)",
        "melting line (90.6941 K) to 625 K, p to 1000 MPa"),
    "h": _Prov(
        SETZMANN_WAGNER_1991, "reference EOS; CoolProp HEOS backend",
        PropertySource.REFERENCE_EOS,
        "NOT_REPORTED as a separate figure; enthalpy is an integral of the "
        "EOS and inherits its density and c_p uncertainties",
        "melting line to 625 K, p to 1000 MPa"),
    "cp": _Prov(
        SETZMANN_WAGNER_1991, "reference EOS; CoolProp HEOS backend",
        PropertySource.REFERENCE_EOS,
        "+/-1% generally (author's own statement)",
        "melting line to 625 K, p to 1000 MPa",
        "the EOS is COMPLETELY ANALYTIC and does NOT reproduce the "
        "theoretical asymptotic scaling behaviour at the critical point. c_p "
        "peaks near the pseudo-critical line are represented well against "
        "data but carry more than the +/-1% headline there."),
    "mu": _Prov(
        QUINONES_DEITERS_2006, "CoolProp's declared methane viscosity model",
        PropertySource.GENERALIZED_MODEL,
        "NOT_REPORTED. CoolProp states NO uncertainties anywhere in its "
        "documentation -- checked the Methane fluid page, the High-Level API "
        "page and the PurePseudoPure page.",
        UNKNOWN_DOMAIN,
        "THIS IS THE WEAKEST LINK IN THE PROPERTY CHAIN. CoolProp does NOT "
        "use a methane-specific reference viscosity correlation; it uses a "
        "GENERALIZED friction-theory model. REFPROP and the 2026 Sotiriadou "
        "reference correlation would both be better. Do not quote a "
        "viscosity uncertainty for this path."),
    "k": _Prov(
        FRIEND_ELY_INGHAM_1989, "CoolProp's declared methane conductivity "
                                "model",
        PropertySource.REFERENCE_TRANSPORT,
        "about 1.5% (author's own statement)",
        "about 91 to 700 K for p < 100 MPa",
        "the CRITICAL ENHANCEMENT term is stated valid only for "
        "185 K < T < 196 K and 7.6 < rho < 12.7 mol/dm^3 (122-204 kg/m^3), "
        "and NO separate uncertainty is stated inside that window. "
        "MEASURED against the reference EOS, that window is reachable "
        "only for p <~ 5.7 MPa (p/p_crit <~ 1.24) -- above it, methane at "
        "185-196 K is denser than 204 kg/m^3. So the conductivity is "
        "least trustworthy in exactly the near-critical band where "
        "Haemisch measured heat-transfer deterioration, and a 9 MPa "
        "channel does NOT pass through it."),
    "Pr": _Prov(
        "derived: Pr = mu c_p / k", "definitional",
        PropertySource.DERIVED,
        "inherits mu (NOT_REPORTED), c_p (+/-1%) and k (1.5%); the mu term "
        "therefore makes the Prandtl uncertainty UNKNOWN",
        "intersection of the three constituents' domains"),
}

#: Friend, Ely & Ingham's critical-enhancement window for thermal
#: conductivity, converted from mol/dm^3 using M = 0.0160428 kg/mol.
_K_CRIT_ENHANCEMENT_T_K = (185.0, 196.0)
_K_CRIT_ENHANCEMENT_RHO_KG_M3 = (7.6e3 * 0.0160428, 12.7e3 * 0.0160428)

#: Setzmann & Wagner's stated EOS domain.
_EOS_T_RANGE_K = (90.6941, 625.0)
_EOS_P_MAX_PA = 1.0e9


@dataclass
class CoolantState:
    """One thermodynamic state, with per-property provenance attached."""

    fluid: str
    pressure_Pa: float
    temperature_K: float
    density_kg_m3: float
    enthalpy_J_kg: float
    specific_heat_J_kgK: float
    viscosity_Pa_s: float
    conductivity_W_mK: float
    prandtl: float
    result: Result
    statuses: frozenset = field(default_factory=frozenset)


def _require_coolprop(equation_id: str, inputs: dict) -> Result | None:
    if not COOLPROP_AVAILABLE:
        return Result(
            float("nan"), "-", equation_id,
            status=Status.INSUFFICIENT_EVIDENCE,
            evidence_level=EvidenceLevel.E0,
            source="NO PROPERTY BACKEND",
            notes="CoolProp is not importable, so real-fluid methane "
                  "properties are unavailable. An ideal-gas substitute is "
                  "NOT provided: near the pseudo-critical line methane's c_p "
                  "varies by a factor of ~3.5 and its density by a factor of "
                  "~4, so ideal-gas methane is not an approximation of this "
                  "state, it is a different fluid.",
            inputs=inputs)
    return None


def coolant_state(*, fluid: str = "Methane", pressure_Pa: float,
                  temperature_K: float) -> CoolantState:
    """Full real-fluid state, with every property's provenance recorded.

    There is deliberately no ideal-gas fallback. Between 150 K and 250 K at
    9 MPa, methane's density falls from 374 to 98 kg/m^3 and its c_p peaks at
    7.7 kJ/(kg.K) against a far-field ~2.2 -- an ideal-gas evaluation there is
    not a rough answer, it is the wrong fluid.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa, "T_K": temperature_K}
    bad = _require_coolprop("CAN-COOL-STATE", inputs)
    if bad is not None:
        n = float("nan")
        return CoolantState(fluid, pressure_Pa, temperature_K, n, n, n, n, n,
                            n, bad, bad.all_statuses)

    if pressure_Pa <= 0 or temperature_K <= 0:
        r = Result(float("nan"), "-", "CAN-COOL-STATE",
                   status=Status.PHYSICALLY_INVALID, inputs=inputs,
                   notes="pressure and temperature must be > 0")
        n = float("nan")
        return CoolantState(fluid, pressure_Pa, temperature_K, n, n, n, n, n,
                            n, r, r.all_statuses)

    statuses: list[Status] = [Status.PASS]
    notes: list[str] = []
    assumptions: list[Assumption] = []

    lo, hi = _EOS_T_RANGE_K
    if not (lo <= temperature_K <= hi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"T = {temperature_K:.1f} K is outside the reference "
                     f"EOS domain {lo:.1f}-{hi:.0f} K")
    if pressure_Pa > _EOS_P_MAX_PA:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append("p exceeds the reference EOS domain of 1000 MPa")

    try:
        rho = _PropsSI("D", "T", temperature_K, "P", pressure_Pa, fluid)
        h = _PropsSI("Hmass", "T", temperature_K, "P", pressure_Pa, fluid)
        cp = _PropsSI("Cpmass", "T", temperature_K, "P", pressure_Pa, fluid)
        mu = _PropsSI("V", "T", temperature_K, "P", pressure_Pa, fluid)
        k = _PropsSI("L", "T", temperature_K, "P", pressure_Pa, fluid)
    except Exception as exc:
        r = Result(float("nan"), "-", "CAN-COOL-STATE",
                   status=Status.PHYSICALLY_INVALID, inputs=inputs,
                   notes=f"the property backend refused this state: {exc}")
        n = float("nan")
        return CoolantState(fluid, pressure_Pa, temperature_K, n, n, n, n, n,
                            n, r, r.all_statuses)

    pr = mu * cp / k

    # The viscosity model is generalized, not methane-specific. That is a
    # standing assumption of every number derived from it, including Pr and
    # therefore every Nusselt number in this module.
    assumptions.append(Assumption(
        "viscosity_model", "Quinones-Cisneros & Deiters friction theory",
        "-", "CoolProp's methane viscosity is a GENERALIZED model, not a "
             "methane-specific reference correlation, and CoolProp states no "
             "uncertainty for it",
        Status.UNVALIDATED_ASSUMPTION))
    statuses.append(Status.UNVALIDATED_ASSUMPTION)

    klo, khi = _K_CRIT_ENHANCEMENT_T_K
    rlo, rhi = _K_CRIT_ENHANCEMENT_RHO_KG_M3
    if klo <= temperature_K <= khi and rlo <= rho <= rhi:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(
            f"INSIDE the thermal-conductivity critical-enhancement window "
            f"({klo:.0f}-{khi:.0f} K and {rlo:.0f}-{rhi:.0f} kg/m^3). Friend, "
            f"Ely & Ingham state NO uncertainty for k there, so the 1.5% "
            f"headline does not apply and Pr inherits that.")

    r = Result(
        temperature_K, "K", "CAN-COOL-STATE", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E2,
        verification=Verification.INDEPENDENT_IMPLEMENTATION,
        validation=Validation.NONE,
        source=f"{SETZMANN_WAGNER_1991} (thermodynamic); "
               f"{FRIEND_ELY_INGHAM_1989} (conductivity); "
               f"{QUINONES_DEITERS_2006} (viscosity)",
        source_locator="all three OPENED; see PROPERTY_PROVENANCE for the "
                       "per-property locator, uncertainty and caveat",
        validity_domain=f"EOS {lo:.1f}-{hi:.0f} K to 1000 MPa; k 91-700 K "
                        f"below 100 MPa",
        uncertainty="rho +/-0.03-0.15%; c_p +/-1%; k 1.5%; mu NOT_REPORTED "
                    "(generalized model). Pr therefore has UNKNOWN "
                    "uncertainty.",
        assumptions=tuple(assumptions), notes="; ".join(notes),
        inputs={**inputs, "rho_kg_m3": rho, "h_J_kg": h, "cp_J_kgK": cp,
                "mu_Pa_s": mu, "k_W_mK": k, "Pr": pr})

    return CoolantState(fluid, pressure_Pa, temperature_K, rho, h, cp, mu, k,
                        pr, r, frozenset(statuses))


def pseudo_critical_temperature(*, fluid: str = "Methane",
                                pressure_Pa: float,
                                search_K: tuple[float, float] = (100.0, 400.0),
                                tolerance_K: float = 0.01) -> Result:
    r"""The pseudo-critical (Widom-line) temperature: where c_p is maximal.

    Above the critical pressure there is no phase boundary, but there IS a
    locus of maximum c_p -- the Widom line -- and crossing it inside the
    boundary layer is the documented condition associated with heat-transfer
    deterioration. It is a property of the fluid at that pressure, located
    here by a golden-section search on c_p rather than taken from a table,
    because no source consulted tabulates T_pc for methane at arbitrary
    pressure.

    Below p_crit this is not defined -- there is a real saturation temperature
    instead -- and the result says so rather than returning the boiling point
    under a misleading name.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa}
    bad = _require_coolprop("CAN-COOL-TPC", inputs)
    if bad is not None:
        return bad
    if pressure_Pa <= 0:
        return Result(float("nan"), "K", "CAN-COOL-TPC",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="pressure must be > 0")
    try:
        p_crit = _PropsSI("Pcrit", fluid)
    except Exception as exc:
        return Result(float("nan"), "K", "CAN-COOL-TPC",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes=f"no critical pressure for {fluid!r}: {exc}")
    if pressure_Pa < p_crit:
        return Result(
            float("nan"), "K", "CAN-COOL-TPC",
            status=Status.INVALID_MODEL_REGIME,
            all_statuses=frozenset({Status.INVALID_MODEL_REGIME}),
            evidence_level=EvidenceLevel.E1,
            source=SETZMANN_WAGNER_1991,
            validity_domain=f"p >= p_crit = {p_crit / 1e6:.4f} MPa",
            notes=f"p = {pressure_Pa / 1e6:.3f} MPa is BELOW p_crit. There is "
                  f"no pseudo-critical temperature below the critical "
                  f"pressure -- there is a saturation temperature, which is a "
                  f"different quantity with a phase boundary attached. This "
                  f"function will not return one under the wrong name.",
            inputs={**inputs, "p_crit_Pa": p_crit})

    def _cp(t: float) -> float:
        try:
            return _PropsSI("Cpmass", "T", t, "P", pressure_Pa, fluid)
        except Exception:
            return -1.0

    a, b = search_K
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = _cp(c), _cp(d)
    for _ in range(200):
        if abs(b - a) < tolerance_K:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = _cp(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = _cp(d)
    t_pc = 0.5 * (a + b)
    cp_peak = _cp(t_pc)
    if cp_peak <= 0:
        return Result(float("nan"), "K", "CAN-COOL-TPC",
                      status=Status.FAIL, inputs=inputs,
                      notes="the c_p search did not find a maximum")

    statuses = [Status.PASS]
    notes = [f"c_p peak = {cp_peak / 1000.0:.2f} kJ/(kg.K)"]
    if not (search_K[0] + tolerance_K < t_pc < search_K[1] - tolerance_K):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append("the maximum sits on the edge of the search interval; "
                     "widen `search_K`")

    return Result(
        t_pc, "K", "CAN-COOL-TPC", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=SETZMANN_WAGNER_1991,
        source_locator="located numerically as the c_p maximum from the "
                       "reference EOS; NOT tabulated by any source consulted",
        validity_domain=f"p >= p_crit = {p_crit / 1e6:.4f} MPa",
        uncertainty="inherits the EOS c_p uncertainty (+/-1% nominal, larger "
                    "near the critical point where the analytic EOS does not "
                    "reproduce asymptotic scaling)",
        notes="; ".join(notes),
        inputs={**inputs, "cp_peak_J_kgK": cp_peak,
                "p_over_p_crit": pressure_Pa / p_crit})


def enthalpy_to_temperature(*, fluid: str = "Methane", pressure_Pa: float,
                            enthalpy_J_kg: float) -> Result:
    r"""Invert h(p,T) -> T(p,h). SOURCE-BACKED, not a constant-c_p shortcut.

    This is what the thermal wave's `energy_balance` needed and did not have.
    The enthalpy form ``h_out = h_in + Q/mdot`` is exact; without this
    inversion the outlet TEMPERATURE was `INSUFFICIENT_EVIDENCE`.

    The inversion is performed by the reference equation of state itself, not
    by a local linearisation, so it is valid across the pseudo-critical line
    where ``dT = Q/(mdot c_p)`` is not merely inaccurate but structurally
    wrong -- c_p there varies by a factor of ~3.5 over a few tens of kelvin.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa, "h_J_kg": enthalpy_J_kg}
    bad = _require_coolprop("CAN-COOL-H2T", inputs)
    if bad is not None:
        return bad
    if pressure_Pa <= 0:
        return Result(float("nan"), "K", "CAN-COOL-H2T",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="pressure must be > 0")
    try:
        t = _PropsSI("T", "P", pressure_Pa, "Hmass", enthalpy_J_kg, fluid)
    except Exception as exc:
        return Result(
            float("nan"), "K", "CAN-COOL-H2T",
            status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
            evidence_level=EvidenceLevel.E0,
            source=SETZMANN_WAGNER_1991,
            notes=f"the reference EOS could not invert this (p, h) state: "
                  f"{exc}. NO fallback linearisation is applied -- a "
                  f"constant-c_p estimate here would be a different equation "
                  f"wearing this function's name.")

    statuses = [Status.PASS]
    notes: list[str] = []
    lo, hi = _EOS_T_RANGE_K
    if not (lo <= t <= hi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"the inverted temperature {t:.1f} K is outside the "
                     f"reference EOS domain {lo:.1f}-{hi:.0f} K")

    return Result(
        t, "K", "CAN-COOL-H2T", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E2,
        verification=Verification.INDEPENDENT_IMPLEMENTATION,
        validation=Validation.NONE,
        source=SETZMANN_WAGNER_1991,
        source_locator="inverted by the reference EOS; the forward and "
                       "inverse are the same equation of state",
        validity_domain=f"{lo:.1f}-{hi:.0f} K, p to 1000 MPa",
        uncertainty="inherits the EOS enthalpy uncertainty; NOT_REPORTED as "
                    "a separate figure by the source",
        notes="; ".join(notes), inputs=inputs)


def coolant_regime_map(*, fluid: str = "Methane", pressure_Pa: float,
                       bulk_temperature_K: float,
                       wall_temperature_K: float | None = None) -> Result:
    """Regime by the FULL state, not by pressure alone.

    The thermal wave's `supercritical_regime` classified on p/p_crit and a
    T_crit crossing. That is not enough: above the critical pressure the
    relevant locus is the PSEUDO-critical temperature, which at 9 MPa sits at
    214 K -- 24 K above T_crit -- and moves with pressure. Classifying on
    T_crit alone mislabels the band between the two.

    Returns one of ``SUBCRITICAL`` / ``LIQUID_LIKE`` / ``PSEUDO_CRITICAL`` /
    ``GAS_LIKE`` / ``NEAR_CRITICAL``, and puts every ratio the literature
    uses into `inputs` so a caller can apply its own criterion.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "T_wall_K": wall_temperature_K}
    bad = _require_coolprop("CAN-COOL-REGIME", inputs)
    if bad is not None:
        return bad
    if pressure_Pa <= 0 or bulk_temperature_K <= 0:
        return Result("INVALID", "-", "CAN-COOL-REGIME",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="pressure and temperature must be > 0")
    try:
        p_c = _PropsSI("Pcrit", fluid)
        t_c = _PropsSI("Tcrit", fluid)
    except Exception as exc:
        return Result("UNKNOWN", "-", "CAN-COOL-REGIME",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      evidence_level=EvidenceLevel.E0,
                      notes=f"no critical point for {fluid!r}: {exc}. The "
                            f"regime is NOT assumed subcritical.")

    p_r = pressure_Pa / p_c
    t_r = bulk_temperature_K / t_c
    statuses = [Status.PASS]
    notes = [f"p/p_crit = {p_r:.3f}, T_bulk/T_crit = {t_r:.3f}"]
    extra: dict[str, Any] = {"p_over_p_crit": p_r, "T_over_T_crit": t_r,
                             "p_crit_Pa": p_c, "T_crit_K": t_c}

    if pressure_Pa < p_c:
        label = "SUBCRITICAL"
        notes.append("below p_crit a phase boundary exists; the "
                     "pseudo-critical concept does not apply and boiling "
                     "models do")
        return Result(
            label, "-", "CAN-COOL-REGIME", status=Status.PASS,
            evidence_level=EvidenceLevel.E1,
            source=SETZMANN_WAGNER_1991, source_locator="critical point",
            notes="; ".join(notes), inputs={**inputs, **extra})

    tpc = pseudo_critical_temperature(fluid=fluid, pressure_Pa=pressure_Pa)
    if not tpc.usable or math.isnan(float(tpc)):
        return Result("UNKNOWN", "-", "CAN-COOL-REGIME",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes=f"the pseudo-critical temperature could not be "
                            f"located: {tpc.notes}")
    t_pc = float(tpc)
    ratio_pc = bulk_temperature_K / t_pc
    extra.update({"T_pc_K": t_pc, "T_bulk_over_T_pc": ratio_pc})
    notes.append(f"T_pc = {t_pc:.2f} K, T_bulk/T_pc = {ratio_pc:.3f}")

    bulk = coolant_state(fluid=fluid, pressure_Pa=pressure_Pa,
                         temperature_K=bulk_temperature_K)
    statuses.extend(bulk.statuses)
    extra.update({"rho_bulk_kg_m3": bulk.density_kg_m3,
                  "cp_bulk_J_kgK": bulk.specific_heat_J_kgK,
                  "Pr_bulk": bulk.prandtl})

    if 0.98 <= p_r <= 1.5:
        label = "NEAR_CRITICAL"
        statuses.append(Status.INVALID_MODEL_REGIME)
        notes.append(
            "0.98 <= p/p_crit <= 1.5. Haemisch et al. observed heat-transfer "
            "deterioration experimentally for p_in below ~70 bar, i.e. "
            "p/p_crit ~ 1.5, and none at p_in >= 100 bar (p/p_crit ~ 2.2).")
    elif bulk_temperature_K < t_pc:
        label = "LIQUID_LIKE"
    else:
        label = "GAS_LIKE"

    if wall_temperature_K is not None and wall_temperature_K > 0:
        wall = coolant_state(fluid=fluid, pressure_Pa=pressure_Pa,
                             temperature_K=wall_temperature_K)
        extra["rho_wall_kg_m3"] = wall.density_kg_m3
        if bulk.density_kg_m3 > 0 and wall.density_kg_m3 > 0:
            extra["rho_wall_over_rho_bulk"] = (wall.density_kg_m3
                                               / bulk.density_kg_m3)
        extra["T_wall_over_T_bulk"] = wall_temperature_K / bulk_temperature_K
        if bulk_temperature_K < t_pc < wall_temperature_K:
            label = "PSEUDO_CRITICAL"
            statuses.append(Status.INVALID_MODEL_REGIME)
            notes.append(
                "THE BOUNDARY LAYER SPANS THE PSEUDO-CRITICAL TEMPERATURE "
                "(T_bulk < T_pc < T_wall). Haemisch describes the mechanism: "
                "a low-density gas-like layer forms at the wall and acts as a "
                "thermal barrier, and the near-wall flow accelerates because "
                "of it. Every constant-property correlation in this module is "
                "outside its assumptions here.")

    return Result(
        label, "-", "CAN-COOL-REGIME", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"{SETZMANN_WAGNER_1991}; regime interpretation from "
               f"{HAEMISCH_2019}",
        source_locator="critical point and c_p locus from the reference EOS; "
                       "the HTD association is Haemisch's, OPENED",
        validity_domain="pure methane; MIXTURES NOT HANDLED (LNG composition "
                        "shifts the pseudo-critical locus and no source "
                        "consulted quantifies that shift)",
        uncertainty="the regime LABEL is categorical; the ratios inherit the "
                    "property uncertainties",
        notes="; ".join(notes), inputs={**inputs, **extra})


# ===========================================================================
# 2. HEAT-TRANSFER DETERIORATION ONSET
# ===========================================================================

#: Urbano & Nasuti (2013), as reproduced by Haemisch et al. (2019) Eq. (5).
#: q/G threshold in J/kg with the inlet pressure in Pa.
_UN_SLOPE_J_KG_PER_PA = 43.2e-6
_UN_INTERCEPT_J_KG = 31.4

#: Haemisch's empirical correction for a rectangular channel heated from ONE
#: side, Eq. (6). See the CONFLICT note in `htd_onset`.
_HAEMISCH_ONE_SIDED_FACTOR = 2.4

#: Haemisch's tested envelope. Outside it the criterion is an extrapolation.
_HAEMISCH_P_RANGE_PA = (58.0e5, 120.0e5)
_HAEMISCH_T_IN_RANGE_K = (135.0, 143.0)
_HAEMISCH_Q_MAX_W_M2 = 18.0e6
_HAEMISCH_DH_RANGE_M = (0.58e-3, 1.5e-3)
_HAEMISCH_AR_RANGE = (1.7, 30.0)


def htd_onset(*, heat_flux_W_m2: float, mass_flux_kg_m2s: float,
              inlet_pressure_Pa: float,
              channel_geometry: str = "CIRCULAR",
              channel_height_m: float | None = None,
              channel_width_m: float | None = None,
              hot_gas_width_per_channel_m: float | None = None,
              inlet_temperature_K: float | None = None,
              hydraulic_diameter_m: float | None = None) -> Result:
    r"""Is the channel above the heat-transfer-deterioration threshold?

    .. math:: \left(\frac{q_w}{G}\right)_{tr} = 43.2\times10^{-6}\,p_{in} + 31.4

    in J/kg with :math:`p_{in}` in Pa, for a **circular, uniformly heated
    tube**. HTD is predicted when the local :math:`q_w/G` EXCEEDS this.

    **THE ATTRIBUTION.** This is Urbano & Nasuti (2013), DOI 10.2514/1.T4001 --
    *not* a Pizzarelli paper, as this project's own task brief stated.
    Confirmed against Crossref, OpenAlex and Haemisch's reference list [8].

    **THE EVIDENCE POSITION IS SECONDHAND AND SAYS SO.** The Urbano & Nasuti
    paper is behind an AIAA 403 and has no open-access copy anywhere. The
    equation above was read in Haemisch et al. (2019), which IS open. So the
    criterion's own derivation range, its data-point count and its stated
    accuracy are all `NOT_REPORTED` to this project.

    **THE RECTANGULAR FORM: READING B IS NOW FALSIFIED.** Haemisch's Eq. (6)
    multiplies by ``(U/L_HG) * 2.4`` for a rectangular channel heated from one
    side, and two readings of the printed equation differ in where the
    parenthesis closes:

        reading A:  [43.2e-6 p + 31.4] * (U/L_HG) * 2.4
        reading B:   43.2e-6 p + 31.4  * (U/L_HG) * 2.4

    The printed bracket could not be read: the PDF text extractor returned
    three mutually inconsistent bracketings across three queries, and the
    thesis that would settle it typographically is a 43 MB file behind a
    30 MB fetch cap.

    **So it was settled by falsification against the source's own Table 6
    instead.** With Q1's geometry (U/L_HG = 2.069, x2.4 = 4.966):

        reading B flags HTD on ALL NINE tabulated test points -- including
        b3 at p_in = 102.64 bar and b4 at 120.84 bar, which the paper states
        in prose show NORMAL heat transfer ("for the higher pressure test
        cases, P_in >= 100 bar, no peak is observable"). A criterion with no
        discriminating power on the dataset it was fitted to is not the
        criterion.

        reading A reproduces b1 (HTD), b2 (HTD) and b4 (no HTD) correctly,
        and is indeterminate on b3 -- q/G = 2.48 kJ/kg against a threshold of
        2.36, a 5% margin well inside the source's own +/-11.4% on q/G.

    Reading A is therefore USED, the result carries `UNVALIDATED_ASSUMPTION`
    (not `PASS`) because the equation was never read typographically, and
    BOTH thresholds stay in `inputs`. See
    `docs/validation/coolant/HTD_CRITERION_VALIDATION.md`.

    **WHAT THIS CANNOT DO.** It is a binary onset flag. It does not give the
    deteriorated heat-transfer coefficient. The correlation that would --
    Pizzarelli, Numer. Heat Transfer A 69(3), 2016 -- is behind a Taylor &
    Francis 403 and **not one of its coefficients is available to this
    project.**
    """
    inputs = {"q_W_m2": heat_flux_W_m2, "G_kg_m2s": mass_flux_kg_m2s,
              "p_in_Pa": inlet_pressure_Pa, "geometry": channel_geometry}
    if heat_flux_W_m2 < 0:
        return Result(False, "-", "CAN-COOL-HTD",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="heat flux into the wall cannot be negative here")
    if mass_flux_kg_m2s <= 0:
        return Result(False, "-", "CAN-COOL-HTD",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="mass flux must be > 0. A zero-flow channel has "
                            "no forced-convection heat transfer at all, which "
                            "is a different and worse problem than HTD.")
    if inlet_pressure_Pa <= 0:
        return Result(False, "-", "CAN-COOL-HTD",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="inlet pressure must be > 0")

    q_over_g = heat_flux_W_m2 / mass_flux_kg_m2s
    base = _UN_SLOPE_J_KG_PER_PA * inlet_pressure_Pa + _UN_INTERCEPT_J_KG

    statuses: list[Status] = [Status.PASS]
    notes: list[str] = []
    assumptions: list[Assumption] = []
    extra: dict[str, Any] = {"q_over_G_J_kg": q_over_g,
                             "threshold_circular_J_kg": base}

    geom = channel_geometry.strip().upper()
    if geom == "CIRCULAR":
        threshold = base
    elif geom == "RECTANGULAR_ONE_SIDED":
        if (channel_height_m is None or channel_width_m is None
                or hot_gas_width_per_channel_m is None):
            return Result(
                False, "-", "CAN-COOL-HTD",
                status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                evidence_level=EvidenceLevel.E0,
                source=HAEMISCH_2019, source_locator="Eq. (6)",
                notes="the rectangular form needs the channel height and "
                      "width (for U = 2h + b) and the hot-gas-side width per "
                      "channel L_HG. Without them the geometric factor cannot "
                      "be formed, and the circular threshold is NOT "
                      "substituted -- Haemisch's factor U/L_HG is ~2.6 for "
                      "his own hardware, so the difference is a factor of "
                      "six, not a rounding.")
        if (channel_height_m <= 0 or channel_width_m <= 0
                or hot_gas_width_per_channel_m <= 0):
            return Result(False, "-", "CAN-COOL-HTD",
                          status=Status.PHYSICALLY_INVALID, inputs=inputs,
                          notes="channel dimensions must be > 0")
        u = 2.0 * channel_height_m + channel_width_m
        geo = u / hot_gas_width_per_channel_m
        thr_a = base * geo * _HAEMISCH_ONE_SIDED_FACTOR
        thr_b = (_UN_SLOPE_J_KG_PER_PA * inlet_pressure_Pa
                 + _UN_INTERCEPT_J_KG * geo * _HAEMISCH_ONE_SIDED_FACTOR)
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
        assumptions.append(Assumption(
            "haemisch_eq6_parenthesization", "reading A", "-",
            "reading B is FALSIFIED by the source's own tabulated data -- it "
            "predicts HTD for all nine Table 6 test points, including the "
            "two at p_in >= 100 bar that the paper states show normal heat "
            "transfer. Reading A is USED but was never read typographically.",
            Status.UNVALIDATED_ASSUMPTION))
        notes.append(
            f"Haemisch Eq. (6) parenthesization: reading A = {thr_a:.1f} J/kg, "
            f"reading B = {thr_b:.1f} J/kg, a factor of "
            f"{max(thr_a, thr_b) / max(min(thr_a, thr_b), 1e-9):.1f}. READING "
            f"B IS FALSIFIED against Table 6 of the source (it flags HTD on "
            f"every one of the nine tabulated cases, including b3 at 102.6 "
            f"bar and b4 at 120.8 bar, which the paper states are NORMAL "
            f"heat transfer). Reading A reproduces 3 of the 4 prose-stated "
            f"cases and is indeterminate on the fourth within the source's "
            f"own +/-11.4% on q/G. The printed bracket was never read -- both "
            f"values remain in `inputs`.")
        threshold = thr_a
        extra.update({"threshold_reading_A_J_kg": thr_a,
                      "threshold_reading_B_J_kg": thr_b,
                      "U_over_L_HG": geo})
    else:
        return Result(
            False, "-", "CAN-COOL-HTD",
            status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
            notes=f"channel_geometry={channel_geometry!r} is not one the "
                  f"sources cover. Urbano & Nasuti derived the threshold for "
                  f"a CIRCULAR uniformly heated tube; Haemisch corrected it "
                  f"for a RECTANGULAR channel heated from ONE side. No other "
                  f"geometry has a source.")

    deteriorated = q_over_g > threshold
    extra["threshold_used_J_kg"] = threshold

    # --- Haemisch's tested envelope -------------------------------------
    plo, phi = _HAEMISCH_P_RANGE_PA
    if not (plo <= inlet_pressure_Pa <= phi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"p_in = {inlet_pressure_Pa / 1e5:.1f} bar is outside "
                     f"Haemisch's tested {plo / 1e5:.0f}-{phi / 1e5:.0f} bar")
    if heat_flux_W_m2 > _HAEMISCH_Q_MAX_W_M2:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"q = {heat_flux_W_m2 / 1e6:.1f} MW/m^2 exceeds "
                     f"Haemisch's tested maximum of "
                     f"{_HAEMISCH_Q_MAX_W_M2 / 1e6:.0f} MW/m^2")
    if inlet_temperature_K is not None:
        tlo, thi = _HAEMISCH_T_IN_RANGE_K
        if not (tlo <= inlet_temperature_K <= thi):
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(f"T_in = {inlet_temperature_K:.1f} K is outside "
                         f"Haemisch's tested {tlo:.0f}-{thi:.0f} K")
    if hydraulic_diameter_m is not None:
        dlo, dhi = _HAEMISCH_DH_RANGE_M
        if not (dlo <= hydraulic_diameter_m <= dhi):
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(f"d_H = {hydraulic_diameter_m * 1e3:.2f} mm is "
                         f"outside Haemisch's tested "
                         f"{dlo * 1e3:.2f}-{dhi * 1e3:.2f} mm")
    if (channel_height_m is not None and channel_width_m is not None
            and channel_width_m > 0):
        ar = channel_height_m / channel_width_m
        extra["aspect_ratio"] = ar
        if ar > 20.0:
            statuses.append(Status.OUT_OF_CORRELATION_RANGE)
            notes.append(
                f"aspect ratio {ar:.1f}: Haemisch states his criterion "
                f"predicts >98% of test points correctly WITH THE EXCEPTION "
                f"of high-aspect-ratio channels, where thermal stratification "
                f"shrinks the susceptible area and his AR=30 channel showed "
                f"no HTD at all.")

    if deteriorated:
        notes.insert(0, f"q/G = {q_over_g:.1f} J/kg EXCEEDS the threshold "
                        f"{threshold:.1f} J/kg -> HTD_DETECTED")
    else:
        notes.insert(0, f"q/G = {q_over_g:.1f} J/kg is below the threshold "
                        f"{threshold:.1f} J/kg -> "
                        f"NORMAL_SUPERCRITICAL_HEAT_TRANSFER "
                        f"(margin {threshold - q_over_g:.1f} J/kg)")

    return Result(
        deteriorated, "-", "CAN-COOL-HTD", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"{URBANO_NASUTI_2013} (threshold; SECONDHAND -- read in "
               f"Haemisch, primary is AIAA 403 with no OA copy); "
               f"{HAEMISCH_2019} (rectangular correction and experimental "
               f"envelope; OPENED at DLR elib)",
        source_locator="Haemisch Eq. (5) reproducing Urbano & Nasuti; "
                       "Haemisch Eq. (6) for the rectangular form",
        validity_domain=f"Haemisch's tested envelope: p_in "
                        f"{plo / 1e5:.0f}-{phi / 1e5:.0f} bar, T_in "
                        f"{_HAEMISCH_T_IN_RANGE_K[0]:.0f}-"
                        f"{_HAEMISCH_T_IN_RANGE_K[1]:.0f} K, q <= "
                        f"{_HAEMISCH_Q_MAX_W_M2 / 1e6:.0f} MW/m^2, d_H "
                        f"{_HAEMISCH_DH_RANGE_M[0] * 1e3:.2f}-"
                        f"{_HAEMISCH_DH_RANGE_M[1] * 1e3:.2f} mm, copper, "
                        f"one-sided heating, k_s ~ 0.2 um (SMOOTH)",
        uncertainty="Haemisch: '>98% of all test points correct', EXCEPT "
                    "high-aspect-ratio channels. The Urbano & Nasuti "
                    "threshold's own accuracy is NOT_REPORTED to this "
                    "project. This is a BINARY FLAG -- it does not quantify "
                    "the deteriorated h_c, and no accessible source does.",
        assumptions=tuple(assumptions),
        notes="; ".join(notes), inputs={**inputs, **extra})


# ===========================================================================
# 3. ROUGHNESS -- the term that actually dominates
# ===========================================================================

class RoughnessRegime(str, Enum):
    HYDRAULICALLY_SMOOTH = "HYDRAULICALLY_SMOOTH"
    TRANSITIONAL = "TRANSITIONAL"
    FULLY_ROUGH = "FULLY_ROUGH"


#: Nikuradse's own boundaries are stated in log10(v* k / nu): 0.55 and 1.83.
#: 10^0.55 = 3.548 and 10^1.83 = 67.61. The "3.5" and "68" universally quoted
#: are those two numbers rounded -- they are not independent conventions.
_NIKURADSE_SMOOTH_LIMIT = 10.0 ** 0.55
_NIKURADSE_ROUGH_LIMIT = 10.0 ** 1.83


def roughness_regime(*, reynolds: float, relative_roughness: float,
                     friction_factor: float) -> Result:
    r"""Classify by the roughness Reynolds number.

    .. math:: h_s^+ = Re\,\frac{\varepsilon}{D}\sqrt{f/8}

    Nikuradse's thresholds, traced to their source: he states the boundaries
    as :math:`\log_{10}(v^* k/\nu) = 0.55` and :math:`1.83`, which are
    :math:`h_s^+ = 3.548` and :math:`67.61` -- the familiar "3.5" and "68"
    are those numbers rounded, not separate conventions. Dipprey & Sabersky's
    ":math:`\varepsilon^* > 67`" is the same boundary rounded down.

    **BUT 68 IS A SAND-GRAIN NUMBER.** Flack & Schultz (2010) document that
    fully-rough onset is strongly surface-dependent and occurs as low as
    :math:`k_s^+ \approx 18`-25 for non-sand-grain surfaces (honed pipe,
    scratched plate, commercial steel). An additively manufactured channel is
    emphatically not uniform sand. So a surface that clears 68 is fully rough
    by every published criterion, and one that does NOT clear it may still be
    fully rough by the surface-specific ones. The result says so.
    """
    inputs = {"Re": reynolds, "eps_over_D": relative_roughness,
              "f": friction_factor}
    if reynolds <= 0 or friction_factor <= 0:
        return Result("INVALID", "-", "CAN-COOL-ROUGHREGIME",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and f must be > 0")
    if relative_roughness < 0:
        return Result("INVALID", "-", "CAN-COOL-ROUGHREGIME",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="relative roughness cannot be negative")

    h_plus = reynolds * relative_roughness * math.sqrt(friction_factor / 8.0)
    statuses = [Status.PASS]
    notes = [f"h_s+ = {h_plus:.2f}"]

    if h_plus < _NIKURADSE_SMOOTH_LIMIT:
        label = RoughnessRegime.HYDRAULICALLY_SMOOTH.value
    elif h_plus < _NIKURADSE_ROUGH_LIMIT:
        label = RoughnessRegime.TRANSITIONAL.value
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(
            "TRANSITIONAL by Nikuradse's sand-grain boundaries. Neither a "
            "smooth-tube correlation nor the fully-rough Dipprey-Sabersky "
            "form is inside its own domain here. Flack & Schultz document "
            "fully-rough onset as low as k_s+ ~ 18-25 for non-sand surfaces, "
            "so an AM channel in this band may already be fully rough -- but "
            "no source consulted gives an AM-specific threshold.")
    else:
        label = RoughnessRegime.FULLY_ROUGH.value
        notes.append("fully rough by Nikuradse's own boundary, which is the "
                     "most conservative of the published thresholds")

    return Result(
        label, "-", "CAN-COOL-ROUGHREGIME", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"{NIKURADSE_1933}; surface-dependence from "
               f"{FLACK_SCHULTZ_2010}",
        source_locator="Nikuradse Eq. group (3.21), boundaries at "
                       "log10(v* k/nu) = 0.55 and 1.83",
        validity_domain="UNIFORM SAND-GRAIN roughness. An LPBF surface is "
                        "adhered partially-sintered powder with positive "
                        "skewness -- not sand.",
        uncertainty="the 3.5/68 boundaries are exact conversions of "
                    "Nikuradse's own log10 values; their APPLICABILITY to a "
                    "non-sand surface is the open question",
        notes="; ".join(notes),
        inputs={**inputs, "h_s_plus": h_plus,
                "nikuradse_smooth_limit": _NIKURADSE_SMOOTH_LIMIT,
                "nikuradse_rough_limit": _NIKURADSE_ROUGH_LIMIT})


def sand_grain_from_Ra(*, Ra_m: float, hydraulic_diameter_m: float | None = None,
                       rms_roughness_m: float | None = None,
                       skewness: float | None = None) -> Result:
    r"""Convert a metrology Ra into an equivalent sand-grain roughness.

    **THIS RETURNS `CONFLICT_UNRESOLVED` AND IT ALWAYS WILL, UNTIL SOMEBODY
    MEASURES A PRESSURE DROP.**

    Three open primary sources give structurally incompatible answers:

        Adams, Grant & Watson (2012)   eps = 5.863 Ra           (a constant)
        Stimpson et al. (2017)         k_s/D_h = 18 Ra/D_h - 0.05
        Flack & Schultz (2010)         k_s = 4.43 k_rms (1+s_k)^1.37

    They are not scattered around a common value -- they are different
    FUNCTIONS. Adams' ratio k_s/Ra is a constant 5.863; Stimpson's varies from
    ~0.4 to ~16 across his own validity range; Flack & Schultz's does not
    depend on Ra at all, but on the RMS height and the SKEWNESS of the height
    distribution. The last is the physically better-motivated form, and it is
    precisely why no pure Ra -> k_s factor can be universal: two surfaces with
    identical Ra and different skewness have different k_s. An LPBF surface,
    which is adhered partially-sintered powder, is the worst case for this.

    **None of the three states an uncertainty for its own conversion.** Adams'
    own demonstration case under-predicted by 33%.

    Consequence for Kryptonis: converting a metrology Ra into a Colebrook
    epsilon introduces an uncertainty **larger than every other uncertainty in
    the coolant chain**. The only reliable route in the literature -- and the
    one Latini/Fiore/Nasuti use -- is to back epsilon out of a MEASURED
    pressure drop, which is not available at design time.
    """
    inputs = {"Ra_m": Ra_m, "D_h_m": hydraulic_diameter_m,
              "k_rms_m": rms_roughness_m, "skewness": skewness}
    if Ra_m <= 0:
        return Result(float("nan"), "m", "CAN-COOL-RA2EPS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Ra must be > 0")

    candidates: dict[str, float] = {}
    notes: list[str] = []

    adams = 5.863 * Ra_m
    candidates["Adams_2012"] = adams
    notes.append(f"Adams (eps = 5.863 Ra): {adams * 1e6:.2f} um")

    if hydraulic_diameter_m is not None and hydraulic_diameter_m > 0:
        ra_over_d = Ra_m / hydraulic_diameter_m
        stimpson = (18.0 * ra_over_d - 0.05) * hydraulic_diameter_m
        if ra_over_d <= 0.0028:
            notes.append(
                f"Stimpson: Ra/D_h = {ra_over_d:.5f} is at or below his "
                f"stated 0.0028 floor, where his affine form returns a "
                f"NEGATIVE k_s. NOT EVALUATED.")
        else:
            candidates["Stimpson_2017"] = stimpson
            notes.append(f"Stimpson (k_s/D_h = 18 Ra/D_h - 0.05): "
                         f"{stimpson * 1e6:.2f} um "
                         f"(k_s/Ra = {stimpson / Ra_m:.2f})")
    else:
        notes.append("Stimpson NOT EVALUATED: needs the hydraulic diameter -- "
                     "his conversion is a function of Ra/D_h, not of Ra")

    if rms_roughness_m is not None and skewness is not None:
        if rms_roughness_m > 0 and (1.0 + skewness) > 0:
            fs = 4.43 * rms_roughness_m * (1.0 + skewness) ** 1.37
            candidates["Flack_Schultz_2010"] = fs
            notes.append(f"Flack & Schultz (k_s = 4.43 k_rms (1+s_k)^1.37): "
                         f"{fs * 1e6:.2f} um")
    else:
        notes.append("Flack & Schultz NOT EVALUATED: needs the RMS height AND "
                     "the skewness. Ra alone cannot reach it, which is their "
                     "point.")

    lo = min(candidates.values())
    hi = max(candidates.values())
    spread = hi / lo if lo > 0 else float("inf")

    return Result(
        float("nan"), "m", "CAN-COOL-RA2EPS",
        status=Status.CONFLICT_UNRESOLVED,
        all_statuses=frozenset({Status.CONFLICT_UNRESOLVED}),
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE, validation=Validation.NONE,
        source=f"{ADAMS_2012}; {STIMPSON_2017}; {FLACK_SCHULTZ_2010}",
        source_locator="Adams Table 1; Stimpson Eq. (9); Flack & Schultz "
                       "Eq. (27) -- all three OPENED",
        validity_domain=UNKNOWN_DOMAIN + " -- the three sources do not share "
                                         "a domain or even a functional form",
        uncertainty=f"NO source states an uncertainty for its own "
                    f"conversion. Candidates here span {spread:.1f}x. Adams' "
                    f"own demonstration case under-predicted by 33%.",
        notes="NO VALUE IS RETURNED. " + " | ".join(notes),
        inputs={**inputs, "candidates_m": candidates,
                "spread_factor": spread})


# ===========================================================================
# 4. FRICTION FACTORS -- the Gnielinski question, resolved
# ===========================================================================

class GnielinskiVintage(str, Enum):
    """Which Gnielinski publication you are citing. They differ.

    This enum exists because "the Gnielinski correlation" is two things. The
    friction factor is not a free choice: it follows from which paper you name.
    """

    #: Forsch. Ing.-Wes. 41(1) 1975 / Int. Chem. Eng. 16(2) 1976 -> Filonenko
    ICE_1976 = "ICE_1976"
    #: VDI Heat Atlas 2nd English ed. 2010, Ch. G1 -> Konakov
    VDI_2010 = "VDI_2010"


def friction_factor_filonenko(*, reynolds: float) -> Result:
    r""":math:`f = (1.82\log_{10}Re - 1.64)^{-2}`. Smooth tube.

    Filonenko's own paper is not digitised anywhere. The equation is
    reproduced verbatim in **Petukhov (1970), Eq. (49)**, which WAS opened --
    a scanned reproduction of the Advances in Heat Transfer chapter -- and the
    Reynolds range below is Petukhov's statement for that equation.

    Note: in the accessible portion of Petukhov's text, Eq. (49) is attributed
    to **Petukhov & Popov**, not to Filonenko. The common claim "Petukhov used
    Filonenko's formula" is therefore `NOT_CONFIRMED` by the primary text.
    """
    inputs = {"Re": reynolds}
    if reynolds <= 0:
        return Result(float("nan"), "-", "CAN-COOL-F-FILONENKO",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Reynolds number must be > 0")
    denom = 1.82 * math.log10(reynolds) - 1.64
    if denom <= 0:
        return Result(float("nan"), "-", "CAN-COOL-F-FILONENKO",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the bracket is non-positive at this Re")
    statuses = [Status.PASS]
    notes: list[str] = []
    if not (1.0e4 <= reynolds <= 5.0e6):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside Petukhov's stated "
                     f"1e4-5e6 for this equation")
    return Result(
        denom ** -2, "-", "CAN-COOL-F-FILONENKO", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"{FILONENKO_1954} (NOT DIGITISED; equation read in "
               f"{PETUKHOV_1970}, Eq. (49), which was OPENED)",
        source_locator="Petukhov (1970) Eq. (49); Filonenko's own printed "
                       "equation is NOT_REPORTED",
        validity_domain="SMOOTH TUBE; Re 1e4-5e6 per Petukhov. The equation "
                        "contains no roughness parameter and cannot "
                        "represent a rough pipe.",
        uncertainty=NOT_REPORTED,
        notes="; ".join(notes), inputs=inputs)


def friction_factor_konakov(*, reynolds: float) -> Result:
    r""":math:`f = (1.8\log_{10}Re - 1.5)^{-2}`. Smooth tube, by its title.

    Konakov's 1946 paper is not digitised. Its TITLE, as printed in
    Gnielinski's own VDI Heat Atlas G1 reference list [17], is *"A new
    correlation for the friction coefficient in smooth tubes"* -- so the
    smooth-tube restriction is established by the source's own title even
    though the paper could not be read.

    **Konakov's own stated Reynolds range is `NOT_REPORTED`.** No accessible
    source gives one. No range check is therefore applied here, because
    inventing one would be worse than having none.
    """
    inputs = {"Re": reynolds}
    if reynolds <= 0:
        return Result(float("nan"), "-", "CAN-COOL-F-KONAKOV",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Reynolds number must be > 0")
    denom = 1.80 * math.log10(reynolds) - 1.50
    if denom <= 0:
        return Result(float("nan"), "-", "CAN-COOL-F-KONAKOV",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the bracket is non-positive at this Re")
    return Result(
        denom ** -2, "-", "CAN-COOL-F-KONAKOV",
        status=Status.UNVALIDATED_ASSUMPTION,
        all_statuses=frozenset({Status.UNVALIDATED_ASSUMPTION}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"{KONAKOV_1946} (NOT DIGITISED; equation read SECONDHAND, "
               f"cited to the VDI-Waermeatlas)",
        source_locator="cited in Gnielinski's own VDI Heat Atlas Ch. G1 "
                       "reference list [17]; Konakov's own printed equation "
                       "is NOT_REPORTED",
        validity_domain="SMOOTH TUBES (established by the source's own "
                        "title). Reynolds range NOT_REPORTED -- no range "
                        "check is applied rather than inventing one.",
        uncertainty=NOT_REPORTED,
        notes="UNVALIDATED_ASSUMPTION because no stated validity range exists "
              "to check against, not because the equation is doubted",
        inputs=inputs)


def gnielinski_friction_factor(*, reynolds: float,
                               vintage: GnielinskiVintage) -> Result:
    """THE resolution of the friction-factor question. It was never a conflict.

    Gnielinski used **Filonenko** in the 1975 German original (and therefore
    in its 1976 English translation) and **Konakov** in the 2010 VDI Heat
    Atlas chapter he wrote himself. The evidence is his own reference lists:

        Forsch. Ing.-Wes. 41(1) 1975, Schrifttum [34]  -> Filonenko,
            Teploenergetika 1(4) 1954, 40/44.
            **Konakov does not appear in that reference list.**

        VDI Heat Atlas 2nd Eng. ed. 2010, Ch. G1, ref [17] -> Konakov,
            Berichte der Akademie der Wissenschaften der UdSSR 51:503-506.
            **Filonenko does not appear in that reference list.**

    So the two camps of secondary sources were each quoting Gnielinski
    correctly, from different papers. There is nothing to resolve except the
    LABEL. `vintage` is required for exactly that reason: naming the paper
    determines the friction factor, and this function will not let a caller
    have the correlation without naming the paper.

    His 2013 IJHMT restatement cites BOTH -- which of them it adopts is
    `NOT_OPENED` (ScienceDirect robots.txt), and there is a Corrigendum
    (IJHMT 81, 2015, 638) that should be read first.

    The practical stake is small: ~1.1% in f at Re = 1e5 and ~0.4% at
    Re = 1e6, sub-1% in Nu. **The smooth-versus-rough question is worth tens
    of percent and is the one that matters.**
    """
    if vintage is GnielinskiVintage.ICE_1976:
        r = friction_factor_filonenko(reynolds=reynolds)
    elif vintage is GnielinskiVintage.VDI_2010:
        r = friction_factor_konakov(reynolds=reynolds)
    else:  # pragma: no cover - enum is exhaustive
        raise ValueError(f"unknown Gnielinski vintage {vintage!r}")

    other = (friction_factor_konakov(reynolds=reynolds)
             if vintage is GnielinskiVintage.ICE_1976
             else friction_factor_filonenko(reynolds=reynolds))
    spread = (abs(float(other) / float(r) - 1.0) * 100.0
              if r.usable and float(r) > 0 and other.usable else float("nan"))

    return Result(
        r.value, "-", "CAN-COOL-F-GNIELINSKI", status=r.status,
        all_statuses=r.all_statuses, evidence_level=r.evidence_level,
        verification=r.verification, validation=r.validation,
        source=r.source, source_locator=r.source_locator,
        validity_domain=r.validity_domain, uncertainty=r.uncertainty,
        notes=f"vintage {vintage.value}: {r.equation_id.split('-')[-1]}. "
              f"The other vintage differs by {spread:.2f}% here. "
              f"{r.notes}".strip(),
        inputs={**r.inputs, "vintage": vintage.value,
                "other_vintage_f": other.value,
                "spread_percent": spread})


# ===========================================================================
# 5. NUSSELT CORRELATIONS NOT ALREADY IN canonical.thermal
# ===========================================================================

def nusselt_petukhov_1970(*, reynolds: float, prandtl: float,
                          friction_factor: float,
                          variable_property_corrections: bool = False) -> Result:
    r"""Petukhov (1970) Eq. (48). **The leading constant is 1.07, not 1.**

    .. math::
        Nu = \frac{(\xi/8)\,Re\,Pr}
                  {1.07 + 12.7\sqrt{\xi/8}\,(Pr^{2/3}-1)}

    This matters because Gnielinski's better-known form uses ``1 +`` and
    ``(Re - 1000)`` instead: **those are Gnielinski's modifications, not
    Petukhov's**, and code that calls this "Petukhov" while using 1.0 is
    citing the wrong paper for the wrong equation.

    Petukhov's own stated envelope, from the opened scan: **Re 1e4-5e6,
    Pr 0.5-2000**, accuracy *"within 1 %"* except in 5e5 < Re < 5e6.

    The variable-property corrections ``K_1(xi) = 1 + 3.45 xi`` and
    ``K_2(Pr) = 11.7 + 1.8 Pr^(-1/3)`` appear in the same equation. They are
    NOT applied by default because the accessible text does not make clear
    how they are normalised, and applying a correction whose normalisation is
    guessed is worse than not applying it.
    """
    inputs = {"Re": reynolds, "Pr": prandtl, "f": friction_factor}
    if reynolds <= 0 or prandtl <= 0 or friction_factor <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-PETUKHOV",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re, Pr and f must be > 0")
    f8 = friction_factor / 8.0
    denom = 1.07 + 12.7 * math.sqrt(f8) * (prandtl ** (2.0 / 3.0) - 1.0)
    if denom <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-PETUKHOV",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="denominator non-positive at this Pr and f")
    nu = f8 * reynolds * prandtl / denom

    statuses = [Status.PASS]
    notes: list[str] = []
    if variable_property_corrections:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append(
            "the K_1/K_2 variable-property corrections are NOT applied: the "
            "accessible scan does not establish their normalisation, and a "
            "guessed normalisation is not a correction.")
    if not (1.0e4 <= reynolds <= 5.0e6):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside Petukhov's stated 1e4-5e6")
    if not (0.5 <= prandtl <= 2000.0):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Pr = {prandtl:.3f} outside Petukhov's stated 0.5-2000")

    return Result(
        nu, "-", "CAN-COOL-NU-PETUKHOV", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=PETUKHOV_1970,
        source_locator="Eq. (48) with the friction factor at Eq. (49); "
                       "OPENED as a scanned reproduction",
        validity_domain="Re 1e4-5e6, Pr 0.5-2000, SMOOTH TUBE",
        uncertainty="'within 1%' per Petukhov, except in 5e5 < Re < 5e6",
        notes="; ".join(notes),
        inputs={**inputs, "leading_constant": 1.07})


#: Dipprey & Sabersky's fully-rough heat-transfer function, JPL TR 32-269
#: Eq. (23), with the Nikuradse fully-rough constant A from p.28.
_DS_KF = 5.19
_DS_EPS_EXPONENT = 0.20
_DS_PR_EXPONENT = 0.44
_DS_A = 8.48
_DS_FULLY_ROUGH_MIN = 67.0
_DS_PR_RANGE = (1.20, 5.94)
_DS_RE_RANGE = (1.4e4, 5.0e5)
_DS_EPS_OVER_D_RANGE = (0.0024, 0.0488)


def nusselt_dipprey_sabersky(*, reynolds: float, prandtl: float,
                             relative_roughness: float,
                             friction_factor: float) -> Result:
    r"""Fully-rough heat transfer. THE correlation LRE cooling groups use.

    .. math::
        St = \frac{f/8}
             {1 + \sqrt{f/8}\,\bigl[5.19\,(\varepsilon^+)^{0.20}Pr^{0.44}
              - 8.48\bigr]},
        \qquad \varepsilon^+ = Re\,\frac{\varepsilon}{D}\sqrt{f/8}

    and :math:`Nu = St\,Re\,Pr`.

    Primary source OPENED: the IJHMT 1963 paper is robots-blocked, but its
    JPL precursor (TR 32-269, 6 June 1962) is fully open and is itself a
    primary document. Coefficients read there: :math:`k_f = 5.19` (p.21 and
    p.28), exponents 0.20 and 0.44 (Eq. 23), :math:`A = 8.48` (p.28, *"the
    value of A was set at 8.48"*, Nikuradse's fully-rough constant).

    **IT WAS FITTED ON DISTILLED WATER.** Prandtl number was varied by
    varying the bulk water temperature over 1.20-5.94, not by changing fluid.
    It has never been fitted to hydrogen, to methane, or to any supercritical
    fluid. Supercritical methane at 9 MPa has Pr ~ 1.2-1.8, which is inside
    the fitted band by coincidence of number, not by similarity of fluid.

    **AND ROCKET PRACTICE EXTRAPOLATES IT 4x.** Dipprey & Sabersky tested
    :math:`\varepsilon/D \le 0.0488`; Latini/Fiore/Nasuti apply it at 0.21,
    and state no heat-transfer error figure for doing so.

    The stated :math:`\pm 5\%` in the source is EXPERIMENTAL uncertainty on
    the measured coefficients. **A goodness-of-fit for Eq. (23) itself is
    NOT_REPORTED**, so +/-5% must not be quoted as the correlation's accuracy.
    """
    inputs = {"Re": reynolds, "Pr": prandtl,
              "eps_over_D": relative_roughness, "f": friction_factor}
    if reynolds <= 0 or prandtl <= 0 or friction_factor <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-DIPPREY",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re, Pr and f must be > 0")
    if relative_roughness <= 0:
        return Result(
            float("nan"), "-", "CAN-COOL-NU-DIPPREY",
            status=Status.INVALID_MODEL_REGIME, inputs=inputs,
            notes="this is a ROUGH-wall correlation; a smooth wall has no "
                  "roughness Reynolds number and the fully-rough form is "
                  "undefined there. Use a smooth-tube correlation instead.")

    f8 = friction_factor / 8.0
    eps_plus = reynolds * relative_roughness * math.sqrt(f8)
    g = _DS_KF * eps_plus ** _DS_EPS_EXPONENT * prandtl ** _DS_PR_EXPONENT
    denom = 1.0 + math.sqrt(f8) * (g - _DS_A)
    if denom <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-DIPPREY",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes=f"denominator non-positive (g = {g:.3f}, "
                            f"A = {_DS_A}); the roughness function has fallen "
                            f"below the additive constant, which is outside "
                            f"the similarity law's meaningful range")
    st = f8 / denom
    nu = st * reynolds * prandtl

    statuses = [Status.PASS]
    notes = [f"eps+ = {eps_plus:.1f}, St = {st:.6f}"]
    if eps_plus < _DS_FULLY_ROUGH_MIN:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"eps+ = {eps_plus:.1f} is below the source's own "
                     f"fully-rough lower limit of {_DS_FULLY_ROUGH_MIN:.0f} "
                     f"(Fig. 18 label). This is the TRANSITIONAL regime and "
                     f"Eq. (23) is not the applicable branch.")
    plo, phi = _DS_PR_RANGE
    if not (plo <= prandtl <= phi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Pr = {prandtl:.3f} outside the fitted {plo}-{phi} "
                     f"(DISTILLED WATER only)")
    rlo, rhi = _DS_RE_RANGE
    if not (rlo <= reynolds <= rhi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside the tested "
                     f"{rlo:.1e}-{rhi:.1e}")
    elo, ehi = _DS_EPS_OVER_D_RANGE
    if not (elo <= relative_roughness <= ehi):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"eps/D = {relative_roughness:.4f} outside the tested "
                     f"{elo}-{ehi}")

    return Result(
        nu, "-", "CAN-COOL-NU-DIPPREY", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=DIPPREY_SABERSKY_1962,
        source_locator="Eq. (23) with k_f = 5.19, p = 0.20, m = 0.44; "
                       "A = 8.48 from p.28; JPL TR 32-269 OPENED in full",
        validity_domain=f"FULLY ROUGH (eps+ > {_DS_FULLY_ROUGH_MIN:.0f}); "
                        f"Pr {plo}-{phi}; Re {rlo:.1e}-{rhi:.1e}; "
                        f"eps/D {elo}-{ehi}; DISTILLED WATER, close-packed "
                        f"granular electroplated-nickel roughness",
        uncertainty="the source's +/-5% is EXPERIMENTAL uncertainty on the "
                    "measured coefficients. A goodness-of-fit for Eq. (23) "
                    "itself is NOT_REPORTED and +/-5% must not be quoted as "
                    "the correlation's accuracy.",
        notes="; ".join(notes),
        inputs={**inputs, "eps_plus": eps_plus, "St": st,
                "roughness_function_g": g})


def nusselt_stimpson_am(*, reynolds: float, prandtl: float,
                        friction_factor: float,
                        ks_over_Dh: float | None = None,
                        friction_factor_measured: bool = False) -> Result:
    r"""Stimpson et al. (2017), fitted on ADDITIVELY MANUFACTURED channels.

    .. math:: Nu = Re^{0.5}\,\frac{29\,Pr\,(f/8)^{0.6}}{(1+Pr)^{2/3}}

    Included because these authors **abandoned the classical route** for AM
    channels -- they state Gnielinski and the smooth-tube family do not
    describe their data -- and fitted this instead on DMLS CoCr and Inconel
    718 channels. It is the only correlation in this project fitted on the
    actual manufacturing process a Kryptonis chamber would use.

    Stated: *"much more capable of predicting flow for relative roughness
    between 0.07 < k_s/D_h < 0.5"*, with **+/-15% when f is measured directly
    and +/-30% when f is predicted.** Their envelope: D_h 406-1275 um,
    Ra 9.27-13.8 um, Re 3 000-30 000, 45-degree build angle, 30 um powder,
    40 um layer.

    **Re 3 000-30 000 is the problem.** A rocket regenerative channel runs at
    Re ~ 1e5-1e6, one to two orders above. Setting `friction_factor_measured`
    tells the result which accuracy band applies; it does not change the
    number.
    """
    inputs = {"Re": reynolds, "Pr": prandtl, "f": friction_factor,
              "ks_over_Dh": ks_over_Dh,
              "f_measured": friction_factor_measured}
    if reynolds <= 0 or prandtl <= 0 or friction_factor <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-STIMPSON",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re, Pr and f must be > 0")
    nu = (reynolds ** 0.5 * 29.0 * prandtl * (friction_factor / 8.0) ** 0.6
          / (1.0 + prandtl) ** (2.0 / 3.0))

    statuses = [Status.PASS]
    notes: list[str] = []
    if not (3.0e3 <= reynolds <= 3.0e4):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"Re = {reynolds:.3e} outside the fitted 3e3-3e4. A "
                     f"rocket regenerative channel runs one to two orders "
                     f"above this and the extrapolation is unbounded.")
    if ks_over_Dh is not None and not (0.07 <= ks_over_Dh <= 0.5):
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"k_s/D_h = {ks_over_Dh:.4f} outside the stated "
                     f"0.07-0.5")
    if ks_over_Dh is None:
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
        notes.append("k_s/D_h not supplied, so the correlation's stated "
                     "roughness range could not be checked")

    return Result(
        nu, "-", "CAN-COOL-NU-STIMPSON", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=STIMPSON_2017,
        source_locator="Eq. (12); roughness conversion at Eq. (9); friction "
                       "closed with Colebrook at Eq. (10)",
        validity_domain="AM (DMLS) channels, 0.07 < k_s/D_h < 0.5, "
                        "Re 3e3-3e4, D_h 406-1275 um, CoCr and Inconel 718, "
                        "AIR",
        uncertainty="+/-15% with f measured directly; +/-30% with f "
                    "predicted. Roughness measurement itself +/-4 um from CT, "
                    "which the authors call semiquantitative."
                    if friction_factor_measured else
                    "+/-30% (f predicted rather than measured). +/-15% is "
                    "available only when f comes from a direct measurement.",
        notes="; ".join(notes), inputs=inputs)


def nusselt_schacht_quentmeyer(*, reynolds: float, prandtl: float,
                               properties_are_integral_averaged: bool = False
                               ) -> Result:
    r"""NASA TN D-7207. ``St* Pr*^0.6 = 0.023 Re*^-0.2``, i.e.
    :math:`Nu = 0.023\,Re^{0.8}Pr^{0.4}` **with integral-averaged properties.**

    The whole content of this correlation is the property averaging, Eq. (8):

    .. math::
        \bar X = \frac{1}{T_{c,w}-T_{c,s}}\int_{T_{c,s}}^{T_{c,w}} x(T)\,dT

    over :math:`x \in \{c_p, \rho, \mu, c_p/k\}` **at static pressure** --
    note it is :math:`c_p/k` that is integrated, not :math:`k` separately.
    Evaluated with bulk properties instead, it is just the 0.023 form and
    carries none of the paper's claim.

    **COEFFICIENT CONFLICT, UNRESOLVED.** Locke & Landrum (2008) are widely
    cited as giving this correlation with **0.025**. TN D-7207 was opened and
    its own new correlation is **0.023**; the string "0.025" appears in that
    report's INTRODUCTION describing the prior practice it set out to improve
    on. Locke & Landrum could not be opened (redirect loop / paywall) to
    determine whether they mis-transcribed or deliberately re-fit.
    **8.7% in Nu is at stake.** The primary's 0.023 is used and the result
    carries `CONFLICT_UNRESOLVED`.

    **Ranges and accuracy are NOT_REPORTED by the primary.** TN D-7207 gives
    no consolidated Re, pressure or temperature envelope -- the conditions
    live in per-station tables -- and states only that the correlation *"does
    a reasonable job"*, noting it UNDER-predicts at the throat. The
    frequently-quoted +/-56% for 95% coverage is **Locke & Landrum's statistic
    on their own database**, not the report's.
    """
    inputs = {"Re": reynolds, "Pr": prandtl,
              "integral_averaged": properties_are_integral_averaged}
    if reynolds <= 0 or prandtl <= 0:
        return Result(float("nan"), "-", "CAN-COOL-NU-SCHACHT",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="Re and Pr must be > 0")
    nu = 0.023 * reynolds ** 0.8 * prandtl ** 0.4

    statuses = [Status.PASS, Status.CONFLICT_UNRESOLVED]
    notes = [f"coefficient CONFLICT_UNRESOLVED: TN D-7207 prints 0.023; "
             f"Locke & Landrum are cited as 0.025 (+8.7% in Nu). The primary "
             f"was opened, the secondary was not, so 0.023 is used."]
    if not properties_are_integral_averaged:
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
        notes.append(
            "PROPERTIES ARE NOT INTEGRAL-AVERAGED. The integral averaging of "
            "{c_p, rho, mu, c_p/k} between the coolant static and wall "
            "temperatures IS this correlation's contribution; without it "
            "this is the bare 0.023 form and Schacht & Quentmeyer should not "
            "be cited for it.")

    return Result(
        nu, "-", "CAN-COOL-NU-SCHACHT", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=SCHACHT_QUENTMEYER_1973,
        source_locator="Eq. (7) with m=0.6, C=0.023, n=-0.2; property "
                       "averaging at Eq. (8), p.7. OPENED via NTRS.",
        validity_domain="HYDROGEN near the pseudo-critical region, in a "
                        "fired H2/O2 thrust chamber. Numerical Re, pressure "
                        "and temperature ranges are NOT_REPORTED by the "
                        "primary -- the report publishes no applicability "
                        "envelope.",
        uncertainty="NOT_REPORTED as a percentage by the primary; it states "
                    "only that the correlation 'does a reasonable job' and "
                    "under-predicts at the throat. The +/-56% for 95% "
                    "coverage is Locke & Landrum's statistic, not this "
                    "report's.",
        notes="; ".join(notes), inputs=inputs)


# ===========================================================================
# 6. THE SELECTOR -- and its right to refuse
# ===========================================================================

@dataclass
class CoolantCorrelationChoice:
    """What the selector chose, or its refusal, with the reasoning attached."""

    correlation: str | None
    result: Result
    considered: dict[str, str] = field(default_factory=dict)


def select_coolant_correlation(
    *, fluid: str, pressure_Pa: float, bulk_temperature_K: float,
    wall_temperature_K: float | None, reynolds: float, prandtl: float,
    relative_roughness: float | None, friction_factor: float | None,
    heat_flux_W_m2: float | None = None,
    mass_flux_kg_m2s: float | None = None,
) -> CoolantCorrelationChoice:
    """Choose a coolant-side correlation, or REFUSE.

    The rule this encodes is the one the thermal wave's audit implies: a
    correlation may be selected only when the state is inside ITS OWN stated
    domain. Applying a hydrogen correlation to methane, a smooth-tube
    correlation to a fully-rough channel, or a water-fitted rough-wall
    correlation four times outside its tested roughness are all failures of
    the same kind, and the selector refuses all three rather than ranking
    them.

    It returns `INSUFFICIENT_EVIDENCE` far more often than it returns a
    correlation. That is the intended behaviour, not a defect: for
    supercritical methane in a rough rocket channel there is currently **no
    correlation in the open literature that is simultaneously inside its
    fluid domain, its roughness domain and its Reynolds domain.**
    """
    considered: dict[str, str] = {}
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "T_wall_K": wall_temperature_K,
              "Re": reynolds, "Pr": prandtl,
              "eps_over_D": relative_roughness, "f": friction_factor}
    statuses: list[Status] = []
    notes: list[str] = []

    # --- regime ----------------------------------------------------------
    regime = coolant_regime_map(fluid=fluid, pressure_Pa=pressure_Pa,
                                bulk_temperature_K=bulk_temperature_K,
                                wall_temperature_K=wall_temperature_K)
    statuses.extend(regime.all_statuses)
    inputs["regime"] = regime.value
    if regime.has(Status.INVALID_MODEL_REGIME):
        notes.append(f"regime {regime.value}: every constant-property "
                     f"correlation is outside its assumptions")

    # --- HTD -------------------------------------------------------------
    htd_state = "NOT_EVALUATED"
    if heat_flux_W_m2 is not None and mass_flux_kg_m2s is not None:
        h = htd_onset(heat_flux_W_m2=heat_flux_W_m2,
                      mass_flux_kg_m2s=mass_flux_kg_m2s,
                      inlet_pressure_Pa=pressure_Pa)
        statuses.extend(h.all_statuses)
        htd_state = "HTD_DETECTED" if h.value is True else \
                    "NORMAL_SUPERCRITICAL_HEAT_TRANSFER"
        if h.value is True:
            statuses.append(Status.INSUFFICIENT_EVIDENCE)
            notes.append(
                "HTD_DETECTED. No accessible source quantifies the "
                "deteriorated heat-transfer coefficient -- the Pizzarelli "
                "correlation that would is behind a 403 and none of its "
                "coefficients is available. NO correlation may be selected "
                "here.")
    else:
        htd_state = "HTD_ONSET_UNKNOWN"
        notes.append("heat flux and mass flux not supplied, so the HTD "
                     "criterion could not be evaluated")
    inputs["HTD_status"] = htd_state

    # --- roughness -------------------------------------------------------
    rough_label = "UNKNOWN"
    if relative_roughness is not None and friction_factor is not None:
        rr = roughness_regime(reynolds=reynolds,
                              relative_roughness=relative_roughness,
                              friction_factor=friction_factor)
        statuses.extend(rr.all_statuses)
        rough_label = str(rr.value)
        inputs["h_s_plus"] = rr.inputs.get("h_s_plus")
    else:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append(
            "roughness or friction factor not supplied. The channel's "
            "roughness regime is the single largest determinant of which "
            "correlation applies -- an AM channel at eps/D = 0.013-0.105 is "
            "FULLY ROUGH, and every smooth-tube correlation is out of domain "
            "there. It is not assumed smooth.")
    inputs["roughness_regime"] = rough_label

    # --- candidates ------------------------------------------------------
    is_methane = fluid.strip().upper().replace("-", "") in (
        "CH4", "LCH4", "GCH4", "METHANE", "LNG")

    if rough_label == RoughnessRegime.FULLY_ROUGH.value:
        considered["Dipprey_Sabersky_1962"] = (
            "fully-rough branch is correct, BUT fitted on DISTILLED WATER at "
            "Pr 1.20-5.94 and eps/D <= 0.0488")
        considered["Stimpson_2017"] = (
            "fitted on AM channels, BUT on AIR at Re 3e3-3e4, one to two "
            "orders below a rocket channel")
        considered["Gnielinski/McAdams/Taylor"] = (
            "SMOOTH-TUBE correlations, out of domain in a fully rough channel")
    elif rough_label == RoughnessRegime.HYDRAULICALLY_SMOOTH.value:
        considered["Gnielinski_1976"] = "smooth branch is correct"
        considered["Petukhov_1970"] = "smooth branch is correct"
        considered["Dipprey_Sabersky_1962"] = "rough-wall form is undefined "
        considered["Taylor_TN_D_4332"] = (
            "fitted on HYDROGEN" if is_methane else "fluid domain ok")
    else:
        considered["all"] = (
            "TRANSITIONAL roughness: neither the smooth-tube family nor the "
            "fully-rough form is inside its own domain")

    # --- decide ----------------------------------------------------------
    choice: str | None = None
    if (rough_label == RoughnessRegime.HYDRAULICALLY_SMOOTH.value
            and not regime.has(Status.INVALID_MODEL_REGIME)
            and htd_state == "NORMAL_SUPERCRITICAL_HEAT_TRANSFER"
            and friction_factor is not None
            and 3.0e3 <= reynolds <= 1.0e6 and 0.5 <= prandtl <= 200.0):
        choice = "GNIELINSKI_1976"
        notes.append("selected Gnielinski (1976) with the Filonenko friction "
                     "factor: smooth channel, no pseudo-critical crossing, "
                     "no HTD, Re and Pr inside the narrowest stated range")
    else:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append(
            "NO SUPPORTED CORRELATION. For supercritical methane in a rough "
            "rocket cooling channel there is currently no correlation in the "
            "open literature simultaneously inside its fluid domain, its "
            "roughness domain and its Reynolds domain. Dipprey-Sabersky is "
            "the right BRANCH and the wrong FLUID; Stimpson is the right "
            "MANUFACTURING PROCESS and the wrong REYNOLDS NUMBER; the "
            "smooth-tube family is the wrong WALL.")

    res = Result(
        choice or "NO_SUPPORTED_CORRELATION", "-", "CAN-COOL-SELECT",
        status=worst(statuses) if statuses else Status.PASS,
        all_statuses=frozenset(statuses) if statuses else frozenset(
            {Status.PASS}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.NONE, validation=Validation.NONE,
        source="selection policy derived from each candidate's OWN stated "
               "domain; no ranking, no preference, no tie-break",
        source_locator="see each correlation's source_locator",
        validity_domain="the intersection of the selected correlation's "
                        "stated domains",
        uncertainty="inherited from the selected correlation; UNKNOWN when "
                    "none is selected",
        notes="; ".join(notes), inputs=inputs)
    return CoolantCorrelationChoice(choice, res, considered)


def coolant_side_h_c(
    *, fluid: str, pressure_Pa: float, bulk_temperature_K: float,
    wall_temperature_K: float | None, mass_flux_kg_m2s: float,
    hydraulic_diameter_m: float, relative_roughness: float | None,
    friction_factor: float | None, heat_flux_W_m2: float | None = None,
) -> Result:
    """h_c from real-fluid properties, with the correlation chosen by evidence.

    This is the function the wall solver should call. It runs the property
    chain, the regime map, the HTD criterion and the selector, and returns
    either a film coefficient with full provenance or `INSUFFICIENT_EVIDENCE`
    with the reason.

    It does NOT fall back to a smooth-tube correlation when the selector
    refuses. That fallback is the entire behaviour this wave exists to remove:
    the previous chain silently used the 0.023 form everywhere, including in
    fully-rough channels where it is out of domain by a factor of ~2.9 against
    the fully-rough branch.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "T_wall_K": wall_temperature_K,
              "G_kg_m2s": mass_flux_kg_m2s, "D_h_m": hydraulic_diameter_m,
              "eps_over_D": relative_roughness, "f": friction_factor}
    if mass_flux_kg_m2s <= 0 or hydraulic_diameter_m <= 0:
        return Result(float("nan"), "W/(m^2.K)", "CAN-COOL-HC",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="mass flux and hydraulic diameter must be > 0")

    st = coolant_state(fluid=fluid, pressure_Pa=pressure_Pa,
                       temperature_K=bulk_temperature_K)
    if not st.result.usable or math.isnan(st.density_kg_m3):
        return Result(float("nan"), "W/(m^2.K)", "CAN-COOL-HC",
                      status=st.result.status,
                      all_statuses=st.result.all_statuses,
                      source=st.result.source, notes=st.result.notes,
                      inputs=inputs)

    reynolds = mass_flux_kg_m2s * hydraulic_diameter_m / st.viscosity_Pa_s
    choice = select_coolant_correlation(
        fluid=fluid, pressure_Pa=pressure_Pa,
        bulk_temperature_K=bulk_temperature_K,
        wall_temperature_K=wall_temperature_K, reynolds=reynolds,
        prandtl=st.prandtl, relative_roughness=relative_roughness,
        friction_factor=friction_factor, heat_flux_W_m2=heat_flux_W_m2,
        mass_flux_kg_m2s=mass_flux_kg_m2s)

    statuses = set(choice.result.all_statuses) | set(st.statuses)
    base = {**inputs, "Re": reynolds, "Pr": st.prandtl,
            "rho_kg_m3": st.density_kg_m3, "cp_J_kgK": st.specific_heat_J_kgK,
            "mu_Pa_s": st.viscosity_Pa_s, "k_W_mK": st.conductivity_W_mK,
            "regime": choice.result.inputs.get("regime"),
            "HTD_status": choice.result.inputs.get("HTD_status"),
            "roughness_regime": choice.result.inputs.get("roughness_regime"),
            "correlation": choice.correlation,
            "considered": choice.considered}

    if choice.correlation is None:
        return Result(
            float("nan"), "W/(m^2.K)", "CAN-COOL-HC",
            status=Status.INSUFFICIENT_EVIDENCE,
            all_statuses=frozenset(statuses | {Status.INSUFFICIENT_EVIDENCE}),
            evidence_level=EvidenceLevel.E0,
            verification=Verification.NONE, validation=Validation.NONE,
            source="no correlation is inside its own stated domain here",
            source_locator="see `considered` for each candidate's reason",
            validity_domain=UNKNOWN_DOMAIN,
            uncertainty=NOT_REPORTED,
            notes="NO FILM COEFFICIENT IS RETURNED. " + choice.result.notes,
            inputs=base)

    if friction_factor is None:
        return Result(float("nan"), "W/(m^2.K)", "CAN-COOL-HC",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=base,
                      notes="the selected correlation needs a friction "
                            "factor and none was supplied")

    from kryptonis.canonical.thermal import nusselt_gnielinski
    nu_res = nusselt_gnielinski(reynolds=reynolds, prandtl=st.prandtl,
                                friction_factor=friction_factor)
    statuses |= set(nu_res.all_statuses)
    if not nu_res.usable or math.isnan(float(nu_res)):
        return Result(float("nan"), "W/(m^2.K)", "CAN-COOL-HC",
                      status=nu_res.status,
                      all_statuses=frozenset(statuses),
                      notes=nu_res.notes, inputs=base)

    h_c = float(nu_res) * st.conductivity_W_mK / hydraulic_diameter_m
    return Result(
        h_c, "W/(m^2.K)", "CAN-COOL-HC", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=nu_res.source,
        source_locator=nu_res.source_locator,
        validity_domain=nu_res.validity_domain,
        uncertainty="inherits the correlation's uncertainty AND the property "
                    "chain's, in which the methane viscosity model is a "
                    "generalized friction-theory model with NO stated "
                    "uncertainty",
        assumptions=st.result.assumptions,
        notes=choice.result.notes, inputs={**base, "Nu": float(nu_res)})


# ===========================================================================
# 7. HYDRAULIC ROUGHNESS FROM A MEASURED PRESSURE DROP
# ===========================================================================
#
# The way OUT of the Ra -> eps deadlock. Instead of converting a metrology
# statistic into a hydraulic parameter through a mapping nobody agrees on,
# MEASURE the hydraulic parameter directly: a pressure drop and a flow rate
# determine the friction factor exactly, and inverting Colebrook then gives
# the equivalent sand-grain roughness the correlations actually want.
#
# This is what Latini/Fiore/Nasuti do, and it is why their epsilon is
# trustworthy while a converted Ra is not.

def friction_factor_from_pressure_drop(
    *, pressure_drop_Pa: float, length_m: float, hydraulic_diameter_m: float,
    density_kg_m3: float, mass_flux_kg_m2s: float) -> Result:
    r"""Darcy friction factor from a measured pressure drop. Definitional.

    .. math:: f = \frac{2\,\Delta p\,D_h\,\rho}{L\,G^2}

    from Darcy-Weisbach with :math:`G = \rho u`. **Exact** for a fully
    developed, adiabatic, constant-area, constant-density duct -- and each of
    those four words is a real restriction that the result carries, because a
    heated rocket cooling channel violates the first and the third.
    """
    inputs = {"dp_Pa": pressure_drop_Pa, "L_m": length_m,
              "D_h_m": hydraulic_diameter_m, "rho_kg_m3": density_kg_m3,
              "G_kg_m2s": mass_flux_kg_m2s}
    if (length_m <= 0 or hydraulic_diameter_m <= 0 or density_kg_m3 <= 0
            or mass_flux_kg_m2s <= 0):
        return Result(float("nan"), "-", "CAN-COOL-F-FROM-DP",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="length, diameter, density and mass flux must "
                            "all be > 0")
    if pressure_drop_Pa <= 0:
        return Result(float("nan"), "-", "CAN-COOL-F-FROM-DP",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="a duct cannot gain pressure through friction; a "
                            "non-positive measured drop means the "
                            "measurement includes something other than "
                            "friction (acceleration, elevation, an error)")
    f = (2.0 * pressure_drop_Pa * hydraulic_diameter_m * density_kg_m3
         / (length_m * mass_flux_kg_m2s ** 2))
    return Result(
        f, "-", "CAN-COOL-F-FROM-DP", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="Darcy-Weisbach, rearranged for f with G = rho u",
        source_locator="definitional",
        validity_domain="fully developed, ADIABATIC, constant-area, "
                        "constant-density duct. A HEATED channel violates the "
                        "adiabatic and constant-density conditions: the "
                        "acceleration pressure drop from the density fall is "
                        "NOT friction and is not separated here.",
        uncertainty="EXACT given the inputs; the measurement uncertainty is "
                    "the caller's. Note d(f)/f = d(dp)/dp + 2 d(G)/G -- the "
                    "mass-flux error enters DOUBLED.",
        inputs=inputs)


def hydraulic_roughness_from_pressure_drop(
    *, pressure_drop_Pa: float, length_m: float, hydraulic_diameter_m: float,
    density_kg_m3: float, mass_flux_kg_m2s: float, viscosity_Pa_s: float,
    tolerance: float = 1e-12, max_iter: int = 200) -> Result:
    r"""Equivalent sand-grain roughness from a MEASURED pressure drop.

    Two exact steps and one inversion:

    1. :math:`f` from Darcy-Weisbach (exact),
    2. :math:`Re = G D_h/\mu`,
    3. Colebrook inverted for :math:`\varepsilon/D` at that :math:`(f, Re)`:

    .. math::
        \frac{\varepsilon}{D} = 3.7\left(10^{-1/(2\sqrt f)}
                                 - \frac{2.51}{Re\sqrt f}\right)

    **THIS IS THE ANSWER TO THE Ra PROBLEM, AND IT IS NOT A CONVERSION.**
    `sand_grain_from_Ra` refuses because Adams, Stimpson and Flack & Schultz
    are different functions spanning 0.4x to 16x. This function does not
    convert anything: it measures the quantity the correlations actually
    take. The cost is that it needs hardware -- a flowed channel with a
    calibrated dP -- so it is available at TEST time, not at design time.

    Two honest limits, both carried in the result:

    * The closure is **Colebrook's**, so the epsilon returned is a
      *Colebrook-equivalent* roughness. Feed it to a different friction model
      and it is no longer the right number. It is a model-conditioned
      parameter, not a property of the surface.
    * A **heated** channel's measured drop contains an acceleration
      contribution from the density fall that is not friction. Using a
      heated-flow dP here overstates f and therefore overstates epsilon.
      Cold-flow calibration is the clean case, and the result says so.

    Returns `INVALID_MODEL_REGIME` when the measured f is at or below the
    hydraulically smooth value for that Re -- there is no real roughness that
    produces less friction than a smooth pipe, so the correct answer is that
    the surface is smooth (or the measurement is wrong), not a negative
    epsilon.
    """
    inputs = {"dp_Pa": pressure_drop_Pa, "L_m": length_m,
              "D_h_m": hydraulic_diameter_m, "rho_kg_m3": density_kg_m3,
              "G_kg_m2s": mass_flux_kg_m2s, "mu_Pa_s": viscosity_Pa_s}
    if viscosity_Pa_s <= 0:
        return Result(float("nan"), "m", "CAN-COOL-EPS-FROM-DP",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="viscosity must be > 0")
    f_res = friction_factor_from_pressure_drop(
        pressure_drop_Pa=pressure_drop_Pa, length_m=length_m,
        hydraulic_diameter_m=hydraulic_diameter_m,
        density_kg_m3=density_kg_m3, mass_flux_kg_m2s=mass_flux_kg_m2s)
    if not f_res.usable or math.isnan(float(f_res)):
        return Result(float("nan"), "m", "CAN-COOL-EPS-FROM-DP",
                      status=f_res.status, all_statuses=f_res.all_statuses,
                      notes=f_res.notes, inputs=inputs)
    f = float(f_res)
    re = mass_flux_kg_m2s * hydraulic_diameter_m / viscosity_Pa_s

    from kryptonis.canonical.thermal import friction_factor_colebrook

    statuses = [Status.PASS]
    notes = [f"measured f = {f:.6f}, Re = {re:.3e}"]
    if re < 4000.0:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append("Re < 4000: Colebrook is a turbulent closure and the "
                     "inversion is outside it")

    # Colebrook solved for eps/D at fixed (f, Re).
    sqrt_f = math.sqrt(f)
    smooth_term = 2.51 / (re * sqrt_f)
    rhs = 10.0 ** (-1.0 / (2.0 * sqrt_f))
    eps_over_d = 3.7 * (rhs - smooth_term)

    if eps_over_d <= 0.0:
        f_smooth = friction_factor_colebrook(reynolds=re,
                                             relative_roughness=0.0)
        return Result(
            float("nan"), "m", "CAN-COOL-EPS-FROM-DP",
            status=Status.INVALID_MODEL_REGIME,
            all_statuses=frozenset({Status.INVALID_MODEL_REGIME}),
            evidence_level=EvidenceLevel.E1,
            source=COLEBROOK_1939,
            source_locator="Colebrook inverted for eps/D; PRIMARY NOT_OPENED",
            validity_domain="turbulent, fully developed",
            notes=f"the measured friction factor {f:.6f} is at or below the "
                  f"hydraulically smooth value {float(f_smooth):.6f} at "
                  f"Re = {re:.3e}. There is no real roughness that produces "
                  f"LESS friction than a smooth pipe, so no epsilon is "
                  f"returned. Either the surface is hydraulically smooth, or "
                  f"the measured drop is not purely frictional.",
            inputs={**inputs, "f_measured": f, "Re": re,
                    "f_smooth_at_Re": float(f_smooth)})

    eps = eps_over_d * hydraulic_diameter_m
    h_plus = re * eps_over_d * math.sqrt(f / 8.0)
    notes.append(f"eps/D = {eps_over_d:.5f}, h_s+ = {h_plus:.1f}")

    # Round-trip check: does Colebrook reproduce the measured f from this eps?
    back = friction_factor_colebrook(reynolds=re,
                                     relative_roughness=eps_over_d)
    residual = abs(float(back) / f - 1.0) if back.usable else float("nan")
    if not (residual < 1e-6):
        statuses.append(Status.FAIL)
        notes.append(f"round trip failed: Colebrook at the inverted eps/D "
                     f"gives f = {float(back):.6f}, {residual:.2%} from the "
                     f"measured value")

    return Result(
        eps, "m", "CAN-COOL-EPS-FROM-DP", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source=f"Darcy-Weisbach (exact) composed with {COLEBROOK_1939} "
               f"inverted for eps/D",
        source_locator="Colebrook PRIMARY NOT_OPENED (ICE 403); the "
                       "inversion is exact algebra on the form this project "
                       "implements",
        validity_domain="turbulent, fully developed, ADIABATIC, "
                        "constant-area. The returned epsilon is a "
                        "COLEBROOK-EQUIVALENT roughness -- a model-conditioned "
                        "parameter, not a property of the surface.",
        uncertainty="d(f)/f = d(dp)/dp + 2 d(G)/G, and d(eps)/eps is a steep "
                    "function of d(f)/f in the fully rough regime. A HEATED "
                    "measurement additionally contains an acceleration "
                    "pressure drop that is not friction and is NOT removed "
                    "here.",
        assumptions=(Assumption(
            "closure_model", "Colebrook-White", "-",
            "the epsilon returned is conditioned on Colebrook; a different "
            "friction model would return a different epsilon from the same "
            "measurement", Status.UNVALIDATED_ASSUMPTION),),
        notes="; ".join(notes),
        inputs={**inputs, "f_measured": f, "Re": re,
                "eps_over_D": eps_over_d, "h_s_plus": h_plus,
                "round_trip_residual": residual})


# ===========================================================================
# 8. NORMAL AND DETERIORATED HEAT TRANSFER -- SEPARATE, NOT MULTIPLIED
# ===========================================================================
#
# The tempting shape is
#     h_c = normal_correlation * deterioration_factor
# and NO SOURCE ESTABLISHES THAT STRUCTURE for supercritical methane.
# Pizzarelli's 2016 abstract describes "a Nusselt number correlation able to
# describe the convective heat-transfer characteristics of supercritical flow
# exhibiting deterioration" -- which reads like a STANDALONE correlation, not
# a multiplier. That reading is an inference from phrasing and the paper is
# behind a 403, so the structural question is genuinely open and the two
# regimes are kept as separate models with separate evidence.

def normal_supercritical_ht(
    *, fluid: str, pressure_Pa: float, bulk_temperature_K: float,
    wall_temperature_K: float | None, mass_flux_kg_m2s: float,
    hydraulic_diameter_m: float, relative_roughness: float | None,
    friction_factor: float | None, heat_flux_W_m2: float | None = None,
) -> Result:
    """Film coefficient for NON-deteriorated supercritical flow.

    Thin wrapper over `coolant_side_h_c` that exists so the two regimes are
    separately nameable. It carries the same refusal: no correlation, no
    number.
    """
    r = coolant_side_h_c(
        fluid=fluid, pressure_Pa=pressure_Pa,
        bulk_temperature_K=bulk_temperature_K,
        wall_temperature_K=wall_temperature_K,
        mass_flux_kg_m2s=mass_flux_kg_m2s,
        hydraulic_diameter_m=hydraulic_diameter_m,
        relative_roughness=relative_roughness,
        friction_factor=friction_factor, heat_flux_W_m2=heat_flux_W_m2)
    return Result(
        r.value, r.units, "CAN-COOL-HC-NORMAL", status=r.status,
        all_statuses=r.all_statuses, evidence_level=r.evidence_level,
        verification=r.verification, validation=r.validation,
        source=r.source, source_locator=r.source_locator,
        validity_domain=r.validity_domain, uncertainty=r.uncertainty,
        assumptions=r.assumptions,
        notes="NORMAL (non-deteriorated) branch. " + r.notes, inputs=r.inputs)


def deteriorated_supercritical_ht(
    *, fluid: str, pressure_Pa: float, bulk_temperature_K: float,
    mass_flux_kg_m2s: float, heat_flux_W_m2: float,
    hydraulic_diameter_m: float | None = None,
) -> Result:
    """Film coefficient DURING heat-transfer deterioration. Returns UNKNOWN.

    **There is no accessible source that quantifies this.** The one paper
    that would -- Pizzarelli, *Numerical Heat Transfer Part A* 69(3), 2016,
    242-264, "A CFD-derived correlation for methane heat transfer
    deterioration" -- was pursued through Unpaywall, OpenAlex, OpenAIRE,
    Crossref, Semantic Scholar, IRIS Sapienza, CORE, BASE, fatcat, EUCASS and
    a survey of every paper found citing it. Result:

        Unpaywall:  is_oa = false, oa_status = "closed",
                    has_repository_copy = false, oa_locations = []
        IRIS:       one attachment, Pizzarelli_a-cfd_2016.pdf, 2.63 MB,
                    access = "Restricted (Archive managers only)"
        secondary:  NOT ONE paper reprints the correlation. Every citing
                    work found describes it in words only.

    What IS known, from the abstract alone: the correlation is fitted on CFD
    of methane in a **heated tube** (circular), the solver was validated on
    near-critical **hydrogen** in heated tubes, and it is for flow "exhibiting
    deterioration and **negligible buoyancy effects**". Whether it is a
    standalone Nusselt correlation or a multiplier on a normal-flow
    correlation is `NOT_OPENED` and is NOT inferred here.

    So this function computes nothing. It returns `INSUFFICIENT_EVIDENCE`
    with the exact document required to close it. Inserting a plausible
    knock-down factor -- 0.7, 0.5, anything -- would be the invented HTD
    multiplier the brief forbids.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "G_kg_m2s": mass_flux_kg_m2s,
              "q_W_m2": heat_flux_W_m2, "D_h_m": hydraulic_diameter_m}
    return Result(
        float("nan"), "W/(m^2.K)", "CAN-COOL-HC-DETERIORATED",
        status=Status.INSUFFICIENT_EVIDENCE,
        all_statuses=frozenset({Status.INSUFFICIENT_EVIDENCE}),
        evidence_level=EvidenceLevel.E0,
        verification=Verification.NONE, validation=Validation.NONE,
        source="NO ACCESSIBLE SOURCE. Required: Pizzarelli, M., 'A "
               "CFD-derived correlation for methane heat transfer "
               "deterioration', Numerical Heat Transfer Part A 69(3), 2016, "
               "pp.242-264, DOI 10.1080/10407782.2015.1080575",
        source_locator="NOT_OPENED. Taylor & Francis 403; Unpaywall reports "
                       "oa_status='closed' with no repository copy; the IRIS "
                       "Sapienza copy is 'Restricted (Archive managers "
                       "only)'. A request-a-copy form exists at "
                       "iris.uniroma1.it/request-item?handle=11573/868669 "
                       "and requires a human to submit it.",
        validity_domain="the source is fitted for methane in a CIRCULAR "
                        "HEATED TUBE with NEGLIGIBLE BUOYANCY (from its "
                        "abstract). Whether it applies to a rectangular "
                        "one-sided-heated rocket channel is unknown.",
        uncertainty=NOT_REPORTED,
        notes="HTD_MAGNITUDE_UNKNOWN. The onset is known and this is not. No "
              "knock-down factor is applied: a plausible multiplier would be "
              "an invented number wearing a correlation's name.",
        inputs=inputs)


@dataclass
class CoolantHeatTransfer:
    """The full coolant answer: regime, both branches, and what is unknown."""

    regime: str
    htd_status: str
    normal_h_W_m2K: float
    deteriorated_h_W_m2K: float
    result: Result
    normal: Result
    deteriorated: Result
    htd: Result


def coolant_heat_transfer(
    *, fluid: str, pressure_Pa: float, bulk_temperature_K: float,
    wall_temperature_K: float | None, mass_flux_kg_m2s: float,
    hydraulic_diameter_m: float, heat_flux_W_m2: float,
    relative_roughness: float | None = None,
    friction_factor: float | None = None,
    channel_geometry: str = "CIRCULAR", **htd_geometry: Any,
) -> CoolantHeatTransfer:
    """THE coolant-side entry point. Never collapses the two regimes.

    Returns `regime`, `htd_status`, a normal-branch film coefficient (or
    refusal) and a deteriorated-branch film coefficient (which is always
    UNKNOWN today). The roll-up `result` is `INSUFFICIENT_EVIDENCE` whenever
    HTD is predicted, because at that point the model knows the design is in
    trouble and cannot say how much -- which is a worse position to be in
    than not knowing either, and must not be reported as a number.
    """
    htd = htd_onset(heat_flux_W_m2=heat_flux_W_m2,
                    mass_flux_kg_m2s=mass_flux_kg_m2s,
                    inlet_pressure_Pa=pressure_Pa,
                    channel_geometry=channel_geometry,
                    hydraulic_diameter_m=hydraulic_diameter_m,
                    **htd_geometry)
    if htd.value is True:
        htd_status = "HTD_DETECTED"
    elif htd.value is False and htd.status is not Status.INSUFFICIENT_EVIDENCE:
        htd_status = "NORMAL_SUPERCRITICAL_HEAT_TRANSFER"
    else:
        htd_status = "HTD_ONSET_UNKNOWN"

    regime = coolant_regime_map(fluid=fluid, pressure_Pa=pressure_Pa,
                                bulk_temperature_K=bulk_temperature_K,
                                wall_temperature_K=wall_temperature_K)
    normal = normal_supercritical_ht(
        fluid=fluid, pressure_Pa=pressure_Pa,
        bulk_temperature_K=bulk_temperature_K,
        wall_temperature_K=wall_temperature_K,
        mass_flux_kg_m2s=mass_flux_kg_m2s,
        hydraulic_diameter_m=hydraulic_diameter_m,
        relative_roughness=relative_roughness,
        friction_factor=friction_factor, heat_flux_W_m2=heat_flux_W_m2)
    det = deteriorated_supercritical_ht(
        fluid=fluid, pressure_Pa=pressure_Pa,
        bulk_temperature_K=bulk_temperature_K,
        mass_flux_kg_m2s=mass_flux_kg_m2s, heat_flux_W_m2=heat_flux_W_m2,
        hydraulic_diameter_m=hydraulic_diameter_m)

    statuses = set(htd.all_statuses) | set(regime.all_statuses) \
        | set(normal.all_statuses)
    notes = [f"regime={regime.value}", f"htd={htd_status}"]
    if htd_status == "HTD_DETECTED":
        statuses |= set(det.all_statuses)
        notes.append("HTD is predicted and its MAGNITUDE IS UNKNOWN, so no "
                     "film coefficient is reported. The normal-branch value "
                     "is retained for inspection only and must not be used "
                     "as the design h_c.")
        value = float("nan")
    else:
        value = normal.value
        notes.append(normal.notes)

    res = Result(
        value, "W/(m^2.K)", "CAN-COOL-HT", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=min(normal.evidence_level, EvidenceLevel.E1,
                           key=lambda e: e.value),
        verification=normal.verification, validation=Validation.NONE,
        source=normal.source, source_locator=normal.source_locator,
        validity_domain=normal.validity_domain,
        uncertainty=normal.uncertainty,
        assumptions=normal.assumptions, notes="; ".join(n for n in notes if n),
        inputs={"regime": regime.value, "HTD_status": htd_status,
                "normal_h_W_m2K": normal.value,
                "deteriorated_h_W_m2K": "UNKNOWN",
                "htd_threshold_J_kg": htd.inputs.get("threshold_used_J_kg"),
                "q_over_G_J_kg": htd.inputs.get("q_over_G_J_kg"),
                **{k: v for k, v in normal.inputs.items()
                   if k in ("Re", "Pr", "correlation", "roughness_regime")}})

    return CoolantHeatTransfer(
        regime=str(regime.value), htd_status=htd_status,
        normal_h_W_m2K=(float(normal) if normal.usable else float("nan")),
        deteriorated_h_W_m2K=float("nan"),
        result=res, normal=normal, deteriorated=det, htd=htd)


# ===========================================================================
# 9. UNCERTAINTY DECOMPOSITION -- computed, not assumed
# ===========================================================================

#: Stated (or, where the source states none, explicitly UNKNOWN) relative
#: uncertainties of every input the coolant film coefficient depends on.
#: `None` means the source reports NO uncertainty, and that is not the same
#: as zero.
_INPUT_UNCERTAINTY: dict[str, tuple[float | None, str]] = {
    "rho": (0.0015, "Setzmann & Wagner: +/-0.03% below 12 MPa and 350 K, "
                    "+/-0.03-0.15% above; the worst case is taken"),
    "cp": (0.01, "Setzmann & Wagner: '+/-1% generally'"),
    "k": (0.015, "Friend, Ely & Ingham: 'about 1.5%'"),
    "mu": (None, "NOT_REPORTED. The backend uses a GENERALIZED "
                 "friction-theory model, not a methane-specific reference "
                 "correlation, and states no uncertainty anywhere."),
    "f": (0.0142, "Haaland's own claim is NOT_REPORTED; this is a "
                  "THIRD-PARTY evaluation (Jaric et al., FME Trans.): max "
                  "+1.420%/-1.314% vs Colebrook over Re 4e3-1e8"),
    "eps": (None, "the Ra -> eps conversion spans 0.4x to 16x across three "
                  "structurally incompatible sources, and none states an "
                  "uncertainty. A MEASURED hydraulic epsilon replaces this "
                  "with the pressure-drop measurement's own uncertainty."),
}


def uncertainty_decomposition(
    *, fluid: str = "Methane", pressure_Pa: float, bulk_temperature_K: float,
    mass_flux_kg_m2s: float, hydraulic_diameter_m: float,
    relative_roughness: float, viscosity_hypotheses: tuple = (0.02, 0.05, 0.10),
    roughness_span: tuple = (0.4, 16.0),
) -> Result:
    r"""Decompose the coolant film coefficient's uncertainty, term by term.

    Each term is computed by PERTURBING that input and re-running the chain,
    not by assuming a sensitivity. The output is the relative change in the
    fully-rough Nusselt number (Dipprey-Sabersky, the correct branch for an
    AM channel) per stated input uncertainty.

    Two inputs have **no stated uncertainty**: the methane viscosity model,
    and the Ra -> epsilon conversion. For those the function reports a
    SENSITIVITY -- how much a 1 % input error costs -- and then propagates a
    range of hypotheses, so the reader can see that the total is unbounded
    because the *input* is unbounded, not because the arithmetic is loose.

    The brief asked for the dominant term to be verified quantitatively
    rather than assumed. It is.
    """
    inputs = {"fluid": fluid, "p_Pa": pressure_Pa,
              "T_bulk_K": bulk_temperature_K, "G_kg_m2s": mass_flux_kg_m2s,
              "D_h_m": hydraulic_diameter_m, "eps_over_D": relative_roughness}
    bad = _require_coolprop("CAN-COOL-UNCERTAINTY", inputs)
    if bad is not None:
        return bad
    if relative_roughness <= 0:
        return Result(float("nan"), "-", "CAN-COOL-UNCERTAINTY",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="this decomposition is for the fully-rough "
                            "branch and needs a positive roughness")

    from kryptonis.canonical.thermal import friction_factor_colebrook

    st = coolant_state(fluid=fluid, pressure_Pa=pressure_Pa,
                       temperature_K=bulk_temperature_K)
    if math.isnan(st.density_kg_m3):
        return Result(float("nan"), "-", "CAN-COOL-UNCERTAINTY",
                      status=st.result.status, inputs=inputs,
                      notes=st.result.notes)

    def _nu(mu=None, cp=None, k=None, eps=None, f_scale=1.0):
        mu = st.viscosity_Pa_s if mu is None else mu
        cp = st.specific_heat_J_kgK if cp is None else cp
        k = st.conductivity_W_mK if k is None else k
        eps = relative_roughness if eps is None else eps
        re = mass_flux_kg_m2s * hydraulic_diameter_m / mu
        pr = mu * cp / k
        f = float(friction_factor_colebrook(reynolds=re,
                                            relative_roughness=eps)) * f_scale
        r = nusselt_dipprey_sabersky(reynolds=re, prandtl=pr,
                                     relative_roughness=eps,
                                     friction_factor=f)
        return float(r) if r.usable else float("nan")

    base = _nu()
    if math.isnan(base) or base <= 0:
        return Result(float("nan"), "-", "CAN-COOL-UNCERTAINTY",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes="the baseline Nusselt number is not computable "
                            "at this state, so no decomposition is possible")

    def _rel(v):
        return abs(v / base - 1.0) if not math.isnan(v) else float("nan")

    terms: dict[str, dict[str, Any]] = {}

    # --- terms with a STATED uncertainty --------------------------------
    u_cp = _INPUT_UNCERTAINTY["cp"][0]
    terms["cp_EOS"] = {"input_uncertainty": u_cp,
                       "d_Nu_over_Nu": _rel(_nu(cp=st.specific_heat_J_kgK
                                                * (1 + u_cp))),
                       "basis": _INPUT_UNCERTAINTY["cp"][1]}
    u_k = _INPUT_UNCERTAINTY["k"][0]
    terms["k_transport"] = {"input_uncertainty": u_k,
                            "d_Nu_over_Nu": _rel(_nu(
                                k=st.conductivity_W_mK * (1 + u_k))),
                            "basis": _INPUT_UNCERTAINTY["k"][1]}
    u_f = _INPUT_UNCERTAINTY["f"][0]
    terms["friction_factor"] = {"input_uncertainty": u_f,
                                "d_Nu_over_Nu": _rel(_nu(f_scale=1 + u_f)),
                                "basis": _INPUT_UNCERTAINTY["f"][1]}
    # density does not enter Nu at fixed mass flux; recorded as exactly that
    terms["rho_EOS"] = {"input_uncertainty": _INPUT_UNCERTAINTY["rho"][0],
                        "d_Nu_over_Nu": 0.0,
                        "basis": _INPUT_UNCERTAINTY["rho"][1]
                                 + ". Density does not enter Nu at fixed MASS "
                                   "FLUX -- it enters the pressure drop and "
                                   "the bulk-temperature march instead."}

    # --- terms with NO stated uncertainty -------------------------------
    sens_mu = _rel(_nu(mu=st.viscosity_Pa_s * 1.01)) / 0.01
    terms["viscosity_model"] = {
        "input_uncertainty": None,
        "sensitivity_per_percent": sens_mu,
        "hypotheses": {f"{h:.0%}": _rel(_nu(mu=st.viscosity_Pa_s * (1 + h)))
                       for h in viscosity_hypotheses},
        "basis": _INPUT_UNCERTAINTY["mu"][1]}

    lo, hi = roughness_span
    terms["roughness_conversion"] = {
        "input_uncertainty": None,
        "sensitivity_per_percent": _rel(_nu(eps=relative_roughness * 1.01))
                                   / 0.01,
        "hypotheses": {"Ra_conversion_low": _rel(_nu(
                           eps=relative_roughness * lo)),
                       "Ra_conversion_high": _rel(_nu(
                           eps=relative_roughness * hi))},
        "basis": _INPUT_UNCERTAINTY["eps"][1]}

    # --- correlation-FORM uncertainty: not an input error at all ---------
    re0 = mass_flux_kg_m2s * hydraulic_diameter_m / st.viscosity_Pa_s
    pr0 = st.prandtl
    f0 = float(friction_factor_colebrook(reynolds=re0,
                                         relative_roughness=relative_roughness))
    nu_smooth = float(nusselt_petukhov_1970(
        reynolds=re0, prandtl=pr0,
        friction_factor=float(friction_factor_filonenko(reynolds=re0))))
    terms["correlation_form"] = {
        "input_uncertainty": None,
        "d_Nu_over_Nu": abs(nu_smooth / base - 1.0),
        "basis": "NOT an input error. This is the gap between the "
                 "fully-rough branch (Dipprey-Sabersky) and the smooth-tube "
                 "family (Petukhov) at the SAME state -- i.e. the cost of "
                 "choosing the wrong branch. It is a MODEL-SELECTION error "
                 "and does not combine in quadrature with the others."}

    # --- rank the terms that HAVE a stated uncertainty -------------------
    stated = {k: v["d_Nu_over_Nu"] for k, v in terms.items()
              if v.get("input_uncertainty") is not None}
    dominant_stated = max(stated, key=stated.get) if stated else None
    rss = math.sqrt(sum(v ** 2 for v in stated.values()))

    unbounded = [k for k, v in terms.items()
                 if v.get("input_uncertainty") is None
                 and k != "correlation_form"]

    notes = [
        f"baseline Nu = {base:.1f} (Dipprey-Sabersky, fully rough)",
        f"largest STATED-uncertainty term: {dominant_stated} at "
        f"{stated[dominant_stated] * 100:.2f}%" if dominant_stated else "",
        f"root-sum-square of all STATED terms: {rss * 100:.2f}%",
        f"UNBOUNDED terms (no source states an uncertainty): "
        f"{', '.join(unbounded)}",
        f"the roughness-conversion hypotheses alone span "
        f"{terms['roughness_conversion']['hypotheses']['Ra_conversion_low'] * 100:.0f}% "
        f"to "
        f"{terms['roughness_conversion']['hypotheses']['Ra_conversion_high'] * 100:.0f}% "
        f"-- an order of magnitude above every stated term combined, which is "
        f"why `sand_grain_from_Ra` refuses and "
        f"`hydraulic_roughness_from_pressure_drop` exists",
        f"correlation-FORM gap (wrong branch): "
        f"{terms['correlation_form']['d_Nu_over_Nu'] * 100:.0f}%",
    ]

    return Result(
        rss, "-", "CAN-COOL-UNCERTAINTY",
        status=Status.INSUFFICIENT_EVIDENCE,
        all_statuses=frozenset({Status.INSUFFICIENT_EVIDENCE}),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="each term's own source; see terms[*]['basis']",
        source_locator="computed by perturbation of the actual chain, not by "
                       "an assumed sensitivity",
        validity_domain="this state only; the sensitivities are local",
        uncertainty="THE TOTAL IS NOT A NUMBER. Two inputs have no stated "
                    "uncertainty, so the root-sum-square of the stated terms "
                    "is a LOWER BOUND and is reported as the value only for "
                    "that reason.",
        notes="; ".join(n for n in notes if n),
        inputs={**inputs, "baseline_Nu": base, "terms": terms,
                "dominant_stated_term": dominant_stated,
                "rss_of_stated_terms": rss, "unbounded_terms": unbounded})
