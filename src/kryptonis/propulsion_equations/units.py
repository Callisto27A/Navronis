"""
Explicit units and explicit status. RULE 5 and RULE 6, made mechanical.
=======================================================================

Two failure modes this module exists to stop, both of which have actually
happened in this repository:

  1. A variable called ``T`` that means kelvin in one module and rankine in
     another, or ``M`` that means g/mol here and kg/mol there. Bartz's own
     viscosity relation needs g/mol and kelvin and NOTHING in the code said so.

  2. A correlation transcribed with the wrong length unit. The contraction
     correlation was implemented with the throat diameter in INCHES when its
     source uses CENTIMETRES -- a 54% error that survived an audit because the
     unit convention lived in a comment rather than in the code.

Every canonical equation therefore takes SI, states its internal convention,
and converts explicitly at the boundary. `Q` is not a units library: it is a
named conversion table so that `to_cm(d)` appears at the call site and cannot
be mistaken for `d * 100`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

__all__ = [
    "INCH_M", "CM_M", "MM_M", "RANKINE_PER_KELVIN", "G_PER_MOL_TO_KG_PER_MOL",
    "to_inch", "to_cm", "to_mm", "to_rankine", "kelvin_from_rankine",
    "molar_mass_g_per_mol", "Status", "EvidenceLevel", "Verification",
    "Validation", "Assumption", "Result", "NOT_REPORTED", "UNKNOWN_DOMAIN",
]

# --- exact conversion factors ----------------------------------------------
INCH_M = 0.0254                      # exact, by definition
CM_M = 0.01
MM_M = 0.001
RANKINE_PER_KELVIN = 1.8             # exact
G_PER_MOL_TO_KG_PER_MOL = 1.0e-3

NOT_REPORTED = "NOT_REPORTED"
UNKNOWN_DOMAIN = "VALIDITY_DOMAIN_UNKNOWN"


def to_inch(metres: float) -> float:
    """Length [m] -> [in]. Named so a call site cannot hide the convention."""
    return metres / INCH_M


def to_cm(metres: float) -> float:
    """Length [m] -> [cm]."""
    return metres / CM_M


def to_mm(metres: float) -> float:
    return metres / MM_M


def to_rankine(kelvin: float) -> float:
    return kelvin * RANKINE_PER_KELVIN


def kelvin_from_rankine(rankine: float) -> float:
    return rankine / RANKINE_PER_KELVIN


def molar_mass_g_per_mol(kg_per_mol: float) -> float:
    """Molar mass [kg/mol] -> [g/mol].

    Bartz's viscosity relation is fitted in g/mol. Passing kg/mol gives an
    error of sqrt(1000) ~ 31.6x in mu, which is 2.0x in h_g.
    """
    return kg_per_mol / G_PER_MOL_TO_KG_PER_MOL


# --- status taxonomy (frozen; mirrors combustor.verification.status) -------

class Status(str, Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    #: The inputs describe a state that cannot exist. There is no number.
    PHYSICALLY_INVALID = "PHYSICALLY_INVALID"
    #: The state is real; THIS correlation does not describe it. A number is
    #: still returned, because the caller may legitimately want the
    #: extrapolation as long as it is labelled.
    OUT_OF_CORRELATION_RANGE = "OUT_OF_CORRELATION_RANGE"
    #: Nothing establishes the answer at these conditions.
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    #: A required input was defaulted rather than supplied.
    UNVALIDATED_ASSUMPTION = "UNVALIDATED_ASSUMPTION"
    #: Sources disagree and the disagreement is not resolved.
    CONFLICT_UNRESOLVED = "CONFLICT_UNRESOLVED"
    #: The model was asked about a regime it does not describe.
    INVALID_MODEL_REGIME = "INVALID_MODEL_REGIME"


_SEVERITY = [Status.PHYSICALLY_INVALID, Status.FAIL,
             Status.INVALID_MODEL_REGIME, Status.CONFLICT_UNRESOLVED,
             Status.OUT_OF_CORRELATION_RANGE, Status.INSUFFICIENT_EVIDENCE,
             Status.UNVALIDATED_ASSUMPTION, Status.PASS]


def worst(statuses) -> Status:
    ss = set(statuses)
    for s in _SEVERITY:
        if s in ss:
            return s
    return Status.PASS


class EvidenceLevel(str, Enum):
    E0 = "E0"
    E1 = "E1"
    E2 = "E2"
    E3 = "E3"
    E4 = "E4"


class Verification(str, Enum):
    NONE = "NONE"
    DIMENSIONAL = "DIMENSIONAL"
    ALGEBRAIC = "ALGEBRAIC"
    INTERNAL_CONSISTENCY = "INTERNAL_CONSISTENCY"
    INDEPENDENT_IMPLEMENTATION = "INDEPENDENT_IMPLEMENTATION"
    #: Agreement with another SOLVER. E1 by definition, always.
    CODE_TO_CODE = "CODE_TO_CODE"


class Validation(str, Enum):
    NONE = "NONE"
    LITERATURE_DATA = "LITERATURE_DATA"
    BLIND_TEST = "BLIND_TEST"
    FALSIFIED = "FALSIFIED"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"


@dataclass(frozen=True)
class Assumption:
    """A value the caller did not supply, and where it came from."""

    name: str
    value: Any
    units: str
    basis: str
    status: Status = Status.UNVALIDATED_ASSUMPTION

    def __str__(self) -> str:
        return f"{self.name}={self.value} {self.units} ({self.basis})"


@dataclass
class Result:
    """One canonical equation's answer, with everything needed to judge it.

    A bare float cannot say "this is 4 % outside the range the correlation was
    fitted over", and every place this repository returned a bare float is a
    place that information was lost.
    """

    value: Any
    units: str
    equation_id: str
    status: Status = Status.PASS
    evidence_level: EvidenceLevel = EvidenceLevel.E0
    verification: Verification = Verification.NONE
    validation: Validation = Validation.NONE
    source: str = ""
    source_locator: str = ""
    validity_domain: str = UNKNOWN_DOMAIN
    uncertainty: str = NOT_REPORTED
    assumptions: tuple[Assumption, ...] = ()
    notes: str = ""
    inputs: dict[str, Any] = field(default_factory=dict)
    #: EVERY condition that fired, not just the worst one. `status` is the
    #: roll-up and is what a gate should read; `all_statuses` is what a caller
    #: reads to ask a specific question. Without this, a CONFLICT_UNRESOLVED
    #: on the curvature exponent MASKS an OUT_OF_CORRELATION_RANGE on the area
    #: ratio -- two independent facts collapsed into one, which is the exact
    #: failure mode the status taxonomy exists to prevent.
    all_statuses: frozenset = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        if not self.all_statuses:
            object.__setattr__(self, "all_statuses", frozenset({self.status}))

    def __float__(self) -> float:
        return float(self.value)

    def has(self, status: "Status") -> bool:
        """Did this specific condition fire, whatever the roll-up says?"""
        return status in self.all_statuses

    @property
    def usable(self) -> bool:
        """False only where there is genuinely no number to use."""
        return self.status is not Status.PHYSICALLY_INVALID

    def as_trace(self) -> dict[str, Any]:
        return {
            "equation_id": self.equation_id,
            "value": self.value,
            "units": self.units,
            "status": self.status.value,
            "evidence_level": self.evidence_level.value,
            "verification": self.verification.value,
            "validation": self.validation.value,
            "source": self.source,
            "source_locator": self.source_locator,
            "validity_domain": self.validity_domain,
            "uncertainty": self.uncertainty,
            "all_statuses": sorted(x.value for x in self.all_statuses),
            "assumptions": [str(a) for a in self.assumptions],
            "inputs": self.inputs,
            "notes": self.notes,
        }
