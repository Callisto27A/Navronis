r"""Wall-conduction measurement-bias model. NOT an engine FEA.

===========================================================================
WHAT THIS MODULE IS FOR, AND WHAT IT IS NOT FOR
===========================================================================

The experimental design reconstructs the wetted-wall temperature from two
embedded sensors by **one-dimensional planar extrapolation**:

.. math:: T_w = T_1 - g x_1,\qquad g = \frac{T_2-T_1}{x_2-x_1}

and flags that reconstruction `UNVALIDATED_ASSUMPTION` because it neglects
channel-corner curvature, rib spreading, lateral conduction and axial
conduction. **This module exists to put a number on that neglect** --
nothing else.

It is therefore the *smallest* numerical conduction model that can answer
the question: a cell-centred finite-volume solver on a rectilinear grid
aligned exactly to the channel geometry, with temperature-dependent
conductivity and harmonic-mean face conductances.

**It is not an engine design tool.** It contains no coolant correlation, no
stress, no fatigue, no life. It cannot predict a flight wall temperature and
does not try to. Where a coolant-side film coefficient is needed it is
supplied **parametrically by the caller** and every result built on one
carries `INSUFFICIENT_EVIDENCE`, because the production coolant model still
returns `NO_SUPPORTED_CORRELATION` and nothing here changes that.

===========================================================================
THE GOVERNING EQUATION, AND WHAT IS ASSUMED ABOUT IT
===========================================================================

Steady conduction with temperature-dependent, isotropic conductivity and no
internal generation:

.. math:: \nabla\cdot\left[k(T)\,\nabla T\right] = 0

**Why q''' = 0 is justified here** and is not merely convenient: the
specimen is heated at an external surface by a separate heater element. If
the facility instead chooses **direct Joule heating of the specimen itself**
-- which is a permitted `DESIGN_CONVENTION` in the facility spec -- then
:math:`q''' = \rho_e J^2 \ne 0` inside the wall and this assumption is
**violated**. `solve_conduction` therefore accepts a volumetric source and
returns `UNVALIDATED_ASSUMPTION` when the caller leaves it at zero without
declaring the heating method. The assumption is recorded, not buried.

**The discretisation.** Cell-centred finite volume on a rectilinear grid.
Face conductivity is the **harmonic mean**, which for two half-cells of
equal thickness is :math:`k_f = 2k_Pk_E/(k_P+k_E)`. Its stated derivation
assumptions are steady state, equal flux through both half-cells, equal
half-thicknesses and one-dimensional flow across the face (Humphrey, Univ.
of Wyoming teaching note, OPENED; the same result is attributed throughout
the literature to Patankar, *Numerical Heat Transfer and Fluid Flow*, 1980,
Ch. 4 §"Steady One-dimensional Conduction" -> "Interface Conductivity",
**which this project could NOT open** -- see
`docs/research/WALL_CONDUCTION_PRIMARY_SOURCE_AUDIT.md`).

**A gap that is stated rather than papered over:** every rigorous accuracy
analysis of harmonic averaging that this project opened -- Kadioglu,
Nourgaliev & Mousseau (INL/EXT-08-13999, 2008) and Pan, Xu & Li
(arXiv:2502.09413, 2025) -- **assumes a uniform grid**. Neither treats
non-uniform spacing. The grids built here are uniform *within* each
geometric block but the block spacings differ, so the second-order claim is
`VALIDITY_DOMAIN_UNKNOWN` at block interfaces. That is why this module
verifies itself by **grid convergence on the actual geometry** and by an
**exact energy balance**, rather than by citing an order of accuracy it
cannot support.

===========================================================================
WHAT THIS MODULE FOUND THAT THE 1-D MODEL CANNOT SEE
===========================================================================

Two distinct biases, which must not be conflated:

**(i) Reconstruction bias.** The 1-D formula applied to the true field's own
sensor readings, versus the true field's wetted-wall temperature *at the
same lateral station*. This is the error of the formula.

**(ii) Definition bias.** The wetted-wall temperature at the sensor station,
versus the area-averaged wetted-wall temperature that the metered flux
:math:`q''=P/(wL)` is actually paired with in :math:`h_c=q''/(T_w-T_b)`.
This is the error of the *definition*, and it exists even if the formula
were exact.

They are reported separately and are never summed in quadrature: a bias is
not a random error and combining them as though they were independent
Gaussians would understate one and misrepresent both.
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

try:                                                        # pragma: no cover
    import numpy as _np
    from scipy.sparse import csr_matrix as _csr
    from scipy.sparse.linalg import spsolve as _spsolve
    NUMERICS_AVAILABLE = True
except Exception:                                           # pragma: no cover
    _np = None
    NUMERICS_AVAILABLE = False

__all__ = [
    "NUMERICS_AVAILABLE",
    "HeatingMethod", "SolidProperties", "GRCOP42_LPBF", "GRCOP84_EXTRUDED",
    "SensorInclusion", "sensor_reading",
    "ANISOTROPY_EVIDENCE",
    "ChannelSection", "ConductionSolution",
    "solve_conduction", "cross_section_solution", "segment_solution",
    "wetted_face_heat_split", "reconstruction_bias",
    "axial_leak_from_field", "guard_matching_requirement",
    "identify_h_c", "allowable_modelling_bias",
    "depth_error_sweep", "optimal_sensor_depths",
    "AnalyticField", "ANALYTIC_FIELDS", "verify_reconstruction_on_field",
    "conductivity_sensitivity",
    "conduction_heat_flux", "fourier_wall_temperature_drop", "cylindrical_wall_correction",
]

def conduction_heat_flux(*, wall_thickness_m: float, T_hot_K: float, T_cold_K: float, thermal_conductivity_W_m_K: float) -> Result:
    r""".. math:: q = \frac{k}{t_{wall}} (T_{hot} - T_{cold})"""
    if wall_thickness_m <= 0 or thermal_conductivity_W_m_K <= 0:
        return Result(float("nan"), "W/m^2", "CAN-Q-COND", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = thermal_conductivity_W_m_K * (T_hot_K - T_cold_K) / wall_thickness_m
    return Result(val, "W/m^2", "CAN-Q-COND", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="Fourier 1D heat conduction law", source_locator="Incropera & DeWitt",
                  validity_domain="steady 1D conduction", uncertainty="exact",
                  inputs={"wall_thickness_m": wall_thickness_m, "T_hot_K": T_hot_K, "T_cold_K": T_cold_K, "k": thermal_conductivity_W_m_K})

def fourier_wall_temperature_drop(*, heat_flux_W_m2: float, wall_thickness_m: float, thermal_conductivity_W_m_K: float) -> Result:
    r""".. math:: \Delta T = \frac{q \cdot t_{wall}}{k}"""
    if wall_thickness_m <= 0 or thermal_conductivity_W_m_K <= 0:
        return Result(float("nan"), "K", "CAN-DT-COND", status=Status.PHYSICALLY_INVALID, notes="Inputs must be > 0")
    val = (heat_flux_W_m2 * wall_thickness_m) / thermal_conductivity_W_m_K
    return Result(val, "K", "CAN-DT-COND", status=Status.PASS, evidence_level=EvidenceLevel.E2,
                  verification=Verification.ALGEBRAIC, validation=Validation.NONE,
                  source="Fourier 1D heat conduction law", source_locator="Incropera & DeWitt",
                  validity_domain="steady 1D conduction", uncertainty="exact",
                  inputs={"heat_flux_W_m2": heat_flux_W_m2, "wall_thickness_m": wall_thickness_m, "k": thermal_conductivity_W_m_K})

cylindrical_wall_correction = fourier_wall_temperature_drop



# ===========================================================================
# 1. MATERIAL -- INCLUDING A CORRECTION TO THIS PROJECT'S OWN EARLIER CLAIM
# ===========================================================================

class HeatingMethod(str, Enum):
    """How the specimen is heated. It decides whether q''' = 0 is legal."""

    #: A separate element bonded to or radiating onto the outer surface. The
    #: wall carries no internal generation and q''' = 0 is EXACT.
    EXTERNAL_SURFACE = "EXTERNAL_SURFACE"
    #: Current passed through the specimen itself. q''' = rho_e J^2 != 0 and
    #: a zero-source solve is WRONG, not approximate.
    DIRECT_JOULE = "DIRECT_JOULE"
    #: The facility has not said. The solve proceeds with q''' = 0 and says
    #: so, because assuming external heating silently is how a modelling
    #: error becomes invisible.
    NOT_DECLARED = "NOT_DECLARED"


@dataclass(frozen=True)
class SolidProperties:
    """k(T) for one material in one condition, with its domain attached."""

    name: str
    k: Callable[[float], float]
    valid_range_K: tuple[float, float]
    source: str
    source_locator: str
    source_access: str
    condition: str
    uncertainty: str = NOT_REPORTED
    #: k_axial / k_through_thickness. `None` means UNKNOWN -- see
    #: ANISOTROPY_EVIDENCE. It is NOT silently 1.0.
    anisotropy_ratio: float | None = None
    notes: str = ""

    def k_at(self, T: float) -> tuple[float, bool]:
        """Conductivity and whether T was inside the source's stated range."""
        lo, hi = self.valid_range_K
        inside = lo <= T <= hi
        Tc = min(max(T, lo), hi)
        return self.k(Tc), inside


#: Chen, Y., Zeng, C., Ding, H., Emanet, S., Gradl, P. R., Ellis, D. L. &
#: Guo, S., 'Thermophysical properties of additively manufactured (AM)
#: GRCop-42 and GRCop-84', Materials Today Communications 36, 106665, 2023,
#: DOI 10.1016/j.mtcomm.2023.106665, **Eq. (3)**.
#:
#: PRINTED VERBATIM (confirmed on three independent fetches of the
#: NSF-PAR open copy, https://par.nsf.gov/servlets/purl/10496372):
#:
#:     k(T) = - 39.86T3 + 4.17T2 + 0.97x + 329.14        (3)
#:     "where T is temperature in K/1000."
#:
#: THE PRINTED LINEAR TERM READS `0.97x` AND `x` IS NEVER DEFINED. The
#: sentence that follows the equation defines only T, and the paper's
#: companion specific-heat equations carry the same slipped variable name,
#: so this is a typesetting error in the source. Reading it as `0.97T` is
#: the only evaluable reading -- but it is a READING, and it is flagged.
#: Its consequence is bounded and tiny: over the whole stated range the
#: linear term spans 0.29-0.94 W/(m.K) out of ~300, i.e. **at most 0.3 %**,
#: so no engineering conclusion in this project can turn on it.
def _grcop42_lpbf_k(T_K: float) -> float:
    t = T_K / 1000.0
    return -39.86 * t ** 3 + 4.17 * t ** 2 + 0.97 * t + 329.14


#: **A CORRECTION TO THIS PROJECT'S OWN EARLIER STATEMENT.** The coolant
#: test-facility wave recorded `conductivity_source = NOT_REPORTED` for LPBF
#: GRCop-42 and asserted that no k(T) curve for it existed in any opened
#: source. **That assertion was wrong.** A deeper source sweep opened Chen
#: et al. (2023) in full through the NSF public-access repository. The
#: earlier statement is corrected here rather than quietly replaced.
#:
#: What the curve does and does not cover:
#:   * eight GRCop-42 samples from six vendors, five machine models
#:   * xenon flash (Netzsch LFA 467 HT), 25-700 degC
#:   * **all samples HIP'd.** An as-built (non-HIP) specimen is OUTSIDE
#:     this curve's condition domain.
#:   * sample-to-sample variation stated as "less than +/- 4 %", against an
#:     instrument uncertainty of 3 %
#:   * published as a FIT plus figures. **No tabulated k values exist in the
#:     paper**, so any point value is DERIVED, not REPORTED.
GRCOP42_LPBF = SolidProperties(
    name="GRCop-42 (L-PBF, HIP'd)",
    k=_grcop42_lpbf_k,
    valid_range_K=(298.15, 973.15),
    source="Chen, Y. et al., 'Thermophysical properties of additively "
           "manufactured (AM) GRCop-42 and GRCop-84', Materials Today "
           "Communications 36, 106665, 2023, DOI 10.1016/j.mtcomm.2023.106665",
    source_locator="Eq. (3), thermal-conductivity section; T in K/1000; "
                   "stated range 25-700 degC",
    source_access="OPENED in full via the NSF public-access repository "
                  "(par.nsf.gov/servlets/purl/10496372). The ScienceDirect "
                  "copy is robots-disallowed and the LSU repository copy "
                  "returns 403.",
    condition="L-PBF, HIP'd. As-built (non-HIP) material is OUT OF DOMAIN.",
    uncertainty="sample-to-sample variation stated as < +/-4 %; instrument "
                "uncertainty 3 %. The fit itself carries NO stated "
                "confidence interval.",
    anisotropy_ratio=None,
    notes="printed linear term reads '0.97x' with x undefined; read as "
          "0.97T, worth at most 0.3 % over the range")


def _grcop84_extruded_k(T_K: float) -> float:
    return 243.8 + 0.1792 * T_K - 1.325e-4 * T_K * T_K


#: NASA/CR-2000-210055 Eq. (15). **Powder metallurgy plus extrusion at a
#: 29.5:1 reduction -- NOT additively manufactured.** Retained only as a
#: contrast case and as the historical entry this project used before Chen
#: et al. was opened. It must not be substituted for an LPBF specimen.
GRCOP84_EXTRUDED = SolidProperties(
    name="GRCop-84 (extruded, powder metallurgy)",
    k=_grcop84_extruded_k,
    valid_range_K=(296.0, 1173.0),
    source="NASA/CR-2000-210055, D. L. Ellis, 'Thermophysical Properties of "
           "GRCop-84'",
    source_locator="Eq. (15), p.6",
    source_access="OPENED",
    condition="extruded, 29.5:1 area reduction. NOT additively manufactured.",
    uncertainty="mean regression over 15 specimens; a lower 95% CI is given "
                "separately as Eq. (16)",
    anisotropy_ratio=None,
    notes="different ALLOY and different PRODUCT FORM from the specimen; "
          "present for contrast only")


#: **Anisotropy of k in LPBF GRCop has never been measured.** What exists:
#:
#:   * the CAUSE is documented: Ellis & Jennings (NTRS 20230007103) report a
#:     strong build-axis crystallographic texture and a columnar Cu grain
#:     structure surviving HIP.
#:   * the EFFECT is not. Chen et al. (2023) never mention build orientation
#:     and never state where on the plate the eight disks were taken from.
#:   * the nearest alloys CONTRADICT each other. Biffi et al. (2024,
#:     DOI 10.1007/s40516-023-00240-7, OPENED) report LPBF CuCrZr k HIGHER
#:     along the growth direction. Xie et al. (JMRT 2023, ABSTRACT_ONLY,
#:     full text robots-disallowed) report heat-treated LPBF CuCrZr at
#:     255 W/(m.K) along the build direction against 307 transverse -- k
#:     17 % LOWER along build. Qu et al. (2022, DOI
#:     10.1016/j.addma.2022.103082, OPENED) measured LPBF PURE COPPER at
#:     0/45/90 degrees and found conductivity "near-isotropic".
#:
#: So the ratio is `CONFLICT_UNRESOLVED`, not merely unknown, and a single
#: value is NOT chosen. The bracket below is used ONLY for a sensitivity
#: sweep and is labelled as coming from a DIFFERENT ALLOY.
ANISOTROPY_EVIDENCE: dict[str, Any] = {
    "grcop_measured": False,
    "status": Status.CONFLICT_UNRESOLVED,
    "sensitivity_bracket": (0.83, 1.16),
    "bracket_basis": "0.83 from Xie et al. (255/307 on heat-treated LPBF "
                     "CuCrZr, ABSTRACT_ONLY, arithmetic on two published "
                     "numbers); 1.16 from Biffi et al. (about 110 vs 95 "
                     "W/(m.K) as-built LPBF CuCrZr at 25 degC, values the "
                     "authors themselves prefix 'about'). DIFFERENT ALLOY, "
                     "different heat treatment, and the two sources "
                     "disagree on the SIGN. This is a sweep range, not a "
                     "property.",
    "closure": "a laser-flash measurement on two coupons cut from the SAME "
               "GRCop-42 build, one with heat flow along the build axis and "
               "one transverse, would close it. Nothing else will.",
}


# ===========================================================================
# 2. GEOMETRY
# ===========================================================================

@dataclass(frozen=True)
class SensorInclusion:
    r"""An embedded thermocouple, represented as what it physically is.

    **The published thermocouple-disturbance literature does not cover this
    case.** Every correction method this project could open -- Pope et al.
    (2021), Beck's lineage, the fire-science literature -- is built for a
    RELATIVELY CONDUCTIVE sensor in a POORLY CONDUCTING solid, with
    :math:`K = k_{solid}/k_{sensor}` in the range 0.003-0.01. In GRCop-42
    at 344 W/(m.K) with a Type K MI sensor (Inconel 600 sheath 14.8, MgO
    ~3-30, chromel 19.2, alumel 29.7 W/(m.K)) that same ratio is **12 to
    115** -- three to four orders of magnitude outside every validated
    method, and on the other side of unity. **The sensor is an insulating
    inclusion, not a thermal bridge**, and no opened source treats that
    regime.

    Four things are therefore modelled explicitly rather than assumed:

    ``k_W_mK``
        the sensor's own conductivity. Not the wall's.
    ``contact_conductance_W_m2K``
        the interface between sensor and wall. **No published measured
        value exists for an embedded thermocouple in metal, by any
        installation method** -- press fit, peened, brazed, soldered,
        spark-welded or AM-embedded. It is a free parameter and it is
        swept, never chosen.
    ``lead_conductance_W_K``
        the fin path out along the wires. Boelter & Lockhart (NACA TN 2427)
        measured 20-200 degF of error from lead routing alone on a thin
        plate, falling to 20-60 degF when the wires were held against the
        plate for ~100 mm first. This is the term that turns a contact
        resistance into a bias.
    ``junction_offset_m``
        where the measuring junction sits inside the inclusion, measured
        from its centre, positive DEEPER. **This matters more than it
        looks**: a passive inclusion centred in a linear gradient reads the
        undisturbed centre temperature exactly, whatever its conductivity
        and whatever its contact resistance. The bias comes from the
        junction being off-centre and from the leads draining heat -- not
        from the contact resistance on its own.
    """

    name: str
    depth_m: float                      # junction depth below the wetted wall
    diameter_m: float
    k_W_mK: float
    contact_conductance_W_m2K: float | None = None
    lead_conductance_W_K: float = 0.0
    lead_exit_temperature_K: float | None = None
    junction_offset_m: float = 0.0
    lateral_position_m: float = 0.0
    installation: str = NOT_REPORTED

    def z_bounds(self, wall_thickness_m: float) -> tuple[float, float]:
        zc = wall_thickness_m - self.depth_m
        return zc - 0.5 * self.diameter_m, zc + 0.5 * self.diameter_m

    def y_bounds(self) -> tuple[float, float]:
        return (max(0.0, self.lateral_position_m - 0.5 * self.diameter_m),
                self.lateral_position_m + 0.5 * self.diameter_m)


@dataclass(frozen=True)
class ChannelSection:
    """Half-pitch cross-section: heated face, wall, channel, rib, backing.

    Lateral symmetry is exploited, so the domain runs from the channel
    centreline (y=0) to the rib centreline (y=L_HG/2). Both are adiabatic
    symmetry planes, which is exact for an infinite periodic array of
    identical channels -- the same idealisation the 1-D model makes, so the
    comparison isolates the conduction dimensionality and nothing else.
    """

    channel_width_m: float
    channel_height_m: float
    wall_thickness_m: float          # heated face -> channel top
    rib_width_m: float
    backing_height_m: float          # channel bottom -> back face

    @property
    def half_pitch_m(self) -> float:
        return 0.5 * (self.channel_width_m + self.rib_width_m)

    @property
    def half_width_m(self) -> float:
        """Half the CHANNEL width -- **not** half the domain width.

        The domain runs to `half_pitch_m`. These differ by the rib, and
        confusing them understates the solid cross-section by 42 % at this
        geometry. Use `solid_cross_section_m2` rather than reassembling it.
        """
        return 0.5 * self.channel_width_m

    @property
    def total_height_m(self) -> float:
        return (self.wall_thickness_m + self.channel_height_m
                + self.backing_height_m)

    @property
    def solid_cross_section_m2(self) -> float:
        """Solid metal area of the half-pitch domain, channel removed.

        **THE ONE CANONICAL IMPLEMENTATION.** It sets the axial conduction
        leak, the segment-length floor and -- for a volumetric heater --
        the generation density. Every one of those is proportional to it,
        so an ad-hoc reassembly is a first-order error.
        """
        return (self.half_pitch_m * self.total_height_m
                - self.half_width_m * self.channel_height_m)

    def validate(self) -> list[str]:
        bad = []
        for n, v in (("channel_width", self.channel_width_m),
                     ("channel_height", self.channel_height_m),
                     ("wall_thickness", self.wall_thickness_m),
                     ("rib_width", self.rib_width_m),
                     ("backing_height", self.backing_height_m)):
            if v <= 0:
                bad.append(f"{n} must be > 0")
        return bad


def _grid(edges: Sequence[float], counts: Sequence[int],
          extra_edges: Sequence[float] = ()) -> "Any":
    """Face coordinates: uniform inside each block, exact at every edge.

    ``extra_edges`` are additional interfaces that must fall exactly on a
    grid line -- the faces of an embedded sensor, for instance. Each parent
    block is subdivided at them and the parent's cell count is distributed
    in proportion to length, so refining for a sensor never coarsens the
    rest of the wall.
    """
    lo, hi = edges[0], edges[-1]
    extra = sorted({float(e) for e in extra_edges
                    if lo + 1e-15 < e < hi - 1e-15})
    out = [lo]
    for i, n in enumerate(counts):
        a, b = edges[i], edges[i + 1]
        inner = [e for e in extra if a + 1e-15 < e < b - 1e-15]
        bounds = [a] + inner + [b]
        span = b - a
        for j in range(len(bounds) - 1):
            p, q = bounds[j], bounds[j + 1]
            m = max(1, int(round(n * (q - p) / span)))
            for t in range(1, m + 1):
                out.append(p + (q - p) * t / m)
    return _np.array(sorted(set(out)), dtype=float)


# ===========================================================================
# 3. THE SOLVER
# ===========================================================================

@dataclass
class ConductionSolution:
    """A solved temperature field, plus everything needed to distrust it."""

    T: Any                      # (nx, ny, nz) with NaN in fluid cells
    xf: Any                     # axial face coordinates   (nx+1,)
    yf: Any                     # lateral face coordinates (ny+1,)
    zf: Any                     # through-thickness faces  (nz+1,)
    solid: Any                  # boolean (nx, ny, nz)
    section: ChannelSection
    heat_in_W: float
    heat_out_W: float
    energy_residual: float
    picard_iterations: int
    converged: bool
    result: Result
    #: per-wetted-face power, keyed TOP / SIDE / BOTTOM
    face_power_W: dict[str, float] = field(default_factory=dict)
    #: cell -> sensor index, or -1. Set when sensors were modelled.
    sensor_id: Any = None
    sensors: tuple = ()

    # -- sampling -------------------------------------------------------
    def _cell(self, xf, x):
        i = int(_np.searchsorted(xf, x) - 1)
        return min(max(i, 0), len(xf) - 2)

    def at(self, *, x: float, y: float, z: float) -> float:
        """Temperature at a point, linearly interpolated **through solid only**.

        A sensor sits at a definite depth, and sampling the nearest cell
        centre instead makes the reading jump as the grid changes -- which
        turns a grid-refinement study into noise. So this interpolates
        linearly in z between the two bracketing solid cell centres in the
        same column, and falls back to the nearest cell whenever a
        bracketing neighbour is fluid. It never interpolates across a
        solid/fluid boundary, because that would manufacture a temperature
        inside the coolant.
        """
        i = self._cell(self.xf, x)
        j = self._cell(self.yf, y)
        zc = 0.5 * (self.zf[:-1] + self.zf[1:])
        k = self._cell(self.zf, z)
        if not self.solid[i, j, k]:
            return float("nan")
        if z >= zc[k]:
            k2 = k + 1
        else:
            k2 = k - 1
        if not (0 <= k2 < len(zc)) or not self.solid[i, j, k2]:
            return float(self.T[i, j, k])
        t1, t2 = float(self.T[i, j, k]), float(self.T[i, j, k2])
        w = (z - zc[k]) / (zc[k2] - zc[k])
        return t1 + w * (t2 - t1)

    def _h_and_bulk(self) -> tuple[Callable[[float], float] | None,
                                    Callable[[float], float]]:
        """The film coefficient AS A FUNCTION OF x, and the bulk callable.

        The film coefficient may vary axially, so this returns the callable
        rather than a scalar. Every wetted-surface reconstruction must
        evaluate it at the station it is reconstructing, or an imposed h_c
        perturbation would be silently reconstructed with the wrong h.
        """
        h = self.result.inputs.get("_h_callable")
        tb = self.result.inputs.get("_bulk_callable")
        return h, tb

    def wetted_surface_temperature(self, *, x: float, y: float) -> float:
        r"""Temperature ON the wetted top surface, not at a cell centre.

        Taken from the converged face flux rather than from the nearest cell,
        so it is grid-independent to the same order as the solve itself:

        .. math::
            q_f = \frac{T_P - T_b}{1/h + \delta/2k},\qquad
            T_s = T_b + \frac{q_f}{h}

        Reading a cell centre instead would report a temperature half a cell
        *inside the metal*, which at 160 K/mm is tens of kelvin and moves
        every time the grid changes.
        """
        i = self._cell(self.xf, x)
        j = self._cell(self.yf, y)
        s = self.section
        k = self._cell(self.zf, s.wall_thickness_m - 1e-12)
        if not self.solid[i, j, k]:
            return float("nan")
        hf, tbf = self._h_and_bulk()
        if hf is None or tbf is None:
            return float("nan")
        xc = 0.5 * (self.xf[i] + self.xf[i + 1])
        h = float(hf(float(xc)))
        if not (h == h) or h <= 0.0:
            return float("nan")
        kP = self._k_of(float(self.T[i, j, k]))
        dz = self.zf[k + 1] - self.zf[k]
        gg = 1.0 / (1.0 / h + 0.5 * dz / max(kP, 1e-30))
        tb = float(tbf(float(xc)))
        q_f = gg * (float(self.T[i, j, k]) - tb)
        return tb + q_f / h

    def _k_of(self, T: float) -> float:
        f = self.result.inputs.get("_k_callable")
        return float(f(T)) if f is not None else float("nan")

    def wetted_top_profile(self, *, x: float) -> tuple[Any, Any]:
        """Lateral profile of the wetted TOP SURFACE at one axial station."""
        ys, ts = [], []
        for j in range(len(self.yf) - 1):
            yc = 0.5 * (self.yf[j] + self.yf[j + 1])
            if yc <= self.section.half_width_m:
                t = self.wetted_surface_temperature(x=x, y=yc)
                if t == t:
                    ys.append(yc)
                    ts.append(t)
        return _np.array(ys), _np.array(ts)


def solve_conduction(
    *, section: ChannelSection, solid: SolidProperties,
    heat_flux_W_m2: Callable[[float], float] | float,
    film_coefficient_W_m2K: Callable[[float], float] | float,
    bulk_temperature_K: Callable[[float], float] | float,
    axial_length_m: float = 0.0,
    n_axial: int = 1,
    n_lateral: tuple[int, int] = (12, 8),
    n_vertical: tuple[int, int, int] = (10, 8, 10),
    volumetric_source_W_m3: float = 0.0,
    heating_method: HeatingMethod = HeatingMethod.NOT_DECLARED,
    anisotropy_ratio: float | None = None,
    max_picard: int = 40, tolerance_K: float = 1.0e-8,
    back_face_h_W_m2K: float = 0.0,
    back_face_T_K: float = 293.15,
    sensors: Sequence["SensorInclusion"] = (),
) -> ConductionSolution:
    r"""Solve :math:`\nabla\cdot[k(T)\nabla T] + q''' = 0` on the section.

    Set ``axial_length_m = 0`` (the default) for the **2-D cross-section**:
    one axial cell with adiabatic end faces, which is the exact
    two-dimensional limit rather than an approximation to it. Give a length
    and ``n_axial > 1`` for the **3-D segment**, where ``heat_flux_W_m2``
    and ``bulk_temperature_K`` may be callables of the axial coordinate so
    that metered segments, guard segments and the coolant enthalpy rise are
    all represented.

    Boundary conditions, all of them:

    ==========================  =====================================
    heated face (z = 0)         specified flux ``heat_flux_W_m2(x)``
    wetted channel faces        ``h (T_face - T_bulk(x))``
    y = 0, y = half pitch       adiabatic (lateral symmetry)
    x = 0, x = L                adiabatic (axial symmetry / 2-D limit)
    back face (z = z_max)       ``back_face_h (T - back_face_T)``,
                                default h = 0 i.e. perfectly insulated
    ==========================  =====================================

    The default perfectly-insulated back face is an **idealisation of the
    facility's vacuum jacket**, not a measurement, and it biases the solve
    toward *higher* wall temperatures. It is exposed as an argument so the
    real insulation can replace it, and a non-zero value is reported.
    """
    inputs: dict[str, Any] = {
        "material": solid.name,
        "h_c_W_m2K": ("AXIALLY VARYING" if callable(film_coefficient_W_m2K)
                      else film_coefficient_W_m2K),
        "axial_length_m": axial_length_m, "n_axial": n_axial,
        "n_lateral": n_lateral, "n_vertical": n_vertical,
        "heating_method": heating_method.value,
        "q_volumetric_W_m3": volumetric_source_W_m3,
        "anisotropy_ratio": anisotropy_ratio}

    def _fail(status: Status, note: str) -> ConductionSolution:
        r = Result(float("nan"), "K", "CAN-WALL-SOLVE", status=status,
                   inputs=inputs, notes=note)
        return ConductionSolution(None, None, None, None, None, section,
                                  float("nan"), float("nan"), float("nan"),
                                  0, False, r)

    if not NUMERICS_AVAILABLE:                              # pragma: no cover
        return _fail(Status.INSUFFICIENT_EVIDENCE,
                     "numpy/scipy unavailable; a conduction bias cannot be "
                     "estimated without solving the field and is NOT assumed")
    bad = section.validate()
    if bad:
        return _fail(Status.PHYSICALLY_INVALID, "; ".join(bad))
    if not callable(film_coefficient_W_m2K) and film_coefficient_W_m2K <= 0:
        return _fail(Status.PHYSICALLY_INVALID,
                     "the coolant-side film coefficient must be > 0. A zero "
                     "means no heat leaves the specimen and the steady "
                     "problem has no solution.")
    if axial_length_m < 0 or n_axial < 1:
        return _fail(Status.PHYSICALLY_INVALID,
                     "axial length must be >= 0 and n_axial >= 1")

    statuses: list[Status] = [Status.PASS]
    notes: list[str] = []
    assumptions: list[Assumption] = []

    # --- the q''' = 0 question, asked out loud --------------------------
    if heating_method is HeatingMethod.DIRECT_JOULE and \
            volumetric_source_W_m3 == 0.0:
        return _fail(Status.PHYSICALLY_INVALID,
                     "DIRECT_JOULE heating with a zero volumetric source is "
                     "not an approximation, it is wrong: q''' = rho_e J^2 is "
                     "the entire heat input in that configuration.")
    if heating_method is HeatingMethod.NOT_DECLARED:
        statuses.append(Status.UNVALIDATED_ASSUMPTION)
        assumptions.append(Assumption(
            "heating_method", "NOT_DECLARED -> q''' = 0 assumed", "-",
            "q''' = 0 is EXACT for a separate surface heater and WRONG for "
            "direct Joule heating of the specimen. The facility has not "
            "declared which it will use, and external heating is not "
            "assumed silently.", Status.UNVALIDATED_ASSUMPTION))

    # --- anisotropy -----------------------------------------------------
    kz_over_kx = 1.0
    if anisotropy_ratio is not None:
        if anisotropy_ratio <= 0:
            return _fail(Status.PHYSICALLY_INVALID,
                         "the anisotropy ratio must be > 0")
        kz_over_kx = anisotropy_ratio
        statuses.append(Status.CONFLICT_UNRESOLVED)
        assumptions.append(Assumption(
            "anisotropy_ratio", anisotropy_ratio, "k_axial/k_normal",
            "no anisotropy of thermal conductivity has EVER been measured "
            "for GRCop in any form. The nearest alloy, LPBF CuCrZr, has two "
            "sources that disagree on the SIGN of the effect, and LPBF pure "
            "copper is reported near-isotropic. Any ratio applied here is a "
            "SENSITIVITY SWEEP, never a property.",
            Status.CONFLICT_UNRESOLVED))
    elif solid.anisotropy_ratio is None:
        assumptions.append(Assumption(
            "anisotropy", "ISOTROPIC assumed", "-",
            "the material's anisotropy_ratio is None, meaning UNKNOWN. An "
            "isotropic solve is used and the assumption is recorded. See "
            "ANISOTROPY_EVIDENCE for what would close it.",
            Status.UNVALIDATED_ASSUMPTION))
        statuses.append(Status.UNVALIDATED_ASSUMPTION)

    # --- grids ----------------------------------------------------------
    s = section
    ny1, ny2 = n_lateral
    nz1, nz2, nz3 = n_vertical
    y_extra: list[float] = []
    z_extra: list[float] = []
    for sn in sensors:
        z_extra.extend(sn.z_bounds(s.wall_thickness_m))
        y_extra.extend(sn.y_bounds())
    yf = _grid([0.0, s.half_width_m, s.half_pitch_m], [ny1, ny2], y_extra)
    zf = _grid([0.0, s.wall_thickness_m,
                s.wall_thickness_m + s.channel_height_m,
                s.total_height_m], [nz1, nz2, nz3], z_extra)
    if axial_length_m > 0.0:
        xf = _np.linspace(0.0, axial_length_m, n_axial + 1)
    else:
        xf = _np.array([0.0, 1.0])         # unit depth, 2-D limit
        n_axial = 1
    nx, ny, nz = len(xf) - 1, len(yf) - 1, len(zf) - 1
    dx = _np.diff(xf)
    dy = _np.diff(yf)
    dz = _np.diff(zf)
    xc = 0.5 * (xf[:-1] + xf[1:])
    yc = 0.5 * (yf[:-1] + yf[1:])
    zc = 0.5 * (zf[:-1] + zf[1:])

    # --- fluid mask (the channel), exact because the grid aligns --------
    solid_mask = _np.ones((nx, ny, nz), dtype=bool)
    in_ch_y = yc < s.half_width_m
    in_ch_z = (zc > s.wall_thickness_m) & (
        zc < s.wall_thickness_m + s.channel_height_m)
    solid_mask[:, _np.ix_(in_ch_y, in_ch_z)[0][:, None],
               _np.ix_(in_ch_y, in_ch_z)[1][None, :]] = False

    # --- sensor inclusions: cells that are solid but NOT the wall --------
    sensor_id = -_np.ones((nx, ny, nz), dtype=int)
    for si, sn in enumerate(sensors):
        zlo, zhi = sn.z_bounds(s.wall_thickness_m)
        ylo, yhi = sn.y_bounds()
        for j in range(ny):
            if not (ylo - 1e-15 <= yc[j] <= yhi + 1e-15):
                continue
            for k in range(nz):
                if zlo - 1e-15 <= zc[k] <= zhi + 1e-15:
                    sensor_id[:, j, k] = si
    for sn in sensors:
        zlo, zhi = sn.z_bounds(s.wall_thickness_m)
        if zlo <= 1e-12:
            return _fail(
                Status.PHYSICALLY_INVALID,
                f"sensor {sn.name!r} of diameter {sn.diameter_m * 1e3:.3f} mm "
                f"at depth {sn.depth_m * 1e3:.3f} mm BREACHES THE HEATED "
                f"FACE of a {s.wall_thickness_m * 1e3:.3f} mm wall (its body "
                f"reaches z = {zlo * 1e3:.3f} mm). A keep-out measured to the "
                f"JUNCTION is not a keep-out: the sensor has a diameter. "
                f"Deepest usable junction is "
                f"{(s.wall_thickness_m - 0.5 * sn.diameter_m) * 1e3:.3f} mm.")
        if zhi > s.wall_thickness_m + 1e-15:
            return _fail(
                Status.PHYSICALLY_INVALID,
                f"sensor {sn.name!r} BREAKS THROUGH INTO THE COOLANT CHANNEL "
                f"(its body reaches z = {zhi * 1e3:.3f} mm against a wall of "
                f"{s.wall_thickness_m * 1e3:.3f} mm). Shallowest usable "
                f"junction is {0.5 * sn.diameter_m * 1e3:.3f} mm.")
    if len(sensors) and _np.all(sensor_id < 0):
        return _fail(Status.PHYSICALLY_INVALID,
                     "no grid cell fell inside a declared sensor: the "
                     "sensor is smaller than one cell, so it is not being "
                     "modelled at all. Refine the grid or enlarge the "
                     "sensor -- do NOT proceed with an unrepresented sensor.")
    if sensors:
        # NUMERICAL HYGIENE. Aligning the grid to a very small sensor
        # creates a cell orders of magnitude thinner than its neighbours.
        # The solve would still converge and still conserve energy, and it
        # would still be untrustworthy, so it is refused rather than
        # reported with a clean-looking residual.
        for lbl, dd in (("lateral", dy), ("through-thickness", dz)):
            ratio = float(_np.max(dd) / max(float(_np.min(dd)), 1e-30))
            if ratio > 1.0e3:
                return _fail(
                    Status.PHYSICALLY_INVALID,
                    f"aligning the grid to the declared sensors forces a "
                    f"{lbl} cell-size ratio of {ratio:.3g}. A grid that "
                    f"anisotropic gives a converged, energy-conserving and "
                    f"untrustworthy field. Use a larger sensor or a coarser "
                    f"alignment -- do not read this solve.")

    idx = -_np.ones((nx, ny, nz), dtype=int)
    n_unk = int(solid_mask.sum())
    idx[solid_mask] = _np.arange(n_unk)

    qflux = (heat_flux_W_m2 if callable(heat_flux_W_m2)
             else (lambda _x, _q=float(heat_flux_W_m2): _q))
    tbulk = (bulk_temperature_K if callable(bulk_temperature_K)
             else (lambda _x, _t=float(bulk_temperature_K): _t))
    # The film coefficient may vary axially. This is what makes an imposed
    # h_c perturbation -- the thing HTD actually is -- representable
    # EXACTLY, rather than emulated through an equivalent bulk temperature.
    hfilm = (film_coefficient_W_m2K if callable(film_coefficient_W_m2K)
             else (lambda _x, _h=float(film_coefficient_W_m2K): _h))
    if callable(film_coefficient_W_m2K):
        _hs = [float(hfilm(float(v))) for v in xc]
        if min(_hs) <= 0.0:
            return _fail(Status.PHYSICALLY_INVALID,
                         f"the axially varying film coefficient reaches "
                         f"{min(_hs):.6g} W/(m^2.K) at some station. It must "
                         f"be > 0 everywhere.")

    # --- Picard ---------------------------------------------------------
    T = _np.full((nx, ny, nz), float(tbulk(xc[0])) + 50.0)
    T[~solid_mask] = _np.nan
    prev = None
    it = 0
    converged = False
    out_of_range = False
    for it in range(1, max_picard + 1):
        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []
        rhs = _np.zeros(n_unk)

        def kcell(i, j, k):
            nonlocal out_of_range
            si = int(sensor_id[i, j, k])
            if si >= 0:
                # the sensor is its OWN material. Using the wall's k here
                # would silently model a hole full of GRCop.
                return max(sensors[si].k_W_mK, 1e-30)
            kk, inside = solid.k_at(float(T[i, j, k]))
            if not inside:
                out_of_range = True
            return kk

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if not solid_mask[i, j, k]:
                        continue
                    p = idx[i, j, k]
                    ap = 0.0
                    kP = kcell(i, j, k)
                    Ax = dy[j] * dz[k]
                    Ay = dx[i] * dz[k]
                    Az = dx[i] * dy[j]

                    for (di, dj, dk, A, d_here, d_there, kmult) in (
                            (-1, 0, 0, Ax, dx[i], dx[i - 1] if i > 0 else 0.0,
                             kz_over_kx),
                            (+1, 0, 0, Ax, dx[i],
                             dx[i + 1] if i < nx - 1 else 0.0, kz_over_kx),
                            (0, -1, 0, Ay, dy[j], dy[j - 1] if j > 0 else 0.0,
                             1.0),
                            (0, +1, 0, Ay, dy[j],
                             dy[j + 1] if j < ny - 1 else 0.0, 1.0),
                            (0, 0, -1, Az, dz[k], dz[k - 1] if k > 0 else 0.0,
                             1.0),
                            (0, 0, +1, Az, dz[k],
                             dz[k + 1] if k < nz - 1 else 0.0, 1.0)):
                        ii, jj, kk_ = i + di, j + dj, k + dk
                        inside_dom = (0 <= ii < nx and 0 <= jj < ny
                                      and 0 <= kk_ < nz)
                        if inside_dom and solid_mask[ii, jj, kk_]:
                            kN = kcell(ii, jj, kk_)
                            # harmonic mean weighted by the two half-widths
                            r_series = (d_here / max(kP, 1e-30)
                                        + d_there / max(kN, 1e-30))
                            # a SENSOR/WALL face carries a contact resistance
                            # in series. It is never omitted: an omitted
                            # contact resistance is an assumed perfect bond.
                            sa = int(sensor_id[i, j, k])
                            sb = int(sensor_id[ii, jj, kk_])
                            if sa != sb:
                                sn_ = sensors[sa if sa >= 0 else sb]
                                hc_ = sn_.contact_conductance_W_m2K
                                if hc_ is None:
                                    return _fail(
                                        Status.INSUFFICIENT_EVIDENCE,
                                        f"sensor {sn_.name!r} declares no "
                                        f"contact conductance. A missing "
                                        f"interface resistance is an ASSUMED "
                                        f"PERFECT BOND, and no published "
                                        f"measured value exists for an "
                                        f"embedded thermocouple in metal by "
                                        f"any installation method. Supply a "
                                        f"value to sweep, or do not model "
                                        f"the sensor.")
                                r_series += (d_here + d_there) / max(
                                    hc_ * (d_here + d_there), 1e-30)
                            kf = (d_here + d_there) / r_series
                            g = kmult * kf * A / (0.5 * (d_here + d_there))
                            ap += g
                            rows.append(p)
                            cols.append(idx[ii, jj, kk_])
                            vals.append(-g)
                            continue
                        # --- a boundary or a fluid face ------------------
                        if dk == -1 and k == 0:
                            rhs[p] += qflux(float(xc[i])) * A
                            continue
                        if dk == +1 and k == nz - 1:
                            if back_face_h_W_m2K > 0.0:
                                gg = 1.0 / (1.0 / back_face_h_W_m2K
                                            + 0.5 * d_here / max(kP, 1e-30))
                                ap += gg * A
                                rhs[p] += gg * A * back_face_T_K
                            continue
                        if di != 0 and (ii < 0 or ii >= nx):
                            continue                    # axial symmetry
                        if dj != 0 and (jj < 0 or jj >= ny):
                            continue                    # lateral symmetry
                        # remaining case: face onto the coolant
                        gg = 1.0 / (1.0 / hfilm(float(xc[i]))
                                    + 0.5 * d_here / max(kP, 1e-30))
                        ap += gg * A
                        rhs[p] += gg * A * tbulk(float(xc[i]))

                    # --- lead conduction out of a sensor cell -----------
                    si_here = int(sensor_id[i, j, k])
                    if si_here >= 0:
                        sn_ = sensors[si_here]
                        if sn_.lead_conductance_W_K > 0.0:
                            if sn_.lead_exit_temperature_K is None:
                                return _fail(
                                    Status.INSUFFICIENT_EVIDENCE,
                                    f"sensor {sn_.name!r} has a lead "
                                    f"conductance but no lead exit "
                                    f"temperature. The fin drains heat TO "
                                    f"somewhere and that somewhere is not "
                                    f"assumed.")
                            n_cells_s = int((sensor_id == si_here).sum())
                            gl = sn_.lead_conductance_W_K / max(n_cells_s, 1)
                            ap += gl
                            rhs[p] += gl * sn_.lead_exit_temperature_K
                    if volumetric_source_W_m3:
                        rhs[p] += volumetric_source_W_m3 * Ax * dx[i]
                    rows.append(p)
                    cols.append(p)
                    vals.append(ap)

        A_mat = _csr((vals, (rows, cols)), shape=(n_unk, n_unk))
        sol = _spsolve(A_mat.tocsc(), rhs)
        newT = _np.full((nx, ny, nz), _np.nan)
        newT[solid_mask] = sol
        if prev is not None:
            delta = float(_np.nanmax(_np.abs(newT - T)))
            T = newT
            if delta < tolerance_K:
                converged = True
                break
        else:
            T = newT
        prev = True
    else:
        converged = False

    # --- energy balance and the wetted-face split -----------------------
    q_in = 0.0
    for i in range(nx):
        for j in range(ny):
            if solid_mask[i, j, 0]:
                q_in += qflux(float(xc[i])) * dx[i] * dy[j]
    # VOLUMETRIC GENERATION IS HEAT IN. Omitting it from the balance would
    # make the residual meaningless for exactly the architecture -- direct
    # Joule or induction heating -- where the balance matters most.
    if volumetric_source_W_m3:
        vol = 0.0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if solid_mask[i, j, k]:
                        vol += dx[i] * dy[j] * dz[k]
        q_in += volumetric_source_W_m3 * vol
    faces = {"TOP": 0.0, "SIDE": 0.0, "BOTTOM": 0.0, "BACK": 0.0}
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not solid_mask[i, j, k]:
                    continue
                kP, _ = solid.k_at(float(T[i, j, k]))
                for (di, dj, dk, A, d_here, label) in (
                        (0, 0, +1, dx[i] * dy[j], dz[k], "TOP"),
                        (0, 0, -1, dx[i] * dy[j], dz[k], "BOTTOM"),
                        (0, +1, 0, dx[i] * dz[k], dy[j], "SIDE"),
                        (0, -1, 0, dx[i] * dz[k], dy[j], "SIDE")):
                    ii, jj, kk_ = i + di, j + dj, k + dk
                    if not (0 <= ii < nx and 0 <= jj < ny and 0 <= kk_ < nz):
                        continue
                    if solid_mask[ii, jj, kk_]:
                        continue
                    gg = 1.0 / (1.0 / hfilm(float(xc[i]))
                                + 0.5 * d_here / max(kP, 1e-30))
                    faces[label] += gg * A * (float(T[i, j, k])
                                              - tbulk(float(xc[i])))
    if back_face_h_W_m2K > 0.0:
        for i in range(nx):
            for j in range(ny):
                k = nz - 1
                if solid_mask[i, j, k]:
                    kP, _ = solid.k_at(float(T[i, j, k]))
                    gg = 1.0 / (1.0 / back_face_h_W_m2K
                                + 0.5 * dz[k] / max(kP, 1e-30))
                    faces["BACK"] += gg * dx[i] * dy[j] * (
                        float(T[i, j, k]) - back_face_T_K)
    q_out = sum(faces.values())
    resid = abs(q_in - q_out) / max(abs(q_in), 1e-30)

    if not converged:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append(f"Picard did not converge in {max_picard} iterations; "
                     f"the field is NOT returned as an answer")
    if resid > 1e-8:
        statuses.append(Status.INSUFFICIENT_EVIDENCE)
        notes.append(f"energy balance residual {resid:.3e} exceeds 1e-8 -- "
                     f"the discrete problem does not conserve energy and no "
                     f"bias computed from it can be trusted")
    if out_of_range:
        statuses.append(Status.OUT_OF_CORRELATION_RANGE)
        notes.append(f"the field reached temperatures outside "
                     f"{solid.name}'s stated {solid.valid_range_K[0]:.0f}-"
                     f"{solid.valid_range_K[1]:.0f} K range; k was HELD at "
                     f"the range edge and NEVER extrapolated")
    notes.append(f"Q_in {q_in:.4f} W, Q_out {q_out:.4f} W, residual "
                 f"{resid:.2e}; wetted split TOP {faces['TOP'] / max(q_out, 1e-30):.1%} "
                 f"SIDE {faces['SIDE'] / max(q_out, 1e-30):.1%} BOTTOM "
                 f"{faces['BOTTOM'] / max(q_out, 1e-30):.1%}")

    res = Result(
        float(_np.nanmax(T)), "K", "CAN-WALL-SOLVE",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="steady conduction div[k(T) grad T] + q''' = 0; cell-centred "
               "finite volume with harmonic-mean face conductivity",
        source_locator="see docs/research/"
                       "WALL_CONDUCTION_PRIMARY_SOURCE_AUDIT.md; k(T) from "
                       + solid.source_locator,
        validity_domain="steady, opaque solid, no radiation inside the "
                        "channel, laterally periodic array of identical "
                        "channels; k(T) domain " + str(solid.valid_range_K),
        uncertainty=f"discrete energy balance closes to {resid:.1e}; the "
                    f"harmonic-mean order of accuracy on a BLOCK-NON-UNIFORM "
                    f"grid is VALIDITY_DOMAIN_UNKNOWN (every rigorous "
                    f"analysis opened assumes a uniform grid), so grid "
                    f"convergence is demonstrated rather than assumed",
        assumptions=tuple(assumptions), notes="; ".join(notes),
        inputs={**inputs, "n_cells": int(n_unk), "Q_in_W": q_in,
                "Q_out_W": q_out, "energy_residual": resid,
                "picard_iterations": it,
                # kept so the solution object can evaluate the wetted
                # SURFACE temperature from the face flux rather than from a
                # cell centre half a cell inside the metal
                "_bulk_callable": tbulk, "_h_callable": hfilm,
                "_k_callable": lambda T: solid.k_at(T)[0]})
    sol = ConductionSolution(T, xf, yf, zf, solid_mask, section, q_in, q_out,
                             resid, it, converged, res, faces)
    sol.sensor_id = sensor_id
    sol.sensors = tuple(sensors)
    return sol


def cross_section_solution(**kw) -> ConductionSolution:
    """The 2-D limit: one axial cell, adiabatic ends. Exact, not approximate."""
    kw.pop("axial_length_m", None)
    kw.pop("n_axial", None)
    return solve_conduction(axial_length_m=0.0, n_axial=1, **kw)


def segment_solution(*, axial_length_m: float, n_axial: int = 24,
                     **kw) -> ConductionSolution:
    """The 3-D case: one metered segment plus its neighbouring guards."""
    return solve_conduction(axial_length_m=axial_length_m, n_axial=n_axial,
                            **kw)


# ===========================================================================
# 4. WHAT THE FIELD SAYS ABOUT THE MEASUREMENT
# ===========================================================================

def sensor_reading(sol: ConductionSolution, index: int = 0) -> Result:
    r"""What the thermocouple actually reads, and where that came from.

    The junction is a point, not the whole inclusion, so this samples the
    cell containing the junction rather than averaging the sensor volume.
    Averaging would quietly assume an isothermal bead, which for a 0.5 mm
    inclusion in a 150 K/mm gradient is a 75 K assumption.
    """
    if sol.T is None or sol.sensor_id is None or index >= len(sol.sensors):
        return Result(float("nan"), "K", "CAN-WALL-SENSORREAD",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="no sensor of that index was modelled")
    sn = sol.sensors[index]
    zc_j = (sol.section.wall_thickness_m - sn.depth_m
            - sn.junction_offset_m)
    i = 0 if len(sol.xf) == 2 else (len(sol.xf) - 1) // 2
    j = sol._cell(sol.yf, sn.lateral_position_m)
    k = sol._cell(sol.zf, zc_j)
    if int(sol.sensor_id[i, j, k]) != index:
        return Result(
            float("nan"), "K", "CAN-WALL-SENSORREAD",
            status=Status.PHYSICALLY_INVALID,
            notes=f"the junction of {sn.name!r} at offset "
                  f"{sn.junction_offset_m * 1e6:.0f} um falls OUTSIDE its "
                  f"own inclusion. A junction cannot sit in the metal and "
                  f"in the sensor at once.")
    return Result(
        float(sol.T[i, j, k]), "K", "CAN-WALL-SENSORREAD",
        status=Status.PASS, evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="the junction cell of the modelled sensor inclusion",
        source_locator="-",
        validity_domain="the solved section",
        uncertainty="inherits the field's energy residual",
        notes=f"{sn.name}: junction at depth "
              f"{(sn.depth_m + sn.junction_offset_m) * 1e3:.3f} mm reads "
              f"{float(sol.T[i, j, k]):.3f} K",
        inputs={"sensor": sn.name, "depth_m": sn.depth_m,
                "junction_offset_m": sn.junction_offset_m,
                "k_sensor_W_mK": sn.k_W_mK,
                "contact_conductance_W_m2K": sn.contact_conductance_W_m2K,
                "lead_conductance_W_K": sn.lead_conductance_W_K})


def wetted_face_heat_split(sol: ConductionSolution) -> Result:
    r"""How much of the metered power actually leaves through the TOP wall.

    **This is the bias the 1-D model cannot even represent.** The facility
    defines :math:`q'' = P/(w L_{seg})` -- power over the *heated width* --
    and pairs it with a top-wall temperature. But the rib is a fin: heat
    conducts laterally and leaves through the channel's **side** and
    **bottom** walls as well. If a fraction :math:`\phi` of the power leaves
    the top, then the true top-wall flux is :math:`\phi P/(wL)` and using
    the full :math:`P/(wL)` overstates :math:`h_c` by :math:`1/\phi`.

    That error is **multiplicative and systematic**. It does not shrink with
    better instruments and it is invisible to any 1-D reconstruction.
    """
    if sol.T is None:
        return Result(float("nan"), "-", "CAN-WALL-SPLIT",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="no field to split")
    total = sum(sol.face_power_W.values())
    frac = {k: v / total for k, v in sol.face_power_W.items() if total}
    phi = frac.get("TOP", float("nan"))
    return Result(
        phi, "-", "CAN-WALL-SPLIT",
        status=Status.PASS, evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="integral of the convective flux over each wetted face of "
               "the converged conduction field",
        source_locator="-",
        validity_domain="the solved section only",
        uncertainty=f"inherits the field's energy residual "
                    f"{sol.energy_residual:.1e}",
        notes=f"TOP {frac.get('TOP', 0):.1%}, SIDE {frac.get('SIDE', 0):.1%}, "
              f"BOTTOM {frac.get('BOTTOM', 0):.1%}, BACK "
              f"{frac.get('BACK', 0):.1%}. Using q'' = P/(w L) with a "
              f"top-wall temperature overstates h_c by a factor "
              f"{1.0 / max(phi, 1e-30):.3f} if phi is not applied.",
        inputs={"fractions": frac, "face_power_W": sol.face_power_W})


def reconstruction_bias(
    sol: ConductionSolution, *, sensor_depths_m: tuple[float, float],
    lateral_position_m: float = 0.0, axial_position_m: float | None = None,
) -> Result:
    r"""Bias of the 1-D two-sensor reconstruction against the true field.

    Reports **two separate numbers, never summed**:

    ``bias_reconstruction_K``
        the 1-D formula applied to the field's OWN sensor readings, minus
        the field's wetted-wall temperature *at the same lateral station*.
        This is the error of the formula.

    ``bias_definition_K``
        that same local wetted-wall temperature, minus the **area-averaged**
        wetted top-wall temperature that :math:`q''=P/(wL)` is actually
        paired with. This is the error of the *definition* and survives even
        a perfect formula.

    A bias is not a random error. They are not combined here, and
    `allowable_modelling_bias` treats them as the systematic terms they are.
    """
    if sol.T is None:
        return Result(float("nan"), "K", "CAN-WALL-BIAS",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="no field to compare against")
    x1, x2 = sensor_depths_m
    s = sol.section
    inputs = {"x1_m": x1, "x2_m": x2, "y_m": lateral_position_m}
    if not (0.0 < x1 < x2):
        return Result(float("nan"), "K", "CAN-WALL-BIAS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="depths must satisfy 0 < x1 < x2; they are "
                            "measured from the WETTED surface into the solid")
    if x2 >= s.wall_thickness_m:
        return Result(
            float("nan"), "K", "CAN-WALL-BIAS",
            status=Status.PHYSICALLY_INVALID, inputs=inputs,
            notes=f"the deeper sensor at {x2 * 1e3:.3f} mm is at or beyond "
                  f"the heated face ({s.wall_thickness_m * 1e3:.3f} mm). A "
                  f"sensor cannot sit in the heater. Either thicken the wall "
                  f"or move the sensor.")
    if lateral_position_m > s.half_width_m:
        return Result(float("nan"), "K", "CAN-WALL-BIAS",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the lateral station must lie over the channel, "
                            "not over the rib")

    xa = (axial_position_m if axial_position_m is not None
          else 0.5 * (sol.xf[0] + sol.xf[-1]))
    z1 = s.wall_thickness_m - x1
    z2 = s.wall_thickness_m - x2
    t1 = sol.at(x=xa, y=lateral_position_m, z=z1)
    t2 = sol.at(x=xa, y=lateral_position_m, z=z2)
    g = (t2 - t1) / (x2 - x1)
    t_1d = t1 - g * x1

    t_true_local = sol.wetted_surface_temperature(x=xa,
                                                  y=lateral_position_m)
    ys, ts = sol.wetted_top_profile(x=xa)
    t_area_avg = float(_np.mean(ts)) if len(ts) else float("nan")

    bias_rec = t_1d - t_true_local
    bias_def = t_true_local - t_area_avg
    lateral_spread = (float(_np.max(ts) - _np.min(ts)) if len(ts)
                      else float("nan"))

    statuses = [Status.PASS, Status.UNVALIDATED_ASSUMPTION]
    return Result(
        bias_rec, "K", "CAN-WALL-BIAS",
        status=worst(statuses), all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="the 1-D reconstruction of facility."
               "wall_temperature_from_embedded_sensors applied to a solved "
               "2-D/3-D conduction field",
        source_locator="-",
        validity_domain="the solved section and the stated sensor station",
        uncertainty="a MODELLING BIAS, not a random error; it must not be "
                    "combined in quadrature with instrument uncertainties",
        assumptions=(Assumption(
            "bias_is_systematic", bias_rec, "K",
            "this is a deterministic offset of the reconstruction, not a "
            "scatter. Adding it in quadrature with thermometry would "
            "understate it whenever it dominates and misrepresent it "
            "always.", Status.UNVALIDATED_ASSUMPTION),),
        notes=f"T_1D {t_1d:.2f} K, T_true(local) {t_true_local:.2f} K, "
              f"T_true(area avg) {t_area_avg:.2f} K; reconstruction bias "
              f"{bias_rec:+.2f} K, definition bias {bias_def:+.2f} K; "
              f"lateral spread across the wetted top wall "
              f"{lateral_spread:.2f} K",
        inputs={**inputs, "T_sensor_1_K": t1, "T_sensor_2_K": t2,
                "gradient_K_per_m": g, "T_1D_K": t_1d,
                "T_true_local_K": t_true_local,
                "T_true_area_avg_K": t_area_avg,
                "bias_reconstruction_K": bias_rec,
                "bias_definition_K": bias_def,
                "lateral_spread_K": lateral_spread,
                "implied_flux_W_m2": abs(g) * sol.solid_k_hint
                if hasattr(sol, "solid_k_hint") else float("nan")})


def axial_leak_from_field(
    sol: ConductionSolution, *, segment_start_m: float, segment_end_m: float,
    solid: SolidProperties, segment_power_W: float) -> Result:
    r"""Inter-segment conduction leak measured on the ACTUAL field.

    The facility spec bounded this with a **bulk-gradient proxy**: it took
    the axial gradient of the coolant bulk temperature and multiplied by the
    solid cross-section. That is a lower bound, because the wall's axial
    gradient is steeper than the bulk's wherever the film coefficient
    varies.

    This function replaces the proxy with the integral of
    :math:`-k\,\partial T/\partial x` over the actual solid cross-section at
    the segment boundary, and reports the ratio between them as
    ``guard_bias_factor``. **That factor is measured on the field, not
    fitted.**
    """
    if sol.T is None or sol.xf is None or len(sol.xf) < 3:
        return Result(float("nan"), "-", "CAN-WALL-AXLEAK",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="an axial leak needs a 3-D field with more than "
                            "one axial cell; the 2-D limit cannot show it")
    if segment_power_W <= 0:
        return Result(float("nan"), "-", "CAN-WALL-AXLEAK",
                      status=Status.PHYSICALLY_INVALID,
                      notes="segment power must be > 0")

    # THE LEAK MUST BE EVALUATED WITH THE AXIAL CONDUCTIVITY. If the field
    # was solved anisotropically, integrating it with k_normal would report
    # a leak the solve never produced -- and would get the SIGN of the
    # anisotropy effect backwards, because a lower k_axial sharpens the
    # boundary gradient at the same time as it reduces the conductance.
    ani = sol.result.inputs.get("anisotropy_ratio")
    k_axial_mult = float(ani) if ani else 1.0

    def _flux_at(xb: float) -> float:
        i = int(_np.searchsorted(sol.xf, xb) - 1)
        i = min(max(i, 0), len(sol.xf) - 3)
        q = 0.0
        dy = _np.diff(sol.yf)
        dz = _np.diff(sol.zf)
        dxa = sol.xf[i + 1] - sol.xf[i]
        dxb = sol.xf[i + 2] - sol.xf[i + 1]
        for j in range(len(dy)):
            for k in range(len(dz)):
                if not (sol.solid[i, j, k] and sol.solid[i + 1, j, k]):
                    continue
                kA, _ = solid.k_at(float(sol.T[i, j, k]))
                kB, _ = solid.k_at(float(sol.T[i + 1, j, k]))
                kf = k_axial_mult * (dxa + dxb) / (dxa / kA + dxb / kB)
                grad = ((float(sol.T[i + 1, j, k]) - float(sol.T[i, j, k]))
                        / (0.5 * (dxa + dxb)))
                q += -kf * grad * dy[j] * dz[k]
        return q

    q_lo = _flux_at(segment_start_m)
    q_hi = _flux_at(segment_end_m)
    net = abs(q_lo) + abs(q_hi)
    frac = net / segment_power_W
    return Result(
        frac, "-", "CAN-WALL-AXLEAK", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="Fourier conduction integrated over the solid cross-section "
               "at each segment boundary of the converged 3-D field",
        source_locator="-",
        validity_domain="the solved segment stack",
        uncertainty=f"inherits the field's energy residual "
                    f"{sol.energy_residual:.1e}",
        notes=f"leak in {q_lo:+.3f} W, out {q_hi:+.3f} W, total magnitude "
              f"{net:.3f} W against {segment_power_W:.1f} W of segment "
              f"power = {frac:.3%}",
        inputs={"Q_in_boundary_W": q_lo, "Q_out_boundary_W": q_hi,
                "net_magnitude_W": net, "P_segment_W": segment_power_W,
                "leak_fraction": frac})


def guard_matching_requirement(
    *, leak_fraction_matched: float, leak_fraction_at_mismatch: float,
    mismatch_used: float, guard_budget_rel: float) -> Result:
    r"""How closely must the guard power track the metered segment?

    The facility spec bounded the inter-segment leak with the **coolant bulk
    axial gradient** and derived a minimum segment length from it. The 3-D
    field shows that bound is reached **only when the guards are perfect**.
    A guard that is 10 % low pushes the leak far past it, and no segment
    length recovers that -- the leak scales with the *mismatch*, not with
    the geometry.

    Linear in the mismatch to first order, so

    .. math::
        m_{max} = \frac{f_{budget} - f_0}{(f_m - f_0)/m}

    is the tolerance on guard power. `PHYSICALLY_INVALID` when even a
    perfectly matched guard already exceeds the budget, because then the
    budget itself is unachievable and the geometry must change.
    """
    inputs = {"f0": leak_fraction_matched, "fm": leak_fraction_at_mismatch,
              "m": mismatch_used, "budget": guard_budget_rel}
    if mismatch_used <= 0 or guard_budget_rel <= 0:
        return Result(float("nan"), "-", "CAN-WALL-GUARD",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the mismatch probe and the budget must be > 0")
    if leak_fraction_at_mismatch <= leak_fraction_matched:
        return Result(float("nan"), "-", "CAN-WALL-GUARD",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="a mismatched guard must leak MORE than a "
                            "matched one; these two numbers are inconsistent")
    if leak_fraction_matched >= guard_budget_rel:
        return Result(
            float("nan"), "-", "CAN-WALL-GUARD",
            status=Status.PHYSICALLY_INVALID, inputs=inputs,
            evidence_level=EvidenceLevel.E1,
            notes=f"even a PERFECTLY matched guard leaks "
                  f"{leak_fraction_matched:.3%}, at or beyond the "
                  f"{guard_budget_rel:.1%} budget. Guard control cannot fix "
                  f"this -- lengthen the segment, thin the solid path, or "
                  f"widen the budget explicitly.")
    slope = (leak_fraction_at_mismatch - leak_fraction_matched) / mismatch_used
    m_max = (guard_budget_rel - leak_fraction_matched) / slope
    return Result(
        m_max, "-", "CAN-WALL-GUARD", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="linear response of the 3-D inter-segment conduction leak to "
               "guard power mismatch, measured on the solved field",
        source_locator="-",
        validity_domain="small mismatches, where the response is linear; "
                        "re-measure for a different geometry",
        uncertainty="inherits the field's energy residual",
        notes=f"a matched guard leaks {leak_fraction_matched:.3%}; a "
              f"{mismatch_used:.0%} mismatch leaks "
              f"{leak_fraction_at_mismatch:.3%}, so the guard must track the "
              f"metered segment to within {m_max:.2%} to hold the "
              f"{guard_budget_rel:.1%} budget. THE BULK-GRADIENT BOUND IS "
              f"THE FLOOR, NOT THE REQUIREMENT: it is reached only by a "
              f"perfect guard, and guard control is what actually decides "
              f"the leak.",
        inputs={**inputs, "slope_per_unit_mismatch": slope,
                "max_mismatch": m_max})


def identify_h_c(
    *, section: ChannelSection, solid: SolidProperties,
    heat_flux_W_m2: float, bulk_temperature_K: float,
    measured_sensor_temperatures_K: tuple[float, float],
    sensor_depths_m: tuple[float, float],
    bracket_W_m2K: tuple[float, float] = (1.0e4, 1.0e6),
    tolerance_rel: float = 1.0e-4, max_iterations: int = 40,
    **solve_kw) -> Result:
    r"""Identify :math:`h_c` by matching the conduction field to the sensors.

    **This replaces the divide.** The obvious reduction,
    :math:`h_c = q''/(T_w-T_b)`, needs a local wetted-wall flux, and the
    2-D field says the two candidate ways of getting one are both badly
    biased:

    * **metered power over the channel width** overstates the local top-wall
      flux by roughly a factor of two, because the rib is a fin and only
      about half the power leaves through the top;
    * **the two-sensor gradient times k** overstates it by 13-27 %, because
      the flux at mid-wall depth has not yet finished spreading sideways.

    Both errors are *computable from the conduction model* -- but the model
    needs :math:`h_c`, which is the unknown. So the reduction is an
    **inverse problem**, not a division: find the :math:`h_c` whose field
    reproduces the measured sensor temperatures.

    It is well conditioned, and unusually so. The sensor temperature moves
    by hundreds of kelvin per e-fold of :math:`h_c`, so a kelvin of
    thermometric error is a fraction of a percent of :math:`h_c` -- and the
    :math:`\pm 4\%` uncertainty on :math:`k(T)` enters at about the same
    small level, because :math:`k` and :math:`h_c` move the field in
    similar ways with very different leverage.

    **It is not free.** The identified value inherits every assumption of
    the conduction model -- the planar-array symmetry, the insulated back
    face, isotropic :math:`k`, no contact resistance at the sensor -- so it
    carries `UNVALIDATED_ASSUMPTION`, and it is a *measurement reduction*,
    never a correlation. `select_coolant_correlation` is unaffected.
    """
    inputs = {"q_W_m2": heat_flux_W_m2, "T_bulk_K": bulk_temperature_K,
              "T_sensors_K": list(measured_sensor_temperatures_K),
              "x_m": list(sensor_depths_m), "bracket": list(bracket_W_m2K)}
    if not NUMERICS_AVAILABLE:                              # pragma: no cover
        return Result(float("nan"), "W/(m^2.K)", "CAN-WALL-IDENT",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes="numerics unavailable")
    x1, x2 = sensor_depths_m
    if not (0.0 < x1 < x2 < section.wall_thickness_m):
        return Result(float("nan"), "W/(m^2.K)", "CAN-WALL-IDENT",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="need 0 < x1 < x2 < wall thickness")
    t1m, t2m = measured_sensor_temperatures_K
    if t2m <= t1m:
        return Result(float("nan"), "W/(m^2.K)", "CAN-WALL-IDENT",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the deeper sensor must be hotter; heat flows "
                            "from the heated face toward the coolant")

    z1 = section.wall_thickness_m - x1
    z2 = section.wall_thickness_m - x2

    def residual(h: float) -> float:
        sol = solve_conduction(
            section=section, solid=solid, heat_flux_W_m2=heat_flux_W_m2,
            film_coefficient_W_m2K=h, bulk_temperature_K=bulk_temperature_K,
            **solve_kw)
        if sol.T is None:
            return float("nan")
        xa = 0.5 * (sol.xf[0] + sol.xf[-1])
        t1 = sol.at(x=xa, y=0.0, z=z1)
        t2 = sol.at(x=xa, y=0.0, z=z2)
        return 0.5 * ((t1 - t1m) + (t2 - t2m))

    lo, hi = bracket_W_m2K
    f_lo, f_hi = residual(lo), residual(hi)
    if not (f_lo == f_lo and f_hi == f_hi):                 # pragma: no cover
        return Result(float("nan"), "W/(m^2.K)", "CAN-WALL-IDENT",
                      status=Status.INSUFFICIENT_EVIDENCE, inputs=inputs,
                      notes="the conduction solve failed at a bracket end")
    if f_lo * f_hi > 0:
        return Result(
            float("nan"), "W/(m^2.K)", "CAN-WALL-IDENT",
            status=Status.PHYSICALLY_INVALID,
            inputs={**inputs, "residual_lo_K": f_lo, "residual_hi_K": f_hi},
            evidence_level=EvidenceLevel.E1,
            notes=f"the measured sensor temperatures are not reproducible by "
                  f"ANY h_c in {lo:.3g}-{hi:.3g} W/(m^2.K): the residual is "
                  f"{f_lo:+.1f} K at the low end and {f_hi:+.1f} K at the "
                  f"high end and does not change sign. Either the "
                  f"measurement or the conduction model is wrong, and this "
                  f"refuses rather than reporting a bracket end.")
    it = 0
    for it in range(1, max_iterations + 1):
        mid = math.sqrt(lo * hi)                # geometric -- h spans decades
        f_mid = residual(mid)
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
        if (hi - lo) / mid < tolerance_rel:
            break
    h_id = math.sqrt(lo * hi)

    # conditioning, measured rather than asserted
    d = 0.02
    r_p = residual(h_id * (1.0 + d))
    r_m = residual(h_id * (1.0 - d))
    dT_dlnh = (r_p - r_m) / (math.log(1.0 + d) - math.log(1.0 - d))
    per_K = abs(1.0 / dT_dlnh) if dT_dlnh else float("nan")

    statuses = [Status.PASS, Status.UNVALIDATED_ASSUMPTION,
                Status.INSUFFICIENT_EVIDENCE]
    return Result(
        h_id, "W/(m^2.K)", "CAN-WALL-IDENT", status=worst(statuses),
        all_statuses=frozenset(statuses),
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="inverse solution of the 2-D conduction field against the "
               "measured embedded-sensor temperatures",
        source_locator="-",
        validity_domain="the modelled section, at the modelled boundary "
                        "conditions; the identified value inherits EVERY "
                        "assumption of the conduction model",
        uncertainty=f"conditioning {abs(dT_dlnh):.0f} K per e-fold of h_c, "
                    f"so 1 K of sensor error is {per_K:.2%} of h_c",
        assumptions=(Assumption(
            "identification_inherits_the_model", "2-D planar-array section, "
            "insulated back face, isotropic k, no sensor contact "
            "resistance", "-",
            "an identified h_c is only as good as the field it was "
            "identified in. These four are the model's load-bearing "
            "assumptions and none of them has been measured.",
            Status.UNVALIDATED_ASSUMPTION),),
        notes=f"h_c = {h_id:.4g} W/(m^2.K) in {it} bisections; "
              f"{abs(dT_dlnh):.0f} K per e-fold means 1 K of sensor error is "
              f"{per_K:.2%} of h_c. THIS IS A MEASUREMENT REDUCTION, NOT A "
              f"CORRELATION: it says what h_c WAS in this test, never what "
              f"h_c WILL BE in a chamber, and select_coolant_correlation is "
              f"untouched by it.",
        inputs={**inputs, "h_c_W_m2K": h_id, "bisections": it,
                "dT_per_e_fold_K": dT_dlnh,
                "rel_h_c_per_K_of_sensor_error": per_K})


def allowable_modelling_bias(
    *, wall_temperature_target_K: float,
    random_uncertainty_K: float) -> Result:
    r"""How much systematic bias fits inside the wall-temperature target?

    **Bias and random uncertainty are not added in quadrature.** A bias is a
    deterministic offset; a random uncertainty is a distribution width. The
    defensible combination for a target that must be met *in the worst case*
    is the linear one:

    .. math:: |b| + k\,u \le U_{target}

    with :math:`k = 1` here because the project's stated targets (+/-25 K
    and +/-50 K) were themselves derived as standard uncertainties, not as
    expanded ones. Applying a coverage factor to one term and not the other
    would be worse than applying none.

    Returns the bias budget remaining. `PHYSICALLY_INVALID` when the random
    term alone already consumes the target -- because then no conduction
    model, however good, rescues the measurement.
    """
    inputs = {"target_K": wall_temperature_target_K,
              "random_K": random_uncertainty_K}
    if wall_temperature_target_K <= 0 or random_uncertainty_K < 0:
        return Result(float("nan"), "K", "CAN-WALL-BIASBUDGET",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="the target must be > 0 and the random term >= 0")
    budget = wall_temperature_target_K - random_uncertainty_K
    if budget <= 0:
        return Result(
            float("nan"), "K", "CAN-WALL-BIASBUDGET",
            status=Status.PHYSICALLY_INVALID, inputs=inputs,
            evidence_level=EvidenceLevel.E1,
            verification=Verification.ALGEBRAIC,
            source="linear combination of a systematic bias and a random "
                   "uncertainty",
            notes=f"the RANDOM term alone ({random_uncertainty_K:.2f} K) "
                  f"meets or exceeds the {wall_temperature_target_K:.1f} K "
                  f"target. No conduction model improves this -- the "
                  f"instrument does.")
    return Result(
        budget, "K", "CAN-WALL-BIASBUDGET", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="linear combination |b| + u <= U_target; a bias is not a "
               "distribution and is not root-sum-squared with one",
        source_locator="targets from docs/research/"
                       "COOLANT_EXPERIMENTAL_CLOSURE_PLAN.md section 1",
        validity_domain="worst-case interpretation of the target",
        uncertainty="exact given the inputs",
        notes=f"{budget:.2f} K of the {wall_temperature_target_K:.1f} K "
              f"target remains for MODELLING BIAS after "
              f"{random_uncertainty_K:.2f} K of random uncertainty",
        inputs={**inputs, "bias_budget_K": budget})


# ===========================================================================
# 5. SENSOR PLACEMENT
# ===========================================================================

def depth_error_sweep(
    *, heat_flux_W_m2: float, conductivity_W_mK: float,
    sensor_uncertainty_K: float, depths_m: tuple[float, float],
    depth_uncertainties_m: Sequence[float],
    depth_errors_correlated: bool = False) -> Result:
    r"""Thermometric versus geometric error, over a sweep of depth accuracy.

    With :math:`a = x_1/(x_2-x_1)` the two-sensor reconstruction gives

    .. math::
        \delta T_{therm} = u_T\sqrt{(1+a)^2+a^2},\qquad
        \delta T_{geom} = |g|\,\delta x\sqrt{1+a^2}

    for **independent** depth errors. If both sensors are placed by the same
    process -- one drill setup, one print, one CT calibration -- the errors
    are **correlated**, and a common-mode shift partly cancels because it
    moves both sensors together:

    .. math:: \delta T_{geom}^{corr} = |g|\,\delta x\,|1 - a + a| = |g|\,\delta x

    Both are reported. Assuming independence when the errors are common-mode
    **overstates** the geometric term; assuming correlation when they are
    independent **understates** it. The facility must say which it is, and
    until it does the independent (conservative) value governs.
    """
    x1, x2 = depths_m
    inputs = {"q_W_m2": heat_flux_W_m2, "k_W_mK": conductivity_W_mK,
              "x1_m": x1, "x2_m": x2, "u_T_K": sensor_uncertainty_K}
    if not (0.0 < x1 < x2) or conductivity_W_mK <= 0 or heat_flux_W_m2 <= 0:
        return Result(float("nan"), "K", "CAN-WALL-DEPTHSWEEP",
                      status=Status.PHYSICALLY_INVALID, inputs=inputs,
                      notes="need 0 < x1 < x2, k > 0 and q > 0")
    a = x1 / (x2 - x1)
    amp = math.sqrt((1.0 + a) ** 2 + a ** 2)
    g = heat_flux_W_m2 / conductivity_W_mK
    rows = []
    crossover = None
    for dx in depth_uncertainties_m:
        d_th = sensor_uncertainty_K * amp
        d_ge_ind = abs(g) * dx * math.hypot(1.0, a)
        d_ge_cor = abs(g) * dx
        tot_ind = math.hypot(d_th, d_ge_ind)
        rows.append({"depth_uncertainty_m": dx,
                     "thermometric_K": d_th,
                     "geometric_independent_K": d_ge_ind,
                     "geometric_correlated_K": d_ge_cor,
                     "total_independent_K": tot_ind,
                     "geometry_dominates": d_ge_ind > d_th})
        if crossover is None and d_ge_ind > d_th:
            crossover = dx
    exact_cross = (sensor_uncertainty_K * amp
                   / (abs(g) * math.hypot(1.0, a)))
    return Result(
        rows, "K", "CAN-WALL-DEPTHSWEEP", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="analytic sensitivity coefficients of the two-sensor "
               "reconstruction, differentiated exactly",
        source_locator="-",
        validity_domain="1-D planar reconstruction; says nothing about the "
                        "2-D/3-D bias, which is a separate and additive term",
        uncertainty="exact given the inputs",
        notes=f"gradient {g / 1e3:.1f} K/mm, amplification {amp:.3f}x; "
              f"geometry overtakes thermometry at a depth uncertainty of "
              f"{exact_cross * 1e6:.2f} um"
              + (f" (first swept point past it: {crossover * 1e6:.0f} um)"
                 if crossover else " (not reached in the swept range)"),
        inputs={**inputs, "a": a, "amplification": amp,
                "gradient_K_per_m": g,
                "crossover_depth_uncertainty_m": exact_cross,
                "correlated_assumed": depth_errors_correlated})


def optimal_sensor_depths(
    *, wall_thickness_m: float, heat_flux_W_m2: float,
    conductivity_W_mK: float, sensor_uncertainty_K: float,
    depth_uncertainty_m: float,
    min_depth_m: float, min_separation_m: float,
    keep_out_from_heated_face_m: float,
    n_grid: int = 200) -> Result:
    r"""Best depths for the two sensors, with each constraint named.

    Minimises the analytic total :math:`\sqrt{\delta T_{therm}^2 +
    \delta T_{geom}^2}` subject to manufacturing and instrumentation
    constraints. **The unconstrained optimum is degenerate and the function
    says so:** both terms fall monotonically as :math:`a\to 0`, i.e. as the
    first sensor approaches the wetted surface and the pair separates as far
    as possible. So the answer is set entirely by the constraints, and the
    honest output is *which constraint binds*, not a number that looks
    optimised.
    """
    inputs = {"t_wall_m": wall_thickness_m, "q_W_m2": heat_flux_W_m2,
              "k_W_mK": conductivity_W_mK, "u_T_K": sensor_uncertainty_K,
              "u_x_m": depth_uncertainty_m, "min_depth_m": min_depth_m,
              "min_sep_m": min_separation_m,
              "keep_out_m": keep_out_from_heated_face_m}
    hi = wall_thickness_m - keep_out_from_heated_face_m
    if hi <= min_depth_m + min_separation_m:
        return Result(
            float("nan"), "m", "CAN-WALL-SENSOROPT",
            status=Status.PHYSICALLY_INVALID, inputs=inputs,
            evidence_level=EvidenceLevel.E1,
            notes=f"the wall is too thin: usable depth band is "
                  f"{min_depth_m * 1e3:.3f}-{hi * 1e3:.3f} mm, which cannot "
                  f"hold two sensors {min_separation_m * 1e3:.3f} mm apart. "
                  f"THICKEN THE WALL -- this is not a placement problem.")
    g = heat_flux_W_m2 / conductivity_W_mK

    def total(x1: float, x2: float) -> float:
        a = x1 / (x2 - x1)
        return math.hypot(sensor_uncertainty_K * math.sqrt(
            (1.0 + a) ** 2 + a ** 2),
            abs(g) * depth_uncertainty_m * math.hypot(1.0, a))

    best = None
    for i in range(n_grid + 1):
        x1 = min_depth_m + (hi - min_separation_m - min_depth_m) * i / n_grid
        for j in range(n_grid + 1):
            lo2 = x1 + min_separation_m
            x2 = lo2 + (hi - lo2) * j / n_grid
            if x2 <= x1:
                continue
            v = total(x1, x2)
            if best is None or v < best[0]:
                best = (v, x1, x2)
    v, x1, x2 = best
    binds = []
    if abs(x1 - min_depth_m) < 1e-9:
        binds.append("MANUFACTURING: minimum depth below the wetted surface")
    if abs(x2 - hi) < 1e-9:
        binds.append("DESIGN: keep-out from the heated face")
    if abs((x2 - x1) - min_separation_m) < 1e-9:
        binds.append("INSTRUMENTATION: minimum sensor separation")
    return Result(
        {"x1_m": x1, "x2_m": x2, "total_uncertainty_K": v},
        "m", "CAN-WALL-SENSOROPT", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="analytic sensitivity of the two-sensor reconstruction, "
               "minimised on a grid over the feasible depth band",
        source_locator="-",
        validity_domain="the 1-D reconstruction only; the 2-D/3-D bias is "
                        "NOT in this objective and is reported separately",
        uncertainty="the objective is exact; the grid resolution is "
                    f"{(hi - min_depth_m) / n_grid * 1e6:.1f} um",
        notes=f"optimum x1 {x1 * 1e6:.0f} um, x2 {x2 * 1e6:.0f} um, total "
              f"{v:.2f} K. THE OPTIMUM IS CONSTRAINT-BOUND, NOT INTERIOR: "
              f"the objective falls monotonically as the first sensor "
              f"approaches the wetted surface and the pair separates, so "
              f"what fixes the answer is "
              + ("; ".join(binds) if binds else "nothing -- check the bounds"),
        inputs={**inputs, "binding_constraints": binds,
                "gradient_K_per_m": g})


# ===========================================================================
# 6. SYNTHETIC VERIFICATION -- E1 ONLY, NEVER VALIDATION
# ===========================================================================

@dataclass(frozen=True)
class AnalyticField:
    """A constructed T(x) used ONLY to verify the reconstruction algebra."""

    name: str
    T: Callable[[float], float]        # depth below the wetted surface -> K
    exact_wall_K: float
    expected: str


def _mk_fields(t_wall: float, g: float, T_w: float) -> dict[str,
                                                            AnalyticField]:
    return {
        "linear": AnalyticField(
            "linear", lambda x: T_w + g * x, T_w,
            "EXACT -- the reconstruction is built for this field"),
        "quadratic": AnalyticField(
            "quadratic", lambda x: T_w + g * x + 0.25 * g / t_wall * x * x,
            T_w, "BIASED -- curvature is what the linear fit cannot see"),
        "cubic": AnalyticField(
            "cubic",
            lambda x: T_w + g * x + 0.15 * g / (t_wall ** 2) * x ** 3, T_w,
            "BIASED, and the sign depends on where the sensors sit"),
        "exponential": AnalyticField(
            "exponential",
            lambda x: T_w + g * t_wall * (math.exp(x / t_wall) - 1.0), T_w,
            "BIASED -- the k(T) variation makes the real profile curve "
            "like this"),
    }


ANALYTIC_FIELDS = _mk_fields


def verify_reconstruction_on_field(
    *, fieldspec: AnalyticField, sensor_depths_m: tuple[float, float]
) -> Result:
    """Feed a known field through the 1-D reconstruction. **E1 only.**

    This is *verification* -- the algebra recovers a constructed answer. It
    is **not** validation and may never be described as such: no measurement
    is involved anywhere in it.
    """
    x1, x2 = sensor_depths_m
    if not (0.0 < x1 < x2):
        return Result(float("nan"), "K", "CAN-WALL-VERIFY",
                      status=Status.PHYSICALLY_INVALID,
                      notes="need 0 < x1 < x2")
    t1, t2 = fieldspec.T(x1), fieldspec.T(x2)
    g = (t2 - t1) / (x2 - x1)
    t_1d = t1 - g * x1
    err = t_1d - fieldspec.exact_wall_K
    exact = abs(err) < 1e-9
    return Result(
        err, "K", "CAN-WALL-VERIFY", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC,
        validation=Validation.NONE,
        source="constructed analytic temperature field; no measurement is "
               "involved",
        source_locator="-",
        validity_domain="algebraic verification of the reconstruction only",
        uncertainty="exact arithmetic",
        notes=f"{fieldspec.name}: T_1D {t_1d:.6f} K against an exact wall "
              f"temperature of {fieldspec.exact_wall_K:.6f} K, error "
              f"{err:+.6f} K -- {'EXACT' if exact else 'BIASED'}. "
              f"{fieldspec.expected}. THIS IS E1 VERIFICATION AND IS NOT "
              f"EXPERIMENTAL VALIDATION.",
        inputs={"field": fieldspec.name, "x1_m": x1, "x2_m": x2,
                "T1_K": t1, "T2_K": t2, "T_1D_K": t_1d,
                "error_K": err, "exact": exact})


def conductivity_sensitivity(
    *, solid: SolidProperties, section: ChannelSection,
    heat_flux_W_m2: float, film_coefficient_W_m2K: float,
    bulk_temperature_K: float, sensor_depths_m: tuple[float, float],
    relative_perturbations: Sequence[float] = (-0.04, 0.0, 0.04),
    **solve_kw) -> Result:
    r"""``dT_wall/dk``: how much does the witness-coupon measurement matter?

    The specimen's k(T) is now available (Chen et al. 2023) but carries a
    stated **+/-4 % sample-to-sample variation** and was measured on **HIP'd**
    material from other vendors. This function perturbs k by that band and
    reports the effect on the reconstructed wall temperature, which is the
    quantitative answer to *"is a laser-flash witness coupon worth its cost
    on this build?"*

    Note the structure of the answer, which is not obvious: the **two-sensor
    reconstruction does not divide by k at all** -- it measures the
    gradient. So k enters the reconstructed wall temperature only through
    the *field*, i.e. through how the real temperature distributes itself,
    not through the formula.
    """
    if not NUMERICS_AVAILABLE:                              # pragma: no cover
        return Result(float("nan"), "K", "CAN-WALL-KSENS",
                      status=Status.INSUFFICIENT_EVIDENCE,
                      notes="numerics unavailable")
    rows = []
    base = None
    for p in relative_perturbations:
        scaled = SolidProperties(
            name=f"{solid.name} x{1 + p:.3f}",
            k=lambda T, _p=p: solid.k(T) * (1.0 + _p),
            valid_range_K=solid.valid_range_K, source=solid.source,
            source_locator=solid.source_locator,
            source_access=solid.source_access, condition=solid.condition,
            uncertainty=solid.uncertainty)
        sol = solve_conduction(
            section=section, solid=scaled, heat_flux_W_m2=heat_flux_W_m2,
            film_coefficient_W_m2K=film_coefficient_W_m2K,
            bulk_temperature_K=bulk_temperature_K, **solve_kw)
        if sol.T is None:
            continue
        b = reconstruction_bias(sol, sensor_depths_m=sensor_depths_m)
        t1d = b.inputs.get("T_1D_K", float("nan"))
        rows.append({"dk_over_k": p, "T_1D_K": t1d,
                     "T_true_local_K": b.inputs.get("T_true_local_K"),
                     "bias_reconstruction_K": b.inputs.get(
                         "bias_reconstruction_K")})
        if p == 0.0:
            base = t1d
    slope = float("nan")
    if len(rows) >= 2 and base is not None and base == base:
        lo, hi = rows[0], rows[-1]
        if hi["dk_over_k"] != lo["dk_over_k"]:
            slope = ((hi["T_1D_K"] - lo["T_1D_K"])
                     / ((hi["dk_over_k"] - lo["dk_over_k"]) * 100.0))
    return Result(
        slope, "K per % of k", "CAN-WALL-KSENS", status=Status.PASS,
        evidence_level=EvidenceLevel.E1,
        verification=Verification.ALGEBRAIC, validation=Validation.NONE,
        source="finite difference of the solved field with respect to a "
               "uniform scaling of k(T)",
        source_locator="perturbation band from the source's own stated "
                       "sample-to-sample variation: " + solid.uncertainty,
        validity_domain="the solved section at the stated condition",
        uncertainty="a numerical derivative; inherits the field residual",
        notes=f"dT_wall/dk = {slope:+.4f} K per 1% of k. The two-sensor "
              f"reconstruction does not divide by k -- it measures the "
              f"gradient -- so k enters only through the field.",
        inputs={"rows": rows, "slope_K_per_percent": slope,
                "material": solid.name})
