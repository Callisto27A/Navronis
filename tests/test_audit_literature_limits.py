"""
Audit & Regression Tests: Primary Literature Review of Combustor & Bartz Limits
==============================================================================
Addresses Navronis Issue #1:
1. Contraction Ratio (eps_c): Humble (1995) correlation across small-scale collegiate
   engines (Dt < 25 mm) vs heavy orbital boosters (Dt > 0.5 m).
2. Bartz Convective Heat Transfer: sigma property variation factor behavior, bounds,
   and supercritical pressure regime (Pc > 100 bar).
3. Cylindrical Barrel Volume Formulation: Algebraic and geometric equivalence
   between frustum integration and Huzel & Huang (NASA SP-125 Eq. 4-5) closed-form identity.

References:
- Humble, R. W., et al., Space Propulsion Analysis and Design, McGraw-Hill, 1995.
- Huzel, D. K., & Huang, D. H., Modern Engineering for Design of Liquid-Propellant
  Rocket Engines, NASA SP-125, 1967.
- Bartz, D. R., 'A Simple Equation for Rapid Estimation of Rocket Nozzle Convective
  Heat Transfer Coefficients', Jet Propulsion 27(1), 1957.
- Sutton, G. P., & Biblarz, O., Rocket Propulsion Elements, 9th ed., Wiley, 2016.
"""

import math
import pytest

from kryptonis.propulsion_equations.chamber import (
    contraction_ratio,
    chamber_diameter,
    chamber_volume,
    convergent_length,
    convergent_volume,
    cylinder_length,
    throat_area,
    throat_diameter,
    EPS_C_TURBOPUMP,
    EPS_C_PRESSURE_FED,
    EPS_C_ABSOLUTE_MIN,
    EPS_C_ABSOLUTE_MAX,
)
from kryptonis.propulsion_equations.bartz import (
    bartz_sigma,
    recovery_temperature,
    bartz_film_coefficient,
    bartz_prandtl,
    bartz_viscosity,
)
from kryptonis.propulsion_equations.units import Status


# ===========================================================================
# 1. CONTRACTION RATIO (Humble 1995) LITERATURE AUDIT ACROSS ENGINE SCALES
# ===========================================================================

class TestContractionRatioLiteratureAudit:
    """Verifies Humble (1995) correlation: eps_c = 1.25 + 8.0 * (Dt[cm])^(-0.6)
    across small-scale collegiate thrusters vs massive orbital boosters."""

    def test_small_scale_collegiate_engine_boundary(self):
        """Small collegiate engines (Dt < 25 mm) yield high contraction ratios.
        At Dt = 10 mm (1.0 cm), eps_c = 1.25 + 8.0*(1.0)^-0.6 = 9.25.
        This exceeds SP-125's empirical envelope (max 6.0), triggering
        OUT_OF_CORRELATION_RANGE as expected for collegiate subscale designs.
        """
        d_throat_10mm = 0.010  # 10 mm
        res_10mm = contraction_ratio(throat_diameter_m=d_throat_10mm)
        assert pytest.approx(res_10mm.value, rel=1e-5) == 9.25
        assert Status.OUT_OF_CORRELATION_RANGE in res_10mm.all_statuses

        d_throat_20mm = 0.020  # 20 mm = 2.0 cm
        res_20mm = contraction_ratio(throat_diameter_m=d_throat_20mm)
        expected_20mm = 1.25 + 8.0 * (2.0 ** -0.6)  # ~6.528
        assert pytest.approx(res_20mm.value, rel=1e-5) == expected_20mm
        assert Status.OUT_OF_CORRELATION_RANGE in res_20mm.all_statuses

    def test_mid_scale_sounding_rocket_regime(self):
        """Mid-scale sounding rockets (e.g. Dt = 50 to 100 mm):
        At Dt = 80 mm (8.0 cm), eps_c = 1.25 + 8.0*(8.0)^-0.6 ~ 3.54.
        This falls squarely inside the general SP-125 envelope (1.3 to 6.0)
        and pressure-fed bounds (2.0 to 5.0).
        """
        d_throat_80mm = 0.080  # 80 mm = 8.0 cm
        res_80mm = contraction_ratio(throat_diameter_m=d_throat_80mm, feed_system="pressure_fed")
        expected_80mm = 1.25 + 8.0 * (8.0 ** -0.6)
        assert pytest.approx(res_80mm.value, rel=1e-5) == expected_80mm
        assert Status.OUT_OF_CORRELATION_RANGE not in res_80mm.all_statuses
        assert EPS_C_PRESSURE_FED[0] <= res_80mm.value <= EPS_C_PRESSURE_FED[1]

    def test_heavy_orbital_booster_asymptotic_limit(self):
        """Heavy orbital booster engines (Dt >= 0.5 m to 1.0 m, e.g. F-1, RD-170, RS-25):
        As Dt grows large, the diameter term 8.0 * Dt^-0.6 asymptotically approaches 0,
        bringing eps_c down to ~ 1.7 to 2.0.
        This matches the SP-125 turbopump range (1.3 to 2.5) where minimizing chamber
        diameter reduces chamber hoop stress and dry engine mass.
        """
        # 500 mm throat (0.5 m = 50 cm)
        d_throat_500mm = 0.500
        res_500mm = contraction_ratio(throat_diameter_m=d_throat_500mm, feed_system="turbopump")
        expected_500mm = 1.25 + 8.0 * (50.0 ** -0.6)  # ~2.017
        assert pytest.approx(res_500mm.value, rel=1e-5) == expected_500mm
        assert Status.OUT_OF_CORRELATION_RANGE not in res_500mm.all_statuses
        assert EPS_C_TURBOPUMP[0] <= res_500mm.value <= EPS_C_TURBOPUMP[1]

        # 1000 mm throat (1.0 m = 100 cm, F-1 scale)
        d_throat_1000mm = 1.000
        res_1000mm = contraction_ratio(throat_diameter_m=d_throat_1000mm, feed_system="turbopump")
        expected_1000mm = 1.25 + 8.0 * (100.0 ** -0.6)  # ~1.755
        assert pytest.approx(res_1000mm.value, rel=1e-5) == expected_1000mm
        assert EPS_C_TURBOPUMP[0] <= res_1000mm.value <= EPS_C_TURBOPUMP[1]

    def test_supplied_designer_value_takes_precedence(self):
        """A human or optimization loop supplying an explicit contraction ratio
        overrides the Humble empirical fit while preserving bounds checks."""
        res_custom = contraction_ratio(throat_diameter_m=0.015, supplied_value=3.5, feed_system="pressure_fed")
        assert res_custom.value == 3.5
        assert Status.UNVALIDATED_ASSUMPTION not in res_custom.all_statuses


# ===========================================================================
# 2. CONVERGENT VOLUME & CYLINDRICAL BARREL ALGEBRAIC IDENTITIES
# ===========================================================================

class TestConvergentVolumeAnalyticalEquivalence:
    """Verifies that the conical frustum volume integration is algebraically
    and numerically identical to Huzel & Huang (NASA SP-125 Eq. 4-5) closed-form:
    V_conv = (pi * Dt^3 / (24 * tan(theta))) * (eps_c^(3/2) - 1)
           = (At * Dt / (6 * tan(theta))) * (eps_c^(3/2) - 1)
    """

    @pytest.mark.parametrize("Dt, eps_c, half_angle_deg", [
        (0.020, 4.0, 20.0),
        (0.080, 2.8, 15.0),
        (0.250, 1.8, 25.0),
        (0.750, 1.5, 30.0),
    ])
    def test_analytical_huzel_huang_equivalence(self, Dt, eps_c, half_angle_deg):
        """Compare standard conical frustum function with Huzel & Huang Eq. 4-5 closed form."""
        Dc = Dt * math.sqrt(eps_c)
        At = 0.25 * math.pi * Dt ** 2

        # Direct Navronis implementation
        res_vconv = convergent_volume(
            throat_diameter_m=Dt,
            chamber_diameter_m=Dc,
            half_angle_deg=half_angle_deg
        )
        assert res_vconv.status == Status.PASS

        # NASA SP-125 Eq. (4-5) closed-form analytical identity:
        # V_conv = (At * Dt / (6 * tan(theta))) * (eps_c^(3/2) - 1)
        theta_rad = math.radians(half_angle_deg)
        expected_vconv = (At * Dt / (6.0 * math.tan(theta_rad))) * (eps_c ** 1.5 - 1.0)

        assert pytest.approx(res_vconv.value, rel=1e-12) == expected_vconv

    def test_chamber_mass_conservation_closure(self):
        """Verify that Vc = V_barrel + V_conv identically equals L* * At."""
        At = 0.010  # m^2
        Dt = math.sqrt(4.0 * At / math.pi)
        eps_c = 2.5
        Ac = At * eps_c
        Dc = math.sqrt(4.0 * Ac / math.pi)
        L_star = 1.10  # 1.10 m

        res_Vc = chamber_volume(l_star_m=L_star, throat_area_m2=At)
        res_Vconv = convergent_volume(throat_diameter_m=Dt, chamber_diameter_m=Dc, half_angle_deg=22.0)
        res_Lcyl = cylinder_length(
            chamber_volume_m3=res_Vc.value,
            convergent_volume_m3=res_Vconv.value,
            chamber_area_m2=Ac
        )

        assert res_Lcyl.status == Status.PASS
        # Reconstruct total chamber volume from cylinder and cone
        reconstructed_Vc = Ac * res_Lcyl.value + res_Vconv.value
        assert pytest.approx(reconstructed_Vc, rel=1e-12) == res_Vc.value


# ===========================================================================
# 3. BARTZ HEAT FLUX & PROPERTY VARIATION FACTOR SIGMA BOUNDARY AUDIT
# ===========================================================================

class TestBartzHeatFluxBoundaryAudit:
    """Audits Bartz (1957) property-variation factor sigma across temperature
    ratios and Mach numbers, documenting supercritical chamber pressure regimes."""

    def test_sigma_physical_monotonicity(self):
        """As wall temperature Tw approaches recovery/stagnation temperature T0,
        the cold-wall density enhancement diminishes, and sigma monotonically decreases towards ~0.8-1.0.
        At cold wall (Tw/T0 = 0.2), near-wall density is high, raising sigma > 1.
        """
        T0 = 3400.0  # K
        gamma = 1.20
        mach = 1.0   # Throat condition

        sigmas = []
        for ratio in [0.2, 0.4, 0.6, 0.8, 1.0]:
            Tw = ratio * T0
            res = bartz_sigma(wall_temperature_K=Tw, stagnation_temperature_K=T0, gamma=gamma, mach=mach)
            assert res.status == Status.PASS
            sigmas.append(res.value)

        # Strictly decreasing as wall temperature increases (cold-wall density effect)
        for i in range(len(sigmas) - 1):
            assert sigmas[i] > sigmas[i + 1]

        # Cold wall raises sigma above 1.0 (physical behavior confirmed by SP-125 Eq. 4-14)
        assert sigmas[0] > 1.0

    def test_bartz_film_coefficient_supercritical_boundary(self):
        """Bartz film coefficient calculation across high chamber pressure (Pc >= 100 bar).
        Ensures stability and positive finite values across supercritical pressure levels.
        """
        Dt = 0.100  # 100 mm throat
        At = 0.25 * math.pi * Dt ** 2
        R_curv = 0.150  # 150 mm curvature radius
        c_star = 1780.0
        gamma = 1.20
        Tc = 3500.0
        M_mol = 0.022  # 22 g/mol
        Tw = 800.0     # 800 K regenerative channel wall

        # Test across subcritical (20 bar), supercritical (100 bar, 200 bar, 300 bar)
        for Pc_bar in [20.0, 100.0, 200.0, 300.0]:
            Pc_Pa = Pc_bar * 1e5
            res_hg = bartz_film_coefficient(
                chamber_pressure_Pa=Pc_Pa,
                c_star_m_s=c_star,
                throat_diameter_m=Dt,
                throat_radius_curvature_m=R_curv,
                specific_heat_J_kgK=1800.0,
                stagnation_temperature_K=Tc,
                wall_temperature_K=Tw,
                mach=1.0,
                area_ratio_local_over_throat=1.0,
                gamma=gamma,
                molar_mass_kg_per_mol=M_mol,
            )
            assert Status.PHYSICALLY_INVALID not in res_hg.all_statuses
            assert res_hg.value > 0.0
            # Heat transfer coefficient scales as ~ Pc^0.8
            assert math.isfinite(res_hg.value)
