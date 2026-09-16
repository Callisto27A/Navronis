"""
Kryptonis Propulsion Equations: Storable Propellants Benchmark & Regression Tests
================================================================================
Integrates and validates the CEA deck benchmarks, SP-125 Table 4-1 mappings,
and multi-regime thermochemical presets audited across PR #6 and PR #8.
"""

import math
import pytest

from kryptonis.propulsion_equations.units import Status
from kryptonis.propulsion_equations.chamber import (
    characteristic_length,
    get_thermochemical_preset,
    THERMOCHEMICAL_PRESETS,
    CEA_BENCHMARK_PRESETS_ALT,
)
from kryptonis.propulsion_equations.combustor import CombustorDesign


def test_storable_propellant_sp125_mappings():
    # MMH and N2H4 map onto N2O4/hydrazine-base (30-35 inches -> midpoint 32.5 in = 0.8255 m)
    res_mmh = characteristic_length(fuel="MMH")
    assert res_mmh.status == Status.UNVALIDATED_ASSUMPTION
    assert pytest.approx(res_mmh.value, rel=1e-4) == 32.5 * 0.0254

    res_n2h4 = characteristic_length(fuel="N2H4")
    assert res_n2h4.status == Status.UNVALIDATED_ASSUMPTION
    assert pytest.approx(res_n2h4.value, rel=1e-4) == 32.5 * 0.0254

    # Alias check
    res_nto_mmh = characteristic_length(fuel="N2O4/MMH")
    assert res_nto_mmh.status == Status.UNVALIDATED_ASSUMPTION
    assert pytest.approx(res_nto_mmh.value, rel=1e-4) == 32.5 * 0.0254

    # Ethanol has no SP-125 Table 4-1 row, must return INSUFFICIENT_EVIDENCE
    # and dynamic notes must reference 'ETHANOL' directly
    res_eth = characteristic_length(fuel="ETHANOL")
    assert res_eth.status == Status.INSUFFICIENT_EVIDENCE
    assert math.isnan(res_eth.value)
    assert "ETHANOL" in res_eth.notes


def test_storable_propellants_combustor_design_matrix():
    # 1. N2O4 / MMH Hypergolic in-space / RCS engine
    des_mmh = CombustorDesign(
        thrust=10000.0,
        chamber_pressure=15.0e5,
        mixture_ratio=1.65,
        propellant="N2O4/MMH",
    )
    res_mmh = des_mmh.solve()
    assert res_mmh.throat_diameter > 0.0
    assert res_mmh.chamber_diameter > res_mmh.throat_diameter
    assert 1600.0 <= res_mmh.c_star_ideal <= 1850.0
    assert pytest.approx(res_mmh.design.mixture_ratio, rel=1e-4) == 1.65

    # 2. Monopropellant Hydrazine catalytic thruster (Tc ~ 1200 K, O/F = 0)
    des_n2h4 = CombustorDesign(
        thrust=500.0,
        chamber_pressure=10.0e5,
        mixture_ratio=0.0,
        propellant="HYDRAZINE",
    )
    res_n2h4 = des_n2h4.solve()
    assert res_n2h4.throat_diameter > 0.0
    assert 1200.0 <= res_n2h4.c_star_ideal <= 1450.0
    assert res_n2h4.oxidizer_mass_flow == 0.0
    assert res_n2h4.fuel_mass_flow > 0.0

    # 3. Green storable N2O / Ethanol thruster
    des_green = CombustorDesign(
        thrust=5000.0,
        chamber_pressure=25.0e5,
        mixture_ratio=4.0,
        propellant="N2O/ETHANOL",
    )
    res_green = des_green.solve()
    assert res_green.throat_diameter > 0.0
    assert 1500.0 <= res_green.c_star_ideal <= 1750.0


def test_pr8_alternative_deck_b_benchmarks():
    assert "N2O4/MMH" in CEA_BENCHMARK_PRESETS_ALT
    assert "Hydrazine" in CEA_BENCHMARK_PRESETS_ALT
    assert "N2O/Ethanol" in CEA_BENCHMARK_PRESETS_ALT

    # Verify alternative preset retrieval via use_alt_deck flag
    mmh_alt = get_thermochemical_preset("N2O4/MMH", use_alt_deck=True)
    assert mmh_alt["tc"] == 3150.0
    assert mmh_alt["isp"] == 285.0

    n2o_alt = get_thermochemical_preset("N2O/Ethanol", use_alt_deck=True)
    assert n2o_alt["of"] == 4.5
    assert n2o_alt["tc"] == 2900.0


def test_cli_and_alias_coverage():
    from kryptonis.propulsion_equations.cli import run_injector_sizing, run_chamber_sizing

    assert get_thermochemical_preset("NTO/MMH")["tc"] == 3120.0
    assert get_thermochemical_preset("AEROZINE50")["gamma"] == 1.24
    assert get_thermochemical_preset("N2O/ETOH")["of"] == 4.0
    assert get_thermochemical_preset("MONOPROPELLANT_HYDRAZINE")["tc"] == 1200.0

    # CLI execution checks
    assert run_chamber_sizing(thrust_n=4000.0, pc_bar=15.0, propellants="NTO/MMH") == 0
    assert run_injector_sizing(thrust_n=4000.0, pc_bar=15.0, propellants="N2O4/MMH", injector_type="swirl") == 0
    assert run_injector_sizing(thrust_n=3000.0, pc_bar=20.0, propellants="N2O/ETHANOL", injector_type="pintle") == 0

    # Monopropellant guard test: Hydrazine with O/F = 0 returns 0 cleanly instead of throwing ValueError
    assert run_injector_sizing(thrust_n=1000.0, pc_bar=10.0, propellants="HYDRAZINE", injector_type="pintle") == 0