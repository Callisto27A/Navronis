"""
Unit tests for the Day 1 CombustorDesign public API.
"""
import os
import pytest
from kryptonis.propulsion_equations import CombustorDesign, CombustorResult

def test_combustor_design_nominal():
    design = CombustorDesign(
        thrust=30000.0,
        chamber_pressure=12.0e6,
        mixture_ratio=2.6,
        propellant="LOX/CH4",
    )
    res = design.solve()
    assert isinstance(res, CombustorResult)
    assert res.total_mass_flow > 0
    assert res.fuel_mass_flow > 0
    assert res.oxidizer_mass_flow > 0
    assert abs(res.fuel_mass_flow + res.oxidizer_mass_flow - res.total_mass_flow) < 1e-6
    assert 0.03 < res.throat_diameter < 0.06
    assert res.chamber_diameter > res.throat_diameter
    assert res.chamber_length > 0.15
    assert res.heat_flux_screen > 1e6
    assert res.acoustic_modes["1T_Hz"] > 1000.0
    assert res.wall_thickness_screen > 0.001

def test_combustor_provenance_explain():
    design = CombustorDesign(
        thrust=35000.0,
        chamber_pressure=7.0e6,
        mixture_ratio=2.5,
        propellant="LOX/RP-1",
    )
    res = design.solve()
    expl = res.explain("throat_diameter")
    assert "NASA SP-125" in expl
    assert "PRELIMINARY_DESIGN" in expl

def test_combustor_exports(tmp_path):
    design = CombustorDesign(
        thrust=25000.0,
        chamber_pressure=5.0e6,
        mixture_ratio=6.0,
        propellant="LOX/LH2",
    )
    res = design.solve()
    json_file = str(tmp_path / "combustor.json")
    csv_file = str(tmp_path / "combustor.csv")
    res.to_json(json_file)
    res.to_csv(csv_file)
    assert os.path.exists(json_file)
    assert os.path.exists(csv_file)
