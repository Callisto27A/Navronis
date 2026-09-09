"""
Engine inner-wall profile generator.
=====================================
Generates the complete (x, r) coordinate arrays for the axisymmetric
inner wall of a liquid rocket thrust chamber, from injector face to
nozzle exit lip.

Sections (left to right):
  1. Injector face (flat wall at x=0)
  2. Cylindrical barrel (constant radius = R_c)
  3. Convergent cone (linear taper from R_c to R_t)
  4. Throat (minimum radius R_t, with upstream/downstream circular arcs)
  5. Divergent nozzle (conical or Rao bell contour, from R_t to R_e)

All dimensions are in SI (metres). The origin x=0 is at the injector face.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Tuple

from kryptonis.propulsion_equations.contour import (
    generate_rao_bell_coordinates,
    calculate_conical_nozzle_length,
)


@dataclass
class ChamberProfile:
    """Complete axisymmetric inner-wall profile of a thrust chamber + nozzle."""

    # Input parameters (all SI)
    throat_diameter_m: float
    chamber_diameter_m: float
    convergent_half_angle_deg: float
    cylindrical_length_m: float
    expansion_ratio: float
    nozzle_type: str  # "conical" or "bell"
    divergent_half_angle_deg: float  # conical: exit half-angle; bell: initial wall angle
    bell_exit_angle_deg: float  # bell nozzles only
    bell_fractional_length: float  # fraction of equivalent 15-deg cone

    # Derived
    throat_radius_m: float = field(init=False)
    chamber_radius_m: float = field(init=False)
    exit_radius_m: float = field(init=False)
    convergent_length_m: float = field(init=False)
    divergent_length_m: float = field(init=False)
    total_length_m: float = field(init=False)

    # Profile points: list of (x, r) in metres
    x: List[float] = field(init=False, default_factory=list)
    r: List[float] = field(init=False, default_factory=list)

    # Section boundaries (x-coordinates)
    x_injector: float = field(init=False, default=0.0)
    x_cyl_end: float = field(init=False)
    x_throat: float = field(init=False)
    x_exit: float = field(init=False)

    def __post_init__(self):
        self.throat_radius_m = self.throat_diameter_m / 2.0
        self.chamber_radius_m = self.chamber_diameter_m / 2.0
        self.exit_radius_m = self.throat_radius_m * math.sqrt(self.expansion_ratio)

        # Convergent cone length
        self.convergent_length_m = (self.chamber_radius_m - self.throat_radius_m) / \
            math.tan(math.radians(self.convergent_half_angle_deg))

        # Section x-boundaries
        self.x_cyl_end = self.cylindrical_length_m
        self.x_throat = self.x_cyl_end + self.convergent_length_m

        # Divergent section
        if self.nozzle_type == "bell":
            cone_len = calculate_conical_nozzle_length(
                self.throat_radius_m, self.exit_radius_m, 15.0)
            self.divergent_length_m = self.bell_fractional_length * cone_len
        else:
            self.divergent_length_m = calculate_conical_nozzle_length(
                self.throat_radius_m, self.exit_radius_m,
                self.divergent_half_angle_deg)

        self.x_exit = self.x_throat + self.divergent_length_m
        self.total_length_m = self.x_exit

        self._generate_points()

    def _generate_points(self, n_cyl: int = 20, n_conv: int = 30,
                          n_div: int = 60):
        """Build the (x, r) arrays for each section."""
        xs: List[float] = []
        rs: List[float] = []

        R_c = self.chamber_radius_m
        R_t = self.throat_radius_m
        R_e = self.exit_radius_m

        # 1. Cylindrical barrel: x = 0 .. x_cyl_end, r = R_c
        for i in range(n_cyl + 1):
            x = self.x_injector + i * self.cylindrical_length_m / n_cyl
            xs.append(x)
            rs.append(R_c)

        # 2. Convergent cone: x = x_cyl_end .. x_throat, r = R_c .. R_t
        for i in range(1, n_conv + 1):
            frac = i / n_conv
            x = self.x_cyl_end + frac * self.convergent_length_m
            r = R_c - frac * (R_c - R_t)
            xs.append(x)
            rs.append(r)

        # 3. Divergent nozzle
        if self.nozzle_type == "bell":
            bell_pts = generate_rao_bell_coordinates(
                R_t=R_t, R_e=R_e,
                fractional_length=self.bell_fractional_length,
                theta_n_deg=self.divergent_half_angle_deg,
                theta_e_deg=self.bell_exit_angle_deg,
                num_points=n_div,
            )
            for bx, by in bell_pts[1:]:  # skip first point (duplicate throat)
                xs.append(self.x_throat + bx)
                rs.append(by)
        else:
            # Conical divergent
            for i in range(1, n_div + 1):
                frac = i / n_div
                x = self.x_throat + frac * self.divergent_length_m
                r = R_t + frac * (R_e - R_t)
                xs.append(x)
                rs.append(r)

        self.x = xs
        self.r = rs

    def points_mm(self) -> List[Tuple[float, float]]:
        """Return profile as (x_mm, r_mm) tuples."""
        return [(xi * 1e3, ri * 1e3) for xi, ri in zip(self.x, self.r)]

    def as_dict(self) -> dict:
        """Serialisable dictionary of all parameters and profile."""
        return {
            "inputs": {
                "throat_diameter_mm": self.throat_diameter_m * 1e3,
                "chamber_diameter_mm": self.chamber_diameter_m * 1e3,
                "convergent_half_angle_deg": self.convergent_half_angle_deg,
                "cylindrical_length_mm": self.cylindrical_length_m * 1e3,
                "expansion_ratio": self.expansion_ratio,
                "nozzle_type": self.nozzle_type,
                "divergent_half_angle_deg": self.divergent_half_angle_deg,
                "bell_exit_angle_deg": self.bell_exit_angle_deg,
                "bell_fractional_length": self.bell_fractional_length,
            },
            "derived": {
                "throat_radius_mm": self.throat_radius_m * 1e3,
                "chamber_radius_mm": self.chamber_radius_m * 1e3,
                "exit_radius_mm": self.exit_radius_m * 1e3,
                "exit_diameter_mm": self.exit_radius_m * 2e3,
                "convergent_length_mm": self.convergent_length_m * 1e3,
                "divergent_length_mm": self.divergent_length_m * 1e3,
                "total_length_mm": self.total_length_m * 1e3,
            },
            "section_boundaries_mm": {
                "injector_face": 0.0,
                "cylinder_end": self.x_cyl_end * 1e3,
                "throat": self.x_throat * 1e3,
                "nozzle_exit": self.x_exit * 1e3,
            },
            "profile_points_mm": self.points_mm(),
            "num_points": len(self.x),
        }


def generate_chamber_profile(
    *,
    throat_diameter_m: float,
    chamber_diameter_m: float,
    cylindrical_length_m: float,
    convergent_half_angle_deg: float = 25.0,
    expansion_ratio: float = 20.0,
    nozzle_type: str = "bell",
    divergent_half_angle_deg: float = 30.0,
    bell_exit_angle_deg: float = 8.0,
    bell_fractional_length: float = 0.80,
) -> ChamberProfile:
    """High-level API to generate a complete chamber + nozzle profile."""
    return ChamberProfile(
        throat_diameter_m=throat_diameter_m,
        chamber_diameter_m=chamber_diameter_m,
        convergent_half_angle_deg=convergent_half_angle_deg,
        cylindrical_length_m=cylindrical_length_m,
        expansion_ratio=expansion_ratio,
        nozzle_type=nozzle_type,
        divergent_half_angle_deg=divergent_half_angle_deg,
        bell_exit_angle_deg=bell_exit_angle_deg,
        bell_fractional_length=bell_fractional_length,
    )
