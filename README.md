# Kryptonis Propulsion — Rocket Propulsion Analytical Toolkit

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-brightgreen.svg)](https://python.org)
[![Release](https://img.shields.io/badge/Release-v0.1.0%20(Day%201:%20Combustor)-orange.svg)](#current-release-v010--combustor-preliminary-design)
[![Physics](https://img.shields.io/badge/Physics-Zero%20Silent%20Defaults-success.svg)](#engineering-sources--pedigree)

**Design a preliminary liquid-rocket combustion chamber from thrust, chamber pressure, mixture ratio and propellants — with every major result traceable to its engineering source.**

---

## Current Release: `v0.1.0` — Combustor Preliminary Design

This is **Day 1** of the Kryptonis open-source propulsion release program. We are releasing one standalone, source-traceable engineering tool at a time:

```text
DAY 1: Combustor Preliminary Design  [RELEASED: v0.1.0]
DAY 2: Nozzle Gas Dynamics & Rao Bell Contour (Next)
DAY 3: Bartz Convective Heat Transfer & Thermal Screening
DAY 4: Regenerative Cooling Channel Fluid Dynamics
DAY 5: Acoustic Stability & Combustion Chamber Modes
DAY 6: Injector Orifice Sizing & Momentum Ratio
DAY 7: Pressure Vessel Structural Hoop Stress
DAY 8: Turbopump Preliminary Sizing
```

---

## What It Calculates

Given four preliminary mission inputs (`thrust`, `chamber_pressure`, `mixture_ratio`, `propellant`):

1. **Propellant Mass Flows:** Total mass flow $\dot{m}$, fuel flow $\dot{m}_f$, oxidizer flow $\dot{m}_{ox}$.
2. **Throat Sizing:** Choked throat area $A_t$ and diameter $D_t$ (NASA SP-125).
3. **Chamber Sizing:** Contraction ratio $\epsilon_c$, chamber area $A_c$, and diameter $D_c$.
4. **Axial Geometry:** Truncated convergent length $L_{conv}$, cylindrical barrel length $L_{cyl}$, total length $L_c$.
5. **Chamber Volume & Stay Time:** Combustion volume $V_c$ from characteristic stay length $L^*$ and stay time $\tau_s$.
6. **Thermal Screening:** Peak throat convective heat flux $q_{throat}$ via Bartz (1957).
7. **Acoustic Mode Screening:** 1st Tangential (1T), 1st Radial (1R), and 1st Longitudinal (1L) acoustic modes (NASA SP-194).
8. **Structural Screening:** Minimum thin-shell wall thickness $t_w$ from material yield strength and safety factor (ASME Sec VIII).

---

## Quick Start (Python API)

```python
from kryptonis.propulsion_equations import CombustorDesign

# 1. Define preliminary engine requirements
engine = CombustorDesign(
    thrust=30000.0,            # 30 kN
    chamber_pressure=12.0e6,    # 120 bar (12.0 MPa)
    mixture_ratio=2.6,         # LOX / LCH4
    propellant="LOX/CH4",
)

# 2. Solve closed-form analytical equations
res = engine.solve()

# 3. Access results
print(f"Total Flow:       {res.total_mass_flow:.2f} kg/s")
print(f"Throat Diameter:  {res.throat_diameter * 1000:.2f} mm")
print(f"Chamber Diameter: {res.chamber_diameter * 1000:.2f} mm")
print(f"Chamber Length:   {res.chamber_length * 1000:.2f} mm")
print(f"1T Acoustic Mode: {res.acoustic_modes['1T_Hz']:.1f} Hz")

# 4. Inspect literature provenance of any quantity
print(res.explain("throat_diameter"))

# 5. Export for CAD / Spreadsheet
res.to_json("combustor.json")
res.to_csv("combustor.csv")
res.plot(save_path="combustor.png")
```

---

## Terminal Command-Line Interface (Instant CLI)

```bash
# Size a 30 kN LOX/CH4 chamber at 120 bar with exports
kryptonis-chamber --thrust 30000 --pc 120 --propellants LOX/CH4 --plot-save chamber.png --export-json chamber.json
```

---

## Example Output (30 kN LOX/CH4 @ 120 bar)

```text
==============================================================================
           KRYPTONIS PROPULSION ENGINE SIZER -- COMPONENT 1: COMBUSTOR
==============================================================================
Inputs:
  Thrust:                30.00 kN (30000 N)
  Chamber Pressure (Pc): 120.00 bar (12.00 MPa)
  Propellants:           LOX/CH4 (gamma=1.20, Tc=3450 K, Mw=22.0 g/mol)
  Characteristic L*:     1.00 m
------------------------------------------------------------------------------
1. THROAT & CHAMBER SIZING (NASA SP-125 / Huzel & Huang)
------------------------------------------------------------------------------
  Throat Area (At):            14.29 cm^2
  Throat Diameter (Dt):        42.65 mm
  Contraction Ratio (eps_c):    4.60
  Chamber Area (Ac):           65.73 cm^2
  Chamber Diameter (Dc):       91.48 mm
  Chamber Volume (Vc):         1.429 liters
------------------------------------------------------------------------------
2. AXIAL PROFILE & STAY TIME
------------------------------------------------------------------------------
  Convergent Half-Angle:        30.0 deg
  Convergent Length (Lconv):   42.29 mm
  Cylindrical Barrel (Lcyl):  193.62 mm
  Total Chamber Length (Lc):  235.91 mm
  Mean Stay / Res. Time:        1.35 ms
------------------------------------------------------------------------------
3. ACOUSTIC STABILITY FREQUENCIES (NASA SP-194)
------------------------------------------------------------------------------
  Chamber Speed of Sound:     1250.9 m/s
  1st Tangential Mode (1T):   8013.7 Hz  [Primary transverse buzz]
  1st Radial Mode (1R):      16677.3 Hz
  1st Longitudinal (1L):      2651.1 Hz  [Chugging / organ-pipe coupling]
------------------------------------------------------------------------------
4. MECHANICAL INTEGRITY (ASME Sec VIII / Thin-Shell Hoop)
------------------------------------------------------------------------------
  Material Yield Strength:     280.0 MPa (SF: 1.5)
  Min Wall Thickness:           2.94 mm
  Recommended Thickness:        3.68 mm (+25% machining margin)
==============================================================================
```

---

## Engineering Sources & Pedigree

Every formula implemented cites its original peer-reviewed or technical monograph source:
- **Chamber Sizing & Geometry:** Huzel & Huang, *Design of Liquid Propellant Rocket Engines* (NASA SP-125).
- **Contraction Ratio Policy:** Traced to Humble, *Space Propulsion Analysis and Design* (1995).
- **Acoustic Mode Eigenvalues:** NASA SP-194, *Liquid Propellant Rocket Combustion Instability*. Exact roots of $J'_m(\alpha) = 0$ via `scipy.special.jnp_zeros`.
- **Hot-Gas Convection:** Bartz, D. R. (1957), *A Simple Equation for Rapid Estimation of Rocket Nozzle Convective Heat Transfer Coefficients*, Jet Propulsion 27(1).

---

## Important Scope & Limitations

> [!IMPORTANT]
> **This library is a preliminary engineering calculation toolkit.**
> - It does **NOT** design manufacturing-ready rocket engines.
> - It does **NOT** replace 3D reacting CFD, conjugate heat transfer, or finite element analysis (FEA).
> - Acoustic modes are **frequencies only**, not a stability assessment (which requires hardware-measured injector response functions $n$ and $\tau$).
> - Advanced capabilities — automatic physics selection, multidisciplinary fixed-point cycle closure, 679-relation turbopump compilers, and manufacturing DFM — belong to the private Kryptonis engineering platform.

---

## Installation

```bash
git clone https://github.com/Ritik27Payak/Navronis.git
cd Navronis
pip install -e .
```

To include plotting dependencies:
```bash
pip install -e ".[dev]"
```

Run test suite:
```bash
pytest tests/ -v
```

---

## Next Planned Release

**Day 2: Nozzle Gas Dynamics & Rao Bell Contour (`v0.2.0`)**  
- Invertible Area-Mach relations
- Rao parabolic bell contour coordinate generation
- Conical and bell divergence efficiency $\lambda$
- Flow separation criteria (Summerfield, Schmucker)
