# Navronis (`navronis-propulsion`)

<div align="center">

[![Live Web App](https://img.shields.io/badge/Live%20Demo-Interactive%20Web%20App-blueviolet.svg?style=for-the-badge&logo=firefox)](https://callisto27a.github.io/Navronis/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-brightgreen.svg)](https://python.org)
[![Build & Test](https://img.shields.io/badge/Tests-27%20Passed-success.svg)](#running-tests--verification)
[![Engineering Pedigree](https://img.shields.io/badge/Physics-Strict%20Provenance-orange.svg)](#primary-literature--pedigree)
[![Community Issues](https://img.shields.io/badge/Community%20Audits-5%20Active%20Issues-informational.svg)](https://github.com/Callisto27A/Navronis/issues)

**An authority-controlled, first-principles liquid rocket preliminary propulsion analytical toolkit with primary literature provenance.**

[🚀 **Live Interactive App**](https://callisto27a.github.io/Navronis/) •
[What Is Launched](#what-is-currently-launched) •
[Step-by-Step Guide](#step-by-step-usage-instructions) •
[4 Canonical Injectors](#the-4-canonical-injector-families) •
[Governing Equations](#governing-equations--literature-provenance) •
[Community Issues](#community-auditing--open-issues) •
[Architecture](#architecture--project-structure)

</div>

---

> 🚀 **Live Interactive Web Application:** [**callisto27a.github.io/Navronis**](https://callisto27a.github.io/Navronis/)
>
> Run combustor and injector sizing, visualize real-time SVG CAD geometry, toggle propellants, inspect droplet Sauter Mean Diameter ($D_{32}$) calculations, and verify test matrices directly in your web browser with zero installation.

---

## Overview

Most preliminary rocket propulsion scripts in circulation rely on undocumented constants, uncalibrated combustion efficiencies, or misattributed empirical formulas.

**Navronis** is built on a strict mathematical foundation:
- **Zero Silent Defaults:** Inputs and intermediate results are strictly validated. Non-physical states (e.g. negative throat areas, sub-unity contraction ratios, negative injector pressure drops) fail explicitly rather than silently propagating errors.
- **Traceable Engineering Provenance:** Every computed parameter carries its governing equation, validity domain, and exact literature citation (NASA SP-125, NASA SP-194, Bartz 1957, Lorenzetto-Lefebvre 1977, Bazarov-Yang 1998, Dressler 2000, Rupe 1956).
- **Self-Explaining Engine:** The calculation engine can explain every derived value upon request via `.explain("parameter_name")`.
- **Zero Heavy CAD/CFD Dependencies:** Pure Python, NumPy, and SciPy core. Runs instantaneously on Linux, macOS, and Windows.

---

## What Is Currently Launched

| Component | Status | Governing Physics & Scope | Primary References |
|---|---|---|---|
| **Day 1: Combustor Assembly** | 🟢 **Released** | Choked sonic throat continuity, Humble contraction ratio, barrel length, stay time, Bartz throat convective heat flux, NASA SP-194 acoustic cavity buzz modes (1T, 1R, 1L), ASME Section VIII hoop wall thickness. | NASA SP-125, Bartz (1957), NASA SP-194, ASME Sec VIII |
| **Day 2: 4 Canonical Injector Families** | 🟢 **Released** | Torricelli hydraulics, chugging feed-system acoustic decoupling ($\Delta P / P_c \ge 15\%$), element packaging pitch, momentum flux ratio $J$, swirl sheet breakup, pintle momentum deflection, droplet atomization Sauter Mean Diameter ($D_{32}$). | NASA SP-194, Yang et al. (2004), Bazarov & Yang (1998), Dressler (2000), Rupe (1956) |
| **Live Browser Web Sizer** | 🟢 **Live** | Interactive client-side UI hosted on GitHub Pages. Real-time SVG CAD cross-section rendering, propellant presets, instant recalculation, and live 27-test verification matrix. | Client-side pure HTML5/SVG/JS |
| **Community Issue Tracker** | 🟢 **Active** | 5 open peer-review issues on GitHub for literature verification, atomization benchmarking, and upcoming regenerative cooling channels. | GitHub Issues #1–#5 |

---

## The 4 Canonical Injector Families

Rather than attempting to model dozens of obscure injector variations, Navronis strictly focuses on the **4 canonical injector families** that power virtually all operational liquid rocket propulsion systems:

```
   1. SHEAR COAXIAL                2. SWIRL COAXIAL                3. CENTRAL PINTLE               4. IMPINGING DOUBLET
  (SSME, RL10, Raptor)           (RD-170, RD-180, NK-33)           (LMDE, Merlin 1D)             (Apollo SPS, Titan II)

      Gas Annulus                     Swirl Outer Gas               Fuel Annular Sheet               Fuel      Oxidizer
     ┌───────────┐                   ┌───────────────┐              ┌────────────────┐                \          /
     │   ┌───┐   │                   │   ┌───────┐   │              │   ┌────────┐   │                 \        /
     │   │Ox │   │                   │   │ Vortex│   │              │   │ Radial │   │                  \      /
     │   │Liq│   │                   │   │ Liquid│   │              │   │ Pintle │   │                   ▼    ▼
     │   └───┘   │                   │   └───────┘   │              │   └────────┘   │                 Collision Fan
     └───────────┘                   └───────────────┘              └────────────────┘                 (Atomization)
```

1. **Shear Coaxial Injectors (`coaxial`):** High-speed annular gas shears a central liquid core (e.g. SSME / RS-25, RL10, Vulcain, Raptor).
   - *Key physics:* Momentum flux ratio $J = \frac{\rho_g V_g^2}{\rho_l V_l^2}$ ($2 \le J \le 20$), Lorenzetto-Lefebvre (1977) droplet atomization ($D_{32}$), recess ratio $R_L$.
2. **Centrifugal Swirl Coaxial Injectors (`swirl`):** Tangential inlet ports create a rotating liquid film with a hollow central gas core, enveloped in an annular gas stream (e.g. RD-170, RD-180, NK-33).
   - *Key physics:* Abramovich centrifugal swirl discharge coefficient ($C_d$) from geometric swirl characteristic ($K$), hollow gas core diameter, annular liquid film thickness, spray half-cone angle $\theta$, and Lefebvre (1989) pressure-swirl droplet SMD ($D_{32}$).
3. **Central Pintle Injectors (`pintle`):** Continuous annular sheet deflected outward by an impinging central radial spray, providing wide throttling margins (e.g. Apollo LMDE, SpaceX Merlin 1D).
   - *Key physics:* Total Momentum Ratio $TMR = \frac{\dot{m}_{rad} V_{rad}}{\dot{m}_{ann} V_{ann}}$, resultant spray cone half-angle $\beta = \arccos\left(\frac{1}{1 + TMR}\right)$, acoustic decoupling stiffness.
4. **Unlike Impinging Doublet Injectors (`impinging`):** Discrete intersecting liquid jets forming an atomizing spray fan (e.g. Apollo SPS, Titan II, Viking).
   - *Key physics:* Rupe (1956) dynamic head balance $\frac{\rho_1 v_1^2 d_1}{\rho_2 v_2^2 d_2} = 1.0$, free jet impingement length, Ingebo (1958) droplet atomization.

---

## Step-by-Step Usage Instructions

### Method 1: Using the Live Web Application (Zero Installation)

The easiest way to explore and calculate rocket engine parameters is directly in your browser:

1. Open [**https://callisto27a.github.io/Navronis/**](https://callisto27a.github.io/Navronis/).
2. Select your subsystem tab:
   - **Combustor Assembly (Day 1):** Adjust Thrust, Chamber Pressure ($P_c$), Mixture Ratio ($O/F$), and Propellant combination. The 2D SVG chamber profile, sonic throat diameter, chamber barrel length, peak Bartz heat flux, and acoustic resonance buzz frequencies update immediately.
   - **4 Canonical Injectors (Day 2):** Select an injector architecture (`coaxial`, `swirl`, `pintle`, or `impinging`). Modify the element count and injector pressure drop ratio ($\Delta P / P_c$). The tool evaluates orifice geometries, velocities, momentum ratios, droplet Sauter Mean Diameters ($D_{32}$), and chugging decoupling status in real time.
3. Expand the **Test Suite Verification Drawer** at the bottom to inspect passing assertions across all 27 unit tests.

---

### Method 2: Command-Line Interface (CLI)

Navronis includes a fast, scriptable terminal CLI:

#### 1. Installation
```bash
git clone https://github.com/Callisto27A/Navronis.git
cd Navronis
python -m pip install -e .
```

#### 2. Combustor Sizing Command
To size a 30 kN Methalox combustor operating at 120 bar:
```bash
navronis --subsystem chamber --thrust 30000 --pc 120 --propellants LOX/CH4
```

**Where to enter custom values:**
- `--thrust <N>`: Vacuum or sea-level thrust in Newtons (e.g. `50000` for 50 kN).
- `--pc <bar>`: Chamber pressure in bar (e.g. `100` for 10 MPa).
- `--propellants <str>`: Supported presets: `LOX/CH4`, `LOX/RP1`, `LOX/LH2`.
- `--lstar <m>`: Characteristic chamber length $L^*$ (default: `1.0` m).
- `--angle <deg>`: Convergent half-angle (default: `30.0` degrees).

#### 3. Injector Sizing Commands (4 Canonical Families)
```bash
# Shear Coaxial (19 elements)
navronis --subsystem injector --injector-type coaxial --thrust 30000 --pc 120 --propellants LOX/CH4 --elements 19

# Centrifugal Swirl Coaxial (19 elements)
navronis --subsystem injector --injector-type swirl --thrust 30000 --pc 120 --propellants LOX/CH4 --elements 19

# Central Throttleable Pintle
navronis --subsystem injector --injector-type pintle --thrust 30000 --pc 120 --propellants LOX/CH4

# Unlike Impinging Doublets (16 pairs)
navronis --subsystem injector --injector-type impinging --thrust 30000 --pc 120 --propellants LOX/CH4 --elements 16
```

**CLI Output Example:**
```text
==============================================================================
NAVRONIS PROPULSION -- INJECTOR SIZING REPORT: COAXIAL
Propellant: LOX/CH4 | Thrust: 30.0 kN | Pc: 120.0 bar
==============================================================================
Mass Flow: Total = 10.370 kg/s (LOX: 8.066 kg/s, Fuel: 2.304 kg/s)
Injector Delta P: 24.00 bar (20.0% Pc)
Chugging Decoupling Margin: PASS (>=15%)
------------------------------------------------------------------------------
Elements:                  19
Liquid Post ID:            2.86 mm
Liquid Post OD:            3.86 mm
Gas Sleeve ID:             5.46 mm
Annular Gap:               0.80 mm
Liquid Ox Velocity:        65.86 m/s
Gas/Fuel Velocity:         104.99 m/s
Momentum Flux Ratio J:     9.35  [Target: 2.0 - 20.0]
Velocity Ratio VR:         1.59
Recess Length:             2.86 mm
Droplet SMD (D32):         24.3 µm
==============================================================================
```

---

### Method 3: Python API

For integration into automated preliminary design pipelines, trajectory optimizers, or parametric trade studies:

```python
from kryptonis.propulsion_equations import CombustorDesign, InjectorDesign

# -------------------------------------------------------------
# 1. Size the Thrust Chamber Assembly
# -------------------------------------------------------------
chamber = CombustorDesign(
    thrust=30000.0,            # 30,000 N (30 kN)
    chamber_pressure=12.0e6,    # 12.0 MPa (120 bar)
    mixture_ratio=2.6,         # O/F ratio
    propellant="LOX/CH4",       # Built-in thermochemical dataset
    convergent_half_angle_deg=30.0,
)
result = chamber.solve()

print(f"Throat Diameter:   {result.throat_diameter * 1000:.2f} mm")
print(f"Chamber Diameter:  {result.chamber_diameter * 1000:.2f} mm")
print(f"Chamber Length:    {result.chamber_length * 1000:.2f} mm")
print(f"Peak Heat Flux:    {result.heat_flux_screen / 1e6:.2f} MW/m²")
print(f"1T Acoustic Buzz:  {result.acoustic_modes['1T_Hz']:.1f} Hz")

# -------------------------------------------------------------
# 2. Size the Injector Head Subsystem
# -------------------------------------------------------------
injector = InjectorDesign(
    injector_type="coaxial",    # "coaxial", "swirl", "pintle", or "impinging"
    chamber_pressure=12.0e6,     # Pa
    mass_flow_ox=result.mass_flow_ox,      # Fed directly from combustor result
    mass_flow_fuel=result.mass_flow_fuel,  # Fed directly from combustor result
    rho_ox=1141.0,               # kg/m³ (LOX)
    rho_fuel=422.0,              # kg/m³ (LCH4)
    n_elements=19,               # Number of elements on faceplate
    delta_p_ratio=0.20,          # 20% Pc injector drop
)
inj_result = injector.solve()

print(f"Post ID:           {inj_result['post_id_mm']:.2f} mm")
print(f"Annular Gap:       {inj_result['annular_gap_mm']:.2f} mm")
print(f"Momentum Ratio J:  {inj_result['momentum_flux_ratio_J']:.2f}")
print(f"Droplet SMD D32:   {inj_result['smd_um']:.1f} µm")
print(f"Chugging Margin:   {inj_result['chugging_margin_adequate']}")

# -------------------------------------------------------------
# 3. Query Engineering Provenance & Equations
# -------------------------------------------------------------
print(result.explain("throat_diameter"))
```

---

## Running Tests & Verification

The test suite validates every analytical formula against independent mathematical limits, monotonicity checks, and published literature values.

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all 27 automated unit tests
pytest tests/ -v
```

**Verification Suite Overview:**
- `test_chamber_and_combustion.py`: Verifies Vandenkerckhove $\Gamma(\gamma)$, sonic throat area, isentropic mass flow, Humble contraction ratio, and stay time.
- `test_bartz.py`: Verifies Bartz recovery temperature bounds, convective heat transfer coefficient, and $\sigma$ boundary-layer correction factor.
- `test_chamber_acoustics.py`: Verifies exact Bessel roots for 1T, 1R, and 1L transverse and longitudinal acoustic cavity frequencies (NASA SP-194).
- `test_injector.py`: Verifies Torricelli hydraulics, chugging decoupling stiffness ($\Delta P / P_c \ge 0.15$), shear coaxial momentum flux ratios, centrifugal swirl core physics, pintle momentum balances, and impinging doublet atomization.
- `test_coolant_and_conduction.py`: Verifies Fourier radial wall conduction, Haaland friction factor, and Gnielinski/Dittus-Boelter coolant heat transfer.

---

## Governing Equations & Literature Provenance

| Category | Mathematical Formulation | Literature Reference | Output Parameters |
|---|---|---|---|
| **Choked Throat Flow** | $\dot{m} = \frac{P_c A_t}{c^*}$, $A_t = \frac{\dot{m} \sqrt{\gamma R T_0}}{P_c \Gamma(\gamma)}$ | NASA SP-125 Eq. (1-9) | Total $\dot{m}$, Throat Area $A_t$, Throat Diameter $D_t$ |
| **Chamber Diameter** | $\epsilon_c = \frac{A_c}{A_t} = 1.25 + 8.0 D_t^{-0.6}$ | Humble (1995) / NASA SP-125 | Contraction ratio $\epsilon_c$, Chamber Diameter $D_c$ |
| **Chamber Barrel & Volume** | $V_c = L^* A_t$, $L_{cyl} = \frac{V_c - V_{conv}}{A_c}$ | NASA SP-125 Eq. (4-4, 4-5) | Chamber Volume $V_c$, Barrel Length $L_{cyl}$, Total Length $L_c$ |
| **Stay Residence Time** | $	au_s = \frac{V_c ho_c}{\dot{m}}$ | NASA SP-125 p. 87 | Residence stay time $	au_s$ (ms) |
| **Throat Convective Heat Flux** | $h_g = \frac{0.026}{D_t^{0.2}} \left(\frac{\mu^{0.2} C_p}{Pr^{0.6}}\right)_0 \left(\frac{P_c}{c^*}\right)^{0.8} \left(\frac{D_t}{R_c}\right)^{0.1} \sigma$ | Bartz (1957) Jet Propulsion | Peak throat convective HTC $h_g$, Heat flux $q$ (MW/m²) |
| **Acoustic Cavity Buzz** | $f_{1T} = \frac{1.8412 a}{\pi D_c}$, $f_{1R} = \frac{3.8317 a}{\pi D_c}$, $f_{1L} = \frac{a}{2 L_c}$ | NASA SP-194 (Harrje & Reardon) | First Tangential (1T), Radial (1R), and Longitudinal (1L) frequencies (Hz) |
| **Mechanical Wall Thickness** | $t_w = \frac{P_c R_c}{\sigma_{yield} / SF - 0.6 P_c}$ | ASME Boiler & Pressure Vessel Sec VIII | Minimum structural wall thickness $t_w$ (mm) |
| **Injector Hydraulics** | $A_o = \frac{\dot{m}}{C_d \sqrt{2 ho \Delta P}}$, $\frac{\Delta P}{P_c} \ge 0.15$ | NASA SP-194 / Huzel & Huang | Orifice areas, injection velocities, feed-system chugging margin |
| **Shear Coaxial Atomization** | $J = \frac{ho_g v_g^2}{ho_l v_l^2}$, $D_{32} \propto \frac{\sigma_l^{0.5}}{\Delta v}$ | Yang et al. (2004) / Lorenzetto-Lefebvre | Post ID/OD, annular gap, momentum flux ratio $J$, droplet SMD $D_{32}$ |
| **Centrifugal Swirl Physics** | $K = \frac{\pi R_{in} r_o}{n A_p}$, $\phi = f(K)$, $C_d = \mu_{swirl}$ | Bazarov & Yang (1998) / Lefebvre | Swirl orifice, gas core diameter, film thickness, spray angle, SMD $D_{32}$ |
| **Pintle Spray Mechanics** | $TMR = \frac{\dot{m}_{rad} v_{rad}}{\dot{m}_{ann} v_{ann}}$, $eta = \arccos\left(\frac{1}{1 + TMR}\right)$ | Dressler (2000) / Heister (2019) | Pintle diameter, slot height, annular gap, $TMR$, spray cone angle $eta$ |
| **Unlike Impinging Doublet** | $\frac{ho_1 v_1^2 d_1}{ho_2 v_2^2 d_2} = 1.0$, $L_{free} = (3	ext{--}5) d_o$ | Rupe (1956) JPL / Ingebo (1958) | Orifice diameters, impingement point, droplet SMD $D_{32}$ |

---

## Community Auditing & Open Issues

To guarantee peer-reviewed rigor and transparency, Navronis hosts public auditing discussions on [GitHub Issues](https://github.com/Callisto27A/Navronis/issues). We invite propulsion engineers, aerospace researchers, and students to audit and contribute:

| Issue | Category | Description | Status |
|---|---|---|---|
| [**#1**](https://github.com/Callisto27A/Navronis/issues/1) | **Chamber / Heat Transfer** | Primary Literature Review: Contraction Ratio & Bartz Heat Flux Boundary Limits | 🟢 Open |
| [**#2**](https://github.com/Callisto27A/Navronis/issues/2) | **Injectors / Atomization** | Benchmarking Droplet SMD ($D_{32}$) Models Against Experimental Cold-Flow Data | 🟢 Open |
| [**#3**](https://github.com/Callisto27A/Navronis/issues/3) | **Thermal / Cooling** | Component 3: Axial Bartz Thermal Marching & Regenerative Cooling Channels | 🟢 Open |
| [**#4**](https://github.com/Callisto27A/Navronis/issues/4) | **Propellants / Presets** | Good First Issue: Add Storable Monopropellants & Hypergolic Presets (N₂O₄/MMH, Hydrazine) | 🟢 Open |
| [**#5**](https://github.com/Callisto27A/Navronis/issues/5) | **Verification & Testing** | Equation Verification Matrix: Cross-Validation with CEA, NASA SP-125 Datasets, and Cold-Flow Benchmarks | 🟢 Open |

### How to Contribute a Verification Test or Equation Correction
1. Fork the repository and create a branch (`feature/verify-equations`).
2. Add or modify tests in `tests/`.
3. Provide the citation (NASA SP, JANNAF, AIAA, ASME, or CEA output deck).
4. Submit a Pull Request. Every PR triggers automated unit tests across Python 3.10–3.12.

---

## CAD Export & Visualization

### 1. 2D Axial Contour Plotting
Generate publication-quality dimensioned cross-sections with cosmetic wall thickness:

```python
result.plot(save_path="combustor_profile.png")
```

<div align="center">
  <img src="docs/assets/combustor_30kn.png" alt="Navronis Combustor Profile Cross-Section" width="850">
</div>

### 2. Multi-Format CAD Export
Export coordinate arrays and design parameters for CAD packages:

```python
result.to_json("combustor.json")   # CadQuery, OpenSCAD, or FreeCAD ingestion
result.to_csv("combustor.csv")     # Spreadsheets, MATLAB, FreeCAD Draft Workbench
```

---

## Architecture & Project Structure

```text
Navronis/
├── src/
│   └── kryptonis/
│       └── propulsion_equations/
│           ├── __init__.py           # Unified top-level API
│           ├── combustor.py          # High-level CombustorDesign & CombustorResult
│           ├── injector.py           # 4 Canonical Injector Families & atomization physics
│           ├── chamber.py            # NASA SP-125 throat, contraction & volume equations
│           ├── combustion.py         # Characteristic velocity c* & stay time
│           ├── bartz.py              # Canonical Bartz 1957 convective film coefficient
│           ├── chamber_acoustics.py  # NASA SP-194 Bessel acoustic cavity modes
│           ├── profile.py            # Axial (x, r) coordinate inner-wall generator
│           ├── plotting.py           # Matplotlib 2D cross-section visualizer
│           ├── export.py             # Structured JSON & CSV profile exporters
│           ├── units.py              # Strongly-typed Result, Status & provenance containers
│           ├── aerodynamics.py       # 1D isentropic gas dynamics Area-Mach solver
│           ├── thermal.py            # Conjugate heat transfer & thermal marching
│           ├── coolant.py            # Supercritical fluid convection & channel friction
│           ├── wall_conduction.py    # 1D radial Fourier heat conduction
│           └── cli.py                # Standalone terminal CLI (navronis / kryptonis-chamber)
├── docs/
│   ├── index.html                    # Interactive client-side Web Application
│   ├── .nojekyll                     # GitHub Pages static asset bypass
│   └── assets/                       # Dimensioned plots and figures
├── examples/
│   ├── combustor_30kn.py             # 30 kN LOX/CH4 chamber sizing example
│   ├── 01_thrust_chamber_sizing.py   # SP-125 analytical chamber sizing
│   ├── 02_bartz_heat_flux.py         # Throat convective heat flux & conduction
│   ├── 03_regenerative_cooling.py    # Cooling channel friction & curvature
│   ├── 04_nozzle_divergence_and_separation.py # Gas dynamics & separation limits
│   └── 05_injector_sizing.py         # Sizing the 4 canonical injector families
├── tests/                            # 27 automated unit tests (100% pass rate)
├── pyproject.toml                    # PEP 621 standard package build metadata
├── LICENSE                           # Apache-2.0 open-source license
└── README.md                         # Project documentation and engineering handbook
```

---

## Scope & Limitations

> [!IMPORTANT]
> **Navronis is a preliminary engineering calculation toolkit.**
> - It does **NOT** design manufacturing-ready rocket engines.
> - It is **NOT** a substitute for 3D reacting CFD, conjugate heat transfer (CHT), or finite element structural analysis (FEA).
> - Acoustic modes are **frequencies only**, not a dynamic stability assessment (which requires hardware-measured injector response functions $n$ and $\tau$).

---

## License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](LICENSE) file for details.
