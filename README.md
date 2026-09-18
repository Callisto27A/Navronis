# Navronis (`navronis-propulsion`)

<div align="center">

[![Live Web App](https://img.shields.io/badge/Live%20Demo-Interactive%20Web%20App-blueviolet.svg?style=for-the-badge&logo=firefox)](https://callisto27a.github.io/Navronis/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-brightgreen.svg)](https://python.org)
[![Build & Test](https://img.shields.io/badge/Tests-77%20Passed-success.svg)](#running-tests--verification)
[![Engineering Pedigree](https://img.shields.io/badge/Physics-Strict%20Provenance-orange.svg)](#primary-literature--pedigree)
[![Release](https://img.shields.io/badge/Release-v0.4.0%20(Day%204)-blue.svg)](https://github.com/Callisto27A/Navronis/releases)

**An authority-controlled, first-principles liquid rocket preliminary propulsion analytical toolkit with primary literature provenance.**

[🚀 **Live Interactive App**](https://callisto27a.github.io/Navronis/) •
[What Is Launched](#what-is-currently-launched) •
[Step-by-Step Guide](#step-by-step-usage-instructions) •
[4 Canonical Injectors](#the-4-canonical-injector-families) •
[Day 3: Cooling Channels](#day-3-regenerative-cooling-channels) •
[Day 4: Supersonic Nozzle](#day-4-supersonic-nozzle-aerodynamics--altitude-engine) •
[Governing Equations](#governing-equations--literature-provenance) •
[Community Issues](#community-auditing--open-issues) •
[Architecture](#architecture--project-structure)

</div>

---

> 🚀 **Live Interactive Web Application:** [**callisto27a.github.io/Navronis**](https://callisto27a.github.io/Navronis/)
>
> Run combustor, injector, and regenerative cooling sizing, visualize real-time SVG CAD geometry, toggle propellants, inspect droplet Sauter Mean Diameter ($D_{32}$) calculations, and verify test matrices directly in your web browser with zero installation.

---

## Overview

Most preliminary rocket propulsion scripts in circulation rely on undocumented constants, uncalibrated combustion efficiencies, or misattributed empirical formulas.

**Navronis** is built on a strict mathematical foundation:
- **Zero Silent Defaults:** Inputs and intermediate results are strictly validated. Non-physical states (e.g. negative throat areas, sub-unity contraction ratios, negative injector pressure drops, or channels exceeding perimeter) fail explicitly rather than silently propagating errors.
- **Traceable Engineering Provenance:** Every computed parameter carries its governing equation, validity domain, and exact literature citation (NASA SP-125, NASA SP-194, NASA SP-8087, Bartz 1957, Gnielinski 1976, Haaland 1983, Lorenzetto-Lefebvre 1977, Bazarov-Yang 1998, Dressler 2000, Rupe 1956).
- **Self-Explaining Engine:** The calculation engine explains every derived value upon request via `.explain("parameter_name")`.
- **Zero Heavy CAD/CFD Dependencies:** Pure Python, NumPy, and SciPy core. Runs instantaneously on Linux, macOS, and Windows.

---

## What Is Currently Launched

| Component | Status | Governing Physics & Scope | Primary References |
|---|---|---|---|
| **Day 1: Combustor Assembly** | 🟢 **Released (`v0.1.0`)** | Choked sonic throat continuity, Humble contraction ratio, barrel length, stay time, Bartz throat convective heat flux, NASA SP-194 acoustic cavity buzz modes (1T, 1R, 1L), ASME Section VIII hoop wall thickness. | NASA SP-125, Bartz (1957), NASA SP-194, ASME Sec VIII |
| **Day 2: 4 Canonical Injectors** | 🟢 **Released (`v0.2.0`)** | Torricelli hydraulics, chugging feed-system acoustic decoupling ($\Delta P / P_c \ge 15\%$), element packaging pitch, momentum flux ratio $J$, swirl sheet breakup, pintle momentum deflection, droplet atomization Sauter Mean Diameter ($D_{32}$). | NASA SP-194, Yang et al. (2004), Bazarov & Yang (1998), Dressler (2000), Rupe (1956) |
| **Day 3: Regen Cooling Channels** | 🟢 **Released (`v0.4.0`)** | Milled rectangular channel aspect ratio, Haaland friction factor, Gnielinski turbulent forced convection, fin efficiency enhancement ($\eta_{fin}$), 3-resistance conjugate wall thermal equilibrium ($T_{wg}$, $T_{wc}$), Darcy-Weisbach $\Delta P$, bulk $\Delta T$, and thermal-structural yield margin. | NASA SP-8087, Gnielinski (1976), Haaland (1983), Bartz (1957), Incropera |
| **Live Browser Web Sizer** | 🟢 **Live** | Interactive client-side UI hosted on GitHub Pages. Real-time SVG CAD cross-section rendering with cooling channels, propellant presets, and live 35-test verification matrix. | Client-side pure HTML5/SVG/JS |
| **Community Issue Tracker** | 🟢 **Active** | 5 open peer-review issues on GitHub for literature verification, atomization benchmarking, and thermal models. | GitHub Issues #1–#5 |

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

## Day 3: Regenerative Cooling Channels & 1D Conjugate Thermal Marching

Day 3 introduces high-pressure milled cooling jacket design. Fuel circulating through milled rectangular passages absorbs heat conducted through the chamber liner:

```
               COOLANT MANIFOLD / CHANNELS (High Speed Coolant vc ~ 50 m/s)
         ┌──────────────────────────────────────────────────────────────────┐
         │   CHANNEL   │   RIB FIN   │   CHANNEL   │   RIB FIN   │  CHANNEL │
         │ (w_c x h_c) │   (t_fin)   │ (w_c x h_c) │   (t_fin)   │          │
         └───────▲────────────▲──────────────▲────────────▲─────────────▲───┘
                 │            │              │            │             │  qc = hc,eff (Twc - Tbulk)
   ══════════════╪════════════╪══════════════╪════════════╪═════════════╪═══════ Coolant Wall (Twc)
                 │            │  RADIAL FOURIER CONDUCTION: q = (Twg - Twc) / (tw / kw)
   ══════════════╪════════════╪══════════════╪════════════╪═════════════╪═══════ Hot-Gas Wall (Twg)
                 │            │              │            │             │  qg = hg (Taw - Twg)
               HOT COMBUSTION GAS & NOZZLE ACCELERATION (Flame Temp Tc ~ 3,400 K)
```

### Governing Equations:
1. **Perimeter Packing:** $w_c = \frac{\pi D_t - N \cdot t_{fin}}{N}$, Aspect Ratio $AR = \frac{h_c}{w_c}$, $D_h = \frac{2 w_c h_c}{w_c + h_c}$
2. **Coolant Convection:** $Nu = \frac{(f/8)(Re - 1000)Pr}{1 + 12.7\sqrt{f/8}(Pr^{2/3} - 1)}$ (Gnielinski 1976 / Filonenko)
3. **Rough Channel Friction:** $\frac{1}{\sqrt{f}} = -1.8 \log_{10}\left[\left(\frac{\epsilon/D_h}{3.7}\right)^{1.11} + \frac{6.9}{Re}\right]$ (Haaland 1983)
4. **Fin Efficiency:** $\eta_{fin} = \frac{\tanh(m h_c)}{m h_c}$, where $m = \sqrt{\frac{2 h_c}{k_w t_{fin}}}$
5. **Conjugate Wall Closure:** $q = \frac{T_{aw} - T_{bulk}}{\frac{1}{h_g} + \frac{t_w}{k_w} + \frac{1}{h_{c,eff}}}$, $T_{wg} = T_{aw} - \frac{q}{h_g}$, $T_{wc} = T_{bulk} + \frac{q}{h_{c,eff}}$
6. **Coolant Pressure Drop:** $\Delta P = f \left(\frac{L_{channel}}{D_h}\right) \left(\frac{\rho_c v_c^2}{2}\right)$
7. **Thermal Stress:** $\sigma_{th} = \frac{E \alpha (T_{wg} - T_{wc})}{2(1 - \nu)}$, Yield Margin $MS = \frac{S_y}{\sigma_{th}} - 1.0$

---

## Step-by-Step Usage Instructions

### Method 1: Using the Live Web Application (Zero Installation)

1. Open [**https://callisto27a.github.io/Navronis/**](https://callisto27a.github.io/Navronis/).
2. Tweak your engine requirements in real time:
   - **🔥 Combustor Inputs (Day 1):** Thrust, $P_c$, Propellants, $L^*$, Convergent Angle.
   - **💧 Injector Inputs (Day 2):** Architecture (Coaxial, Swirl, Pintle, Impinging), Elements, $\Delta P / P_c$.
   - **❄️ Regen Cooling Inputs (Day 3):** Channel Count $N$, Depth $h_c$, Fin Width $t_{fin}$, Wall $t_w$, Alloy (CuCrZr, GRCop-42).
3. Observe instantaneous recalculation and inspect the live **2D Dynamic CAD Assembly Cross-Section** showing the combustor barrel, throat, injector faceplate, and green cooling channels with flow vectors.

---

### Method 2: Command-Line Interface (CLI)

#### 1. Installation
```bash
git clone https://github.com/Callisto27A/Navronis.git
cd Navronis
python -m pip install -e .
```

#### 2. Combustor Sizing Command (Day 1)
```bash
navronis --subsystem chamber --thrust 30000 --pc 120 --propellants LOX/CH4
```

#### 3. Injector Sizing Command (Day 2)
```bash
navronis --subsystem injector --injector-type coaxial --thrust 30000 --pc 120 --propellants LOX/CH4 --elements 19
```

#### 4. Regenerative Cooling Jacket Sizing (Day 3)
```bash
navronis --subsystem cooling --thrust 30000 --pc 120 --propellants LOX/CH4 --channels 80 --channel-height 1.8 --fin-thickness 0.8
```

**CLI Output Example (Day 3):**
```text
==============================================================================
NAVRONIS PROPULSION -- REGENERATIVE COOLING REPORT (DAY 3)
Propellant: LOX/CH4 | Thrust: 30.0 kN | Pc: 120.0 bar | Liner: CuCrZr
==============================================================================
Coolant Mass Flow:         2.304 kg/s (CH4)
Channels Count:            80 milled channels
------------------------------------------------------------------------------
Throat Channel Width (wc): 0.875 mm
Channel Height (hc):       1.800 mm
Channel Aspect Ratio (AR): 2.06
Hydraulic Diameter (Dh):   1.177 mm
Liner Wall Thickness (tw): 1.500 mm
------------------------------------------------------------------------------
Coolant Velocity (vc):     43.3 m/s
Reynolds Number (Re):      195803 (Fully Turbulent)
Friction Factor (f):       0.0415 (Haaland 1983)
Coolant Base HTC (hc):     91055.1 W/m²-K (Gnielinski 1976)
Fin Efficiency:            59.8%
Enhanced Effective HTC:    164669.6 W/m²-K
------------------------------------------------------------------------------
Hot-Gas Wall Temp (Twg):   139.8 K (-133.3 °C)
Coolant Wall Temp (Twc):   131.2 K (-141.9 °C)
Peak Throat Heat Flux (q): 1.84 MW/m²
Coolant Delta P:           44.30 bar
Coolant Bulk Temp Rise:    5.5 K
Thermal Compressive Stress:13.1 MPa
Yield Safety Margin (MS):  20.37 (PASS)
Thermal/Boiling Margin:    PASS (Below thermal limit)
==============================================================================
```

---

### Method 3: Python API

```python
from kryptonis.propulsion_equations import (
    CombustorDesign,
    InjectorDesign,
    RegenCoolingJacket,
)

# 1. Size Combustor (Day 1)
combustor = CombustorDesign(
    thrust=30000.0,
    chamber_pressure=12.0e6,
    mixture_ratio=2.6,
    propellant="LOX/CH4",
).solve()

# 2. Size Injector Head (Day 2)
injector = InjectorDesign(
    injector_type="coaxial",
    chamber_pressure=12.0e6,
    mass_flow_ox=combustor.total_mass_flow * 2.6 / 3.6,
    mass_flow_fuel=combustor.total_mass_flow / 3.6,
    rho_ox=1141.0,
    rho_fuel=422.0,
    n_elements=19,
).solve()

# 3. Size Regenerative Cooling Jacket (Day 3)
cooling = RegenCoolingJacket(
    throat_diameter_m=combustor.throat_diameter,
    chamber_diameter_m=combustor.chamber_diameter,
    chamber_length_m=combustor.chamber_length,
    mass_flow_coolant_kg_s=combustor.total_mass_flow / 3.6,
    chamber_pressure_pa=12.0e6,
    gas_recovery_temp_k=3400.0,
    gas_throat_htc_w_m2k=12500.0,
    n_channels=80,
    channel_height_m=0.0018,
    fin_thickness_m=0.0008,
    liner_material="CuCrZr",
).solve()

print(f"Throat Width:  {cooling.channel_width_mm:.3f} mm")
print(f"Coolant Speed: {cooling.coolant_velocity_m_s:.1f} m/s")
print(f"Hot-Wall Twg:  {cooling.hot_gas_wall_temp_K:.1f} K")
print(f"Coolant ΔP:    {cooling.coolant_pressure_drop_bar:.2f} bar")

# Query Literature Citation
print(cooling.explain("hot_gas_wall_temp_K"))
```

---

## Running Tests & Verification

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all 35 automated unit tests
pytest tests/ -v
```

**Verification Suite Overview (35/35 Passed):**
- `test_chamber_and_combustion.py`: Sonic throat area, isentropic mass flow, Humble contraction ratio, and stay time.
- `test_bartz.py`: Bartz recovery temperature bounds, convective heat transfer coefficient, and $\sigma$ factor.
- `test_chamber_acoustics.py`: Exact Bessel roots for 1T, 1R, and 1L acoustic cavity frequencies (NASA SP-194).
- `test_injector.py`: Torricelli hydraulics, chugging decoupling stiffness ($\Delta P / P_c \ge 0.15$), coaxial, swirl, pintle, and doublet atomization.
- `test_coolant_and_conduction.py`: Fourier radial wall conduction, Haaland friction factor, and Gnielinski/Dittus-Boelter coolant heat transfer.
- `test_regen_channel.py` **(Day 3)**: Milled channel geometry, perimeter limits, Gnielinski HTC, Haaland friction, fin efficiency, 3-resistance conjugate closure, Darcy-Weisbach $\Delta P$, thermal stress, and provenance explainers.

---

## Community Auditing & Open Issues

| Issue | Category | Description | Status |
|---|---|---|---|
| [**#1**](https://github.com/Callisto27A/Navronis/issues/1) | **Chamber / Heat Transfer** | Primary Literature Review: Contraction Ratio & Bartz Heat Flux Boundary Limits | 🟢 Open |
| [**#2**](https://github.com/Callisto27A/Navronis/issues/2) | **Injectors / Atomization** | Benchmarking Droplet SMD ($D_{32}$) Models Against Experimental Cold-Flow Data | 🟢 Open |
| [**#3**](https://github.com/Callisto27A/Navronis/issues/3) | **Thermal / Cooling** | Component 3: Axial Bartz Thermal Marching & Regenerative Cooling Channels | 🟢 **Implemented (`v0.4.0`)** |
| [**#4**](https://github.com/Callisto27A/Navronis/issues/4) | **Propellants / Presets** | Good First Issue: Add Storable Monopropellants & Hypergolic Presets (N₂O₄/MMH, Hydrazine) | 🟢 Open |
| [**#5**](https://github.com/Callisto27A/Navronis/issues/5) | **Verification & Testing** | Equation Verification Matrix: Cross-Validation with CEA, NASA SP-125 Datasets, and Cold-Flow Benchmarks | 🟢 Open |

---

## Architecture & Project Structure

```text
Navronis/
├── src/
│   └── kryptonis/
│       └── propulsion_equations/
│           ├── __init__.py           # Unified top-level API (v0.4.0)
│           ├── combustor.py          # Day 1: CombustorDesign & CombustorResult
│           ├── injector.py           # Day 2: 4 Canonical Injector Families
│           ├── regen_channel.py      # Day 3: RegenCoolingJacket & Conjugate Wall Solver
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
│           └── cli.py                # Standalone terminal CLI (chamber, injector, cooling)
├── docs/
│   ├── index.html                    # Interactive client-side Web Application (Day 1, 2, 3)
│   ├── .nojekyll                     # GitHub Pages static asset bypass
│   └── assets/                       # Dimensioned plots and figures
├── examples/
│   ├── combustor_30kn.py             # 30 kN LOX/CH4 chamber sizing example
│   ├── 01_thrust_chamber_sizing.py   # SP-125 analytical chamber sizing
│   ├── 02_bartz_heat_flux.py         # Throat convective heat flux & conduction
│   ├── 03_regenerative_cooling.py    # Cooling channel friction & curvature
│   ├── 04_nozzle_divergence_and_separation.py # Gas dynamics & separation limits
│   ├── 05_injector_sizing.py         # Sizing the 4 canonical injector families
│   └── 06_regenerative_cooling_channels.py # Day 3 regenerative cooling jacket sizing
├── tests/                            # 77 automated unit tests (100% pass rate)
├── pyproject.toml                    # PEP 621 standard package build metadata
├── LICENSE                           # Apache-2.0 open-source license
└── README.md                         # Project documentation and engineering handbook
```

---

## License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](LICENSE) file for details.
