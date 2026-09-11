# Changelog

All notable changes to the Kryptonis Propulsion open-source toolkit will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.2.0] — 2026-09-11 (Injector Head & Atomization Elements)

### Added
- **`InjectorDesign` Facade API:** Unified multi-architecture entry point supporting shear coaxial, pintle, and impinging doublet configurations with strict validation.
- **Hydraulic Orifice Sizing & Chugging Decoupling:**
  - Standard incompressible orifice discharge equation ($A_o = \dot{m} / (C_d \sqrt{2 \rho \Delta P})$)
  - NASA SP-194 low-frequency chugging acoustic stability criterion ($\Delta P / P_c \ge 15\% - 25\%$)
- **Shear Coaxial Injector Elements (Yang 2004, Huzel & Huang):**
  - Annular gas / liquid post hydraulic sizing
  - Momentum flux ratio ($J = \frac{\rho_g V_g^2}{\rho_l V_l^2}$) screening against the stable combustion window ($2 \le J \le 20$)
  - Lorenzetto & Lefebvre (1977) Sauter Mean Diameter (SMD $D_{32}$) droplet atomization model
  - Recess ratio ($R_L$) calculation for internal mixing and flame anchoring
- **Pintle Injector Elements (Dressler 2000, Heister 2019):**
  - Total Momentum Ratio ($TMR = \frac{\dot{m}_{rad} V_{rad}}{\dot{m}_{ann} V_{ann}}$)
  - Resultant spray cone half-angle ($\beta = \arccos\left(\frac{1}{1 + TMR}\right)$)
  - Annular fuel gap and discrete radial LOX slot sizing
- **Unlike Impinging Doublet Elements (Rupe 1956, Ingebo 1958):**
  - Rupe momentum ratio balance for uniform mixture ratio distribution
  - Free jet length before impingement ($L_{jet} = d_o / \tan(\theta/2)$)
  - Ingebo high-velocity liquid jet breakup droplet SMD ($D_{32}$)
- **CLI Extension:** `navronis --subsystem injector` with `--injector-type`, `--elements`, and `--delta-p-ratio` flags.
- **Automated Verification:** 6 new pytest unit tests (26 total unit tests passing in under 1.5s).
- **Comprehensive Sizing Example:** `examples/05_injector_sizing.py` demonstrating coaxial, pintle, and impinging doublet sizing for a 30 kN Methalox engine.

---

## [0.1.0] — 2026-09-09 (Combustor Preliminary Design)

### Added
- **`CombustorDesign` High-Level API:** User-facing class taking thrust, chamber pressure, mixture ratio, and propellant pair to solve preliminary chamber geometry.
- **`CombustorResult` with Mathematical Provenance:** Every quantity exposes `.explain("<parameter>")` citing primary literature (NASA SP-125, Bartz 1957, NASA SP-194).
- **Throat & Chamber Sizing:** Exact choked throat area ($A_t$), diameter ($D_t$), contraction ratio ($\epsilon_c$), chamber diameter ($D_c$), volume ($V_c$), and cylindrical barrel length ($L_{cyl}$).
- **Preliminary Physical Screening:**
  - Peak convective heat flux via Bartz (1957)
  - 1T, 1R, and 1L cylindrical acoustic cavity mode frequencies via Bessel eigenvalues (NASA SP-194)
  - Thin-shell hoop stress wall thickness screening via ASME Section VIII
- **Visual Cross-Section Plotting:** `plot_chamber_profile()` generating publication-quality 2D axial contour plots with dimension callouts.
- **Multi-Format Export:** `.to_json()` and `.to_csv()` exporting dimensions for CAD tools (CadQuery, FreeCAD, OpenSCAD).
- **Command-Line Interface:** `kryptonis-chamber` console script providing zero-code terminal sizing.
- **Automated Verification:** 20 automated pytest unit tests passing in under 2 seconds.
