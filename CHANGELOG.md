# Changelog

All notable changes to the Kryptonis Propulsion open-source toolkit will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.0] — 2026-09-09 (Day 1: Combustor Preliminary Design)

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

---

### Upcoming
- **Day 2 (`v0.2.0`):** Nozzle Gas Dynamics & Rao Bell Contour
- **Day 3 (`v0.3.0`):** Bartz Convective Heat Transfer & Thermal Axial Profile
- **Day 4 (`v0.4.0`):** Regenerative Cooling Channels & Pressure Drop
