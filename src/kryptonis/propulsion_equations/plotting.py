"""
Matplotlib engine cross-section visualization.
===============================================
Draws a publication-quality 2D axial cross-section of the thrust chamber
and nozzle, with dimension annotations and section labels.

Requires matplotlib (optional dependency: ``pip install kryptonis-propulsion-equations[plot]``).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kryptonis.propulsion_equations.profile import ChamberProfile

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend by default
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import PolyCollection
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def plot_chamber_profile(
    profile: "ChamberProfile",
    *,
    title: str | None = None,
    show_dimensions: bool = True,
    show_wall_thickness: bool = True,
    wall_thickness_mm: float = 2.5,
    save_path: str | None = None,
    show: bool = True,
    dpi: int = 150,
    figsize: tuple[float, float] = (14, 6),
) -> None:
    """Draw a 2D axial cross-section of the engine.

    Parameters
    ----------
    profile
        A ``ChamberProfile`` from ``generate_chamber_profile()``.
    title
        Plot title. Auto-generated if None.
    show_dimensions
        Annotate key dimensions with arrows.
    show_wall_thickness
        Draw the outer wall contour (cosmetic, not structural FEA).
    wall_thickness_mm
        Visual wall thickness in mm (cosmetic only).
    save_path
        If set, save the figure to this file path (PNG, PDF, SVG).
    show
        If True, display the plot interactively (calls plt.show()).
    """
    if not HAS_MPL:
        raise ImportError(
            "matplotlib is required for plotting. Install with:\n"
            "  pip install kryptonis-propulsion-equations[plot]"
        )

    # Convert to mm for display
    x_mm = [xi * 1e3 for xi in profile.x]
    r_mm = [ri * 1e3 for ri in profile.r]

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    # Inner wall — upper half
    ax.plot(x_mm, r_mm, color="#1a5276", linewidth=2.0, label="Inner wall")
    # Inner wall — lower half (mirror)
    ax.plot(x_mm, [-r for r in r_mm], color="#1a5276", linewidth=2.0)

    # Outer wall (cosmetic)
    if show_wall_thickness:
        wt = wall_thickness_mm
        r_outer = [r + wt for r in r_mm]
        ax.plot(x_mm, r_outer, color="#5d6d7e", linewidth=1.2,
                linestyle="--", label=f"Outer wall (+{wt:.1f} mm)")
        ax.plot(x_mm, [-r for r in r_outer], color="#5d6d7e",
                linewidth=1.2, linestyle="--")

        # Fill wall region
        ax.fill_between(x_mm, r_mm, r_outer, alpha=0.15, color="#aab7b8")
        ax.fill_between(x_mm, [-r for r in r_mm], [-r for r in r_outer],
                        alpha=0.15, color="#aab7b8")

    # Fill combustion gas region (light orange)
    ax.fill_between(x_mm, r_mm, [-r for r in r_mm], alpha=0.08,
                     color="#e74c3c", label="Gas flow region")

    # Injector face
    R_c_mm = profile.chamber_radius_m * 1e3
    wt_draw = wall_thickness_mm if show_wall_thickness else 0
    ax.plot([0, 0], [-(R_c_mm + wt_draw), (R_c_mm + wt_draw)],
            color="#2c3e50", linewidth=3)

    # Nozzle exit plane
    x_exit_mm = profile.x_exit * 1e3
    R_e_mm = profile.exit_radius_m * 1e3
    ax.plot([x_exit_mm, x_exit_mm], [-(R_e_mm + wt_draw), -(R_e_mm)],
            color="#2c3e50", linewidth=2)
    ax.plot([x_exit_mm, x_exit_mm], [(R_e_mm), (R_e_mm + wt_draw)],
            color="#2c3e50", linewidth=2)

    # Centreline
    ax.axhline(0, color="#95a5a6", linewidth=0.5, linestyle="-.")

    # Section dividers
    x_cyl_mm = profile.x_cyl_end * 1e3
    x_thr_mm = profile.x_throat * 1e3
    for xv, label in [(x_cyl_mm, ""), (x_thr_mm, "Throat")]:
        ax.axvline(xv, color="#d5d8dc", linewidth=0.8, linestyle=":")

    # Section labels
    R_c = profile.chamber_radius_m * 1e3
    R_t = profile.throat_radius_m * 1e3
    y_label = max(r_mm) * 1.15
    ax.text(x_cyl_mm / 2, y_label, "Cylindrical\nBarrel",
            ha="center", va="bottom", fontsize=8, color="#2c3e50")
    ax.text((x_cyl_mm + x_thr_mm) / 2, y_label, "Convergent\nCone",
            ha="center", va="bottom", fontsize=8, color="#2c3e50")
    ax.text((x_thr_mm + x_exit_mm) / 2, y_label,
            f"{'Bell' if profile.nozzle_type == 'bell' else 'Conical'}\nNozzle",
            ha="center", va="bottom", fontsize=8, color="#2c3e50")

    # Dimension annotations
    if show_dimensions:
        ann_kw = dict(fontsize=7.5, color="#2874a6",
                      arrowprops=dict(arrowstyle="<->", color="#2874a6",
                                      lw=1.0))

        # Throat diameter
        ax.annotate("", xy=(x_thr_mm, R_t), xytext=(x_thr_mm, -R_t),
                     arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
        ax.text(x_thr_mm + 3, 0, f"Dt={R_t*2:.1f}",
                fontsize=7.5, color="#c0392b", va="center")

        # Chamber diameter
        ax.annotate("", xy=(x_cyl_mm * 0.3, R_c),
                     xytext=(x_cyl_mm * 0.3, -R_c),
                     arrowprops=dict(arrowstyle="<->", color="#27ae60", lw=1.2))
        ax.text(x_cyl_mm * 0.3 + 3, 0, f"Dc={R_c*2:.1f}",
                fontsize=7.5, color="#27ae60", va="center")

        # Exit diameter
        ax.annotate("", xy=(x_exit_mm - 2, R_e_mm),
                     xytext=(x_exit_mm - 2, -R_e_mm),
                     arrowprops=dict(arrowstyle="<->", color="#8e44ad", lw=1.0))
        ax.text(x_exit_mm - 15, -R_e_mm * 0.5,
                f"De={R_e_mm*2:.1f}", fontsize=7, color="#8e44ad", va="center")

        # Total length (along bottom)
        y_bot = -(max(r_mm) * 1.25)
        ax.annotate("", xy=(0, y_bot), xytext=(x_exit_mm, y_bot),
                     arrowprops=dict(arrowstyle="<->", color="#2c3e50", lw=1.0))
        ax.text(x_exit_mm / 2, y_bot - max(r_mm) * 0.08,
                f"Total Length = {x_exit_mm:.1f} mm",
                ha="center", va="top", fontsize=8, color="#2c3e50",
                fontweight="bold")

    # Title
    if title is None:
        title = (f"Kryptonis Chamber Profile  |  "
                 f"Dt={profile.throat_diameter_m*1e3:.1f} mm  "
                 f"Dc={profile.chamber_diameter_m*1e3:.1f} mm  "
                 f"eps={profile.expansion_ratio:.0f}  "
                 f"({profile.nozzle_type} nozzle)")
    ax.set_title(title, fontsize=11, fontweight="bold", pad=12)

    ax.set_xlabel("Axial Position [mm]", fontsize=10)
    ax.set_ylabel("Radius [mm]", fontsize=10)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"  [SAVED] {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig
