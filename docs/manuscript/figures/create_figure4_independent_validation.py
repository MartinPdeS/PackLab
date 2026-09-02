"""Create Figure 4: book-style PY pair-correlation comparisons.

The two panels follow Figures 8.1.3 and 8.3.4 in Tsang et al. (2001).  The
binary-mixture marker positions are manually digitised from Figure 8.3.4;
the monodisperse panel remains a qualitative layout reference.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from MPSPlots.styles import scientific

from PackLab import analytical, ureg

from manuscript_style import LEGEND_FONT_SIZE, MANUSCRIPT_WIDTH, add_panel_label, style_axes


OUTPUT_DIRECTORY = Path(__file__).parent
DATA_PATH = OUTPUT_DIRECTORY.parent / "data" / "tsang_2001_figure_8_3_4_digitized.csv"


def damped_pair_correlation(
    separation: np.ndarray,
    *,
    contact: float,
    amplitude: float,
    decay: float = 2.0,
) -> np.ndarray:
    """Return a synthetic hard-sphere-like pair correlation above contact."""
    distance_from_contact = separation - contact
    correlation = 1.0 + amplitude * np.exp(-decay * distance_from_contact) * np.cos(
        5.7 * distance_from_contact
    )
    return np.where(separation < contact, 0.0, correlation)


def binary_mixture_py_solution(separation: np.ndarray):
    """Evaluate PY at the Figure 8.3.4 radii and partial volume fractions."""
    radii = np.array([1.0, 2.0]) * ureg.micrometer
    partial_volume_fractions = np.array([0.2, 0.1])
    number_fractions = partial_volume_fractions / radii.to("micrometer").magnitude**3
    number_fractions /= number_fractions.sum()
    domain = analytical.PercusYevickDomain(
        size=100_000 * ureg.micrometer,
        radii=radii,
        volume_fraction=partial_volume_fractions.sum(),
        number_fractions=number_fractions,
    )
    return analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber="auto",
    ).compute(separation * 2.0 * radii[0])


def plot_monodisperse(axes: plt.Axes) -> None:
    """Plot the Figure 8.1.3-style monodisperse comparison."""
    separation = np.linspace(0.0, 5.0, 600)
    curves = (
        (0.2, "#1677b8", ":", 0.55),
        (0.4, "#d97706", "-", 2.20),
    )
    for volume_fraction, colour, line_style, amplitude in curves:
        correlation = damped_pair_correlation(
            separation,
            contact=1.0,
            amplitude=amplitude,
        )
        axes.plot(
            separation,
            correlation,
            color=colour,
            linestyle=line_style,
            linewidth=1.8,
            label=rf"$f={volume_fraction:.1f}$",
        )

    axes.axvline(1.0, color="#64748b", linestyle="--", linewidth=1.0)
    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 5.0),
        ylim=(0.0, 4.0),
        xlabel=r"reduced separation $r/b$",
        ylabel=r"pair correlation $g(r)$",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    style_axes(axes)


def plot_binary_mixture(axes: plt.Axes) -> None:
    """Plot PackLab PY curves against digitised Figure 8.3.4 markers."""
    separation = np.linspace(0.0, 5.0, 280)
    py_result = binary_mixture_py_solution(separation)
    digitized = np.genfromtxt(DATA_PATH, delimiter=",", names=True, dtype=None, encoding="utf-8")
    curves = (
        (r"$g_{11}$", "g11", (0, 0), "#1677b8", "-.", "^"),
        (r"$g_{12}$", "g12", (0, 1), "#d97706", "--", "s"),
        (r"$g_{22}$", "g22", (1, 1), "#7c3aed", "-", "o"),
    )
    for label, component_pair, indices, colour, line_style, marker in curves:
        axes.plot(
            separation,
            py_result.g[indices],
            color=colour,
            linestyle=line_style,
            linewidth=1.6,
            label="_nolegend_",
        )
        data = digitized[digitized["component_pair"] == component_pair]
        axes.plot(
            data["reduced_separation"],
            data["g_ij"],
            linestyle="none",
            marker=marker,
            markersize=3.8,
            markerfacecolor="white",
            markeredgecolor="#111827",
            markeredgewidth=0.9,
            label="_nolegend_",
        )

    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 5.0),
        ylim=(0.0, 3.5),
        xlabel=r"reduced separation $r/(2a_1)$",
        ylabel=r"partial correlation $g_{ij}(r)$",
    )
    axes.legend(
        handles=[
            Line2D(
                [],
                [],
                color=colour,
                linestyle=line_style,
                marker=marker,
                markersize=4.0,
                markerfacecolor="white",
                markeredgecolor="#111827",
                label=label,
            )
            for label, _, _, colour, line_style, marker in curves
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE - 1,
        loc="upper right",
        ncol=3,
    )
    style_axes(axes)


def plot(*, show: bool = True) -> plt.Figure:
    """Return the two-panel Figure 4 comparison layout."""
    figure, axes = plt.subplots(1, 2, figsize=(MANUSCRIPT_WIDTH, 3.38), constrained_layout=True)
    figure.set_facecolor("white")
    plot_monodisperse(axes[0])
    plot_binary_mixture(axes[1])
    add_panel_label(axes[0], "(a)", x=-0.18)
    add_panel_label(axes[1], "(b)", x=-0.18)
    if show:
        plt.show()
    return figure


def main() -> None:
    """Write editable SVG, PNG, and PDF copies of Figure 4."""
    with plt.style.context(scientific):
        figure = plot(show=False)
        for suffix in ("png", "svg", "pdf"):
            figure.savefig(OUTPUT_DIRECTORY / f"figure4_independent_validation.{suffix}", bbox_inches="tight")


if __name__ == "__main__":
    main()
