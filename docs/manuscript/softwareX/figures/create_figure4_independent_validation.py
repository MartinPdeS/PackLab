"""Create the independent Percus--Yevick validation figure.

The upper panels reproduce the monodisperse and binary-mixture PY presentations
in Figures 8.1.3 and 8.3.4 of Tsang et al. (2001) using digitised curves. The
lower panel compares reciprocal-space partial structure factors with mctpy.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from MPSPlots.styles import scientific

from PackLab import analytical, ureg

from manuscript_style import LEGEND_FONT_SIZE, MANUSCRIPT_WIDTH, add_panel_label, style_axes


OUTPUT_DIRECTORY = Path(__file__).parent
REPOSITORY_ROOT = OUTPUT_DIRECTORY.parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT))

from development.compare_mctpy import MIXTURES, compare_mixture  # noqa: E402


DATA_PATH = OUTPUT_DIRECTORY.parent / "data" / "tsang_2001_figure_8_3_4_digitized.csv"
MONODISPERSE_DATA_PATH = OUTPUT_DIRECTORY.parent / "data" / "tsang_2001_figure_8_1_3_digitized.csv"


def monodisperse_py_solution(separation: np.ndarray, volume_fraction: float):
    """Evaluate the monodisperse PY reference in the book coordinate r/b."""
    radii = np.array([1.0]) * ureg.micrometer
    domain = analytical.PercusYevickDomain(
        size=100_000 * ureg.micrometer,
        radii=radii,
        volume_fraction=volume_fraction,
        number_fractions=np.array([1.0]),
    )
    return analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber="auto",
    ).compute(separation * 2.0 * radii[0])


def binary_mixture_py_solution(separation: np.ndarray):
    """Evaluate PY at the Figure 8.3.4 radii and partial volume fractions."""
    radii = np.array([1.0, 2.0]) * ureg.micrometer
    partial_volume_fractions = np.array([0.2, 0.1])
    number_fractions = partial_volume_fractions / radii.to("micrometer").magnitude ** 3
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
    """Plot PackLab PY curves against digitised Figure 8.1.3 curves."""
    separation = np.linspace(0.0, 5.0, 600)
    digitized = np.genfromtxt(
        MONODISPERSE_DATA_PATH, delimiter=",", names=True, dtype=None, encoding="utf-8"
    )
    curves = ((0.2, "#1677b8", "o"), (0.3, "#d97706", "s"))
    for volume_fraction, colour, marker in curves:
        py_result = monodisperse_py_solution(separation, volume_fraction)
        reference = digitized[np.isclose(digitized["volume_fraction"], volume_fraction)]
        axes.plot(
            reference["reduced_separation"],
            reference["g_r"],
            linestyle="--",
            color="#111827",
            label="_nolegend_",
            zorder=5,
        )
        axes.plot(
            separation,
            py_result.g[0, 0],
            color=colour,
            linestyle="-",
            linewidth=4.2,
            label="_nolegend_",
            zorder=4,
        )

    axes.axvline(1.0, color="#64748b", linestyle="--", linewidth=1.0)
    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 3.5),
        ylim=(0.6, 2.8),
        xlabel=r"reduced separation $r/b$",
        ylabel=r"pair correlation $g(r)$",
    )
    axes.legend(
        handles=[
            Line2D(
                [],
                [],
                color=colour,
                linewidth=4.2,
                linestyle="-",
                label=rf"$f={volume_fraction:.1f}$",
            )
            for volume_fraction, colour, marker in curves
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper right",
    )
    style_axes(axes)


def plot_binary_mixture(axes: plt.Axes) -> None:
    """Plot PackLab PY curves against digitised Figure 8.3.4 curves."""
    separation = np.linspace(0.0, 5.0, 280)
    py_result = binary_mixture_py_solution(separation)
    digitized = np.genfromtxt(DATA_PATH, delimiter=",", names=True, dtype=None, encoding="utf-8")
    curves = (
        (r"$g_{11}$", "g11", (0, 0), "#1677b8", "^"),
        (r"$g_{12}$", "g12", (0, 1), "#d97706", "s"),
        (r"$g_{22}$", "g22", (1, 1), "#16a34a", "o"),
    )
    for label, component_pair, indices, colour, marker in curves:
        data = digitized[digitized["component_pair"] == component_pair]
        axes.plot(
            data["reduced_separation"],
            data["g_ij"],
            linestyle="--",
            color="#111827",
            label="_nolegend_",
            zorder=5,
        )
        axes.plot(
            separation,
            py_result.g[indices],
            color=colour,
            linestyle="-",
            linewidth=4.2,
            label="_nolegend_",
            zorder=4,
        )

    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 5.0),
        ylim=(0.6, 3.5),
        xlabel=r"reduced separation $r/(2a_1)$",
        ylabel=r"partial correlation $g_{ij}(r)$",
    )
    axes.legend(
        handles=[
            Line2D(
                [],
                [],
                color=colour,
                linewidth=4.2,
                linestyle="-",
                label=label,
            )
            for label, _, _, colour, marker in curves
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE - 1,
        loc="upper right",
        ncol=3,
    )
    style_axes(axes)


def mctpy_binary_comparison() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return the shared-grid binary structure factors and maximum error."""
    mixture = next(item for item in MIXTURES if item.name == "binary")
    wavenumber, packlab, mctpy, maximum_absolute, _ = compare_mixture(mixture)
    return wavenumber, packlab, mctpy, maximum_absolute


def plot_mctpy_structure_factors(
    axes: plt.Axes,
    wavenumber: np.ndarray,
    packlab: np.ndarray,
    mctpy: np.ndarray,
) -> None:
    """Plot binary partial structure factors from both implementations."""
    curves = (
        (r"$S_{11}$", (0, 0), "#1677b8"),
        (r"$S_{12}$", (0, 1), "#d97706"),
        (r"$S_{22}$", (1, 1), "#16a34a"),
    )
    solid_linewidth = 4.2
    dashed_linewidth = 1.5
    for _, indices, colour in curves:
        axes.plot(
            wavenumber,
            packlab[indices],
            color=colour,
            linewidth=solid_linewidth,
            solid_capstyle="round",
            zorder=2,
        )
        axes.plot(
            wavenumber,
            mctpy[indices],
            color="#111827",
            linestyle="--",
            linewidth=dashed_linewidth,
            dash_capstyle="butt",
            zorder=3,
        )

    axes.set(
        xlim=(0.0, 30.0),
        xlabel=r"wavenumber $k$ ($\mathrm{\mu m}^{-1}$)",
        ylabel=r"partial structure factor $S_{ij}(k)$",
    )
    axes.legend(
        handles=[
            *(
                Line2D([], [], color=colour, linewidth=solid_linewidth, label=rf"PackLab {label}")
                for label, _, colour in curves
            ),
            Line2D(
                [], [], color="#111827", linestyle="--", linewidth=dashed_linewidth, label="mctpy"
            ),
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE - 1,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=4,
    )
    style_axes(axes)


def plot(*, show: bool = True) -> plt.Figure:
    """Return the three-panel independent-validation layout."""
    figure = plt.figure(figsize=(MANUSCRIPT_WIDTH, 6.35), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 1.05))
    axes = (
        figure.add_subplot(grid[0, 0]),
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[1, :]),
    )
    figure.set_facecolor("white")
    plot_monodisperse(axes[0])
    plot_binary_mixture(axes[1])
    wavenumber, packlab, mctpy, _ = mctpy_binary_comparison()
    plot_mctpy_structure_factors(axes[2], wavenumber, packlab, mctpy)
    for panel_axes, label, label_x in zip(axes, ("(a)", "(b)", "(c)"), (-0.30, -0.30, -0.12)):
        add_panel_label(panel_axes, label, x=label_x, y_offset=20 if label == "(c)" else 6)
    if show:
        plt.show()
    return figure


def main() -> None:
    """Write editable SVG, PNG, and PDF copies of Figure 4."""
    with plt.style.context(scientific):
        figure = plot(show=False)
        for suffix in ("png", "svg", "pdf"):
            figure.savefig(
                OUTPUT_DIRECTORY / f"figure4_independent_validation.{suffix}", bbox_inches="tight"
            )
    _, _, _, maximum_absolute = mctpy_binary_comparison()
    print(f"Binary mctpy comparison maximum absolute difference: {maximum_absolute:.5e}")


if __name__ == "__main__":
    main()
