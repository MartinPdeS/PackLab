"""Create the reproducible RSA example used as manuscript Figure 2.

The figure is generated from a fixed-seed PackLab simulation. Set
``INCLUDE_COMPOSITION_PANEL`` to ``False`` for a compact two-panel version.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Circle
import numpy as np
from MPSPlots.styles import scientific

from PackLab import monte_carlo, samplers, ureg
from manuscript_style import (
    LEGEND_FONT_SIZE,
    MANUSCRIPT_WIDTH,
    TICK_LABEL_SIZE,
    add_panel_label,
    style_axes,
)


OUTPUT_DIRECTORY = Path(__file__).parent

# --- Simulation parameters -------------------------------------------------
SEED = 20260828
BOX_LENGTH = 10.0 * ureg.micrometer
RADII = np.array([100.0, 200.0]) * ureg.nanometer
NUMBER_FRACTIONS = np.array([0.65, 0.35])
TARGET_PACKING_FRACTION = 0.20

# --- Figure parameters -----------------------------------------------------
INCLUDE_COMPOSITION_PANEL = False
SLICE_THICKNESS = 0.55  # micrometres
PAIR_CORRELATION_BINS = 500
PAIR_CORRELATION_MAXIMUM_RADIUS = 2.0  # micrometres
ESTIMATOR_SAMPLES = 64
COLOURS = np.array(["#1677b8", "#d97706"])
DEFAULT_TITLE = "Periodic random sequential adsorption of a polydisperse hard-sphere mixture"
DEFAULT_TITLE = ""


def run_simulation():
    """Run the fixed-seed, periodic, two-radius RSA simulation."""
    domain = monte_carlo.PackingDomain(
        length_x=BOX_LENGTH,
        length_y=BOX_LENGTH,
        length_z=BOX_LENGTH,
        use_periodic_boundaries=True,
    )
    sampler = samplers.DiscreteRadiusSampler(radii=RADII, weights=NUMBER_FRACTIONS)
    options = monte_carlo.RSAOptions()
    options.random_seed = SEED
    options.maximum_attempts = 100_000
    options.maximum_consecutive_rejections = 20_000
    options.target_packing_fraction = TARGET_PACKING_FRACTION
    options.enforce_radii_distribution = True
    return monte_carlo.RSASimulator(domain, sampler, options).run()


def estimate_pair_correlation():
    """Estimate the small-sphere partial correlation over an RSA ensemble."""
    domain = monte_carlo.PackingDomain(
        length_x=BOX_LENGTH,
        length_y=BOX_LENGTH,
        length_z=BOX_LENGTH,
        use_periodic_boundaries=True,
    )
    sampler = samplers.DiscreteRadiusSampler(radii=RADII, weights=NUMBER_FRACTIONS)
    options = monte_carlo.RSAOptions()
    options.random_seed = SEED
    options.maximum_attempts = 100_000
    options.maximum_consecutive_rejections = 20_000
    options.target_packing_fraction = TARGET_PACKING_FRACTION
    options.enforce_radii_distribution = True
    estimator = monte_carlo.PackingEstimator(
        domain,
        sampler,
        options,
        number_of_bins=PAIR_CORRELATION_BINS,
    )
    estimate = estimator.estimate(number_of_samples=ESTIMATOR_SAMPLES)
    return estimate, estimator.statistics


def plot_cross_section(
    axes: Axes,
    positions_um: np.ndarray,
    radii_um: np.ndarray,
    class_index: np.ndarray,
    box_length_um: float,
) -> None:
    """Plot the physical intersections of spheres with the central z-slice."""
    slice_center = box_length_um / 2.0
    displacement = positions_um[:, 2] - slice_center
    displacement -= box_length_um * np.round(displacement / box_length_um)
    in_slice = np.abs(displacement) <= np.minimum(radii_um, SLICE_THICKNESS / 2.0)

    for position, radius, particle_class, offset in zip(
        positions_um[in_slice],
        radii_um[in_slice],
        class_index[in_slice],
        displacement[in_slice],
    ):
        cross_section_radius = np.sqrt(max(radius**2 - offset**2, 0.0))
        axes.add_patch(
            Circle(
                (position[0], position[1]),
                cross_section_radius,
                facecolor=COLOURS[particle_class],
                edgecolor="white",
                linewidth=0.45,
                alpha=0.83,
            )
        )

    axes.set(
        aspect="equal",
        xlim=(0, box_length_um),
        ylim=(0, box_length_um),
        xlabel=r"$x$ (µm)",
        ylabel=r"$y$ (µm)",
    )
    style_axes(axes)
    # axes.text(
    #     0.03,
    #     0.03,
    #     f"{in_slice.sum()} intersected spheres\nperiodic domain",
    #     transform=axes.transAxes,
    #     fontsize=10,
    #     color="#334e68",
    #     va="bottom",
    # )


def plot_composition(
    axes: Axes, class_index: np.ndarray, particle_count: int, packing_fraction: float
) -> None:
    """Compare the requested and accepted number fractions for each radius."""
    locations = np.arange(len(RADII))
    accepted_counts = np.bincount(class_index, minlength=len(RADII))
    axes.bar(
        locations - 0.18,
        NUMBER_FRACTIONS,
        width=0.36,
        color="#b7ddf6",
        label="requested",
    )
    axes.bar(
        locations + 0.18,
        accepted_counts / particle_count,
        width=0.36,
        color=COLOURS,
        label="accepted",
    )
    axes.set(
        xticks=locations,
        xticklabels=[f"{radius:.0f} nm" for radius in RADII.to("nanometer").magnitude],
        ylim=(0, 0.8),
        xlabel="sphere radius",
        ylabel="number fraction",
        title="Radius composition",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    axes.text(
        0.04,
        0.92,
        f"N = {particle_count}\n$\\phi$ = {packing_fraction:.3f}",
        transform=axes.transAxes,
        fontsize=11,
        color="#334e68",
        va="top",
    )


def plot_pair_correlation(axes: Axes, estimate, radii_um: np.ndarray) -> None:
    """Plot all ordered partial correlations with estimator uncertainty."""
    distances_um = estimate.centers.to("micrometer").magnitude
    mean_g = np.asarray(estimate.mean_g)
    std_g = np.asarray(estimate.std_g)

    pairs = ((0, 0), (0, 1), (1, 1))
    labels = (r"$g_{00}$", r"$g_{01}/g_{10}$", r"$g_{11}$")
    colours = ("#1677b8", "#d97706", "#7c3aed")
    line_styles = ("-", "-", "-")
    contact_distances_um = []

    for (i, j), label, colour, line_style in zip(pairs, labels, colours, line_styles, strict=True):
        axes.plot(
            distances_um,
            mean_g[i, j],
            color=colour,
            linestyle=line_style,
            linewidth=2.0,
            label=rf"{label} ($\pm 1\sigma$)",
        )
        axes.fill_between(
            distances_um,
            mean_g[i, j] - std_g[i, j],
            mean_g[i, j] + std_g[i, j],
            color=colour,
            alpha=0.14,
            label="_nolegend_",
        )
        contact_distances_um.append(radii_um[i] + radii_um[j])

    axes.axhline(1.0, color="#94a3b8", linewidth=1.1, linestyle="--")
    axes.axvspan(
        0.0,
        min(contact_distances_um),
        color="#dc2626",
        alpha=0.12,
        zorder=0,
        label="_nolegend_",
    )
    for index, contact_distance in enumerate(sorted(set(contact_distances_um))):
        axes.axvline(
            contact_distance,
            color="#dc2626",
            linestyle=":",
            linewidth=1.2,
            label="_nolegend_",
        )
    axes.set(
        xlim=(0, PAIR_CORRELATION_MAXIMUM_RADIUS / 2),
        ylim=(0, 1.8),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"$g_{ij}(r)$",
    )
    axes.grid(alpha=0.18)
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    style_axes(axes)


def add_rsa_statistics(axes: Axes, statistics) -> None:
    """Display compact statistics from the fixed-seed RSA realisation."""
    acceptance_rate = statistics.accepted_insertions / statistics.attempted_insertions
    axes.text(
        0.04,
        0.04,
        (
            f"$N = {statistics.sphere_count}$\n"
            f"$\\phi = {statistics.packing_fraction_geometry:.3f}$\n"
            f"acceptance = {acceptance_rate:.1%}"
        ),
        transform=axes.transAxes,
        fontsize=TICK_LABEL_SIZE,
        va="bottom",
        ha="left",
        color="#102a43",
        bbox={"facecolor": "white", "edgecolor": "#94a3b8", "alpha": 0.88, "pad": 2.0},
    )


def plot(
    *,
    show: bool = True,
    figsize: tuple[float, float] | None = None,
    title: str | None = DEFAULT_TITLE,
    include_composition_panel: bool = INCLUDE_COMPOSITION_PANEL,
) -> plt.Figure:
    """Create the Figure 2 RSA visualisation.

    Parameters
    ----------
    show : bool, default=True
        Display the figure after it has been created.
    figsize : tuple of float, optional
        Figure width and height in inches. By default, width scales with the
        number of displayed panels.
    title : str or None, default=DEFAULT_TITLE
        Figure-level title. Pass ``None`` to omit it.
    include_composition_panel : bool, default=INCLUDE_COMPOSITION_PANEL
        Include the requested-versus-accepted radius-composition panel.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure.
    """
    result = run_simulation()
    estimate, _ = estimate_pair_correlation()
    positions_um = np.asarray(result.positions.to("micrometer").magnitude)
    radii_um = np.asarray(result.radii.to("micrometer").magnitude)
    class_index = np.asarray(result.sphere_configuration.classes_index)
    box_length_um = BOX_LENGTH.to("micrometer").magnitude

    panel_count = 3 if include_composition_panel else 2
    width_ratios = (1.0, 1.0, 1.65) if include_composition_panel else (1.0, 1.8)
    figure, axes = plt.subplots(
        1,
        panel_count,
        figsize=figsize or (4.65 * panel_count, 4.45),
        constrained_layout=True,
        gridspec_kw={"width_ratios": width_ratios},
    )
    axes = np.atleast_1d(axes)
    slice_axes, correlation_axes = axes[0], axes[-1]
    figure.set_facecolor("white")

    plot_cross_section(slice_axes, positions_um, radii_um, class_index, box_length_um)
    add_rsa_statistics(slice_axes, result.statistics)
    add_panel_label(slice_axes, "(a)", x=-0.23)
    if include_composition_panel:
        plot_composition(
            axes[1],
            class_index,
            particle_count=len(radii_um),
            packing_fraction=result.statistics.packing_fraction_geometry,
        )
    plot_pair_correlation(correlation_axes, estimate, radii_um)
    add_panel_label(correlation_axes, "(b)", x=-0.23)

    if title is not None:
        figure.suptitle(title, fontsize=16, fontweight="bold", color="#102a43")

    if show:
        plt.show()

    return figure


def main() -> None:
    """Write the default Figure 2 as editable SVG and manuscript PDF files."""
    with plt.style.context(scientific):
        figure = plot(show=False, figsize=(MANUSCRIPT_WIDTH, 3.38))
        figure.savefig(OUTPUT_DIRECTORY / "figure2_rsa_example.svg", bbox_inches="tight")
        figure.savefig(OUTPUT_DIRECTORY / "figure2_rsa_example.pdf", bbox_inches="tight")
        figure.savefig(OUTPUT_DIRECTORY / "figure2_rsa_example.png", bbox_inches="tight")


if __name__ == "__main__":
    main()
