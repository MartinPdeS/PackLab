"""Create manuscript Figure 3: RSA--PY comparison and grid convergence.

RSA and Percus--Yevick use matched particle inputs but represent distinct
physics: RSA is irreversible and history-dependent, while Percus--Yevick is
an equilibrium hard-sphere reference. The comparison panels do not imply that
the curves should coincide; the final panel is a separate numerical check.
"""

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
import numpy as np
from MPSPlots.styles import scientific

from PackLab import analytical, monte_carlo, samplers, ureg
from manuscript_style import (
    AXIS_LABEL_SIZE,
    LEGEND_FONT_SIZE,
    MANUSCRIPT_WIDTH,
    TICK_LABEL_SIZE,
    add_panel_label,
    style_axes,
)
LEGEND_FONT_SIZE -= 1  # Slightly compact legend for the two-panel layout.
AXIS_LABEL_SIZE -= 2  # Slightly compact axis labels for the two-panel layout.

OUTPUT_DIRECTORY = Path(__file__).parent

# Matched binary-mixture inputs.
PARTICLE_RADII = np.array([0.5 / 1.4, 1.0 / 1.4]) * ureg.micrometer  # Two hard-sphere radii shared by RSA and PY.
NUMBER_FRACTIONS = np.array([0.50, 0.50])  # Equal requested counts of the two radius classes.
VOLUME_FRACTION = 0.25  # Target total occupied volume fraction for both workflows.
RSA_BOX_LENGTH = 30.0 * ureg.micrometer  # Side length of the periodic RSA cube.
RSA_ENSEMBLE_SIZE = 50 * 2  # Independent RSA realisations used for the mean and one-standard-error band.
RSA_NUMBER_OF_BINS = 150 * 10  # Radial histogram bins used by each RSA pair-correlation estimate.

# Numerical-convergence controls.
MAXIMUM_DISTANCE = 10.0 * ureg.micrometer  # Largest separation included in the PY comparison and convergence test.
NUMBER_OF_RADIAL_POINTS = 1000  # Evaluation points on the common real-space distance grid.
SAMPLES_PER_OSCILLATION = (4, 6, 8, 10, 12, 16, 24, 36, 48, 64, 96)  # Resolutions tested in panel (b).
REFERENCE_SAMPLES_PER_OSCILLATION = 512  # Refined PY grid used only as the convergence reference.
WARNING_THRESHOLD = 8  # PackLab warns below this many samples per fastest sinc-kernel oscillation.
AUTOMATIC_GRID_DEFAULT = 12  # Samples per fastest oscillation selected by wavenumber="auto".
DEFAULT_TITLE = None  # No figure-level title; the manuscript caption supplies context.

PAIR_INDICES = ((0, 0), (0, 1), (1, 1))  # The three unique partial correlations of a binary mixture.
PAIR_LABELS = (r"$g_{00}(r)$", r"$g_{01}(r)$", r"$g_{11}(r)$")  # Legend labels for those partial correlations.
PAIR_COLOURS = ("#1677b8", "#d97706", "#7c3aed")  # Consistent colours for the three partial correlations.



def make_rsa_estimate():
    """Return a deterministic RSA ensemble for the binary mixture."""
    domain = monte_carlo.PackingDomain(
        RSA_BOX_LENGTH,
        RSA_BOX_LENGTH,
        RSA_BOX_LENGTH,
        use_periodic_boundaries=True,
    )
    sampler = samplers.DiscreteRadiusSampler(
        radii=PARTICLE_RADII,
        weights=NUMBER_FRACTIONS,
    )
    options = monte_carlo.RSAOptions()
    options.random_seed = 2026
    options.maximum_attempts = 250_000
    options.maximum_consecutive_rejections = 50_000
    options.target_packing_fraction = VOLUME_FRACTION
    options.enforce_radii_distribution = True

    estimator = monte_carlo.PackingEstimator(
        domain,
        sampler,
        options,
        number_of_bins=RSA_NUMBER_OF_BINS,
    )

    result = estimator.estimate(number_of_samples=RSA_ENSEMBLE_SIZE, progress=True), estimator.statistics
    estimator.print_statistics()
    return result

def make_py_domain() -> analytical.PercusYevickDomain:
    """Return an equilibrium PY domain with the same mixture inputs as RSA."""
    return analytical.PercusYevickDomain(
        size=100_000.0 * ureg.micrometer,
        radii=PARTICLE_RADII,
        volume_fraction=VOLUME_FRACTION,
        number_fractions=NUMBER_FRACTIONS,
    )


def solve_py(distances, *, wavenumber="auto"):
    """Evaluate the matched PY mixture on ``distances``."""
    domain = make_py_domain()
    return analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber=wavenumber,
    ).compute(distances)


def grid_convergence() -> tuple[np.ndarray, np.ndarray]:
    """Return PY errors relative to an independently refined discretisation."""
    distances = np.linspace(
        0.0,
        MAXIMUM_DISTANCE.to("micrometer").magnitude,
        NUMBER_OF_RADIAL_POINTS,
    ) * ureg.micrometer
    reference_grid = analytical.make_wavenumber_grid(
        radial_resolution=PARTICLE_RADII.min() / 20,
        maximum_distance=distances.max(),
        samples_per_oscillation=REFERENCE_SAMPLES_PER_OSCILLATION,
    )
    reference = solve_py(distances, wavenumber=reference_grid)

    errors = []
    for samples in SAMPLES_PER_OSCILLATION:
        grid = analytical.make_wavenumber_grid(
            radial_resolution=PARTICLE_RADII.min() / 20,
            maximum_distance=distances.max(),
            samples_per_oscillation=samples,
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The wavenumber grid is likely too coarse.*",
                category=RuntimeWarning,
            )
            result = solve_py(distances, wavenumber=grid)
        errors.append(np.sqrt(np.mean((result.g - reference.g) ** 2)))
    return np.asarray(SAMPLES_PER_OSCILLATION), np.asarray(errors)


def plot_pair_comparisons(
    axes: Axes,
    *,
    estimate,
    completed_samples: int,
    py_result,
) -> None:
    """Compare all unique RSA partial correlations with their PY references."""
    centers = estimate.centers.to("micrometer").magnitude
    for pair, label, colour in zip(PAIR_INDICES, PAIR_LABELS, PAIR_COLOURS, strict=True):
        i, j = pair
        rsa_mean = np.asarray(estimate.mean_g)[i, j]
        rsa_standard_error = np.asarray(estimate.std_g)[i, j] / np.sqrt(completed_samples)
        axes.plot(centers, rsa_mean, color=colour, linewidth=1.8)
        axes.fill_between(
            centers,
            rsa_mean - rsa_standard_error,
            rsa_mean + rsa_standard_error,
            color=colour,
            alpha=0.24,
        )
        axes.plot(centers, py_result.g[i, j], color="#111827", linestyle="--", linewidth=1.5)
        contact = (PARTICLE_RADII[i] + PARTICLE_RADII[j]).to("micrometer").magnitude
        axes.axvline(contact, color=colour, linestyle=":", linewidth=1.0)

    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 5.0),
        ylim=(0.6, 2.1),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"$g_{ij}(r)$",
    )
    axes.legend(
        handles=[
            *(Line2D([], [], color=colour, linewidth=1.8, label=label) for label, colour in zip(PAIR_LABELS, PAIR_COLOURS, strict=True)),
            Line2D([], [], color="#111827", linestyle="--", linewidth=1.5, label="PY"),
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper right",
        ncol=4,
    )
    style_axes(axes)


def plot_convergence(axes: Axes) -> None:
    """Plot PY numerical convergence, independently of the RSA comparison."""
    samples, errors = grid_convergence()
    axes.semilogy(samples, errors, "o-", color="#d97706", linewidth=1.8, markersize=4.5)
    # axes.axvline(WARNING_THRESHOLD, color="#b91c1c", linestyle="--", linewidth=1.1, label="warning threshold")
    # axes.axvline(AUTOMATIC_GRID_DEFAULT, color="#1677b8", linestyle="--", linewidth=1.1, label="automatic default")
    axes.set(
        xlabel="samples per sinc-kernel\noscillation",
        ylabel=r"PY RMS difference in $g_{ij}(r)$",
    )
    # axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    axes.grid(alpha=0.18, which="both")
    axes.tick_params(labelsize=TICK_LABEL_SIZE)
    axes.xaxis.label.set_size(AXIS_LABEL_SIZE)
    axes.yaxis.label.set_size(AXIS_LABEL_SIZE)


def plot(
    *,
    show: bool = True,
    figsize: tuple[float, float] = (MANUSCRIPT_WIDTH, 4.38),
    title: str | None = DEFAULT_TITLE,
) -> plt.Figure:
    """Create Figure 3 with model comparison and numerical convergence.

    Returns
    -------
    matplotlib.figure.Figure
        A combined RSA--PY pair-correlation comparison and a separate PY
        wavenumber-grid convergence panel.
    """
    estimate, statistics = make_rsa_estimate()
    py_result = solve_py(estimate.centers)

    figure = plt.figure(figsize=figsize, constrained_layout=True)
    figure.set_facecolor("white")
    grid = figure.add_gridspec(1, 2, width_ratios=(1.5, 1.0))
    comparison_axes = figure.add_subplot(grid[0, 0])
    convergence_axes = figure.add_subplot(grid[0, 1])

    plot_pair_comparisons(
        comparison_axes,
        estimate=estimate,
        completed_samples=statistics.completed_samples,
        py_result=py_result,
    )
    add_panel_label(comparison_axes, "(a)", x=-0.12)

    plot_convergence(convergence_axes)
    add_panel_label(convergence_axes, "(b)", x=-0.12)
    if title is not None:
        figure.suptitle(title, fontsize=15, fontweight="bold", color="#102a43")
    if show:
        plt.show()
    return figure


def main() -> None:
    """Write Figure 3 as editable SVG and manuscript PDF files."""
    with plt.style.context(scientific):
        figure = plot(show=True)
        figure.savefig(OUTPUT_DIRECTORY / "figure3_validation.png", bbox_inches="tight")
        figure.savefig(OUTPUT_DIRECTORY / "figure3_validation.svg", bbox_inches="tight")
        figure.savefig(OUTPUT_DIRECTORY / "figure3_validation.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
