"""Create rendered Figure 4: RSA--PY, Metropolis--PY, and grid checks.

RSA and Percus--Yevick use matched particle inputs but represent distinct
physics: RSA is irreversible and history-dependent, while Percus--Yevick is
an equilibrium hard-sphere reference. Metropolis sampling evolves independent
valid configurations toward the fixed-volume equilibrium ensemble. The final
panel is a separate numerical check of the analytical Fourier grid.
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
PARTICLE_RADII = (
    np.array([0.5 / 1.4, 1.0 / 1.4]) * ureg.micrometer
)  # Two hard-sphere radii shared by RSA and PY.
NUMBER_FRACTIONS = np.array([0.50, 0.50])  # Equal requested counts of the two radius classes.
VOLUME_FRACTION = 0.25  # Target total occupied volume fraction for both workflows.
RSA_BOX_LENGTH = 30.0 * ureg.micrometer  # Side length of the periodic RSA cube.
RSA_ENSEMBLE_SIZE = (
    50 * 2
)  # Independent RSA realisations used for the mean and one-standard-error band.
RSA_NUMBER_OF_BINS = 150 * 10  # Radial histogram bins used by each RSA pair-correlation estimate.

# Metropolis comparison controls. Independent RSA configurations are used only
# as valid starting points, and one post-burn-in sample is retained per chain.
MC_BOX_LENGTH = 15.0 * ureg.micrometer
MC_ENSEMBLE_SIZE = 12
MC_BURN_IN_SWEEPS = 600
MC_MAXIMUM_DISPLACEMENT = 0.12 * ureg.micrometer
MC_NUMBER_OF_BINS = 250
MC_MAXIMUM_PAIRS = 1_000_000

# Numerical-convergence controls.
MAXIMUM_DISTANCE = (
    10.0 * ureg.micrometer
)  # Largest separation included in the PY comparison and convergence test.
NUMBER_OF_RADIAL_POINTS = 1000  # Evaluation points on the common real-space distance grid.
SAMPLES_PER_OSCILLATION = (
    4,
    6,
    8,
    10,
    12,
    16,
    24,
    36,
    48,
    64,
    96,
)  # Resolutions tested in panel (c).
REFERENCE_SAMPLES_PER_OSCILLATION = 512  # Refined PY grid used only as the convergence reference.
WARNING_THRESHOLD = 8  # PackLab warns below this many samples per fastest sinc-kernel oscillation.
AUTOMATIC_GRID_DEFAULT = 12  # Samples per fastest oscillation selected by wavenumber="auto".
DEFAULT_TITLE = None  # No figure-level title; the manuscript caption supplies context.

PAIR_INDICES = (
    (0, 0),
    (0, 1),
    (1, 1),
)  # The three unique partial correlations of a binary mixture.
PAIR_LABELS = (
    r"$g_{00}(r)$",
    r"$g_{01}(r)$",
    r"$g_{11}(r)$",
)  # Legend labels for those partial correlations.
PAIR_COLOURS = (
    "#1677b8",
    "#d97706",
    "#16a34a",
)  # Consistent colours for the three partial correlations.


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

    result = (
        estimator.estimate(number_of_samples=RSA_ENSEMBLE_SIZE, progress=True),
        estimator.statistics,
    )
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


def make_metropolis_estimate():
    """Return independent post-burn-in Metropolis correlation estimates."""
    domain = monte_carlo.PackingDomain(
        MC_BOX_LENGTH,
        MC_BOX_LENGTH,
        MC_BOX_LENGTH,
        use_periodic_boundaries=True,
    )
    sampler = samplers.DiscreteRadiusSampler(
        radii=PARTICLE_RADII,
        weights=NUMBER_FRACTIONS,
    )
    correlation_samples = []
    midpoint_samples = []
    acceptance_rates = []
    packing_fractions = []
    centers = None

    for chain_index in range(MC_ENSEMBLE_SIZE):
        rsa_options = monte_carlo.RSAOptions()
        rsa_options.random_seed = 4100 + chain_index
        rsa_options.maximum_attempts = 250_000
        rsa_options.maximum_consecutive_rejections = 50_000
        rsa_options.target_packing_fraction = VOLUME_FRACTION
        rsa_options.enforce_radii_distribution = True
        initial_result = monte_carlo.RSASimulator(domain, sampler, rsa_options).run()

        metropolis_options = monte_carlo.MetropolisOptions()
        metropolis_options.random_seed = 5100 + chain_index
        metropolis_options.number_of_sweeps = MC_BURN_IN_SWEEPS // 2
        metropolis_options.maximum_displacement = MC_MAXIMUM_DISPLACEMENT
        simulator = monte_carlo.MetropolisSimulator(
            domain,
            initial_result.sphere_configuration,
            metropolis_options,
        )
        midpoint_sample = simulator.run()
        _, midpoint_g = midpoint_sample.compute_partial_pair_correlation_function(
            n_bins=MC_NUMBER_OF_BINS,
            maximum_pairs=MC_MAXIMUM_PAIRS,
        )
        midpoint_samples.append(np.asarray(midpoint_g))

        metropolis_options.number_of_sweeps = MC_BURN_IN_SWEEPS // 2
        sample = simulator.run()
        sample_centers, partial_g = sample.compute_partial_pair_correlation_function(
            n_bins=MC_NUMBER_OF_BINS,
            maximum_pairs=MC_MAXIMUM_PAIRS,
        )
        if centers is None:
            centers = sample_centers
        correlation_samples.append(np.asarray(partial_g))
        acceptance_rates.append(simulator.statistics.acceptance_rate)
        packing_fractions.append(initial_result.statistics.packing_fraction_geometry)

    samples = np.asarray(correlation_samples)
    midpoint_samples = np.asarray(midpoint_samples)
    midpoint_mean = midpoint_samples.mean(axis=0)
    mean_g = samples.mean(axis=0)
    midpoint_standard_error = midpoint_samples.std(axis=0, ddof=1) / np.sqrt(MC_ENSEMBLE_SIZE)
    standard_error = samples.std(axis=0, ddof=1) / np.sqrt(MC_ENSEMBLE_SIZE)
    burn_in_change = np.sqrt(np.mean((mean_g - midpoint_mean) ** 2))
    combined_standard_error = np.sqrt(midpoint_standard_error**2 + standard_error**2)
    resolved = combined_standard_error > np.finfo(float).eps
    stability_fraction = np.mean(
        np.abs(mean_g[resolved] - midpoint_mean[resolved])
        <= 1.96 * combined_standard_error[resolved]
    )
    print(
        "Metropolis validation: "
        f"{MC_ENSEMBLE_SIZE} independent chains, "
        f"{MC_BURN_IN_SWEEPS} burn-in sweeps, "
        f"mean acceptance={np.mean(acceptance_rates):.3f}, "
        f"mean packing fraction={np.mean(packing_fractions):.6f}, "
        f"{MC_BURN_IN_SWEEPS // 2}--{MC_BURN_IN_SWEEPS}-sweep "
        f"RMS change={burn_in_change:.5f}, "
        f"within combined 95% uncertainty={stability_fraction:.3f}"
    )
    return (
        centers,
        mean_g,
        standard_error,
        np.mean(acceptance_rates),
        burn_in_change,
        stability_fraction,
    )


def post_contact_rms(centers, observed: np.ndarray, reference) -> tuple[float, ...]:
    """Return per-pair RMS differences beyond hard-core contact and below 5 µm."""
    center_values = centers.to("micrometer").magnitude
    errors = []
    for i, j in PAIR_INDICES:
        contact = (PARTICLE_RADII[i] + PARTICLE_RADII[j]).to("micrometer").magnitude
        mask = (center_values >= contact) & (center_values <= 5.0)
        errors.append(np.sqrt(np.mean((observed[i, j, mask] - reference.g[i, j, mask]) ** 2)))
    return tuple(errors)


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
    distances = (
        np.linspace(
            0.0,
            MAXIMUM_DISTANCE.to("micrometer").magnitude,
            NUMBER_OF_RADIAL_POINTS,
        )
        * ureg.micrometer
    )
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
        axes.plot(centers, rsa_mean, color=colour, linewidth=4.2)
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
        ylim=(0.6, 2.3),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"$g_{ij}(r)$",
    )
    axes.legend(
        handles=[
            *(
                Line2D([], [], color=colour, linewidth=4.2, label=label)
                for label, colour in zip(PAIR_LABELS, PAIR_COLOURS, strict=True)
            ),
            Line2D([], [], color="#111827", linestyle="--", linewidth=1.5, label="PY"),
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper right",
        ncol=2,
    )
    style_axes(axes)


def plot_metropolis_comparison(
    axes: Axes,
    *,
    centers,
    mean_g: np.ndarray,
    standard_error: np.ndarray,
    py_result,
) -> None:
    """Compare independent post-burn-in Metropolis samples with PY."""
    center_values = centers.to("micrometer").magnitude
    for pair, label, colour in zip(PAIR_INDICES, PAIR_LABELS, PAIR_COLOURS, strict=True):
        i, j = pair
        axes.plot(center_values, mean_g[i, j], color=colour, linewidth=4.2)
        axes.fill_between(
            center_values,
            mean_g[i, j] - standard_error[i, j],
            mean_g[i, j] + standard_error[i, j],
            color=colour,
            alpha=0.24,
        )
        axes.plot(
            center_values,
            py_result.g[i, j],
            color="#111827",
            linestyle="--",
            linewidth=1.5,
            zorder=4,
        )
        contact = (PARTICLE_RADII[i] + PARTICLE_RADII[j]).to("micrometer").magnitude
        axes.axvline(contact, color=colour, linestyle=":", linewidth=1.0)

    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, 5.0),
        ylim=(0.6, 2.3),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"$g_{ij}(r)$",
    )
    axes.legend(
        handles=[
            *(
                Line2D([], [], color=colour, linewidth=4.2, label=label)
                for label, colour in zip(PAIR_LABELS, PAIR_COLOURS, strict=True)
            ),
            Line2D([], [], color="#111827", linestyle="--", linewidth=1.5, label="PY"),
        ],
        frameon=True,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper right",
        ncol=2,
    )
    style_axes(axes)


def plot_convergence(axes: Axes) -> None:
    """Plot PY numerical convergence, independently of the RSA comparison."""
    samples, errors = grid_convergence()
    axes.semilogy(samples, errors, "o-", color="#d97706", linewidth=1.8, markersize=4.5)
    axes.axvline(
        WARNING_THRESHOLD, color="#b91c1c", linestyle="--", linewidth=2.0, label="warning threshold"
    )
    axes.axvline(
        AUTOMATIC_GRID_DEFAULT,
        color="#1677b8",
        linestyle="--",
        linewidth=2.0,
        label="automatic default",
    )
    axes.set(
        xlabel="samples per sinc-kernel oscillation",
        ylabel=r"PY RMS difference in $g_{ij}(r)$",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    axes.grid(alpha=0.18, which="both")
    axes.tick_params(labelsize=TICK_LABEL_SIZE)
    axes.xaxis.label.set_size(AXIS_LABEL_SIZE)
    axes.yaxis.label.set_size(AXIS_LABEL_SIZE)


def plot(
    *,
    show: bool = True,
    figsize: tuple[float, float] = (MANUSCRIPT_WIDTH, 6.35),
    title: str | None = DEFAULT_TITLE,
) -> plt.Figure:
    """Create the three-panel model-comparison and convergence figure.

    Returns
    -------
    matplotlib.figure.Figure
        RSA--PY and Metropolis--PY comparisons with a separate PY
        wavenumber-grid convergence panel.
    """
    estimate, statistics = make_rsa_estimate()
    rsa_py_result = solve_py(estimate.centers)
    (
        mc_centers,
        mc_mean_g,
        mc_standard_error,
        mc_acceptance,
        burn_in_change,
        stability_fraction,
    ) = make_metropolis_estimate()
    mc_py_result = solve_py(mc_centers)
    mc_rms = post_contact_rms(mc_centers, mc_mean_g, mc_py_result)
    print(
        "Metropolis--PY post-contact RMS differences: "
        + ", ".join(f"{label}={error:.5f}" for label, error in zip(PAIR_LABELS, mc_rms))
    )
    print(
        f"Metropolis summary: acceptance={mc_acceptance:.3f}, "
        f"burn-in RMS change={burn_in_change:.5f}, "
        f"stability fraction={stability_fraction:.3f}"
    )

    figure = plt.figure(figsize=figsize, constrained_layout=True)
    figure.set_facecolor("white")
    grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 1.05))
    rsa_axes = figure.add_subplot(grid[0, 0])
    metropolis_axes = figure.add_subplot(grid[0, 1])
    convergence_axes = figure.add_subplot(grid[1, :])

    plot_pair_comparisons(
        rsa_axes,
        estimate=estimate,
        completed_samples=statistics.completed_samples,
        py_result=rsa_py_result,
    )
    add_panel_label(rsa_axes, "(a)", x=-0.30)

    plot_metropolis_comparison(
        metropolis_axes,
        centers=mc_centers,
        mean_g=mc_mean_g,
        standard_error=mc_standard_error,
        py_result=mc_py_result,
    )
    add_panel_label(metropolis_axes, "(b)", x=-0.30)

    plot_convergence(convergence_axes)
    add_panel_label(convergence_axes, "(c)", x=-0.12, y_offset=20)
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
