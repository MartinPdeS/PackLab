"""Create an independent numerical-validation figure for the manuscript."""

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from MPSPlots.styles import scientific

from PackLab import analytical, ureg
from manuscript_style import LEGEND_FONT_SIZE, MANUSCRIPT_WIDTH, style_axes


OUTPUT_DIRECTORY = Path(__file__).parent

# --- Binary hard-sphere mixture -------------------------------------------
PARTICLE_RADII = np.array([150.0, 300.0]) * ureg.nanometer
NUMBER_FRACTIONS = np.array([0.70, 0.30])
VOLUME_FRACTION = 0.20
MAXIMUM_DISTANCE = 3.0 * ureg.micrometer
NUMBER_OF_RADIAL_POINTS = 700

# --- Validation controls ---------------------------------------------------
SAMPLES_PER_OSCILLATION = (4, 6, 8, 10, 12, 16, 24, 36, 48, 64, 96, 128, 192, 256)
REFERENCE_SAMPLES_PER_OSCILLATION = 512
WARNING_THRESHOLD = 8
AUTOMATIC_GRID_DEFAULT = 12
DEFAULT_TITLE = ""

PAIR_LABELS = (r"$g_{00}$", r"$g_{01}$", r"$g_{11}$")
PAIR_INDICES = ((0, 0), (0, 1), (1, 1))
PAIR_COLOURS = ("#1677b8", "#d97706", "#7c3aed")


def make_mixture_domain() -> analytical.PercusYevickDomain:
    """Return the binary mixture used in every validation calculation."""
    return analytical.PercusYevickDomain(
        size=100.0 * ureg.micrometer,
        radii=PARTICLE_RADII,
        volume_fraction=VOLUME_FRACTION,
        number_fractions=NUMBER_FRACTIONS,
    )


def distances():
    """Return the common radial evaluation grid."""
    return (
        np.linspace(0.0, MAXIMUM_DISTANCE.to("micrometer").magnitude, NUMBER_OF_RADIAL_POINTS)
        * ureg.micrometer
    )


def solve(distances, *, wavenumber="auto"):
    """Evaluate the binary Percus--Yevick solution on a specified grid."""
    mixture = make_mixture_domain()
    return analytical.PercusYevickSolver(
        densities=mixture.particle_densities_per_radius,
        radii=mixture.radii,
        wavenumber=wavenumber,
    ).compute(distances)


def refined_result(radial_distances):
    """Return the unplotted high-resolution reference solution."""
    wavenumber = analytical.make_wavenumber_grid(
        radial_resolution=PARTICLE_RADII.min() / 20,
        maximum_distance=radial_distances.max(),
        samples_per_oscillation=REFERENCE_SAMPLES_PER_OSCILLATION,
    )
    return solve(radial_distances, wavenumber=wavenumber)


def grid_convergence(radial_distances, reference) -> tuple[np.ndarray, np.ndarray]:
    """Return mixture-wide RMS differences from ``reference``."""
    errors = []
    for samples in SAMPLES_PER_OSCILLATION:
        wavenumber = analytical.make_wavenumber_grid(
            radial_resolution=PARTICLE_RADII.min() / 20,
            maximum_distance=radial_distances.max(),
            samples_per_oscillation=samples,
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The wavenumber grid is likely too coarse.*",
                category=RuntimeWarning,
            )
            result = solve(radial_distances, wavenumber=wavenumber)
        errors.append(np.sqrt(np.mean((result.g - reference.g) ** 2)))
    return np.asarray(SAMPLES_PER_OSCILLATION), np.asarray(errors)


def plot_binary_mixture(axes: Axes, radial_distances, result) -> None:
    """Show the three unique partial correlations of the binary mixture."""
    distances_um = radial_distances.to("micrometer").magnitude
    for (i, j), label, colour in zip(PAIR_INDICES, PAIR_LABELS, PAIR_COLOURS, strict=True):
        axes.plot(distances_um, result.g[i, j], color=colour, linewidth=1.8, label=label)

    for contact, colour in zip(
        (2 * PARTICLE_RADII[0], PARTICLE_RADII.sum(), 2 * PARTICLE_RADII[1]),
        PAIR_COLOURS,
        strict=True,
    ):
        axes.axvline(contact.to("micrometer").magnitude, color=colour, linestyle=":", linewidth=1.0)

    axes.axhline(1.0, color="#94a3b8", linestyle="--", linewidth=0.9)
    axes.set(
        xlim=(0.0, MAXIMUM_DISTANCE.to("micrometer").magnitude),
        ylim=(0.0, 2.8),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"partial $g_{ij}(r)$",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right", ncols=3)
    style_axes(axes)


def plot_auto_residual(axes: Axes, radial_distances, automatic, reference) -> None:
    """Show a readable local error measure for the automatic solution."""
    distances_um = radial_distances.to("micrometer").magnitude
    bin_size = 10
    usable_size = (len(distances_um) // bin_size) * bin_size
    binned_distances = distances_um[:usable_size].reshape(-1, bin_size).mean(axis=1)
    for (i, j), label, colour in zip(PAIR_INDICES, PAIR_LABELS, PAIR_COLOURS, strict=True):
        difference = automatic.g[i, j] - reference.g[i, j]
        local_rms = np.sqrt(np.mean(difference[:usable_size].reshape(-1, bin_size) ** 2, axis=1))
        # The solver enforces the hard core exactly, so residuals there are
        # identically zero and cannot be represented on a logarithmic axis.
        nonzero = local_rms > 1e-14
        axes.semilogy(
            binned_distances[nonzero],
            local_rms[nonzero],
            color=colour,
            linewidth=1.5,
            label=label,
        )

    axes.set(
        xlim=(0.0, MAXIMUM_DISTANCE.to("micrometer").magnitude),
        xlabel=r"separation $r$ (µm)",
        ylabel=r"local RMS difference in $g_{ij}(r)$",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right", ncols=3)
    style_axes(axes)
    axes.grid(alpha=0.18, which="both")


def plot_convergence(axes: Axes, radial_distances, reference) -> None:
    """Show convergence to the independently specified refined grid."""
    samples, errors = grid_convergence(radial_distances, reference)
    axes.semilogy(samples, errors, "o-", color="#d97706", linewidth=1.8, markersize=4.5)
    axes.axvline(
        WARNING_THRESHOLD, color="#b91c1c", linestyle="--", linewidth=1.1, label="warning threshold"
    )
    axes.axvline(
        AUTOMATIC_GRID_DEFAULT,
        color="#1677b8",
        linestyle="--",
        linewidth=1.1,
        label="automatic default",
    )
    axes.set(
        xlabel="samples per sinc-kernel oscillation",
        ylabel=r"RMS difference in $g_{ij}(r)$",
    )
    axes.legend(frameon=True, fontsize=LEGEND_FONT_SIZE, loc="upper right")
    style_axes(axes)
    axes.grid(alpha=0.18, which="both")


def plot(
    *,
    show: bool = True,
    figsize: tuple[float, float] = (MANUSCRIPT_WIDTH, 3.38),
    title: str | None = DEFAULT_TITLE,
) -> plt.Figure:
    """Create the complete manuscript validation figure.

    The panels respectively demonstrate binary-mixture output, local agreement
    of the automatic grid with a refined reference, and global grid convergence.
    """
    radial_distances = distances()
    automatic = solve(radial_distances)
    reference = refined_result(radial_distances)

    figure, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    figure.set_facecolor("white")
    plot_binary_mixture(axes[0], radial_distances, automatic)
    plot_auto_residual(axes[1], radial_distances, automatic, reference)
    plot_convergence(axes[2], radial_distances, reference)
    if title is not None:
        figure.suptitle(title, fontsize=15, fontweight="bold", color="#102a43")
    if show:
        plt.show()
    return figure


def main() -> None:
    """Write the validation figure as editable SVG and manuscript PDF files."""
    with plt.style.context(scientific):
        figure = plot(show=True)
        figure.savefig(OUTPUT_DIRECTORY / "validation_figure.png", bbox_inches="tight")
        # figure.savefig(OUTPUT_DIRECTORY / "validation_figure.svg", bbox_inches="tight")
        # figure.savefig(OUTPUT_DIRECTORY / "validation_figure.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
