"""Generate the curated figures displayed in the project README."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

from PackLab import analytical, monte_carlo, samplers, ureg


OUTPUT_DIRECTORY = Path(__file__).parent


def create_rsa_packing_figure() -> None:
    """Save a deterministic two-dimensional slice through an RSA packing."""
    domain = monte_carlo.PackingDomain(
        5 * ureg.micrometer,
        5 * ureg.micrometer,
        5 * ureg.micrometer,
        use_periodic_boundaries=True,
    )
    sampler = samplers.UniformRadiusSampler(
        90 * ureg.nanometer,
        170 * ureg.nanometer,
        bins=8,
    )
    options = monte_carlo.RSAOptions()
    options.random_seed = 42
    options.maximum_attempts = 40_000
    options.maximum_consecutive_rejections = 8_000
    options.target_packing_fraction = 0.12

    result = monte_carlo.RSASimulator(domain, sampler, options).run()
    figure = result.plot_slice_2d(show=False, slice_thickness_fraction=0.12)
    axis = figure.axes[0]
    for circle in axis.patches:
        circle.set_edgecolor("#1864ab")
        circle.set_facecolor("#d0ebff")
        circle.set_alpha(0.75)

    micrometer_formatter = FuncFormatter(lambda value, _: f"{value / 1e-6:g}")
    axis.xaxis.set_major_formatter(micrometer_formatter)
    axis.yaxis.set_major_formatter(micrometer_formatter)
    axis.set_title("Periodic RSA packing")
    axis.set_xlabel(r"$x$ ($\mu$m)")
    axis.set_ylabel(r"$y$ ($\mu$m)")
    figure.set_size_inches(5.2, 4.4)
    figure.tight_layout()
    figure.savefig(OUTPUT_DIRECTORY / "readme_rsa_packing.png", dpi=180, bbox_inches="tight")
    plt.close(figure)


def create_percus_yevick_figure() -> None:
    """Save partial pair correlations for a binary PY hard-sphere mixture."""
    radii = np.array([75, 140]) * ureg.nanometer
    domain = analytical.PercusYevickDomain(
        size=50 * ureg.micrometer,
        radii=radii,
        volume_fraction=0.25,
        number_fractions=np.array([0.7, 0.3]),
    )
    distances = np.linspace(0.0, 1.5, 300) * ureg.micrometer
    result = analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber="auto",
    ).compute(distances)

    figure, axis = plt.subplots(figsize=(5.2, 4.4))
    for i, j, label in ((0, 0, r"$g_{11}$"), (0, 1, r"$g_{12}$"), (1, 1, r"$g_{22}$")):
        _ = axis.plot(distances.to("micrometer").magnitude, result.g[i, j], label=label)

    _ = axis.set(
        xlabel=r"separation $r$ ($\mu$m)",
        ylabel=r"partial correlation $g_{ij}(r)$",
        title="Binary Percus--Yevick mixture",
    )
    _ = axis.legend()
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(OUTPUT_DIRECTORY / "readme_percus_yevick.png", dpi=180, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    create_rsa_packing_figure()
    create_percus_yevick_figure()
