"""Compare PackLab with stored mctpy Percus--Yevick reference data.

This is an independent-development validation script, not a runtime feature.
It does not import mctpy: the reference values produced with mctpy 0.0.2.post1
are versioned as CSV data in the repository. Run the comparison with::

    python development/compare_mctpy.py

Both implementations use the density-scaled mixture convention

    S_ij(k) = delta_ij + sqrt(rho_i rho_j) h_tilde_ij(k),

for which the structure-factor matrix approaches the identity at large
wavenumber. The stored reference uses densities in micrometre**-3, diameters
in micrometres, and wavenumbers in micrometre**-1. See
``docs/manuscript/softwareX/data/README.md`` for its provenance and
regeneration command.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import PackLab
from PackLab import analytical, ureg


@dataclass(frozen=True)
class Mixture:
    """Physical inputs for one cross-implementation comparison."""

    name: str
    radii_um: tuple[float, ...]
    number_fractions: tuple[float, ...]
    volume_fraction: float


MIXTURES = (
    Mixture("monodisperse", (0.5,), (1.0,), 0.20),
    Mixture("binary", (0.35, 0.70), (0.65, 0.35), 0.25),
    Mixture("ternary", (0.30, 0.50, 0.80), (0.50, 0.30, 0.20), 0.30),
)

REFERENCE_DATA_PATH = (
    Path(__file__).parents[1]
    / "docs"
    / "manuscript"
    / "softwareX"
    / "data"
    / "mctpy_0_0_2_structure_factors.csv"
)


def partial_number_densities(mixture: Mixture) -> np.ndarray:
    """Return exact partial number densities in micrometre**-3."""
    radii = np.asarray(mixture.radii_um)
    fractions = np.asarray(mixture.number_fractions)
    mean_particle_volume = np.sum(fractions * (4.0 / 3.0) * np.pi * radii**3)
    return fractions * mixture.volume_fraction / mean_particle_volume


def packlab_structure_factor(
    mixture: Mixture, wavenumber_um_inverse: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return PackLab densities and partial structure factors."""
    radii = np.asarray(mixture.radii_um) * ureg.micrometer
    density_values = partial_number_densities(mixture)
    densities = density_values / ureg.micrometer**3
    solver = analytical.PercusYevickSolver(
        densities=densities,
        radii=radii,
        wavenumber=wavenumber_um_inverse / ureg.micrometer,
    )

    # compute() also evaluates g(r), although only its reciprocal-space result
    # is used here. The short distance interval comfortably satisfies the
    # solver's Fourier-grid sampling diagnostic.
    result = solver.compute(np.linspace(0.0, 2.0, 64) * ureg.micrometer)
    identity = np.eye(len(radii))[:, :, np.newaxis]
    structure_factor = identity + result.H
    return density_values, structure_factor


def load_mctpy_structure_factor(
    mixture: Mixture,
    reference_path: Path = REFERENCE_DATA_PATH,
) -> tuple[np.ndarray, np.ndarray]:
    """Load one mctpy reference and return it in PackLab axis order."""
    data = np.genfromtxt(
        reference_path,
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )
    selected = data[data["mixture"] == mixture.name]
    if selected.size == 0:
        raise RuntimeError(f"No reference rows found for mixture {mixture.name!r}.")

    wavenumber = np.unique(selected["wavenumber_um_inverse"])
    species_count = len(mixture.radii_um)
    structure_factor = np.empty((species_count, species_count, wavenumber.size))
    structure_factor.fill(np.nan)
    for row in range(species_count):
        for column in range(row, species_count):
            component = selected[(selected["row"] == row) & (selected["column"] == column)]
            if component.size != wavenumber.size:
                raise RuntimeError(
                    f"Reference component ({row}, {column}) for {mixture.name} "
                    f"contains {component.size} rows; expected {wavenumber.size}."
                )
            order = np.argsort(component["wavenumber_um_inverse"])
            values = component["S_ij"][order]
            structure_factor[row, column] = values
            structure_factor[column, row] = values

    if not np.all(np.isfinite(structure_factor)):
        raise RuntimeError(f"Reference data for {mixture.name} contain missing values.")
    return wavenumber, structure_factor


def compare_mixture(mixture: Mixture) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Evaluate PackLab on the stored reference grid and return error statistics."""
    wavenumber, mctpy = load_mctpy_structure_factor(mixture)
    _, packlab = packlab_structure_factor(mixture, wavenumber)

    if packlab.shape != mctpy.shape:
        raise RuntimeError(
            f"Unexpected array shapes for {mixture.name}: "
            f"PackLab {packlab.shape}, mctpy {mctpy.shape}."
        )

    difference = packlab - mctpy
    maximum_absolute = float(np.max(np.abs(difference)))
    root_mean_square = float(np.sqrt(np.mean(difference**2)))
    return wavenumber, packlab, mctpy, maximum_absolute, root_mean_square


def plot_comparison(
    mixture: Mixture,
    wavenumber_um_inverse: np.ndarray,
    packlab: np.ndarray,
    mctpy: np.ndarray,
    output_directory: Path,
) -> None:
    """Save one curve-and-residual comparison figure."""
    species_count = packlab.shape[0]
    figure, (curve_axis, residual_axis) = plt.subplots(
        2,
        1,
        figsize=(7.0, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": (3, 1)},
        constrained_layout=True,
    )

    for row in range(species_count):
        for column in range(row, species_count):
            label = rf"$S_{{{row + 1}{column + 1}}}$"
            (line,) = curve_axis.plot(
                wavenumber_um_inverse,
                packlab[row, column],
                label=f"PackLab {label}",
            )
            curve_axis.plot(
                wavenumber_um_inverse,
                mctpy[row, column],
                "--",
                color=line.get_color(),
                label=f"mctpy {label}",
            )
            residual_axis.plot(
                wavenumber_um_inverse,
                packlab[row, column] - mctpy[row, column],
                color=line.get_color(),
                label=label,
            )

    curve_axis.set_title(
        f"{mixture.name.capitalize()} PY structure-factor comparison "
        f"($\\phi={mixture.volume_fraction:.2f}$)"
    )
    curve_axis.set_ylabel(r"$S_{ij}(k)$")
    curve_axis.legend(ncol=2, fontsize="small")
    curve_axis.grid(alpha=0.2)
    residual_axis.axhline(0.0, color="black", linewidth=0.8)
    residual_axis.set_xlabel(r"Wavenumber $k$ ($\mathrm{\mu m}^{-1}$)")
    residual_axis.set_ylabel("Difference")
    residual_axis.grid(alpha=0.2)

    output_directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_directory / f"mctpy_{mixture.name}.png", dpi=180)
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("development/mctpy_comparison"),
        help="Directory for comparison figures.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-8,
        help="Maximum permitted absolute difference; use a negative value to disable.",
    )
    return parser.parse_args()


def main() -> None:
    """Run all comparisons, save figures, and enforce the error tolerance."""
    arguments = parse_arguments()
    failed = False
    print(f"PackLab {PackLab.__version__}; stored reference: mctpy 0.0.2.post1")
    print("mixture       species    max |delta S|       RMS delta S")
    print("----------------------------------------------------------")
    for mixture in MIXTURES:
        wavenumber, packlab, mctpy_values, maximum_absolute, root_mean_square = compare_mixture(
            mixture
        )
        print(
            f"{mixture.name:<13} {packlab.shape[0]:>3d}       "
            f"{maximum_absolute:>12.5e}    {root_mean_square:>12.5e}"
        )
        plot_comparison(
            mixture,
            wavenumber,
            packlab,
            mctpy_values,
            arguments.output_directory,
        )
        failed |= arguments.tolerance >= 0.0 and maximum_absolute > arguments.tolerance

    print(f"\nFigures written to {arguments.output_directory.resolve()}")
    if failed:
        raise SystemExit(
            f"At least one comparison exceeded the absolute tolerance {arguments.tolerance:.5e}."
        )


if __name__ == "__main__":
    main()
