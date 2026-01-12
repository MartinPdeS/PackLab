from dataclasses import dataclass
from typing import Callable, Optional, Tuple, Literal
from abc import ABC, abstractmethod

import numpy as np


BinSpacing = Literal["linear", "log"]

__all__ = [
    "Delta",
    "Uniform",
    "Gaussian",
    "LogNormal",
    "Custom",
    "Discrete"
]


def _as_quantity_1d(value):
    """Best effort conversion to a 1D Pint quantity array (no unit changes)."""
    value = np.asarray(value.magnitude) * value.units
    if value.ndim != 1:
        raise ValueError("Expected a one dimensional quantity array.")
    return value


def _get_unit_reference(*candidates):
    """Return the first Pint unit found among candidates."""
    for candidate in candidates:
        if candidate is None:
            continue
        return candidate.units
    raise ValueError("Could not infer units. Provide at least one Pint quantity parameter.")


def _make_bin_edges(
    radius_min,
    radius_max,
    number_of_bins: int,
    bin_spacing: BinSpacing,
):
    if number_of_bins < 1:
        raise ValueError("number_of_bins must be >= 1.")

    unit_reference = _get_unit_reference(radius_min, radius_max)
    radius_min = radius_min.to(unit_reference)
    radius_max = radius_max.to(unit_reference)

    if radius_max.magnitude <= radius_min.magnitude:
        raise ValueError("radius_max must be strictly larger than radius_min.")

    if bin_spacing == "linear":
        edges_magnitude = np.linspace(radius_min.magnitude, radius_max.magnitude, number_of_bins + 1)
        return edges_magnitude * unit_reference

    if bin_spacing == "log":
        if radius_min.magnitude <= 0:
            raise ValueError("radius_min must be > 0 for log spacing.")
        edges_magnitude = np.logspace(
            np.log10(radius_min.magnitude),
            np.log10(radius_max.magnitude),
            number_of_bins + 1,
        )
        return edges_magnitude * unit_reference

    raise ValueError("bin_spacing must be 'linear' or 'log'.")


def _edges_to_centers_and_widths(edges, bin_spacing: BinSpacing):
    unit_reference = edges.units
    edges_magnitude = np.asarray(edges.magnitude)

    if bin_spacing == "linear":
        centers_magnitude = 0.5 * (edges_magnitude[:-1] + edges_magnitude[1:])
        widths_magnitude = edges_magnitude[1:] - edges_magnitude[:-1]
        return centers_magnitude * unit_reference, widths_magnitude * unit_reference

    if bin_spacing == "log":
        centers_magnitude = np.sqrt(edges_magnitude[:-1] * edges_magnitude[1:])
        widths_magnitude = edges_magnitude[1:] - edges_magnitude[:-1]
        return centers_magnitude * unit_reference, widths_magnitude * unit_reference

    raise ValueError("bin_spacing must be 'linear' or 'log'.")


def _normalize_weights(weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    weights[~np.isfinite(weights)] = 0.0
    weights[weights < 0.0] = 0.0

    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError("Distribution produced zero total weight. Check parameters or range.")
    return weights / total


class RadiusDistribution(ABC):
    """Base class for parametric particle radius distributions.

    A distribution converts its parameters into a discrete mixture representation
    suitable for solvers that operate on bins.

    The discrete representation is:
    - particle_radii: Pint array of bin centers
    - weights: float array that sums to 1
    """

    @abstractmethod
    def to_bins(self) -> Tuple[object, np.ndarray]:
        """Discretize the distribution into bin centers and number fractions.

        Returns
        -------
        particle_radii : object
            Pint quantity array of bin centers, shape (number_of_bins,).
        weights : np.ndarray
            Array of weights, shape (number_of_bins,).
        """
        raise NotImplementedError


@dataclass(slots=True)
class Delta(RadiusDistribution):
    """Delta distribution (monodisperse).

    Parameters
    ----------
    radius : object
        The single radius value.
    """

    radius: object

    def to_bins(self) -> Tuple[object, np.ndarray]:
        particle_radii = np.asarray([self.radius.to(self.radius.units).magnitude]) * self.radius.units
        weights = np.asarray([1.0], dtype=float)
        return particle_radii, weights


@dataclass(slots=True)
class Uniform(RadiusDistribution):
    """Uniform distribution over a finite interval.

    The probability mass is approximated by integrating a constant pdf over each bin,
    which reduces to weights proportional to bin widths.

    Parameters
    ----------
    radius_min : object
        Lower bound of the support.
    radius_max : object
        Upper bound of the support.
    number_of_bins : int
        Number of bins.
    bin_spacing : BinSpacing
        "linear" or "log".
    """

    radius_min: object
    radius_max: object
    number_of_bins: int
    bin_spacing: BinSpacing = "linear"

    def to_bins(self) -> Tuple[object, np.ndarray]:
        edges = _make_bin_edges(self.radius_min, self.radius_max, self.number_of_bins, self.bin_spacing)
        particle_radii, bin_widths = _edges_to_centers_and_widths(edges, self.bin_spacing)

        weights = bin_widths.to(particle_radii.units).magnitude
        weights = _normalize_weights(weights)
        return particle_radii, weights


@dataclass(slots=True)
class Gaussian(RadiusDistribution):
    """
    Gaussian distribution truncated to a finite interval.

    The number fractions are obtained by approximating probability mass in each bin as: weight_k = pdf(center_k) * width_k followed by normalization.

    Parameters
    ----------
    mean_radius : object
        Mean radius.
    standard_deviation : object
        Standard deviation.
    radius_min : object
        Lower bound of truncation.
    radius_max : object
        Upper bound of truncation.
    number_of_bins : int
        Number of bins.
    bin_spacing : BinSpacing
        "linear" or "log".

    """

    mean_radius: object
    standard_deviation: object
    radius_min: object
    radius_max: object
    number_of_bins: int
    bin_spacing: BinSpacing = "linear"

    def to_bins(self) -> Tuple[object, np.ndarray]:
        unit_reference = _get_unit_reference(self.mean_radius, self.standard_deviation, self.radius_min, self.radius_max)

        mean_radius = self.mean_radius.to(unit_reference)
        standard_deviation = self.standard_deviation.to(unit_reference)

        if standard_deviation.magnitude <= 0:
            raise ValueError("standard_deviation must be > 0.")

        edges = _make_bin_edges(self.radius_min, self.radius_max, self.number_of_bins, self.bin_spacing)
        particle_radii, bin_widths = _edges_to_centers_and_widths(edges, self.bin_spacing)

        x = particle_radii.to(unit_reference).magnitude
        mu = mean_radius.magnitude
        sigma = standard_deviation.magnitude

        pdf = np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))
        weights = pdf * bin_widths.to(unit_reference).magnitude

        weights = _normalize_weights(weights)
        return particle_radii, weights


@dataclass(slots=True)
class LogNormal(RadiusDistribution):
    """Log normal distribution truncated to a finite interval.

    Parameterization uses median radius and geometric standard deviation.

    The pdf is evaluated in a dimensionless log space using magnitudes expressed
    in a consistent unit (the unit of median_radius).

    Parameters
    ----------
    median_radius : object
        Median radius.
    geometric_standard_deviation : object
        Geometric standard deviation (must be > 1).
    radius_min : object
        Lower bound of truncation.
    radius_max : object
        Upper bound of truncation.
    number_of_bins : int
        Number of bins.
    bin_spacing : BinSpacing
        "linear" or "log".
    """

    median_radius: object
    geometric_standard_deviation: object
    radius_min: object
    radius_max: object
    number_of_bins: int
    bin_spacing: BinSpacing = "log"

    def to_bins(self) -> Tuple[object, np.ndarray]:
        unit_reference = _get_unit_reference(self.median_radius, self.radius_min, self.radius_max)

        median_radius = self.median_radius.to(unit_reference)
        geometric_standard_deviation = self.geometric_standard_deviation

        if geometric_standard_deviation <= 0:
            raise ValueError("geometric_standard_deviation must be > 0.")

        sigma = np.log(float(geometric_standard_deviation))
        if sigma <= 0:
            raise ValueError("geometric_standard_deviation must be > 1 to represent a non degenerate distribution.")

        edges = _make_bin_edges(self.radius_min, self.radius_max, self.number_of_bins, self.bin_spacing)
        particle_radii, bin_widths = _edges_to_centers_and_widths(edges, self.bin_spacing)

        x = particle_radii.to(unit_reference).magnitude
        m = median_radius.magnitude

        if np.any(x <= 0) or m <= 0:
            raise ValueError("Log normal requires strictly positive radii and median_radius.")

        mu = np.log(m)

        pdf = np.exp(-((np.log(x) - mu) ** 2) / (2.0 * sigma**2)) / (x * sigma * np.sqrt(2.0 * np.pi))
        weights = pdf * bin_widths.to(unit_reference).magnitude

        weights = _normalize_weights(weights)
        return particle_radii, weights


@dataclass(slots=True)
class Custom(RadiusDistribution):
    """User defined pdf over a finite interval.

    The callable is evaluated at bin centers in magnitude space, expressed in the
    unit of radius_min. Probability mass is approximated as: weight_k = pdf(center_k) * width_k then normalized.

    Parameters
    ----------
    pdf : Callable[[np.ndarray], np.ndarray]
        User defined probability density function. Accepts and returns arrays of shape (number_of_bins,).
    radius_min : object
        Lower bound of the support.
    radius_max : object
        Upper bound of the support.
    number_of_bins : int
        Number of bins.
    bin_spacing : BinSpacing
        "linear" or "log".
    """

    pdf: Callable[[np.ndarray], np.ndarray]
    radius_min: object
    radius_max: object
    number_of_bins: int
    bin_spacing: BinSpacing = "linear"

    def to_bins(self) -> Tuple[object, np.ndarray]:
        unit_reference = _get_unit_reference(self.radius_min, self.radius_max)

        edges = _make_bin_edges(self.radius_min, self.radius_max, self.number_of_bins, self.bin_spacing)
        particle_radii, bin_widths = _edges_to_centers_and_widths(edges, self.bin_spacing)

        x = particle_radii.to(unit_reference).magnitude
        pdf_values = np.asarray(self.pdf(x), dtype=float)

        if pdf_values.shape != x.shape:
            raise ValueError("Custom pdf must return an array with the same shape as its input.")

        weights = pdf_values * bin_widths.to(unit_reference).magnitude
        weights = _normalize_weights(weights)
        return particle_radii, weights


@dataclass(slots=True)
class Discrete(RadiusDistribution):
    """Explicit discrete radius distribution.

    This distribution is already discretized: the user provides bin centers and
    associated probability masses. The class validates and normalizes the masses.

    Parameters
    ----------
    particle_radii : object
        Pint quantity array of bin centers.
    weights : np.ndarray
        Array of number fractions (not necessarily normalized).
    """

    particle_radii: object
    weights: np.ndarray

    def to_bins(self) -> Tuple[object, np.ndarray]:
        particle_radii = _as_quantity_1d(self.particle_radii)

        weights = np.asarray(self.weights, dtype=float)
        if weights.ndim != 1:
            raise ValueError("weights must be a one dimensional array.")

        if particle_radii.magnitude.shape[0] != weights.shape[0]:
            raise ValueError(
                "particle_radii and weights must have the same length."
            )

        if weights.size < 1:
            raise ValueError("Discrete distribution must contain at least one radius value.")

        weights = _normalize_weights(weights)
        return particle_radii, weights