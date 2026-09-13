"""Empirical finite-configuration structure calculations."""

from __future__ import annotations

from dataclasses import dataclass
import warnings
from typing import Any, Literal

import numpy as np

from PackLab.analytical.grid import make_wavenumber_grid
from PackLab.units import ureg


@dataclass(frozen=True, slots=True)
class EmpiricalStructure:
    """
    Pair structure estimated from one explicit finite packing configuration.

    This is empirical finite-configuration data, not a Percus--Yevick result.
    ``H`` follows the existing scattering convention
    ``sqrt(n_i n_j) Fourier[g_ij - 1]`` and is dimensionless; it should not be
    interpreted as a PY closure solution.
    """

    source: str
    densities: Any
    radii: Any
    classes: np.ndarray
    distances: Any
    g: np.ndarray
    h: np.ndarray
    wavenumber: Any
    H: np.ndarray
    S: np.ndarray

    def scattering_inputs(self) -> tuple[Any, np.ndarray, Any]:
        """
        Return ``(densities, H, wavenumber)`` for ``ScatteringDataset`` methods.

        Returns
        -------
        tuple
            Species densities, dimensionless correlation tensor, and reciprocal
            wavenumber grid in PackLab's existing mixture convention.
        """
        return self.densities, self.H, self.wavenumber


def _validate_wavenumber(wavenumber: Any) -> np.ndarray:
    try:
        values = np.asarray(wavenumber.to("1 / meter").magnitude, dtype=float)
    except AttributeError as error:
        raise TypeError(
            "wavenumber must be a unit-bearing quantity convertible to 1 / meter."
        ) from error
    if values.ndim != 1 or values.size < 2:
        raise ValueError("wavenumber must be one-dimensional with at least two points.")
    if (
        not np.all(np.isfinite(values))
        or np.any(values < 0.0)
        or np.any(np.diff(values) <= 0.0)
    ):
        raise ValueError("wavenumber must be finite, non-negative, and strictly increasing.")
    return values


def _packing_data(packing: Any, domain: Any | None) -> tuple[Any, Any, np.ndarray, Any]:
    if hasattr(packing, "sphere_configuration") and hasattr(packing, "domain"):
        return (
            packing.sphere_configuration,
            packing.domain,
            np.asarray(packing.sphere_configuration.classes_index),
            packing,
        )
    if domain is None:
        raise ValueError("domain is required for a bare PackingConfiguration.")
    if not hasattr(packing, "classes_index"):
        raise TypeError(
            "packing must be a PackingResult, loaded configuration, or PackingConfiguration."
        )
    return packing, domain, np.asarray(packing.classes_index), packing


def empirical_structure(
    packing: Any,
    *,
    domain: Any | None = None,
    n_bins: int = 64,
    maximum_pairs: int = 1_000_000,
    wavenumber: Any | Literal["auto"] = "auto",
    samples_per_oscillation: int = 12,
) -> EmpiricalStructure:
    """
    Estimate empirical ``g_ij(r)`` and reciprocal correlations from a packing.

    Parameters
    ----------
    packing : PackingResult, LoadedPackingConfiguration, or PackingConfiguration
        Explicit hard-sphere configuration. This workflow is empirical and is
        not a Percus--Yevick calculation.
    domain : PackingDomain, optional
        Required only with a bare native ``PackingConfiguration``.
    n_bins : int, default=64
        Number of radial bins used for partial pair correlations.
    maximum_pairs : int, default=1_000_000
        Maximum native pair count; zero requests all available pairs.
    wavenumber : pint.Quantity or "auto", default="auto"
        Reciprocal grid. Automatic selection uses the radial bin width and
        finite correlation range.
    samples_per_oscillation : int, default=12
        Reciprocal samples per sinc-kernel period for automatic selection.

    Returns
    -------
    EmpiricalStructure
        Finite-configuration densities, class radii, ``g``, ``h=g-1``,
        dimensionless scattering-compatible ``H``, and ``S=I+H``.

    Warns
    -----
    RuntimeWarning
        Pair counting and truncating ``h(r)`` at half the shortest box length
        produce finite-size and reciprocal-resolution artifacts.
    """
    if (
        isinstance(n_bins, bool)
        or not isinstance(n_bins, (int, np.integer))
        or n_bins < 2
    ):
        raise ValueError("n_bins must be an integer of at least 2.")
    if (
        isinstance(maximum_pairs, bool)
        or not isinstance(maximum_pairs, (int, np.integer))
        or maximum_pairs < 0
    ):
        raise ValueError("maximum_pairs must be a non-negative integer.")
    configuration, packing_domain, classes, owner = _packing_data(packing, domain)
    if configuration.count < 2:
        raise ValueError("at least two particles are required for empirical structure estimation.")
    if classes.shape != (configuration.count,) or np.any(classes < 0):
        raise ValueError(
            "configuration class labels must be non-negative and match particle count."
        )

    if hasattr(owner, "sphere_configuration"):
        distances, g = owner.compute_partial_pair_correlation_function(
            n_bins=n_bins, maximum_pairs=maximum_pairs
        )
    else:
        distances, g = configuration.compute_partial_pair_correlation_function(
            packing_domain, n_bins=n_bins, maximum_pairs=maximum_pairs
        )
        distances = np.asarray(distances) * ureg.meter
    distances_m = np.asarray(distances.to("meter").magnitude, dtype=float)
    g = np.asarray(g, dtype=float)
    if not np.all(np.isfinite(g)):
        raise ValueError("empirical pair correlation contains non-finite values.")
    active_classes = np.unique(classes)
    g = g[active_classes][:, active_classes]
    radial_resolution = float(np.min(np.diff(distances_m)))
    if radial_resolution <= 0.0:
        raise ValueError("pair-correlation distances must be strictly increasing.")
    if isinstance(wavenumber, str) and wavenumber == "auto":
        if samples_per_oscillation < 2:
            raise ValueError("samples_per_oscillation must be at least 2.")
        wavenumber_values = np.asarray(
            make_wavenumber_grid(
                radial_resolution * ureg.meter,
                distances_m[-1] * ureg.meter,
                samples_per_oscillation,
            )
            .to("1 / meter")
            .magnitude
        )
    elif isinstance(wavenumber, str):
        raise ValueError("wavenumber must be a unit-bearing quantity or 'auto'.")
    else:
        wavenumber_values = _validate_wavenumber(wavenumber)
        samples = 2.0 * np.pi / (np.max(np.diff(wavenumber_values)) * distances_m[-1])
        if samples < 8.0:
            warnings.warn(
                "The supplied wavenumber grid is likely too coarse for the "
                "empirical correlation range; finite-size and transform "
                "artifacts may be substantial.",
                RuntimeWarning,
                stacklevel=2,
            )

    radii_values = np.asarray(configuration.radii.to("meter").magnitude)
    number_of_classes = active_classes.size
    counts = np.asarray([(classes == class_index).sum() for class_index in active_classes])
    volume = float(packing_domain.volume)
    densities = counts / volume
    class_radii = np.array(
        [
            radii_values[classes == class_index].mean() for class_index in active_classes
        ]
    )
    h = g - 1.0
    kernel = np.sinc(np.outer(wavenumber_values, distances_m) / np.pi)
    transform = 4.0 * np.pi * np.trapezoid(
        h[:, :, None, :] * (distances_m**2 * kernel)[None, None, :, :],
        x=distances_m,
        axis=-1,
    )
    H = np.sqrt(densities[:, None] * densities[None, :])[:, :, None] * transform
    S = H + np.eye(number_of_classes)[:, :, None]
    warnings.warn(
        "Empirical finite-configuration correlations are truncated at the sampled box range; "
        "finite-size, binning, and reciprocal-transform artifacts are expected.",
        RuntimeWarning,
        stacklevel=2,
    )
    return EmpiricalStructure(
        source="empirical finite-configuration",
        densities=densities / ureg.meter**3,
        radii=class_radii * ureg.meter,
        classes=active_classes,
        distances=distances_m * ureg.meter,
        g=g,
        h=h,
        wavenumber=wavenumber_values / ureg.meter,
        H=H,
        S=S,
    )
