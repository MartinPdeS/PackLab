"""Utilities for constructing reciprocal-space grids."""

from typing import Any

import numpy as np

from PackLab.units import ureg


def make_wavenumber_grid(
    radial_resolution: Any,
    maximum_distance: Any,
    samples_per_oscillation: int = 12,
) -> Any:
    """Create a reciprocal-space grid for the radial inverse transform.

    Parameters
    ----------
    radial_resolution : pint.Quantity
        Smallest real-space feature to resolve. It determines the maximum
        wavenumber through ``pi / radial_resolution``.
    maximum_distance : pint.Quantity
        Largest real-space distance at which the correlation function will be
        evaluated. It determines the wavenumber spacing.
    samples_per_oscillation : int, default=12
        Samples used for one period of the sinc kernel at ``maximum_distance``.
        Values below 8 are not recommended.

    Returns
    -------
    pint.Quantity
        A uniformly spaced, zero-inclusive wavenumber grid with units of
        inverse length.
    """
    if samples_per_oscillation < 2:
        raise ValueError("samples_per_oscillation must be at least 2.")

    resolution_m = radial_resolution.to("meter")
    maximum_distance_m = maximum_distance.to("meter")
    if resolution_m.magnitude <= 0:
        raise ValueError("radial_resolution must be positive.")
    if maximum_distance_m.magnitude <= 0:
        raise ValueError("maximum_distance must be positive.")

    maximum_wavenumber = np.pi / resolution_m
    wavenumber_step = 2 * np.pi / (samples_per_oscillation * maximum_distance_m)
    number_of_points = int(
        np.ceil((maximum_wavenumber / wavenumber_step).to_base_units().magnitude)
    ) + 1

    return np.linspace(
        0.0,
        maximum_wavenumber.to("1 / meter").magnitude,
        number_of_points,
    ) / ureg.meter
