"""Tests for analytical reciprocal-space grid handling."""

import numpy as np
import pytest

from PackLab import analytical
from PackLab.units import ureg


def _solver(wavenumber):
    return analytical.PercusYevickSolver(
        densities=np.array([0.1]) / ureg.meter**3,
        radii=np.array([0.1]) * ureg.meter,
        wavenumber=wavenumber,
    )


def test_make_wavenumber_grid_resolves_requested_distance():
    wavenumber = analytical.make_wavenumber_grid(
        radial_resolution=0.05 * ureg.meter,
        maximum_distance=2.0 * ureg.meter,
    )

    assert wavenumber[0] == 0 / ureg.meter
    assert len(wavenumber) > 2
    _solver(wavenumber).compute(np.linspace(0.2, 2.0, 20) * ureg.meter)


def test_solver_warns_for_a_coarse_wavenumber_grid():
    wavenumber = np.linspace(0.0, 10.0, 3) / ureg.meter

    with pytest.warns(RuntimeWarning, match="wavenumber grid is likely too coarse"):
        _solver(wavenumber).compute(np.linspace(0.2, 10.0, 20) * ureg.meter)


def test_solver_rejects_non_increasing_wavenumber_grid():
    with pytest.raises(ValueError, match="strictly increasing"):
        _solver(np.array([0.0, 2.0, 1.0]) / ureg.meter)


def test_solver_does_not_accept_the_removed_p_keyword():
    with pytest.raises(TypeError):
        analytical.PercusYevickSolver(
            densities=np.array([0.1]) / ureg.meter**3,
            radii=np.array([0.1]) * ureg.meter,
            p=np.linspace(0.0, 10.0, 20) / ureg.meter,
        )


def test_solver_automatically_selects_a_wavenumber_grid():
    solver = analytical.PercusYevickSolver(
        densities=np.array([0.1]) / ureg.meter**3,
        radii=np.array([0.1]) * ureg.meter,
        wavenumber="auto",
    )
    result = solver.compute(np.linspace(0.2, 2.0, 20) * ureg.meter)

    assert "wavenumber='auto'" in repr(solver)
    assert result.wavenumber[0] == 0 / ureg.meter
    assert len(result.wavenumber) > 2


def test_automatic_wavenumber_grid_accepts_resolution_controls():
    solver = analytical.PercusYevickSolver(
        densities=np.array([0.1]) / ureg.meter**3,
        radii=np.array([0.1]) * ureg.meter,
        wavenumber="auto",
        radial_resolution=0.025 * ureg.meter,
        samples_per_oscillation=16,
    )
    result = solver.compute(np.linspace(0.2, 2.0, 20) * ureg.meter)

    assert np.isclose(
        result.wavenumber.max().to("1 / meter").magnitude,
        np.pi / 0.025,
    )
