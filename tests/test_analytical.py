"""Tests for analytical reciprocal-space grid handling."""

from concurrent.futures import ThreadPoolExecutor

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


def test_solver_enforces_hard_core_exclusion():
    """Fourier ringing must not produce non-zero in-core correlations."""
    radius = 0.1 * ureg.meter
    solver = analytical.PercusYevickSolver(
        densities=np.array([0.1]) / ureg.meter**3,
        radii=np.array([radius.magnitude]) * radius.units,
        wavenumber="auto",
    )
    distances = np.linspace(0.0, 0.199, 20) * ureg.meter

    result = solver.compute(distances)

    np.testing.assert_array_equal(result.g[0, 0], 0.0)


def test_automatic_wavenumber_grid_converges_against_refined_grid():
    """The default automatic grid should agree with a twice-refined integral."""
    radius = 0.1 * ureg.meter
    volume_fraction = 0.20
    density = volume_fraction / ((4.0 / 3.0) * np.pi * radius**3)
    radii = np.array([radius.magnitude]) * radius.units
    densities = np.array([density.magnitude]) * density.units
    distances = np.linspace(0.2, 0.8, 120) * ureg.meter

    automatic = analytical.PercusYevickSolver(
        densities=densities,
        radii=radii,
        wavenumber="auto",
    ).compute(distances)

    refined_wavenumber = analytical.make_wavenumber_grid(
        radial_resolution=radius / 20,
        maximum_distance=distances.max(),
        samples_per_oscillation=24,
    )
    refined = analytical.PercusYevickSolver(
        densities=densities,
        radii=radii,
        wavenumber=refined_wavenumber,
    ).compute(distances)

    assert len(refined.wavenumber) > len(automatic.wavenumber)
    np.testing.assert_allclose(automatic.g, refined.g, rtol=5e-3, atol=5e-4)


@pytest.mark.parametrize("volume_fraction,size_ratio", [(0.05, 1.5), (0.20, 2.0), (0.35, 3.0)])
def test_binary_automatic_grid_converges_across_regimes(volume_fraction, size_ratio):
    """Resolve all partial correlations away from their hard-core contacts."""
    radii = np.array([0.1, 0.1 * size_ratio]) * ureg.micrometer
    domain = analytical.PercusYevickDomain(
        size=100 * ureg.micrometer,
        radii=radii,
        volume_fraction=volume_fraction,
        number_fractions=[0.5, 0.5],
    )
    distances = np.linspace(0, 1.5, 160) * ureg.micrometer
    kwargs = dict(densities=domain.particle_densities_per_radius, radii=domain.radii)
    automatic = analytical.PercusYevickSolver(**kwargs, wavenumber="auto").compute(distances)
    grid = analytical.make_wavenumber_grid(
        radial_resolution=radii.min() / 20,
        maximum_distance=distances.max(),
        samples_per_oscillation=24,
    )
    refined = analytical.PercusYevickSolver(**kwargs, wavenumber=grid).compute(distances)

    for i, j in ((0, 0), (0, 1), (1, 1)):
        contact = (radii[i] + radii[j]).to("micrometer").magnitude
        mask = distances.magnitude > contact + 0.04
        np.testing.assert_allclose(
            automatic.g[i, j, mask], refined.g[i, j, mask], atol=5e-4, rtol=5e-3
        )


@pytest.mark.parametrize("start", [0.0, 0.15])
def test_inverse_transform_matches_direct_nonuniform_quadrature(start):
    """Independently integrate H for all pairs, including zero density and r=0."""
    wavenumber = start + np.linspace(0.0, 1.0, 401) ** 1.3 * 80.0
    distances = np.array([0.0, 0.05, 0.2, 0.4, 0.8])
    densities = np.array([0.7, 1.2, 0.0])
    radii = np.array([0.1, 0.2, 0.3])
    result = analytical.PercusYevickSolver(
        densities=densities / ureg.meter**3,
        radii=radii * ureg.meter,
        wavenumber=wavenumber / ureg.meter,
    ).compute(distances * ureg.meter)

    for i in range(3):
        for j in range(3):
            if densities[i] == 0.0 or densities[j] == 0.0:
                np.testing.assert_array_equal(result.h[i, j], 0.0)
                continue
            for index, distance in enumerate(distances):
                integrand = (
                    result.H[i, j] * wavenumber**2
                    * np.sinc(wavenumber * distance / np.pi)
                )
                integral = np.sum(
                    0.5 * np.diff(wavenumber) * (integrand[:-1] + integrand[1:])
                )
                expected = integral / (2 * np.pi**2 * np.sqrt(densities[i] * densities[j]))
                assert result.h[i, j, index] == pytest.approx(expected, rel=1e-11, abs=1e-12)
    np.testing.assert_allclose(result.Cpy, result.Cpy.swapaxes(0, 1), atol=1e-14)
    np.testing.assert_allclose(result.H, result.H.swapaxes(0, 1), atol=1e-14)


def test_solver_supports_concurrent_calls():
    solver = _solver("auto")
    distances = np.linspace(0.0, 1.0, 80) * ureg.meter
    expected = solver.compute(distances).g
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(solver.compute, [distances] * 4))
    for result in results:
        np.testing.assert_allclose(result.g, expected, rtol=1e-13, atol=1e-13)
