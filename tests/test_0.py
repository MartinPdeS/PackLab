"""
PackLab test.
"""

import pytest
import numpy as np
import pytest
from unittest.mock import patch

from PackLab.monte_carlo import (
    Domain,
    Options,
    Simulator,
    UniformRadiusSampler,
)


# ----------------------------------------------------------
# Utility: minimum image distance to test overlaps
# ----------------------------------------------------------
def _min_distance_periodic(positions, box_lengths):
    """Compute the minimum center-to-center distance with
    minimum-image convention. Used only for testing overlap."""
    if len(positions) < 2:
        return np.inf

    positions = np.asarray(positions)
    Lx, Ly, Lz = box_lengths

    dx = positions[:, None, 0] - positions[None, :, 0]
    dy = positions[:, None, 1] - positions[None, :, 1]
    dz = positions[:, None, 2] - positions[None, :, 2]

    dx -= Lx * np.round(dx / Lx)
    dy -= Ly * np.round(dy / Ly)
    dz -= Lz * np.round(dz / Lz)

    distances = np.sqrt(dx**2 + dy**2 + dz**2)

    i_upper = np.triu_indices_from(distances, k=1)
    return float(np.min(distances[i_upper]))


# ==========================================================
#  Basic sanity tests
# ==========================================================
def test_basic_rsa_run():
    domain = Domain(
        length_x=5.0,
        length_y=5.0,
        length_z=5.0,
        use_periodic_boundaries=True,
    )

    radius_sampler = UniformRadiusSampler(
        minimum_radius=0.2,
        maximum_radius=0.2,
        bins=10
    )

    options = Options()
    options.random_seed = 123
    options.maximum_attempts = 200_000
    options.maximum_consecutive_rejections = 50_000
    options.target_packing_fraction = 0.10

    simulator = Simulator(
        domain=domain,
        radius_sampler=radius_sampler,
        options=options
    )
    result = simulator.run()

    positions = result.positions
    radii = result.radii
    stats = result.statistics

    # Shapes
    assert positions.ndim == 2
    assert positions.shape[1] == 3
    assert radii.ndim == 1
    assert radii.shape[0] == positions.shape[0]

    # Statistics exist and are finite
    assert stats.sphere_count == positions.shape[0]
    assert np.isfinite(stats.packing_fraction_simulator or stats.packing_fraction_geometry)


# ==========================================================
#  Test NO OVERLAP in periodic mode
# ==========================================================
def test_no_overlap_periodic():
    domain = Domain(6.0, 6.0, 6.0, use_periodic_boundaries=True)

    radius_sampler = UniformRadiusSampler(0.15, 0.15, bins=10)

    options = Options()
    options.random_seed = 99
    options.maximum_attempts = 300_000
    options.maximum_consecutive_rejections = 30_000
    options.target_packing_fraction = 0.15

    simulator = Simulator(domain=domain, radius_sampler=radius_sampler, options=options)
    result = simulator.run()

    positions = result.positions
    radii = result.radii

    if len(positions) < 2:
        pytest.skip("Too few spheres placed to test overlap.")

    min_dist = _min_distance_periodic(
        positions, (domain.length_x, domain.length_y, domain.length_z)
    )
    cutoff = 2.0 * np.min(radii)

    assert min_dist + 1e-12 >= cutoff


# ==========================================================
#  Test that packing fraction matches geometry
# ==========================================================
def test_packing_fraction_consistency():
    domain = Domain(6.0, 6.0, 6.0, use_periodic_boundaries=False)

    radius_sampler = UniformRadiusSampler(0.2, 0.2, bins=10)

    options = Options()
    options.random_seed = 42
    options.maximum_attempts = 200_000
    options.maximum_consecutive_rejections = 50_000
    options.target_packing_fraction = 0.12

    simulator = Simulator(domain=domain, radius_sampler=radius_sampler, options=options)
    result = simulator.run()

    positions = result.positions
    radii = result.radii
    stats = result.statistics

    if len(radii) == 0:
        pytest.skip("No spheres placed.")

    volume = domain.length_x * domain.length_y * domain.length_z
    sphere_volume = (4.0 / 3.0) * np.pi * np.sum(radii**3)
    pf_geom = sphere_volume / volume

    pf_sim = stats.packing_fraction_simulator or stats.packing_fraction_geometry

    assert abs(pf_geom - pf_sim) < 5e-3


# ==========================================================
#  Test stopping based on maximum spheres
# ==========================================================
def test_stop_by_maximum_spheres():
    domain = Domain(10.0, 10.0, 10.0, use_periodic_boundaries=True)

    radius_sampler = UniformRadiusSampler(0.1, 0.1, bins=10)

    options = Options()
    options.random_seed = 11
    options.maximum_attempts = 1_000_000
    options.maximum_consecutive_rejections = 1_000_000
    options.maximum_spheres = 25

    simulator = Simulator(domain=domain, radius_sampler=radius_sampler, options=options)
    result = simulator.run()

    assert result.positions.shape[0] <= 25


# ==========================================================
#  Test that plotting functions run without errors
# ==========================================================
@patch('matplotlib.pyplot.show')
def test_plot_slice_runs(patch):
    domain = Domain(4.0, 4.0, 4.0, use_periodic_boundaries=True)
    radius_sampler = UniformRadiusSampler(0.15, 0.15, bins=10)

    options = Options()
    options.random_seed = 2
    options.maximum_attempts = 80_000

    result = Simulator(domain, radius_sampler, options).run()

    fig = result.plot_slice_2d(show=False)


@patch('matplotlib.pyplot.show')
def test_plot_pair_correlation_runs(patch):
    domain = Domain(4.0, 4.0, 4.0, use_periodic_boundaries=True)
    radius_sampler = UniformRadiusSampler(0.15, 0.15, bins=10)

    options = Options()
    options.random_seed = 3
    options.maximum_attempts = 80_000

    result = Simulator(domain, radius_sampler, options).run()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
