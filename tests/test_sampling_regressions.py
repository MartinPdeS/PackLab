"""Accuracy and boundary regressions through PackLab's public workflows."""

from statistics import NormalDist

import numpy as np
import pytest

from PackLab import monte_carlo, samplers
from PackLab.units import ureg


def _configuration(positions, radii, classes=None):
    if classes is None:
        classes = np.zeros(len(radii), dtype=int)
    return monte_carlo.PackingConfiguration.from_arrays(
        np.asarray(positions) * ureg.meter,
        np.asarray(radii) * ureg.meter,
        np.asarray(classes, dtype=int),
    )


def _metropolis(domain, configuration, displacement=0.1, seed=37):
    options = monte_carlo.MetropolisOptions()
    options.random_seed = seed
    options.maximum_displacement = displacement * ureg.meter
    return monte_carlo.MetropolisSimulator(domain, configuration, options)


def test_capped_pair_correlations_converge_for_both_particle_orders():
    """A class-sorted prefix must not bias the sampled mixture correlations."""
    domain = monte_carlo.PackingDomain(*(8 * ureg.meter for _ in range(3)), True)
    positions = np.column_stack((np.linspace(0.2, 7.8, 40), np.ones(40), np.ones(40)))
    classes = np.repeat([0, 1], 20)
    forward = _configuration(positions, np.full(40, 0.001), classes)
    reverse = _configuration(positions[::-1], np.full(40, 0.001), classes[::-1])
    _, exact = forward.compute_partial_pair_correlation_function(domain, n_bins=8, maximum_pairs=0)
    _, reversed_exact = reverse.compute_partial_pair_correlation_function(
        domain, n_bins=8, maximum_pairs=0
    )
    np.testing.assert_array_equal(exact, reversed_exact)
    for configuration in (forward, reverse):
        values = [configuration.compute_partial_pair_correlation_function(
            domain, n_bins=8, maximum_pairs=120, random_seed=seed
        )[1] for seed in range(256)]
        averaged = np.mean(values, axis=0)
        assert np.linalg.norm(averaged - exact) / np.linalg.norm(exact) < 0.04
        first = configuration.compute_partial_pair_correlation_function(
            domain, n_bins=8, maximum_pairs=120, random_seed=7
        )[1]
        second = configuration.compute_partial_pair_correlation_function(
            domain, n_bins=8, maximum_pairs=120, random_seed=7
        )[1]
        np.testing.assert_array_equal(first, second)


def test_distance_sampling_is_uniform_without_replacement_and_respects_cap():
    domain = monte_carlo.PackingDomain(*(20 * ureg.meter for _ in range(3)), True)
    positions = np.column_stack((np.arange(8) * 0.1, np.ones(8), np.ones(8)))
    configuration = _configuration(positions, np.full(8, 0.01), np.arange(8))
    result = _metropolis(domain, configuration).run_sweeps(0)
    counts = np.zeros((8, 8), dtype=int)
    for seed in range(400):
        distances, first, second = result.binding.compute_partial_pair_distances(5, seed)
        assert len(distances) == 5
        assert len(set(zip(first, second))) == 5
        assert np.all(first < second)
        counts[first, second] += 1
    observed = counts[np.triu_indices(8, 1)]
    expected = 400 * 5 / 28
    assert np.all(np.abs(observed - expected) < 0.5 * expected)
    all_distances = result.binding.compute_partial_pair_distances(0, 7)[0]
    assert len(all_distances) == 28
    np.testing.assert_array_equal(
        all_distances, result.binding.compute_partial_pair_distances(100, 99)[0]
    )


def test_distance_cap_counts_examined_pairs_before_radial_filtering():
    domain = monte_carlo.PackingDomain(*(4 * ureg.meter for _ in range(3)), False)
    configuration = _configuration([[0.1, 1, 1], [3.9, 1, 1], [0.2, 1, 1]], [0.01] * 3)
    result = _metropolis(domain, configuration).run_sweeps(0)
    for seed in range(20):
        distances = result.binding.compute_partial_pair_distances(1, seed)[0]
        assert len(distances) <= 1
        assert np.all(distances < 2)
    # Non-periodic separations must not be shortened across opposite walls.
    np.testing.assert_allclose(result.binding.compute_partial_pair_distances(0)[0], [0.1])


@pytest.mark.parametrize("enforce", [False, True])
@pytest.mark.parametrize("radius, expected", [(1.01, 0), (1.0, 1), (0.99, 1)])
def test_periodic_rsa_rejects_self_image_overlap_but_allows_contact(enforce, radius, expected):
    domain = monte_carlo.PackingDomain(2 * ureg.meter, 3 * ureg.meter, 4 * ureg.meter, True)
    options = monte_carlo.RSAOptions()
    options.random_seed = 5
    options.maximum_attempts = 10
    options.maximum_consecutive_rejections = 10
    options.maximum_spheres = 1
    options.enforce_radii_distribution = enforce
    result = monte_carlo.RSASimulator(
        domain, samplers.ConstantRadiusSampler(radius * ureg.meter), options
    ).run()
    assert result.sphere_configuration.count == expected


def test_metropolis_and_archives_reject_periodic_self_image_overlap(tmp_path):
    domain = monte_carlo.PackingDomain(2 * ureg.meter, 3 * ureg.meter, 4 * ureg.meter, True)
    oversized = _configuration([[1, 1, 1]], [1.01])
    with pytest.raises(ValueError, match="own periodic image"):
        _metropolis(domain, oversized)
    path = tmp_path / "packing.npz"
    with pytest.raises(ValueError, match="own periodic image"):
        monte_carlo.save_packing(oversized, path, domain=domain)
    touching = _configuration([[1, 1, 1]], [1.0])
    _metropolis(domain, touching).run_sweeps(0)
    monte_carlo.save_packing(touching, path, domain=domain)
    assert monte_carlo.load_packing(path).sphere_configuration.count == 1
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    arrays["radii_m"] = np.array([1.01])
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match="own periodic image"):
        monte_carlo.load_packing(path)


def test_metropolis_unwrapped_msd_includes_multiple_crossings_and_call_origin():
    domain = monte_carlo.PackingDomain(*(1 * ureg.meter for _ in range(3)), True)
    configuration = _configuration([[0.5, 0.5, 0.5]], [0.01])
    simulator = _metropolis(domain, configuration, displacement=3.0)
    simulator.run_sweeps(5)
    origin = simulator.cumulative_displacements.to("meter").magnitude.copy()
    report = monte_carlo.run_metropolis_diagnostics(
        simulator, 8, sample_interval=2, block_size=1, observable=lambda _: 1.0
    )
    total = simulator.cumulative_displacements.to("meter").magnitude
    assert report.mean_squared_displacements.magnitude[-1] == pytest.approx(
        np.sum((total - origin) ** 2)
    )
    assert np.max(report.mean_squared_displacements.magnitude) > 0.75  # wrapped-MSD bound
    np.testing.assert_allclose(
        simulator.sphere_configuration.positions.magnitude,
        (configuration.positions.magnitude + total) % 1.0,
        atol=1e-13,
    )
    # Property returns independent storage, and reset clears the accumulated travel.
    simulator.cumulative_displacements.magnitude[:] = 999
    np.testing.assert_array_equal(simulator.cumulative_displacements.magnitude, total)
    simulator.reset()
    np.testing.assert_array_equal(simulator.cumulative_displacements.magnitude, np.zeros((1, 3)))


def test_unwrapped_displacements_exclude_rejected_moves():
    domain = monte_carlo.PackingDomain(*(1 * ureg.meter for _ in range(3)), False)
    configuration = _configuration([[0.5, 0.5, 0.5]], [0.1])
    simulator = _metropolis(domain, configuration, displacement=0.6)
    simulator.run_sweeps(100)
    assert simulator.statistics.rejected_moves > 0
    np.testing.assert_allclose(
        simulator.cumulative_displacements.magnitude,
        simulator.sphere_configuration.positions.magnitude - configuration.positions.magnitude,
        atol=1e-14,
    )


def test_dilute_periodic_msd_matches_proposal_variance():
    """Independent walkers provide an analytical MSD accuracy reference."""
    domain = monte_carlo.PackingDomain(*(1000 * ureg.meter for _ in range(3)), True)
    rng = np.random.default_rng(17)
    configuration = _configuration(rng.uniform(0, 1000, (1000, 3)), np.full(1000, 1e-6))
    simulator = _metropolis(domain, configuration, displacement=2.0)
    report = monte_carlo.run_metropolis_diagnostics(
        simulator, 8, sample_interval=2, block_size=1, observable=lambda _: 1.0
    )
    # Each component has variance d^2/3 per accepted sweep; three components give t*d^2.
    expected = report.sample_sweeps * 2.0**2
    np.testing.assert_allclose(report.mean_squared_displacements.magnitude, expected, rtol=0.1)
    assert simulator.statistics.rejected_moves == 0


@pytest.mark.parametrize("bins", [0, 8])
@pytest.mark.parametrize("mean", [0.0, 0.5])
def test_normal_samples_match_positive_conditioned_upper_clipped_probabilities(
    bins, mean, tmp_path
):
    sigma, cutoff = 1.0, 2.0
    normal = NormalDist(mean, sigma)
    edges = np.linspace(0, cutoff, 9)
    expected = np.diff([normal.cdf(edge) for edge in edges]) / (1 - normal.cdf(0))
    expected[-1] += (1 - normal.cdf(cutoff)) / (1 - normal.cdf(0))
    sampler = samplers.NormalRadiusSampler(
        mean * ureg.meter, sigma * ureg.meter, cutoff * ureg.meter, bins=bins
    )
    domain = monte_carlo.PackingDomain(*(1000 * ureg.meter for _ in range(3)), True)
    options = monte_carlo.RSAOptions()
    options.random_seed = 19
    options.maximum_spheres = 6000
    options.maximum_attempts = 10000
    result = monte_carlo.RSASimulator(domain, sampler, options).run()
    radii = result.radii.to("meter").magnitude
    assert len(radii) == 6000
    assert np.all((radii > 0) & (radii <= cutoff))
    observed = np.histogram(radii, bins=edges)[0] / radii.size
    np.testing.assert_allclose(observed, expected, atol=0.02, rtol=0)
    if bins:
        _, weights = sampler.to_bins()
        np.testing.assert_allclose(weights, expected, atol=1e-14)
    monte_carlo.save_packing(result, tmp_path / "normal.npz")


@pytest.mark.parametrize("mean,sigma", [(0, 0), (-1, 1), (-1, 0)])
def test_normal_rejects_invalid_location_or_zero_radius_degenerate_distribution(mean, sigma):
    with pytest.raises(ValueError, match="mean"):
        samplers.NormalRadiusSampler(mean * ureg.meter, sigma * ureg.meter, 2 * ureg.meter)


def test_normal_zero_variance_is_positive_and_upper_clipped():
    sampler = samplers.NormalRadiusSampler(3 * ureg.meter, 0 * ureg.meter, 2 * ureg.meter, 4)
    radii, weights = sampler.to_bins()
    np.testing.assert_array_equal(weights, [0, 0, 0, 1])
    assert radii.magnitude[-1] == pytest.approx(1.75)


@pytest.mark.parametrize("second_x", [0.9, 3.9])
def test_metropolis_checks_overlap_with_translated_periodic_images(second_x):
    domain = monte_carlo.PackingDomain(*(1 * ureg.meter for _ in range(3)), True)
    configuration = _configuration([[0.1, 0.5, 0.5], [second_x, 0.5, 0.5]], [0.15, 0.15])
    with pytest.raises(ValueError, match="overlapping spheres"):
        _metropolis(domain, configuration)


def test_normal_zero_variance_without_bins_runs_as_positive_constant():
    domain = monte_carlo.PackingDomain(*(10 * ureg.meter for _ in range(3)), True)
    options = monte_carlo.RSAOptions()
    options.maximum_spheres = 1
    options.random_seed = 7
    result = monte_carlo.RSASimulator(
        domain,
        samplers.NormalRadiusSampler(3 * ureg.meter, 0 * ureg.meter, 2 * ureg.meter),
        options,
    ).run()
    np.testing.assert_array_equal(result.radii.magnitude, [2.0])
    np.testing.assert_array_equal(result.sphere_configuration.classes_index, [0])
