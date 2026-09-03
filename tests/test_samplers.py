# test_samplers.py

import math
import numpy as np
import pytest

from PackLab.units import ureg

from PackLab import samplers


def _as_meters_1d(quantity) -> np.ndarray:
    values = quantity.to("meter").magnitude
    values = np.asarray(values, dtype=float)
    assert values.ndim == 1
    return values


def _assert_weights_normalized(weights, tol=1e-12):
    weights = np.asarray(weights, dtype=float)
    assert weights.ndim == 1
    assert np.all(np.isfinite(weights))
    assert np.all(weights >= 0.0)
    total = float(np.sum(weights))
    assert math.isfinite(total)
    assert total > 0.0
    assert abs(total - 1.0) <= tol


def test_constant_to_bins_returns_single_center_and_unit_weight():
    radius = 2.0 * ureg.micrometer
    sampler = samplers.ConstantRadiusSampler(radius=2.0 * ureg.micrometer)

    particle_radii, weights = sampler.to_bins()

    radii_m = _as_meters_1d(particle_radii)
    weights = np.asarray(weights, dtype=float)

    assert radii_m.shape == (1,)
    assert weights.shape == (1,)
    assert radii_m[0] == pytest.approx((2.0 * ureg.micrometer).to("meter").magnitude)
    assert weights[0] == pytest.approx(1.0)


def test_plot_histogram_returns_axes_with_radius_classes():
    import matplotlib.pyplot as plt

    sampler = samplers.UniformRadiusSampler(
        minimum_radius=100.0 * ureg.nanometer,
        maximum_radius=200.0 * ureg.nanometer,
        bins=4,
    )
    figure, axis = plt.subplots()

    returned_axis = sampler.plot_histogram(ax=axis, color="tab:blue")

    assert returned_axis is axis
    assert len(axis.patches) == 4
    assert axis.get_xlabel() == "particle radius (nanometer)"
    assert axis.get_ylabel() == "number fraction"
    plt.close(figure)


def test_constant_sampler_runs_rsa_with_a_valid_single_class_index():
    """A default constant sampler must not emit the invalid class index -1."""
    from PackLab import monte_carlo

    domain = monte_carlo.PackingDomain(
        3.0 * ureg.micrometer,
        3.0 * ureg.micrometer,
        3.0 * ureg.micrometer,
        use_periodic_boundaries=True,
    )
    options = monte_carlo.RSAOptions()
    options.random_seed = 7
    options.maximum_attempts = 1_000
    options.maximum_spheres = 12

    result = monte_carlo.RSASimulator(
        domain,
        samplers.ConstantRadiusSampler(150.0 * ureg.nanometer),
        options,
    ).run()

    class_index = np.asarray(result.sphere_configuration.classes_index)
    assert result.sphere_configuration.number_of_classes == 1
    assert np.all(class_index == 0)


def test_uniform_to_bins_requires_bins_positive_if_that_is_your_policy():
    sampler = samplers.UniformRadiusSampler(
        minimum_radius=1.0 * ureg.micrometer,
        maximum_radius=3.0 * ureg.micrometer,
        bins=0,
    )

    # If your C++ implementation requires bins>0 for to_bins, keep this test.
    # If you changed policy to allow bins=0, remove this and add an expected behavior test.
    with pytest.raises(RuntimeError):
        sampler.to_bins()


def test_uniform_to_bins_linear_edges_gives_uniform_weights():
    bins = 4
    rmin = 1.0 * ureg.micrometer
    rmax = 3.0 * ureg.micrometer

    sampler = samplers.UniformRadiusSampler(minimum_radius=rmin, maximum_radius=rmax, bins=bins)
    particle_radii, weights = sampler.to_bins()

    radii_m = _as_meters_1d(particle_radii)
    weights = np.asarray(weights, dtype=float)

    assert radii_m.shape == (bins,)
    assert weights.shape == (bins,)

    _assert_weights_normalized(weights)

    # Uniform distribution on linear bins -> weights should be equal (or extremely close)
    assert np.allclose(weights, np.ones(bins) / bins, rtol=0.0, atol=1e-12)

    # Centers should be midpoints of equally spaced edges
    rmin_m = rmin.to("meter").magnitude
    rmax_m = rmax.to("meter").magnitude
    edges = np.linspace(rmin_m, rmax_m, bins + 1)
    centers_expected = 0.5 * (edges[:-1] + edges[1:])
    assert np.allclose(radii_m, centers_expected, rtol=0.0, atol=1e-15)


def test_lognormal_to_bins_positive_weights_and_normalized():
    bins = 40


    sampler = samplers.LogNormalRadiusSampler(
        median_radius=1.0 * ureg.micrometer,
        geometric_standard_deviation=1.4,
        maximum_radius_clip=100.0 * ureg.micrometer,
        bins=bins,
    )

    particle_radii, weights = sampler.to_bins()
    radii_m = _as_meters_1d(particle_radii)
    weights = np.asarray(weights, dtype=float)

    assert radii_m.shape == (bins,)
    assert weights.shape == (bins,)
    _assert_weights_normalized(weights)
    assert np.all(radii_m > 0.0)



def test_discrete_to_bins_returns_normalized_weights_and_original_radii():
    particle_radii = np.array([1.0, 2.0]) * ureg.micrometer
    weights_in = np.array([1.0, 1.0], dtype=float)

    sampler = samplers.DiscreteRadiusSampler(radii=particle_radii, weights=weights_in)

    particle_radii_out, weights_out = sampler.to_bins()

    radii_m = _as_meters_1d(particle_radii_out)
    weights_out = np.asarray(weights_out, dtype=float)

    assert radii_m.shape == (2,)
    assert np.allclose(radii_m, particle_radii.to("meter").magnitude, rtol=0.0, atol=0.0)

    _assert_weights_normalized(weights_out)
    assert np.allclose(weights_out, np.array([0.5, 0.5]), rtol=0.0, atol=1e-12)


def test_discrete_rejects_mismatched_lengths():
    particle_radii = np.array([1.0, 2.0]) * ureg.micrometer
    weights_in = np.array([1.0], dtype=float)

    with pytest.raises(Exception):
        samplers.DiscreteRadiusSampler(radii=particle_radii, weights=weights_in)


def test_discrete_rejects_negative_weights():
    particle_radii = np.array([1.0, 2.0]) * ureg.micrometer
    weights_in = np.array([1.0, -1.0], dtype=float)

    with pytest.raises(Exception):
        samplers.DiscreteRadiusSampler(radii=particle_radii, weights=weights_in)


def test_uniform_rejects_invalid_range():
    with pytest.raises(Exception):
        samplers.UniformRadiusSampler(
            minimum_radius=2.0 * ureg.micrometer,
            maximum_radius=1.0 * ureg.micrometer,
            bins=10,
        )

if __name__ == "__main__":
    pytest.main(["-W error", __file__])
