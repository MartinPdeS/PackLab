"""Coverage for validation and plotting helpers implemented in Python."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from PackLab.monte_carlo.diagnostics import run_metropolis_diagnostics
from PackLab.scattering import ScatteringData, ScatteringDataset
from PackLab.scattering.plottings import (
    _coerce_phase_function_axes,
    plot_phase_function_2d_projection,
    plot_phase_function_3d,
)
from PackLab.units import ureg


def _dataset(phi=None):
    phi = np.array([0.0, np.pi / 4, np.pi / 2]) * ureg.radian if phi is None else phi
    dataset = ScatteringDataset()
    dataset.k = 2.0 / ureg.meter
    dataset.phi = phi
    dataset.append(
        ScatteringData(
            S1=np.array([1.0, 2.0, 3.0]) * ureg.dimensionless,
            S2=np.array([3.0, 2.0, 1.0]) * ureg.dimensionless,
            k=dataset.k,
            Csca=2.0 * ureg.meter**2,
            phi=phi,
        )
    )
    return dataset


def test_scattering_dataset_validates_processing_inputs():
    with pytest.raises(ValueError, match="at least one"):
        ScatteringDataset().process()

    empty_phi = ScatteringDataset()
    empty_phi.append(ScatteringData([], [], 1 / ureg.meter, 1 * ureg.meter**2, []))
    with pytest.raises(ValueError, match="phi must be set"):
        empty_phi.process()

    bad_item = ScatteringDataset()
    bad_item.phi = np.array([0.0, 1.0]) * ureg.radian
    bad_item.append(object())
    with pytest.raises(TypeError, match="ScatteringData"):
        bad_item.process()

    bad_grid = _dataset(np.array([0.0, 0.5, 0.4]) * ureg.radian)
    with pytest.raises(ValueError, match="strictly increasing"):
        bad_grid.process()


def test_scattering_dataset_rejects_invalid_mixture_values():
    dataset = _dataset()
    wavenumber = np.array([0.0, 1.0, 2.0, 3.0]) / ureg.meter
    H = np.zeros((1, 1, 4))

    with pytest.raises(ValueError, match="non-negative"):
        dataset.get_mu_independant(np.array([-1.0]) / ureg.meter**3)
    with pytest.raises(ValueError, match="at least 2"):
        dataset.get_phase_function(
            np.array([1.0]) / ureg.meter**3,
            H,
            wavenumber,
            theta_points=1,
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        dataset.get_phase_function(
            np.array([1.0]) / ureg.meter**3,
            H,
            np.array([0.0, 2.0, 1.0]) / ureg.meter,
        )


def test_phase_function_plot_helpers_support_both_layouts_and_modes():
    phi = np.array([0.0, np.pi / 2])
    theta = np.linspace(0.0, 2.0 * np.pi, 3)
    values = np.arange(6, dtype=float).reshape(2, 3)

    np.testing.assert_array_equal(_coerce_phase_function_axes(values.T, phi, theta), values)
    with pytest.raises(ValueError, match="incompatible shape"):
        _coerce_phase_function_axes(np.zeros((4, 4)), phi, theta)

    figure_3d = plot_phase_function_3d(phi, theta, values, mode="surface", normalize=False)
    assert figure_3d.axes[0].get_zlabel() == "P"
    figure_spherical = plot_phase_function_3d(
        phi, theta, values.T, mode="spherical", use_magnitude=False
    )
    assert figure_spherical.axes[0].get_title().startswith("Phase function mapped")

    figure_average = plot_phase_function_2d_projection(
        phi, theta, values.T, projection="azimuth_average", normalize=True
    )
    assert len(figure_average.axes[0].lines) == 1
    figure_heatmap = plot_phase_function_2d_projection(
        phi, theta, values, projection="heatmap", use_magnitude=False
    )
    assert len(figure_heatmap.axes[0].images) == 1

    with pytest.raises(ValueError, match='projection must be'):
        plot_phase_function_2d_projection(phi, theta, values, projection="invalid")
    with pytest.raises(ValueError, match='mode must be'):
        plot_phase_function_3d(phi, theta, values, mode="invalid")


class _FakeStatistics:
    accepted_moves = 0
    rejected_moves = 0


class _FakeDomain:
    length_x = 4.0 * ureg.meter
    length_y = 4.0 * ureg.meter
    length_z = 4.0 * ureg.meter
    use_periodic_boundaries = True


class _FakeConfiguration:
    positions = np.array([[0.5, 0.5, 0.5]]) * ureg.meter


class _FakeSimulator:
    def __init__(self):
        self.domain = _FakeDomain()
        self.sphere_configuration = _FakeConfiguration()
        self.statistics = _FakeStatistics()
        self.cumulative_displacements = np.zeros((1, 3)) * ureg.meter

    def run_sweeps(self, sweeps):
        self.statistics.accepted_moves += sweeps
        self.sphere_configuration.positions += 0.1 * sweeps * ureg.meter
        self.cumulative_displacements += 0.1 * sweeps * ureg.meter
        return object()


def test_metropolis_diagnostics_support_custom_observables():
    simulator = _FakeSimulator()
    report = run_metropolis_diagnostics(
        simulator,
        4,
        sample_interval=1,
        block_size=2,
        observable=lambda configuration: float(
            configuration.positions.to("meter").magnitude[0, 0]
        ),
        observable_name="x_position_m",
    )

    assert report.observable_name == "x_position_m"
    assert report.observable.tolist() == pytest.approx([0.6, 0.7, 0.8, 0.9])
    assert report.acceptance_rates.tolist() == [1.0] * 4
    assert len(report.block_means) == 2


def test_metropolis_diagnostics_correlated_stationary_observable():
    """An AR(1) process has tau=(1+rho)/(2*(1-rho)) and known mean error."""
    rng = np.random.default_rng(2026)
    rho = 0.8
    values = np.empty(8192)
    values[0] = rng.normal()
    for index in range(1, len(values)):
        values[index] = rho * values[index - 1] + np.sqrt(1 - rho**2) * rng.normal()
    samples = iter(values)
    report = run_metropolis_diagnostics(
        _FakeSimulator(), len(values), sample_interval=1, block_size=64,
        observable=lambda _: next(samples), observable_name="AR(1)",
    )
    expected_tau = (1 + rho) / (2 * (1 - rho))
    expected_error = np.sqrt(2 * expected_tau / len(values))
    assert report.autocorrelation_time == pytest.approx(expected_tau, rel=0.25)
    assert report.effective_sample_size == pytest.approx(len(values) / (2 * expected_tau), rel=0.3)
    assert report.standard_error == pytest.approx(expected_error, rel=0.3)
    assert report.confidence_interval[0] < 0 < report.confidence_interval[1]


def test_metropolis_diagnostics_excludes_short_final_interval_from_statistics():
    values = iter([1.0, 2.0, 3.0, 100.0])
    report = run_metropolis_diagnostics(
        _FakeSimulator(), 7, sample_interval=2, block_size=1,
        observable=lambda _: next(values),
    )
    assert report.sample_sweeps.tolist() == [2, 4, 6, 7]
    assert report.observable.tolist() == [1.0, 2.0, 3.0, 100.0]
    np.testing.assert_array_equal(report.block_means, [1.0, 2.0, 3.0])


def test_metropolis_diagnostics_constant_observable_has_no_uncertainty_estimate():
    report = run_metropolis_diagnostics(
        _FakeSimulator(), 4, sample_interval=1, observable=lambda _: 1.0,
    )
    np.testing.assert_array_equal(report.block_means, [1.0] * 4)
    assert report.autocorrelation_time is None
    assert report.effective_sample_size is None
    assert report.standard_error is None
    assert report.confidence_interval is None
    assert report.confidence_level is None


def test_metropolis_diagnostics_validates_known_errors_before_advancing():
    simulator = _FakeSimulator()
    with pytest.raises(ValueError, match="block_size requires an observable"):
        run_metropolis_diagnostics(simulator, 10, block_size=1)
    with pytest.raises(ValueError, match="two complete blocks"):
        run_metropolis_diagnostics(simulator, 10, sample_interval=2, block_size=3,
                                   observable=lambda _: 1.0)
    assert simulator.statistics.accepted_moves == 0


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"number_of_sweeps": 0}, "number_of_sweeps"),
        ({"number_of_sweeps": 2, "sample_interval": 0}, "sample_interval"),
        ({"number_of_sweeps": 2, "burn_in_sweeps": 2}, "burn_in_sweeps"),
        ({"number_of_sweeps": 2, "confidence_level": 1.0}, "confidence_level"),
    ],
)
def test_metropolis_diagnostics_rejects_invalid_controls(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        run_metropolis_diagnostics(_FakeSimulator(), **kwargs)
