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

    def run_sweeps(self, sweeps):
        self.statistics.accepted_moves += sweeps
        self.sphere_configuration.positions += 0.1 * sweeps * ureg.meter
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
