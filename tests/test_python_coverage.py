"""Additional branch coverage for the Python-facing workflow helpers."""

import json

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from PackLab.analytical.grid import make_wavenumber_grid
from PackLab.monte_carlo.persistence import (
    _extract_packing,
    _quantity_values,
    _validate_arrays,
    load_packing,
)
from PackLab.monte_carlo.results import PackingResult
from PackLab.monte_carlo.structure import _validate_wavenumber, empirical_structure
from PackLab.scattering import ScatteringData, ScatteringDataset
from PackLab.scattering.model import compute_scattering_amplitudes
from PackLab.units import ureg


class _Domain:
    length_x = 4.0 * ureg.meter
    length_y = 5.0 * ureg.meter
    length_z = 6.0 * ureg.meter
    use_periodic_boundaries = True


class _Configuration:
    positions = np.array([[0.5, 0.5, 0.5], [2.0, 2.0, 2.0]]) * ureg.meter
    radii = np.array([0.1, 0.2]) * ureg.meter
    classes_index = np.array([0, 1], dtype=np.int64)
    count = 2

    def compute_partial_pair_correlation_function(self, *args, **kwargs):
        return np.linspace(0.1, 0.4, 4), np.ones((2, 2, 4))


class _Binding:
    sphere_configuration = _Configuration()
    domain = _Domain()
    statistics = {"accepted": 2}
    partial_volume_fractions = [0.2, 0.3]
    partial_volumes = [1.0, 2.0]
    pair_correlation_centers = [0.1, 0.2]
    pair_correlation_values = [1.0, 0.8]

    def compute_partial_pair_correlation_function(self, **kwargs):
        return np.linspace(0.1, 0.4, 4), np.ones((2, 2, 4))

    def compute_pair_correlation_function(self, **kwargs):
        self.pair_correlation_called = kwargs
        return "computed"


def test_packing_result_properties_delegation_and_plots():
    result = PackingResult(_Binding(), source="rsa", run_metadata={"seed": 1})

    assert result.source == "rsa"
    assert result.run_metadata == {"seed": 1}
    assert result.positions is result.positions
    assert result.radii is result.radii
    assert result.statistics == {"accepted": 2}
    assert result.sphere_configuration is result.sphere_configuration
    assert result.domain is result.domain
    assert result.partial_volume_fractions.tolist() == [0.2, 0.3]
    assert result.partial_volumes.tolist() == [1.0, 2.0]
    centers, values = result.compute_partial_pair_correlation_function(n_bins=4)
    assert centers.units == ureg.meter
    assert values.shape == (2, 2, 4)
    assert result.compute_pair_correlation_function(n_bins=2) == "computed"
    assert result.pair_correlation_centers.shape == (2,)
    assert result.pair_correlation_values.shape == (2,)

    figures = [
        result.plot_centers_3d(show=False),
        result.plot_radius_distribution(show=False, density=False),
        result.plot_slice_2d(slice_axis="x", show=False),
        result.plot_slice_2d(slice_axis="y", show=False),
        result.plot_pair_correlation(show=False),
    ]
    assert all(len(figure.axes) for figure in figures)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"slice_center_fraction": -0.1}, "slice_center_fraction"),
        ({"slice_thickness_fraction": 1.1}, "slice_thickness_fraction"),
        ({"slice_axis": "q"}, "q"),
    ],
)
def test_packing_result_slice_validates_controls(kwargs, message):
    with pytest.raises((ValueError, KeyError), match=message):
        PackingResult(_Binding()).plot_slice_2d(show=False, **kwargs)


def test_result_slice_mask_handles_periodic_and_nonperiodic_domains():
    result = PackingResult(_Binding())
    coord = np.array([0.1, 3.9])
    periodic = result._compute_slice_mask(coord, 0.0, 0.4, 4.0)
    assert periodic.tolist() == [True, True]
    result.domain.use_periodic_boundaries = False
    nonperiodic = result._compute_slice_mask(coord, 0.0, 0.4, 4.0)
    assert nonperiodic.tolist() == [True, False]
    result.domain.use_periodic_boundaries = True


def test_analytical_and_structure_grid_validation_errors():
    grid = make_wavenumber_grid(0.1 * ureg.meter, 1.0 * ureg.meter)
    assert grid[0] == 0 / ureg.meter
    assert grid[-1].to("1 / meter").magnitude == pytest.approx(10 * np.pi)

    for kwargs, message in [
        ({"samples_per_oscillation": 1}, "at least 2"),
        ({"radial_resolution": -1 * ureg.meter}, "radial_resolution"),
        ({"maximum_distance": -1 * ureg.meter}, "maximum_distance"),
    ]:
        values = {
            "radial_resolution": 0.1 * ureg.meter,
            "maximum_distance": 1.0 * ureg.meter,
            **kwargs,
        }
        with pytest.raises(ValueError, match=message):
            make_wavenumber_grid(**values)

    with pytest.raises(TypeError, match="unit-bearing"):
        _validate_wavenumber([0.0, 1.0])
    with pytest.raises(ValueError, match="at least two"):
        _validate_wavenumber(np.array([0.0]) / ureg.meter)
    with pytest.raises(ValueError, match="non-negative"):
        _validate_wavenumber(np.array([-1.0, 1.0]) / ureg.meter)


class _StructureDomain:
    volume = 125.0


class _StructureConfiguration:
    count = 3
    classes_index = np.array([0, 2, 2], dtype=np.int64)
    radii = np.array([0.1, 0.2, 0.2]) * ureg.meter

    def compute_partial_pair_correlation_function(self, **kwargs):
        distances = np.linspace(0.1, 1.0, 8) * ureg.meter
        values = np.ones((3, 3, 8))
        return distances, values


class _StructureOwner:
    sphere_configuration = _StructureConfiguration()
    domain = _StructureDomain()

    def compute_partial_pair_correlation_function(self, **kwargs):
        return self.sphere_configuration.compute_partial_pair_correlation_function(**kwargs)


def test_empirical_structure_supports_automatic_and_manual_grids():
    owner = _StructureOwner()
    with pytest.warns(RuntimeWarning, match="finite-configuration"):
        automatic = empirical_structure(owner, n_bins=8, samples_per_oscillation=8)
    assert automatic.source == "empirical finite-configuration"
    assert automatic.classes.tolist() == [0, 2]
    assert automatic.g.shape == (2, 2, 8)
    assert automatic.H.shape == (2, 2, automatic.wavenumber.size)
    assert automatic.S.shape == automatic.H.shape
    assert automatic.scattering_inputs()[0].units == ureg.meter**-3

    manual_grid = np.linspace(0.0, 10.0, 8) / ureg.meter
    with pytest.warns(RuntimeWarning, match="likely too coarse") as warning:
        with pytest.warns(RuntimeWarning, match="finite-configuration"):
            manual = empirical_structure(owner, n_bins=8, wavenumber=manual_grid)
    assert len(warning) == 1
    assert manual.wavenumber.shape == manual_grid.shape

    with pytest.raises(ValueError, match="wavenumber must be"):
        empirical_structure(owner, n_bins=8, wavenumber="bad")
    with pytest.raises(ValueError, match="samples_per_oscillation"):
        empirical_structure(owner, n_bins=8, samples_per_oscillation=1)


def test_persistence_array_and_quantity_validation():
    valid_positions = np.array([[0.5, 0.5, 0.5], [2.0, 2.0, 2.0]])
    valid_radii = np.array([0.1, 0.1])
    valid_classes = np.array([0, 1], dtype=np.int64)
    valid_lengths = np.array([4.0, 4.0, 4.0])
    _validate_arrays(valid_positions, valid_radii, valid_classes, valid_lengths, True)
    assert _quantity_values(2 * ureg.meter, "meter", "length") == 2.0

    invalid_cases = [
        (valid_positions[:, :2], valid_radii, valid_classes, valid_lengths, True, "shape"),
        (valid_positions, np.array([0.1]), valid_classes, valid_lengths, True, "radii"),
        (valid_positions, valid_radii, np.array([-1, 0]), valid_lengths, True, "non-negative"),
        (valid_positions, valid_radii, valid_classes, np.array([4.0, 4.0]), True, "lengths"),
        (
            np.array([[0.0, 0.5, 0.5], [2.0, 2.0, 2.0]]),
            valid_radii,
            valid_classes,
            valid_lengths,
            False,
            "outside",
        ),
    ]
    for positions, radii, classes, lengths, periodic, message in invalid_cases:
        with pytest.raises(ValueError, match=message):
            _validate_arrays(positions, radii, classes, lengths, periodic)

    with pytest.raises(TypeError, match="unit-bearing"):
        _quantity_values(2.0, "meter", "length")


def test_persistence_extracts_supported_shapes_and_rejects_unknown():
    result = PackingResult(_Binding(), source="rsa", run_metadata={"seed": 2})
    configuration, domain, source, metadata = _extract_packing(result, None)
    assert configuration is result.sphere_configuration
    assert domain is result.domain
    assert (source, metadata) == ("rsa", {"seed": 2})

    with pytest.raises(TypeError, match="packing must be"):
        _extract_packing(object(), None)


def test_load_packing_rejects_bad_archive_metadata(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(path, unexpected=np.array([1]))
    with pytest.raises(ValueError, match="missing"):
        load_packing(path)

    metadata = {
        "schema": "wrong",
        "version": 1,
        "source": "unknown",
        "run_metadata": {},
        "units": {"positions_m": "meter", "radii_m": "meter", "box_lengths_m": "meter"},
    }
    np.savez(
        path,
        positions_m=np.array([[0.5, 0.5, 0.5]]),
        radii_m=np.array([0.1]),
        classes_index=np.array([0], dtype=np.int64),
        box_lengths_m=np.array([2.0, 2.0, 2.0]),
        periodic=np.asarray(True),
        metadata_json=np.asarray(json.dumps(metadata)),
    )
    with pytest.raises(ValueError, match="unsupported packing schema"):
        load_packing(path)


def test_scattering_dataset_interpolation_and_mixture_calculations():
    phi = np.array([0.0, np.pi / 4, np.pi / 2]) * ureg.radian
    dataset = ScatteringDataset(
        [
            ScatteringData(
                np.array([1.0, 2.0, 3.0]) * ureg.dimensionless,
                np.array([3.0, 2.0, 1.0]) * ureg.dimensionless,
                2 / ureg.meter,
                2 * ureg.meter**2,
                phi,
            )
        ]
    )
    dataset.k = 2 / ureg.meter
    dataset.phi = phi
    densities = np.array([2.0]) / ureg.meter**3
    wavenumber = np.linspace(0.0, 4.0, 5) / ureg.meter
    H = np.ones((1, 1, 5))

    assert dataset.get_mu_independant(densities).to("1/meter").magnitude == pytest.approx(4.0)
    dependent = dataset.get_mu_dependant(densities, H, wavenumber, theta_points=5)
    assert np.isfinite(dependent.magnitude)
    total = dataset.get_mu(densities, H, wavenumber, theta_points=5)
    assert np.isfinite(total.magnitude)
    assert dataset.interpolate_last_axis_linear(
        np.array([[0.0, 2.0]]), np.array([0.0, 2.0]), np.array([1.0])
    ).item() == 1.0


def test_scattering_model_debug_and_plot_branches(monkeypatch, capsys):
    class Gaussian:
        def __init__(self, **kwargs):
            self.wavenumber_vacuum = 2 / ureg.meter

    class Sphere:
        Csca = 1 * ureg.meter**2

        def __init__(self, **kwargs):
            pass

    class PolarizationState:
        def __init__(self, **kwargs):
            pass

    class Setup:
        def __init__(self, **kwargs):
            pass

        def get_s1s2(self, *, angles):
            return angles * 0 + 1, angles * 0 + 2

        def get(self, name):
            return 1 * ureg.meter**2

    from PackLab.scattering import model

    monkeypatch.setattr(model, "Gaussian", Gaussian)
    monkeypatch.setattr(model, "Sphere", Sphere)
    monkeypatch.setattr(model, "PolarizationState", PolarizationState)
    monkeypatch.setattr(model, "Setup", Setup)
    result = compute_scattering_amplitudes(
        500 * ureg.nanometer,
        [100 * ureg.nanometer],
        1.5,
        1.0,
        np.array([0.0, 0.5]) * ureg.radian,
        plot=True,
        debug_mode=True,
    )
    assert len(result) == 1
    assert "Diameter:" in capsys.readouterr().out
