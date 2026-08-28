"""Tests for the public scattering data model and helper."""

import numpy as np
import pytest

from PackLab.scattering import ScatteringData, compute_scattering_amplitudes
from PackLab.scattering import model
from PackLab.units import ureg


class _FakeGaussian:
    def __init__(self, **kwargs):
        self.wavenumber_vacuum = 2.0 / ureg.meter


class _FakeSphere:
    def __init__(self, **kwargs):
        pass


class _FakePolarizationState:
    def __init__(self, **kwargs):
        pass


class _FakeSetup:
    def __init__(self, **kwargs):
        pass

    def get_s1s2(self, *, angles):
        values = np.asarray(angles.to("radian").magnitude)
        return values * ureg.dimensionless, (values + 1) * ureg.dimensionless

    def get(self, name):
        assert name == "Csca"
        return 3.0 * ureg.meter**2


def test_get_s1s2_rejects_empty_diameters():
    with pytest.raises(ValueError, match="at least one"):
        compute_scattering_amplitudes(
            wavelength=500 * ureg.nanometer,
            diameters=[],
            material=1.5,
            medium=1.0,
            phi=np.linspace(0, np.pi, 3) * ureg.radian,
        )


def test_domain_apis_match_documented_namespaces():
    from PackLab import analytical, monte_carlo, samplers

    assert monte_carlo.PackingDomain.__module__ == "PackLab.monte_carlo.domain"
    assert monte_carlo.RSAOptions.__module__ == "PackLab.monte_carlo.simulator"
    assert monte_carlo.RSASimulator.__module__ == "PackLab.monte_carlo.simulator"
    assert samplers.UniformRadiusSampler.__module__ == "PackLab.samplers"
    assert analytical.PercusYevickDomain.__module__ == "PackLab.analytical.domain"


def test_public_binding_classes_have_docstrings_and_representations():
    from PackLab import monte_carlo, samplers

    domain = monte_carlo.PackingDomain(
        1.0 * ureg.meter,
        2.0 * ureg.meter,
        3.0 * ureg.meter,
        use_periodic_boundaries=True,
    )
    options = monte_carlo.RSAOptions()
    sampler = samplers.UniformRadiusSampler(0.1 * ureg.meter, 0.2 * ureg.meter)

    assert "Parameters" in monte_carlo.PackingDomain.__doc__
    assert "Parameters" in monte_carlo.RSASimulator.__doc__
    assert "Parameters" in samplers.UniformRadiusSampler.__doc__
    assert "PackingDomain" in repr(domain)
    assert "RSAOptions" in repr(options)
    assert "UniformRadiusSampler" in repr(sampler)


def test_get_s1s2_returns_typed_data(monkeypatch):
    monkeypatch.setattr(model, "Gaussian", _FakeGaussian)
    monkeypatch.setattr(model, "Sphere", _FakeSphere)
    monkeypatch.setattr(model, "Setup", _FakeSetup)
    monkeypatch.setattr(model, "PolarizationState", _FakePolarizationState)

    phi = np.linspace(0, np.pi, 3) * ureg.radian
    datas = compute_scattering_amplitudes(
        wavelength=500 * ureg.nanometer,
        diameters=[100, 200] * ureg.nanometer,
        material=1.5,
        medium=1.0,
        phi=phi,
    )

    assert len(datas) == 2
    assert all(isinstance(data, ScatteringData) for data in datas)
    assert datas.k == 2.0 / ureg.meter
    assert np.allclose(datas[0].S1.magnitude, phi.magnitude + np.pi / 2)


def test_scattering_data_plot_returns_figure():
    data = ScatteringData(
        S1=np.array([1.0, 2.0]) * ureg.dimensionless,
        S2=np.array([2.0, 1.0]) * ureg.dimensionless,
        k=1.0 / ureg.meter,
        Csca=1.0 * ureg.meter**2,
        phi=np.array([0.0, np.pi]) * ureg.radian,
    )

    figure = data.plot()
    assert len(figure.axes) == 1
