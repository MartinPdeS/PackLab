from TypedUnit.units import ureg, Length, RefractiveIndex, Angle
from PackLab.scattering.data import ScatteringDataset, ScatteringData
import numpy as np

Gaussian = Sphere = Setup = PolarizationState = None


def compute_scattering_amplitudes(
        wavelength: Length,
        diameters: Length,
        material: RefractiveIndex,
        medium: RefractiveIndex,
        phi: Angle,
        plot: bool = False,
        polarization: float = 0 * ureg.degree,
        debug_mode: bool = False
    ) -> ScatteringDataset:
    """
    Compute far field amplitude scattering functions S1 and S2 for a set of sphere diameters.

    This helper constructs a `Gaussian` source and a `Sphere` scatterer for each diameter,
    calls `Setup.get_s1s2(...)`, and stores the resulting objects in a `ScatteringDataset` container.

    Parameters
    ----------
    diameters
        Iterable of particle diameters. Each element is expected to be a Pint Quantity
        with length units (for example `6 * ureg.nanometer`).
    plot
        If True, call `data.plot(...)` for each diameter.

    Returns
    -------
    ScatteringDataset
        A list like container where each element is the return value of
        `Setup.get_s1s2(angles=...)`. Additional attributes are attached:

        * `datas.k`: optical wavenumber of the source (1/length)
        * `datas.phi`: angular sampling grid taken from the first result (radians)
        * after calling `datas.process()`, arrays `S1`, `S2`, `Csca` are available

    Notes
    -----
    This function sets two extra attributes on each `data` element:

    * `data.k` is set to `source.wavenumber`
    * `data.Csca` is set to `scatterer.Csca`

    The `ScatteringDataset` instance also receives `k` and `phi` from the last and first element
    respectively. If you want stricter correctness, you should assert that all returned
    `data.phi` grids are identical.
    """
    diameters = tuple(diameters)
    if not diameters:
        raise ValueError("diameters must contain at least one diameter")

    # Keep importing the data containers lightweight; PyMieSim is only needed when
    # this numerical helper is actually called.
    global Gaussian, Sphere, Setup, PolarizationState
    if Gaussian is None:
        from PyMieSim.single.scatterer import Sphere as _Sphere
        from PyMieSim.single.source import Gaussian as _Gaussian
        from PyMieSim.single import Setup as _Setup
        from PyMieSim.polarization import PolarizationState as _PolarizationState

        Gaussian = _Gaussian
        Sphere = _Sphere
        Setup = _Setup
        PolarizationState = _PolarizationState

    datas = ScatteringDataset()

    for diameter in diameters:
        source = Gaussian(
            wavelength=wavelength,
            polarization=PolarizationState(angle=polarization),
            optical_power=1 * ureg.watt,
            numerical_aperture=0.3,
        )

        scatterer = Sphere(
            diameter=diameter,
            medium=medium,
            material=material,
        )

        setup = Setup(source=source, scatterer=scatterer)

        s1, s2 = setup.get_s1s2(
            angles=phi + np.pi / 2 * ureg.radian
        )

        data = ScatteringData(
            S1=s1,
            S2=s2,
            k=source.wavenumber_vacuum * medium,
            Csca=setup.get("Csca"),
            phi=phi,
        )

        if plot:
            data.plot(tight_layout=False)

        datas.append(data)

        if debug_mode:
            print(f"[compute_scattering_amplitudes] Diameter: {diameter}, material: {material}, medium: {medium}, Csca: {data.Csca}")

    datas.k = data.k
    datas.phi = phi

    return datas
