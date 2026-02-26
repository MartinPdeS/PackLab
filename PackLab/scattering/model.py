from TypedUnit.units import ureg, Length, RefractiveIndex, Angle
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

from PackLab.scattering.data import Datas
import numpy as np

def get_s1s2(
        wavelength: Length,
        diameters: Length,
        refractive_index: RefractiveIndex,
        medium_refractive_index: RefractiveIndex,
        phi: Angle,
        plot: bool = False,
        polarization: float = 0 * ureg.degree
    ) -> Datas:
    """
    Compute far field amplitude scattering functions S1 and S2 for a set of sphere diameters.

    This helper constructs a `Gaussian` source and a `Sphere` scatterer for each diameter,
    calls `scatterer.get_s1s2(...)`, and stores the resulting objects in a `Datas` container.

    Parameters
    ----------
    diameters
        Iterable of particle diameters. Each element is expected to be a Pint Quantity
        with length units (for example `6 * ureg.nanometer`).
    plot
        If True, call `data.plot(...)` for each diameter.

    Returns
    -------
    Datas
        A list like container where each element is the return value of
        `Sphere.get_s1s2(sampling=...)`. Additional attributes are attached:

        * `datas.k`: optical wavenumber of the source (1/length)
        * `datas.phi`: angular sampling grid taken from the first result (radians)
        * after calling `datas.process()`, arrays `S1`, `S2`, `Csca` are available

    Notes
    -----
    This function sets two extra attributes on each `data` element:

    * `data.k` is set to `source.wavenumber`
    * `data.Csca` is set to `scatterer.Csca`

    The `Datas` instance also receives `k` and `phi` from the last and first element
    respectively. If you want stricter correctness, you should assert that all returned
    `data.phi` grids are identical.
    """
    class Temp():
        pass

    datas = Datas()

    for diameter in diameters:
        source = Gaussian(
            wavelength=wavelength,
            polarization=polarization,
            optical_power=1 * ureg.watt,
            NA=0.3 * ureg.AU,
        )

        scatterer = Sphere(
            diameter=diameter,
            source=source,
            medium_refractive_index=medium_refractive_index,
            refractive_index=refractive_index,
        )

        data = Temp()
        data.S1, data.S2 = scatterer.get_s1s2(
            phi=phi + np.pi / 2 * ureg.radian
        )

        if plot:
            data.plot(tight_layout=False)

        data.k = source.wavenumber_vacuum * medium_refractive_index
        data.Csca = scatterer.Csca

        datas.append(data)

    datas.k = data.k
    datas.phi = phi

    return datas
