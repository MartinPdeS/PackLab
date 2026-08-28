Scattering calculations
=======================

Scattering support is an optional integration with PyMieSim. Install it with:

.. code-block:: console

   $ pip install "packlab[scattering]"

``compute_scattering_amplitudes`` evaluates the single-particle amplitudes for
a refractive-index and diameter grid. ``ScatteringDataset`` then stores those
amplitudes with their coordinates and derives scattering quantities such as a
phase function.

.. code-block:: python

   import numpy as np

   from PackLab.scattering import compute_scattering_amplitudes
   from PackLab.units import ureg

   dataset = compute_scattering_amplitudes(
       wavelength=532 * ureg.nanometer,
       diameters=[100, 200] * ureg.nanometer,
       material=1.59,
       medium=1.33,
       phi=np.linspace(-90, 90, 181) * ureg.degree,
   )
   dataset.process()

The gallery's scattering example demonstrates a complete, unit-aware workflow
and the available plotting helpers.
