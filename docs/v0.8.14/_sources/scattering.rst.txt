Scattering calculations
=======================

Scattering support is an optional integration with PyMieSim. Install it with:

.. code-block:: console

   $ pip install "packlab[scattering]"

``compute_scattering_amplitudes`` evaluates the single-particle amplitudes for
a refractive-index and diameter grid. ``ScatteringDataset`` then stores those
amplitudes with their coordinates and derives scattering quantities such as a
phase function.
Input angles are measured relative to the transverse plane; the returned
dataset stores the corresponding physical polar angles from 0 to :math:`\pi`.

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

The scattering gallery provides separate examples for single-particle
amplitudes and cross sections, and for a mixture phase function using a
Percus--Yevick structure factor.

Examples
--------

* :ref:`sphx_glr_gallery_scattering_01_single_particle_amplitudes.py`
  evaluates amplitudes and cross sections for individual sphere diameters.
* :ref:`sphx_glr_gallery_scattering_02_structure_factor_phase_function.py`
  combines a binary Percus--Yevick mixture with optical amplitudes to obtain a
  structure-corrected phase function.
