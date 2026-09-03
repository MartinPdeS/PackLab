|logo|

PackLab
=======

**PackLab computes structure in three-dimensional hard-sphere systems.** Its
primary analytical workflow evaluates the Percus--Yevick (PY) approximation
for equilibrium mixtures. It also generates explicit random sequential
adsorption (RSA) configurations and samples fixed-volume equilibrium
configurations with Metropolis Monte Carlo (MC).

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Project
     - Status
   * - Package
     - |PyPI| |anaconda| |python|
   * - Documentation and tests
     - |docs| |ci/cd| |coverage|
   * - Citation
     - |doi|

The workflows share physical inputs such as particle radii, number fractions,
and volume fraction, but they answer different questions:

* **PY** is a fast analytical equilibrium reference for pair correlations,
  structure factors, and structure-aware scattering calculations.
* **RSA** creates an explicit, non-overlapping deposition configuration. It is
  irreversible and retains the history of accepted particles.
* **Metropolis MC** moves particles in a valid configuration to sample an
  equilibrium hard-sphere system at fixed volume, particle count, and radii.

Installation
------------

Install the core package from PyPI:

.. code-block:: bash

   pip install packlab

For scattering calculations, install the optional PyMieSim integration:

.. code-block:: bash

   pip install "packlab[scattering]"

Or install the Conda package:

.. code-block:: bash

   conda install -c martinpdes packlab

Verify the compiled package with:

.. code-block:: bash

   python -c "import PackLab; print(PackLab.__version__)"

1. Compute a Percus--Yevick equilibrium reference
--------------------------------------------------

Use PY when you need equilibrium mixture correlations without generating an
explicit packing. PackLab computes the partial pair correlations
:math:`g_{ij}(r)` and reciprocal-space correlations on an automatically
resolved wavenumber grid.

.. image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/master/docs/images/readme_percus_yevick.png
   :alt: Partial pair correlations of a binary Percus--Yevick hard-sphere mixture.
   :width: 62%
   :align: center

.. code-block:: python

   import numpy as np

   from PackLab import analytical, ureg

   radii = np.array([75, 140]) * ureg.nanometer
   domain = analytical.PercusYevickDomain(
       size=50 * ureg.micrometer,
       radii=radii,
       volume_fraction=0.25,
       number_fractions=np.array([0.7, 0.3]),
   )
   distances = np.linspace(0.0, 1.5, 300) * ureg.micrometer
   result = analytical.PercusYevickSolver(
       densities=domain.particle_densities_per_radius,
       radii=domain.radii,
       wavenumber="auto",
   ).compute(distances)

   g_12 = result.g[0, 1]
   wavenumber = result.wavenumber

The automatic grid is a useful default. For a resolution study, use
``analytical.make_wavenumber_grid(...)`` and compare the resulting curves.
PY is an analytical approximation to an equilibrium hard-sphere mixture; it
does not create particle centres or reproduce the irreversible RSA process.

2. Generate an explicit RSA packing
------------------------------------

RSA proposes particles one at a time and keeps only non-overlapping proposals.
Accepted particles never move, so the final configuration is physically valid
but history-dependent rather than an equilibrium sample.

.. image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/master/docs/images/readme_rsa_packing.png
   :alt: Two-dimensional slice through a periodic random sequential adsorption packing.
   :width: 62%
   :align: center

.. code-block:: python

   from PackLab import monte_carlo, samplers, ureg

   domain = monte_carlo.PackingDomain(
       5 * ureg.micrometer,
       5 * ureg.micrometer,
       5 * ureg.micrometer,
       use_periodic_boundaries=True,
   )
   sampler = samplers.UniformRadiusSampler(
       90 * ureg.nanometer,
       170 * ureg.nanometer,
       bins=8,
   )
   options = monte_carlo.RSAOptions()
   options.random_seed = 42
   options.maximum_attempts = 40_000
   options.target_packing_fraction = 0.12

   rsa_result = monte_carlo.RSASimulator(domain, sampler, options).run()
   print(rsa_result.statistics.packing_fraction_geometry)
   figure = rsa_result.plot_slice_2d(show=False)

``PackingResult`` provides accepted centres, sampled radii, packing
statistics, pair-correlation estimators, and plotting helpers. Radius samplers
support constant, uniform, normal, log-normal, and discrete distributions.

3. Equilibrate hard spheres with Metropolis MC
----------------------------------------------

Use Metropolis MC when an equilibrium configuration is needed. It can start
from the valid RSA configuration above, but then proposes particle
displacements; particle count, radii, and class labels remain fixed.

.. image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/master/docs/images/readme_metropolis.png
   :alt: Two-dimensional slice of a hard-sphere configuration after Metropolis Monte Carlo moves.
   :width: 62%
   :align: center

.. code-block:: python

   options = monte_carlo.MetropolisOptions()
   options.random_seed = 34
   options.number_of_sweeps = 500
   options.maximum_displacement = 50 * ureg.nanometer

   simulator = monte_carlo.MetropolisSimulator(
       domain,
       rsa_result.sphere_configuration,
       options,
   )
   mc_result = simulator.run()
   print(simulator.statistics.acceptance_rate)
   figure = mc_result.plot_slice_2d(show=False)

The number of sweeps alone does not establish equilibration. Discard an
initial burn-in interval, assess autocorrelation for the quantity of interest,
and compare larger systems when finite-size effects may matter.

Choosing the right workflow
---------------------------

.. list-table::
   :widths: 22 39 39
   :header-rows: 1

   * - Workflow
     - Use it when you need
     - Important limitation
   * - PY analytical
     - Fast equilibrium pair correlations, structure factors, parameter
       sweeps, or scattering inputs.
     - It is an equilibrium approximation, not an explicit packing.
   * - RSA
     - Particle centres, radius-sampling effects, deposition history, or a
       finite non-overlapping configuration.
     - It is irreversible and is not an equilibrium sampler.
   * - Metropolis MC
     - An explicit equilibrium hard-sphere configuration at fixed volume and
       composition.
     - Equilibration, autocorrelation, and finite-size effects require checks.

Scattering
----------

The optional ``PackLab.scattering`` workflow computes optical amplitudes with
PyMieSim. Combine them with the PY correlation tensor when you need a
structure-corrected mixture phase function. See the
`scattering examples <https://martinpdes.github.io/PackLab/docs/latest/scattering.html>`_
for runnable single-particle and mixture calculations.

Documentation, validation, and citation
----------------------------------------

The `online documentation <https://martinpdes.github.io/PackLab/docs/latest/index.html>`_
contains theory, API reference, output conventions, assumptions, and
executable galleries for PY, RSA, Metropolis MC, scattering, validation, and
benchmarks.

For development:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PackLab.git
   cd PackLab
   pip install -e ".[testing,documentation]"
   pytest

If PackLab contributes to academic work, cite the archived Zenodo release you
used. Release metadata is included in ``.zenodo.json``.

.. |logo| image:: https://github.com/MartinPdeS/PackLab/raw/master/docs/images/logo.png
   :alt: PackLab logo - sphere packing and correlation curve.
.. |python| image:: https://img.shields.io/pypi/pyversions/packlab.svg
   :alt: Supported Python versions
   :target: https://www.python.org/
.. |docs| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_documentation.yml/badge.svg
   :alt: Documentation status
   :target: https://martinpdes.github.io/PackLab/
.. |ci/cd| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_coverage.yml/badge.svg
   :alt: Continuous integration status
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/python-coverage-comment-action-data/badge.svg
   :alt: Test coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |PyPI| image:: https://badge.fury.io/py/packlab.svg
   :alt: PyPI version
   :target: https://badge.fury.io/py/PackLab
.. |anaconda| image:: https://anaconda.org/martinpdes/packlab/badges/version.svg
   :alt: Anaconda version
   :target: https://anaconda.org/martinpdes/packlab
.. |doi| image:: https://zenodo.org/badge/1105416581.svg
   :alt: Cite PackLab on Zenodo
   :target: https://doi.org/10.5281/zenodo.20207570
