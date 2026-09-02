|logo|

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Badge
     - Status
   * - Python versions
     - |python|
   * - Documentation
     - |docs|
   * - Continuous integration
     - |ci/cd|
   * - Test coverage
     - |coverage|
   * - PyPI package
     - |PyPI|
   * - PyPI downloads
     - |PyPI_download|
   * - Anaconda package
     - |anaconda|
   * - Anaconda downloads
     - |anaconda_download|
   * - Latest Anaconda release
     - |anaconda_date|

PackLab
=======

**PackLab** is an open-source Python package for generating and analysing
three-dimensional hard-sphere packings. Its C++ core provides fast random
sequential adsorption (RSA), while its analytical tools implement a
Percus--Yevick model for mixture correlations and structure factors.

Use PackLab when you need an explicit non-overlapping configuration, a
reproducible packing statistic, or a fast analytical reference for validating
an RSA result.

Features
--------

* Random sequential adsorption of mono- and polydisperse spheres.
* Periodic or finite box domains with configurable stopping criteria.
* Radius samplers for constant, uniform, normal, log-normal, and discrete
  distributions.
* Pair-correlation estimates, packing statistics, and Matplotlib plots.
* Percus--Yevick mixture solver with automatic, resolution-aware wavenumber
  grids.
* Optional PyMieSim integration for scattering and phase-function workflows.
* Unit-aware quantities throughout, via ``TypedUnit``.

Installation
------------

Install the core package from PyPI:

.. code-block:: bash

   pip install packlab

Install optional scattering support:

.. code-block:: bash

   pip install "packlab[scattering]"

The conda package is also available:

.. code-block:: bash

   conda install -c martinpdes packlab

Verify that the compiled extensions are available with the interpreter you
will use for simulations:

.. code-block:: bash

   python -c "import PackLab; print(PackLab.__version__)"

First RSA packing
-----------------

Create a periodic domain, choose a radius distribution, configure the RSA
stopping conditions, and run the simulation. Dimensional inputs carry units.

.. code-block:: python

   from PackLab import monte_carlo, samplers, ureg

   domain = monte_carlo.PackingDomain(
       length_x=6 * ureg.micrometer,
       length_y=6 * ureg.micrometer,
       length_z=6 * ureg.micrometer,
       use_periodic_boundaries=True,
   )
   radii = samplers.UniformRadiusSampler(
       minimum_radius=100 * ureg.nanometer,
       maximum_radius=200 * ureg.nanometer,
       bins=12,
   )
   options = monte_carlo.RSAOptions()
   options.random_seed = 42
   options.maximum_attempts = 100_000
   options.target_packing_fraction = 0.15

   result = monte_carlo.RSASimulator(domain, radii, options).run()
   print(result.statistics.packing_fraction_geometry)
   result.plot_slice_2d()

Analytical reference
--------------------

Use the analytical solver for a fast Percus--Yevick reference. With
``wavenumber="auto"``, PackLab chooses a zero-inclusive wavenumber grid from
the particle radii and requested distance range.

.. code-block:: python

   import numpy as np

   from PackLab import analytical, ureg

   domain = analytical.PercusYevickDomain(
       size=10 * ureg.micrometer,
       radii=[100, 150] * ureg.nanometer,
       volume_fraction=0.15,
       number_fractions=[0.7, 0.3],
   )
   distances = np.linspace(0.2, 1.5, 300) * ureg.micrometer
   solver = analytical.PercusYevickSolver(
       densities=domain.particle_densities_per_radius,
       radii=domain.radii,
       wavenumber="auto",
   )
   result = solver.compute(distances)
   print(result.wavenumber)

For an explicit grid, use ``analytical.make_wavenumber_grid(...)``. PackLab
warns when a manually supplied grid is too coarse for the requested distances.

Choosing a workflow
-------------------

* Use ``PackLab.monte_carlo`` when individual centres, sampled radii, box
  boundaries, or finite-size effects are important.
* Use ``PackLab.analytical`` for fast parameter sweeps and an analytical
  correlation reference.
* Use the validation gallery examples to compare a matching RSA configuration
  against the analytical model.

Documentation and examples
--------------------------

The `online documentation <https://martinpdes.github.io/PackLab/docs/latest/index.html>`_
contains theory, API reference, and executable examples organised into
Monte-Carlo, analytical, and validation workflows.

Building from source
--------------------

For development, clone the repository and install it in editable mode. A C++20
compiler and CMake are required to build the native extensions.

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PackLab.git
   cd PackLab
   pip install -e ".[testing,documentation]"

Testing
-------

Run the test suite with:

.. code-block:: bash

   pytest

Citing PackLab
--------------

If PackLab contributes to academic work, cite the archived Zenodo release you
used. Release metadata is included in ``.zenodo.json``.

Contributing and contact
------------------------

Issues and pull requests are welcome. For questions or collaborations, contact
`Martin Poinsinet de Sivry-Houle <mailto:martin.poinsinet.de.sivry@gmail.com>`_.

.. |logo| image:: https://github.com/MartinPdeS/PackLab/raw/master/docs/images/logo.png
   :alt: PackLab logo - sphere packing and correlation curve.
.. |python| image:: https://img.shields.io/pypi/pyversions/packlab.svg
   :alt: Supported Python versions
   :target: https://www.python.org/
.. |docs| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_documentation.yml/badge.svg
   :alt: Documentation status
   :target: https://martinpdes.github.io/PackLab/
.. |PyPI| image:: https://badge.fury.io/py/packlab.svg
   :alt: PyPI version
   :target: https://badge.fury.io/py/PackLab
.. |PyPI_download| image:: https://img.shields.io/pypi/dm/PackLab?label=PyPI%20downloads
   :alt: PyPI downloads
   :target: https://pypistats.org/packages/packlab
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/python-coverage-comment-action-data/badge.svg
   :alt: Test coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |ci/cd| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_coverage.yml/badge.svg
   :alt: Continuous integration status
.. |anaconda| image:: https://anaconda.org/martinpdes/packlab/badges/version.svg
   :alt: Anaconda version
   :target: https://anaconda.org/martinpdes/packlab
.. |anaconda_download| image:: https://anaconda.org/martinpdes/packlab/badges/downloads.svg
   :alt: Anaconda downloads
   :target: https://anaconda.org/martinpdes/packlab
.. |anaconda_date| image:: https://anaconda.org/martinpdes/packlab/badges/latest_release_relative_date.svg
   :alt: Latest Anaconda release date
   :target: https://anaconda.org/martinpdes/packlab
