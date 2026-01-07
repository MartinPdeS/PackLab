|logo|

.. list-table::
   :widths: 10 25 25 25
   :header-rows: 0

   * - Meta
     - |python|
     - |docs|
     -
   * - Testing
     - |ci/cd|
     - |coverage|
     - |colab|
   * - PyPI
     - |PyPI|
     - |PyPI_download|
     -
   * - Anaconda
     - |anaconda|
     - |anaconda_download|
     - |anaconda_date|

PackLab
=======

PackLab is a fast Random Sequential Addition (RSA) particle packing library with a C plus plus core and Python bindings.
It helps you generate non overlapping sphere configurations inside a 3D box with optional periodic boundaries.
It supports mono disperse and poly disperse radius sampling through a flexible sampling interface.

Features
********

* Fast C++ implementation with Python bindings
* 3D box domain with optional periodic boundaries
* Configurable stopping criteria: maximum attempts, maximum consecutive rejections, target packing fraction, maximum spheres
* Radius samplers designed for future extension to poly disperse systems
* Python level convenience API with plotting utilities and summary statistics

Installation
************

From PyPI:

.. code-block:: bash

   pip install packlab

From conda:

.. code-block:: bash

   conda install -c martinpdes packlab

From source:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PackLab.git
   cd PackLab
   pip install -e .[testing]

Quick start
***********

.. code-block:: python

    from PackLab import Domain, Options, Simulator, UniformRadiusSampler

    domain = Domain(
        length_x=6.0,
        length_y=6.0,
        length_z=6.0,
        use_periodic_boundaries=True
    )

    radius_sampler = UniformRadiusSampler(
        minimum_radius=0.1,
        maximum_radius=0.4
    )

    options = Options()
    options.random_seed = 42
    options.maximum_attempts = 4_000_000
    options.maximum_consecutive_rejections = 80_000
    options.target_packing_fraction = 0.55
    options.minimum_center_separation_addition = 0.0

    rsa_simulator = Simulator(
        domain=domain,
        radius_sampler=radius_sampler,
        options=options
    )

    result = rsa_simulator.run()

    result.statistics.print()

    result.plot_slice_2d()

    result.plot_pair_correlation(maximum_number_of_pairs=3_000_000)

   import PackLab


Testing
*******

To test locally, clone the repository, install dependencies, and run pytest.

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PackLab.git
   cd PackLab
   pip install -e .[testing]
   pytest

Documentation
*************

The documentation includes tutorials, API reference, and gallery examples.

See |docs|.

Contributing
************

Issues and pull requests are welcome.
If you are using PackLab in research, citations and links to your work are appreciated.

Contact Information
*******************

As of 2025, the project is still under development.
If you want to collaborate, it would be a pleasure. Feel free to contact me.

PackLab was written by `Martin Poinsinet de Sivry Houle <https://github.com/MartinPdS>`_.

Email: `martin.poinsinet.de.sivry@gmail.com <mailto:martin.poinsinet.de.sivry@gmail.com?subject=PackLab>`_.

.. |logo| image:: https://github.com/MartinPdeS/PackLab/raw/master/docs/images/logo.png
    :alt: PackLab logo
.. |python| image:: https://img.shields.io/pypi/pyversions/packlab.svg
    :alt: Python
    :target: https://www.python.org/
.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Google Colab
    :target: https://colab.research.google.com/github/MartinPdeS/PackLab/blob/master/notebook.ipynb
.. |docs| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_documentation.yml/badge.svg
    :target: https://martinpdes.github.io/PackLab/
    :alt: Documentation Status
.. |PyPI| image:: https://badge.fury.io/py/packlab.svg
    :alt: PyPI version
    :target: https://badge.fury.io/py/PackLab
.. |PyPI_download| image:: https://img.shields.io/pypi/dm/PackLab?style=plastic&label=PyPI%20downloads&labelColor=hex&color=hex
    :alt: PyPI downloads
    :target: https://pypistats.org/packages/packlab
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PackLab/python-coverage-comment-action-data/badge.svg
    :alt: Unittest coverage
    :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |ci/cd| image:: https://github.com/martinpdes/packlab/actions/workflows/deploy_coverage.yml/badge.svg
    :alt: Unittest Status
.. |anaconda| image:: https://anaconda.org/martinpdes/packlab/badges/version.svg
    :alt: Anaconda version
    :target: https://anaconda.org/martinpdes/packlab
.. |anaconda_download| image:: https://anaconda.org/martinpdes/packlab/badges/downloads.svg
    :alt: Anaconda downloads
    :target: https://anaconda.org/martinpdes/packlab
.. |anaconda_date| image:: https://anaconda.org/martinpdes/packlab/badges/latest_release_relative_date.svg
    :alt: Latest release date
    :target: https://anaconda.org/martinpdes/packlab
