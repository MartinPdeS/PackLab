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
=========


Testing
*******

To test localy (with cloning the GitHub repository) you'll need to install the dependencies and run the coverage command as

.. code:: python

   >>> git clone https://github.com/MartinPdeS/PackLab.git
   >>> cd PackLab
   >>> pip install -r requirements/requirements.txt
   >>> pytest

----


Coding example
**************

.. code-block:: python

   import PackLab




----


Contact Information
************************
As of 2025, the project is still under development. If you want to collaborate, it would be a pleasure! I encourage you to contact me.

PackLab was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet.de.sivry@gmail.com?subject=PackLab>`_ .

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
