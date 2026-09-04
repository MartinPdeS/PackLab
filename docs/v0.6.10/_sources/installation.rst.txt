Installation and reproducibility
================================

Installation
------------

Install the base package for radius samplers, RSA, Metropolis Monte Carlo, and
Percus--Yevick calculations:

.. code-block:: console

   $ pip install packlab

The optional scattering workflow requires PyMieSim:

.. code-block:: console

   $ pip install "packlab[scattering]"

When using the PackLab Conda channel, install the compiled package with its
dependencies from conda-forge:

.. code-block:: console

   $ conda install -c martinpdes -c conda-forge packlab

Verify the installation and record the version used for a calculation:

.. code-block:: console

   $ python -c "import PackLab; print(PackLab.__version__)"

What to report
--------------

For a reproducible packing, correlation, or scattering calculation, report:

* PackLab version and installation method;
* Python version, operating system, and—for native builds—the compiler;
* radii or radius-distribution parameters, number fractions, and volume
  fraction;
* domain dimensions and whether periodic boundaries were used;
* random seeds and RSA or Metropolis stopping and proposal settings;
* wavenumber-grid settings for Percus--Yevick calculations; and
* optical wavelength, refractive indices, angular grid, and PyMieSim version
  for scattering calculations.

Creating a release
------------------

For maintainers, the release helper updates the version in ``pyproject.toml``,
``meta.yaml``, ``.zenodo.json``, and the generated package version file before
creating an annotated Git tag:

.. code-block:: console

   $ make tag v0.7.0
   $ git push origin HEAD --tags

Run it from a clean worktree. The helper accepts only release tags in the form
``vMAJOR.MINOR.PATCH``.

Citing PackLab
--------------

Use the repository's ``CITATION.cff`` for software citation and retain the
version number in the methods section of a scientific report. The references
page lists the underlying hard-sphere literature and the release metadata.
