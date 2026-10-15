Understanding results
=====================

PackLab keeps physical inputs, generated configurations, and derived
correlations separate. This page maps the result objects returned by each
public workflow.

RSA and Metropolis results
--------------------------

``RSASimulator.run()`` and ``MetropolisSimulator.run()`` return a
``PackingResult`` representing one explicit sphere configuration:

.. code-block:: python

   result = monte_carlo.RSASimulator(domain, sampler, options).run()
   positions = result.positions       # shape: (number_of_spheres, 3)
   radii = result.radii               # shape: (number_of_spheres,)
   statistics = result.statistics
   centers, g_ij = result.compute_partial_pair_correlation_function(n_bins=80)

``positions`` and ``radii`` identify the actual packing. ``statistics``
contains summary quantities such as sphere count and geometric packing
fraction. ``centers`` carries length units and ``g_ij`` has shape
``(number_of_classes, number_of_classes, number_of_bins)``.

The plotting helpers return Matplotlib figures, so they can be labelled,
saved, or embedded in a larger figure:

.. code-block:: python

   figure = result.plot_slice_2d(show=False)
   figure.savefig("packing-slice.png", dpi=200)

Percus--Yevick results
----------------------

``PercusYevickSolver.compute(...)`` returns the analytical result evaluated on
the requested distance and wavenumber grids:

.. code-block:: python

   py_result = solver.compute(distances)
   distances = py_result.distances
   wavenumber = py_result.wavenumber
   partial_g = py_result.g
   total_correlation = py_result.H

For a mixture with :math:`K` size classes, ``g`` is indexed as
``g[i, j, distance_index]`` and ``H`` as ``H[i, j, wavenumber_index]``. Match
the class ordering to the radii and number fractions supplied to the domain.
The result is an analytical equilibrium reference, not a generated sphere
configuration.

Scattering results
------------------

``compute_scattering_amplitudes`` returns a ``ScatteringDataset`` containing
one item per requested diameter. Call ``process()`` before using mixture-level
arrays:

.. code-block:: python

   dataset = scattering.compute_scattering_amplitudes(...)
   dataset.process()
   cross_sections = dataset.Csca
   phi, theta, phase_function = dataset.get_phase_function(
       densities=py_result.densities,
       H=py_result.H,
       wavenumber=py_result.wavenumber,
   )

The phase function combines optical amplitudes with the supplied analytical
correlation tensor. It therefore inherits the assumptions of both the optical
model and the chosen hard-sphere structure model.
