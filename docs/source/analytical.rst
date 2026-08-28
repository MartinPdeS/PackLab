Analytical structure factors
============================

The analytical namespace implements a Percus--Yevick hard-sphere model. It is
useful when you want a fast structure-factor calculation without generating an
explicit packing configuration.

.. code-block:: python

   import numpy as np

   from PackLab import analytical
   from PackLab.units import ureg

   radii = [100, 150] * ureg.nanometer
   number_fractions = [0.7, 0.3]
   domain = analytical.PercusYevickDomain(
       size=10 * ureg.micrometer,
       radii=radii,
       volume_fraction=0.15,
       number_fractions=number_fractions,
   )
   distances = np.linspace(0.2, 1.5, 300) * ureg.micrometer
   solver = analytical.PercusYevickSolver(
       densities=domain.particle_densities_per_radius,
       radii=domain.radii,
       wavenumber="auto",
   )
   result = solver.compute(distances=distances)

The solver result contains the evaluated wavenumber grid and the associated
structure factor. With ``wavenumber="auto"``, PackLab uses zero as the minimum,
one twentieth of the smallest particle radius as the real-space resolution, and
12 samples per sinc-kernel oscillation at the largest requested distance. Pass
``radial_resolution=...`` or ``samples_per_oscillation=...`` to tune that
selection. ``make_wavenumber_grid`` remains available when you need to provide
the grid explicitly. The solver emits a ``RuntimeWarning`` when an explicit
grid has fewer than eight samples per sinc-kernel oscillation.

The model assumes an idealised hard-sphere fluid. It is therefore a useful
reference for Monte-Carlo results, rather than a replacement for an RSA
configuration with finite size, boundaries, and a chosen radius sampler.
