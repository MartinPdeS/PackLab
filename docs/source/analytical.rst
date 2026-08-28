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
   wave_vector = np.linspace(0.0, 80.0, 500) / ureg.micrometer
   solver = analytical.PercusYevickSolver(
       densities=domain.particle_densities_per_radius,
       radii=domain.radii,
       p=wave_vector,
   )
   result = solver.compute(distances=np.linspace(0.2, 1.5, 300) * ureg.micrometer)

The solver result contains the evaluated wave-vector grid and the associated
structure factor. Keep the grid compact while exploring parameters, then use
a denser grid where a plotted feature needs resolution.

The model assumes an idealised hard-sphere fluid. It is therefore a useful
reference for Monte-Carlo results, rather than a replacement for an RSA
configuration with finite size, boundaries, and a chosen radius sampler.
