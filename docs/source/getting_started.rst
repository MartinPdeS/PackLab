Getting started
===============

Install the base package for radius sampling, RSA simulations, and analytical
structure-factor calculations:

.. code-block:: console

   $ pip install packlab

Scattering calculations additionally require PyMieSim:

.. code-block:: console

   $ pip install "packlab[scattering]"

PackLab uses Pint-compatible quantities. Give all dimensional inputs units;
results retain them as well.

Minimal RSA simulation
----------------------

The following creates a small periodic packing of spheres with uniformly
sampled radii. It is deliberately small so it can serve as a first smoke test;
increase the domain size and proposal budget for production work.

.. code-block:: python

   from PackLab import monte_carlo, samplers
   from PackLab.units import ureg

   domain = monte_carlo.PackingDomain(
       length_x=4 * ureg.micrometer,
       length_y=4 * ureg.micrometer,
       length_z=4 * ureg.micrometer,
       use_periodic_boundaries=True,
   )
   sampler = samplers.UniformRadiusSampler(
       minimum_radius=80 * ureg.nanometer,
       maximum_radius=120 * ureg.nanometer,
   )
   options = monte_carlo.RSAOptions()
   options.maximum_attempts = 20_000
   options.target_packing_fraction = 0.08

   result = monte_carlo.RSASimulator(domain, sampler, options).run()
   print(result.statistics.packing_fraction_geometry)

``PackingResult`` stores the accepted centres, sampled radii, and derived
statistics. See :doc:`monte_carlo` for inspecting and plotting a result, and
the :doc:`gallery/index` for complete runnable examples.

Choosing a workflow
-------------------

Use :mod:`PackLab.monte_carlo` when you need an explicit configuration or a
specific radius distribution. Use :mod:`PackLab.analytical` for fast
Percus--Yevick predictions over a wave-vector grid. The two workflows can be
compared through their pair correlations or structure factors, but they are not
interchangeable models.
