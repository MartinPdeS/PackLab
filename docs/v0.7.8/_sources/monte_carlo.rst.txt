Monte-Carlo packings
====================

The Monte-Carlo namespace implements random sequential adsorption (RSA) of
non-overlapping spheres. Its workflow is intentionally explicit: define the
box, define how radii are drawn, choose stopping conditions, then run.

Configure the domain and sampler
--------------------------------

``PackingDomain`` takes the three box widths and a ``periodic`` flag. Radius
samplers accept dimensional radii and include constant, uniform, normal,
log-normal, and discrete distributions.

.. code-block:: python

   from PackLab import monte_carlo, samplers
   from PackLab.units import ureg

   domain = monte_carlo.PackingDomain(
       length_x=10 * ureg.micrometer,
       length_y=10 * ureg.micrometer,
       length_z=10 * ureg.micrometer,
       use_periodic_boundaries=True,
   )
   radii = samplers.LogNormalRadiusSampler(
       median_radius=150 * ureg.nanometer,
       geometric_standard_deviation=1.15,
       maximum_radius_clip=250 * ureg.nanometer,
   )

Run RSA
-------

``maximum_attempts`` controls how long RSA continues after unsuccessful
proposals. ``target_packing_fraction`` is optional and stops the simulation as
soon as the requested fraction is reached.

.. code-block:: python

   options = monte_carlo.RSAOptions()
   options.maximum_attempts = 100_000
   options.target_packing_fraction = 0.15
   result = monte_carlo.RSASimulator(domain, radii, options).run()

Inspecting results
------------------

``PackingResult`` exposes the accepted configuration and helper methods for
visualisation and statistics. For example:

.. code-block:: python

   print(result)
   figure = result.plot_slice_2d()

For a radial pair-correlation estimate, build a ``PackingStatistics`` instance
from the result and choose bins appropriate to the box dimensions. Periodic
domains are generally the right choice for bulk correlation estimates because
they avoid treating the box faces as physical boundaries.

Equilibrating a configuration with Metropolis moves
---------------------------------------------------

RSA is an irreversible deposition process. To sample a fixed-volume
equilibrium hard-sphere system, use its valid configuration only as the
initial state of ``MetropolisSimulator``. The sampler proposes symmetric,
single-particle displacements and accepts a proposal when it creates no hard
sphere overlap. It preserves the number of particles, radii, and class labels.

.. code-block:: python

   metropolis_options = monte_carlo.MetropolisOptions()
   metropolis_options.random_seed = 42
   metropolis_options.number_of_sweeps = 1_000
   metropolis_options.maximum_displacement = 50 * ureg.nanometer

   simulator = monte_carlo.MetropolisSimulator(
       domain,
       result.sphere_configuration,
       metropolis_options,
   )
   equilibrated_result = simulator.run()
   print(simulator.statistics.acceptance_rate)

The displacement controls the acceptance rate and should be tuned for the
density and particle sizes. A finite number of sweeps alone does not establish
equilibration: assess stationarity and collect independent samples for the
observable of interest.
