"""
A full workflow example with PackLab
====================================

This example shows the complete workflow of using PackLab to perform a
Random Sequential Addition (RSA) simulation in three dimensions.

It demonstrates the following steps:

1. Create a simulation domain
2. Define a radius sampler for particle sizes
3. Configure simulation options
4. Construct and run the RSA simulator
5. Access the simulation statistics
6. Visualize the resulting configuration and pair correlation function

This is the recommended starting point when learning how to use PackLab.
"""
from PackLab import monte_carlo
from TypedUnit.units import ureg
# %%
# Simulation domain
# -----------------
# The domain defines the physical volume of the simulation.
# Here we use periodic boundary conditions on a cubic box.

domain = monte_carlo.Domain(
    length_x=6.0 * ureg.millimeter,
    length_y=6.0 * ureg.millimeter,
    length_z=6.0 * ureg.millimeter,
    use_periodic_boundaries=True
)

radius_sampler = monte_carlo.samplers.Discrete(
    radii=[0.1, 0.2] * ureg.millimeter,
    weights=[0.5, 0.5],
)

options = monte_carlo.Options()
options.random_seed = 123
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 50_000
options.target_packing_fraction = 0.50
options.minimum_center_separation_addition = 0.0

rsa_simulator = monte_carlo.Simulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()

result.statistics.print()

result.plot_slice_2d(
    slice_axis="z",
    slice_center_fraction=0.5,
    slice_thickness_fraction=0.08,
    maximum_circles_in_slice=2500,
)

result.plot_pair_correlation(
    n_bins=150,
    maximum_pairs=20_000_000
)
