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

# %%
# Simulation domain
# -----------------
# The domain defines the physical volume of the simulation.
# Here we use periodic boundary conditions on a cubic box.

from PackLab import Domain, Options, Simulator, UniformRadiusSampler

domain = Domain(
    length_x=6.0,
    length_y=6.0,
    length_z=6.0,
    use_periodic_boundaries=True
)

# %%
# Radius sampler
# --------------
# All spheres will have radius 0.2 in this example.

radius_sampler = UniformRadiusSampler(
    minimum_radius=0.15,
    maximum_radius=0.15
)

# %%
# Simulation options
# ------------------
# These control the stochastic behaviour of the RSA algorithm.

options = Options()
options.random_seed = 123
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 50_000
options.target_packing_fraction = 0.50
options.minimum_center_separation_addition = 0.0

# %%
# Run the simulator
# -----------------
# Construct the simulator and run the RSA process.

rsa_simulator = Simulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()

# %%
# Print statistics
# ----------------
# The statistics object reports attempts, accepted insertions,
# packing fraction and other useful diagnostic information.

result.statistics.print()

# %%
# Plot a 2D slice of the configuration
# ------------------------------------
# This provides a visual view of the spatial arrangement of spheres.

result.plot_slice_2d(
    slice_axis="z",
    slice_center_fraction=0.5,
    slice_thickness_fraction=0.08,
    maximum_circles_in_slice=2500,
)

# %%
# Pair correlation function
# -------------------------
# We finally compute and plot the pair correlation function g(r)
# using Monte Carlo sampling of particle pairs.

result.plot_pair_correlation(
    pair_correlation_bins=300,
    maximum_number_of_pairs=20_000_000
)
