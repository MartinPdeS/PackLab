"""
Polydisperse RSA simulation
===========================

This example demonstrates a Random Sequential Addition (RSA) simulation
with polydisperse spheres. In a polydisperse system, each sphere radius
is sampled from a continuous distribution rather than using a fixed value.

The workflow illustrates:

1. Creating a simulation domain
2. Sampling radii from a uniform distribution
3. Running the simulator with polydispersity
4. Inspecting statistics
5. Plotting a slice of the configuration
6. Computing and displaying the pair correlation function

Polydisperse systems are useful when modelling realistic physical packings,
microstructures, or biological particle systems with non uniform sizes.
"""

# %%
# Simulation domain
# -----------------
# We create a cube with periodic boundaries.

from PackLab import Domain, Options, Simulator, UniformRadiusSampler

domain = Domain(
    length_x=6.0,
    length_y=6.0,
    length_z=6.0,
    use_periodic_boundaries=True
)

# %%
# Polydisperse radius sampler
# ---------------------------
# Radii are drawn from a uniform distribution between 0.1 and 0.4.
# This gives a wide range of sphere sizes.

radius_sampler = UniformRadiusSampler(
    minimum_radius=0.1,
    maximum_radius=0.4
)

# %%
# Simulation options
# ------------------
# We increase the number of attempts since small spheres will continue
# to find space even when the large ones are already rejected.

options = Options()
options.random_seed = 42
options.maximum_attempts = 4_000_000
options.maximum_consecutive_rejections = 80_000
options.target_packing_fraction = 0.55
options.minimum_center_separation_addition = 0.0

# %%
# Run the polydisperse RSA simulation
# -----------------------------------

rsa_simulator = Simulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()

# %%
# Inspect simulation statistics
# -----------------------------

result.statistics.print()

# %%
# Plot a slice of the configuration
# ---------------------------------
# This gives a quick visual impression of the spatial structure,
# which is more complex than in the monodisperse case.

result.plot_slice_2d()

# %%
# Pair correlation function
# -------------------------
# Compute and plot the pair correlation function g(r)
# for the polydisperse configuration.

result.plot_pair_correlation(maximum_number_of_pairs=3_000_000)
