"""
Monte Carlo RSA vs Percus Yevick
================================

This example compares the partial pair correlation functions :math:`g_{ij}(r)` obtained from:

1. A Monte Carlo Random Sequential Addition (RSA) packing in 3D (PackLab.monte_carlo)
2. The analytical Percus Yevick solution for a polydisperse hard sphere mixture (PackLab.analytical)

We use a two class discrete radius distribution with radii :math:`r = 1` and :math:`r = 2`,
sampled with equal number fractions. The Monte Carlo simulation is configured to
enforce the target radii distribution during packing, which is useful when larger
particles are rejected more often.

The workflow is:

1. Configure a periodic cubic domain, a discrete radius sampler, and RSA options
2. Run RSA and compute partial :math:`g_{ij}(r)` from the resulting configuration
3. Build the matching analytical polydisperse domain and solve Percus Yevick
4. Plot Monte Carlo curves and overlay the analytical solution
"""

import numpy as np
import matplotlib.pyplot as plt

import PackLab
from PackLab.monte_carlo import Options, Simulator, DiscreteRadiusSampler
from PackLab.analytical import PercusYevickSolver, PolydisperseDomain
from PackLab.analytical.distributions import DiscreteRadiusDistribution
from TypedUnit import ureg


# %%
# Monte Carlo RSA setup
# ---------------------
# We use a periodic cubic domain and a two radius discrete sampler.

domain = PackLab.monte_carlo.Domain(
    length_x=40.0,
    length_y=40.0,
    length_z=40.0,
    use_periodic_boundaries=True,
)

radius_sampler = DiscreteRadiusSampler(
    radii=[1.0, 2.0],
    weights=[0.5, 0.5],
)

options = Options()
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 500_000
options.target_packing_fraction = 0.24
options.minimum_center_separation_addition = 0.0
options.enforce_radii_distribution = True

rsa_simulator = Simulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options,
)

result = rsa_simulator.run()
result.statistics.print()


mc_centers, mc_g_ij = result.compute_partial_pair_correlation_function(
    n_bins=1_000,
    maximum_pairs=0,
)


# %%
# Analytical Percus Yevick setup
# ------------------------------
# We construct an analytical polydisperse domain matching the Monte Carlo mixture.
# The analytical domain uses Pint quantities.

distribution = DiscreteRadiusDistribution(
    particle_radii=[1.0, 2.0] * ureg.micrometer,
    number_fractions=[1.0, 1.0],
)

particle_radii, number_fractions = distribution.to_bins()

py_domain = PolydisperseDomain(
    size=100_000 * ureg.micrometer,
    particle_radii=particle_radii,
    volume_fraction=0.24,
    number_fractions=number_fractions,
)

# Percus Yevick solver radial frequency grid
p_max = 1e3 / py_domain.particle_radii.min()
p = np.linspace(0, p_max * 5, 30_000)

solver = PercusYevickSolver(
    densities=py_domain.particle_densities_per_radius,
    radii=py_domain.particle_radii,
    p=p,
)

distances = np.linspace(
    py_domain.particle_radii.min() * 2,
    py_domain.particle_radii.max() * 10,
    400,
)

py_result = solver.compute(distances=distances)


# %%
# Compare Monte Carlo and analytical g_ij(r)
# -----------------------------------------
# We plot all partial curves on the same axes and overlay the Percus Yevick result in black.
fig, ax = plt.subplots(1, 1)

K = 2
for i in range(K):
    for j in range(K):
        ax.plot(mc_centers, mc_g_ij[i, j], label=rf"RSA $g_{{{i}{j}}}(r)$")
        ax.plot(py_result.distances, py_result.g[i, j], color="black", linewidth=1.5)

ax.set_xlabel("r")
ax.set_ylabel(r"$g_{ij}(r)$")
ax.set_title("Partial pair correlation: RSA vs Percus Yevick")
ax.legend()
plt.show()
