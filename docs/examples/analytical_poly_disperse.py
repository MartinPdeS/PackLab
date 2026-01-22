"""
Percus Yevick mixture solver workflow
=====================================

This example demonstrates the complete workflow for computing the radial
distribution function :math:`g_{ij}(r)` of a polydisperse hard sphere mixture
using a Percus Yevick style solver.

The example covers the following steps:

1. Define a polydisperse domain (radii, volume fraction, number fractions)
2. Build the Fourier grid :math:`p`
3. Construct the Percus Yevick solver
4. Compute :math:`C_{ij}(p)`, :math:`H_{ij}(p)`, :math:`h_{ij}(r)`, and :math:`g_{ij}(r)`
5. Plot all :math:`g_{ij}(r)` curves on a single figure

The main output is a figure showing all pair correlations :math:`g_{ij}(r)` for
each species pair :math:`(i, j)`.

"""

import numpy as np
from PackLab import ureg
from PackLab import analytical, samplers
import matplotlib.pyplot as plt


distribution = samplers.Discrete(
    radii=[1.5, 1.6] * ureg.micrometer,
    weights=[0.5, 0.5],
    bins=4,
)

particle_radii, number_fractions = distribution.to_bins()

domain = analytical.PYDomain(
    size=100 * ureg.micrometer,
    radii=particle_radii,
    volume_fraction=0.3,
    number_fractions=number_fractions,
)

domain.print_bins()

p_max = 1e3 / domain.radii.min()

p = np.linspace(p_max / 1e8, p_max / 2, 20_000)

solver = analytical.Solver(
    densities=domain.particle_densities_per_radius,
    radii=domain.radii,
    p=p,
)

distances = np.linspace(domain.radii.min() * 2, domain.radii.max() * 2, 1500)

result = solver.compute(distances=distances)


figure, ax = plt.subplots(1, 1)

r = result.distances.magnitude
n = result.g.shape[0]


for i in range(n):
    for j in range(n):
        ax.plot(
            result.distances.magnitude,
            result.g[i, j, :],
            label=f"g[{i},{j}]"
        )

ax.set_xlabel("r")
ax.set_ylabel("g(r)")
ax.set_title("Radial distribution functions")

ax.legend()

plt.show()
