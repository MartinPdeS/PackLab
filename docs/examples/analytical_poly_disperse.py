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

distribution = samplers.LogNormal(
    median_radius=1.5 * ureg.micrometer,
    geometric_standard_deviation=1.2,
    maximum_radius_clip=2.5 * ureg.micrometer,
    bins=4,
)

particle_radii, number_fractions = distribution.to_bins()

domain = analytical.Domain(
    size=100 * ureg.micrometer,
    particle_radii=particle_radii,
    volume_fraction=0.2,
    number_fractions=number_fractions,
)

domain.print_bins()

p_max = 1e3 / domain.particle_radii.min()

p = np.linspace(p_max / 1e8, p_max / 2, 10_000)

solver = analytical.Solver(
    densities=domain.particle_densities_per_radius,
    radii=domain.particle_radii,
    p=p,
)

distances = np.linspace(domain.particle_radii.min() * 2, domain.particle_radii.max() * 10, 1500)

result = solver.compute(distances=distances)

result.plot_pair_correlation()