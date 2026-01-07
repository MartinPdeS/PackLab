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
from TypedUnit.units import ureg
from PackLab import analytical

distribution = analytical.GaussianRadiusDistribution(
    mean_radius=1.5 * ureg.micrometer,
    standard_deviation=0.2 * ureg.micrometer,
    radius_min=0.7 * ureg.micrometer,
    radius_max=2.5 * ureg.micrometer,
    number_of_bins=4,
    bin_spacing="linear",
)

particle_radii, number_fractions = distribution.to_bins()

domain = analytical.Domain(
    size=100 * ureg.micrometer,
    particle_radii=particle_radii,
    volume_fraction=0.2,
    number_fractions=number_fractions,
)

domain.print_bins()

domain.plot_radius_distribution()

p_max = 1e3 / domain.particle_radii.min()

p = np.linspace(0, p_max, 5_000)


solver = analytical.Solver(
    densities=domain.particle_densities_per_radius,
    radii=domain.particle_radii,
    p=p,
)

distances = np.linspace(domain.particle_radii.min() * 2, domain.particle_radii.max() * 10, 1500)

result = solver.compute(distances=distances)

result.plot_pair_correlation()