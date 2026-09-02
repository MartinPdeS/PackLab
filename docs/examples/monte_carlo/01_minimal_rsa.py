"""
Minimal random sequential adsorption
====================================

This example creates a small periodic RSA packing from a uniform radius
distribution. It then visualises a central two-dimensional slice of the
accepted sphere configuration.
"""

from PackLab import monte_carlo, samplers
from PackLab.units import ureg


domain = monte_carlo.PackingDomain(
    length_x=4 * ureg.micrometer,
    length_y=4 * ureg.micrometer,
    length_z=4 * ureg.micrometer,
    use_periodic_boundaries=True,
)

radius_sampler = samplers.UniformRadiusSampler(
    minimum_radius=80 * ureg.nanometer,
    maximum_radius=120 * ureg.nanometer,
    bins=8,
)

options = monte_carlo.RSAOptions()
options.random_seed = 12
options.maximum_attempts = 20_000
options.maximum_consecutive_rejections = 5_000
options.target_packing_fraction = 0.08

result = monte_carlo.RSASimulator(domain, radius_sampler, options).run()

print(result)
print(f"Packing fraction: {result.statistics.packing_fraction_geometry:.3f}")

figure = result.plot_slice_2d(show=False)
figure.suptitle("Central slice of a periodic RSA packing")
