import numpy
import numpy as np

from PackLab import monte_carlo, samplers


domain = monte_carlo.PackingDomain(
    length_x=6.0,
    length_y=6.0,
    length_z=6.0,
    use_periodic_boundaries=True
)

domain.scale(10)

radius_sampler = samplers.DiscreteRadiusSampler(
    radii=[1, 2],
    weights=[0.5, 0.5],
)

options = monte_carlo.RSAOptions()
options.random_seed = 123
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 500_000
options.target_packing_fraction = 0.2
options.minimum_center_separation_addition = 0.0
options.enforce_radii_distribution = True

rsa_simulator = monte_carlo.RSASimulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()

result.statistics.print()


centers, g_matrix = result.binding.compute_partial_pair_correlation_function(
    number_of_distance_bins=300,
    maximum_pairs=200_000,
)

print("Centers shape:", centers.shape)
