
from PackLab import monte_carlo, samplers
from PackLab import ureg

domain = monte_carlo.PackingDomain(
    length_x=6.0 * ureg.millimeter,
    length_y=6.0 * ureg.millimeter,
    length_z=6.0 * ureg.millimeter,
    use_periodic_boundaries=True
)

radius_sampler = samplers.DiscreteRadiusSampler(
    radii=[0.1, 0.2] * ureg.millimeter,
    weights=[0.5, 0.5],
)

options = monte_carlo.RSAOptions()
options.random_seed = 123
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 50_000
options.target_packing_fraction = 0.3
options.minimum_center_separation_addition = 0.0

rsa_simulator = monte_carlo.RSASimulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()
