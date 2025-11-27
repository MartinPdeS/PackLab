"""
Workflow example
----------------

Docstring for PackLab.docs.examples.workflow
"""



from PackLab import Domain, Options

from PackLab import Simulator, UniformRadiusSampler


domain = Domain(
    length_x=3.0,
    length_y=3.0,
    length_z=3.0,
    use_periodic_boundaries=True
)

radius_sampler = UniformRadiusSampler(minimum_radius=0.1, maximum_radius=0.1)

options = Options()
options.random_seed = 123
options.maximum_attempts = 1500_000
options.maximum_consecutive_rejections = 50_000
options.target_packing_fraction = 0.80
options.minimum_center_separation_addition = 0.0


rsa_simulator = Simulator(
    domain=domain,
    radius_sampler=radius_sampler,
    options=options
)

result = rsa_simulator.run()

result.statistics.print()

result.plot_slice_2d()