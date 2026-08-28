from PackLab import monte_carlo, samplers, ureg

import time
from contextlib import contextmanager

@contextmanager
def tic_toc(label: str = "Block"):
    start_time = time.perf_counter()
    try:
        yield
    finally:
        elapsed_time_seconds = time.perf_counter() - start_time
        print(f"{label}: {elapsed_time_seconds:.6f} s")

with tic_toc("RSA simulation"):
    domain = monte_carlo.PackingDomain(
        length_x=4.0 * ureg.meter,
        length_y=4.0 * ureg.meter,
        length_z=4.0 * ureg.meter,
        use_periodic_boundaries=True
    )

    radius_sampler = samplers.UniformRadiusSampler(
        minimum_radius=0.05 * ureg.meter,
        maximum_radius=0.05 * ureg.meter,
    )

    options = monte_carlo.RSAOptions()
    options.random_seed = 123
    options.maximum_attempts = 500_000
    options.maximum_consecutive_rejections = 50_000
    options.target_packing_fraction = 0.20
    options.minimum_center_separation_addition = 0.0

    rsa_simulator = monte_carlo.RSASimulator(
        domain=domain,
        radius_sampler=radius_sampler,
        options=options
    )

    rsa_result = rsa_simulator.run()

    rsa_result.statistics.print()
