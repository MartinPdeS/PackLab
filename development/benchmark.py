from PackLab.binary.interface_rsa import SimulationDomainBox, RandomSequentialAdditionOptions, RandomSequentialAdditionSimulator

from PackLab.binary.interface_radius_sampler import UniformRadiusSampler
from PackLab.PackLab.simulator import RSASimulator



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
    simulation_domain_box = SimulationDomainBox(
        length_x=4.0,
        length_y=4.0,
        length_z=4.0,
        use_periodic_boundaries=True
    )

    radius_sampler = UniformRadiusSampler(minimum_radius=0.05, maximum_radius=0.05)

    options = RandomSequentialAdditionOptions()
    options.random_seed = 123
    options.maximum_attempts = 500_000
    options.maximum_consecutive_rejections = 50_000
    options.target_packing_fraction = 0.20
    options.minimum_center_separation_addition = 0.0

    rsa_simulator = RSASimulator(
        simulation_domain_box=simulation_domain_box,
        radius_sampler=radius_sampler,
        options=options
    )

    rsa_result = rsa_simulator.run()

    print(rsa_result.summary())