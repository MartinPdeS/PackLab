"""
RSA cost at different target fractions
=======================================

Compare insertion times at fixed box size and radius. Each point is the
median of three seeded independent runs; print achieved fractions and particle
counts because stopping limits can prevent a target from being reached.
Timings depend on the machine, build, and random seeds. RSA is irreversible.
"""

from time import perf_counter
import matplotlib.pyplot as plt
import numpy as np
from PackLab import monte_carlo, samplers, ureg

domain = monte_carlo.PackingDomain(
    3 * ureg.micrometer, 3 * ureg.micrometer, 3 * ureg.micrometer, use_periodic_boundaries=True
)
sampler = samplers.ConstantRadiusSampler(0.15 * ureg.micrometer)
targets = [0.05, 0.10, 0.15, 0.20]
timings = []
for target in targets:
    runs = []
    for seed in (11, 22, 33):
        options = monte_carlo.RSAOptions()
        options.random_seed = seed
        options.maximum_attempts = 100000
        options.maximum_consecutive_rejections = 20000
        options.target_packing_fraction = target
        start = perf_counter()
        result = monte_carlo.RSASimulator(domain, sampler, options).run()
        runs.append(perf_counter() - start)
        print(
            f"target {target:.2f}, seed {seed}: "
            f"achieved {result.statistics.packing_fraction_geometry:.4f}, "
            f"N={result.statistics.sphere_count}, time={runs[-1]:.4f} s"
        )
    timings.append(np.median(runs))
figure, axis = plt.subplots()
axis.plot(targets, timings, "o-")
axis.set(
    xlabel="target packing fraction", ylabel="median wall time (s)", title="RSA insertion cost"
)
axis.grid(alpha=0.2)
figure.tight_layout()
