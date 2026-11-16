"""
RSA generation scaling
======================

This benchmark measures the time to generate periodic random sequential
adsorption (RSA) configurations as the box volume, and therefore the requested
particle count, increases. The radius, target packing fraction, and stopping
criteria remain fixed.

The timings are specific to the machine and build configuration. They describe
RSA insertion cost, not equilibrium sampling or agreement with an analytical
hard-sphere model.
"""

from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from PackLab import monte_carlo, samplers, ureg


# %%
# Define a fixed RSA workload
# ---------------------------

radius = 100 * ureg.nanometer
target_packing_fraction = 0.08
box_lengths = np.linspace(3.0, 5.0, 10) * ureg.micrometer
sampler = samplers.ConstantRadiusSampler(radius=radius)


# %%
# Time configuration generation at increasing box sizes
# ------------------------------------------------------

sphere_counts = []
elapsed_seconds = []

for index, box_length in enumerate(box_lengths):
    domain = monte_carlo.PackingDomain(
        box_length,
        box_length,
        box_length,
        use_periodic_boundaries=True,
    )
    options = monte_carlo.RSAOptions()
    options.random_seed = 2026 + index
    options.maximum_attempts = 50_000
    options.maximum_consecutive_rejections = 10_000
    options.target_packing_fraction = target_packing_fraction

    started = perf_counter()
    result = monte_carlo.RSASimulator(domain, sampler, options).run()
    elapsed_seconds.append(perf_counter() - started)
    sphere_counts.append(result.statistics.sphere_count)
    print(
        f"{box_length.to('micrometer').magnitude:.1f} µm box: "
        f"{sphere_counts[-1]:4d} spheres in {elapsed_seconds[-1]:.3f} s"
    )


# %%
# Display generation time against the realised workload
# ------------------------------------------------------

figure, axis = plt.subplots(figsize=(6, 4))
_ = axis.plot(sphere_counts, elapsed_seconds, "o-", color="#d97706")
_ = axis.set(
    xlabel="accepted sphere count",
    ylabel="wall-clock time (s)",
    title="Periodic RSA generation scaling",
)
axis.grid(alpha=0.25)
figure.tight_layout()
