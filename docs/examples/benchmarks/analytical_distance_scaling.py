"""
Analytical evaluation scaling
=============================

This benchmark measures the wall-clock time required to evaluate a fixed
binary Percus--Yevick mixture at increasingly dense real-space distance grids.
It measures computational cost only; it is not a numerical-convergence or
physical-validation test.

Run it on the target machine and report the Python environment, processor,
and the median of several repetitions when comparing results between systems.
"""

from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from PackLab import analytical, ureg


# %%
# Define a fixed physical mixture
# -------------------------------

radii = np.array([0.75, 1.5]) * ureg.micrometer
domain = analytical.PercusYevickDomain(
    size=100_000 * ureg.micrometer,
    radii=radii,
    volume_fraction=0.20,
    number_fractions=np.array([0.5, 0.5]),
)
solver = analytical.PercusYevickSolver(
    densities=domain.particle_densities_per_radius,
    radii=domain.radii,
    wavenumber="auto",
)


# %%
# Time evaluations at different distance-grid sizes
# -------------------------------------------------

point_counts = np.array([120, 180, 240, 360, 480, 720, 960, 1440, 1920, 2880])
elapsed_seconds = []

for point_count in point_counts:
    distances = np.linspace(0.0, 12.0, point_count) * ureg.micrometer
    started = perf_counter()
    result = solver.compute(distances)
    elapsed_seconds.append(perf_counter() - started)
    print(
        f"{point_count:4d} distances: {elapsed_seconds[-1]:.3f} s "
        f"({len(result.wavenumber)} wavenumber points)"
    )


# %%
# Display the workload and timing result
# ---------------------------------------

figure, axis = plt.subplots(figsize=(6, 4))
axis.loglog(point_counts, elapsed_seconds, "o-", color="#1677b8")
axis.set(
    xlabel="requested real-space distance points",
    ylabel="wall-clock time (s)",
    title="Percus--Yevick evaluation scaling",
)
axis.grid(alpha=0.25, which="both")
figure.tight_layout()
