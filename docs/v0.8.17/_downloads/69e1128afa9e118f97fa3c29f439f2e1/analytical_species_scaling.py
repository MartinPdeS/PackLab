"""
Percus--Yevick mixture-size scaling
====================================

Measure the cost of evaluating mixtures with different numbers of radius
classes on fixed radial and wavenumber grids. The physical size distribution
changes with the number of classes, so this benchmark measures computational
scaling rather than a numerical convergence sequence. The manual grid uses
12 samples per longest-distance oscillation and a 0.01 µm radial resolution;
verify convergence separately before using it in a scientific calculation.

The reported median of seven runs includes result construction. Timings depend
on hardware, compiler, package versions and thread settings.
"""
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from PackLab import analytical, ureg


distances = np.linspace(0.0, 1.5, 150) * ureg.micrometer
wavenumber = analytical.make_wavenumber_grid(
    radial_resolution=0.01 * ureg.micrometer,
    maximum_distance=distances.max(),
    samples_per_oscillation=12,
)
class_counts = (1, 2, 4, 8)
medians = []
for count in class_counts:
    radii = np.linspace(0.1, 0.2, count) * ureg.micrometer
    domain = analytical.PercusYevickDomain(
        size=100 * ureg.micrometer,
        radii=radii,
        volume_fraction=0.2,
        number_fractions=np.full(count, 1.0 / count),
    )
    solver = analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber=wavenumber,
    )
    solver.compute(distances)  # warm-up
    timings = []
    for _ in range(7):
        start = perf_counter()
        solver.compute(distances)
        timings.append(perf_counter() - start)
    medians.append(np.median(timings))
    print(f"{count} classes, {len(wavenumber)} wavenumbers: {medians[-1]:.4f} s")

figure, axis = plt.subplots()
axis.plot(class_counts, medians, "o-")
axis.set(
    xlabel="radius classes",
    ylabel="median wall time (s)",
    title="Percus--Yevick mixture-size scaling",
)
axis.grid(alpha=0.2)
figure.tight_layout()
