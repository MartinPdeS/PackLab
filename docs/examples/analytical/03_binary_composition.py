"""
Binary-mixture composition at fixed volume fraction
===================================================

Vary the number fraction of small spheres while keeping both radii and total
volume fraction fixed. Number fraction is not volume fraction: large spheres
contribute more volume per particle. These curves are equilibrium
Percus--Yevick references, not RSA configurations.
"""

import matplotlib.pyplot as plt
import numpy as np
from PackLab import analytical, ureg

radii = np.array([0.1, 0.2]) * ureg.micrometer
distances = np.linspace(0, 1.5, 200) * ureg.micrometer
figure, axis = plt.subplots()
for small_fraction in (0.25, 0.5, 0.75):
    domain = analytical.PercusYevickDomain(
        size=100 * ureg.micrometer,
        radii=radii,
        volume_fraction=0.2,
        number_fractions=[small_fraction, 1 - small_fraction],
    )
    result = analytical.PercusYevickSolver(
        densities=domain.particle_densities_per_radius,
        radii=domain.radii,
        wavenumber="auto",
    ).compute(distances)
    axis.plot(distances.magnitude, result.g[0, 1], label=f"small fraction = {small_fraction:.2f}")
axis.set(xlabel="separation (µm)", ylabel=r"$g_{12}(r)$", title="Effect of mixture composition")
axis.legend()
axis.grid(alpha=0.2)
figure.tight_layout()
