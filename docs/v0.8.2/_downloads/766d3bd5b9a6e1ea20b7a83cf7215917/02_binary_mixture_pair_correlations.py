"""
Binary-mixture Percus--Yevick pair correlations
================================================

Evaluate the three distinct partial pair correlations of a binary hard-sphere
mixture.  The hard-core cutoffs occur at the corresponding contact distances
:math:`2a_1`, :math:`a_1+a_2`, and :math:`2a_2`.
"""

import matplotlib.pyplot as plt
import numpy as np

from PackLab import analytical, ureg


# %%
# Define and solve a binary equilibrium reference
# ------------------------------------------------

radii = np.array([75, 140]) * ureg.nanometer
domain = analytical.PercusYevickDomain(
    size=50 * ureg.micrometer,
    radii=radii,
    volume_fraction=0.25,
    number_fractions=np.array([0.7, 0.3]),
)
distances = np.linspace(0.0, 1.5, 300) * ureg.micrometer
result = analytical.PercusYevickSolver(
    densities=domain.particle_densities_per_radius,
    radii=domain.radii,
    wavenumber="auto",
).compute(distances)


# %%
# Compare like- and unlike-species correlations
# ----------------------------------------------

figure, axis = plt.subplots(figsize=(6.5, 4))
for i, j, label in ((0, 0, r"$g_{11}$"), (0, 1, r"$g_{12}$"), (1, 1, r"$g_{22}$")):
    _ = axis.plot(distances.to("micrometer").magnitude, result.g[i, j], label=label)

_ = axis.set(
    xlabel=r"separation $r$ ($\mu$m)",
    ylabel=r"partial pair correlation $g_{ij}(r)$",
    title="Binary Percus--Yevick hard-sphere mixture",
)
_ = axis.legend()
axis.grid(alpha=0.2)
figure.tight_layout()
