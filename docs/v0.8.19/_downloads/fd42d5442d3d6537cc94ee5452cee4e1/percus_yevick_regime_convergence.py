"""
Grid refinement across mixture regimes
======================================

Check how the automatic Percus--Yevick inverse-transform grid responds to
changes in volume fraction and size ratio. The reported maximum difference is
measured away from hard-core contacts, where the discontinuity makes pointwise
Fourier convergence slower. This assesses numerical resolution of one model;
it does not imply that irreversible RSA is equivalent to the equilibrium
Percus--Yevick reference.
"""

import matplotlib.pyplot as plt
import numpy as np
from PackLab import analytical, ureg

small_radius = 0.1
ratios = (1.5, 2.0, 3.0)
volume_fractions = (0.05, 0.2, 0.35)
distances = np.linspace(0, 1.5, 240) * ureg.micrometer
errors = np.empty((len(volume_fractions), len(ratios)))
for row, fraction in enumerate(volume_fractions):
    for column, ratio in enumerate(ratios):
        radii = np.array([small_radius, small_radius * ratio]) * ureg.micrometer
        domain = analytical.PercusYevickDomain(
            size=100 * ureg.micrometer,
            radii=radii,
            volume_fraction=fraction,
            number_fractions=[0.5, 0.5],
        )
        kwargs = dict(densities=domain.particle_densities_per_radius, radii=domain.radii)
        automatic = analytical.PercusYevickSolver(**kwargs, wavenumber="auto").compute(distances)
        grid = analytical.make_wavenumber_grid(
            radial_resolution=radii.min() / 20,
            maximum_distance=distances.max(),
            samples_per_oscillation=24,
        )
        refined = analytical.PercusYevickSolver(**kwargs, wavenumber=grid).compute(distances)
        differences = []
        for i, j in ((0, 0), (0, 1), (1, 1)):
            contact = (radii[i] + radii[j]).to("micrometer").magnitude
            mask = distances.magnitude > contact + 0.04
            differences.append(np.max(np.abs(automatic.g[i, j, mask] - refined.g[i, j, mask])))
        errors[row, column] = max(differences)
        print(
            f"volume fraction={fraction:.2f}, size ratio={ratio:.1f}: "
            f"max change={errors[row, column]:.3g}"
        )

figure, axis = plt.subplots()
image = axis.imshow(errors, origin="lower", aspect="auto")
axis.set(
    xticks=range(len(ratios)),
    xticklabels=ratios,
    yticks=range(len(volume_fractions)),
    yticklabels=volume_fractions,
    xlabel="large/small radius",
    ylabel="volume fraction",
    title="Maximum change under grid refinement",
)
figure.colorbar(image, ax=axis, label="maximum absolute change in partial g")
figure.tight_layout()
