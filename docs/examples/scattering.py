"""
Scattering Example: Percus Yevick Structure Factor and Phase Function
=====================================================================

Example of Percus Yevick structure factor computation for a polydisperse
mixture and corresponding phase function calculation.

"""

import matplotlib.pyplot as plt
import numpy as np

from PackLab import analytical, samplers, scattering
from TypedUnit import ureg
from PackLab.units import ureg


particle_radii = np.asarray([0.1, 5.0 / 10]) * ureg.micrometer

sampler = samplers.Discrete(
    radii=particle_radii,
    weights=[1.0 * 1000, 1.0],
)


particle_radii, number_fractions = sampler.to_bins()

py_domain = analytical.PYDomain(
    size=100 * ureg.micrometer,
    radii=particle_radii,
    volume_fraction=0.24,
    number_fractions=number_fractions,
)

py_domain.print_bins()

# Percus Yevick solver radial frequency grid
p_max = 1e3 / py_domain.radii.min()
p = np.linspace(0, p_max * 1, 2 * 60_000)

solver = analytical.Solver(
    densities=py_domain.particle_densities_per_radius,
    radii=py_domain.radii,
    p=p,
)

distances = np.linspace(
    py_domain.radii.min() * 2,
    py_domain.radii.max() * 10,
    400,
)

py_result = solver.compute(distances=distances)

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

K = 2
for i in range(K):
    for j in range(K):
        ax.plot(py_result.distances.to('micrometer'), py_result.g[i, j], linewidth=1.5, label=f'{i}-{j}')


ax.set_xlabel("r")
ax.set_ylabel(r"$g_{ij}(r)$")
ax.set_title("Partial pair correlation: RSA vs Percus Yevick")
ax.legend()
plt.show()


datas = scattering.get_s1s2(
    wavelength=1_450 * ureg.nanometer,
    diameters=py_result.radii,
    phi=np.linspace(-np.pi / 2, np.pi / 2, 400),
    polarization=0 * ureg.degree,
)

datas.process()


phi, theta, phase_function = datas.get_phase_function(
    densities=py_result.densities,
    H=py_result.H,
    p=py_result.p,
    theta_points=150
)


fig1 = scattering.plottings.plot_phase_function_3d(
    phi=phi,
    theta=theta,
    phase_function=phase_function.to('1 / meter').magnitude,
    mode="spherical"
)

plt.show()